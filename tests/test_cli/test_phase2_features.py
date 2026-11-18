"""Tests des Phases 2 et 3 : concentration par amorce, JSON de analyze, redesign, et chargement unique ViennaRNA.
Phase 2 and 3 tests: per-primer concentration, analyze JSON, redesign ranking/F2/locus, and single-load ViennaRNA.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

from labcraft.cli.main import app
from labcraft.cli.parsers import (
    ParseError,
    build_config_from_cli,
    parse_primer_conc_args,
    read_multi_fasta,
)
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.optimize.redesign import (
    Candidate,
    Locus,
    RedesignError,
    compute_f2_comparison,
    compute_max_3p_complementarity,
    equilibrium_sort_key,
    generate_candidates,
    geometry_satisfied,
    prescreen_sort_key,
    select_top_candidates,
)

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
GT_DIR = REPO_ROOT / "tests" / "fixtures" / "ground_truth"
DUPLEX_PRIMERS = GT_DIR / "duplex_YFV_DENV4_primers.txt"
DUPLEX_TARGETS = GT_DIR / "duplex_YFV_DENV4_targets.fasta"
HALF_CONC_ARGS = [
    "--temperature", "63.0",
    "--conc-fip-bip", "0.4",
    "--conc-loop", "0.2",
    "--conc-f3-b3", "0.1",
]

GT_TARGETS_MAP = {
    "d3": GT_DIR / "DEN3_M93130.fasta",
    "d4": GT_DIR / "Dengue_4_S1028.fasta",
    "yfv": GT_DIR / "YFV_17D_17DD.fasta",
    "dup": DUPLEX_TARGETS,
}


def _sum_amplifiable_families_nm(html: str) -> float:
    """Additionne les concentrations (en nM) des familles amplifiables affichées dans le HTML.
    Sums amplifiable family concentrations (in nM) from the HTML decomposition table.
    """
    vals = re.findall(r"<td>([0-9]+\.[0-9]{4})\s+nM</td>", html)
    return sum(float(v) for v in vals)


# --- 2.1 Concentration par amorce / Per-primer concentration ---------------------------

def test_parse_primer_conc_args_valid_and_invalid() -> None:
    """Analyse de --conc NOM=VALEUR et rejet des formats invalides.
    Parsing of --conc NAME=VALUE and rejection of invalid formats.
    """
    assert parse_primer_conc_args(["Nunes_YFV_FIP=0.1", "ParidaD4_F3=0.05"]) == {
        "Nunes_YFV_FIP": 0.1,
        "ParidaD4_F3": 0.05,
    }
    assert parse_primer_conc_args(None) == {}
    for bad in (["Nunes_YFV_FIP"], ["Nunes_YFV_FIP=abc"], ["=0.1"], ["X_F3=-1"]):
        with pytest.raises(ParseError):
            parse_primer_conc_args(bad)


def _write_panel(tmp_path: Path, fip_header_extra: str = "") -> Path:
    """Écrit le duplex avec un attribut conc= optionnel sur la FIP fièvre jaune.
    Writes the duplex with an optional conc= attribute on the yellow fever FIP.
    """
    out = []
    for h, s in read_multi_fasta(str(DUPLEX_PRIMERS)):
        name = h.split()[0]
        extra = f" {fip_header_extra}" if (name == "Nunes_YFV_FIP" and fip_header_extra) else ""
        out.append(f">{name}{extra}\n{s}\n")
    p = tmp_path / "duplex_conc.txt"
    p.write_text("".join(out), encoding="utf-8")
    return p


def _conc_map(cfg) -> dict:
    return {
        f"{role}_{ps.panel_name}": pc.conc_uM for ps in cfg.primer_sets for role, pc in ps.primers.items()
    }


def test_fasta_header_conc_and_cli_priority(tmp_path: Path) -> None:
    """FASTA conc= s'applique à l'amorce seule, la CLI --conc est prioritaire, les autres gardent leur classe.
    FASTA conc= applies to that primer only, CLI --conc wins, others keep their class value.
    """
    p_file = _write_panel(tmp_path, "conc=0.1")
    cfg = build_config_from_cli(
        primers_path=str(p_file), targets_path=str(DUPLEX_TARGETS),
        conc_fip_bip=0.4, conc_f3_b3=0.1, conc_loop=0.2,
    )
    m = _conc_map(cfg)
    assert m["FIP_Nunes_YFV"] == pytest.approx(0.1)
    assert m["BIP_Nunes_YFV"] == pytest.approx(0.4)
    assert m["FIP_ParidaD4"] == pytest.approx(0.4)
    assert m["F3_ParidaD4"] == pytest.approx(0.1)

    cfg2 = build_config_from_cli(
        primers_path=str(p_file), targets_path=str(DUPLEX_TARGETS),
        conc_fip_bip=0.4, conc_f3_b3=0.1, conc_loop=0.2,
        primer_concs={"Nunes_YFV_FIP": 0.05, "Nunes_YFV_BIP": 0.1},
    )
    m2 = _conc_map(cfg2)
    assert m2["FIP_Nunes_YFV"] == pytest.approx(0.05)
    assert m2["BIP_Nunes_YFV"] == pytest.approx(0.1)
    assert m2["BIP_ParidaD4"] == pytest.approx(0.4)

    with pytest.raises(ParseError):
        build_config_from_cli(
            primers_path=str(p_file), targets_path=str(DUPLEX_TARGETS),
            primer_concs={"Inconnu_FIP": 0.1},
        )


@pytest.fixture(scope="module")
def quarter_fip_bip_run(tmp_path_factory: pytest.TempPathFactory) -> tuple[int, str, dict]:
    """Duplex demi-concentration avec Nunes_YFV_FIP et Nunes_YFV_BIP à 0,1 µM sur cibles complètes.
    Half-concentration duplex with Nunes_YFV_FIP and Nunes_YFV_BIP at 0.1 uM on full targets.
    """
    out = tmp_path_factory.mktemp("quarter") / "dup_quarter.html"
    res = CliRunner().invoke(app, [
        "analyze", "-p", str(DUPLEX_PRIMERS), "-t", str(DUPLEX_TARGETS), *HALF_CONC_ARGS,
        "--conc", "Nunes_YFV_FIP=0.1", "--conc", "Nunes_YFV_BIP=0.1", "-o", str(out),
    ])
    js = json.loads(out.with_suffix(".json").read_text(encoding="utf-8")) if out.with_suffix(".json").exists() else {}
    html = out.read_text(encoding="utf-8") if out.exists() else ""
    return res.exit_code, res.stdout + html, js


def test_quarter_fip_bip_reduces_dominant_family(quarter_fip_bip_run: tuple[int, str, dict]) -> None:
    """La famille FIP_Nunes_YFV × F3_ParidaD4 diminue sous 0,29 nM ; concentrations effectives affichées.
    The FIP_Nunes_YFV x F3_ParidaD4 family decreases below 0.29 nM; effective concentrations displayed.
    """
    exit_code, text, js = quarter_fip_bip_run
    assert exit_code == 0, text
    assert js["primer_concentrations_uM"]["FIP_Nunes_YFV"] == pytest.approx(0.1)
    assert js["primer_concentrations_uM"]["BIP_Nunes_YFV"] == pytest.approx(0.1)
    assert js["primer_concentrations_uM"]["FIP_ParidaD4"] == pytest.approx(0.4)
    assert "FIP_Nunes_YFV=0.1" in text
    assert 'id="primer-concentrations"' in text

    fam = {f["label"]: f["conc_nm"] for f in js["totals"]["families"]}
    key = next(k for k in fam if "FIP_Nunes_YFV" in k and "F3_ParidaD4" in k)
    assert fam[key] < 0.29
    assert js["totals"]["amplifiable_dimers_total_nm"] < 0.3514
    assert js["totals"]["amplifiable_dimers_inter_nm"] < 0.3097
    assert set(js["initiation_access_by_panel"]).issuperset({"U17066.1", "MG272272.1"})
    assert set(js["sites_by_target"]).issuperset({"U17066.1", "MG272272.1"})


# --- 2.2 JSON de analyze & non-régression au centième / analyze JSON & regression -------

def test_analyze_json_schema(quarter_fip_bip_run: tuple[int, str, dict]) -> None:
    """Le JSON contient toutes les rubriques demandées.
    The JSON contains all requested sections.
    """
    _exit, _text, js = quarter_fip_bip_run
    for key in (
        "conditions", "version", "inputs", "primer_concentrations_uM", "free_fractions",
        "sites_by_target", "dimers", "totals", "mispriming_risks", "verdict", "verdict_cause",
    ):
        assert key in js, key
    assert js["inputs"]["primers"]["sha256"]
    site = next(iter(next(iter(js["sites_by_target"].values())).values()))
    for k in ("position", "strand", "mismatches", "dg_hyb", "dg_open", "occupation"):
        assert k in site
    d = js["dimers"][0]
    for k in ("name", "conc_nm", "dg_global", "dg_3prime", "amplifiable", "scope"):
        assert k in d
    amp_sum = sum(x["conc_nm"] for x in js["dimers"] if x["amplifiable"])
    assert amp_sum == pytest.approx(js["totals"]["amplifiable_dimers_total_nm"], rel=1e-6, abs=1e-9)


@pytest.mark.parametrize(
    ("case_id", "primers_file", "target_key", "c_fip", "c_loop", "c_f3", "ref_nm", "expect_failure"),
    [
        ("singleplex_denv3", "primer_dengue_3_panel3.txt", "d3", "0.8", "0.4", "0.2", 0.0, False),
        ("singleplex_denv4", "parida_DENV4.txt", "d4", "0.8", "0.4", "0.2", 0.1540, False),
        ("singleplex_yfv", "nunes_YFV.txt", "yfv", "0.8", "0.4", "0.2", 0.0011, False),
        ("duplex_half", "duplex_YFV_DENV4_primers.txt", "dup", "0.4", "0.2", "0.1", 0.3514, True),
        ("duplex_full", "duplex_YFV_DENV4_primers.txt", "dup", "0.8", "0.4", "0.2", 1.3498, True),
    ],
)
def test_analyze_json_matches_html_five_cases(
    case_id: str,
    primers_file: str,
    target_key: str,
    c_fip: str,
    c_loop: str,
    c_f3: str,
    ref_nm: float,
    expect_failure: bool,
    tmp_path: Path,
) -> None:
    """Sur les 5 cas de non-régression (cibles complètes), JSON == HTML et valeurs inchangées au centième.
    Across all 5 regression cases (full targets), JSON matches HTML and reference values hold within 0.01.
    """
    out_html = tmp_path / f"{case_id}.html"
    res = CliRunner().invoke(
        app,
        [
            "analyze",
            "-p", str(GT_DIR / primers_file),
            "-t", str(GT_TARGETS_MAP[target_key]),
            "--temperature", "63.0",
            "--conc-fip-bip", c_fip,
            "--conc-loop", c_loop,
            "--conc-f3-b3", c_f3,
            "-o", str(out_html),
        ],
    )
    assert res.exit_code == 0, res.stdout
    html = out_html.read_text(encoding="utf-8")
    js = json.loads(out_html.with_suffix(".json").read_text(encoding="utf-8"))

    html_total_nm = _sum_amplifiable_families_nm(html)
    js_total_nm = js["totals"]["amplifiable_dimers_total_nm"]
    # Borne exacte d'arrondi pour la somme de K familles affichées à 10^-4 : K * 0,5e-4
    # Exact rounding bound when summing K family cells formatted at 1e-4: K * 0.5e-4
    rounding_tol = max(1e-4, len(js["totals"]["families"]) * 5e-5)
    assert abs(js_total_nm - html_total_nm) <= rounding_tol, (
        f"{case_id}: JSON {js_total_nm} vs HTML {html_total_nm} (tol={rounding_tol})"
    )
    assert (js["verdict"] == "FAILURE") == expect_failure
    if ref_nm == 0.0:
        assert js_total_nm == 0.0
    else:
        expected_banner = f"{js_total_nm:.4f}" if js_total_nm < 0.01 else f"{js_total_nm:.2f}"
        m_banner = re.search(r"Total dimères amplifiables :</strong>\s*([0-9.]+)\s+nM", html)
        assert m_banner is not None and m_banner.group(1) == expected_banner
        assert abs(js_total_nm - ref_nm) < 0.01
        assert abs(js_total_nm - ref_nm) <= max(0.01 * ref_nm, 5e-5)

    if js["totals"]["is_multiplex"]:
        html_inter = re.search(r"Part inter-panels</strong>.*?<strong>([0-9.]+) nM</strong>", html, re.S)
        assert html_inter is not None
        assert f"{js['totals']['amplifiable_dimers_inter_nm']:.2f}" == html_inter.group(1)


# --- 2.3 Redessin : classement, F2, complémentarité 3', coordonnées réelles -------------

def test_redesign_refuses_fip_and_bip(tmp_path: Path) -> None:
    """Une demande sur une amorce interne (FIP ou BIP) renvoie le message de refus.
    A request on an inner primer (FIP or BIP) returns the refusal message.
    """
    for pname in ("Nunes_YFV_FIP", "ParidaD4_BIP"):
        res = CliRunner().invoke(app, [
            "redesign", "-p", str(DUPLEX_PRIMERS), "-t", str(DUPLEX_TARGETS),
            "--primer", pname, "-o", str(tmp_path / f"{pname}.html"),
        ])
        assert res.exit_code == 1
        assert "Redessin refusé" in res.stdout
        assert "FIP/BIP" in res.stdout


def test_redesign_geometry_and_ranking_unit() -> None:
    """Contraintes géométriques, comparaison F2, complémentarité 3' avec matrice 5' et départage.
    Geometric constraints, F2 stability rule, 5'-templated 3' complementarity, and tie-breaking.
    """
    geo = {
        "F2": Locus(120, 140, "+"),
        "F1": Locus(160, 180, "+"),
        "B1": Locus(200, 220, "-"),
        "B2": Locus(240, 260, "-"),
    }
    assert geometry_satisfied("F3", Locus(95, 115, "+"), geo)
    assert not geometry_satisfied("F3", Locus(105, 125, "+"), geo)
    assert geometry_satisfied("B3", Locus(265, 285, "-"), geo)
    assert not geometry_satisfied("B3", Locus(250, 270, "-"), geo)
    assert geometry_satisfied("LF", Locus(142, 158, "-"), geo)
    assert not geometry_satisfied("LF", Locus(150, 165, "-"), geo)
    assert geometry_satisfied("LB", Locus(222, 238, "+"), geo)
    assert not geometry_satisfied("LB", Locus(215, 235, "+"), geo)

    with pytest.raises(RedesignError):
        geometry_satisfied("F3", Locus(95, 115, "+"), {"F2": None, "F1": None, "B2": None, "B1": None})

    target = "A" * 10 + "CCCCCCCCGG" + "N" * 5 + "TTTTTTTT"
    cands = generate_candidates("CCCCCCCCGG", Locus(10, 20, "+"), target, shift=3, length_min=8, length_max=10)
    assert cands[0].is_reference and cands[0].shift == 0 and cands[0].sequence == "CCCCCCCCGG"
    assert all(c.locus.strand == "+" for c in cands)
    assert all("N" not in c.sequence for c in cands)

    # Vérifie l'exigence de matrice en 5' (position >= 1 dans le partenaire)
    # Checks the 5' template requirement (match at index >= 1 in partner, not flush at 5' end)
    p_flush_5p = PhysicalPrimer.from_simple("P1", "CCGGGGTTTT", PrimerRole.FIP)   # revcomp(CCCCGG)=CCGGGG en pos 0 -> pas de matrice 5'
    p_with_5p = PhysicalPrimer.from_simple("P2", "ATCCGGTTTT", PrimerRole.BIP)    # revcomp(CCGG)=CCGG en pos 2 -> matrice 5' 'AT'
    c_test = Candidate(0, 10, "AAAACCCCGG", Locus(0, 10, "+"))
    compute_max_3p_complementarity(c_test, [p_flush_5p, p_with_5p])
    assert c_test.max_3p_comp == 4
    assert c_test.max_3p_partner == "P2"

    # Ordre de départage : non >F2, puis comp 3', puis |L - L0|, puis GC dans [40, 60]
    # Tie-breaker order: not >F2, then 3' comp, then |L - L0|, then GC in [40, 60]
    dg_f2 = -12.91
    c_too_stable = Candidate(+2, 22, "GGCTATTGAAGTCAGGCCACTT", Locus(0, 22, "+"), dg_duplex=-13.87)
    c_comp4 = Candidate(-9, 18, "CCAAAGAGGCTATTGAAG", Locus(0, 18, "+"), dg_duplex=-9.36)
    c_len19_gc_ok = Candidate(-7, 19, "CAAAGAGGCTATTGAAGTC", Locus(0, 19, "+"), dg_duplex=-9.34)
    c_len19_gc_bad = Candidate(-8, 17, "AAAGAGGCTATTGAAGT", Locus(0, 17, "+"), dg_duplex=-8.78)
    c_best_len18 = Candidate(-7, 18, "AAAGAGGCTATTGAAGTC", Locus(0, 18, "+"), dg_duplex=-9.04)

    for c_item, k3, gc in [
        (c_too_stable, 3, 50.0),
        (c_comp4, 4, 50.0),
        (c_len19_gc_ok, 3, 42.1),
        (c_len19_gc_bad, 3, 35.3),
        (c_best_len18, 3, 38.9),
    ]:
        compute_f2_comparison(c_item, dg_f2)
        c_item.max_3p_comp = k3
        c_item.gc_percent = gc
        c_item.equilibrium = {"amplifiable_total_nm": 0.058569, "amplifiable_inter_nm": 0.0168}

    pool = [c_too_stable, c_comp4, c_len19_gc_bad, c_len19_gc_ok, c_best_len18]
    pool.sort(key=lambda c: prescreen_sort_key(c, orig_len=18))
    assert pool == [c_best_len18, c_len19_gc_ok, c_len19_gc_bad, c_comp4, c_too_stable]

    eq_sorted = sorted(pool, key=lambda c: equilibrium_sort_key(c, orig_len=18))
    assert eq_sorted == [c_best_len18, c_len19_gc_ok, c_len19_gc_bad, c_comp4, c_too_stable]

    kept, n_disc = select_top_candidates(pool, top=4, role="F3")
    assert n_disc == 1
    assert c_too_stable not in kept


def test_redesign_parida_f3_full_exploration(tmp_path: Path) -> None:
    """Exploration complète ParidaD4_F3 (107 candidats) : locus réel 10291-10309, filtre F2 et départage.
    Full ParidaD4_F3 exploration (107 candidates): true locus 10291-10309, F2 rule, and tie-breaking.
    """
    out = tmp_path / "redesign.html"
    res = CliRunner().invoke(app, [
        "redesign", "-p", str(DUPLEX_PRIMERS), "-t", str(DUPLEX_TARGETS), *HALF_CONC_ARGS,
        "--primer", "ParidaD4_F3", "--shift", "12", "--length", "17-22", "--top", "5", "-o", str(out),
    ])
    assert res.exit_code == 0, res.stdout
    js = json.loads(out.with_suffix(".json").read_text(encoding="utf-8"))
    html = out.read_text(encoding="utf-8")

    assert js["n_candidates_geometry_ok"] == 107
    assert js["original_locus"] == {"start": 10291, "end": 10309, "strand": "+"}
    assert js["dg_f2"] == pytest.approx(-12.9052, abs=0.01)
    assert js["n_discarded_more_stable_than_f2"] > 0

    ref = js["ranked"][0]
    assert ref["is_reference"] and ref["shift"] == 0 and ref["length"] == 18
    assert ref["start"] == 10291 and ref["end"] == 10309
    assert ref["max_3p_complementarity"] == 6
    assert abs(ref["equilibrium"]["amplifiable_total_nm"] - 0.3514) < 0.01
    assert any("FIP_Nunes_YFV" in f and "F3_ParidaD4" in f for f in ref["equilibrium"]["families"])

    recs = js["ranked"][1:]
    assert len(recs) == 5
    for c in recs:
        assert c["geometry_ok"]
        assert c["strand"] == "+"
        assert not c["more_stable_than_f2"]
        assert c["dg_duplex"] >= js["dg_f2"]
        assert c["ddg_vs_f2"] >= 0.0
        assert c["max_3p_complementarity"] == 3
        assert "site_occupations" in c["equilibrium"]
        assert not any("FIP_Nunes_YFV" in f and "F3_ParidaD4" in f for f in c["equilibrium"]["families"])

    # Départage : d'abord L=18 (|L - 18| == 0), puis |L - 18| == 1 avec GC entre 40 et 60 %
    # Tie-breaking: first L=18 (|L - 18| == 0), then |L - 18| == 1 with GC within [40%, 60%]
    assert [c["length"] for c in recs[:2]] == [18, 18]
    for c in recs[2:]:
        assert abs(c["length"] - 18) == 1
        assert 40.0 <= c["gc_percent"] <= 60.0

    assert 'id="redesign-ranking"' in html
    assert "10291-10309 (+)" in html
    assert "100-118" not in html
    assert "ΔG vs F2" in html
    assert "Comp. 3' max (nt)" in html


# --- Phase 3 : chargement unique des paramètres ADN ViennaRNA / Single DNA load ---------

def test_phase3_single_dna_params_load_per_execution(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Une exécution complète de `analyze` ne charge Mathews 2004 qu'une seule fois et restaure Turner 2004.
    A full `analyze` execution loads Mathews 2004 DNA params exactly once and restores Turner 2004 RNA at exit.
    """
    import RNA
    from labcraft.thermo.backends.vienna import _SHARED_VIENNA_CACHE
    from labcraft.thermo.backends.vienna_salt import _SHARED_SALT_CACHE
    from labcraft.target.cache import unfolding_cache

    _SHARED_VIENNA_CACHE.clear()
    _SHARED_SALT_CACHE.clear()
    unfolding_cache.clear()

    dna_loads = 0
    rna_loads = 0
    orig_dna = RNA.params_load_DNA_Mathews2004
    orig_rna = RNA.params_load_RNA_Turner2004

    def counted_dna() -> None:
        nonlocal dna_loads
        dna_loads += 1
        orig_dna()

    def counted_rna() -> None:
        nonlocal rna_loads
        rna_loads += 1
        orig_rna()

    orig_rna()
    RNA.cvar.temperature = 37.0
    ref_rna_energy = RNA.fold_compound("GGGGCCAAAAAAGGCCCC").pf()[1]

    monkeypatch.setattr(RNA, "params_load_DNA_Mathews2004", counted_dna)
    monkeypatch.setattr(RNA, "params_load_RNA_Turner2004", counted_rna)

    out_html = tmp_path / "d4_single_load.html"
    res = CliRunner().invoke(app, [
        "analyze",
        "-p", str(GT_DIR / "parida_DENV4.txt"),
        "-t", str(GT_DIR / "Dengue_4_S1028.fasta"),
        "--temperature", "63.0",
        "--conc-fip-bip", "0.8",
        "--conc-loop", "0.4",
        "--conc-f3-b3", "0.2",
        "-o", str(out_html),
    ])
    assert res.exit_code == 0, res.stdout
    assert dna_loads == 1
    assert rna_loads == 1

    # Vérifie que l'état ARN et la température ont bien été restaurés en fin d'exécution
    # Verify that RNA parameters and temperature are properly restored after execution
    assert RNA.cvar.temperature == 37.0
    post_rna_energy = RNA.fold_compound("GGGGCCAAAAAAGGCCCC").pf()[1]
    assert abs(post_rna_energy - ref_rna_energy) < 1e-4
