"""Tests de la Phase 2 : concentration par amorce, sortie JSON de analyze, commande redesign.
Phase 2 tests: per-primer concentration, analyze JSON output, redesign command.
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
    reverse_complement,
)
from labcraft.lamp.domains import _find_iupac_substring
from labcraft.optimize.redesign import (
    Locus,
    RedesignError,
    generate_candidates,
    geometry_satisfied,
)

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
GT_DIR = REPO_ROOT / "tests" / "fixtures" / "ground_truth"
DUPLEX_PRIMERS = GT_DIR / "duplex_YFV_DENV4_primers.txt"
HALF_CONC_ARGS = [
    "--temperature", "63.0",
    "--conc-fip-bip", "0.4",
    "--conc-loop", "0.2",
    "--conc-f3-b3", "0.1",
]


def _slice_around_amplicon(target_seq: str, primer_file: Path, pad: int = 100) -> str:
    """Extrait la fenêtre cible autour de l'amplicon pour accélérer les tests CLI.
    Slices target sequence tightly around the LAMP amplicon to speed up CLI tests.
    """
    recs = dict(read_multi_fasta(str(primer_file)))
    spans = []
    for name, seq in recs.items():
        if name.endswith(("FIP", "BIP")):
            continue
        for cand in (seq, reverse_complement(seq)):
            pos = _find_iupac_substring(cand, target_seq)
            if pos >= 0:
                spans.append((pos, pos + len(cand)))
    assert spans, f"Aucun domaine trouvé pour {primer_file.name}"
    lo = max(0, min(a for a, _ in spans) - pad)
    hi = min(len(target_seq), max(b for _, b in spans) + pad)
    return target_seq[lo:hi]


def _sum_amplifiable_families_nm(html: str) -> float:
    """Additionne les concentrations (en nM) des familles amplifiables affichées dans le HTML.
    Sums amplifiable family concentrations (in nM) from the HTML decomposition table.
    """
    vals = re.findall(r"<td>([0-9]+\.[0-9]{4})\s+nM</td>", html)
    return sum(float(v) for v in vals)


@pytest.fixture(scope="module")
def gt_targets_truncated(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """Cibles de vérité terrain tronquées autour des amplicons (singleplex et duplex).
    Truncated ground-truth target FASTA files around each amplicon (singleplex and duplex).
    """
    tdir = tmp_path_factory.mktemp("p2_targets")
    d3_id, d3_seq = read_multi_fasta(str(GT_DIR / "DEN3_M93130.fasta"))[0]
    d4_id, d4_seq = read_multi_fasta(str(GT_DIR / "Dengue_4_S1028.fasta"))[0]
    y_id, y_seq = read_multi_fasta(str(GT_DIR / "YFV_17D_17DD.fasta"))[0]

    d3w = _slice_around_amplicon(d3_seq, GT_DIR / "primer_dengue_3_panel3.txt")
    d4w = _slice_around_amplicon(d4_seq, GT_DIR / "parida_DENV4.txt")
    yw = _slice_around_amplicon(y_seq, GT_DIR / "nunes_YFV.txt")

    p_d3 = tdir / "DEN3.fasta"
    p_d4 = tdir / "DENV4.fasta"
    p_y = tdir / "YFV.fasta"
    p_dup = tdir / "DUPLEX.fasta"

    p_d3.write_text(f">{d3_id}\n{d3w}\n", encoding="utf-8")
    p_d4.write_text(f">{d4_id}\n{d4w}\n", encoding="utf-8")
    p_y.write_text(f">{y_id}\n{yw}\n", encoding="utf-8")
    p_dup.write_text(f">{y_id}\n{yw}\n>{d4_id}\n{d4w}\n", encoding="utf-8")
    return {"d3": p_d3, "d4": p_d4, "yfv": p_y, "dup": p_dup}


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


def test_fasta_header_conc_and_cli_priority(tmp_path: Path, gt_targets_truncated: dict[str, Path]) -> None:
    """FASTA conc= s'applique à l'amorce seule, la CLI --conc est prioritaire, les autres gardent leur classe.
    FASTA conc= applies to that primer only, CLI --conc wins, others keep their class value.
    """
    p_file = _write_panel(tmp_path, "conc=0.1")
    dup_targets = gt_targets_truncated["dup"]
    cfg = build_config_from_cli(
        primers_path=str(p_file), targets_path=str(dup_targets),
        conc_fip_bip=0.4, conc_f3_b3=0.1, conc_loop=0.2,
    )
    m = _conc_map(cfg)
    assert m["FIP_Nunes_YFV"] == pytest.approx(0.1)
    assert m["BIP_Nunes_YFV"] == pytest.approx(0.4)
    assert m["FIP_ParidaD4"] == pytest.approx(0.4)
    assert m["F3_ParidaD4"] == pytest.approx(0.1)

    cfg2 = build_config_from_cli(
        primers_path=str(p_file), targets_path=str(dup_targets),
        conc_fip_bip=0.4, conc_f3_b3=0.1, conc_loop=0.2,
        primer_concs={"Nunes_YFV_FIP": 0.05, "Nunes_YFV_BIP": 0.1},
    )
    m2 = _conc_map(cfg2)
    assert m2["FIP_Nunes_YFV"] == pytest.approx(0.05)
    assert m2["BIP_Nunes_YFV"] == pytest.approx(0.1)
    assert m2["BIP_ParidaD4"] == pytest.approx(0.4)

    with pytest.raises(ParseError):
        build_config_from_cli(
            primers_path=str(p_file), targets_path=str(dup_targets),
            primer_concs={"Inconnu_FIP": 0.1},
        )


@pytest.fixture(scope="module")
def quarter_fip_bip_run(
    gt_targets_truncated: dict[str, Path],
    tmp_path_factory: pytest.TempPathFactory,
) -> tuple[int, str, dict]:
    """Duplex demi-concentration avec Nunes_YFV_FIP et Nunes_YFV_BIP à 0,1 µM.
    Half-concentration duplex with Nunes_YFV_FIP and Nunes_YFV_BIP at 0.1 uM.
    """
    out = tmp_path_factory.mktemp("quarter") / "dup_quarter.html"
    res = CliRunner().invoke(app, [
        "analyze", "-p", str(DUPLEX_PRIMERS), "-t", str(gt_targets_truncated["dup"]), *HALF_CONC_ARGS,
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
    # Référence demi-concentration : 0,29 nM (Phase 1) / Half-conc reference: 0.29 nM (Phase 1)
    assert fam[key] < 0.29
    assert js["totals"]["amplifiable_dimers_total_nm"] < 0.3514
    assert js["totals"]["amplifiable_dimers_inter_nm"] < 0.3097
    assert set(js["initiation_access_by_panel"]).issuperset({"U17066.1", "MG272272.1"})
    assert set(js["sites_by_target"]).issuperset({"U17066.1", "MG272272.1"})


# --- 2.2 JSON de analyze / analyze JSON -------------------------------------------------

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
    gt_targets_truncated: dict[str, Path],
    tmp_path: Path,
) -> None:
    """Sur les 5 cas de non-régression, les totaux du JSON égalent ceux du HTML et les valeurs restent inchangées.
    Across all 5 regression cases, JSON totals match HTML totals and reference values remain unchanged.
    """
    out_html = tmp_path / f"{case_id}.html"
    res = CliRunner().invoke(
        app,
        [
            "analyze",
            "-p", str(GT_DIR / primers_file),
            "-t", str(gt_targets_truncated[target_key]),
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
    assert abs(js_total_nm - html_total_nm) < 1e-4, f"{case_id}: JSON {js_total_nm} vs HTML {html_total_nm}"
    assert (js["verdict"] == "FAILURE") == expect_failure
    if ref_nm == 0.0:
        assert js_total_nm == 0.0
    else:
        expected_banner = f"{js_total_nm:.4f}" if js_total_nm < 0.01 else f"{js_total_nm:.2f}"
        m_banner = re.search(r"Total dimères amplifiables :</strong>\s*([0-9.]+)\s+nM", html)
        assert m_banner is not None and m_banner.group(1) == expected_banner
        assert abs(js_total_nm - ref_nm) <= max(0.01 * ref_nm, 5e-5)

    if js["totals"]["is_multiplex"]:
        html_inter = re.search(r"Part inter-panels</strong>.*?<strong>([0-9.]+) nM</strong>", html, re.S)
        assert html_inter is not None
        assert f"{js['totals']['amplifiable_dimers_inter_nm']:.2f}" == html_inter.group(1)


# --- 2.3 Redessin / Redesign ------------------------------------------------------------

def test_redesign_refuses_fip_and_bip(gt_targets_truncated: dict[str, Path], tmp_path: Path) -> None:
    """Une demande sur une amorce interne (FIP ou BIP) renvoie le message de refus.
    A request on an inner primer (FIP or BIP) returns the refusal message.
    """
    for pname in ("Nunes_YFV_FIP", "ParidaD4_BIP"):
        res = CliRunner().invoke(app, [
            "redesign", "-p", str(DUPLEX_PRIMERS), "-t", str(gt_targets_truncated["dup"]),
            "--primer", pname, "-o", str(tmp_path / f"{pname}.html"),
        ])
        assert res.exit_code == 1
        assert "Redessin refusé" in res.stdout
        assert "FIP/BIP" in res.stdout


def test_redesign_geometry_and_candidates_unit() -> None:
    """Contraintes géométriques F3/B3/LF/LB et génération de candidats sur le même brin.
    Geometric constraints for F3/B3/LF/LB and candidate generation on the same strand.
    """
    geo = {
        "F2": Locus(120, 140, "+"),
        "F1": Locus(160, 180, "+"),
        "B1": Locus(200, 220, "-"),
        "B2": Locus(240, 260, "-"),
    }
    assert geometry_satisfied("F3", Locus(95, 115, "+"), geo)
    assert not geometry_satisfied("F3", Locus(105, 125, "+"), geo)  # chevauche F2 / overlaps F2
    assert geometry_satisfied("B3", Locus(265, 285, "-"), geo)
    assert not geometry_satisfied("B3", Locus(250, 270, "-"), geo)  # chevauche B2 / overlaps B2
    assert geometry_satisfied("LF", Locus(142, 158, "-"), geo)
    assert not geometry_satisfied("LF", Locus(150, 165, "-"), geo)  # chevauche F1 / overlaps F1
    assert geometry_satisfied("LB", Locus(222, 238, "+"), geo)
    assert not geometry_satisfied("LB", Locus(215, 235, "+"), geo)  # chevauche B1 / overlaps B1

    with pytest.raises(RedesignError):
        geometry_satisfied("F3", Locus(95, 115, "+"), {"F2": None, "F1": None, "B2": None, "B1": None})

    target = "A" * 10 + "CCCCCCCCGG" + "N" * 5 + "TTTTTTTT"
    cands = generate_candidates("CCCCCCCCGG", Locus(10, 20, "+"), target, shift=3, length_min=8, length_max=10)
    assert cands[0].is_reference and cands[0].shift == 0 and cands[0].sequence == "CCCCCCCCGG"
    assert all(c.locus.strand == "+" for c in cands)
    assert all("N" not in c.sequence for c in cands)


def test_redesign_parida_f3_reference_and_family(gt_targets_truncated: dict[str, Path], tmp_path: Path) -> None:
    """La ligne de référence reproduit 0,351 nM et un candidat supprime la famille FIP_Nunes_YFV × F3_ParidaD4.
    The reference row reproduces 0.351 nM and a candidate eliminates the FIP_Nunes_YFV x F3_ParidaD4 family.

    Plage réduite (--shift 4, --top 2) pour borner la durée en CI.
    Reduced range (--shift 4, --top 2) to bound CI runtime.
    """
    out = tmp_path / "redesign.html"
    res = CliRunner().invoke(app, [
        "redesign", "-p", str(DUPLEX_PRIMERS), "-t", str(gt_targets_truncated["dup"]), *HALF_CONC_ARGS,
        "--primer", "ParidaD4_F3", "--shift", "4", "--top", "2", "-o", str(out),
    ])
    assert res.exit_code == 0, res.stdout
    js = json.loads(out.with_suffix(".json").read_text(encoding="utf-8"))
    ref = js["ranked"][0]
    assert ref["is_reference"] and ref["shift"] == 0 and ref["length"] == len(js["original_sequence"])
    assert abs(ref["equilibrium"]["amplifiable_total_nm"] - 0.3514) / 0.3514 < 0.01
    assert any("FIP_Nunes_YFV" in f and "F3_ParidaD4" in f for f in ref["equilibrium"]["families"])
    for c in js["ranked"][1:]:
        assert c["geometry_ok"]
        assert c["strand"] == js["original_locus"]["strand"]
        assert "site_occupations" in c["equilibrium"]
    assert any(
        not any("FIP_Nunes_YFV" in f and "F3_ParidaD4" in f for f in c["equilibrium"]["families"])
        for c in js["ranked"][1:]
    )
    assert 'id="redesign-ranking"' in out.read_text(encoding="utf-8")
