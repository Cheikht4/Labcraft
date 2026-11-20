"""Tests unitaires et d'intégration de l'Étape 5 : boucles internes/renflements et probabilités 3' de Boltzmann.

Unit and integration tests for Step 5: bulge/internal-loop structure thermodynamics in
`ViennaSaltShiftBackend` and Boltzmann 3' base-pairing probability (`pf_cofold`) in `amplifiable_dimer`.
"""
from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from labcraft.cli.main import app
from labcraft.diagnostics.amplifiable_dimer import (
    calc_boltzmann_3p_pairing,
    evaluate_pair_amplifiable,
)
from labcraft.diagnostics.enzyme import BST_2_0
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.thermo.backends.vienna_salt import (
    ViennaSaltShiftBackend,
    estimate_helix_thermo,
    estimate_structure_thermo,
)
from labcraft.thermo.salt import UnifiedSaltModel
from labcraft.thermo.vienna import dna_session


def test_bulge_and_internal_loop_thermo_gibbs_helmholtz_and_tm() -> None:
    """Vérifie l'auto-cohérence Gibbs-Helmholtz et la baisse de Tm sur duplexes à bulge et boucle interne.
    Verifies Gibbs-Helmholtz self-consistency and Tm drop on bulged and internal-looped duplexes.
    """
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())

    # Duplexe parfait de 16 pb vs même duplexe avec renflement central d'1 nt ('T')
    # Perfect 16-bp duplex vs same duplex with a central 1-nt ('T') bulge
    seq_top_perf = "ATGCATGCATGCATGC"
    seq_bot_perf = "GCATGCATGCATGCAT"
    seq_top_bulge1 = "ATGCATGCTATGCATGC"
    seq_top_bulge3 = "ATGCATGCTTTATGCATGC"
    seq_top_int2x2 = "ATGCATGCCAATGCATGC"
    seq_bot_int2x2 = "GCATGCAACGGCATGCAT"

    with dna_session(temp_celsius=65.0):
        res_perf = backend.calc_duplex(
            seq_top_perf, seq_bot_perf, temp_celsius=65.0, na_mm=50.0, mg_mm=8.0, dntp_mm=1.4
        )
        res_bulge1 = backend.calc_duplex(
            seq_top_bulge1, seq_bot_perf, temp_celsius=65.0, na_mm=50.0, mg_mm=8.0, dntp_mm=1.4
        )
        res_bulge3 = backend.calc_duplex(
            seq_top_bulge3, seq_bot_perf, temp_celsius=65.0, na_mm=50.0, mg_mm=8.0, dntp_mm=1.4
        )
        res_int2x2 = backend.calc_duplex(
            seq_top_int2x2, seq_bot_int2x2, temp_celsius=65.0, na_mm=50.0, mg_mm=8.0, dntp_mm=1.4
        )
        res_bulge1_no_mg = backend.calc_duplex(
            seq_top_bulge1, seq_bot_perf, temp_celsius=65.0, na_mm=50.0, mg_mm=0.0, dntp_mm=0.0
        )

    # 1. Le renflement déstabilise le duplexe (dG moins négatif et Tm plus bas que le duplexe parfait)
    # 1. The bulge destabilizes the duplex (less negative dG and lower Tm than perfect duplex)
    assert res_bulge1.dg_kcal > res_perf.dg_kcal
    assert res_bulge1.tm_celsius < res_perf.tm_celsius
    assert res_bulge3.tm_celsius < res_perf.tm_celsius

    # 2. Cohérence thermodynamique Gibbs-Helmholtz sur structures à boucles : dG(T) == dH - T * dS
    # 2. Gibbs-Helmholtz thermodynamic consistency on loop/bulge structures: dG(T) == dH - T * dS
    temp_k = 65.0 + 273.15
    for r in (res_bulge1, res_bulge3, res_int2x2):
        dg_reconstructed = r.dh_kcal - temp_k * (r.ds_cal_per_k / 1000.0)
        assert abs(dg_reconstructed - r.dg_kcal) < 0.03
    assert abs((res_perf.dh_kcal - temp_k * (res_perf.ds_cal_per_k / 1000.0)) - res_perf.dg_kcal) < 0.25

    # 3. Stabilisation monotone par Mg2+ sur un duplexe à renflement
    # 3. Monotonic Mg2+ stabilization on a bulged duplex
    assert res_bulge1.dg_kcal < res_bulge1_no_mg.dg_kcal
    assert res_bulge1.tm_celsius > res_bulge1_no_mg.tm_celsius

    # 4. Vérification directe de `estimate_structure_thermo` vs `estimate_helix_thermo` sur bulge 1-nt
    # 4. Direct check of `estimate_structure_thermo` vs `estimate_helix_thermo` on 1-nt bulge
    dh_wc, ds_wc, n_bp, _, _, _ = estimate_helix_thermo(
        f"{seq_top_bulge1}&{seq_bot_perf}", res_bulge1.structure
    )
    dh_st, ds_st, n_bp_st, _, _, _ = estimate_structure_thermo(
        f"{seq_top_bulge1}&{seq_bot_perf}", res_bulge1.structure
    )
    assert n_bp == n_bp_st == 16
    assert (dh_st, ds_st) != (dh_wc, ds_wc)


def test_boltzmann_3p_pairing_and_suboptimal_masked_dimer() -> None:
    """Vérifie `calc_boltzmann_3p_pairing` et la détection d'une structure 3' sous-optimale masquée par le MFE.
    Verifies `calc_boltzmann_3p_pairing` and detection of a suboptimal 3' dimer masked by the MFE.
    """
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    salt_kw = {"na_mm": 50.0, "mg_mm": 8.0, "dntp_mm": 1.4}

    with dna_session(temp_celsius=63.0):
        # 1. Dimère fortement apparié en 3' vs dimère uniquement apparié en 5'
        # 1. Strongly 3'-paired dimer vs 5'-only paired dimer
        pf_strong = calc_boltzmann_3p_pairing("TTTTTTTGCATGCAT", "AAAAAAAATGCATGC", temp_celsius=63.0)
        pf_5p_only = calc_boltzmann_3p_pairing("GCATGCATAAAAAAA", "ATGCATGCTTTTTTT", temp_celsius=63.0)
        assert pf_strong["boltzmann_3p_prob"] > 0.70
        assert pf_5p_only["boltzmann_3p_prob"] < 0.05

        # 2. Le blocage chimique 3' annule la probabilité d'extension sans annuler l'appariement physique
        # 2. Chemical 3' blocking zeroes extensible probability while keeping physical pairing probability
        pf_blocked = calc_boltzmann_3p_pairing(
            "TTTTTTTGCATGCAT", "AAAAAAAATGCATGC", temp_celsius=63.0, blocked_a=True, blocked_b=True
        )
        assert pf_blocked["prob_3p_paired_a"] > 0.50
        assert pf_blocked["prob_3p_extensible_a"] == 0.0
        assert pf_blocked["boltzmann_3p_prob"] == 0.0

        # 3. Conformation 3' sous-optimale masquée par une tige 5' GC légèrement plus stable au MFE
        # 3. Suboptimal 3' conformation masked by a slightly more stable 5' GC stem in the single MFE
        p_sub_a = PhysicalPrimer("Sub_A", "CGCGCGCAAAAAATGCATGC", PrimerRole.F3, "CGCGCGCAAAAAATGCATGC")
        p_sub_b = PhysicalPrimer("Sub_B", "CGCGCGCTTTTTGCATGCAT", PrimerRole.B3, "CGCGCGCTTTTTGCATGCAT")
        is_amp, _dg3p, det = evaluate_pair_amplifiable(
            p_sub_a,
            p_sub_b,
            backend,
            BST_2_0,
            63.0,
            include_boltzmann_3p=True,
            suboptimal_3p_prob_threshold=0.10,
            **salt_kw,
        )
        assert is_amp is False
        assert det["boltzmann_3p_prob"] > 0.12
        assert det["suboptimal_3p_masked"] is True


def test_analyze_json_exports_boltzmann_3p_fields(tmp_path: Path) -> None:
    """Vérifie que `labcraft analyze` exporte `boltzmann_3p_prob` et `suboptimal_3p_masked` dans le JSON.
    Verifies `labcraft analyze` exports `boltzmann_3p_prob` and `suboptimal_3p_masked` in the JSON output.
    """
    primers_file = tmp_path / "primers.txt"
    target_file = tmp_path / "target.fasta"
    out_html = tmp_path / "report.html"

    target_seq = (
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "CGCGCGCAAAAAATGCATGC"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "ATGCATGCAAAAAGCGCGCG"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    )
    target_file.write_text(f">T1\n{target_seq}\n", encoding="utf-8")
    primers_file.write_text(
        ">SubPanel_F3\nCGCGCGCAAAAAATGCATGC\n"
        ">SubPanel_B3\nCGCGCGCTTTTTGCATGCAT\n",
        encoding="utf-8",
    )

    runner = CliRunner()
    res = runner.invoke(
        app,
        [
            "analyze",
            "-p", str(primers_file),
            "-t", str(target_file),
            "--temperature", "63.0",
            "--conc-f3-b3", "0.8",
            "-o", str(out_html),
        ],
    )
    assert res.exit_code == 0, res.stdout
    data = json.loads((tmp_path / "report.json").read_text(encoding="utf-8"))
    assert len(data["dimers"]) >= 1
    hetero = next(d for d in data["dimers"] if "F3" in d["name"] and "B3" in d["name"])
    assert hetero["boltzmann_3p_prob"] is not None
    assert hetero["boltzmann_3p_prob"] > 0.10
    assert hetero["suboptimal_3p_masked"] is True
    assert hetero["prob_3p_intra_a"] is not None
    assert hetero["prob_3p_intra_b"] is not None
    html_text = out_html.read_text(encoding="utf-8")
    assert "heuristique non calibré" in html_text
    assert "sans effet sur le verdict" in html_text
    assert "validé expérimentalement" in html_text
    readme_text = (Path(__file__).resolve().parent.parent.parent / "README.md").read_text(encoding="utf-8")
    assert "ne doit pas entrer dans le verdict tant qu'il n'est pas validé expérimentalement" in readme_text


def test_point3_boltzmann_3p_inter_vs_intra_contiguous_and_connected() -> None:
    """Vérifie la définition corrigée de `calc_boltzmann_3p_pairing` (Point 3) :
    (a) rejet des contacts isolés d'1 seule paire en 3' (`min_3p_bp >= 2`, défaut 3) et respect du veto ARMS ;
    (b) séparation stricte entre épingle auto-amorçante intramoléculaire (`prob_3p_intra_*`) et dimère A:B ;
    (c) conditionnement exact sur l'existence d'au moins une paire intermoléculaire entre A et B (`FAB`).

    Verifies the corrected definition of `calc_boltzmann_3p_pairing` (Point 3):
    (a) single-base 3' contacts rejected (`min_3p_bp >= 2`, default 3) and ARMS veto enforced;
    (b) strict separation of intramolecular self-priming hairpins (`prob_3p_intra_*`) from A:B dimers;
    (c) exact conditioning on >= 1 intermolecular pair between A and B (`FAB`).
    """
    import math
    import pytest
    import RNA
    from labcraft.diagnostics.enzyme import get_enzyme
    from labcraft.thermo.constants import R_GAS_KCAL_MOL_K

    with pytest.raises(ValueError, match="min_3p_bp"):
        calc_boltzmann_3p_pairing("ATGCATGC", "GCATGCAT", min_3p_bp=1)

    with dna_session(temp_celsius=63.0):
        # (a1) Contact 3' d'une seule base malgré une forte tige 5' : rejeté (prob == 0.0)
        # (a1) Lone 1-bp 3' contact alongside a strong 5' stem: rejected (prob == 0.0)
        pf_1bp = calc_boltzmann_3p_pairing("GCATGCATAAAAAAC", "ATGCATGCTGTGAAA", temp_celsius=63.0)
        assert pf_1bp["boltzmann_3p_prob"] == 0.0

        # (a2) Wobble G·T en position pénultième (-2) : bloqué en ARMS strict, autorisé si activé sur l'enzyme
        # (a2) Penultimate (-2) G·T wobble: blocked under strict ARMS, allowed when enabled on enzyme
        # A 3' = ... G(-3) T(-2) C(-1) ; B 5'->3' = AAAA G(+1 vs C) G(+2 wobble vs T) C(+3 vs G) AAAA
        pf_wob_strict = calc_boltzmann_3p_pairing(
            "AAAAAAAAAGTC", "AAAAGGCAAAAA", temp_celsius=50.0, enzyme=BST_2_0
        )
        bst_wob = get_enzyme({"name": "bst2.0", "allow_penultimate_gt_wobble": True})
        pf_wob_loose = calc_boltzmann_3p_pairing(
            "AAAAAAAAAGTC", "AAAAGGCAAAAA", temp_celsius=50.0, enzyme=bst_wob
        )
        assert pf_wob_strict["prob_3p_inter_a"] == 0.0
        assert pf_wob_loose["prob_3p_inter_a"] > 0.10

        # (b) Épingle auto-amorçante sur A face à un brin B sans complémentarité 3' :
        #     `boltzmann_3p_prob` (intermoléculaire) vaut 0.0 tandis que `prob_3p_intra_a` > 0.80
        # (b) Self-priming hairpin on A paired with strand B lacking 3' complementarity:
        #     intermolecular `boltzmann_3p_prob` is 0.0 while `prob_3p_intra_a` > 0.80
        pf_hp = calc_boltzmann_3p_pairing("AAAAAGCATGCTTTTGCATGC", "CCCCCCCCCCCCCCC", temp_celsius=63.0)
        assert pf_hp["prob_3p_inter_a"] == 0.0
        assert pf_hp["prob_3p_inter_b"] == 0.0
        assert pf_hp["boltzmann_3p_prob"] == 0.0
        assert pf_hp["prob_3p_intra_a"] > 0.80
        assert pf_hp["prob_3p_intra_b"] == 0.0

        # (c1) Deux brins sans aucune paire intermoléculaire possible : `connected_dimer_dg_kcal` est None
        # (c1) Two strands with zero possible intermolecular pairs: `connected_dimer_dg_kcal` is None
        pf_unb = calc_boltzmann_3p_pairing("AAATATAATTTTTATTATAAA", "CCCCCCCCCCCCCCC", temp_celsius=37.0)
        assert pf_unb["connected_dimer_dg_kcal"] is None
        assert pf_unb["boltzmann_3p_prob"] == 0.0

        # (c2) Accord analytique exact avec l'énumération exhaustive des structures connectées sur un 8-mer
        # (c2) Exact analytical agreement against brute-force enumeration of connected structures on an 8-mer
        sa, sb = "AAAAAGGC", "AAAGCCAA"
        pf_8 = calc_boltzmann_3p_pairing(sa, sb, temp_celsius=63.0, min_3p_bp=3)
        fc_ref = RNA.fold_compound(f"{sa}&{sb}")
        rt = R_GAS_KCAL_MOL_K * (63.0 + 273.15)
        seq_cat = sa + sb
        allowed_pairs = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"), ("G", "T"), ("T", "G")}
        memo: dict[tuple[int, int], list[str]] = {}

        def _enum(i: int, j: int) -> list[str]:
            if i > j:
                return [""]
            if (i, j) in memo:
                return memo[(i, j)]
            out = ["." + s for s in _enum(i + 1, j)]
            for k in range(i + 1, j + 1):
                same = (i < 8 and k < 8) or (i >= 8 and k >= 8)
                if same and (k - i - 1) < 3:
                    continue
                if (seq_cat[i], seq_cat[k]) in allowed_pairs:
                    for inside in _enum(i + 1, k - 1):
                        for right in _enum(k + 1, j):
                            out.append("(" + inside + ")" + right)
            memo[(i, j)] = out
            return out

        z_conn = 0.0
        z_3p = 0.0
        for st in _enum(0, len(seq_cat) - 1):
            stk: list[int] = []
            prs: dict[int, int] = {}
            for idx1, ch in enumerate(st, 1):
                if ch == "(":
                    stk.append(idx1)
                elif ch == ")":
                    j1 = stk.pop()
                    prs[j1] = idx1
                    prs[idx1] = j1
            if any(u <= 8 < v for u, v in prs.items()):
                w = math.exp(-fc_ref.eval_structure(st) / rt)
                z_conn += w
                if prs.get(8) == 12 and prs.get(7) == 13 and prs.get(6) == 14:
                    z_3p += w
        assert pf_8["prob_3p_inter_a"] == pytest.approx(z_3p / z_conn, abs=0.005)


def test_point2_santalucia_hicks_2004_table4_exact_values_and_default_equilibrium_dg() -> None:
    """Vérifie sans dépendance externe les tables SantaLucia & Hicks (2004) Table 4 & Table 1 et la
    non-modification du ΔG d'équilibre par défaut (`couple_loops_in_salt_dg=False`).

    Verifies without external GPL dependencies the SantaLucia & Hicks (2004) Table 4 & Table 1 constants
    and that default `couple_loops_in_salt_dg=False` leaves equilibrium ΔG untouched.
    """
    import pytest

    from labcraft.cli.config import build_engine_from_config
    from labcraft.cli.parsers import build_config_from_cli
    from labcraft.thermo.backends.vienna_salt import (
        _BULGE_LOOP_TABLE4_DG37_KCAL,
        _TERMINAL_AT_PENALTY,
    )
    from labcraft.thermo.mismatch import (
        _INTERNAL_LOOP_DG37_KCAL,
        _INTERNAL_LOOP_TABLE4_DG37_KCAL,
        _JACOBSON_STOCKMAYER_EXP,
    )

    expected_internal_table4 = {
        3: 3.2, 4: 3.6, 5: 4.0, 6: 4.4, 7: 4.6, 8: 4.8, 9: 4.9, 10: 4.9,
        12: 5.2, 14: 5.4, 15: 5.5, 16: 5.6, 18: 5.8, 20: 5.9, 25: 6.3, 30: 6.6,
    }
    expected_bulge_table4 = {
        1: 4.0, 2: 2.9, 3: 3.1, 4: 3.2, 5: 3.3, 6: 3.5, 7: 3.7, 8: 3.9,
        9: 4.1, 10: 4.3, 12: 4.5, 14: 4.8, 15: 4.9, 16: 5.0, 18: 5.2,
        20: 5.3, 25: 5.6, 30: 5.9,
    }
    assert _INTERNAL_LOOP_TABLE4_DG37_KCAL == expected_internal_table4
    assert _BULGE_LOOP_TABLE4_DG37_KCAL == expected_bulge_table4
    assert _JACOBSON_STOCKMAYER_EXP == pytest.approx(2.44)
    assert _TERMINAL_AT_PENALTY == (2.2, 6.9)
    for sz, val in expected_internal_table4.items():
        assert _INTERNAL_LOOP_DG37_KCAL[sz] == pytest.approx(val, abs=1e-9)

    cfg = build_config_from_cli(config_path="examples/denv2_lopez.yaml")
    _, _, eng_backend, *_ = build_engine_from_config(cfg, {"DENV2": "A" * 300})
    assert isinstance(eng_backend, ViennaSaltShiftBackend)
    assert eng_backend.couple_loops_in_salt_dg is False

    salt_model = UnifiedSaltModel()
    backend_default = ViennaSaltShiftBackend(salt_model)
    with dna_session(temp_celsius=65.0):
        raw_v = backend_default.vienna.calc_heterodimer("ATGCATGCTTTATGCATGC", "GCATGCATGCATGCAT", temp_celsius=65.0)
        shifted_def = backend_default.calc_heterodimer(
            "ATGCATGCTTTATGCATGC", "GCATGCATGCATGCAT", temp_celsius=65.0, na_mm=50.0, mg_mm=8.0, dntp_mm=1.4
        )
    dh_stem, ds_stem, n_bp_stem, f_gc_stem, _, _ = estimate_helix_thermo(
        "ATGCATGCTTTATGCATGC&GCATGCATGCATGCAT", raw_v.structure, include_loops=False
    )
    _, ds_stem_corr = salt_model.corrected_thermodynamics(
        dh_stem, ds_stem, f_gc_stem, n_bp_stem, 0.050, 0.0, 0.0, 0.008, 0.0014
    )
    expected_dg_stem_only = raw_v.dg_kcal - (65.0 + 273.15) * (ds_stem_corr - ds_stem) / 1000.0
    assert shifted_def.dg_kcal == pytest.approx(expected_dg_stem_only, abs=1e-12)
    assert shifted_def.dh_kcal != round(dh_stem, 2)


def test_point2_santalucia_hicks_2004_matches_primer3_reference() -> None:
    """Recoupe ligne par ligne les boucles internes et bulges contre Primer3 `thal` / `loops.ds`
    (proprement ignoré via `pytest.importorskip` si `primer3-py` n'est pas installé).

    Cross-validates internal-loop and bulge tables line by line against Primer3 `thal` / `loops.ds`
    (cleanly skipped via `pytest.importorskip` when optional `primer3-py` is not installed).
    """
    import os
    import pytest

    primer3 = pytest.importorskip("primer3")

    from labcraft.thermo.backends.vienna_salt import _bulge_loop_entropy_cal
    from labcraft.thermo.mismatch import _internal_loop_entropy_cal

    loops_ds_path = os.path.join(
        os.path.dirname(primer3.__file__), "src", "libprimer3", "primer3_config", "loops.ds"
    )
    with open(loops_ds_path, encoding="utf-8") as fh:
        for line in fh:
            n_s, ds_int_s, ds_bul_s, _ = line.split()
            n_sz = int(n_s)
            ds_int_ref = float(ds_int_s)
            ds_bul_ref = float(ds_bul_s)
            assert _bulge_loop_entropy_cal(n_sz) == pytest.approx(ds_bul_ref, abs=0.02)
            if n_sz >= 3:
                assert _internal_loop_entropy_cal(n_sz) == pytest.approx(ds_int_ref, abs=0.02)

    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}

    def rc(s: str) -> str:
        return "".join(comp[b] for b in reversed(s))

    for k_b in (1, 2, 4, 6):
        for left, bul_ch, right in (("CGCGCG", "T" * k_b, "GCGCGC"), ("CGCGCA", "C" * k_b, "ACGCGC")):
            top_p = left + right
            top_b = left + bul_ch + right
            bot = rc(top_p)
            struct_p = "(" * len(top_p) + "&" + ")" * len(bot)
            struct_b = "(" * len(left) + "." * k_b + "(" * len(right) + "&" + ")" * len(bot)
            dh0, ds0, *_ = estimate_structure_thermo(f"{top_p}&{bot}", struct_p)
            dhb, dsb, *_ = estimate_structure_thermo(f"{top_b}&{bot}", struct_b)
            p0 = primer3.calc_heterodimer(top_p, bot, mv_conc=1000.0, dv_conc=0.0, dntp_conc=0.0, temp_c=37.0)
            pb = primer3.calc_heterodimer(top_b, bot, mv_conc=1000.0, dv_conc=0.0, dntp_conc=0.0, temp_c=37.0)
            assert (dhb - dh0) == pytest.approx((pb.dh - p0.dh) / 1000.0, abs=0.02)
            assert (dsb - ds0) == pytest.approx(pb.ds - p0.ds, abs=0.05)

    for top_loop, bot_loop in (("AA", "A"), ("AA", "AA"), ("AAA", "AAA")):
        left, right = "CGCGCG", "GCGCGC"
        top_p = left + right
        bot_p = rc(top_p)
        top_i = left + top_loop + right
        bot_i = rc(right) + bot_loop + rc(left)
        struct_p = "(" * len(top_p) + "&" + ")" * len(bot_p)
        struct_i = (
            "(" * len(left) + "." * len(top_loop) + "(" * len(right)
            + "&"
            + ")" * len(right) + "." * len(bot_loop) + ")" * len(left)
        )
        dh0, ds0, *_ = estimate_structure_thermo(f"{top_p}&{bot_p}", struct_p)
        dhi, dsi, *_ = estimate_structure_thermo(f"{top_i}&{bot_i}", struct_i)
        p0 = primer3.calc_heterodimer(top_p, bot_p, mv_conc=1000.0, dv_conc=0.0, dntp_conc=0.0, temp_c=37.0)
        pi = primer3.calc_heterodimer(top_i, bot_i, mv_conc=1000.0, dv_conc=0.0, dntp_conc=0.0, temp_c=37.0)
        assert (dhi - dh0) == pytest.approx((pi.dh - p0.dh) / 1000.0, abs=0.02)
        assert (dsi - ds0) == pytest.approx(pi.ds - p0.ds, abs=0.05)


def test_relative_boltzmann_3p_enrichment_floor_from_dimer_conc_max() -> None:
    """Vérifie que l'enrichissement Boltzmann 3' des dimères non amplifiables suit le seuil relatif
    `conc >= 0.1 * dimer_conc_max` (en mol/L : `0.1 * dimer_conc_max_nm * 1e-9`) dérivé de `--dimer-conc-max`,
    et qu'un dimère sub-nanomolaire (`0.050 nM < 1 nM`) avec `suboptimal_3p_masked = True` est conservé
    dans `evaluate_risks` et le JSON sans gonfler `blocking_dimers_total_nm` (`> 1 nM`).

    Verifies that Boltzmann 3' enrichment of non-amplifiable dimers uses the relative floor
    `conc >= 0.1 * dimer_conc_max` (`0.1 * dimer_conc_max_nm * 1e-9` M) derived from `--dimer-conc-max`,
    and that a sub-nanomolar dimer (`0.050 nM < 1 nM`) with `suboptimal_3p_masked = True` is retained in
    `evaluate_risks` and JSON output without inflating `blocking_dimers_total_nm` (`> 1 nM`).
    """
    import numpy as np
    import pytest

    from labcraft.diagnostics.enzyme import BST_2_0
    from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
    from labcraft.metrics.fractions import PrimerFractions
    from labcraft.metrics.risk import evaluate_risks
    from labcraft.metrics.verdict import PanelVerdict
    from labcraft.pipeline.analysis import evaluate_dimer_complexes
    from labcraft.report.json_export import build_analyze_json_payload
    from labcraft.solver.types import EquilibriumProblem, EquilibriumResult, SolverMethod
    from labcraft.thermo.backends.vienna import ViennaRNABackend

    p_a = PhysicalPrimer.from_simple("SubPanel_F3", "CGCGCGCAAAAAATGCATGC", PrimerRole.F3, nominal_concentration=2e-7)
    p_b = PhysicalPrimer.from_simple("SubPanel_B3", "CGCGCGCTTTTTGCATGCAT", PrimerRole.B3, nominal_concentration=2e-7)
    prob = EquilibriumProblem(
        n_strands=2,
        n_complexes=3,
        stoichiometry=np.array([[1, 0], [0, 1], [1, 1]], dtype=int),
        delta_g=np.array([0.0, 0.0, -6.5], dtype=float),
        total_concentrations=np.array([2e-7, 2e-7], dtype=float),
        temperature_kelvin=63.0 + 273.15,
    )
    # 0,050 nM = 5,0e-11 M : entre 0,023 nM (0,1 * 0,230 nM) et 1,0 nM (ancien filtre codé en dur)
    # 0.050 nM = 5.0e-11 M: between 0.023 nM (0.1 * 0.230 nM) and 1.0 nM (previous hardcoded cutoff)
    res_eq = EquilibriumResult(
        concentrations=np.array([1.9995e-7, 1.9995e-7, 5.0e-11], dtype=float),
        free_concentrations=np.array([1.9995e-7, 1.9995e-7], dtype=float),
        log_free_concentrations=np.log(np.array([1.9995e-7, 1.9995e-7], dtype=float)),
        residuals=np.zeros(2, dtype=float),
        max_residual=0.0,
        n_iterations=1,
        method=SolverMethod.DUAL_NEWTON,
        converged=True,
    )
    c_names = ["SubPanel_F3_free", "SubPanel_B3_free", "SubPanel_F3_SubPanel_B3"]
    backend = ViennaRNABackend()
    pmap = {"SubPanel_F3": "SubPanel", "SubPanel_B3": "SubPanel"}

    # 1. Avec dimer_conc_max_nm = 0.230 (plancher = 0.023 nM <= 0.050 nM) : enrichi et retenu
    amp_flags, concs_out, details_list = evaluate_dimer_complexes(
        prob,
        res_eq,
        c_names,
        [p_a, p_b],
        backend,
        BST_2_0,
        63.0,
        {},
        pmap,
        dimer_conc_max_nm=0.23,
    )
    assert amp_flags == [False, False, False]
    assert details_list[2].get("boltzmann_3p_prob") is not None
    assert details_list[2]["suboptimal_3p_masked"] is True

    risks = evaluate_risks(c_names, concs_out, amp_flags, dimer_details=details_list)
    assert len(risks) == 1
    assert risks[0].suboptimal_3p_masked is True

    payload = build_analyze_json_payload(
        version="0.0.1",
        timestamp="2026-11-19T00:00:00Z",
        inputs={},
        conditions={"temperature_celsius": 63.0, "enzyme": "Bst 2.0", "dimer_dg_threshold": -2.0},
        primer_concentrations_uM={"SubPanel_F3": 0.2, "SubPanel_B3": 0.2},
        fractions={
            "SubPanel_F3": PrimerFractions(0.99, 0.0, 0.0, 0.01, 0.0, 0.0),
            "SubPanel_B3": PrimerFractions(0.99, 0.0, 0.0, 0.01, 0.0, 0.0),
        },
        unfolding_penalties={},
        target_occupations={},
        risks=risks,
        dimer_decomposition={},
        mispriming_risks=[],
        verdict=PanelVerdict(status="OK", issues=[], global_cause="OK"),
    )
    assert len(payload["dimers"]) == 1
    assert payload["dimers"][0]["suboptimal_3p_masked"] is True
    assert payload["totals"]["blocking_dimers_total_nm"] == pytest.approx(0.0, abs=1e-12)

    # 2. Avec dimer_conc_max_nm = 1.000 (plancher = 0.100 nM > 0.050 nM) : non enrichi
    _, _, details_high = evaluate_dimer_complexes(
        prob,
        res_eq,
        c_names,
        [p_a, p_b],
        backend,
        BST_2_0,
        63.0,
        {},
        pmap,
        dimer_conc_max_nm=1.0,
    )
    assert "boltzmann_3p_prob" not in details_high[2]


def test_boltzmann_max_dimers_cap_and_json_skipped_export(tmp_path: Path) -> None:
    """Vérifie le dédoublonnage par couple parent et le plafond `--boltzmann-max-dimers` (défaut 25) :
    1. 100 % des dimères amplifiables sont évalués même si `boltzmann_max_dimers == 2` ;
    2. Les variants dégénérés d'un même couple parent (`P_LF#1 × P_LB` à `0.090 nM` et `P_LF#2 × P_LB`
       à `0.085 nM`) sont dédoublonnés en conservant uniquement le plus concentré (`P_LF#1 × P_LB`),
       de sorte que le plafond `K = 2` s'applique aux couples parents (`(P_LB, P_LF)` puis `(P_LF, P_LF)`) ;
    3. Le JSON exporte `boltzmann_eligible_parent_pairs` et `boltzmann_skipped_parent_pairs` dans `totals`.

    Verifies parent-pair deduplication and `--boltzmann-max-dimers` (default 25):
    1. 100% of amplifiable dimers are evaluated;
    2. Degenerate variants sharing a parent pair keep only the highest-concentration variant before applying `K`;
    3. JSON `totals` exports `boltzmann_eligible_parent_pairs` and `boltzmann_skipped_parent_pairs`.
    """
    import numpy as np

    from labcraft.diagnostics.enzyme import BST_2_0
    from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
    from labcraft.pipeline.analysis import evaluate_dimer_complexes
    from labcraft.solver.types import EquilibriumProblem, EquilibriumResult, SolverMethod
    from labcraft.thermo.backends.vienna import ViennaRNABackend

    p_amp1 = PhysicalPrimer.from_simple("P_F3", "AAAAAAATGCATGC", PrimerRole.F3, nominal_concentration=2e-7)
    p_amp2 = PhysicalPrimer.from_simple("P_B3", "TTTTTTTGCATGCA", PrimerRole.B3, nominal_concentration=2e-7)
    p_lf1 = PhysicalPrimer.from_simple(
        "P_LF#1", "CGCGCGCAAAAAATGCATGC", PrimerRole.LF, nominal_concentration=1e-7, parent_name="P_LF"
    )
    p_lf2 = PhysicalPrimer.from_simple(
        "P_LF#2", "CGCGCGCTAAAAATGCATGC", PrimerRole.LF, nominal_concentration=1e-7, parent_name="P_LF"
    )
    p_lb = PhysicalPrimer.from_simple("P_LB", "CGCGCGCTTTTTGCATGCAT", PrimerRole.LB, nominal_concentration=2e-7)

    primers = [p_amp1, p_amp2, p_lf1, p_lf2, p_lb]
    # 5 monomères + 1 amplifiable + 4 non-amp éligibles répartis sur 3 couples parents :
    # - idx 6 : P_LF#1 × P_LB (0.090 nM, couple parent (P_LB, P_LF) -> représentant retenu)
    # - idx 7 : P_LF#2 × P_LB (0.085 nM, même couple parent (P_LB, P_LF) -> écarté par dédoublonnage !)
    # - idx 8 : P_LF#1_homo   (0.070 nM, couple parent (P_LF, P_LF) -> rang 2 <= K=2, évalué !)
    # - idx 9 : P_LB_homo     (0.050 nM, couple parent (P_LB, P_LB) -> rang 3 > K=2, non évalué)
    stoich = np.array(
        [
            [1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0],
            [0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1],
            [1, 1, 0, 0, 0],  # idx 5 : amplifiable (0.010 nM)
            [0, 0, 1, 0, 1],  # idx 6 : P_LF#1_P_LB (0.090 nM)
            [0, 0, 0, 1, 1],  # idx 7 : P_LF#2_P_LB (0.085 nM)
            [0, 0, 2, 0, 0],  # idx 8 : P_LF#1_homo (0.070 nM)
            [0, 0, 0, 0, 2],  # idx 9 : P_LB_homo   (0.050 nM)
        ],
        dtype=int,
    )
    prob = EquilibriumProblem(
        n_strands=5,
        n_complexes=10,
        stoichiometry=stoich,
        delta_g=np.array([0.0, 0.0, 0.0, 0.0, 0.0, -7.0, -6.5, -6.4, -6.2, -6.0], dtype=float),
        total_concentrations=np.array([2e-7, 2e-7, 1e-7, 1e-7, 2e-7], dtype=float),
        temperature_kelvin=63.0 + 273.15,
    )
    res_eq = EquilibriumResult(
        concentrations=np.array(
            [1.9e-7, 1.9e-7, 0.9e-7, 0.9e-7, 1.9e-7, 1.0e-11, 9.0e-11, 8.5e-11, 7.0e-11, 5.0e-11]
        ),
        free_concentrations=np.array([1.9e-7, 1.9e-7, 0.9e-7, 0.9e-7, 1.9e-7], dtype=float),
        log_free_concentrations=np.log(np.array([1.9e-7, 1.9e-7, 0.9e-7, 0.9e-7, 1.9e-7])),
        residuals=np.zeros(5, dtype=float),
        max_residual=0.0,
        n_iterations=1,
        method=SolverMethod.DUAL_NEWTON,
        converged=True,
    )
    c_names = [
        "P_F3_free",
        "P_B3_free",
        "P_LF#1_free",
        "P_LF#2_free",
        "P_LB_free",
        "P_F3_P_B3",
        "P_LF#1_P_LB",
        "P_LF#2_P_LB",
        "P_LF#1_homo",
        "P_LB_homo",
    ]
    stats: dict[str, object] = {}
    amp_flags, _, details_list = evaluate_dimer_complexes(
        prob,
        res_eq,
        c_names,
        primers,
        ViennaRNABackend(),
        BST_2_0,
        63.0,
        {},
        {"P_F3": "P", "P_B3": "P", "P_LF": "P", "P_LB": "P"},
        dimer_conc_max_nm=0.23,
        boltzmann_max_dimers=2,
        boltzmann_stats_out=stats,
    )
    assert amp_flags[5] is True
    assert "boltzmann_3p_prob" in details_list[5]      # 100 % des amplifiables évalués
    assert "boltzmann_3p_prob" in details_list[6]      # couple #1 (P_LB, P_LF) : variant le plus concentré évalué
    assert "boltzmann_3p_prob" not in details_list[7]  # même couple (P_LB, P_LF) : doublon dégénéré ignoré
    assert "boltzmann_3p_prob" in details_list[8]      # couple #2 (P_LF, P_LF) : évalué grâce au dédoublonnage !
    assert "boltzmann_3p_prob" not in details_list[9]  # couple #3 (P_LB, P_LB) : au-delà de K=2
    assert stats["eligible_non_amp_dimers"] == 4
    assert stats["eligible_parent_pairs"] == 3
    assert stats["evaluated_non_amp_dimers"] == 2
    assert stats["skipped_parent_pairs"] == 1

    # Vérification CLI de bout en bout sur DENV3 (5 couples parents éligibles) avec --boltzmann-max-dimers 2
    gt_dir = Path(__file__).resolve().parent.parent / "fixtures" / "ground_truth"
    out_html = tmp_path / "d3_cap.html"
    runner = CliRunner()
    res_cli = runner.invoke(
        app,
        [
            "analyze",
            "-p", str(gt_dir / "primer_dengue_3_panel3.txt"),
            "-t", str(gt_dir / "DEN3_M93130.fasta"),
            "--temperature", "63.0",
            "--conc-fip-bip", "0.8",
            "--conc-loop", "0.4",
            "--conc-f3-b3", "0.2",
            "--boltzmann-max-dimers", "2",
            "-o", str(out_html),
        ],
    )
    assert res_cli.exit_code == 0, res_cli.stdout
    data = json.loads((tmp_path / "d3_cap.json").read_text(encoding="utf-8"))
    assert data["conditions"]["boltzmann_max_dimers"] == 2
    assert data["totals"]["boltzmann_eligible_parent_pairs"] == 5
    assert data["totals"]["boltzmann_skipped_parent_pairs"] == 3


