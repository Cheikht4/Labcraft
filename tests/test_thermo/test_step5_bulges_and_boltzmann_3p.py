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
