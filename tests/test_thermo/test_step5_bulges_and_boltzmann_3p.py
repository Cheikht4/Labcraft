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
