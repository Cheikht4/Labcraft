"""Tests unitaires et d'intégration de l'Étape 6 : épingles intramoléculaires et connexité intermoléculaire.

Unit and integration tests for Step 6: monomolecular `{p.name}_hairpin` enumeration,
mass-balance conservation (`PrimerFractions.sum == 1.0`), JSON `"hairpin"` export, and
intermolecular connectivity enforcement on `RNA.cofold` dimers via `RNA.duplexfold`.
"""
from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from labcraft.cli.main import app
from labcraft.lamp.complex_enumeration import enumerate_complexes
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.lamp.stoichiometry import ConcentrationProfile
from labcraft.metrics.fractions import compute_fractions
from labcraft.solver.dual import solve_dual
from labcraft.thermo.backends.native import NativeBackend
from labcraft.thermo.backends.vienna import ViennaRNABackend, _has_intermolecular_pair
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.thermo.salt import UnifiedSaltModel


def test_hairpin_enumeration_fractions_sum_and_native_fallback() -> None:
    """Vérifie l'énumération `{p}_hairpin`, `PrimerFractions.sum == 1.0` et le repli propre de `NativeBackend`.
    Verifies `{p}_hairpin` enumeration, `PrimerFractions.sum == 1.0`, and graceful `NativeBackend` fallback.
    """
    p_hp = PhysicalPrimer(
        name="P_hp",
        sequence="CGCGCGCAAAAGCGCGCG",
        role=PrimerRole.F3,
        binding_domain="CGCGCGCAAAAGCGCGCG",
        nominal_concentration=0.2e-6,
    )
    p_poly_a = PhysicalPrimer(
        name="P_polyA",
        sequence="AAAAAAAAAAAAAAAAAA",
        role=PrimerRole.B3,
        binding_domain="AAAAAAAAAAAAAAAAAA",
        nominal_concentration=0.2e-6,
    )
    profile = ConcentrationProfile(target=1e-12, fip_bip=0.8e-6, f3_b3=0.2e-6, lf_lb=0.4e-6)
    salt_backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    buf = {"na_mM": 50.0, "mg_mM": 8.0, "dntp_mM": 1.4}

    prob, strands, complexes, _ = enumerate_complexes(
        [p_hp, p_poly_a],
        None,
        salt_backend,
        profile=profile,
        temp_celsius=63.0,
        mon_molar=0.192,
        buffer=buf,
    )
    assert "P_hp_hairpin" in complexes
    assert "P_polyA_hairpin" not in complexes
    # Aucun faux hétérodimère entre C/G/A pur et poly-A pur malgré l'épingle de P_hp
    # No spurious heterodimer between pure C/G/A and pure poly-A despite P_hp's strong hairpin
    assert "P_hp_P_polyA" not in complexes

    res = solve_dual(prob)
    fracs = compute_fractions(
        strands, complexes, prob.stoichiometry, res.free_concentrations, prob.delta_g, 63.0
    )
    assert fracs["P_hp"].hairpin > 0.50
    assert abs(fracs["P_hp"].sum - 1.0) < 1e-9
    assert abs(fracs["P_polyA"].free - 1.0) < 1e-9
    assert abs(fracs["P_polyA"].sum - 1.0) < 1e-9

    # NativeBackend lève NotImplementedError sur calc_hairpin : doit être géré sans erreur
    # NativeBackend raises NotImplementedError on calc_hairpin: must be handled gracefully
    prob_nat, _st_n, comp_nat, _ = enumerate_complexes(
        [p_hp], None, NativeBackend(), profile=profile, temp_celsius=63.0
    )
    assert "P_hp_hairpin" not in comp_nat
    assert prob_nat.n_complexes >= 1


def test_intermolecular_connectivity_enforced_via_duplexfold() -> None:
    """Vérifie que `calc_heterodimer` rejette les co-repliements disjoints et retourne un duplex intermoléculaire.
    Verifies `calc_heterodimer` rejects disconnected cofolds and returns a true intermolecular duplex.
    """
    vienna = ViennaRNABackend()
    # s1 possède une épingle très stable en 5' ET une queue 3' complémentaire à s2
    # s1 has a strong 5' hairpin AND a 3' tail complementary to s2
    s1 = "CGCGCGCAAAAGCGCGCGATGCATGC"
    s2 = "GCATGCATAAAAAAAAAA"
    res = vienna.calc_heterodimer(s1, s2, temp_celsius=63.0)
    assert len(res.structure) == len(s1) + len(s2)
    if res.dg_kcal < 0.0:
        assert _has_intermolecular_pair(res.structure, len(s1))

    # Cas sans aucune complémentarité inter-brins : dg_kcal == 0.0 au lieu de l'épingle seule (-6,66 kcal/mol)
    # Case with zero cross-strand complementarity: dg_kcal == 0.0 instead of single-strand hairpin (-6.66 kcal/mol)
    res_disjoint = vienna.calc_heterodimer("CGCGCGCAAAAGCGCGCG", "AAAAAAAAAAAAAAAAAA", temp_celsius=63.0)
    assert res_disjoint.dg_kcal >= 0.0


def test_cli_analyze_json_exports_hairpin_fraction_without_fake_dimer(tmp_path: Path) -> None:
    """Vérifie que `labcraft analyze` exporte `fractions[p].hairpin > 0` dans le JSON sans faux dimère bloquant.
    Verifies `labcraft analyze` exports `fractions[p].hairpin > 0` in JSON without fake blocking dimers.
    """
    primers_file = tmp_path / "hp_primers.txt"
    target_file = tmp_path / "hp_target.fasta"
    out_html = tmp_path / "hp_report.html"

    target_seq = (
        "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
        "CGCGCGCAAAAGCGCGCG"
        "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
        "TTTTTTTTTTTTTTTTTT"
        "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    )
    target_file.write_text(f">T1\n{target_seq}\n", encoding="utf-8")
    primers_file.write_text(
        ">HPTest_F3\nCGCGCGCAAAAGCGCGCG\n"
        ">HPTest_B3\nAAAAAAAAAAAAAAAAAA\n",
        encoding="utf-8",
    )

    res = CliRunner().invoke(
        app,
        [
            "analyze",
            "-p", str(primers_file),
            "-t", str(target_file),
            "--temperature", "63.0",
            "-o", str(out_html),
        ],
    )
    assert res.exit_code == 0, res.stdout
    js = json.loads((tmp_path / "hp_report.json").read_text(encoding="utf-8"))
    assert "hairpin" in js["fractions"]["F3_HPTest"]
    assert js["fractions"]["F3_HPTest"]["hairpin"] > 0.50
    assert not any(d["name"].endswith("_hairpin") for d in js["dimers"])
    assert not any("F3_HPTest" in d["name"] and "B3_HPTest" in d["name"] for d in js["dimers"])
