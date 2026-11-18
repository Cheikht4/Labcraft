"""Tests unitaires de non-régression pour l'Étape 1 (correctifs de fiabilité).
Non-regression unit tests for Step 1 (reliability bug fixes).
"""
from __future__ import annotations

import warnings
import pytest
import RNA

from labcraft.cli.config import BufferConfig, ExperimentConfig, PanelConfig, PrimerConfig, PrimerSetConfig, TargetConfig, build_engine_from_config
from labcraft.diagnostics.amplifiable_dimer import is_amplifiable_dimer
from labcraft.diagnostics.enzyme import BST_2_0
from labcraft.lamp.complex_enumeration import enumerate_complexes
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole, _revcomp
from labcraft.optimize.concentrations import optimize_concentrations
from labcraft.optimize.constraints import evaluate_variant_dimers, filter_no_new_artefacts
from labcraft.pipeline.titration import simulate_titration
from labcraft.report.alignment import dotbracket_to_alignment
from labcraft.target.seeding import find_candidate_sites, primer_matches_sequence
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.thermo.salt import UnifiedSaltModel
from labcraft.thermo.vienna import dna_params, dna_session


def test_amplifiable_dimer_strand_b_3prime_veto_direction():
    """Vérifie que `is_amplifiable_dimer` détecte une extrémité 3' extensible sur le brin B (second brin).
    Verifies that `is_amplifiable_dimer` detects an extensible 3' end on strand B (second strand).
    """
    # Miroir exact : A offre la matrice 5' ("TTTT..."), B s'étend en 3' ("...CGCG")
    # Exact mirror: A provides the 5' template ("TTTT..."), B extends at its 3' end ("...CGCG")
    primer_a = "TTTTCGCGCATGCAT"
    primer_b = "ATGCATGCGCG"

    with dna_session(65.0):
        struct, mfe = RNA.cofold(f"{primer_a}&{primer_b}")
        struct_clean = struct.replace("&", "")
        is_amp, dg_3p, ext_strand, blocked = is_amplifiable_dimer(
            primer_a, primer_b, struct_clean, mfe, BST_2_0, temp_celsius=65.0
        )

    assert is_amp is True
    assert ext_strand == "b"
    assert dg_3p <= BST_2_0.dimer_dg_threshold
    assert blocked is False


def test_seeding_antisense_not_masked_by_noisy_sense_hit():
    """Vérifie qu'un faux positif bruité sur '+' ne masque pas un site parfait sur '-' (antisens).
    Verifies that a noisy false positive on '+' does not mask a perfect site on '-' (antisense).
    """
    b3_seq = "CTGGCTTGAGACATCTTCT"  # 19 nt
    # Leurre sur '+' : 2 substitutions en 5' + 1 substitution terminale en 3' (3 erreurs au total)
    # Decoy on '+': 2 substitutions in 5' + 1 terminal 3' substitution (3 errors total)
    decoy_plus = "CAGGATTGAGACATCTTCA"
    # Vrai site parfait sur '-' : complément inverse exact sur le brin '+'
    # True perfect site on '-': exact reverse complement on the '+' strand
    perfect_minus_on_sense = _revcomp(b3_seq)

    genome = "A" * 40 + decoy_plus + "T" * 60 + perfect_minus_on_sense + "A" * 40
    expected_minus_pos = 40 + len(decoy_plus) + 60

    # 1. primer_matches_sequence doit retourner le site parfait sur '-'
    # 1. primer_matches_sequence must return the perfect '-' locus
    hit = primer_matches_sequence(genome, b3_seq, max_errors=2, strict_3prime_len=3, strict_3prime_tolerate=2)
    assert hit is not None
    assert hit == (expected_minus_pos, expected_minus_pos + len(b3_seq), "-")

    # 2. find_candidate_sites doit retenir le site parfait sur '-' (0 mésappariement)
    # 2. find_candidate_sites must retain the perfect '-' site (0 mismatches)
    b3_primer = PhysicalPrimer("P_B3", b3_seq, PrimerRole.B3, b3_seq)
    cands = find_candidate_sites({"S1": genome}, [b3_primer], max_errors=2, strict_3prime_len=3, strict_3prime_tolerate=2)
    assert len(cands) == 1
    assert cands[0]["strand"] == "-"
    assert cands[0]["position"] == expected_minus_pos
    assert cands[0]["n_mismatches"] == 0


def test_seeding_skips_5prime_n_padding():
    """Vérifie qu'un bloc de padding 5' ('NNNN...') ne fait pas échouer la détection du vrai locus en aval.
    Verifies that a 5' assembly gap ('NNNN...') does not abort detection of the downstream true locus.
    """
    f3_seq = "CTCGTGTAGGAATGGGAG"
    genome = ("N" * 250) + ("ACGT" * 10) + f3_seq + ("TGCA" * 10)
    expected_pos = 250 + 40

    f3_primer = PhysicalPrimer("P_F3", f3_seq, PrimerRole.F3, f3_seq)
    cands = find_candidate_sites({"GISAID_padded": genome}, [f3_primer], max_errors=2, strict_3prime_len=3)
    assert len(cands) == 1
    assert cands[0]["position"] == expected_pos
    assert cands[0]["strand"] == "+"
    assert cands[0]["n_mismatches"] == 0

    hit = primer_matches_sequence(genome, f3_seq, max_errors=2, strict_3prime_len=3)
    assert hit == (expected_pos, expected_pos + len(f3_seq), "+")


def test_optimize_concentrations_tracks_degenerate_and_multitarget_initiations(capsys):
    """Vérifie le suivi des occupations d'initiation pour amorces dégénérées et sites `@t_id_site`, sans `print`.
    Verifies initiation occupation tracking for degenerate primers and `@t_id_site` loci, with no stray `print`.
    """
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    t1_seq = "AAAACCCC" + "ATCGTACGATCGATC" + "GGGGTTTT"
    t2_seq = "TTTTGGGG" + "GCTAGCTAGCTAGCT" + "CCCCAAAA"

    # Panel 1 avec FIP dégénérée (#1 matchant T1, #2 mutée) + LB formant un dimère amplifiable
    # Panel 1 with degenerate FIP (#1 matching T1, #2 mutated) + LB forming an amplifiable dimer
    fip_v1 = PhysicalPrimer(
        "FIP_P1#1", "TTTTATCGTACGATCGATC", PrimerRole.FIP, "ATCGTACGATCGATC",
        nominal_concentration=0.8e-6, parent_name="FIP_P1", parent_binding_domain="ATCGTACGATCGATY",
    )
    fip_v2 = PhysicalPrimer(
        "FIP_P1#2", "TTTTATCGTACGATCGATT", PrimerRole.FIP, "ATCGTACGATCGATT",
        nominal_concentration=0.8e-6, parent_name="FIP_P1", parent_binding_domain="ATCGTACGATCGATY",
    )
    lb_p1 = PhysicalPrimer(
        "LB_P1", "AAACCCGGGGATCGTACGAT", PrimerRole.LB, "AAACCCGGGGATCGTACGAT",
        nominal_concentration=0.8e-6, parent_name="LB_P1",
    )
    fip_p2 = PhysicalPrimer(
        "FIP_P2", "AAAAGCTAGCTAGCTAGCT", PrimerRole.FIP, "GCTAGCTAGCTAGCT",
        nominal_concentration=1.6e-6, parent_name="FIP_P2",
    )

    primers = [fip_v1, fip_v2, lb_p1, fip_p2]
    p2t = {
        "FIP_P1": "T1", "FIP_P1#1": "T1", "FIP_P1#2": "T1", "LB_P1": "T1",
        "FIP_P2": "T2",
    }
    prob, strands, complexes, _ = enumerate_complexes(
        primers, {"T1": t1_seq, "T2": t2_seq}, backend, temp_celsius=65.0,
        enzyme=BST_2_0, primer_to_target=p2t,
    )
    assert "FIP_P1@T1_site" in strands
    assert "FIP_P2@T2_site" in strands

    # Même avec max_iter=1, aucun print parasite sur stdout, et les occupations d'initiation sont peuplées
    # Even with max_iter=1, no stray stdout print, and initiation occupations are populated
    res = optimize_concentrations(
        prob_template=prob,
        species_names=strands,
        primers=primers,
        target_dict={"T1": t1_seq, "T2": t2_seq},
        primer_to_panel=p2t,
        original_free_fractions={},
        original_target_occupations={},
        complex_names=complexes,
        temp_celsius=65.0,
        backend=backend,
        enzyme=BST_2_0,
        max_iter=1,
    )
    captured = capsys.readouterr()
    assert "Optimize concentrations stopped" not in captured.out
    if isinstance(res, dict):
        assert "FIP_P1" in res["metrics"]["initiations_before"]
        assert "FIP_P2" in res["metrics"]["initiations_before"]
        assert "LB_P1" not in res["metrics"]["initiations_before"]


def test_config_mg_dntp_complete_chelation_warning_and_conc_map():
    """Vérifie la distinction entre chélation complète (Mg <= dNTP) et ratio faible (< 3.0), ainsi que `conc_map`.
    Verifies distinction between complete chelation (Mg <= dNTP) and low ratio (< 3.0), plus `conc_map`.
    """
    # 1. Mg <= dNTP avec dNTP > 0 -> avertissement de chélation complète
    # 1. Mg <= dNTP with dNTP > 0 -> complete chelation warning
    with pytest.warns(UserWarning, match="Chélation complète du Mg2\\+"):
        PanelConfig(experiment=ExperimentConfig(buffer=BufferConfig(na_mM=50.0, mg_mM=1.4, dntp_mM=1.4)))

    # 2. dNTP < Mg < 3 * dNTP -> avertissement de rapport faible
    # 2. dNTP < Mg < 3 * dNTP -> low Mg/dNTP ratio warning
    with pytest.warns(UserWarning, match="Rapport Mg/dNTP faible"):
        PanelConfig(experiment=ExperimentConfig(buffer=BufferConfig(na_mM=50.0, mg_mM=2.8, dntp_mM=1.4)))

    # 3. Tampon purement monovalent (Mg=0, dNTP=0) -> aucun avertissement de chélation
    # 3. Purely monovalent buffer (Mg=0, dNTP=0) -> no chelation warning
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        cfg = PanelConfig(
            experiment=ExperimentConfig(buffer=BufferConfig(na_mM=50.0, mg_mM=0.0, dntp_mM=0.0)),
            targets=[TargetConfig(id="T1", sequence_file="dummy.fa")],
            primer_sets=[
                PrimerSetConfig(
                    target="T1",
                    panel_name="P1",
                    primers={
                        "F3": PrimerConfig(seq="ATGCATGCATGCATGC", conc_uM=0.35),
                        "FIP": PrimerConfig(seq="GCATGCATGCATGCATATGCATGCATGCATGC", conc_uM=1.2),
                    },
                )
            ],
        )
    assert not any("Mg" in str(w.message) for w in rec)

    _, _, _, _, _, _, _, profiles = build_engine_from_config(cfg, {"T1": "ATGC" * 20})
    assert profiles["T1"].f3_b3 == pytest.approx(0.35e-6)
    assert profiles["T1"].fip_bip == pytest.approx(1.2e-6)


def test_titration_constraints_and_alignment_fixes():
    """Vérifie `simulate_titration` (_p.name), `evaluate_variant_dimers` (4-tuple + 2 ordres) et `dotbracket_to_alignment`.
    Verifies `simulate_titration` (_p.name), `evaluate_variant_dimers` (4-tuple + 2 orders), and `dotbracket_to_alignment`.
    """
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())

    # 1. Titration : P1 et P2 fortement complémentaires (< 20 % libres), P3 inerte en dernière position
    # 1. Titration: P1 and P2 strongly complementary (< 20% free), inert P3 placed last
    p1 = PhysicalPrimer("P1", "GCGCGCGCGCGCGCGCGC", PrimerRole.FIP, "GCGCGCGCGCGCGCGCGC")
    p2 = PhysicalPrimer("P2", "GCGCGCGCGCGCGCGCGC", PrimerRole.BIP, "GCGCGCGCGCGCGCGCGC")
    p3 = PhysicalPrimer("P3_inert", "TTTTTTTTTTTTTTTTTT", PrimerRole.F3, "TTTTTTTTTTTTTTTTTT")
    tit = simulate_titration(
        [p1, p2, p3],
        backend,
        {"P1": 2.0e-6, "P2": 2.0e-6, "P3_inert": 0.2e-6},
        dilutions=[1.0, 0.5],
        temp_celsius=65.0,
    )
    assert "P1" in tit.dominant_complexes
    assert "P2" in tit.dominant_complexes
    assert "P3_inert" not in tit.dominant_complexes

    # 2. evaluate_variant_dimers & filter_no_new_artefacts (dépaquetage 4 valeurs + ordre inverse)
    # 2. evaluate_variant_dimers & filter_no_new_artefacts (4-value unpack + reverse order)
    panel_dicts = [
        {"name": "VAR", "sequence": "TTTTCGCGCATGCAT"},
        {"name": "PARTNER", "sequence": "ATGCATGCGCG"},
    ]
    with dna_params(65.0):
        amp_list = evaluate_variant_dimers("TTTTCGCGCATGCAT", panel_dicts, 0, BST_2_0, 65.0)
        is_ok, _ = filter_no_new_artefacts("TTTTCGCGCATGCAT", panel_dicts, 0, BST_2_0, 65.0)
    assert len(amp_list) == 1
    assert amp_list[0][0] == 1
    assert is_ok is False

    # 3. dotbracket_to_alignment avec et sans '&'
    # 3. dotbracket_to_alignment with and without '&'
    a_with = dotbracket_to_alignment("ATGC", "GCAT", "((((&))))")
    a_without = dotbracket_to_alignment("ATGC", "GCAT", "(((())))")
    assert a_with == a_without
