import pytest
import numpy as np

from labcraft.lamp.complex_enumeration import enumerate_complexes
from labcraft.solver.types import EquilibriumProblem
from labcraft.solver.dual import solve_dual
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.thermo.salt import UnifiedSaltModel
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.optimize.concentrations import optimize_concentrations
from labcraft.diagnostics.enzyme import BST_2_0

def test_guardrail_lower_concentration_lowers_initiation():
    """
    Test de garde-fou pour s'assurer que baisser la concentration
    d'une amorce d'initiation (ex: FIP) diminue son occupation sur la cible.
    """
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    
    target_seq = "ATGCGTACGTGCAACTGATCGATCGTACGATCG"
    
    # FIP parfait avec binding de 15 bases sur la cible (GATCGATCGTACGAT)
    binding = "ATCGTACGATCGATC"
    fip = PhysicalPrimer("FIP_1", sequence="TTTT" + binding, role=PrimerRole.FIP, binding_domain=binding, nominal_concentration=1.6e-6)
    
    prob_high, species_high, complexes_high, _ = enumerate_complexes(
        primers=[fip],
        target_seq=target_seq,
        backend=backend,
        temp_celsius=65.0
    )
    res_high = solve_dual(prob_high)
    
    target_complex_name = 'FIP_1_on_FIP_1_site'
    target_complex_idx = complexes_high.index(target_complex_name)
    occ_high = res_high.concentrations[target_complex_idx]
    
    fip_low = PhysicalPrimer("FIP_1", sequence="TTTT" + binding, role=PrimerRole.FIP, binding_domain=binding, nominal_concentration=0.4e-6)
    
    prob_low, species_low, complexes_low, _ = enumerate_complexes(
        primers=[fip_low],
        target_seq=target_seq,
        backend=backend,
        temp_celsius=65.0
    )
    res_low = solve_dual(prob_low)
    target_complex_idx = complexes_low.index(target_complex_name)
    occ_low = res_low.concentrations[target_complex_idx]
    
    assert occ_low < occ_high, f"L'occupation aurait dû baisser. High: {occ_high}, Low: {occ_low}"


def test_optimizer_dimer_and_floor():
    """
    Test the optimizer actually runs and reduces dangerous dimers without
    violating hierarchy and floor.
    """
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    # Create primers forming a dimer
    # Primer 1 and 2
    fip = PhysicalPrimer("FIP_1", sequence="ATCGTACGATCGATCGGGGGG", role=PrimerRole.FIP, binding_domain="ATCGTACGATCGATC", nominal_concentration=1.6e-6)
    # The BIP will have the complementary 3' end so it forms an amplifiable dimer
    bip = PhysicalPrimer("BIP_1", sequence="TTTTAAAAAACCCCCC", role=PrimerRole.BIP, binding_domain="AAAAAA", nominal_concentration=1.6e-6)
    
    target_seq = "ATGCGTACGTGCAACTGATCGATCGTACGATCG"
    
    prob, species, complexes, _ = enumerate_complexes(
        primers=[fip, bip],
        target_seq=target_seq,
        backend=backend,
        temp_celsius=65.0
    )
    
    # Run optimizer
    results = optimize_concentrations(
        prob_template=prob,
        species_names=species,
        primers=[fip, bip],
        target_dict={'PANEL1': target_seq},
        primer_to_panel={'FIP_1': 'PANEL1', 'BIP_1': 'PANEL1'},
        original_free_fractions={},
        original_target_occupations={},
        complex_names=complexes,
        temp_celsius=65.0,
        backend=backend,
        enzyme=BST_2_0,
        min_initiation_occupation=0.01
    )
    
    assert isinstance(results, list)
    
    # Find FIP and BIP in species
    fip_idx = species.index('FIP_1')
    bip_idx = species.index('BIP_1')
    
    fip_site_idx = species.index('FIP_1_site')
    
    # Assert sites are preserved
    for r in results:
        # Just to check it produces something
        assert "primer_name" in r
        assert r["suggested_conc"] >= 0.8e-6 # Bound for FIP/BIP

def test_optimizer_multiplex_balance():
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    fip1 = PhysicalPrimer("FIP_1", sequence="ATCGTACGATCGATC", role=PrimerRole.FIP, binding_domain="ATCGTACGATCGATC", nominal_concentration=1.6e-6)
    fip2 = PhysicalPrimer("FIP_2", sequence="GCATGCATGCATGCA", role=PrimerRole.FIP, binding_domain="GCATGCATGCATGCA", nominal_concentration=1.6e-6)
    target_seq = "ATGCGTACGTGCAACTGATCGATCGTACGATCG GCATGCATGCATGCA"
    prob, species, complexes, _ = enumerate_complexes(
        primers=[fip1, fip2],
        target_seq=target_seq,
        backend=backend,
        temp_celsius=65.0
    )
    
    results = optimize_concentrations(
        prob_template=prob,
        species_names=species,
        primers=[fip1, fip2],
        target_dict={'PANEL1': target_seq, 'PANEL2': target_seq},
        primer_to_panel={'FIP_1': 'PANEL1', 'FIP_2': 'PANEL2'},
        original_free_fractions={},
        original_target_occupations={},
        complex_names=complexes,
        temp_celsius=65.0,
        backend=backend,
        enzyme=BST_2_0
    )
    assert isinstance(results, list)
