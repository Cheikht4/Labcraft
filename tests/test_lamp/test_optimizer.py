import pytest
import numpy as np

from labcraft.lamp.complex_enumeration import enumerate_complexes
from labcraft.solver.types import EquilibriumProblem
from labcraft.solver.dual import solve_dual
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.thermo.backends.native import NativeBackend
from labcraft.thermo.salt import UnifiedSaltModel, SaltCorrectedBackend

def test_guardrail_lower_concentration_lowers_initiation():
    """
    Test de garde-fou pour s'assurer que baisser la concentration
    d'une amorce d'initiation (ex: FIP) diminue son occupation sur la cible.
    """
    backend = SaltCorrectedBackend(NativeBackend(), UnifiedSaltModel())
    
    target_seq = "ATGCGTACGTGCAACTGATCGATCGTACGATCG"
    
    # FIP parfait avec binding de 15 bases sur la cible (GATCGATCGTACGAT)
    # L'inverse complémentaire de ce binding est ATCGTACGATCGATC
    binding = "ATCGTACGATCGATC"
    fip = PhysicalPrimer("FIP_1", sequence="TTTT" + binding, role=PrimerRole.FIP, binding_domain=binding, nominal_concentration=1.6e-6)
    
    prob_high, species_high, complexes_high, _ = enumerate_complexes(
        primers=[fip],
        target_seq=target_seq,
        backend=backend,
        temp_celsius=65.0
    )
    res_high = solve_dual(prob_high)
    
    # Trouver l'index de FIP_1_target dans complexes
    target_complex_name = 'FIP_1_on_FIP_1_site'
    target_complex_idx = complexes_high.index(target_complex_name)
    occ_high = res_high.concentrations[target_complex_idx]
    
    # On baisse la concentration de FIP par 4
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
