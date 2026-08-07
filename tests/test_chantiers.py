import pytest
from labcraft.diagnostics.enzyme import get_enzyme

def test_enzyme_registry():
    # bst3.0 must have a less negative threshold than bst2.0
    bst2 = get_enzyme("bst2.0")
    bst3 = get_enzyme("bst3.0")
    assert bst3.dimer_dg_threshold > bst2.dimer_dg_threshold
    assert bst3.dimer_dg_threshold == -1.5

def test_enzyme_override():
    spec = {"name": "bst2.0", "dimer_dg_threshold": -1.0, "three_prime_window": 5}
    enz = get_enzyme(spec)
    assert enz.name == "Bst 2.0"
    assert enz.dimer_dg_threshold == -1.0
    assert enz.three_prime_window == 5

# Tests for Chantier 2
from labcraft.thermo.mismatch import nn_duplex_energy, three_prime_extensible

def test_mismatch_perfect():
    # Non-régression parfait : nn_duplex_energy("CGTGACGTCACG", "GCACTGCAGTGC") donne ΔG = -8.449 kcal/mol
    dh, ds, dg = nn_duplex_energy("CGTGACGTCACG", "GCACTGCAGTGC")
    assert abs(dg - -8.45) < 0.05
    assert abs(dh - -98.6) < 0.1
    assert abs(ds - -266.6) < 0.1

def test_mismatch_internal():
    # Mésappariement interne unique, déstabilisation attendue
    # pos 3 (G face à A), brin bas GCAATGCAGTGC : ΔG ≈ -5.97
    dh, ds, dg = nn_duplex_energy("CGTGACGTCACG", "GCAATGCAGTGC")
    assert abs(dg - -5.97) < 0.05
    
    # pos 5 (C face à A), brin bas GCACTACAGTGC : ΔG ≈ -5.21
    dh, ds, dg = nn_duplex_energy("CGTGACGTCACG", "GCACTACAGTGC")
    assert abs(dg - -5.21) < 0.05
    
    # pos 7 (T face à C), brin bas GCACTGCCGTGC : ΔG ≈ -5.82
    dh, ds, dg = nn_duplex_energy("CGTGACGTCACG", "GCACTGCCGTGC")
    assert abs(dg - -5.82) < 0.05

def test_three_prime_extensible():
    enz = get_enzyme("bst")
    enz.three_prime_window = 3
    
    # OK
    ext, bad = three_prime_extensible("ATGC", "TACG", enz)
    assert ext is True
    
    # Veto terminal (pos 1)
    ext, bad = three_prime_extensible("ATGC", "TACA", enz)
    assert ext is False
    assert bad == 1
    
    # Veto pénultième (pos 2)
    ext, bad = three_prime_extensible("ATGC", "TATG", enz)
    assert ext is False
    assert bad == 2
    
    # Veto antépénultième (pos 3)
    ext, bad = three_prime_extensible("ATGC", "TTCG", enz)
    assert ext is False
    assert bad == 3
    
    # Loin du 3'
    ext, bad = three_prime_extensible("ATGCATGC", "AACGTACG", enz)
    assert ext is True
