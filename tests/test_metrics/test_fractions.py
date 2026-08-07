import pytest
import numpy as np
from labcraft.metrics.fractions import compute_fractions

def test_fractions_sum_to_one():
    # Simulation d'un cas jouet
    primer_names = ["P1", "P2", "T_site"]
    complex_names = [
        "P1_free", "P2_free", "T_site_free", 
        "P1_homo", "P1_P2", "P1_on_T"
    ]
    
    # 6 complexes, 3 espèces
    stoich = np.array([
        [1, 0, 0], # P1_free
        [0, 1, 0], # P2_free
        [0, 0, 1], # T_site_free
        [2, 0, 0], # P1_homo (attention au 2)
        [1, 1, 0], # P1_P2
        [1, 0, 1], # P1_on_T
    ])
    
    # Dg arbitraires
    dg = np.array([0.0, 0.0, 0.0, -10.0, -5.0, -15.0])
    
    # R * (273.15 + 65) = 0.672
    RT = 0.001987 * (273.15 + 65.0)
    
    # On fixe les free_concentrations pour voir si la conservation est bonne
    free_conc = np.array([1e-6, 1e-6, 1e-12])
    
    fractions = compute_fractions(
        primer_names, complex_names, stoich, free_conc, dg, 65.0
    )
    
    # T_site ne doit pas être dans les fractions
    assert "T_site" not in fractions
    
    # La somme des fractions de P1 et P2 doit faire exactement 1.0
    assert abs(fractions["P1"].sum - 1.0) < 1e-6
    assert abs(fractions["P2"].sum - 1.0) < 1e-6

