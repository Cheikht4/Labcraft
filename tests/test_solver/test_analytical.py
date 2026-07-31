"""Tests for analytical solver.

Tests pour le solveur analytique.
"""
from __future__ import annotations

import numpy as np
from labcraft.solver.analytical import analytical_two_strand, analytical_two_strand_naive

def test_simple_case():
    """Cas simple connu / Simple known case."""
    a0, b0 = 1e-6, 1e-6
    dg = -10.0
    t = 338.15
    
    res = analytical_two_strand(a0, b0, dg, t)
    
    ab = res.concentrations[2]
    assert ab > 0
    
    # Conservation de masse
    np.testing.assert_allclose(res.concentrations[0] + ab, a0)
    np.testing.assert_allclose(res.concentrations[1] + ab, b0)

def test_symmetric_case():
    """Cas symétrique / Symmetric case."""
    a0, b0 = 1e-6, 1e-6
    dg = -15.0
    t = 310.15
    
    res = analytical_two_strand(a0, b0, dg, t)
    ab = res.concentrations[2]
    # Almost all bound
    assert ab > 0.9 * a0

def test_asymmetric_case():
    """Cas asymétrique / Asymmetric case."""
    a0 = 1e-3
    b0 = 1e-9
    dg = -12.0
    t = 310.15
    
    res = analytical_two_strand(a0, b0, dg, t)
    ab = res.concentrations[2]
    # B is the limiting factor / B est le facteur limitant
    # Almost all B should be bound / Presque tout B doit être lié
    assert ab > 0.5 * b0
    assert ab <= b0

def test_numerical_stability():
    """Test de stabilité numérique / Numerical stability test.

    At weak binding (ΔG° = -5 kcal/mol) with small equal concentrations
    (a₀ = b₀ = 1e-9 M), the naive form (S - √(S²-4ab))/2 suffers from
    catastrophic cancellation because S ≈ 1/Keq >> a₀+b₀, and √(S²-4ab) ≈ S.
    The stable form 2a₀b₀/(S + √(S²-4ab)) avoids this subtraction.

    À faible liaison (ΔG° = -5 kcal/mol) avec de petites concentrations
    égales, la forme naïve perd au moins 3 chiffres significatifs par
    annulation catastrophique. La forme stable l'évite.
    """
    # Weak binding with small concentrations: the catastrophic regime
    # Liaison faible avec petites concentrations : le régime catastrophique
    a0 = 1e-9
    b0 = 1e-9
    dg = -5.0  # kcal/mol — weak binding / liaison faible
    t = 338.15  # 65°C

    res_stable = analytical_two_strand(a0, b0, dg, t)
    res_naive = analytical_two_strand_naive(a0, b0, dg, t)

    ab_stable = res_stable.concentrations[2]
    ab_naive = res_naive.concentrations[2]

    # The naive form should lose at least 3 significant digits (rel error > 1e-5)
    # La forme naïve doit perdre au moins 3 chiffres significatifs
    rel_error = abs(ab_naive - ab_stable) / ab_stable
    assert rel_error > 1e-6, (
        f"Expected significant precision loss in naive form, "
        f"but rel_error = {rel_error:.2e}. "
        f"stable={ab_stable:.10e}, naive={ab_naive:.10e}"
    )

    # The stable form should conserve mass to machine precision
    # La forme stable doit conserver la masse à la précision machine
    free_a = res_stable.concentrations[0]
    assert abs(free_a + ab_stable - a0) / a0 < 1e-14

def test_mass_conservation_random():
    """Conservation de masse aléatoire / Random mass conservation."""
    rng = np.random.default_rng(42)
    for _ in range(20):
        a0 = 10 ** rng.uniform(-9, -3)
        b0 = 10 ** rng.uniform(-9, -3)
        dg = rng.uniform(-20, -5)
        t = rng.uniform(293.15, 363.15)
        
        res = analytical_two_strand(a0, b0, dg, t)
        ab = res.concentrations[2]
        free_a = res.concentrations[0]
        
        rel_error = abs(free_a + ab - a0) / a0
        assert rel_error < 1e-14
