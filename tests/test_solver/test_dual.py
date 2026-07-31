"""Tests for dual solver.

Tests pour le solveur dual.
"""
from __future__ import annotations

import numpy as np
import pytest
from labcraft.solver.types import EquilibriumProblem, SolverMethod, ConvergenceError
from labcraft.solver.dual import solve_dual
from labcraft.solver.analytical import analytical_two_strand

def test_dual_two_strand():
    """Cas à 2 brins / 2-strand case."""
    a0, b0 = 1e-6, 1e-6
    dg = -10.0
    t = 338.15
    
    ref = analytical_two_strand(a0, b0, dg, t)
    
    prob = EquilibriumProblem(
        n_strands=2,
        n_complexes=3,
        stoichiometry=np.array([[1, 0], [0, 1], [1, 1]], dtype=np.float64),
        delta_g=np.array([0.0, 0.0, dg]),
        total_concentrations=np.array([a0, b0]),
        temperature_kelvin=t
    )
    
    res = solve_dual(prob)
    
    np.testing.assert_allclose(res.concentrations, ref.concentrations, rtol=1e-12, atol=1e-15)
    assert res.method == SolverMethod.DUAL_NEWTON
    assert res.converged
    assert res.max_residual < 1e-10

def test_dual_three_strand():
    """Cas à 3 brins / 3-strand case."""
    # A+B ⇌ AB, B+C ⇌ BC
    # Complexes: A, B, C, AB, BC
    # Stoichiometry:
    # A:  1 0 0
    # B:  0 1 0
    # C:  0 0 1
    # AB: 1 1 0
    # BC: 0 1 1
    prob = EquilibriumProblem(
        n_strands=3,
        n_complexes=5,
        stoichiometry=np.array([
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 1, 0],
            [0, 1, 1]
        ], dtype=np.float64),
        delta_g=np.array([0.0, 0.0, 0.0, -10.0, -12.0]),
        total_concentrations=np.array([1e-6, 1e-6, 1e-6]),
        temperature_kelvin=310.15
    )
    
    res = solve_dual(prob)
    assert res.converged
    assert res.max_residual < 1e-10

def test_extreme_case():
    """Cas extrême / Extreme case."""
    prob = EquilibriumProblem(
        n_strands=2,
        n_complexes=3,
        stoichiometry=np.array([[1, 0], [0, 1], [1, 1]], dtype=np.float64),
        delta_g=np.array([0.0, 0.0, -25.0]),
        total_concentrations=np.array([1e-3, 1e-12]),
        temperature_kelvin=310.15
    )
    res = solve_dual(prob)
    assert res.converged
    assert res.max_residual < 1e-10

def test_pathological_zero_dg():
    """Cas pathologique avec ΔG° = 0 / Pathological case."""
    a0, b0 = 1e-6, 1e-6
    prob = EquilibriumProblem(
        n_strands=2,
        n_complexes=3,
        stoichiometry=np.array([[1, 0], [0, 1], [1, 1]], dtype=np.float64),
        delta_g=np.array([0.0, 0.0, 0.0]),
        total_concentrations=np.array([a0, b0]),
        temperature_kelvin=310.15
    )
    res = solve_dual(prob)
    assert res.converged
    
def test_convergence_error():
    """Vérifie ConvergenceError / Check ConvergenceError."""
    prob = EquilibriumProblem(
        n_strands=2,
        n_complexes=3,
        stoichiometry=np.array([[1, 0], [0, 1], [1, 1]], dtype=np.float64),
        delta_g=np.array([0.0, 0.0, -10.0]),
        total_concentrations=np.array([1e-6, 1e-6]),
        temperature_kelvin=310.15
    )
    
    with pytest.raises(ConvergenceError):
        # Unrealistic threshold and 1 max iter
        solve_dual(prob, convergence_threshold=1e-30, max_iterations=1)
