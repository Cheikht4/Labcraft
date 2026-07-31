"""Tests for primal solver.

Tests pour le solveur primal.
"""
from __future__ import annotations

import numpy as np
from labcraft.solver.types import EquilibriumProblem, SolverMethod
from labcraft.solver.primal import solve_primal
from labcraft.solver.dual import solve_dual

def test_primal_two_strand():
    """Cas à 2 brins pour primal / 2-strand case."""
    a0, b0 = 1e-6, 1e-6
    dg = -10.0
    t = 338.15
    
    prob = EquilibriumProblem(
        n_strands=2,
        n_complexes=3,
        stoichiometry=np.array([[1, 0], [0, 1], [1, 1]], dtype=np.float64),
        delta_g=np.array([0.0, 0.0, dg]),
        total_concentrations=np.array([a0, b0]),
        temperature_kelvin=t
    )
    
    res_primal = solve_primal(prob, convergence_threshold=1e-10)
    res_dual = solve_dual(prob)
    
    np.testing.assert_allclose(res_primal.concentrations, res_dual.concentrations, rtol=1e-4, atol=1e-15)
    assert res_primal.method == SolverMethod.TRUST_CONSTR
    assert res_primal.converged

def test_primal_three_strand():
    """Cas à 3 brins / 3-strand case."""
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
    
    res_primal = solve_primal(prob)
    res_dual = solve_dual(prob)
    
    np.testing.assert_allclose(res_primal.concentrations, res_dual.concentrations, rtol=1e-4, atol=1e-15)
    assert res_primal.converged
