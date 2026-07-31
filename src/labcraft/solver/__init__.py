"""Equilibrium solver package / Package du solveur d'équilibre."""
from __future__ import annotations

from .types import (
    EquilibriumProblem,
    EquilibriumResult,
    SolverMethod,
    ConvergenceError
)
from .solver import solve_equilibrium
from .analytical import analytical_two_strand, analytical_two_strand_naive

__all__ = [
    "EquilibriumProblem",
    "EquilibriumResult",
    "SolverMethod",
    "ConvergenceError",
    "solve_equilibrium",
    "analytical_two_strand",
    "analytical_two_strand_naive"
]
