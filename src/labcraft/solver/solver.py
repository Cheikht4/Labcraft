"""Main solver facade for equilibrium problems.

Façade principale pour les problèmes d'équilibre.
"""
from __future__ import annotations

import logging
from .types import EquilibriumProblem, EquilibriumResult, SolverMethod, ConvergenceError
from .dual import solve_dual
from .primal import solve_primal
from .extended import solve_extended

logger = logging.getLogger(__name__)


def solve_equilibrium(
    problem: EquilibriumProblem,
    *,
    method: SolverMethod | None = None,
    convergence_threshold: float = 1e-10,
) -> EquilibriumResult:
    """Solve the multi-state equilibrium problem.
    
    Résout le problème d'équilibre multi-états. Tente plusieurs méthodes en cascade.
    
    Args:
        problem: EquilibriumProblem to solve.
        method: Specific method to use, or None for automatic cascade.
        convergence_threshold: Threshold for relative mass conservation residual.
        
    Returns:
        EquilibriumResult object.
        
    Raises:
        ConvergenceError: If all attempted methods fail.
    """
    problem.validate()
    
    errors = []
    
    if method is None or method == SolverMethod.DUAL_NEWTON:
        try:
            logger.info("Attempting DUAL_NEWTON solver / Tentative du solveur DUAL_NEWTON")
            return solve_dual(problem, convergence_threshold=convergence_threshold)
        except Exception as e:
            logger.warning(f"DUAL_NEWTON failed: {type(e).__name__} - {e}")
            errors.append(f"DUAL_NEWTON: {type(e).__name__} - {e}")
            if method == SolverMethod.DUAL_NEWTON:
                raise
                
    if method is None or method == SolverMethod.TRUST_CONSTR:
        try:
            logger.info("Attempting TRUST_CONSTR solver / Tentative du solveur TRUST_CONSTR")
            return solve_primal(problem, convergence_threshold=convergence_threshold)
        except Exception as e:
            logger.warning(f"TRUST_CONSTR failed: {type(e).__name__} - {e}")
            errors.append(f"TRUST_CONSTR: {type(e).__name__} - {e}")
            if method == SolverMethod.TRUST_CONSTR:
                raise
                
    if method is None or method == SolverMethod.EXTENDED_PRECISION:
        try:
            logger.info("Attempting EXTENDED_PRECISION solver / Tentative de EXTENDED_PRECISION")
            return solve_extended(problem, convergence_threshold=convergence_threshold)
        except ImportError as e:
            logger.warning(f"EXTENDED_PRECISION unavailable: {e}")
            errors.append(f"EXTENDED_PRECISION: {e}")
            if method == SolverMethod.EXTENDED_PRECISION:
                raise
        except Exception as e:
            logger.warning(f"EXTENDED_PRECISION failed: {type(e).__name__} - {e}")
            errors.append(f"EXTENDED_PRECISION: {type(e).__name__} - {e}")
            if method == SolverMethod.EXTENDED_PRECISION:
                raise
                
    # If we get here, all attempted methods failed
    # Si nous arrivons ici, toutes les méthodes tentées ont échoué
    raise ConvergenceError(
        "All solvers failed to converge. Details / Détails : \n" + "\n".join(errors)
    )
