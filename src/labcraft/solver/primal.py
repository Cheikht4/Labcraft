"""Primal solver for equilibrium problems using scipy trust-constr.

Solveur primal pour les problèmes d'équilibre utilisant scipy trust-constr.
Repli (fallback) pour le solveur dual de Newton.
Ref: Dirks et al., SIAM Review 2007, section 4.

Primal solver for equilibrium problems using scipy trust-constr.
Fallback for the dual Newton solver.
"""
from __future__ import annotations

import logging

import numpy as np
import scipy.optimize
from scipy.optimize import Bounds, LinearConstraint

from .types import ConvergenceError, EquilibriumProblem, EquilibriumResult, SolverMethod

logger = logging.getLogger(__name__)


def solve_primal(
    problem: EquilibriumProblem,
    *,
    convergence_threshold: float = 1e-10,
) -> EquilibriumResult:
    """Solve the equilibrium problem by minimizing the total free energy.

    Résout le problème d'équilibre en minimisant l'énergie libre totale.
    Formulation primale : minimise G(x) sous contraintes de conservation
    de masse et de positivité.

    Parameters
    ----------
    problem : EquilibriumProblem
        The equilibrium problem to solve / Le problème d'équilibre à résoudre.
    convergence_threshold : float
        Maximum relative mass conservation residual / Résidu relatif maximal.

    Returns
    -------
    EquilibriumResult
        Solution at equilibrium / Solution à l'équilibre.

    Raises
    ------
    ConvergenceError
        If the solver fails to converge / Si le solveur ne converge pas.
    """
    # Gas constant / Constante des gaz
    r_kcal_mol_k = 1.987e-3
    rt = r_kcal_mol_k * problem.temperature_kelvin

    a_mat = problem.stoichiometry
    dg = problem.delta_g
    xtot = problem.total_concentrations

    n_c = problem.n_complexes
    n_s = problem.n_strands

    # Precomputed ΔG°/RT / ΔG°/RT précalculé
    dg_rt = dg / rt

    # Minimum bound for concentrations (avoid log(0))
    # Borne minimale pour les concentrations (éviter log(0))
    x_min = 1e-30

    # --- Objective: G(x) = Σ_c x_c · (ΔG°_c/RT + ln(x_c) - 1) ---
    # --- Objectif : G(x) = Σ_c x_c · (ΔG°_c/RT + ln(x_c) - 1) ---
    def obj_fun(x: np.ndarray) -> float:
        x_safe = np.maximum(x, x_min)
        return float(np.sum(x_safe * (dg_rt + np.log(x_safe) - 1.0)))

    def obj_jac(x: np.ndarray) -> np.ndarray:
        x_safe = np.maximum(x, x_min)
        return (dg_rt + np.log(x_safe)).astype(np.float64)

    def obj_hess(x: np.ndarray) -> np.ndarray:
        x_safe = np.maximum(x, x_min)
        return np.diag((1.0 / x_safe).astype(np.float64))

    # --- Constraints: A^T · x = x_total ---
    # --- Contraintes : A^T · x = x_total ---
    def mass_eq(x: np.ndarray) -> np.ndarray:
        return a_mat.T @ x - xtot

    def mass_jac(x: np.ndarray) -> np.ndarray:
        return a_mat.T

    mass_constraint = {'type': 'eq', 'fun': mass_eq, 'jac': mass_jac}
    positivity_bounds = [(x_min, None) for _ in range(n_c)]

    # --- Initialization: distribute mass among complexes containing each strand ---
    # --- Initialisation : distribue la masse parmi les complexes contenant chaque brin ---
    # Start with most mass in free strands, small amount in formed complexes
    # On met l'essentiel de la masse dans les brins libres, peu dans les complexes
    x0 = np.full(n_c, x_min)
    for i in range(n_s):
        # Identify the free-strand complex for strand i
        # Identifie le complexe "brin libre" pour le brin i
        for c_idx in range(n_c):
            if (a_mat[c_idx, i] == 1
                    and np.sum(a_mat[c_idx, :]) == 1
                    and dg[c_idx] == 0.0):
                x0[c_idx] = xtot[i] * 0.9  # 90% free / 90% libre
                break
        # Distribute remaining 10% among other complexes containing i
        # Distribue les 10% restants parmi les autres complexes contenant i
        other_complexes = [
            c_idx for c_idx in range(n_c)
            if a_mat[c_idx, i] > 0 and not (
                a_mat[c_idx, i] == 1 and np.sum(a_mat[c_idx, :]) == 1
                and dg[c_idx] == 0.0
            )
        ]
        if other_complexes:
            share = xtot[i] * 0.1 / len(other_complexes)
            for c_idx in other_complexes:
                x0[c_idx] = max(x0[c_idx], share)

    # Solve / Résolution
    try:
        res = scipy.optimize.minimize(
            fun=obj_fun,
            x0=x0,
            method="SLSQP",
            jac=obj_jac,
            constraints=[mass_constraint],
            bounds=positivity_bounds,
            options={
                "maxiter": 1000,
                "ftol": 1e-15,
            },
        )
    except Exception as exc:
        raise ConvergenceError(
            f"Primal solver (SLSQP) raised an exception: {exc}"
        ) from exc

    c_conc = np.maximum(res.x, 0.0)

    # Compute mass conservation residuals / Calcule les résidus de conservation de masse
    g_resid = xtot - a_mat.T @ c_conc
    rel_residuals = np.abs(g_resid) / xtot
    max_res = float(np.max(rel_residuals))

    if max_res > convergence_threshold * 100:
        raise ConvergenceError(
            f"Primal solver failed: max relative residual {max_res:.2e} "
            f"exceeds threshold {convergence_threshold * 100:.2e}. "
            f"scipy message: {res.message}"
        )

    # Extract free strand concentrations / Extrait les concentrations des brins libres
    # Free strands are monomeric complexes with ΔG° = 0
    # Les brins libres sont les complexes monomériques avec ΔG° = 0
    free_conc = np.zeros(n_s)
    for i in range(n_s):
        for c_idx in range(n_c):
            if (dg[c_idx] == 0.0
                    and a_mat[c_idx, i] == 1
                    and np.sum(a_mat[c_idx, :]) == 1):
                free_conc[i] = c_conc[c_idx]
                break

    free_conc = np.maximum(free_conc, x_min)
    u = np.log(free_conc)

    return EquilibriumResult(
        concentrations=c_conc,
        free_concentrations=free_conc,
        log_free_concentrations=u,
        residuals=g_resid,
        max_residual=max_res,
        n_iterations=res.nit,
        method=SolverMethod.TRUST_CONSTR,
        converged=True,
    )
