"""Dual Newton solver for multi-state equilibrium problems.

Solveur Newton sur le dual pour les problèmes d'équilibre multi-états.

This is the primary solver. It maximizes the concave dual function D(u)
whose variables are u_i = ln([x_i]_free), the log-concentrations of free
strands. At the optimum, the gradient equals the mass conservation residual.

C'est le solveur principal. Il maximise la fonction duale concave D(u)
dont les variables sont u_i = ln([x_i]_libre). À l'optimum, le gradient
est le résidu de conservation de masse.

References
----------
Dirks R.M., Bois J.S., Schaeffer J.M., Winfree E., Pierce N.A. (2007),
Thermodynamic analysis of interacting nucleic acid strands,
SIAM Review 49(1):65-88, Section 4.
"""
from __future__ import annotations
from labcraft.thermo.constants import R_GAS_KCAL_MOL_K

import logging

import numpy as np
import scipy.linalg

from .types import ConvergenceError, EquilibriumProblem, EquilibriumResult, SolverMethod

logger = logging.getLogger(__name__)

# Gas constant in kcal/(mol·K) / Constante des gaz en kcal/(mol·K)
_R_KCAL_MOL_K: float = R_GAS_KCAL_MOL_K


def _compute_complex_concentrations(
    u: np.ndarray,
    a_mat: np.ndarray,
    dg_over_rt: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, bool]:
    """Compute complex concentrations from free-strand log-concentrations.

    Calcule les concentrations de complexes à partir des log-concentrations
    des brins libres. Ne masque plus silencieusement les débordements,
    mais signale si un débordement a eu lieu (overflow).

    Parameters
    ----------
    u : ndarray, shape (n_strands,)
        Log free-strand concentrations / Log-concentrations des brins libres.
    a_mat : ndarray, shape (n_complexes, n_strands)
        Stoichiometry matrix / Matrice stœchiométrique.
    dg_over_rt : ndarray, shape (n_complexes,)
        ΔG°/RT for each complex / ΔG°/RT pour chaque complexe.

    Returns
    -------
    log_c : ndarray, shape (n_complexes,)
        Log-concentrations of complexes / Log-concentrations des complexes.
    c_conc : ndarray, shape (n_complexes,)
        Concentrations of complexes / Concentrations des complexes.
    overflow : bool
        True if any log_c > 500 (would cause exp overflow).
    """
    log_c = -dg_over_rt + a_mat @ u
    
    # Check for overflow
    overflow = bool(np.any(log_c > 500.0))
    
    if overflow:
        # Cap only to prevent actual floating point exception,
        # but the caller is informed that the point is invalid.
        log_c_safe = np.clip(log_c, a_min=None, a_max=500.0)
        c_conc = np.exp(log_c_safe)
    else:
        c_conc = np.exp(log_c)
        
    return log_c, c_conc, overflow


def _build_result(
    u: np.ndarray,
    c_conc: np.ndarray,
    g: np.ndarray,
    max_res: float,
    n_iterations: int,
) -> EquilibriumResult:
    """Build an EquilibriumResult from the solver state.

    Construit un EquilibriumResult à partir de l'état du solveur.
    """
    free_conc = np.exp(u)
    return EquilibriumResult(
        concentrations=c_conc,
        free_concentrations=free_conc,
        log_free_concentrations=u.copy(),
        residuals=g.copy(),
        max_residual=max_res,
        n_iterations=n_iterations,
        method=SolverMethod.DUAL_NEWTON,
        converged=True,
    )


def solve_dual(
    problem: EquilibriumProblem,
    *,
    max_iterations: int = 500,
    convergence_threshold: float = 1e-10,
    armijo_c: float = 1e-4,
    armijo_rho: float = 0.5,
    initial_step_size: float = 1.0,
    precondition: bool = True,
) -> EquilibriumResult:
    """Solve the equilibrium problem using the damped dual Newton method.

    Résout le problème d'équilibre par la méthode de Newton duale amortie,
    avec régularisation de Levenberg-Marquardt dans l'espace préconditionné.
    Utilise une globalisation "Newton Complet d'Abord" pour garantir la 
    convergence quadratique finale.

    Parameters
    ----------
    problem : EquilibriumProblem
        Problem specification / Spécification du problème.
    max_iterations : int
        Maximum Newton iterations / Itérations Newton maximales.
    convergence_threshold : float
        Max relative residual for convergence / Résidu relatif max pour convergence.
    armijo_c : float
        Armijo sufficient decrease parameter / Paramètre de décroissance suffisante.
    armijo_rho : float
        Step size reduction factor / Facteur de réduction du pas.
    initial_step_size : float
        Initial step size α₀ / Pas initial α₀.
    precondition : bool
        If True, apply diagonal preconditioning / Si True, préconditionnement diagonal.

    Returns
    -------
    EquilibriumResult
        Equilibrium concentrations / Concentrations à l'équilibre.

    Raises
    ------
    ConvergenceError
        If convergence is not achieved / Si la convergence n'est pas atteinte.
    """
    rt = _R_KCAL_MOL_K * problem.temperature_kelvin
    a_mat = problem.stoichiometry
    dg_over_rt = problem.delta_g / rt
    xtot = problem.total_concentrations

    # Initialization: assume half of each strand is free
    u = np.log(xtot) - np.log(2.0)
    
    lambda_lm = 1e-12

    for it in range(max_iterations):
        # Compute complex concentrations
        _, c_conc, overflow = _compute_complex_concentrations(u, a_mat, dg_over_rt)
        
        if overflow:
            raise ConvergenceError(
                f"Overflow encountered outside line search at iteration {it}."
            )

        # Gradient (= mass conservation residual)
        g = xtot - a_mat.T @ c_conc
        g_norm_sq = float(np.sum(g**2))

        # Relative residual
        rel_residuals = np.abs(g) / xtot
        max_res = float(np.max(rel_residuals))

        # Check convergence
        if max_res < convergence_threshold:
            return _build_result(u, c_conc, g, max_res, it + 1)

        # Hessian: H[i,j] = -Σ_c A[c,i]·A[c,j]·[c]  (negative definite)
        # -H = A^T · diag(c) · A  (positive definite)
        neg_h = (a_mat.T * c_conc) @ a_mat

        if precondition:
            # Diagonal preconditioning: D = diag(1/sqrt(diag(-H)))
            diag_neg_h = np.diag(neg_h)
            diag_neg_h = np.maximum(diag_neg_h, 1e-300)
            d_vec = 1.0 / np.sqrt(diag_neg_h)
            
            h_precond = d_vec[:, np.newaxis] * neg_h * d_vec[np.newaxis, :]
            g_precond = d_vec * g
        else:
            h_precond = neg_h
            g_precond = g
            d_vec = np.ones(problem.n_strands)

        # ---------------------------------------------------------
        # FAST PATH: Try pure full Newton step first (lambda=0, alpha=1)
        # ---------------------------------------------------------
        pure_success = False
        try:
            cho, lower = scipy.linalg.cho_factor(h_precond)
            delta_precond = scipy.linalg.cho_solve((cho, lower), g_precond)
            delta_u_pure = d_vec * delta_precond
            
            # Evaluate full step
            u_next_pure = u + delta_u_pure
            _, c_next_pure, next_overflow_pure = _compute_complex_concentrations(
                u_next_pure, a_mat, dg_over_rt
            )
            
            if not next_overflow_pure:
                with np.errstate(over='ignore', invalid='ignore'):
                    g_next_pure = xtot - a_mat.T @ c_next_pure
                    g_next_norm_sq = float(np.sum(g_next_pure**2))
                    if not np.isfinite(g_next_norm_sq):
                        g_next_norm_sq = np.inf
                
                if g_next_norm_sq < g_norm_sq:
                    # Pure step accepted!
                    u = u_next_pure
                    lambda_lm = 0.0
                    pure_success = True
                    
        except scipy.linalg.LinAlgError:
            pass
            
        if pure_success:
            continue
            
        # ---------------------------------------------------------
        # FALLBACK PATH: Levenberg-Marquardt with Cholesky retries
        # ---------------------------------------------------------
        lambda_lm = max(lambda_lm, 1e-12)
        cholesky_success = False
        
        for _ in range(20):  # max retries for LM
            if precondition:
                h_damped = h_precond + lambda_lm * np.eye(problem.n_strands)
            else:
                h_damped = neg_h + lambda_lm * np.eye(problem.n_strands)

            try:
                cho, lower = scipy.linalg.cho_factor(h_damped)
                delta_precond = scipy.linalg.cho_solve((cho, lower), g_precond)
                delta_u = d_vec * delta_precond
                cholesky_success = True
                break
            except scipy.linalg.LinAlgError:
                lambda_lm = max(lambda_lm * 10, 1e-12)
        
        if not cholesky_success:
            raise ConvergenceError(
                f"Cholesky factorization failed even with LM damping at iteration {it}."
            )

        # Pre-compute values for defensive Dual Armijo if needed
        g_dot_du = float(np.dot(g, delta_u))
        if g_dot_du <= 0:
            logger.debug(f"Newton direction is not ascent: g·Δu = {g_dot_du:.2e}")
            break
            
        d_current = float(np.sum(xtot * u) - np.sum(c_conc))

        alpha = initial_step_size
        min_alpha = 1e-15
        step_accepted = False
        
        # 1. Line search on Residual Norm ||g||^2
        while alpha > min_alpha:
            u_next = u + alpha * delta_u
            _, c_next, next_overflow = _compute_complex_concentrations(u_next, a_mat, dg_over_rt)
            
            if next_overflow:
                alpha *= armijo_rho
                continue
                
            with np.errstate(over='ignore', invalid='ignore'):
                g_next = xtot - a_mat.T @ c_next
                g_next_norm_sq = float(np.sum(g_next**2))
                if not np.isfinite(g_next_norm_sq):
                    g_next_norm_sq = np.inf
            
            if g_next_norm_sq < g_norm_sq:
                step_accepted = True
                break
                
            alpha *= armijo_rho

        # 2. Defensive fallback: Line search on Dual Objective D(u)
        if not step_accepted:
            alpha = initial_step_size
            while alpha > min_alpha:
                u_next = u + alpha * delta_u
                _, c_next, next_overflow = _compute_complex_concentrations(u_next, a_mat, dg_over_rt)
                
                if next_overflow:
                    alpha *= armijo_rho
                    continue
                    
                d_next = float(np.sum(xtot * u_next) - np.sum(c_next))
                if d_next >= d_current + armijo_c * alpha * g_dot_du:
                    step_accepted = True
                    break
                    
                alpha *= armijo_rho

        if not step_accepted:
            logger.debug(f"Armijo step too small at iteration {it}")
            # If the step fails, we might just be blocked. Increase LM damping.
            lambda_lm = max(lambda_lm * 10, 1e-12)
        else:
            # Step accepted: update u
            u = u + alpha * delta_u
            
            # Update LM damping based on step size
            if alpha >= 0.1:
                if lambda_lm > 1e-12:
                    lambda_lm /= 10
                else:
                    lambda_lm = 0.0
            elif alpha < 1e-3:
                lambda_lm = max(lambda_lm * 10, 1e-12)

    raise ConvergenceError(
        f"Dual Newton failed to converge after {it + 1} iterations. "
        f"Max relative residual: {max_res:.2e} "
        f"(threshold: {convergence_threshold:.0e})"
    )
