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

import logging

import numpy as np
import scipy.linalg

from .types import ConvergenceError, EquilibriumProblem, EquilibriumResult, SolverMethod

logger = logging.getLogger(__name__)

# Gas constant in kcal/(mol·K) / Constante des gaz en kcal/(mol·K)
_R_KCAL_MOL_K: float = 1.987e-3

# Maximum log-concentration to prevent exp overflow / Log-concentration max pour éviter overflow
_LOG_CONC_CAP: float = 500.0


def _compute_complex_concentrations(
    u: np.ndarray,
    a_mat: np.ndarray,
    dg_over_rt: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute complex concentrations from free-strand log-concentrations.

    Calcule les concentrations de complexes à partir des log-concentrations
    des brins libres.

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
        Log-concentrations of complexes (capped) / Log-concentrations des complexes.
    c_conc : ndarray, shape (n_complexes,)
        Concentrations of complexes / Concentrations des complexes.
    """
    log_c = -dg_over_rt + a_mat @ u
    log_c_capped = np.clip(log_c, a_min=-_LOG_CONC_CAP, a_max=_LOG_CONC_CAP)
    c_conc = np.exp(log_c_capped)
    return log_c_capped, c_conc


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
    max_iterations: int = 200,
    convergence_threshold: float = 1e-10,
    armijo_c: float = 1e-4,
    armijo_rho: float = 0.5,
    initial_step_size: float = 1.0,
    precondition: bool = True,
    stagnation_rtol: float = 1e-14,
) -> EquilibriumResult:
    """Solve the equilibrium problem using the damped dual Newton method.

    Résout le problème d'équilibre par la méthode de Newton duale amortie.

    The dual function D(u) = Σ_i x_i^total · u_i - Σ_c exp(-ΔG°_c/RT + Σ_i A[c,i]·u_i)
    is concave. Its gradient equals the mass conservation residual, and its
    Hessian is negative definite. We maximize D(u) using damped Newton with
    Armijo backtracking line search.

    La fonction duale D(u) est concave. Son gradient est le résidu de
    conservation de masse, et sa hessienne est définie négative. On maximise
    D(u) par Newton amorti avec recherche linéaire d'Armijo.

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
    stagnation_rtol : float
        Relative tolerance for detecting stagnation. If the residual doesn't
        improve by more than this factor, accept the current solution if
        the residual is within 10x of the threshold.

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
    # Initialisation : on suppose la moitié de chaque brin libre
    u = np.log(xtot) - np.log(2.0)

    prev_max_res = np.inf
    stagnation_count = 0
    max_stagnation = 10  # Accept after N stagnant iterations / Accepte après N itérations stagnantes

    for it in range(max_iterations):
        # Compute complex concentrations / Calcule les concentrations de complexes
        _, c_conc = _compute_complex_concentrations(u, a_mat, dg_over_rt)

        # Gradient (= mass conservation residual)
        # Gradient (= résidu de conservation de masse)
        g = xtot - a_mat.T @ c_conc

        # Relative residual / Résidu relatif
        rel_residuals = np.abs(g) / xtot
        max_res = float(np.max(rel_residuals))

        # Check convergence / Vérifie la convergence
        if max_res < convergence_threshold:
            return _build_result(u, c_conc, g, max_res, it + 1)

        # Detect stagnation: if residual barely improves, the solver has
        # reached machine precision for this problem. Accept if close enough.
        # Détecte la stagnation : si le résidu ne s'améliore plus, le solveur
        # a atteint la précision machine. Accepte si suffisamment proche.
        improvement = abs(prev_max_res - max_res) / max(prev_max_res, 1e-300)
        if improvement < stagnation_rtol:
            stagnation_count += 1
            if stagnation_count >= max_stagnation:
                # Accept if within 10x of threshold
                # Accepte si dans un facteur 10 du seuil
                if max_res < convergence_threshold * 10:
                    logger.info(
                        f"Dual Newton stagnated at residual {max_res:.2e} "
                        f"(threshold={convergence_threshold:.0e}). "
                        f"Accepting as converged. / Stagnation, accepté."
                    )
                    return _build_result(u, c_conc, g, max_res, it + 1)
                break  # Stagnated but too far — let fallback handle it
        else:
            stagnation_count = 0
        prev_max_res = max_res

        # Hessian: H[i,j] = -Σ_c A[c,i]·A[c,j]·[c]  (negative definite)
        # Hessienne : H[i,j] = -Σ_c A[c,i]·A[c,j]·[c]  (définie négative)
        # -H = A^T · diag(c) · A  (positive definite)
        neg_h = (a_mat.T * c_conc) @ a_mat

        # Solve Newton system: (-H) · Δu = g
        # Résout le système Newton : (-H) · Δu = g
        if precondition:
            # Diagonal preconditioning: D = diag(1/sqrt(diag(-H)))
            # Préconditionnement diagonal : D = diag(1/sqrt(diag(-H)))
            diag_neg_h = np.diag(neg_h)
            diag_neg_h = np.maximum(diag_neg_h, 1e-300)
            d_vec = 1.0 / np.sqrt(diag_neg_h)

            # Preconditioned system: (D·(-H)·D) · δ = D·g, then Δu = D·δ
            h_precond = d_vec[:, np.newaxis] * neg_h * d_vec[np.newaxis, :]
            g_precond = d_vec * g

            try:
                cho, lower = scipy.linalg.cho_factor(h_precond)
                delta_u = d_vec * scipy.linalg.cho_solve((cho, lower), g_precond)
            except scipy.linalg.LinAlgError:
                # Fallback to general solve / Repli sur résolution générale
                logger.debug("Cholesky failed, falling back to general solve")
                delta_u = np.linalg.solve(neg_h, g)
        else:
            try:
                cho, lower = scipy.linalg.cho_factor(neg_h)
                delta_u = scipy.linalg.cho_solve((cho, lower), g)
            except scipy.linalg.LinAlgError:
                delta_u = np.linalg.solve(neg_h, g)

        # Armijo backtracking line search / Recherche linéaire d'Armijo
        # Find α such that D(u + α·Δu) ≥ D(u) + c·α·g^T·Δu
        g_dot_du = float(np.dot(g, delta_u))

        # Skip if Newton direction is not ascent / Passe si pas direction de montée
        if g_dot_du <= 0:
            logger.debug(f"Newton direction is not ascent: g·Δu = {g_dot_du:.2e}")
            break

        # Evaluate dual at current point / Évalue le dual au point courant
        d_current = float(np.sum(xtot * u) - np.sum(c_conc))

        alpha = initial_step_size
        min_alpha = 1e-15
        while alpha > min_alpha:
            u_next = u + alpha * delta_u
            _, c_next = _compute_complex_concentrations(u_next, a_mat, dg_over_rt)
            d_next = float(np.sum(xtot * u_next) - np.sum(c_next))
            if d_next >= d_current + armijo_c * alpha * g_dot_du:
                break
            alpha *= armijo_rho

        if alpha <= min_alpha:
            logger.debug(f"Armijo step too small at iteration {it}")

        u = u + alpha * delta_u

    raise ConvergenceError(
        f"Dual Newton failed to converge after {it + 1} iterations. "
        f"Max relative residual: {max_res:.2e} "
        f"(threshold: {convergence_threshold:.0e})"
    )
