"""Extended precision solver for equilibrium problems.

Solveur à précision étendue pour les problèmes d'équilibre.
"""
from __future__ import annotations

import numpy as np
from .types import EquilibriumProblem, EquilibriumResult, SolverMethod, ConvergenceError

try:
    import mpmath
except ImportError:
    mpmath = None  # type: ignore[assignment]


def solve_extended(
    problem: EquilibriumProblem,
    *,
    precision_digits: int = 50,
    convergence_threshold: float = 1e-10,
    max_iterations: int = 200,
) -> EquilibriumResult:
    """Solve the equilibrium problem using extended precision (mpmath).
    
    Résout le problème d'équilibre en utilisant la précision étendue (mpmath).
    """
    if mpmath is None:
        raise ImportError(
            "mpmath is required for the extended precision solver. "
            "Please install it with: pip install mpmath"
        )
        
    mpmath.mp.dps = precision_digits
    
    r_kcal_mol_k = mpmath.mpf('1.987e-3')
    rt = r_kcal_mol_k * mpmath.mpf(problem.temperature_kelvin)
    
    a_mat = problem.stoichiometry
    dg = [mpmath.mpf(x) for x in problem.delta_g]
    xtot = [mpmath.mpf(x) for x in problem.total_concentrations]
    
    n_s = problem.n_strands
    n_c = problem.n_complexes
    
    # Initialisation
    u = [mpmath.log(xtot[i]) - mpmath.log(2.0) for i in range(n_s)]
    
    for it in range(max_iterations):
        # log_c = -ΔG/RT + A * u
        log_c = []
        for c in range(n_c):
            sum_au = mpmath.mpf(0.0)
            for i in range(n_s):
                if a_mat[c, i] > 0:
                    sum_au += mpmath.mpf(a_mat[c, i]) * u[i]
            log_c.append(-dg[c] / rt + sum_au)
            
        # Concentrations
        c_conc = [mpmath.exp(lc) for lc in log_c]
        
        # Gradient
        g = []
        for i in range(n_s):
            sum_ac = mpmath.mpf(0.0)
            for c in range(n_c):
                if a_mat[c, i] > 0:
                    sum_ac += mpmath.mpf(a_mat[c, i]) * c_conc[c]
            g.append(xtot[i] - sum_ac)
            
        # Max relative residual
        max_res = max([abs(g[i]) / xtot[i] for i in range(n_s)])
        
        if float(max_res) < convergence_threshold:
            free_conc_np = np.array([float(mpmath.exp(ui)) for ui in u], dtype=np.float64)
            c_conc_np = np.array([float(cc) for cc in c_conc], dtype=np.float64)
            g_np = np.array([float(gi) for gi in g], dtype=np.float64)
            u_np = np.array([float(ui) for ui in u], dtype=np.float64)
            
            return EquilibriumResult(
                concentrations=c_conc_np,
                free_concentrations=free_conc_np,
                log_free_concentrations=u_np,
                residuals=g_np,
                max_residual=float(max_res),
                n_iterations=it + 1,
                method=SolverMethod.EXTENDED_PRECISION,
                converged=True
            )
            
        # Hessian H = - A^T * diag(c) * A
        h = mpmath.matrix(n_s, n_s)
        for i in range(n_s):
            for j in range(n_s):
                sum_h = mpmath.mpf(0.0)
                for c in range(n_c):
                    if a_mat[c, i] > 0 and a_mat[c, j] > 0:
                        sum_h -= mpmath.mpf(a_mat[c, i] * a_mat[c, j]) * c_conc[c]
                h[i, j] = sum_h
                
        # Newton step (-H) * du = g
        neg_h = -h
        g_vec = mpmath.matrix(n_s, 1)
        for i in range(n_s):
            g_vec[i, 0] = g[i]
            
        try:
            delta_u = mpmath.lu_solve(neg_h, g_vec)
        except ZeroDivisionError:
            raise ConvergenceError("Singular Hessian in extended precision solver")
            
        for i in range(n_s):
            u[i] += delta_u[i, 0]
            
    raise ConvergenceError(
        f"Extended precision solver failed to converge after {max_iterations} iterations. "
        f"Max residual: {float(max_res):.2e}"
    )
