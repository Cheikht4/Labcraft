"""Analytical solutions for simple equilibrium systems.

Solutions analytiques pour les systèmes d'équilibre simples.
"""
from __future__ import annotations
from labcraft.thermo.constants import R_GAS_KCAL_MOL_K

import numpy as np
from .types import EquilibriumResult, SolverMethod


def _get_result(
    a0: float, 
    b0: float, 
    ab_conc: float, 
    method_name: str = "analytical"
) -> EquilibriumResult:
    """Helper to build the EquilibriumResult for a 2-strand system.
    
    Aide pour construire l'EquilibriumResult pour un système à 2 brins.
    """
    free_a = max(0.0, a0 - ab_conc)
    free_b = max(0.0, b0 - ab_conc)
    
    concentrations = np.array([free_a, free_b, ab_conc], dtype=np.float64)
    free_concentrations = np.array([free_a, free_b], dtype=np.float64)
    
    # Avoid log(0)
    log_free_concentrations = np.log(np.maximum(free_concentrations, 1e-300))
    
    residuals = np.array([
        (free_a + ab_conc) - a0,
        (free_b + ab_conc) - b0
    ], dtype=np.float64)
    
    max_residual = max(abs(residuals[0]) / a0, abs(residuals[1]) / b0)
    
    return EquilibriumResult(
        concentrations=concentrations,
        free_concentrations=free_concentrations,
        log_free_concentrations=log_free_concentrations,
        residuals=residuals,
        max_residual=float(max_residual),
        n_iterations=0,
        method=SolverMethod.ANALYTICAL,
        converged=True,
    )


def analytical_two_strand(
    a0: float, 
    b0: float, 
    delta_g_kcal: float, 
    temperature_kelvin: float,
    *, 
    use_stable_form: bool = True,
) -> EquilibriumResult:
    """Analytical solution for a 2-strand system A + B ⇌ AB.
    
    Solution analytique pour un système à 2 brins A + B ⇌ AB.
    
    Args:
        a0: Total concentration of A (mol/L).
        b0: Total concentration of B (mol/L).
        delta_g_kcal: Standard free energy of binding (kcal/mol).
        temperature_kelvin: Temperature in Kelvin.
        use_stable_form: If True, uses the numerically stable form.
        
    Returns:
        EquilibriumResult containing the concentrations.
    """
    r = R_GAS_KCAL_MOL_K  # kcal/(mol K)
    rt = r * temperature_kelvin
    k_eq = np.exp(-delta_g_kcal / rt)
    
    s = a0 + b0 + 1.0 / k_eq
    
    if use_stable_form:
        # Stable form: 2*a0*b0 / (S + sqrt(S^2 - 4*a0*b0))
        # Forme stable
        denominator = s + np.sqrt(max(0.0, s**2 - 4 * a0 * b0))
        ab_conc = (2 * a0 * b0) / denominator if denominator > 0 else 0.0
    else:
        # Naive form: (S - sqrt(S^2 - 4*a0*b0)) / 2
        # Forme naïve
        ab_conc = (s - np.sqrt(max(0.0, s**2 - 4 * a0 * b0))) / 2.0
        
    return _get_result(a0, b0, ab_conc)


def analytical_two_strand_naive(
    a0: float, 
    b0: float, 
    delta_g_kcal: float, 
    temperature_kelvin: float,
) -> EquilibriumResult:
    """Naive analytical solution for a 2-strand system.
    
    Solution analytique naïve pour un système à 2 brins.
    """
    return analytical_two_strand(
        a0, b0, delta_g_kcal, temperature_kelvin, use_stable_form=False
    )
