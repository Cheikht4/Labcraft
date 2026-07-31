"""Types for the equilibrium solver / Types pour le solveur d'équilibre."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import numpy as np


class ConvergenceError(Exception):
    """Raised when the solver fails to converge.
    
    Levée quand le solveur ne parvient pas à converger.
    """
    pass


class SolverMethod(Enum):
    """Methods for solving equilibrium problems / Méthodes de résolution."""
    ANALYTICAL = "analytical"
    DUAL_NEWTON = "dual_newton"
    TRUST_CONSTR = "trust_constr"
    EXTENDED_PRECISION = "extended_precision"


@dataclass(frozen=True)
class EquilibriumProblem:
    """Specification of a multi-state equilibrium problem.
    
    Spécification d'un problème d'équilibre multi-états.
    
    Attributes:
        n_strands: Number of distinct strand species / Nombre d'espèces de brins.
        n_complexes: Number of complexes (including free strands) / Nombre de complexes.
        stoichiometry: (n_complexes, n_strands) matrix A[c,i] / Matrice stœchiométrique.
        delta_g: (n_complexes,) vector of ΔG° in kcal/mol / Vecteur des ΔG°.
        total_concentrations: (n_strands,) vector of total concentrations in M /
            Concentrations totales en M.
        temperature_kelvin: Temperature in Kelvin / Température en Kelvin.
        max_complex_size: Maximum number of strands per complex / Taille max des complexes.
    """
    n_strands: int
    n_complexes: int
    stoichiometry: np.ndarray  # shape (n_complexes, n_strands)
    delta_g: np.ndarray  # shape (n_complexes,), kcal/mol
    total_concentrations: np.ndarray  # shape (n_strands,), mol/L
    temperature_kelvin: float
    max_complex_size: int = 2

    def validate(self) -> None:
        """Validate problem dimensions and constraints.
        
        Valide les dimensions et les contraintes du problème.
        
        Raises:
            ValueError: If constraints are violated / Si des contraintes sont violées.
        """
        if self.stoichiometry.shape != (self.n_complexes, self.n_strands):
            raise ValueError(
                f"Stoichiometry shape mismatch: expected ({self.n_complexes}, {self.n_strands}), "
                f"got {self.stoichiometry.shape}"
            )
        if self.delta_g.shape != (self.n_complexes,):
            raise ValueError(
                f"Delta G shape mismatch: expected ({self.n_complexes},), got {self.delta_g.shape}"
            )
        if self.total_concentrations.shape != (self.n_strands,):
            raise ValueError(
                f"Total concentrations shape mismatch: expected ({self.n_strands},), "
                f"got {self.total_concentrations.shape}"
            )
            
        if np.any(self.total_concentrations <= 0):
            raise ValueError("All total concentrations must be > 0")
            
        if self.temperature_kelvin <= 0:
            raise ValueError("Temperature must be > 0 K")
            
        if np.any(self.stoichiometry < 0):
            raise ValueError("Stoichiometry matrix must be non-negative")
            
        if not np.all(self.stoichiometry == np.floor(self.stoichiometry)):
            raise ValueError("Stoichiometry matrix must contain integers")
            
        # Check that each strand has at least one complex containing it
        # Vérifie que chaque brin a au moins un complexe le contenant
        for i in range(self.n_strands):
            if np.all(self.stoichiometry[:, i] == 0):
                raise ValueError(f"Strand {i} is not present in any complex")


@dataclass(frozen=True)
class EquilibriumResult:
    """Result of equilibrium calculation.
    
    Résultat du calcul d'équilibre.
    """
    concentrations: np.ndarray  # shape (n_complexes,), mol/L
    free_concentrations: np.ndarray  # shape (n_strands,), mol/L
    log_free_concentrations: np.ndarray  # shape (n_strands,), ln([x_i]_libre)
    residuals: np.ndarray  # shape (n_strands,), mass conservation residual
    max_residual: float  # max absolute relative residual
    n_iterations: int
    method: SolverMethod
    converged: bool
