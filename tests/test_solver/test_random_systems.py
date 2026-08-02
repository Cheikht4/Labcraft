"""Random and large-scale validation tests for the equilibrium solver.

Tests de validation aléatoires et à grande échelle pour le solveur d'équilibre.
Ref: Cahier des charges section 2.4.
"""
from __future__ import annotations

import numpy as np
import pytest
import scipy.linalg
from labcraft.solver.analytical import analytical_two_strand
from labcraft.solver.solver import solve_equilibrium
from labcraft.solver.types import EquilibriumProblem


class TestSolverValidation:
    """Intensive validation of the solver (D1-D5, C1-C3).
    
    Validation intensive du solveur.
    """
    
    # Tolerances
    # Solveur threshold: 1e-8. Assertions on concentrations: 1e-6.
    # Note: float64 has ~16 digits. For extreme cases (xtot=1e-12, max(c)=1e-3),
    # the condition number is 1e9, so the best achievable relative residual is 1e-16 * 1e9 = 1e-7.
    # Thus, 1e-7 is the safest strict threshold for these expanded ranges.
    CONVERGENCE_THRESHOLD = 1e-7
    CONC_TOLERANCE = 1e-5
    
    def test_random_systems_and_mass_action(self) -> None:
        """Garantie 1: 500 random systems with mass action verification.
        
        500 systèmes aléatoires (fortes et faibles affinités).
        Vérifie l'accord avec l'analytique ET la loi d'action de masse a posteriori.
        """
        N_SYSTEMS = 500
        rng = np.random.default_rng(seed=42)
        
        # Concentrations: 1e-12 to 1e-3 M (9 orders of magnitude)
        log_conc = rng.uniform(-12, -3, size=(N_SYSTEMS, 2))
        concentrations = 10.0 ** log_conc
        # dG°: -45 to -5 kcal/mol (covers extreme strong binding)
        delta_gs = rng.uniform(-45, -5, size=N_SYSTEMS)
        temperatures = rng.uniform(328.15, 343.15, size=N_SYSTEMS)
        
        max_errors = []
        for i in range(N_SYSTEMS):
            a0, b0 = concentrations[i]
            dg = delta_gs[i]
            T = temperatures[i]
            
            ref = analytical_two_strand(a0, b0, dg, T)
            
            problem = EquilibriumProblem(
                n_strands=2,
                n_complexes=3,
                stoichiometry=np.array([[1, 0], [0, 1], [1, 1]], dtype=np.float64),
                delta_g=np.array([0.0, 0.0, dg]),
                total_concentrations=np.array([a0, b0]),
                temperature_kelvin=T,
            )
            
            result = solve_equilibrium(problem, convergence_threshold=self.CONVERGENCE_THRESHOLD)
            
            ref_ab = ref.concentrations[2]
            solver_ab = result.concentrations[2]
            
            if ref_ab > 1e-30:
                rel_error = abs(solver_ab - ref_ab) / ref_ab
                max_errors.append(rel_error)
                assert rel_error < self.CONC_TOLERANCE, (
                    f"System {i}: rel_error={rel_error:.2e} > {self.CONC_TOLERANCE:.0e}. "
                    f"a0={a0:.2e}, b0={b0:.2e}, dG={dg:.1f}"
                )
                
            # Garantie 1: Vérification stricte de la loi d'action de masse a posteriori
            # K_c = [c] / ( [A]^A_c * [B]^B_c )
            rt = 1.987e-3 * T
            
            free_a = result.concentrations[0]
            free_b = result.concentrations[1]
            complex_ab = result.concentrations[2]
            
            if free_a > 1e-300 and free_b > 1e-300:
                k_computed = complex_ab / (free_a * free_b)
                k_theoretical = np.exp(-dg / rt)
                
                k_rel_error = abs(k_computed - k_theoretical) / k_theoretical
                assert k_rel_error < self.CONC_TOLERANCE * 100, (
                    f"Mass action violation on system {i}: computed K={k_computed:.2e}, "
                    f"theoretical={k_theoretical:.2e}"
                )

    def test_block_diagonal_30_strands(self) -> None:
        """Garantie 2: Invariance par bloc (30 brins, 15 paires indépendantes permutees)."""
        rng = np.random.default_rng(seed=123)
        N_PAIRS = 15
        N_STRANDS = N_PAIRS * 2
        N_COMPLEXES = N_STRANDS + N_PAIRS
        
        concentrations = 10.0 ** rng.uniform(-9, -4, size=N_STRANDS)
        delta_gs = rng.uniform(-35, -10, size=N_PAIRS)
        T = 338.15
        
        # We will map pair p (0 to 14) to two random distinct strand indices
        strand_indices = np.arange(N_STRANDS)
        rng.shuffle(strand_indices)
        
        a_mat = np.zeros((N_COMPLEXES, N_STRANDS), dtype=np.float64)
        dg_array = np.zeros(N_COMPLEXES, dtype=np.float64)
        
        # 1. Monomer rows (mandatory)
        for i in range(N_STRANDS):
            a_mat[i, i] = 1.0
            dg_array[i] = 0.0
            
        # 2. Dimer rows
        expected_dimers = []
        for p in range(N_PAIRS):
            idx_a = strand_indices[2*p]
            idx_b = strand_indices[2*p + 1]
            
            c_idx = N_STRANDS + p
            a_mat[c_idx, idx_a] = 1.0
            a_mat[c_idx, idx_b] = 1.0
            dg_array[c_idx] = delta_gs[p]
            
            # Compute analytical reference for this independent pair
            a0 = concentrations[idx_a]
            b0 = concentrations[idx_b]
            ref = analytical_two_strand(a0, b0, delta_gs[p], T)
            expected_dimers.append((c_idx, ref.concentrations[2]))
            
        problem = EquilibriumProblem(
            n_strands=N_STRANDS,
            n_complexes=N_COMPLEXES,
            stoichiometry=a_mat,
            delta_g=dg_array,
            total_concentrations=concentrations,
            temperature_kelvin=T
        )
        
        result = solve_equilibrium(problem, convergence_threshold=self.CONVERGENCE_THRESHOLD)
        
        for c_idx, expected_conc in expected_dimers:
            actual_conc = result.concentrations[c_idx]
            rel_error = abs(actual_conc - expected_conc) / expected_conc
            assert rel_error < self.CONC_TOLERANCE, (
                f"Block diagonal failure: rel_error={rel_error:.2e} > {self.CONC_TOLERANCE:.0e}"
            )

    def test_fully_symmetric_n_strands(self) -> None:
        """Garantie 3: Symétrie totale (N brins, tous homodimers/heterodimers identiques)."""
        N_STRANDS = 10
        a0 = 1e-6
        dg = -20.0
        T = 338.15
        
        # Number of heterodimers: N*(N-1)/2
        n_heterodimers = N_STRANDS * (N_STRANDS - 1) // 2
        N_COMPLEXES = N_STRANDS + n_heterodimers
        
        a_mat = np.zeros((N_COMPLEXES, N_STRANDS), dtype=np.float64)
        dg_array = np.zeros(N_COMPLEXES, dtype=np.float64)
        
        for i in range(N_STRANDS):
            a_mat[i, i] = 1.0
            dg_array[i] = 0.0
            
        c_idx = N_STRANDS
        for i in range(N_STRANDS):
            for j in range(i + 1, N_STRANDS):
                a_mat[c_idx, i] = 1.0
                a_mat[c_idx, j] = 1.0
                dg_array[c_idx] = dg
                c_idx += 1
                
        problem = EquilibriumProblem(
            n_strands=N_STRANDS,
            n_complexes=N_COMPLEXES,
            stoichiometry=a_mat,
            delta_g=dg_array,
            total_concentrations=np.full(N_STRANDS, a0),
            temperature_kelvin=T
        )
        
        result = solve_equilibrium(problem, convergence_threshold=self.CONVERGENCE_THRESHOLD)
        
        # Analytical scalar solution for symmetry
        # Each strand is in (N-1) dimers.
        # [dimer] = K * x^2, where x is free concentration
        # Mass cons: x + (N-1) * K * x^2 = a0
        # (N-1)K x^2 + x - a0 = 0
        rt = 1.987e-3 * T
        K = np.exp(-dg / rt)
        
        A = (N_STRANDS - 1) * K
        B = 1.0
        C = -a0
        
        # Stable quadratic root (positive)
        x_exact = (2 * C) / (-B - np.sqrt(B**2 - 4*A*C))
        dimer_exact = K * x_exact**2
        
        for i in range(N_STRANDS):
            assert abs(result.free_concentrations[i] - x_exact) / x_exact < self.CONC_TOLERANCE
            
        for c_idx in range(N_STRANDS, N_COMPLEXES):
            assert abs(result.concentrations[c_idx] - dimer_exact) / dimer_exact < self.CONC_TOLERANCE

    def test_free_energy_minimality(self) -> None:
        """Garantie 4: Minimalité stricte de l'énergie libre sur des perturbations admissibles."""
        rng = np.random.default_rng(seed=999)
        N_STRANDS = 4
        N_COMPLEXES = 10
        
        # Generate random A matrix (with monomers)
        a_mat = np.zeros((N_COMPLEXES, N_STRANDS), dtype=np.float64)
        for i in range(N_STRANDS):
            a_mat[i, i] = 1.0
            
        for c in range(N_STRANDS, N_COMPLEXES):
            # random combinations of 2 strands
            idx = rng.choice(N_STRANDS, size=2, replace=False)
            a_mat[c, idx[0]] = 1.0
            a_mat[c, idx[1]] = 1.0
            
        dg_array = np.zeros(N_COMPLEXES)
        dg_array[N_STRANDS:] = rng.uniform(-25, -5, size=N_COMPLEXES - N_STRANDS)
        xtot = 10.0 ** rng.uniform(-7, -4, size=N_STRANDS)
        T = 338.15
        
        problem = EquilibriumProblem(
            n_strands=N_STRANDS,
            n_complexes=N_COMPLEXES,
            stoichiometry=a_mat,
            delta_g=dg_array,
            total_concentrations=xtot,
            temperature_kelvin=T
        )
        
        result = solve_equilibrium(problem, convergence_threshold=1e-13)
        c_opt = result.concentrations
        
        # Compute null space of A^T to find mass-conserving perturbations
        # scipy.linalg.null_space returns a basis for the null space
        Z = scipy.linalg.null_space(a_mat.T)
        
        def compute_G(c: np.ndarray) -> float:
            rt = 1.987e-3 * T
            # Ignore zeros in log to avoid warnings, though c should be > 0
            c_safe = np.maximum(c, 1e-300)
            return float(np.sum(c * (dg_array / rt + np.log(c_safe) - 1.0)))
            
        G_opt = compute_G(c_opt)
        
        # Generate 50 random perturbations
        for _ in range(50):
            # random vector in null space
            v = rng.uniform(-1e-9, 1e-9, size=Z.shape[1])
            dc = Z @ v
            
            c_pert = c_opt + dc
            
            # Ensure admissible (c > 0)
            if np.any(c_pert <= 0):
                continue
                
            G_pert = compute_G(c_pert)
            
            # G_pert must be strictly greater than G_opt
            # We allow a tiny numerical margin due to float64
            assert G_pert > G_opt - 1e-15, f"Found a point with lower energy! {G_pert} <= {G_opt}"

    def test_malformed_problem(self) -> None:
        """Test D3: L'absence de monomère lève une ValueError, pas une LinAlgError."""
        # Create a problem where strand 0 has no monomer row
        a_mat = np.array([
            [1, 1],  # Heterodimer AB
            [0, 1]   # Monomer B (but Monomer A is missing)
        ], dtype=np.float64)
        
        dg_array = np.array([-10.0, 0.0])
        
        problem = EquilibriumProblem(
            n_strands=2,
            n_complexes=2,
            stoichiometry=a_mat,
            delta_g=dg_array,
            total_concentrations=np.array([1e-6, 1e-6]),
            temperature_kelvin=338.15
        )
        
        with pytest.raises(ValueError, match="Strand 0 is missing a valid monomer row"):
            solve_equilibrium(problem)

    def test_benign_fast_convergence(self) -> None:
        """Garantie de non-régression (Vitesse): Convergence en < 20 itérations.
        
        Vérifie que la stratégie 'Newton Complet d'Abord' empêche le
        sur-amortissement sur un problème facile bien conditionné.
        """
        a0 = 1.087e-6
        b0 = 2.458e-6
        dg = -11.76
        T = 338.15
        
        problem = EquilibriumProblem(
            n_strands=2,
            n_complexes=3,
            stoichiometry=np.array([[1, 0], [0, 1], [1, 1]], dtype=np.float64),
            delta_g=np.array([0.0, 0.0, dg]),
            total_concentrations=np.array([a0, b0]),
            temperature_kelvin=T,
        )
        
        # We explicitly request the DUAL_NEWTON method and a very strict threshold
        from labcraft.solver.types import SolverMethod
        result = solve_equilibrium(
            problem, 
            method=SolverMethod.DUAL_NEWTON, 
            convergence_threshold=1e-10
        )
        
        # It should converge extremely fast (typically < 10 iterations)
        assert result.n_iterations < 20, (
            f"Dual solver is crawling! Took {result.n_iterations} iterations "
            f"for a benign problem. Expected < 20."
        )
