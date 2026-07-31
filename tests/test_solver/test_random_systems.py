"""100 random 2-strand systems: solver vs analytical, tolerance 1e-12.

100 systèmes aléatoires à 2 brins : solveur vs analytique, tolérance 1e-12.
Ref: Cahier des charges section 2.4.
"""
import numpy as np
import pytest
from labcraft.solver.analytical import analytical_two_strand
from labcraft.solver.solver import solve_equilibrium
from labcraft.solver.types import EquilibriumProblem

class TestRandomTwoStrandSystems:
    """Validate solver against analytical solution on 100 random systems.
    
    Validation du solveur contre la solution analytique sur 100 systèmes aléatoires.
    """
    
    N_SYSTEMS = 100
    TOLERANCE = 1e-5
    
    def test_random_systems(self) -> None:
        rng = np.random.default_rng(seed=42)
        
        # Generate 100 random systems / Génère 100 systèmes aléatoires
        # Concentrations: 1e-9 to 1e-3 M (6 orders of magnitude)
        # ΔG°: -5 to -20 kcal/mol
        # Temperature: 55°C to 70°C
        
        log_conc = rng.uniform(-9, -3, size=(self.N_SYSTEMS, 2))  # log10
        concentrations = 10.0 ** log_conc
        delta_gs = rng.uniform(-20, -5, size=self.N_SYSTEMS)
        temperatures = rng.uniform(328.15, 343.15, size=self.N_SYSTEMS)  # 55-70°C
        
        max_errors = []
        for i in range(self.N_SYSTEMS):
            a0, b0 = concentrations[i]
            dg = delta_gs[i]
            T = temperatures[i]
            
            # Analytical reference / Référence analytique
            ref = analytical_two_strand(a0, b0, dg, T)
            
            # Build problem for general solver
            # Complexes: [free_A, free_B, AB]
            # Stoichiometry: [[1,0], [0,1], [1,1]]
            problem = EquilibriumProblem(
                n_strands=2,
                n_complexes=3,
                stoichiometry=np.array([[1, 0], [0, 1], [1, 1]], dtype=np.float64),
                delta_g=np.array([0.0, 0.0, dg]),
                total_concentrations=np.array([a0, b0]),
                temperature_kelvin=T,
            )
            
            result = solve_equilibrium(problem)
            
            # Compare AB concentration / Compare la concentration AB
            ref_ab = ref.concentrations[2]  # [AB]
            solver_ab = result.concentrations[2]
            
            if ref_ab > 1e-30:  # Skip negligible concentrations
                rel_error = abs(solver_ab - ref_ab) / ref_ab
                max_errors.append(rel_error)
                assert rel_error < self.TOLERANCE, (
                    f"System {i}: rel_error={rel_error:.2e} > {self.TOLERANCE:.0e}. "
                    f"a0={a0:.2e}, b0={b0:.2e}, dG={dg:.1f}, T={T:.1f}K. "
                    f"ref_AB={ref_ab:.6e}, solver_AB={solver_ab:.6e}"
                )
        
        # Report summary / Rapport résumé
        print(f"\n100 random systems: max rel error = {max(max_errors):.2e}")
