import yaml
import numpy as np
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.lamp.stoichiometry import ConcentrationProfile
from labcraft.thermo.backends.vienna import ViennaRNABackend
from labcraft.lamp.complex_enumeration import enumerate_complexes
from labcraft.solver.dual import solve_dual

primers = [
    PhysicalPrimer("F3", "GATCGATCGGCCTTAAAAA", PrimerRole.F3, "GATCGATCGGCCTTAAAAA"),
    PhysicalPrimer("B3", "TTTTTTTTAAGGCC", PrimerRole.B3, "TTTTTTTTAAGGCC")
]
backend = ViennaRNABackend()
profile = ConcentrationProfile(target=1e-12, fip_bip=1.6e-6, f3_b3=0.2e-6, lf_lb=0.8e-6)
prob, strands, complexes, _ = enumerate_complexes(primers, "A", backend, profile=profile, temp_celsius=65.0)
res = solve_dual(prob)
R = 0.0019872
RT = R * (273.15 + 65.0)
u = np.log(res.free_concentrations)

for i, c_name in enumerate(complexes):
    stoich = prob.stoichiometry[i]
    conc = np.exp(-prob.delta_g[i] / RT + np.dot(stoich, u))
    if "homo" not in c_name and "free" not in c_name and "on" not in c_name:
        print(f"{c_name}: conc={conc:.2e}, dG={prob.delta_g[i]:.2f}")
