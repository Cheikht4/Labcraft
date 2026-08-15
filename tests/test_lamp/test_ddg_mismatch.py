import pytest
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.thermo.salt import UnifiedSaltModel

def test_ddg_mismatch_zero_for_perfect():
    from labcraft.thermo.mismatch import nn_duplex_energy
    primers = [
        PhysicalPrimer('F3_DEN3', 'GCCACCTTAAGCCACAGTA', PrimerRole.F3, 'GCCACCTTAAGCCACAGTA'),
        PhysicalPrimer('B3_DEN3', 'GTTGTGTCATGGGAGGG', PrimerRole.B3, 'GTTGTGTCATGGGAGGG'),
        PhysicalPrimer('FIP_DEN3', 'TGGCTTTTGGGCCTGACTTCTTTTTTGAAGAAGCTGTGCAGCCTG', PrimerRole.FIP, 'GAAGAAGCTGTGCAGCCTG', 'TGGCTTTTGGGCCTGACTTC'),
        PhysicalPrimer('BIP_DEN3', 'CTGTAGCTCCGTCGTGGGGATTTTCTAGTCTGCTACACCGTGC', PrimerRole.BIP, 'CTAGTCTGCTACACCGTGC', 'CTGTAGCTCCGTCGTGGGGA'),
        PhysicalPrimer('LF_DEN3', 'CCTTGGACGGGGCT', PrimerRole.LF, 'CCTTGGACGGGGCT'),
        PhysicalPrimer('LB_DEN3', 'GGAGGCTGCAAACCGTG', PrimerRole.LB, 'GGAGGCTGCAAACCGTG'),
    ]
    temp_celsius = 63.0
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    for p in primers:
        perfect_bottom = "".join(comp.get(c, c) for c in p.binding_domain)
        _, _, dg_perfect_nn = nn_duplex_energy(p.binding_domain, perfect_bottom, temp_celsius)
        _, _, dg_mismatched_nn = nn_duplex_energy(p.binding_domain, perfect_bottom, temp_celsius)
        
        ddg_mismatch = dg_mismatched_nn - dg_perfect_nn
        assert abs(ddg_mismatch) < 1e-6, f"ddg_mismatch is not zero for perfect match {p.name}: {ddg_mismatch}"
