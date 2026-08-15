import pytest
import os
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.lamp.stoichiometry import ConcentrationProfile
from labcraft.lamp.complex_enumeration import enumerate_complexes
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.thermo.salt import UnifiedSaltModel
from labcraft.solver.dual import solve_dual

def get_den3_primers():
    return [
        PhysicalPrimer('F3', 'GCCACCTTAAGCCACAGTA', PrimerRole.F3, 'GCCACCTTAAGCCACAGTA', '', None, 0.2e-6),
        PhysicalPrimer('B3', 'GTTGTGTCATGGGAGGG', PrimerRole.B3, 'GTTGTGTCATGGGAGGG', '', None, 0.2e-6),
        PhysicalPrimer('FIP', 'TGGCTTTTGGGCCTGACTTCTTTTTTGAAGAAGCTGTGCAGCCTG', PrimerRole.FIP, 'GAAGAAGCTGTGCAGCCTG', 'TGGCTTTTGGGCCTGACTTC', None, 1.6e-6),
        PhysicalPrimer('BIP', 'CTGTAGCTCCGTCGTGGGGATTTTCTAGTCTGCTACACCGTGC', PrimerRole.BIP, 'CTAGTCTGCTACACCGTGC', 'CTGTAGCTCCGTCGTGGGGA', None, 1.6e-6),
        PhysicalPrimer('LF', 'CCTTGGACGGGGCT', PrimerRole.LF, 'CCTTGGACGGGGCT', '', None, 0.8e-6),
        PhysicalPrimer('LB', 'GGAGGCTGCAAACCGTG', PrimerRole.LB, 'GGAGGCTGCAAACCGTG', '', None, 0.8e-6),
    ]

def _mock_calc_unfolding(window, local_start, local_end, temp_celsius, mon_molar=None, buffer=None):
    if window[local_start:local_end].upper() == 'GCCACCTTAAGCCACAGTA': return 4.35
    elif window[local_start:local_end].upper() == 'CCCTCCCATGACACAAC': return 3.18
    elif window[local_start:local_end].upper() == 'GAAGAAGCTGTGCAGCCTG': return 3.62
    elif window[local_start:local_end].upper() == 'GCACGGTGTAGCAGACTAG': return 4.27
    return 0.0

def test_non_regression_den3_full_genome():
    fasta_path = os.path.join(os.path.dirname(__file__), '../../DEN3_M93130.fasta')
    with open(fasta_path, 'r') as f:
        target_seq = f.read().split('\n', 1)[1].replace('\n', '').upper()
        
    primers = get_den3_primers()
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    profile = ConcentrationProfile(target=1.66e-16, fip_bip=1.6e-6, f3_b3=0.2e-6, lf_lb=0.8e-6)

    from unittest.mock import patch
    with patch('labcraft.lamp.complex_enumeration.calc_unfolding_penalty', side_effect=_mock_calc_unfolding):
        prob, strands, complexes, penalties = enumerate_complexes(
            primers, target_seq, backend, profile=profile, temp_celsius=63.0, mon_molar=0.05, buffer={'na_mM': 50.0, 'mg_mM': 8.0, 'dntp_mM': 1.4}
        )
    res = solve_dual(prob)

    def get_occ(name):
        c_name = name.replace('_site', '') + '_on_' + name
        if c_name not in complexes: return 0.0
        c_idx = complexes.index(c_name)
        t_idx = strands.index(name)
        return res.concentrations[c_idx] / prob.total_concentrations[t_idx]

    f3_occ = get_occ('F3_site')
    b3_occ = get_occ('B3_site')
    fip_occ = get_occ('FIP_site')
    bip_occ = get_occ('BIP_site')
    
    assert abs(penalties['F3_site']['dg_unfold'] - 4.35) < 0.05
    assert abs(penalties['B3_site']['dg_unfold'] - 3.18) < 0.05
    assert abs(penalties['FIP_site']['dg_unfold'] - 3.62) < 0.05
    assert abs(penalties['BIP_site']['dg_unfold'] - 4.27) < 0.05

    assert abs(f3_occ - 0.026) < 0.05, f"F3 expected 2.6%, got {f3_occ*100}%"
    assert abs(b3_occ - 0.013) < 0.05, f"B3 expected 1.3%, got {b3_occ*100}%"
    assert abs(fip_occ - 0.725) < 0.05, f"FIP expected 72.5%, got {fip_occ*100}%"
    assert abs(bip_occ - 0.239) < 0.05, f"BIP expected 23.9%, got {bip_occ*100}%"
