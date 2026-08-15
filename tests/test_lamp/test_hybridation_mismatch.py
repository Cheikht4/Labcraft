import pytest
from unittest.mock import patch
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.lamp.stoichiometry import ConcentrationProfile
from labcraft.lamp.complex_enumeration import enumerate_complexes
from labcraft.thermo.backends.native import NativeBackend
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.thermo.salt import UnifiedSaltModel
class DummyBackend(NativeBackend):
    def calc_heterodimer(self, *args, **kwargs):
        from labcraft.thermo.backends.base import DuplexResult
        return DuplexResult(0, 0, 10.0, 0, '', 65.0)
    def calc_homodimer(self, *args, **kwargs):
        from labcraft.thermo.backends.base import DuplexResult
        return DuplexResult(0, 0, 10.0, 0, '', 65.0)
    def calc_hairpin(self, *args, **kwargs):
        from labcraft.thermo.backends.base import DuplexResult
        return DuplexResult(0, 0, 10.0, 0, '', 65.0)
from labcraft.solver.dual import solve_dual

def get_den3_primers():
    return [
        PhysicalPrimer('F3', 'GCCACCTTAAGCCACAGTA', PrimerRole.F3, 'GCCACCTTAAGCCACAGTA', '', None, 0.2e-6),
        PhysicalPrimer('B3', 'GTTGTGTCATGGGAGGG', PrimerRole.B3, 'GTTGTGTCATGGGAGGG', '', None, 0.2e-6),
        PhysicalPrimer('FIP', 'TGGCTTTTGGGCCTGACTTCTTTTTTGAAGAAGCTGTGCAGCCTG', PrimerRole.FIP, 'GAAGAAGCTGTGCAGCCTG', 'TGGCTTTTGGGCCTGACTTC', None, 1.6e-6),
        PhysicalPrimer('BIP', 'CTGTAGCTCCGTCGTGGGGATTTTCTAGTCTGCTACACCGTGC', PrimerRole.BIP, 'CTAGTCTGCTACACCGTGC', 'CTGTAGCTCCGTCGTGGGGA', None, 1.6e-6),
    ]

def _mock_calc_unfolding(window, local_start, local_end, temp_celsius, mon_molar):
    length = local_end - local_start
    if length == len('GCCACCTTAAGCCACAGTA'): return 4.354
    elif length == len('GTTGTGTCATGGGAGGG'): return 3.178
    elif length == len('GAAGAAGCTGTGCAGCCTG'): return 3.620
    elif length == len('CTAGTCTGCTACACCGTGC'): return 4.270
    return 0.0

@pytest.fixture
def target_seq():
    from labcraft.lamp.complex_enumeration import _revcomp
    f3 = 'GCCACCTTAAGCCACAGTA'
    f2 = 'GAAGAAGCTGTGCAGCCTG'
    b2 = 'CTAGTCTGCTACACCGTGC'
    b3 = 'GTTGTGTCATGGGAGGG'
    return f3 + 'A'*10 + f2 + 'A'*10 + _revcomp(b2) + 'A'*10 + _revcomp(b3)

def test_mismatch_internal_degrades_dg(target_seq):
    f2_mismatch = 'GAAGCAGCTGTGCAGCCTG'
    from labcraft.lamp.complex_enumeration import _revcomp
    target_mismatch = 'GCCACCTTAAGCCACAGTA' + 'A'*10 + f2_mismatch + 'A'*10 + _revcomp('CTAGTCTGCTACACCGTGC') + 'A'*10 + _revcomp('GTTGTGTCATGGGAGGG')
    primers = get_den3_primers()
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    profile = ConcentrationProfile(target=1.66e-16, fip_bip=1.6e-6, f3_b3=0.2e-6, lf_lb=0.8e-6)
    def _mock_find_iupac(query, target):
        if query == 'GAAGAAGCTGTGCAGCCTG': return target.find(f2_mismatch)
        from labcraft.lamp.domains import _find_iupac_substring as orig
        return orig(query, target)
        
    with patch('labcraft.lamp.complex_enumeration._find_iupac_substring', side_effect=_mock_find_iupac), patch('labcraft.lamp.complex_enumeration.calc_unfolding_penalty', side_effect=_mock_calc_unfolding):
        prob, strands, complexes, penalties = enumerate_complexes(primers, target_mismatch, backend, profile=profile, temp_celsius=63.0, mon_molar=0.05, buffer={'na_mM': 50.0, 'mg_mM': 8.0, 'dntp_mM': 1.4})
        
    assert penalties['FIP_site']['mismatches'] == 1
    assert penalties['FIP_site']['extensible'] == True
    fip_c_idx = complexes.index('FIP_on_FIP_site')
    assert prob.delta_g[fip_c_idx] > -13.5

def test_mismatch_3prime_veto(target_seq):
    f2_mismatch = 'GAAGAAGCTGTGCAGCCTC'
    from labcraft.lamp.complex_enumeration import _revcomp
    target_mismatch = 'GCCACCTTAAGCCACAGTA' + 'A'*10 + f2_mismatch + 'A'*10 + _revcomp('CTAGTCTGCTACACCGTGC') + 'A'*10 + _revcomp('GTTGTGTCATGGGAGGG')
    primers = get_den3_primers()
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    profile = ConcentrationProfile(target=1.66e-16, fip_bip=1.6e-6, f3_b3=0.2e-6, lf_lb=0.8e-6)
    def _mock_find_iupac(query, target):
        if query == 'GAAGAAGCTGTGCAGCCTG': return target.find(f2_mismatch)
        from labcraft.lamp.domains import _find_iupac_substring as orig
        return orig(query, target)
        
    with patch('labcraft.lamp.complex_enumeration._find_iupac_substring', side_effect=_mock_find_iupac), patch('labcraft.lamp.complex_enumeration.calc_unfolding_penalty', side_effect=_mock_calc_unfolding):
        prob, strands, complexes, penalties = enumerate_complexes(primers, target_mismatch, backend, profile=profile, temp_celsius=63.0, mon_molar=0.05, buffer={'na_mM': 50.0, 'mg_mM': 8.0, 'dntp_mM': 1.4})
        
    assert penalties['FIP_site']['mismatches'] == 1
    assert penalties['FIP_site']['extensible'] == False
    assert 'FIP_on_FIP_site' not in complexes
