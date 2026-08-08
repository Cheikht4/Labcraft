import pytest
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.thermo.backends.vienna import ViennaRNABackend
from labcraft.diagnostics.enzyme import BST
from labcraft.diagnostics.amplifiable_dimer import evaluate_pair_amplifiable

def test_blocked_3prime_veto():
    backend = ViennaRNABackend()
    
    # Séquences d'un dimère très fort qui est amplifiable par defaut
    # On va prendre un exemple artificiel qui s'apparie parfaitement sur l'extrémité 3'
    seq_a = "ATCGATCGATCGATCG"
    seq_b = "CGATCGATCGATCGAT"  # reverse complement of seq_a is CGATCGATCGATCGAT
    
    p_a = PhysicalPrimer.from_simple("A", seq_a, PrimerRole.FIP)
    p_b = PhysicalPrimer.from_simple("B", seq_b, PrimerRole.BIP)
    
    # Sans blocage
    is_amp, dg_3p, details = evaluate_pair_amplifiable(p_a, p_b, backend, BST, 65.0)
    assert is_amp == True
    assert details["is_blocked_veto"] == False
    
    # Avec blocage sur A et B
    p_a.blocked_3prime = True
    p_b.blocked_3prime = True
    is_amp, dg_3p, details = evaluate_pair_amplifiable(p_a, p_b, backend, BST, 65.0)
    assert is_amp == False
    assert details["is_blocked_veto"] == True

def test_probe_tm_check():
    from labcraft.diagnostics.probe_tm import check_probes_tm
    backend = ViennaRNABackend()
    
    p_a = PhysicalPrimer.from_simple("FIP", "ATCGATCGATCGATCGATCG", PrimerRole.FIP, nominal_concentration=0.8e-6)
    p_b = PhysicalPrimer.from_simple("BIP", "ATCGATCGATCGATCGATCG", PrimerRole.BIP, nominal_concentration=0.8e-6)
    
    # Tm très haut (fort %GC et longueur)
    probe_seq = "GCGCGCGCGCGCGCGCGCGCGCGC"
    probe = PhysicalPrimer.from_simple("PROBE1", probe_seq, PrimerRole.PROBE, nominal_concentration=0.2e-6)
    
    res = check_probes_tm([p_a, p_b, probe], backend, 65.0)
    assert len(res) == 1
    assert res[0]["probe_name"] == "PROBE1"
    assert res[0]["is_ok"] == True # Probe is very GC rich, Tm should be much higher than AT-rich primers
    
    # Sonde faible
    probe_seq2 = "ATATATATATATATATATAT"
    probe2 = PhysicalPrimer.from_simple("PROBE2", probe_seq2, PrimerRole.PROBE, nominal_concentration=0.2e-6)
    res2 = check_probes_tm([p_a, p_b, probe2], backend, 65.0)
    assert res2[0]["is_ok"] == False
