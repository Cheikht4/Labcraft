import pytest
from labcraft.diagnostics.mispriming import detect_inter_target_mispriming, _revcomp
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.thermo.salt import UnifiedSaltModel
from labcraft.diagnostics.enzyme import ENZYME_REGISTRY

@pytest.fixture
def backend():
    return ViennaSaltShiftBackend(UnifiedSaltModel())

@pytest.fixture
def enzyme():
    return ENZYME_REGISTRY["bst2.0"]

def test_mispriming_perfect_match(backend, enzyme):
    # L'amorce s'apparie parfaitement sur la cible YFV
    p = PhysicalPrimer.from_simple("F3_X", "ATCGATCGATCGATCG", PrimerRole.F3)
    target_seq = "NNNN" + _revcomp("ATCGATCGATCGATCG") + "NNNN"
    
    primer_to_panel = {"F3_X": "X"}
    targets = {"YFV": target_seq}
    
    risks = detect_inter_target_mispriming(
        [p], primer_to_panel, targets, backend, enzyme, temp_celsius=65.0
    )
    
    assert len(risks) > 0
    assert risks[0].primer_name == "F3_X"
    assert risks[0].target_id == "YFV"

def test_mispriming_internal_mismatch(backend, enzyme):
    # L'extrémité 3' (disons les 3 derniers nt) est parfaite.
    # On introduit un mismatch à la position 5 (en partant du 3').
    seq = "ATCGATCGATCGATCG"
    p = PhysicalPrimer.from_simple("F3_X", seq, PrimerRole.F3)
    # revcomp: CGATCGATCGATCGAT
    # 3' end de p = TCG. Son revcomp = CGA.
    # Position 5 du 3' = C (index -5 de seq, seq[-5] = G, revcomp = C).
    # On change ce C en A dans la cible pour faire un mismatch
    rc_seq = _revcomp(seq)
    # L'index de seq[-5] dans rc_seq est 4 (0-indexed).
    target_rc = rc_seq[:4] + "A" + rc_seq[5:]
    
    target_seq = "NNNN" + target_rc + "NNNN"
    
    primer_to_panel = {"F3_X": "X"}
    targets = {"YFV": target_seq}
    
    risks = detect_inter_target_mispriming(
        [p], primer_to_panel, targets, backend, enzyme, temp_celsius=65.0
    )
    
    # Doit être détecté car l'ancrage est sur les 3 derniers, et le mismatch n'est pas dans le veto 3' terminal (pos 1-2).
    assert len(risks) > 0

def test_mispriming_terminal_mismatch(backend, enzyme):
    # L'extrémité 3' (position 1 terminale) a un mismatch.
    # L'ancrage (K_LEN=3) sur l'extrémité 3' EXACTE ne trouvera pas la cible.
    # Ce cas n'est PAS remonté (veto/ancrage).
    seq = "ATCGATCGATCGATCG"
    p = PhysicalPrimer.from_simple("F3_X", seq, PrimerRole.F3)
    
    # 3' terminal est G. Dans le revcomp, c'est le 1er nt (C).
    rc_seq = _revcomp(seq)
    target_rc = "A" + rc_seq[1:]
    
    target_seq = "NNNN" + target_rc + "NNNN"
    
    primer_to_panel = {"F3_X": "X"}
    targets = {"YFV": target_seq}
    
    risks = detect_inter_target_mispriming(
        [p], primer_to_panel, targets, backend, enzyme, temp_celsius=65.0
    )
    
    # Ne doit pas être détecté (veto / pas ancré)
    assert len(risks) == 0

def test_mispriming_no_homology(backend, enzyme):
    p = PhysicalPrimer.from_simple("F3_X", "ATCGATCGATCGATCG", PrimerRole.F3)
    target_seq = "AAAAAAAAAAAAAAAAAAAAAAA"
    
    primer_to_panel = {"F3_X": "X"}
    targets = {"YFV": target_seq}
    
    risks = detect_inter_target_mispriming(
        [p], primer_to_panel, targets, backend, enzyme, temp_celsius=65.0
    )
    
    assert len(risks) == 0
