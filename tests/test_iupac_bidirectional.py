from labcraft.lamp.domains import _match_iupac_substring, _find_iupac_substring

def test_iupac_match():
    # 1. Cible contenant un code IUPAC là où l'amorce a une base fixe
    assert _match_iupac_substring("A", "R") is True  # R = A ou G
    assert _match_iupac_substring("G", "R") is True
    assert _match_iupac_substring("C", "R") is False
    
    # 2. Cible contenant un N
    assert _match_iupac_substring("A", "N") is True
    assert _match_iupac_substring("T", "N") is True
    
    # 3. Amorce dégénérée sur une cible propre
    assert _match_iupac_substring("R", "A") is True
    assert _match_iupac_substring("R", "G") is True
    assert _match_iupac_substring("R", "C") is False
    
    # 4. Multiples positions
    assert _match_iupac_substring("ATGC", "ATGC") is True
    assert _find_iupac_substring("TGC", "AATGC") == 2
    assert _find_iupac_substring("TRC", "AATGC") == 2
    assert _find_iupac_substring("TGC", "AATRC") == 2
