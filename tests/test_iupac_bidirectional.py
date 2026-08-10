"""Tests de correspondance IUPAC bilatérale / Bidirectional IUPAC matching tests.

Vérifie que :
- Les codes d'ambiguïté vraie (R, Y, S, W, K, M, B, D, H, V) s'apparient
  par intersection d'ensembles, dans les deux sens.
- N (base inconnue) ne produit JAMAIS de correspondance.
"""

from labcraft.lamp.domains import _match_iupac_substring, _find_iupac_substring


# --- Codes d'ambiguïté vraie : doivent s'apparier ---

def test_true_ambiguity_target_has_code():
    """Cible contenant un code d'ambiguïté vraie, amorce a une base fixe."""
    assert _match_iupac_substring("A", "R") is True   # R = A ou G, A intersecte
    assert _match_iupac_substring("G", "R") is True   # G intersecte
    assert _match_iupac_substring("C", "R") is False  # C ne fait pas partie de R

def test_true_ambiguity_primer_has_code():
    """Amorce dégénérée sur cible propre."""
    assert _match_iupac_substring("R", "A") is True
    assert _match_iupac_substring("R", "G") is True
    assert _match_iupac_substring("R", "C") is False

def test_true_ambiguity_both_sides():
    """Les deux côtés portent des codes d'ambiguïté vraie."""
    assert _match_iupac_substring("R", "S") is True   # {A,G} ∩ {G,C} = {G}
    assert _match_iupac_substring("Y", "W") is True   # {C,T} ∩ {A,T} = {T}
    assert _match_iupac_substring("Y", "R") is False  # {C,T} ∩ {A,G} = {}

def test_multiposition_with_ambiguity():
    """Recherche multi-positions avec des codes d'ambiguïté."""
    assert _find_iupac_substring("TGC", "AATGC") == 2
    assert _find_iupac_substring("TRC", "AATGC") == 2   # R={A,G}, G matche
    assert _find_iupac_substring("TGC", "AATRC") == 2   # R côté cible, G matche


# --- N (base inconnue) : ne doit JAMAIS s'apparier ---

def test_n_in_target_no_match():
    """Une position N dans la cible ne produit pas de correspondance."""
    assert _match_iupac_substring("A", "N") is False
    assert _match_iupac_substring("T", "N") is False
    assert _match_iupac_substring("ATGC", "NNNNN") is False

def test_n_in_primer_no_match():
    """Une position N dans l'amorce ne produit pas de correspondance."""
    assert _match_iupac_substring("N", "A") is False
    assert _match_iupac_substring("N", "T") is False

def test_all_n_target_no_site():
    """Cible entièrement de N : aucun site trouvé, quelle que soit l'amorce."""
    target_all_n = "N" * 100
    primer = "GCCACCTTAAGCCACAGTA"  # 19 nt d'amorce quelconque
    assert _find_iupac_substring(primer, target_all_n) == -1

def test_partial_n_overlap_no_site():
    """Amorce chevauchant partiellement une plage de N : pas de site."""
    # Les 10 premières bases de l'amorce sont dans la vraie séquence,
    # les 9 suivantes tombent dans des N
    primer = "GCCACCTTAAGCCACAGTA"  # 19 nt
    target = "A" * 30 + primer[:10] + "N" * 9 + "A" * 30
    assert _find_iupac_substring(primer, target) == -1

def test_real_site_with_n_elsewhere():
    """Le vrai site est trouvé même si la cible contient des N ailleurs."""
    primer = "GCCACCTTAAGCCACAGTA"
    target = "N" * 20 + "AAAA" + primer + "AAAA" + "N" * 20
    idx = _find_iupac_substring(primer, target)
    assert idx == 24  # 20 N + 4 A


# --- Combinaisons : ambiguïté vraie ET N coexistent ---

def test_ambiguity_near_n():
    """Vrai code d'ambiguïté côté cible, avec des N ailleurs."""
    # Cible : ...R...N...  L'amorce doit matcher sur le R mais pas sur le N
    primer = "AG"
    target = "RNAG"  # pos 0: R=non (A et R matchent, mais pas G+R au pair)... 
    # On veut que AG matche à la position 2 (A,G)
    assert _find_iupac_substring(primer, target) == 2
    # Mais pas à la position 0 (R,N -> N bloque)
    target2 = "ANAG"
    assert _find_iupac_substring(primer, target2) == 2  # pas pos 0 car pos 1 = N


# --- Non-régression : comportement sur cible propre inchangé ---

def test_clean_target_unchanged():
    """Sur cible et amorce sans aucun code dégénéré, comportement exact standard."""
    assert _match_iupac_substring("ATGC", "ATGC") is True
    assert _match_iupac_substring("ATGC", "AAGC") is False
    assert _find_iupac_substring("ATGC", "CCATGCCC") == 2
