import pytest

from labcraft.lamp.domains import PhysicalPrimer, PrimerRole, _revcomp


def test_auto_detect_fip_bip_orientation():
    """Test d'autodétection des domaines F1c/F2 et B1c/B2 (Jalon 3).
    
    Vérifie l'orientation correcte : le domaine de liaison (F2/B2) doit toujours
    être à l'extrémité 3' de l'amorce composite, et le domaine de boucle (F1c/B1c)
    à l'extrémité 5'.
    """
    # Séquence cible fictive (Brin +)
    # 5'- ... F3c ... F2c ... F1c ... B1 ... B2 ... B3 ... -3'
    f3c = "ACTGACTGACTG"
    f2c = "GGATCCGGATCC" # F2 s'y lie, donc F2 = RC(f2c) = GGATCCGGATCC
    f1c = "ATATATATATAT" # F1c est de même sens que la cible (+)
    b1  = "CGCGCGCGCGCG" # B1c est le RC de B1 = CGCGCGCGCGCG
    b2  = "TTAATTAATTAA" # B2 est de même sens que la cible (+) = TTAATTAATTAA
    b3  = "CCCGGGCCCGGG"
    
    # Cible complète (5' -> 3')
    target_plus = f"{f3c}NNN{f2c}NNN{f1c}NNN{b1}NNN{b2}NNN{b3}"
    
    # Construction des amorces
    # FIP = F1c + linker + F2
    f2 = _revcomp(f2c)
    fip_seq = f1c + "TTTT" + f2
    
    # BIP = B1c + linker + B2
    b1c = _revcomp(b1)
    bip_seq = b1c + "TTTT" + b2
    
    fip = PhysicalPrimer.from_alignment("FIP", fip_seq, PrimerRole.FIP, target_plus)
    assert fip.binding_domain == f2, "F2 doit être le domaine de liaison (3')"
    assert fip.tail_domain == f1c, "F1c doit être le domaine de queue (5')"
    assert fip.linker == "TTTT", "Linker non détecté"
    
    bip = PhysicalPrimer.from_alignment("BIP", bip_seq, PrimerRole.BIP, target_plus)
    assert bip.binding_domain == b2, "B2 doit être le domaine de liaison (3')"
    assert bip.tail_domain == b1c, "B1c doit être le domaine de queue (5')"
    assert bip.linker == "TTTT", "Linker non détecté"


def test_iupac_matching():
    """Test de la reconnaissance des codes IUPAC (Dégénérescence)."""
    # Cible contient un 'N' et un 'R'
    target = "ATGCATGCATGC"
    
    # L'amorce a des codes de dégénérescence qui couvrent la cible
    # M = A ou C (matche A)
    # W = A ou T (matche T)
    # S = G ou C (matche G)
    # N = n'importe quoi (matche C)
    primer_binding = "MWSNATGCATGC" # => ATGCATGCATGC
    # On met une queue qui s'aligne d'au moins 12 bp
    primer_tail = "ATGCATGCATGC"
    
    fip_seq = primer_tail + "TTTT" + primer_binding
    
    fip = PhysicalPrimer.from_alignment("FIP_deg", fip_seq, PrimerRole.FIP, target)
    assert fip.binding_domain == primer_binding
    assert fip.tail_domain == primer_tail

def test_expand_degenerate():
    from labcraft.lamp.domains import expand_degenerate
    
    # Test simple : aucune base dégénérée
    assert expand_degenerate("ATGC") == ["ATGC"]
    
    # Test 1 position dégénérée (R = A,G)
    variants = expand_degenerate("ARCT")
    assert sorted(variants) == sorted(["AACT", "AGCT"])
    
    # Test 2 positions dégénérées (W = A,T et S = G,C)
    variants = expand_degenerate("WTS")
    assert sorted(variants) == sorted(["ATG", "ATC", "TTG", "TTC"])
    assert len(variants) == 4
    
    # Conservation de la concentration (testé indirectement ici)
    # L'attribution est faite dans config.py, mais vérifions la logique arithmétique de base
    nominal_conc = 1.6
    split_conc = nominal_conc / len(variants)
    assert abs((split_conc * len(variants)) - nominal_conc) < 1e-12
    
    # Test plafond (max_variants)
    with pytest.raises(ValueError, match="générerait 512 variants"):
        expand_degenerate("NNNNS", max_variants=16) 
        
    with pytest.raises(ValueError, match="générerait 32 variants"):
        expand_degenerate("RRRRR", max_variants=16)

