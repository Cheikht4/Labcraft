def test_tris_mm_case():
    from labcraft.buffer.monovalent import get_total_monovalent
    
    # Sans Tris
    mon_sans = get_total_monovalent(50.0, 0.0, 0.0)
    
    # Avec Tris
    mon_avec = get_total_monovalent(50.0, 0.0, 20.0)
    
    assert mon_avec > mon_sans
