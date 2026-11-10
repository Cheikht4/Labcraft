# Note : Ce test à -6.27 kcal/mol est un test de non-régression numérique d'un dimère
# (conditions standards 1 M Na+, sans correction saline) et ne valide pas la biologie.
# Note: This test at -6.27 kcal/mol is a numerical non-regression test of a dimer
# (standard 1 M Na+ conditions, without salt correction) and does not validate biology.
def test_evaluate_pair_amplifiable_order():
    from labcraft.diagnostics.amplifiable_dimer import evaluate_pair_amplifiable
    from labcraft.thermo.backends.vienna import ViennaRNABackend
    from labcraft.diagnostics.enzyme import BST
    
    # Parida DEN-3: FLP et BIP
    flp = "CCTTGGACGGGGCT"
    bip = "CTGTAGCTCCGTCGTGGGGATTTTCTAGTCTGCTACACCGTGC"
    
    backend = ViennaRNABackend()
    from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
    p_flp = PhysicalPrimer.from_simple("FLP", flp, PrimerRole.FIP)
    p_bip = PhysicalPrimer.from_simple("BIP", bip, PrimerRole.BIP)
    
    is_amp, dg, details = evaluate_pair_amplifiable(p_flp, p_bip, backend, BST, 63.0)
    
    assert is_amp is True
    assert abs(dg - -6.27) < 0.05
    assert details["order"] in ["a&b", "b&a"]
    
    # Reverse order should give the exact same output
    is_amp2, dg2, details2 = evaluate_pair_amplifiable(p_bip, p_flp, backend, BST, 63.0)
    assert is_amp2 is True
    assert abs(dg2 - -6.27) < 0.05
