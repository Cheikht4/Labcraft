def test_evaluate_pair_amplifiable_order():
    from labcraft.diagnostics.amplifiable_dimer import evaluate_pair_amplifiable
    from labcraft.thermo.backends.vienna import ViennaRNABackend
    from labcraft.diagnostics.enzyme import BST
    
    # Parida DEN-3: FLP et BIP
    flp = "CCTTGGACGGGGCT"
    bip = "CTGTAGCTCCGTCGTGGGGATTTTCTAGTCTGCTACACCGTGC"
    
    backend = ViennaRNABackend()
    is_amp, dg, details = evaluate_pair_amplifiable(flp, bip, backend, BST, 63.0)
    
    assert is_amp is True
    assert abs(dg - -6.27) < 0.05
    assert details["order"] in ["a&b", "b&a"]
    
    # Reverse order should give the exact same output
    is_amp2, dg2, details2 = evaluate_pair_amplifiable(bip, flp, backend, BST, 63.0)
    assert is_amp2 is True
    assert abs(dg2 - -6.27) < 0.05
