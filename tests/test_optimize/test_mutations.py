import pytest
from labcraft.optimize.mutations import enumerate_variants

def test_enumerate_variants_counts():
    seq = "CCTTGGACGGGGCT"
    variants = enumerate_variants(seq, max_mutations=2, window_3p=6)
    
    # original = 1
    # 1 mut on 6 bases = 6 * 3 = 18
    # 2 mut on 6 bases = (6 * 5 / 2) * 3 * 3 = 15 * 9 = 135
    # total = 1 + 18 + 135 = 154
    assert len(variants) == 154
    
    # Check that original is the first
    assert variants[0][0] == seq
    assert len(variants[0][1]) == 0
    
    # Check a 1-mut variant
    mut1 = [v for v in variants if len(v[1]) == 1]
    assert len(mut1) == 18
    
    # Check a 2-mut variant
    mut2 = [v for v in variants if len(v[1]) == 2]
    assert len(mut2) == 135
