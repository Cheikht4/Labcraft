import pytest
from labcraft.thermo.lna import parse_lna_sequence, load_lna_params, _LNA_PARAMS_MXL, _LNA_PARAMS_XLN
from labcraft.thermo.backends.native import NativeBackend

def test_lna_parsing():
    # Parsing basic
    bare, pos = parse_lna_sequence("AG+C+GTA")
    assert bare == "AGCGTA"
    assert pos == [2, 3]
    
    bare, pos = parse_lna_sequence("+A")
    assert bare == "A"
    assert pos == [0]
    
    bare, pos = parse_lna_sequence("ACGT")
    assert bare == "ACGT"
    assert pos == []
    
    with pytest.raises(ValueError, match="Dangling"):
        parse_lna_sequence("ACGT+")

def test_lna_table_loaded():
    load_lna_params()
    assert len(_LNA_PARAMS_MXL) == 16
    assert len(_LNA_PARAMS_XLN) == 16
    
    # Check some known values from the prompt
    # XLN(G,C) = (-0.360, -0.251, -0.414)
    # Wait, XLN table has key (base_lna, base3) -> ("G", "C")
    dh, ds = _LNA_PARAMS_XLN[("G", "C")]
    assert abs(dh - -0.360) < 0.001
    assert abs(ds - -0.251) < 0.001
    
    # MXL(G,C) = (-0.925, -1.111, -0.535)
    dh, ds = _LNA_PARAMS_MXL[("G", "C")]
    assert abs(dh - -0.925) < 0.001
    assert abs(ds - -1.111) < 0.001

def test_golden_anchor():
    """
    Test d'ancrage golden de l'algorithme : pour une LNA-C INTERNE avec voisin 5' = G et
    voisin 3' = T, la somme des deux corrections (MXL(G,C) + XLN(C,T)) vaut
    ΔΔH = -0,217 kcal/mol, ΔΔS = 3,064 cal/mol·K, ΔΔG37 = -1,201 kcal/mol.
    """
    load_lna_params()
    # Sequence with G +C T.
    seq = "AG+CTA"
    backend = NativeBackend()
    
    bare, pos = parse_lna_sequence(seq)
    # Get native without LNA
    res_native = backend._calc_perfect_duplex(bare, 37.0, 1e-6)
    # Get with LNA
    res_lna = backend._calc_perfect_duplex(bare, 37.0, 1e-6, lna_positions=tuple(pos))
    
    ddH = res_lna.dh_kcal - res_native.dh_kcal
    ddS = res_lna.ds_cal_per_k - res_native.ds_cal_per_k
    
    # Note: dg_kcal is at 37°C
    # Let's just check the ddH and ddS
    # Expected: ddH = -0.217, ddS = 3.064
    assert abs(ddH - -0.217) < 0.01, f"Expected -0.217, got {ddH}"
    assert abs(ddS - 3.064) < 0.01, f"Expected 3.064, got {ddS}"
    
def test_global_effect():
    """un duplexe avec une ou plusieurs LNA a un Tm PLUS élevé et un ΔG
    plus négatif que le même duplexe tout-ADN."""
    backend = NativeBackend()
    seq_dna = "CCTTGGACGGG"
    seq_lna = "CCTTGG+ACGGG"
    
    def _revcomp(seq: str) -> str:
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return "".join(complement[c] for c in reversed(seq))
        
    res_dna = backend.calc_duplex(seq_dna, _revcomp(seq_dna), temp_celsius=65.0)
    bare, pos = parse_lna_sequence(seq_lna)
    res_lna = backend.calc_duplex(bare, _revcomp(bare), temp_celsius=65.0, lna_positions=tuple(pos))
    
    assert res_lna.tm_celsius > res_dna.tm_celsius
    assert res_lna.dg_kcal < res_dna.dg_kcal
