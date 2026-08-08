import pytest
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.thermo.salt import UnifiedSaltModel

def test_vienna_salt_shift_reference_conditions():
    """Conditions de référence (mon 1 M, Mg 0) : ΔG_final == ΔG_ref à 1e-6 près."""
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    # FIP_DEN-3 vs FLP_DEN-3 approx
    seq1 = "GCTGCGTTGTGTCTTGGGAGGTTTTCTGTACGCATGGGGTAGC"
    seq2 = "CTCCTCTAACCACTAGTC"
    
    # Conditions de ref (na_mm = 1000)
    res_ref = backend.calc_heterodimer(seq1, seq2, temp_celsius=65.0, na_mm=1000.0, mg_mm=0.0)
    
    # The native ViennaRNA backend gives a specific dG at 65°C.
    # Since we're at reference conditions, the ddg_salt should be 0.
    # So dG should match exactly what ViennaRNA gives natively.
    from labcraft.thermo.backends.vienna import ViennaRNABackend
    vienna = ViennaRNABackend()
    res_vienna = vienna.calc_heterodimer(seq1, seq2, temp_celsius=65.0)
    
    assert abs(res_ref.dg_kcal - res_vienna.dg_kcal) < 1e-6

def test_vienna_salt_shift_stabilizes_with_mg():
    """Sur un duplexe parfait, augmenter le Mg rend le ΔG PLUS négatif et le Tm PLUS haut."""
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    
    seq1 = "ATGCATGCATGCATGCATGC"
    seq2 = "GCATGCATGCATGCATGCAT"
    
    res_no_mg = backend.calc_heterodimer(seq1, seq2, temp_celsius=65.0, na_mm=50.0, mg_mm=0.0)
    res_mg = backend.calc_heterodimer(seq1, seq2, temp_celsius=65.0, na_mm=50.0, mg_mm=6.0)
    
    assert res_mg.dg_kcal < res_no_mg.dg_kcal
    assert res_mg.tm_celsius > res_no_mg.tm_celsius

def test_vienna_salt_shift_non_zero_dimer():
    """Un dimère connu sous tampon donne un ΔG NON nul et fini."""
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    
    # FLP et BIP from the user prompt
    seq1 = "CTCCTCTAACCACTAGTC"
    seq2 = "CCCAACACCAGGGGAAGCTGTTTTTTTGTTGTTGTGCGGGGG"
    
    res = backend.calc_heterodimer(seq1, seq2, temp_celsius=65.0, na_mm=50.0, mg_mm=6.0)
    
    assert res.dg_kcal != 0.0
    assert not (abs(res.dg_kcal) < 1e-3)
    assert res.dg_kcal < -2.0 # It should be around -6.86 or so

def test_vienna_salt_shift_fgc_nbp():
    """Vérifie que f_gc et n_bp sont corrects sur un duplexe parfait."""
    from labcraft.thermo.backends.vienna_salt import estimate_helix_thermo
    # Séquence de 20 bases avec 10 GC -> f_gc = 0.5
    seq = "ATGCATGCATGCATGCATGC&GCATGCATGCATGCATGCAT"
    structure = "((((((((((((((((((((&))))))))))))))))))))"
    structure = structure.replace('&', '')
    dh, ds, n_bp, f_gc = estimate_helix_thermo(seq, structure)
    assert n_bp == 20
    assert abs(f_gc - 0.5) < 0.01

def test_vienna_salt_shift_amplitude():
    """L'amplitude du décalage salin doit coïncider avec SaltCorrectedBackend(NativeBackend) à 0.3 kcal/mol près."""
    from labcraft.thermo.salt import UnifiedSaltModel, SaltCorrectedBackend
    from labcraft.thermo.backends.native import NativeBackend
    from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
    
    native_salt = SaltCorrectedBackend(NativeBackend(), UnifiedSaltModel())
    vienna_salt = ViennaSaltShiftBackend(UnifiedSaltModel())
    
    seq1 = "ATGCATGCATGCATGCATGC"
    seq2 = "GCATGCATGCATGCATGCAT"
    
    res_native_mg0 = native_salt.calc_heterodimer(seq1, seq2, temp_celsius=65.0, na_mm=50.0, mg_mm=0.0)
    res_native_mg6 = native_salt.calc_heterodimer(seq1, seq2, temp_celsius=65.0, na_mm=50.0, mg_mm=6.0)
    native_shift = res_native_mg6.dg_kcal - res_native_mg0.dg_kcal
    
    res_vienna_mg0 = vienna_salt.calc_heterodimer(seq1, seq2, temp_celsius=65.0, na_mm=50.0, mg_mm=0.0)
    res_vienna_mg6 = vienna_salt.calc_heterodimer(seq1, seq2, temp_celsius=65.0, na_mm=50.0, mg_mm=6.0)
    vienna_shift = res_vienna_mg6.dg_kcal - res_vienna_mg0.dg_kcal
    
    assert abs(native_shift - vienna_shift) < 0.3
