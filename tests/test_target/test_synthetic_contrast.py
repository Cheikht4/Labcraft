import pytest
from labcraft.target.unfolding import calc_unfolding_penalty

def test_synthetic_target_contrast():
    """Vérifie que SynthA présente un contraste d'accessibilité net.
    
    Le site F3_A doit être fermement bloqué dans une tige-boucle (> 5 kcal/mol),
    tandis que les régions flanquantes (contrôles) doivent rester très 
    accessibles (< 1 kcal/mol).
    """
    # Séquence de SynthA telle que générée et optimisée
    seq = "CGCCCGAAATTTCAACAAAAGGCATGCCTAGCTAGCAATCGTACGCATAGTAGGCATGCCAAAAAAGGAGTTAGCGAACG"
    
    # 1. Contrôle 5' (20 nt)
    dg_c1 = calc_unfolding_penalty(seq, 0, 20, temp_celsius=65.0)
    
    # 2. Site F3_A (20 nt) : index 30 à 50
    dg_f3 = calc_unfolding_penalty(seq, 30, 50, temp_celsius=65.0)
    
    # 3. Contrôle 3' (20 nt) : index 60 à 80
    dg_c2 = calc_unfolding_penalty(seq, 60, 80, temp_celsius=65.0)
    
    # Assertions selon les spécifications
    assert dg_f3 > 5.0, f"Le site F3_A n'est pas assez bloqué (ΔG = {dg_f3:.2f} kcal/mol)"
    assert dg_c1 < 1.0, f"Le contrôle 5' est trop structuré (ΔG = {dg_c1:.2f} kcal/mol)"
    assert dg_c2 < 1.0, f"Le contrôle 3' est trop structuré (ΔG = {dg_c2:.2f} kcal/mol)"
