import math
import pathlib

import numpy as np
import pandas as pd
import pytest

from labcraft.thermo.backends.native import NativeBackend
from labcraft.thermo.salt import UnifiedSaltModel, SaltCorrectedBackend
from labcraft.buffer.monovalent import get_total_monovalent


def test_owczarzy2008_golden_anchor():
    """Reproduit strictement l'ancre dorée d'Owczarzy 2008."""
    salt_model = UnifiedSaltModel()
    
    # Entrées
    f_gc = 0.6
    n_bp = 20
    tm_1m_celsius = 76.3
    mg_molar = 0.0015
    mon_molar = 0.0
    
    # Delta(1/Tm) par accès privé pour vérif intermédiaire
    delta_inv_tm = salt_model._calc_delta_inv_tm(
        f_gc=f_gc,
        n_bp=n_bp,
        mon_total=mon_molar,
        mg_free=mg_molar,
        na_ref_molar=1.0
    )
    
    # Attendus
    expected_delta_inv = 7.0537e-5
    # Tolérance assez fine pour valider le calcul
    assert abs(delta_inv_tm - expected_delta_inv) < 1e-7, f"Delta(1/Tm) diverge : {delta_inv_tm} vs {expected_delta_inv}"
    
    # Tm prédit
    tm_pred = salt_model.correct_tm(
        tm_ref_celsius=tm_1m_celsius,
        f_gc=f_gc,
        n_bp=n_bp,
        na_molar=0.0,
        mg_molar=mg_molar,
        na_ref_molar=1.0
    )
    expected_tm = 67.9
    assert abs(tm_pred - expected_tm) < 0.1, f"Tm prédit diverge : {tm_pred} vs {expected_tm}"


def test_owczarzy2008_mg_pure():
    """Validation de l'équation 16 pure (Mg2+) sur owczarzy2008_mg.csv."""
    ref_mg_path = pathlib.Path(__file__).parents[2] / "validation/reference_data/owczarzy2008_mg.csv"
    ref_na_path = pathlib.Path(__file__).parents[2] / "validation/reference_data/owczarzy2004_salt.csv"
    
    df_mg = pd.read_csv(ref_mg_path)
    df_na = pd.read_csv(ref_na_path)
    
    salt_model = UnifiedSaltModel()
    
    mg_cols = [
        ('Tm_0.5mM_Mg_C', 0.5e-3),
        ('Tm_1.5mM_Mg_C', 1.5e-3),
        ('Tm_3.0mM_Mg_C', 3.0e-3),
        ('Tm_10mM_Mg_C', 10e-3),
        ('Tm_20mM_Mg_C', 20e-3),
        ('Tm_50mM_Mg_C', 50e-3),
        ('Tm_125mM_Mg_C', 125e-3),
        ('Tm_300mM_Mg_C', 300e-3),
        ('Tm_600mM_Mg_C', 600e-3)
    ]
    
    errors_mg = {col: [] for col, _ in mg_cols}
    
    for _, row in df_mg.iterrows():
        seq = row['sequence']
        f_gc = row['f_GC']
        n_bp = row['length']
        
        # Référence à 1 M (colonne 1.02 M Na+ du dataset 2004)
        na_match = df_na[df_na['sequence'] == seq]
        if na_match.empty:
            continue
        tm_ref_1020 = na_match.iloc[0]['Tm_1020mM_Na_C']
        
        # La référence expérimentale Na+ était à 1.02 M
        # Le SaltModel utilise na_ref_molar=1.02 pour corriger depuis cette référence
        for col_name, mg_molar in mg_cols:
            tm_pred = salt_model.correct_tm(
                tm_ref_celsius=tm_ref_1020,
                f_gc=f_gc,
                n_bp=n_bp,
                na_molar=0.0,
                mg_molar=mg_molar,
                na_ref_molar=1.02
            )
            tm_exp = row[col_name]
            errors_mg[col_name].append(abs(tm_pred - tm_exp))

    print("\\n\\n=== VALIDATION OWCZARZY 2008 (Magnésium Pur) ===")
    for col_name, mg_molar in mg_cols:
        mean_err = np.mean(errors_mg[col_name])
        max_err = np.max(errors_mg[col_name])
        print(f"[{mg_molar*1000:^5.1f} mM Mg] Moyenne: {mean_err:.2f} C | Max: {max_err:.2f} C")
        # Attendu: Moyenne 0.6 à 1.0, Max vers 2.0
        assert mean_err < 1.5, f"Erreur Mg trop élevée pour {col_name}: {mean_err:.2f} C"


def test_owczarzy2008_mixed_regime():
    """Validation du régime mixte et du sélecteur."""
    ref_mixed_path = pathlib.Path(__file__).parents[2] / "validation/reference_data/owczarzy2008_mixed.csv"
    ref_na_path = pathlib.Path(__file__).parents[2] / "validation/reference_data/owczarzy2004_salt.csv"
    
    df_mixed = pd.read_csv(ref_mixed_path)
    df_na = pd.read_csv(ref_na_path)
    
    salt_model = UnifiedSaltModel()
    
    mg_cols = [
        ('Tm_0.5mM_Mg', 0.5e-3),
        ('Tm_1.5mM_Mg', 1.5e-3),
        ('Tm_3.0mM_Mg', 3.0e-3),
        ('Tm_10mM_Mg', 10e-3),
        ('Tm_20mM_Mg', 20e-3),
        ('Tm_50mM_Mg', 50e-3),
        ('Tm_125mM_Mg', 125e-3)
    ]
    
    errors_mixed = []
    errors_by_r = {"na_dom": [], "mixed": [], "mg_dom": []}
    
    for _, row in df_mixed.iterrows():
        seq = row['sequence']
        f_gc = row['f_GC']
        n_bp = row['length']
        mon_molar = row['mon_mM'] / 1000.0
        
        na_match = df_na[df_na['sequence'] == seq]
        if na_match.empty:
            continue
        tm_ref_1020 = na_match.iloc[0]['Tm_1020mM_Na_C']
        
        for col_name, mg_molar in mg_cols:
            tm_exp = row[col_name]
            if pd.isna(tm_exp):
                continue
                
            tm_pred = salt_model.correct_tm(
                tm_ref_celsius=tm_ref_1020,
                f_gc=f_gc,
                n_bp=n_bp,
                na_molar=mon_molar,
                mg_molar=mg_molar,
                na_ref_molar=1.02
            )
            
            err = abs(tm_pred - tm_exp)
            errors_mixed.append(err)
            
            # Catégorisation pour l'affichage
            r_ratio = math.sqrt(mg_molar) / mon_molar
            if r_ratio < 0.22:
                errors_by_r["na_dom"].append(err)
            elif r_ratio < 6.0:
                errors_by_r["mixed"].append(err)
            else:
                errors_by_r["mg_dom"].append(err)
                
    print("\\n=== VALIDATION OWCZARZY 2008 (Régime Mixte) ===")
    mean_all = np.mean(errors_mixed)
    max_all = np.max(errors_mixed)
    print(f"Global: Moyenne: {mean_all:.2f} C | Max: {max_all:.2f} C")
    
    print("Ventilation par Régime (Ratio de Compétition R):")
    if errors_by_r["na_dom"]:
        print(f"R < 0.22 (Na+ dominant)  : {np.mean(errors_by_r['na_dom']):.2f} C (n={len(errors_by_r['na_dom'])})")
    if errors_by_r["mixed"]:
        print(f"Mixte (0.22 <= R < 6.0)  : {np.mean(errors_by_r['mixed']):.2f} C (n={len(errors_by_r['mixed'])})")
    if errors_by_r["mg_dom"]:
        print(f"R >= 6.0 (Mg2+ dominant) : {np.mean(errors_by_r['mg_dom']):.2f} C (n={len(errors_by_r['mg_dom'])})")
        
    assert mean_all < 1.0, f"Moyenne globale {mean_all:.2f} °C > 1.0 °C (attendu < 1.0 °C)"
    assert max_all < 3.5, f"Erreur maximale {max_all:.2f} °C > 3.5 °C (attendu < 3.5 °C)"


def test_solver_path_consistency():
    """Vérifie que la correction entropique (Solveur) produit le même Tm que la formule empirique."""
    salt_model = UnifiedSaltModel()
    native_backend = NativeBackend(default_ct_molar=2e-6)
    corrected_backend = SaltCorrectedBackend(backend=native_backend, salt_model=salt_model)
    
    seq = "GATGCGCTCG"
    f_gc = 0.7
    n_bp = 10
    
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    seq_comp = "".join(complement[b] for b in reversed(seq))
    
    # 1. Native à 1 M
    res_1m = native_backend.calc_heterodimer(seq, seq_comp, ct_molar=2e-6)
    tm_1m = res_1m.tm_celsius
    
    # Conditions de test
    conditions = [
        (100.0, 0.0),   # 100 mM Na
        (0.0, 3.0),     # 3 mM Mg
        (55.0, 1.5)     # Mixte: 55 mM Na, 1.5 mM Mg
    ]
    
    for na_mm, mg_mm in conditions:
        # Tm empirique
        tm_empiric = salt_model.correct_tm(
            tm_ref_celsius=tm_1m,
            f_gc=f_gc,
            n_bp=n_bp,
            na_molar=na_mm/1000.0,
            mg_molar=mg_mm/1000.0,
            na_ref_molar=1.0
        )
        
        # Tm reconstruit via dH et dS corrigés (Solveur)
        res_corr = corrected_backend.calc_heterodimer(
            seq, seq_comp, ct_molar=2e-6, na_mm=na_mm, mg_mm=mg_mm
        )
        tm_entropic = res_corr.tm_celsius
        
        assert abs(tm_empiric - tm_entropic) < 1e-4, (
            f"Divergence de Tm entre empirique ({tm_empiric}) et entropique ({tm_entropic}) "
            f"pour Na={na_mm}, Mg={mg_mm}"
        )
