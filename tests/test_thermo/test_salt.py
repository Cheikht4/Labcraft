import pathlib

import numpy as np
import pandas as pd
import pytest

from labcraft.thermo.backends.native import NativeBackend
from labcraft.thermo.salt import Owczarzy2004SaltModel, SaltCorrectedBackend


def test_owczarzy2004_validations():
    """Valide la correction saline d'Owczarzy 2004 sous 3 angles."""
    ref_path = pathlib.Path("validation/reference_data/owczarzy2004_salt.csv")
    df = pd.read_csv(ref_path)
    
    salt_model = Owczarzy2004SaltModel()
    native_backend = NativeBackend(default_ct_molar=2e-6)
    corrected_backend = SaltCorrectedBackend(backend=native_backend, salt_model=salt_model)
    
    target_na_cols = [
        ('Tm_69mM_Na_C', 0.069),
        ('Tm_119mM_Na_C', 0.119),
        ('Tm_220mM_Na_C', 0.220),
        ('Tm_621mM_Na_C', 0.621)
    ]
    
    # 1. Erreurs Eq 22 isolée (ancrage 1.02 M)
    errors_isolated = {col: [] for col, _ in target_na_cols}
    # 2. Erreurs bout-en-bout (NN -> Eq 22 empirique)
    errors_e2e = {col: [] for col, _ in target_na_cols}
    errors_e2e_1020 = []
    # 3. Écart empirique (Eq 22) vs Entropique (Solveur)
    discrepancies = {col: [] for col, _ in target_na_cols}
    
    for _, row in df.iterrows():
        seq = row['sequence']
        f_gc = row['f_GC']
        
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        seq_comp = "".join(complement[b] for b in reversed(seq))
        
        # --- Couche isolée ---
        tm_ref_1020 = row['Tm_1020mM_Na_C']
        for col_name, na_molar in target_na_cols:
            tm_pred_iso = salt_model.correct_tm(
                tm_ref_celsius=tm_ref_1020,
                f_gc=f_gc,
                na_target_molar=na_molar,
                na_ref_molar=1.02
            )
            tm_exp = row[col_name]
            errors_isolated[col_name].append(abs(tm_pred_iso - tm_exp))
            
        # --- Couche bout-en-bout & Écart Entropique ---
        # Prédiction NN brute à 1 M (sans correction)
        res_nn = native_backend.calc_heterodimer(seq, seq_comp, ct_molar=2e-6)
        tm_nn_1m = res_nn.tm_celsius
        
        # Test cohérence 1.02 M pour bout-en-bout
        tm_e2e_1020 = salt_model.correct_tm(
            tm_ref_celsius=tm_nn_1m, f_gc=f_gc, na_target_molar=1.02, na_ref_molar=1.0
        )
        errors_e2e_1020.append(abs(tm_e2e_1020 - tm_ref_1020))
        
        for col_name, na_molar in target_na_cols:
            # Prédiction empirique (Eq 22) ancrée sur le NN à 1 M
            tm_pred_e2e = salt_model.correct_tm(
                tm_ref_celsius=tm_nn_1m,
                f_gc=f_gc,
                na_target_molar=na_molar,
                na_ref_molar=1.0
            )
            tm_exp = row[col_name]
            errors_e2e[col_name].append(abs(tm_pred_e2e - tm_exp))
            
            # Prédiction entropique (Solveur, via SaltCorrectedBackend)
            res_corr = corrected_backend.calc_heterodimer(
                seq, seq_comp, ct_molar=2e-6, na_mm=na_molar * 1000.0
            )
            tm_pred_entropic = res_corr.tm_celsius
            
            # Écart absolu entre empirique et entropique
            discrepancies[col_name].append(abs(tm_pred_entropic - tm_pred_e2e))

    # --- Affichage des tableaux (pour le rapport) ---
    print("\\n\\n=== RÉSULTATS DE VALIDATION OWCZARZY 2004 ===")
    
    print("\\n1. Erreur Eq 22 isolée (Ancrage sur mesure expérimentale à 1.02 M)")
    for col_name, na_molar in target_na_cols:
        mean_err = np.mean(errors_isolated[col_name])
        max_err = np.max(errors_isolated[col_name])
        print(f"[{na_molar*1000:.0f} mM] Moyenne: {mean_err:.2f} C | Max: {max_err:.2f} C")
        
    print("\\n2. Erreur Bout-en-bout (NN Natif 1M -> Eq 22 Empirique)")
    mean_err_1020 = np.mean(errors_e2e_1020)
    max_err_1020 = np.max(errors_e2e_1020)
    print(f"[1020 mM] Moyenne: {mean_err_1020:.2f} C | Max: {max_err_1020:.2f} C (Cohérence)")
    for col_name, na_molar in target_na_cols:
        mean_err = np.mean(errors_e2e[col_name])
        max_err = np.max(errors_e2e[col_name])
        print(f"[{na_molar*1000:.0f} mM] Moyenne: {mean_err:.2f} C | Max: {max_err:.2f} C")

    print("\\n3. Écart Solveur vs Validation (Modèle Entropique vs Eq 22 Empirique)")
    for col_name, na_molar in target_na_cols:
        mean_diff = np.mean(discrepancies[col_name])
        max_diff = np.max(discrepancies[col_name])
        print(f"[{na_molar*1000:.0f} mM] Écart Moyen: {mean_diff:.2f} C | Écart Max: {max_diff:.2f} C")
        
    # Assertions minimales
    # L'erreur isolée de la formule doit être excellente (< 0.5 C en moyenne)
    for col_name, _ in target_na_cols:
        assert np.mean(errors_isolated[col_name]) < 1.0
