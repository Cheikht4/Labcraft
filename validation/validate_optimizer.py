"""Validation script for the sequence optimizer.

Script de démonstration de l'optimiseur sur le jeu Parida DEN-3 (FLP).
"""
import csv
from typing import List, Dict, Any

from labcraft.optimize.optimizer import optimize_primer
from labcraft.optimize.recommendations import generate_recommendations
from labcraft.metrics.verdict import PanelVerdict, PrimerIssue

def load_parida_csv(path: str) -> Dict[str, List[Dict[str, Any]]]:
    from collections import defaultdict
    serotypes = defaultdict(list)
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            serotypes[row['serotype']].append(row)
    return serotypes

def main():
    path = "validation/reference_data/parida2005_lamp_primers.csv"
    serotypes = load_parida_csv(path)
    den3_primers = serotypes["DEN-3"]
    
    print("--- Démonstration de l'Optimiseur sur Parida DEN-3 ---")
    
    # 1. État initial
    target_primer = "FLP_DEN-3"
    original_seq = ""
    for p in den3_primers:
        if p['name'] == target_primer:
            original_seq = p['sequence']
            break
            
    print(f"Amorce cible : {target_primer}")
    print(f"Séquence originale : {original_seq}")
    print("\nAVERTISSEMENT : FLP est une amorce de boucle dont la cible est la région boucle de l'amplicon.")
    print("La séquence cible n'étant pas fournie dans ce test unitaire, la préservation de la fixation")
    print("à la boucle ne peut pas être formellement vérifiée. En production, fournissez la séquence cible")
    print("pour activer le filtre de préservation d'affinité.\n")
    
    # 2. Recommandations
    # On simule un verdict avec des dimères amplifiables pour générer les recommandations
    issues = [
        PrimerIssue(primer_name="FLP_DEN-3", target_site="N/A", occupation=0.0, is_critical=True, cause="Dimère amplifiable"),
        PrimerIssue(primer_name="BIP_DEN-3", target_site="N/A", occupation=0.0, is_critical=True, cause="Dimère bloquant")
    ]
    verdict = PanelVerdict("FAILURE", issues, "Échec par dimères multiples.")
    recs = generate_recommendations(verdict)
    
    print("Recommandations non séquentielles :")
    for r in recs:
        print(f"  - {r}")
    print()
    
    # 3. Optimisation
    print(f"Lancement de l'optimisation pour {target_primer} (fenêtre 3' de 6 nucléotides)...")
    top_variants = optimize_primer(
        primer_name=target_primer,
        primers=den3_primers,
        temp_celsius=63.0,
        max_mutations=2,
        window_3p=6
    )
    
    # 4. Affichage Avant/Après
    print("\n--- Tableau Avant/Après (Top 3 Variants) ---")
    print(f"{'Rang':<4} | {'Séquence':<20} | {'Muts':<4} | {'A/T':<3} | {'Pire ΔG_3\' (kcal/mol)'}")
    print("-" * 65)
    
    print(f"{'Orig':<4} | {original_seq:<20} | {'0':<4} | {'0':<3} | {-3.65:.2f} (Échec)")
    
    for i, var in enumerate(top_variants):
        seq = var['sequence']
        muts = var['num_mutations']
        at = var['num_at']
        dg = var['worst_dg_3p']
        
        # Surligner la mutation dans la séquence pour l'affichage (optionnel, on laisse brut)
        print(f"{i+1:<4} | {seq:<20} | {muts:<4} | {at:<3} | {dg:>.2f}")
        
    if top_variants:
        best_seq = top_variants[0]['sequence']
        
        # 5. Contrôle complet sur le panel
        print(f"\n--- Contrôle du Panel Modifié (Variant 1) ---")
        
        # On remplace l'amorce dans le panel
        modified_primers = list(den3_primers)
        for i, p in enumerate(modified_primers):
            if p['name'] == target_primer:
                p_copy = dict(p)
                p_copy['sequence'] = best_seq
                modified_primers[i] = p_copy
                
        # On utilise le script validate_parida pour réévaluer
        import sys
        import os
        sys.path.append(os.path.join(os.path.dirname(__file__)))
        from validate_parida import run_analysis
        dimers_after = run_analysis("DEN-3_Optimized", modified_primers, skip_flp=False)
        print(f"Nombre de dimères amplifiables dans le panel complet après optimisation : {len(dimers_after)}")
        
        if len(dimers_after) == 0:
            print("SUCCÈS : Le panel DEN-3 optimisé ne génère plus aucun dimère amplifiable.")
        else:
            print("ÉCHEC : Des dimères amplifiables persistent.")
            for d in dimers_after:
                print(f"  -> {d['primer1']} / {d['primer2']} : dG_3' = {d['dg_3p']:.2f}")

if __name__ == "__main__":
    main()
