"""Main optimizer logic.

Logique principale de l'optimiseur.
"""
from typing import List, Dict, Any, Tuple

from labcraft.optimize.mutations import enumerate_variants
from labcraft.optimize.constraints import filter_no_new_artefacts, calc_intrinsic_3p_dg
from labcraft.optimize.ranking import score_variants
from labcraft.diagnostics.enzyme import BST
from labcraft.thermo.vienna import dna_params

def optimize_primer(
    primer_name: str,
    primers: List[Dict[str, Any]],
    temp_celsius: float = 63.0,
    max_mutations: int = 2,
    window_3p: int = 6
) -> List[Dict[str, Any]]:
    """
    Optimise une amorce ciblée dans un panel pour supprimer les dimères amplifiables.
    
    Args:
        primer_name: Nom de l'amorce à optimiser.
        primers: Liste des amorces du panel (dict avec 'name', 'sequence', etc.).
        temp_celsius: Température de la réaction.
        max_mutations: Nombre max de mutations (1 ou 2).
        window_3p: Taille de la fenêtre 3' à explorer.
        
    Returns:
        Une liste des meilleurs variants (top 3) avec leurs métadonnées.
    """
    # Trouver l'amorce
    modified_idx = -1
    for i, p in enumerate(primers):
        if p['name'] == primer_name:
            modified_idx = i
            break
            
    if modified_idx == -1:
        raise ValueError(f"Amorce '{primer_name}' introuvable dans le panel.")
        
    original_seq = primers[modified_idx]['sequence']
    
    # 1. Énumération des variants
    variants = enumerate_variants(original_seq, max_mutations, window_3p)
    
    # Calcul de la stabilité intrinsèque de l'original
    orig_intrinsic_dg = calc_intrinsic_3p_dg(original_seq, window_3p, temp_celsius)
    
    # 2. Application des filtres durs (évaluation thermodynamique)
    valid_variants = []
    
    # L'évaluation se fait dans le contexte ViennaRNA
    with dna_params(temp_celsius):
        for seq, muts in variants:
            # Filtre : la stabilité intrinsèque ne doit pas être pire (plus négative) que l'original
            # On autorise un léger delta dû aux arrondis (e.g., -0.1) ou on force <= orig ?
            # Le but est "rejette ou pénalise fortement tout variant qui augmente la stabilité 3' intrinsèque"
            # On rejette si c'est strictement plus négatif.
            intrinsic_dg = calc_intrinsic_3p_dg(seq, window_3p, temp_celsius)
            if intrinsic_dg < orig_intrinsic_dg:
                continue
                
            # Filtre strict : aucun dimère amplifiable impliquant le variant (hétéro et homo)
            is_valid, worst_dg_3p = filter_no_new_artefacts(
                seq, primers, modified_idx, BST, temp_celsius
            )
            
            if is_valid:
                valid_variants.append((seq, muts, worst_dg_3p, intrinsic_dg))
                
    # 3. Classement des survivants
    ranked = score_variants(valid_variants)
    
    # Retourne les 3 meilleurs
    return ranked[:3]
