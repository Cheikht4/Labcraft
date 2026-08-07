"""Constraints and filters for sequence optimization.

Filtres et contraintes pour l'optimisation des séquences.
"""
from typing import List, Dict, Tuple, Any
import RNA

from labcraft.diagnostics.amplifiable_dimer import is_amplifiable_dimer
from labcraft.diagnostics.enzyme import PolymeraseProfile

def evaluate_variant_dimers(
    variant_seq: str,
    original_primers: List[Dict[str, Any]],
    modified_idx: int,
    enzyme: PolymeraseProfile,
    temp_celsius: float
) -> List[Tuple[int, float]]:
    """
    Évalue tous les dimères impliquant le variant (hétérodimères et homodimère).
    
    Args:
        variant_seq: La séquence mutée.
        original_primers: Liste des amorces du panel (dict avec 'sequence').
        modified_idx: Index de l'amorce modifiée dans original_primers.
        enzyme: Profil de polymérase.
        temp_celsius: Température d'évaluation.
        
    Returns:
        Liste de (partner_idx, dg_3p) pour chaque dimère amplifiable trouvé.
    """
    amplifiable = []
    
    # Doit être exécuté dans un contexte dna_params(temp_celsius)
    for j, p in enumerate(original_primers):
        seq_partner = variant_seq if j == modified_idx else p['sequence']
        
        # cofold
        seq_concat = f"{variant_seq}&{seq_partner}"
        structure, mfe = RNA.cofold(seq_concat)
        struct_clean = structure.replace('&', '')
        
        is_amp, min_dg_3p = is_amplifiable_dimer(
            variant_seq, seq_partner, struct_clean, mfe, enzyme, temp_celsius=temp_celsius
        )
        
        if is_amp:
            amplifiable.append((j, min_dg_3p))
            
    return amplifiable

def filter_no_new_artefacts(
    variant_seq: str,
    original_primers: List[Dict[str, Any]],
    modified_idx: int,
    enzyme: PolymeraseProfile,
    temp_celsius: float
) -> Tuple[bool, float]:
    """
    Filtre un variant s'il possède un quelconque dimère amplifiable.
    Retourne (is_valid, worst_dg_3p).
    Si is_valid est True, worst_dg_3p est le pire dG_3p (le plus petit/négatif)
    parmi les dimères (non amplifiables) évalués.
    """
    worst_dg_3p = 0.0
    
    for j, p in enumerate(original_primers):
        seq_partner = variant_seq if j == modified_idx else p['sequence']
        
        seq_concat = f"{variant_seq}&{seq_partner}"
        structure, mfe = RNA.cofold(seq_concat)
        struct_clean = structure.replace('&', '')
        
        is_amp, min_dg_3p = is_amplifiable_dimer(
            variant_seq, seq_partner, struct_clean, mfe, enzyme, temp_celsius=temp_celsius
        )
        
        if is_amp:
            return False, min_dg_3p
            
        if min_dg_3p < worst_dg_3p:
            worst_dg_3p = min_dg_3p
            
    return True, worst_dg_3p
