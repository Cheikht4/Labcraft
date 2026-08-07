"""Ranking of optimized variants.

Classement des variants optimisés.
"""
from typing import List, Tuple, Dict, Any

def score_variants(
    valid_variants: List[Tuple[str, List[Tuple[int, str, str]], float]]
) -> List[Dict[str, Any]]:
    """
    Trie les variants survivants.
    
    Args:
        valid_variants: Liste de tuples (variant_seq, mutations, worst_dg_3p).
        
    Returns:
        Liste de dictionnaires décrivant les variants, triée du meilleur au pire.
    """
    scored = []
    
    for variant_seq, mutations, worst_dg_3p in valid_variants:
        num_mutations = len(mutations)
        
        # Compte combien de mutations introduisent un A ou un T
        num_at = sum(1 for _, _, new_b in mutations if new_b in ('A', 'T'))
        
        scored.append({
            'sequence': variant_seq,
            'mutations': mutations,
            'worst_dg_3p': worst_dg_3p,
            'num_mutations': num_mutations,
            'num_at': num_at
        })
        
    # Tri : 
    # 1. Moins de mutations (croissant)
    # 2. Pire dG_3p le plus haut (décroissant)
    # 3. Plus de mutations vers A/T (décroissant)
    scored.sort(key=lambda x: (x['num_mutations'], -x['worst_dg_3p'], -x['num_at']))
    
    return scored
