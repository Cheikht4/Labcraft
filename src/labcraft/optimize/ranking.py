"""Ranking of optimized variants.

Classement des variants optimisés.
"""
from typing import List, Tuple, Dict, Any

def score_variants(
    valid_variants: List[Tuple[str, List[Tuple[int, str, str]], float, float]]
) -> List[Dict[str, Any]]:
    """
    Trie les variants survivants.
    
    Args:
        valid_variants: Liste de tuples (variant_seq, mutations, worst_dg_3p, intrinsic_dg).
        
    Returns:
        Liste de dictionnaires décrivant les variants, triée du meilleur au pire.
    """
    scored = []
    
    for variant_seq, mutations, worst_dg_3p, intrinsic_dg in valid_variants:
        num_mutations = len(mutations)
        
        scored.append({
            'sequence': variant_seq,
            'mutations': mutations,
            'worst_dg_3p': worst_dg_3p,
            'intrinsic_dg': intrinsic_dg,
            'num_mutations': num_mutations
        })
        
    # Tri : 
    # 1. Moins de mutations (croissant)
    # 2. Stabilité intrinsèque la moins collante (le dG le plus positif/haut) (décroissant)
    # 3. Pire dG_3p le plus haut (décroissant)
    scored.sort(key=lambda x: (x['num_mutations'], -x['intrinsic_dg'], -x['worst_dg_3p']))
    
    return scored
