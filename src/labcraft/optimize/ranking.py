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
    scored: List[Dict[str, Any]] = []
    
    for variant_seq, mutations, worst_dg_3p, intrinsic_dg in valid_variants:
        num_mutations = len(mutations)
        
        scored.append({
            'sequence': variant_seq,
            'mutations': mutations,
            'worst_dg_3p': worst_dg_3p,
            'intrinsic_dg': intrinsic_dg,
            'num_mutations': num_mutations
        })
        
    # Tri / Sorting:
    # 1. Moins de mutations (croissant) / Fewer mutations (ascending)
    # 2. Stabilité intrinsèque la moins collante (dG le plus haut) / Least sticky intrinsic stability (highest dG)
    # 3. Pire dG_3p le plus haut / Highest worst dG_3p
    scored.sort(
        key=lambda x: (
            int(x['num_mutations']),
            -float(x['intrinsic_dg']),
            -float(x['worst_dg_3p']),
        )
    )
    
    return scored
