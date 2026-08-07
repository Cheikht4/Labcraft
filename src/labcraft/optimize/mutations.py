"""Mutation enumeration for sequence optimization.

Énumération des mutations pour l'optimisation des séquences.
"""
from typing import List, Tuple
import itertools

def enumerate_variants(sequence: str, max_mutations: int = 2, window_3p: int = 6) -> List[Tuple[str, List[Tuple[int, str, str]]]]:
    """
    Enumerate all substitution variants of a sequence within the 3' window.
    
    Args:
        sequence: Original primer sequence (5' to 3').
        max_mutations: Maximum number of simultaneous substitutions (1 or 2).
        window_3p: Number of bases at the 3' end to consider for mutations.
        
    Returns:
        List of tuples (mutated_sequence, list_of_mutations), where list_of_mutations
        is a list of (index, old_base, new_base).
        Returns the original sequence as well (0 mutations).
    """
    if window_3p > len(sequence):
        window_3p = len(sequence)
        
    start_idx = len(sequence) - window_3p
    bases = ['A', 'C', 'G', 'T']
    
    variants = []
    
    # 0 mutations (original)
    variants.append((sequence, []))
    
    # 1 mutation
    if max_mutations >= 1:
        for i in range(start_idx, len(sequence)):
            old_b = sequence[i]
            for new_b in bases:
                if new_b != old_b:
                    mut_seq = sequence[:i] + new_b + sequence[i+1:]
                    variants.append((mut_seq, [(i, old_b, new_b)]))
                    
    # 2 mutations
    if max_mutations >= 2:
        for i in range(start_idx, len(sequence)):
            for j in range(i + 1, len(sequence)):
                old_b1 = sequence[i]
                old_b2 = sequence[j]
                for new_b1 in bases:
                    if new_b1 != old_b1:
                        for new_b2 in bases:
                            if new_b2 != old_b2:
                                mut_list = list(sequence)
                                mut_list[i] = new_b1
                                mut_list[j] = new_b2
                                mut_seq = "".join(mut_list)
                                variants.append((mut_seq, [(i, old_b1, new_b1), (j, old_b2, new_b2)]))
                                
    return variants
