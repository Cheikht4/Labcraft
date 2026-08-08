import sys
import os
sys.path.insert(0, os.path.abspath('src'))
from labcraft.thermo.backends.native import _NN_PARAMS, _SEQ_TO_NN_KEY, _INIT_GC, _INIT_AT

def extract_pairs(structure: str) -> dict[int, int]:
    pairs = {}
    stack = []
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                j = stack.pop()
                pairs[j] = i
                pairs[i] = j
    return pairs

def estimate_helix_thermo(seq: str, structure: str):
    seq = seq.replace('&', '+').upper()
    pairs = extract_pairs(structure)
    
    dh_total = 0.0
    ds_total = 0.0
    n_bp = 0
    gc_count = 0
    
    for i in range(len(seq) - 1):
        if i in pairs and (i+1) in pairs:
            j = pairs[i]
            j_minus_1 = pairs[i+1]
            if j_minus_1 == j - 1 and i < j and (i+1) < (j-1):
                dinuc_top = seq[i:i+2]
                try:
                    key = _SEQ_TO_NN_KEY[dinuc_top]
                    dh, ds = _NN_PARAMS[key]
                    dh_total += dh
                    ds_total += ds
                except KeyError:
                    pass
                
    for i in pairs:
        if i < pairs[i]:
            j = pairs[i]
            is_left_end = True
            if (i-1) in pairs and pairs[i-1] == j+1:
                is_left_end = False
                
            is_right_end = True
            if (i+1) in pairs and pairs[i+1] == j-1:
                is_right_end = False
                
            base = seq[i]
            # When looking at base pair (i, j), if it's an end, its contribution depends on the pair.
            # But wait, initiation penalty is per helix.
            # A helix is a contiguous stretch of base pairs.
            # So a helix has exactly two ends. 
            # Our code assigns an initiation penalty to EACH end (left and right).
            # Wait, SantaLucia states initiation penalty is per DUPLEX, not per end!
            # Wait, no. "For duplexes with terminal A-T pairs...".
            # The penalties _INIT_AT (2.3, 4.1) and _INIT_GC (0.1, -2.8) are per helix end!
            # See native.py: `for end in (seq[0], seq[-1]):`
            # Yes, each end contributes. So for a helix, we penalize both ends.
            if is_left_end:
                pair = (seq[i], seq[j])
                if pair in [('G', 'C'), ('C', 'G')]:
                    dh_total += _INIT_GC[0]; ds_total += _INIT_GC[1]
                else:
                    dh_total += _INIT_AT[0]; ds_total += _INIT_AT[1]
            if is_right_end:
                pair = (seq[i], seq[j])
                if pair in [('G', 'C'), ('C', 'G')]:
                    dh_total += _INIT_GC[0]; ds_total += _INIT_GC[1]
                else:
                    dh_total += _INIT_AT[0]; ds_total += _INIT_AT[1]
                    
            n_bp += 1
            if seq[i] in ('G', 'C'):
                gc_count += 1
                
    f_gc = gc_count / n_bp if n_bp > 0 else 0.0
    return dh_total, ds_total, n_bp, f_gc

print(estimate_helix_thermo("ATGC&GCAT", "((((+))))"))
