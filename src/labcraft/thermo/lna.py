import csv
from pathlib import Path
from typing import Tuple, List

def parse_lna_sequence(seq: str) -> Tuple[str, List[int]]:
    """
    Parse a sequence that may contain LNA bases indicated by a '+' prefix.
    For example: 'AG+C+GTA' -> ('AGCGTA', [2, 3])
    
    Args:
        seq: The input sequence which may contain '+' signs.
        
    Returns:
        A tuple of (bare_sequence, list_of_lna_indices)
    """
    bare_seq = []
    lna_positions = []
    
    i = 0
    bare_idx = 0
    while i < len(seq):
        if seq[i] == '+':
            # The next base is an LNA
            if i + 1 < len(seq):
                bare_seq.append(seq[i+1])
                lna_positions.append(bare_idx)
                bare_idx += 1
                i += 2
            else:
                raise ValueError("Dangling '+' at the end of the sequence")
        else:
            bare_seq.append(seq[i])
            bare_idx += 1
            i += 1
            
    return "".join(bare_seq), lna_positions

_LNA_PARAMS_MXL = {}
_LNA_PARAMS_XLN = {}
_LNA_LOADED = False

def load_lna_params():
    global _LNA_LOADED
    if _LNA_LOADED:
        return
        
    csv_path = Path(__file__).parent.parent / "data" / "thermo" / "mctigue2004_lna_nn.csv"
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            half = row["half"].strip()
            base5 = row["base5"].strip()
            base3 = row["base3_or_lna"].strip()
            dh = float(row["ddH_kcal"])
            ds = float(row["ddS_cal"])
            dg = float(row["ddG37_kcal"])
            
            if half == "MXL":
                _LNA_PARAMS_MXL[(base5, base3)] = (dh, ds)
            elif half == "XLN":
                _LNA_PARAMS_XLN[(base5, base3)] = (dh, ds)
                
    _LNA_LOADED = True

