def dotbracket_to_alignment(seq_a: str, seq_b: str, struct: str) -> str:
    """
    Convertit une structure dot-bracket (sans le '&') de deux brins
    en un alignement ASCII 2D antiparallèle.
    Brin du haut (A) : 5' -> 3'
    Brin du bas (B)  : 3' -> 5'
    """
    # 1. Extraction des paires
    stack = []
    pairs = {}
    for i, c in enumerate(struct):
        if c == '(': 
            stack.append(i)
        elif c == ')':
            if stack:
                j = stack.pop()
                pairs[i] = j
                pairs[j] = i
                
    # 2. Alignement
    out_a = ""
    out_lines = ""
    out_b = ""
    
    idx_a = 0
    idx_b = len(seq_a) + len(seq_b) - 1
    
    while idx_a < len(seq_a) or idx_b >= len(seq_a):
        # Sont-ils appariés l'un à l'autre ?
        if idx_a < len(seq_a) and idx_b >= len(seq_a) and pairs.get(idx_a) == idx_b:
            out_a += seq_a[idx_a]
            out_b += seq_b[idx_b - len(seq_a)]
            out_lines += "|"
            idx_a += 1
            idx_b -= 1
        else:
            # Sinon, on détermine qui doit avancer en cherchant la prochaine paire inter-moléculaire
            next_a_pair = None
            next_b_pair = None
            
            for i in range(idx_a, len(seq_a)):
                if pairs.get(i) is not None and pairs[i] >= len(seq_a):
                    next_a_pair = i
                    break
            for j in range(idx_b, len(seq_a) - 1, -1):
                if pairs.get(j) is not None and pairs[j] < len(seq_a):
                    next_b_pair = j
                    break
                    
            if next_a_pair is not None and next_b_pair is not None:
                if idx_a < next_a_pair:
                    out_a += seq_a[idx_a]
                    out_b += "-"
                    out_lines += " "
                    idx_a += 1
                elif idx_b > next_b_pair:
                    out_a += "-"
                    out_b += seq_b[idx_b - len(seq_a)]
                    out_lines += " "
                    idx_b -= 1
                else:
                    # Cas théorique de croisement (non géré par dot-bracket standard)
                    out_a += seq_a[idx_a]
                    out_b += "-"
                    out_lines += " "
                    idx_a += 1
            else:
                # Fin des paires inter-moléculaires
                if idx_a < len(seq_a):
                    out_a += seq_a[idx_a]
                    out_b += " "
                    out_lines += " "
                    idx_a += 1
                elif idx_b >= len(seq_a):
                    out_a += " "
                    out_b += seq_b[idx_b - len(seq_a)]
                    out_lines += " "
                    idx_b -= 1

    # On trim les espaces à droite pour la propreté
    return f"5' {out_a} 3'\n   {out_lines}\n3' {out_b} 5'"

from typing import List, Dict, Optional

def get_alignment_columns(seq_a: str, seq_b: str, struct: str, extensible_strand: Optional[str] = None) -> List[Dict]:
    """
    Construit une liste de colonnes pour le rendu graphique de l'alignement.
    Si extensible_strand ('a' ou 'b') est fourni, annote le 3' extensible et la matrice.
    """
    # 1. Extraction des paires
    stack = []
    pairs = {}
    struct_clean = struct.replace('&', '')
    for i, c in enumerate(struct_clean):
        if c == '(': 
            stack.append(i)
        elif c == ')':
            if stack:
                j = stack.pop()
                pairs[i] = j
                pairs[j] = i
                
    columns = []
    l_a = len(seq_a)
    l_b = len(seq_b)
    idx_a = 0
    idx_b = l_a + l_b - 1
    
    # 2. Alignement
    while idx_a < l_a or idx_b >= l_a:
        col = {"top": "-", "bottom": "-", "paired": False, "role": "", "idx_a": None, "idx_b": None, "is_truncation": False}
        if idx_a < l_a and idx_b >= l_a and pairs.get(idx_a) == idx_b:
            col["top"] = seq_a[idx_a]
            col["bottom"] = seq_b[idx_b - l_a]
            col["paired"] = True
            col["idx_a"] = idx_a
            col["idx_b"] = idx_b
            idx_a += 1
            idx_b -= 1
        else:
            next_a_pair = None
            next_b_pair = None
            
            for i in range(idx_a, l_a):
                if pairs.get(i) is not None and pairs[i] >= l_a:
                    next_a_pair = i
                    break
            for j in range(idx_b, l_a - 1, -1):
                if pairs.get(j) is not None and pairs[j] < l_a:
                    next_b_pair = j
                    break
                    
            if next_a_pair is not None and next_b_pair is not None:
                if idx_a < next_a_pair:
                    col["top"] = seq_a[idx_a]
                    col["idx_a"] = idx_a
                    idx_a += 1
                elif idx_b > next_b_pair:
                    col["bottom"] = seq_b[idx_b - l_a]
                    col["idx_b"] = idx_b
                    idx_b -= 1
                else:
                    col["top"] = seq_a[idx_a]
                    col["idx_a"] = idx_a
                    idx_a += 1
            else:
                if idx_a < l_a:
                    col["top"] = seq_a[idx_a]
                    col["idx_a"] = idx_a
                    idx_a += 1
                elif idx_b >= l_a:
                    col["bottom"] = seq_b[idx_b - l_a]
                    col["idx_b"] = idx_b
                    idx_b -= 1
        columns.append(col)
        
    # 3. Annotation des rôles
    if extensible_strand == 'a':
        for col in columns:
            if col["idx_a"] == l_a - 1:
                col["role"] = "three_prime"
            elif col["idx_b"] is not None:
                partner_of_3p = pairs.get(l_a - 1)
                if partner_of_3p is not None and col["idx_b"] > partner_of_3p:
                    col["role"] = "template"
    elif extensible_strand == 'b':
        for col in columns:
            if col["idx_b"] == l_a + l_b - 1:
                col["role"] = "three_prime"
            elif col["idx_a"] is not None:
                partner_of_3p = pairs.get(l_a + l_b - 1)
                if partner_of_3p is not None and col["idx_a"] > partner_of_3p:
                    col["role"] = "template"
                    
    # 4. Compaction
    first_important = -1
    for i, col in enumerate(columns):
        if col["paired"] or col["role"]:
            first_important = i
            break
            
    if first_important > 3:
        new_cols = [{"top": "-", "bottom": "-", "paired": False, "role": "", "is_truncation": True}]
        new_cols.extend(columns[first_important - 2:])
        columns = new_cols
        
    return columns
