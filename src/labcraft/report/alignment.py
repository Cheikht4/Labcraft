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
