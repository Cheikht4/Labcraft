import RNA
from typing import Tuple, List, Dict
from labcraft.diagnostics.enzyme import PolymeraseProfile

def get_pairs(structure: str) -> Dict[int, int]:
    """Parse la structure dot-bracket et retourne un dictionnaire de paires."""
    stack = []
    pairs = {}
    for i, c in enumerate(structure):
        if c == '(':
            stack.append(i)
        elif c == ')':
            if stack:
                j = stack.pop()
                pairs[i] = j
                pairs[j] = i
    return pairs

def is_amplifiable_dimer(
    primer_a: str, 
    primer_b: str, 
    structure: str, 
    mfe: float, 
    enzyme: PolymeraseProfile,
    temp_celsius: float = 65.0
) -> Tuple[bool, float]:
    """
    Détermine si un dimère (A + B) est amplifiable.
    Retourne (amplifiable, dg_3p).
    
    Args:
        primer_a: Séquence de la première amorce (5' -> 3')
        primer_b: Séquence de la deuxième amorce (5' -> 3', ou identique à A pour homodimère)
        structure: Structure dot-bracket retournée par cofold (sans '&')
        mfe: Energie libre globale du dimère
        enzyme: Profil de la polymérase
        temp_celsius: Température
        
    Returns:
        (is_amplifiable, min_dg_3p)
    """
    l_a = len(primer_a)
    l_b = len(primer_b)
    
    if len(structure) != l_a + l_b:
        raise ValueError(f"La longueur de la structure ({len(structure)}) ne correspond pas à la somme des amorces ({l_a + l_b})")
        
    pairs = get_pairs(structure)
    
    def check_3p_end(is_primer_a: bool) -> Tuple[bool, float]:
        """Vérifie si l'extrémité 3' d'une amorce est extensible."""
        # Index du 3' dans la structure concaténée
        idx_3p = (l_a - 1) if is_primer_a else (l_a + l_b - 1)
        
        # 1. Le 3' est-il apparié ?
        if idx_3p not in pairs:
            return False, 0.0
            
        partner_idx = pairs[idx_3p]
        
        # 2. Y a-t-il une matrice à copier ? (Le partenaire dépasse-t-il côté 5' ?)
        if is_primer_a:
            # A s'étend sur le partenaire.
            # Si le partenaire est B (partner_idx >= l_a), le 5' de B est à l'index l_a.
            # A s'étend vers le 5' de B. Il faut que partner_idx > l_a.
            if partner_idx >= l_a:
                has_template = (partner_idx > l_a)
            else:
                # Auto-complémentarité terminale sur A lui-même.
                # A s'étend vers le 5' de A (index 0).
                has_template = (partner_idx > 0)
        else:
            # B s'étend sur le partenaire.
            # Si le partenaire est A (partner_idx < l_a), le 5' de A est à l'index 0.
            # B s'étend vers le 5' de A. Il faut que partner_idx > 0.
            if partner_idx < l_a:
                has_template = (partner_idx > 0)
            else:
                # Auto-complémentarité terminale sur B lui-même.
                # Le 5' de B est à l'index l_a.
                has_template = (partner_idx > l_a)
                
        if not has_template:
            return False, 0.0
            
        # 3. Calcul de l'énergie locale du 3'
        # Parcourt jusqu'à 6 paires contiguës depuis le 3' end
        paired_indices = []
        c_idx = idx_3p
        c_partner = partner_idx
        
        while len(paired_indices) < 6:
            paired_indices.append(c_idx)
            next_idx = c_idx - 1
            expected_partner = c_partner + 1
            
            # Limite 5' de l'amorce courante
            if is_primer_a and next_idx < 0:
                break
            if not is_primer_a and next_idx < l_a:
                break
                
            # Vérifie la contiguïté
            if next_idx not in pairs or pairs[next_idx] != expected_partner:
                break
                
            c_idx = next_idx
            c_partner = expected_partner
            
        paired_indices.reverse() # Remettre en 5' -> 3'
        
        seq_concat = primer_a + primer_b
        segment_seq = "".join(seq_concat[i] for i in paired_indices)
        
        # Calcul NN SantaLucia 1998 (sans pénalités d'initiation)
        dh_total = 0.0
        ds_total = 0.0
        
        from labcraft.thermo.backends.native import _NN_PARAMS, _SEQ_TO_NN_KEY
        
        for i in range(len(segment_seq) - 1):
            dinuc = segment_seq[i:i+2]
            key = _SEQ_TO_NN_KEY.get(dinuc)
            if key:
                dh, ds = _NN_PARAMS[key]
                dh_total += dh
                ds_total += ds
                
        # dG = dH - T*dS
        temp_k = temp_celsius + 273.15
        local_dg = dh_total - temp_k * (ds_total / 1000.0)
        
        is_amp = local_dg <= enzyme.dimer_dg_threshold
        return is_amp, local_dg
        
    amp_a, dg_a = check_3p_end(True)
    amp_b, dg_b = check_3p_end(False)
    
    is_amplifiable = amp_a or amp_b
    min_dg = min(dg_a, dg_b) if (amp_a or amp_b) else 0.0
    
    return is_amplifiable, min_dg

