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
        # On extrait les 6 derniers nucléotides de l'amorce
        seq_primer = primer_a if is_primer_a else primer_b
        sub_primer = seq_primer[-6:]
        
        # On extrait la région correspondante sur le partenaire
        # On cherche les index min et max du partenaire pour ces 6 bases
        partner_indices = []
        start_idx = (l_a - 6) if is_primer_a else (l_a + l_b - 6)
        end_idx = idx_3p
        
        for i in range(max(0 if is_primer_a else l_a, start_idx), end_idx + 1):
            if i in pairs:
                p = pairs[i]
                # S'assurer que le partenaire est bien sur l'autre brin ou plus loin dans le même brin
                partner_indices.append(p)
                
        if not partner_indices:
            return False, 0.0
            
        min_p = min(partner_indices)
        max_p = max(partner_indices)
        
        # Extraction de la séquence partenaire
        concat_seq = primer_a + primer_b
        sub_partner = concat_seq[min_p:max_p+1]
        
        # Calcul du dG local avec RNA.cofold
        # RNA.cvar.temperature est supposé être déjà géré par le context manager dna_params
        # Mais on le set par précaution si on est appelé hors du context manager
        saved_temp = RNA.cvar.temperature
        RNA.cvar.temperature = temp_celsius
        
        try:
            # Évaluation locale
            _, local_dg = RNA.cofold(f"{sub_primer}&{sub_partner}")
        finally:
            RNA.cvar.temperature = saved_temp
            
        is_amp = local_dg <= enzyme.dimer_dg_threshold
        return is_amp, local_dg
        
    amp_a, dg_a = check_3p_end(True)
    amp_b, dg_b = check_3p_end(False)
    
    is_amplifiable = amp_a or amp_b
    min_dg = min(dg_a, dg_b) if (amp_a or amp_b) else 0.0
    
    return is_amplifiable, min_dg

