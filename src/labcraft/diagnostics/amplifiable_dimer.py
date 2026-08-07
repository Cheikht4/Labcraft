import RNA
from typing import Tuple, List, Dict
from labcraft.diagnostics.enzyme import PolymeraseProfile
from labcraft.thermo.backends.base import DuplexEnergyBackend

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
        
        # Application du veto 3'
        # Le veto pour les dimères n'est absolu que sur les positions critiques (1-2 selon ARMS)
        window = enzyme.three_prime_absolute_window
        passed_veto = True
        for k in range(window):
            idx_k = idx_3p - k
            expected_partner_k = partner_idx + k if is_primer_a else partner_idx - k
            
            # Limites
            if is_primer_a and idx_k < 0:
                passed_veto = False
                break
            if not is_primer_a and idx_k < l_a:
                passed_veto = False
                break
                
            if idx_k not in pairs or pairs[idx_k] != expected_partner_k:
                passed_veto = False
                break
                
        is_amp = (local_dg <= enzyme.dimer_dg_threshold) and passed_veto
        return is_amp, local_dg
        
    amp_a, dg_a = check_3p_end(True)
    amp_b, dg_b = check_3p_end(False)
    
    is_amplifiable = amp_a or amp_b
    
    # Pour le dG_3p, on veut remonter la pire (plus négative) énergie parmi celles 
    # qui ont une template extensible, même si ce n'est pas "amplifiable" (< seuil).
    valid_dgs = []
    if dg_a != 0.0:
        valid_dgs.append(dg_a)
    if dg_b != 0.0:
        valid_dgs.append(dg_b)
        
    min_dg = min(valid_dgs) if valid_dgs else 0.0
    
    return is_amplifiable, min_dg

def evaluate_pair_amplifiable(
    seq_a: str, 
    seq_b: str, 
    backend: DuplexEnergyBackend, 
    enzyme: PolymeraseProfile, 
    temp_celsius: float,
    **backend_kwargs
) -> Tuple[bool, float, Dict]:
    """Évalue un couple d'amorces dans les deux ordres de concaténation et
    renvoie le verdict combiné. Un dimère est amplifiable si l'un OU l'autre
    ordre révèle un 3' extensible.
    
    Note: tester deux ordres n'échantillonne que deux structures de l'ensemble dégénéré.
    La solution vraiment complète utiliserait les probabilités d'appariement de la 
    fonction de partition (RNA.pf_fold / co-fold partition function) pour juger si 
    l'extrémité 3' s'apparie avec une probabilité supérieure à un seuil, indépendamment 
    d'une structure MFE arbitraire.
    
    Returns:
        (is_amplifiable, worst_dg_3p, details)
        où details indique l'ordre et l'extrémité 3' responsables.
    """
    if seq_a == seq_b:
        # Homodimère: un seul ordre suffit car la concaténation est symétrique
        res = backend.calc_homodimer(seq_a, temp_celsius=temp_celsius, **backend_kwargs)
        struct, mfe = res.structure.replace('&', ''), res.dg_kcal
        is_amp, dg_3p = is_amplifiable_dimer(seq_a, seq_a, struct, mfe, enzyme, temp_celsius)
        details = {
            "seq_a": seq_a, "seq_b": seq_a, "structure": res.structure, 
            "delta_g": mfe, "delta_g_3p": dg_3p,
            "order": "a&a"
        }
        return is_amp, dg_3p, details
        
    # Hétérodimère: on évalue a&b ET b&a
    res_ab = backend.calc_heterodimer(seq_a, seq_b, temp_celsius=temp_celsius, **backend_kwargs)
    struct_ab, mfe_ab = res_ab.structure.replace('&', ''), res_ab.dg_kcal
    is_amp_ab, dg_3p_ab = is_amplifiable_dimer(seq_a, seq_b, struct_ab, mfe_ab, enzyme, temp_celsius)
    
    res_ba = backend.calc_heterodimer(seq_b, seq_a, temp_celsius=temp_celsius, **backend_kwargs)
    struct_ba, mfe_ba = res_ba.structure.replace('&', ''), res_ba.dg_kcal
    is_amp_ba, dg_3p_ba = is_amplifiable_dimer(seq_b, seq_a, struct_ba, mfe_ba, enzyme, temp_celsius)
    
    is_amplifiable = is_amp_ab or is_amp_ba
    
    # Pour le dg_3p, on prend le pire (le plus négatif) parmi ceux qui sont extensibles,
    # ou 0.0 si aucun n'a un template. 
    # S'il y a un is_amplifiable = True, on prend le dG du is_amplifiable qui l'est.
    if is_amp_ab and not is_amp_ba:
        worst_dg = dg_3p_ab
        chosen_res = res_ab
        chosen_order = "a&b"
        seq_1, seq_2 = seq_a, seq_b
    elif is_amp_ba and not is_amp_ab:
        worst_dg = dg_3p_ba
        chosen_res = res_ba
        chosen_order = "b&a"
        seq_1, seq_2 = seq_b, seq_a
    elif is_amp_ab and is_amp_ba:
        if dg_3p_ab <= dg_3p_ba: # Plus négatif = pire
            worst_dg = dg_3p_ab
            chosen_res = res_ab
            chosen_order = "a&b"
            seq_1, seq_2 = seq_a, seq_b
        else:
            worst_dg = dg_3p_ba
            chosen_res = res_ba
            chosen_order = "b&a"
            seq_1, seq_2 = seq_b, seq_a
    else:
        # Aucun amplifiable, on prend le pire dG qui a quand même une template, s'il y en a.
        valid_dgs = []
        if dg_3p_ab != 0.0: valid_dgs.append((dg_3p_ab, res_ab, "a&b", seq_a, seq_b))
        if dg_3p_ba != 0.0: valid_dgs.append((dg_3p_ba, res_ba, "b&a", seq_b, seq_a))
        
        if valid_dgs:
            best_tuple = min(valid_dgs, key=lambda x: x[0])
            worst_dg = best_tuple[0]
            chosen_res = best_tuple[1]
            chosen_order = best_tuple[2]
            seq_1, seq_2 = best_tuple[3], best_tuple[4]
        else:
            worst_dg = 0.0
            chosen_res = res_ab
            chosen_order = "a&b"
            seq_1, seq_2 = seq_a, seq_b
            
    details = {
        "seq_a": seq_1, "seq_b": seq_2, "structure": chosen_res.structure, 
        "delta_g": chosen_res.dg_kcal, "delta_g_3p": worst_dg,
        "order": chosen_order
    }
    return is_amplifiable, worst_dg, details

