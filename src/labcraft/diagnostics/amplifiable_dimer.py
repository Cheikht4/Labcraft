import math
from typing import Any, Tuple, Dict, Optional
import RNA
from labcraft.diagnostics.enzyme import PolymeraseProfile
from labcraft.thermo.backends.base import DuplexEnergyBackend
from labcraft.thermo.salt import sodium_equivalent_for_folding
from labcraft.thermo.vienna import dna_params
from labcraft.buffer.monovalent import get_total_monovalent

from labcraft.thermo.constants import R_GAS_KCAL_MOL_K

_MAX_PF_3P_CACHE_SIZE = 16384
_SHARED_PF_3P_CACHE: Dict[Tuple[Any, ...], Dict[str, Any]] = {}

# Seuil heuristique non calibré (10 %) signalant dans le rapport HTML/JSON qu'un dimère classé
# non amplifiable au MFE possède néanmoins une sous-population Boltzmann 3' intermoléculaire notable.
# N'intervient JAMAIS dans `is_amplifiable`, `total_amp_dimer_conc_nm`, ni dans le verdict du panel.
# Uncalibrated heuristic threshold (10%) flagging in HTML/JSON reports that a dimer classified as
# non-amplifiable at MFE still carries a notable intermolecular 3' Boltzmann sub-population.
# NEVER enters `is_amplifiable`, `total_amp_dimer_conc_nm`, or the panel verdict.
SUBOPTIMAL_3P_PROB_THRESHOLD_HEURISTIC = 0.10

_WC_STRICT_PAIRS = frozenset({("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")})
_WC_OR_WOBBLE_PAIRS = frozenset({("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"), ("G", "T"), ("T", "G")})


def _is_valid_3p_stem_pair(
    b_3p: str,
    b_tpl: str,
    offset_from_3p: int,
    *,
    abs_window: int = 2,
    allow_penultimate_wobble: bool = False,
) -> bool:
    """Vérifie la compatibilité d'appariement et l'absence de veto ARMS à `offset_from_3p` (0 = base 3' terminale).
    Checks base-pair compatibility and absence of ARMS block at `offset_from_3p` (0 = strict 3' terminal base).
    """
    pair = (b_3p.upper(), b_tpl.upper())
    pos_1based = offset_from_3p + 1
    if pos_1based <= abs_window:
        if allow_penultimate_wobble and pos_1based == 2 and pair in (("G", "T"), ("T", "G")):
            return True
        return pair in _WC_STRICT_PAIRS
    return pair in _WC_OR_WOBBLE_PAIRS


def calc_boltzmann_3p_pairing(
    primer_a: str,
    primer_b: str,
    temp_celsius: float = 65.0,
    blocked_a: bool = False,
    blocked_b: bool = False,
    mon_molar: Optional[float] = None,
    *,
    min_3p_bp: int = 3,
    enzyme: Optional[PolymeraseProfile] = None,
) -> Dict[str, Any]:
    """Calcule les probabilités de Boltzmann d'appariement 3' en séparant strictement :
    1. La probabilité INTERMOLÉCULAIRE conditionnelle P(3' >= `min_3p_bp` pb contiguës avec matrice 5' | >= 1 paire A:B),
       évaluée analytiquement via `pf_dimer()` (`FAB` connecté) et fonctions de partition contraintes (`hc_add_bp`),
       en excluant les états dissociés sans paire entre A et B ainsi que tout auto-repliement intramoléculaire.
    2. La probabilité INTRAMOLÉCULAIRE `prob_3p_intra_{a,b}` d'épingle auto-amorçante sur chaque monomère isolé.

    Note : cet indicateur est purement informatif (rapport HTML et JSON) et n'entre pas dans le verdict.

    Computes Boltzmann 3'-end pairing probabilities while strictly separating:
    1. Conditional INTERMOLECULAR probability P(3' >= `min_3p_bp` contiguous bp with 5' template | >= 1 A:B pair),
       evaluated analytically via `pf_dimer()` (connected `FAB`) and hard-constrained partition functions (`hc_add_bp`),
       excluding unbound A&B states and intramolecular self-hairpins.
    2. Separate INTRAMOLECULAR self-priming hairpin probabilities `prob_3p_intra_{a,b}` on each isolated monomer.

    Note: this indicator is purely informational (HTML and JSON reports) and never enters the panel verdict.
    """
    if min_3p_bp < 2:
        raise ValueError("min_3p_bp doit être >= 2 pour exclure les contacts isolés d'une seule base en 3'.")

    l_a = len(primer_a)
    l_b = len(primer_b)
    if l_a == 0 or l_b == 0:
        return {
            "prob_3p_paired_a": 0.0,
            "prob_3p_paired_b": 0.0,
            "prob_3p_inter_a": 0.0,
            "prob_3p_inter_b": 0.0,
            "prob_3p_extensible_a": 0.0,
            "prob_3p_extensible_b": 0.0,
            "prob_3p_intra_a": 0.0,
            "prob_3p_intra_b": 0.0,
            "boltzmann_3p_prob": 0.0,
            "connected_dimer_dg_kcal": None,
            "ensemble_dg_kcal": 0.0,
            "min_3p_bp": min_3p_bp,
        }

    abs_win = enzyme.three_prime_absolute_window if enzyme is not None else 2
    allow_wob = bool(getattr(enzyme, "allow_penultimate_gt_wobble", False)) if enzyme is not None else False
    salt_key = round(mon_molar, 6) if mon_molar is not None else None
    cache_key: Tuple[Any, ...] = (
        primer_a,
        primer_b,
        round(temp_celsius, 3),
        blocked_a,
        blocked_b,
        salt_key,
        min_3p_bp,
        abs_win,
        allow_wob,
    )
    if cache_key in _SHARED_PF_3P_CACHE:
        return dict(_SHARED_PF_3P_CACHE[cache_key])

    rt = R_GAS_KCAL_MOL_K * (temp_celsius + 273.15)
    hc_flags = RNA.CONSTRAINT_CONTEXT_ALL_LOOPS | RNA.CONSTRAINT_CONTEXT_ENFORCE

    def _inter_conditional_3p(s_ext: str, s_tpl: str) -> Tuple[float, float, float]:
        """Retourne (P_inter(3' de s_ext | s_ext:s_tpl connecté), FAB_conn, F0AB_total) dans l'ordre `s_ext&s_tpl`.
        Placer `s_ext` juste avant '&' garantit que l'hélice 3' `(l1 - r, k + r)` est sans pseudonœud
        vis-à-vis de toute tige intermoléculaire additionnelle en amont.
        """
        l1, l2 = len(s_ext), len(s_tpl)
        fc = RNA.fold_compound(f"{s_ext}&{s_tpl}")
        _st, _fa, _fb, fab_conn_raw, f0ab_raw = fc.pf_dimer()
        fab_conn = float(fab_conn_raw)
        f0ab = float(f0ab_raw)
        if fab_conn >= 900.0 or l1 < min_3p_bp or l2 < min_3p_bp + 1 or rt <= 0.0:
            return 0.0, fab_conn, f0ab

        prob = 0.0
        # k (1-based dans s_ext&s_tpl) apparié à la base 3' terminale `l1` de s_ext :
        # k >= l1 + 2 impose au moins 1 nt de matrice 5' sur s_tpl (dont le 5' est à l1 + 1).
        for k in range(l1 + 2, l1 + l2 - min_3p_bp + 2):
            if not all(
                _is_valid_3p_stem_pair(
                    s_ext[l1 - 1 - r],
                    s_tpl[(k - l1 - 1) + r],
                    r,
                    abs_window=abs_win,
                    allow_penultimate_wobble=allow_wob,
                )
                for r in range(min_3p_bp)
            ):
                continue
            fc_k = RNA.fold_compound(f"{s_ext}&{s_tpl}")
            for r in range(min_3p_bp):
                fc_k.hc_add_bp(l1 - r, k + r, hc_flags)
            _stk, _fak, _fbk, fab_k_raw, _f0k = fc_k.pf_dimer()
            fab_k = float(fab_k_raw)
            if fab_k < 900.0:
                prob += math.exp(-(fab_k - fab_conn) / rt)
        return min(1.0, max(0.0, prob)), fab_conn, f0ab

    def _intra_monomer_3p(seq_m: str) -> float:
        """Retourne P_intra(3' >= min_3p_bp pb contiguës avec matrice 5') sur le monomère isolé `seq_m`."""
        lm = len(seq_m)
        # Boucle d'épingle interne >= 3 nt entre `k + min_3p_bp - 1` et `lm - min_3p_bp + 1`, et `k >= 2`
        if lm < 2 * min_3p_bp + 4 or rt <= 0.0:
            return 0.0
        fc_m = RNA.fold_compound(seq_m)
        _, gm_raw = fc_m.pf()
        gm = float(gm_raw)
        prob = 0.0
        for k in range(2, lm - 2 * min_3p_bp):
            if not all(
                _is_valid_3p_stem_pair(
                    seq_m[lm - 1 - r],
                    seq_m[(k - 1) + r],
                    r,
                    abs_window=abs_win,
                    allow_penultimate_wobble=allow_wob,
                )
                for r in range(min_3p_bp)
            ):
                continue
            fc_mk = RNA.fold_compound(seq_m)
            for r in range(min_3p_bp):
                fc_mk.hc_add_bp(k + r, lm - r, hc_flags)
            _, gmk_raw = fc_mk.pf()
            gmk = float(gmk_raw)
            if gmk < 900.0:
                prob += math.exp(-(gmk - gm) / rt)
        return min(1.0, max(0.0, prob))

    with dna_params(temp_celsius=temp_celsius, mon_molar=mon_molar):
        prob_inter_a, fab_ab, f0_ab = _inter_conditional_3p(primer_a, primer_b)
        if primer_a == primer_b:
            prob_inter_b, fab_ba, f0_ba = prob_inter_a, fab_ab, f0_ab
        else:
            prob_inter_b, fab_ba, f0_ba = _inter_conditional_3p(primer_b, primer_a)
        prob_intra_a = _intra_monomer_3p(primer_a)
        prob_intra_b = prob_intra_a if primer_a == primer_b else _intra_monomer_3p(primer_b)

    prob_ext_a = 0.0 if blocked_a else prob_inter_a
    prob_ext_b = 0.0 if blocked_b else prob_inter_b
    fab_best = min(fab_ab, fab_ba)
    f0_best = min(f0_ab, f0_ba)

    result: Dict[str, Any] = {
        "prob_3p_paired_a": prob_inter_a,
        "prob_3p_paired_b": prob_inter_b,
        "prob_3p_inter_a": prob_inter_a,
        "prob_3p_inter_b": prob_inter_b,
        "prob_3p_extensible_a": prob_ext_a,
        "prob_3p_extensible_b": prob_ext_b,
        "prob_3p_intra_a": prob_intra_a,
        "prob_3p_intra_b": prob_intra_b,
        "boltzmann_3p_prob": max(prob_ext_a, prob_ext_b),
        "connected_dimer_dg_kcal": fab_best if fab_best < 900.0 else None,
        "ensemble_dg_kcal": f0_best,
        "min_3p_bp": min_3p_bp,
    }
    if len(_SHARED_PF_3P_CACHE) >= _MAX_PF_3P_CACHE_SIZE:
        _SHARED_PF_3P_CACHE.clear()
    _SHARED_PF_3P_CACHE[cache_key] = result
    return dict(result)

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
    temp_celsius: float = 65.0,
    blocked_a: bool = False,
    blocked_b: bool = False,
    lna_positions: tuple[int, ...] = (),
    na_mm: Optional[float] = None,
    k_mm: float = 0.0,
    tris_mm: float = 0.0,
    mg_mm: float = 0.0,
    dntp_mm: float = 0.0
) -> Tuple[bool, float, Optional[str], bool]:
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
        (is_amplifiable, min_dg_3p, extensible_strand, is_blocked_veto)
    """
    l_a = len(primer_a)
    l_b = len(primer_b)
    
    if len(structure) != l_a + l_b:
        raise ValueError(f"La longueur de la structure ({len(structure)}) ne correspond pas à la somme des amorces ({l_a + l_b})")
        
    pairs = get_pairs(structure)
    
    def check_3p_end(is_primer_a: bool) -> Tuple[bool, float, bool]:
        """Vérifie si l'extrémité 3' d'une amorce est extensible."""
        # Index du 3' dans la structure concaténée
        idx_3p = (l_a - 1) if is_primer_a else (l_a + l_b - 1)
        
        # 1. Le 3' est-il apparié ?
        if idx_3p not in pairs:
            return False, 0.0, False
            
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
            return False, 0.0, False
            
        # 3. Calcul de l'énergie locale du 3'
        # Parcourt jusqu'à 6 paires contiguës depuis le 3' end
        paired_indices: list[int] = []
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
        
        from labcraft.thermo.lna import load_lna_params, _LNA_PARAMS_MXL, _LNA_PARAMS_XLN
        load_lna_params()
        
        for i in range(len(segment_seq) - 1):
            dinuc = segment_seq[i:i+2]
            key = _SEQ_TO_NN_KEY.get(dinuc)
            if key:
                dh, ds = _NN_PARAMS[key]
                dh_total += dh
                ds_total += ds
                
            # LNA correction for this step
            actual_idx = paired_indices[i]
            next_actual_idx = paired_indices[i+1]
            # Since paired_indices are contiguous 5' to 3' along the strand, 
            # actual_idx and next_actual_idx are adjacent indices in the concatenated string.
            if actual_idx in lna_positions:
                key_xln = (segment_seq[i], segment_seq[i+1])
                if key_xln in _LNA_PARAMS_XLN:
                    dh_total += _LNA_PARAMS_XLN[key_xln][0]
                    ds_total += _LNA_PARAMS_XLN[key_xln][1]
            if next_actual_idx in lna_positions:
                key_mxl = (segment_seq[i], segment_seq[i+1])
                if key_mxl in _LNA_PARAMS_MXL:
                    dh_total += _LNA_PARAMS_MXL[key_mxl][0]
                    ds_total += _LNA_PARAMS_MXL[key_mxl][1]
                
        # Correction saline sur le segment 3' court via l'équivalent sodium (inversion Owczarzy)
        # et la correction entropique de SantaLucia 1998 : dS += 0.368 * (N - 1) * ln[Na_eq].
        # Salt correction on the short 3' segment using sodium equivalent (Owczarzy inversion)
        # and SantaLucia 1998 entropic correction: dS += 0.368 * (N - 1) * ln[Na_eq].
        n_bp = len(segment_seq)
        if n_bp > 1 and (na_mm is not None or k_mm > 0.0 or tris_mm > 0.0 or mg_mm > 0.0):
            mon_molar = get_total_monovalent(
                na_mm if na_mm is not None else 0.0, k_mm, tris_mm
            ) / 1000.0
            mg_molar = mg_mm / 1000.0
            dntp_molar = (dntp_mm or 0.0) / 1000.0
            na_eq = sodium_equivalent_for_folding(mon_molar, mg_molar, dntp_molar)
            if na_eq > 0.0:
                ds_total += 0.368 * (n_bp - 1) * math.log(na_eq)

        # dG = dH - T*dS
        temp_k = temp_celsius + 273.15
        local_dg = dh_total - temp_k * (ds_total / 1000.0)
        
        # Application du veto 3'
        # Le veto pour les dimères n'est absolu que sur les positions critiques (1-2 selon ARMS)
        window = enzyme.three_prime_absolute_window
        passed_veto = True
        for k in range(window):
            idx_k = idx_3p - k
            # Dans la concaténation 5'->3' (A + B), reculer depuis le 3' (idx_3p - k)
            # impose d'avancer en 5'->3' sur le partenaire antiparallèle (partner_idx + k) pour A comme pour B.
            # In the 5'->3' concatenation (A + B), stepping upstream from 3' (idx_3p - k)
            # requires stepping 5'->3' along the antiparallel partner (partner_idx + k) for both A and B.
            expected_partner_k = partner_idx + k
            
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
        is_blocked = False
        
        if is_amp:
            # Application du veto de blocage 3'
            if is_primer_a and blocked_a:
                is_amp = False
                is_blocked = True
            elif not is_primer_a and blocked_b:
                is_amp = False
                is_blocked = True
                
        return is_amp, local_dg, is_blocked
        
    amp_a, dg_a, blocked_veto_a = check_3p_end(True)
    amp_b, dg_b, blocked_veto_b = check_3p_end(False)
    
    is_amplifiable = amp_a or amp_b
    
    # Pour le dG_3p, on privilégie les brins effectivement amplifiables s'il y en a,
    # sinon la pire (plus négative) énergie parmi celles qui ont une matrice extensible.
    # For dG_3p, prefer amplifiable strands if any, otherwise the worst extensible template.
    amp_dgs = []
    if amp_a:
        amp_dgs.append((dg_a, 'a'))
    if amp_b:
        amp_dgs.append((dg_b, 'b'))

    valid_dgs = []
    if dg_a != 0.0:
        valid_dgs.append((dg_a, 'a'))
    if dg_b != 0.0:
        valid_dgs.append((dg_b, 'b'))
        
    candidate_dgs = amp_dgs if amp_dgs else valid_dgs
    if candidate_dgs:
        best_tuple = min(candidate_dgs, key=lambda x: x[0])
        min_dg = best_tuple[0]
        ext_strand = best_tuple[1]
    else:
        min_dg = 0.0
        ext_strand = None
    
    return is_amplifiable, min_dg, ext_strand, (blocked_veto_a or blocked_veto_b)

from labcraft.lamp.domains import PhysicalPrimer

def enrich_dimer_boltzmann_3p(
    details: Dict[str, Any],
    p_a: PhysicalPrimer,
    p_b: PhysicalPrimer,
    is_amplifiable: bool,
    temp_celsius: float,
    *,
    suboptimal_3p_prob_threshold: float = SUBOPTIMAL_3P_PROB_THRESHOLD_HEURISTIC,
    min_3p_bp: int = 3,
    enzyme: Optional[PolymeraseProfile] = None,
    **backend_kwargs: Any,
) -> Dict[str, Any]:
    """Enrichit les métadonnées de rapport d'un dimère avec les probabilités de Boltzmann 3'
    intermoléculaires conditionnelles et intramoléculaires séparées (sans modifier le verdict).
    Le seuil `suboptimal_3p_prob_threshold = 0.10` (10 %) est une heuristique non calibrée.

    Enriches dimer report metadata with conditional intermolecular and separated intramolecular
    Boltzmann 3' probabilities (never alters the verdict). The `0.10` (10%) cutoff for
    `suboptimal_3p_masked` is an uncalibrated heuristic.
    """
    na_mm = backend_kwargs.get("na_mm")
    k_mm = float(backend_kwargs.get("k_mm", 0.0) or 0.0)
    tris_mm = float(backend_kwargs.get("tris_mm", 0.0) or 0.0)
    mg_mm = float(backend_kwargs.get("mg_mm", 0.0) or 0.0)
    dntp_mm = float(backend_kwargs.get("dntp_mm", 0.0) or 0.0)

    mon_molar: Optional[float] = None
    if na_mm is not None or k_mm > 0.0 or tris_mm > 0.0 or mg_mm > 0.0:
        mon_true = get_total_monovalent(float(na_mm or 0.0), k_mm, tris_mm) / 1000.0
        mon_molar = sodium_equivalent_for_folding(mon_true, mg_mm / 1000.0, dntp_mm / 1000.0)

    order = details.get("order", "a&b")
    if order == "b&a":
        s1, s2 = p_b.sequence, p_a.sequence
        blk1, blk2 = p_b.blocked_3prime, p_a.blocked_3prime
    else:
        s1, s2 = p_a.sequence, p_b.sequence
        blk1, blk2 = p_a.blocked_3prime, p_b.blocked_3prime

    pf_info = calc_boltzmann_3p_pairing(
        s1,
        s2,
        temp_celsius=temp_celsius,
        blocked_a=blk1,
        blocked_b=blk2,
        mon_molar=mon_molar,
        min_3p_bp=min_3p_bp,
        enzyme=enzyme,
    )
    if order == "b&a":
        prob_a = pf_info["prob_3p_extensible_b"]
        prob_b = pf_info["prob_3p_extensible_a"]
        intra_a = pf_info["prob_3p_intra_b"]
        intra_b = pf_info["prob_3p_intra_a"]
    else:
        prob_a = pf_info["prob_3p_extensible_a"]
        prob_b = pf_info["prob_3p_extensible_b"]
        intra_a = pf_info["prob_3p_intra_a"]
        intra_b = pf_info["prob_3p_intra_b"]

    max_prob = float(pf_info["boltzmann_3p_prob"])
    details["prob_3p_extensible_a"] = float(prob_a)
    details["prob_3p_extensible_b"] = float(prob_b)
    details["prob_3p_intra_a"] = float(intra_a)
    details["prob_3p_intra_b"] = float(intra_b)
    details["boltzmann_3p_prob"] = max_prob
    details["connected_dimer_dg_kcal"] = pf_info["connected_dimer_dg_kcal"]
    details["ensemble_dg_kcal"] = float(pf_info["ensemble_dg_kcal"])
    # Heuristique non calibrée (seuil 10 %) : informatif uniquement, n'entre pas dans le verdict
    # Uncalibrated heuristic (10% threshold): informational only, does not enter the verdict
    details["suboptimal_3p_masked"] = bool(
        (not is_amplifiable) and (max_prob >= suboptimal_3p_prob_threshold)
    )
    return details


def evaluate_pair_amplifiable(
    p_a: PhysicalPrimer, 
    p_b: PhysicalPrimer, 
    backend: DuplexEnergyBackend, 
    enzyme: PolymeraseProfile, 
    temp_celsius: float,
    *,
    include_boltzmann_3p: bool = False,
    suboptimal_3p_prob_threshold: float = 0.10,
    **backend_kwargs: Any
) -> Tuple[bool, float, Dict[str, Any]]:
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
    seq_a = p_a.sequence
    seq_b = p_b.sequence
    blocked_a = p_a.blocked_3prime
    blocked_b = p_b.blocked_3prime
    
    lna_a = p_a.lna_positions
    lna_b = p_b.lna_positions

    salt_kwargs = {
        "na_mm": backend_kwargs.get("na_mm"),
        "k_mm": backend_kwargs.get("k_mm", 0.0),
        "tris_mm": backend_kwargs.get("tris_mm", 0.0),
        "mg_mm": backend_kwargs.get("mg_mm", 0.0),
        "dntp_mm": backend_kwargs.get("dntp_mm", 0.0),
    }

    pa_parent = p_a.parent_name or p_a.name.split("#")[0]
    pb_parent = p_b.parent_name or p_b.name.split("#")[0]

    def _resolve_ext_and_template(order_str: str, ext_strand: Optional[str]) -> Tuple[str, str]:
        """Résout les noms d'amorces parentes étendue et matrice selon l'ordre de repliement.
        Resolves extended and template parent primer names accounting for fold order.
        """
        if order_str == "a&a":
            return (pa_parent, pa_parent) if ext_strand else ("", "")
        strand_a_parent = pa_parent if order_str == "a&b" else pb_parent
        strand_b_parent = pb_parent if order_str == "a&b" else pa_parent
        if ext_strand == "a":
            return strand_a_parent, strand_b_parent
        if ext_strand == "b":
            return strand_b_parent, strand_a_parent
        return "", ""

    if seq_a == seq_b:
        # Homodimère: un seul ordre suffit car la concaténation est symétrique
        mapped_lna = tuple(list(lna_a) + [pos + len(seq_a) for pos in lna_a])
        backend_kwargs['lna_positions'] = lna_a
        res = backend.calc_homodimer(seq_a, temp_celsius=temp_celsius, **backend_kwargs)
        struct, mfe = res.structure.replace('&', ''), res.dg_kcal
        is_amp, dg_3p, ext_strand, blocked_veto = is_amplifiable_dimer(seq_a, seq_a, struct, mfe, enzyme, temp_celsius, blocked_a, blocked_b, lna_positions=mapped_lna, **salt_kwargs)
        ext_par, tpl_par = _resolve_ext_and_template("a&a", ext_strand)
        details = {
            "seq_a": seq_a, "seq_b": seq_a, "structure": res.structure, 
            "delta_g": mfe, "delta_g_3p": dg_3p,
            "order": "a&a",
            "extensible_strand": ext_strand,
            "parent_a": pa_parent,
            "parent_b": pb_parent,
            "extensible_primer_parent": ext_par,
            "template_primer_parent": tpl_par,
            "is_blocked_veto": blocked_veto
        }
        if include_boltzmann_3p:
            enrich_dimer_boltzmann_3p(
                details,
                p_a,
                p_b,
                is_amp,
                temp_celsius,
                suboptimal_3p_prob_threshold=suboptimal_3p_prob_threshold,
                enzyme=enzyme,
                **salt_kwargs,
            )
        return is_amp, dg_3p, details
        
    # Hétérodimère: on évalue a&b ET b&a
    backend_kwargs['lna_positions_a'] = lna_a
    backend_kwargs['lna_positions_b'] = lna_b
    
    res_ab = backend.calc_heterodimer(seq_a, seq_b, temp_celsius=temp_celsius, **backend_kwargs)
    struct_ab, mfe_ab = res_ab.structure.replace('&', ''), res_ab.dg_kcal
    mapped_lna_ab = tuple(list(lna_a) + [pos + len(seq_a) for pos in lna_b])
    is_amp_ab, dg_3p_ab, ext_strand_ab, blocked_veto_ab = is_amplifiable_dimer(seq_a, seq_b, struct_ab, mfe_ab, enzyme, temp_celsius, blocked_a, blocked_b, lna_positions=mapped_lna_ab, **salt_kwargs)
    
    # Pour b&a, on intervertit les rôles a et b
    backend_kwargs['lna_positions_a'] = lna_b
    backend_kwargs['lna_positions_b'] = lna_a
    res_ba = backend.calc_heterodimer(seq_b, seq_a, temp_celsius=temp_celsius, **backend_kwargs)
    struct_ba, mfe_ba = res_ba.structure.replace('&', ''), res_ba.dg_kcal
    mapped_lna_ba = tuple(list(lna_b) + [pos + len(seq_b) for pos in lna_a])
    is_amp_ba, dg_3p_ba, ext_strand_ba, blocked_veto_ba = is_amplifiable_dimer(seq_b, seq_a, struct_ba, mfe_ba, enzyme, temp_celsius, blocked_b, blocked_a, lna_positions=mapped_lna_ba, **salt_kwargs)
    
    is_amplifiable = is_amp_ab or is_amp_ba
    
    # Pour le dg_3p, on prend le pire (le plus négatif) parmi ceux qui sont extensibles,
    # ou 0.0 si aucun n'a un template. 
    # S'il y a un is_amplifiable = True, on prend le dG du is_amplifiable qui l'est.
    if is_amp_ab and not is_amp_ba:
        worst_dg = dg_3p_ab
        chosen_res = res_ab
        chosen_order = "a&b"
        seq_1, seq_2 = seq_a, seq_b
        worst_ext_strand = ext_strand_ab
        worst_blocked = blocked_veto_ab
    elif is_amp_ba and not is_amp_ab:
        worst_dg = dg_3p_ba
        chosen_res = res_ba
        chosen_order = "b&a"
        seq_1, seq_2 = seq_b, seq_a
        worst_ext_strand = ext_strand_ba
        worst_blocked = blocked_veto_ba
    elif is_amp_ab and is_amp_ba:
        if dg_3p_ab <= dg_3p_ba: # Plus négatif = pire
            worst_dg = dg_3p_ab
            chosen_res = res_ab
            chosen_order = "a&b"
            seq_1, seq_2 = seq_a, seq_b
            worst_ext_strand = ext_strand_ab
            worst_blocked = blocked_veto_ab
        else:
            worst_dg = dg_3p_ba
            chosen_res = res_ba
            chosen_order = "b&a"
            seq_1, seq_2 = seq_b, seq_a
            worst_ext_strand = ext_strand_ba
            worst_blocked = blocked_veto_ba
    else:
        # Aucun amplifiable : on prend la structure qui a le 3' le plus favorable (même si > seuil),
        # ou à défaut la structure la plus stable globalement.
        # Attention : si un dimère n'est pas amplifiable UNIQUEMENT à cause du veto,
        # on veut le remonter pour l'afficher !
        if blocked_veto_ab and not blocked_veto_ba:
            worst_dg = dg_3p_ab
            chosen_res = res_ab
            chosen_order = "a&b"
            seq_1, seq_2 = seq_a, seq_b
            worst_ext_strand = ext_strand_ab
            worst_blocked = True
        elif blocked_veto_ba and not blocked_veto_ab:
            worst_dg = dg_3p_ba
            chosen_res = res_ba
            chosen_order = "b&a"
            seq_1, seq_2 = seq_b, seq_a
            worst_ext_strand = ext_strand_ba
            worst_blocked = True
        elif dg_3p_ab <= dg_3p_ba and dg_3p_ab < 0.0:
            worst_dg = dg_3p_ab
            chosen_res = res_ab
            chosen_order = "a&b"
            seq_1, seq_2 = seq_a, seq_b
            worst_ext_strand = ext_strand_ab
            worst_blocked = blocked_veto_ab
        elif dg_3p_ba < 0.0:
            worst_dg = dg_3p_ba
            chosen_res = res_ba
            chosen_order = "b&a"
            seq_1, seq_2 = seq_b, seq_a
            worst_ext_strand = ext_strand_ba
            worst_blocked = blocked_veto_ba
        elif mfe_ab <= mfe_ba:
            worst_dg = dg_3p_ab
            chosen_res = res_ab
            chosen_order = "a&b"
            seq_1, seq_2 = seq_a, seq_b
            worst_ext_strand = ext_strand_ab
            worst_blocked = blocked_veto_ab
        else:
            worst_dg = dg_3p_ba
            chosen_res = res_ba
            chosen_order = "b&a"
            seq_1, seq_2 = seq_b, seq_a
            worst_ext_strand = ext_strand_ba
            worst_blocked = blocked_veto_ba
            
    ext_par, tpl_par = _resolve_ext_and_template(chosen_order, worst_ext_strand)
    details = {
        "seq_a": seq_1,
        "seq_b": seq_2,
        "structure": chosen_res.structure,
        "delta_g": chosen_res.dg_kcal,
        "delta_g_3p": worst_dg,
        "order": chosen_order,
        "extensible_strand": worst_ext_strand,
        "parent_a": pa_parent,
        "parent_b": pb_parent,
        "extensible_primer_parent": ext_par,
        "template_primer_parent": tpl_par,
        "is_blocked_veto": worst_blocked
    }
    if include_boltzmann_3p:
        enrich_dimer_boltzmann_3p(
            details,
            p_a,
            p_b,
            is_amplifiable,
            temp_celsius,
            suboptimal_3p_prob_threshold=suboptimal_3p_prob_threshold,
            enzyme=enzyme,
            **salt_kwargs,
        )
    return is_amplifiable, worst_dg, details
