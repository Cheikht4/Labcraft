"""Thermodynamique des mésappariements et règle du 3'.
Mismatch thermodynamics and 3' enzymatic extension rules.

Fournit le calcul d'énergie pour des duplexes avec mésappariements
(internes isolés, terminaux et boucles internes symétriques consécutives)
et l'heuristique d'extensibilité du 3'.
Provides energy calculation for duplexes with mismatches (single internal,
terminal, and consecutive symmetric internal loops) and 3' extensibility heuristics.
"""
import csv
import math
import os
from typing import Tuple, Optional, Dict
from labcraft.diagnostics.enzyme import PolymeraseProfile
from labcraft.thermo.backends.native import _NN_PARAMS, _SEQ_TO_NN_KEY, _INIT_AT, _INIT_GC
from labcraft.thermo.constants import R_GAS_KCAL_MOL_K

# (dH, dS) par clé / (dH, dS) per dinucleotide step key
_MISMATCH_INTERNAL: Dict[str, Tuple[float, float]] = {}
_MISMATCH_TERMINAL: Dict[str, Tuple[float, float]] = {}

# Initiation entropique des boucles internes ADN à 37 °C en kcal/mol (SantaLucia & Hicks 2004, Table 5 ;
# identique à la table `internal_loop` de Mathews 2004 dans ViennaRNA).
# Entropic DNA internal loop initiation at 37 °C in kcal/mol (SantaLucia & Hicks 2004, Table 5;
# identical to the Mathews 2004 `internal_loop` table in ViennaRNA).
_INTERNAL_LOOP_DG37_KCAL: Dict[int, float] = {
    4: 3.1, 5: 3.5, 6: 3.9, 7: 4.1, 8: 4.2, 9: 4.3, 10: 4.5,
    11: 4.6, 12: 4.6, 13: 4.7, 14: 4.8, 15: 4.9, 16: 5.0, 17: 5.0,
    18: 5.1, 19: 5.1, 20: 5.2, 21: 5.3, 22: 5.3, 23: 5.3, 24: 5.4,
    25: 5.4, 26: 5.5, 27: 5.5, 28: 5.6, 29: 5.6, 30: 5.6,
}
_T37_KELVIN = 310.15
_JACOBSON_STOCKMAYER_EXP = 1.75

def _load_tables():
    data_dir = os.path.join(os.path.dirname(__file__), "..", "data", "thermo")
    
    int_file = os.path.join(data_dir, "mismatch_nn_internal.csv")
    if os.path.exists(int_file):
        with open(int_file, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                key = row["nn_step_top_bottom"]
                if 'I' in key:
                    continue
                _MISMATCH_INTERNAL[key] = (float(row["dH_kcal_per_mol"]), float(row["dS_cal_per_molK"]))
                
    term_file = os.path.join(data_dir, "mismatch_nn_terminal.csv")
    if os.path.exists(term_file):
        with open(term_file, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                key = row["nn_step_top_bottom"]
                if 'I' in key:
                    continue
                _MISMATCH_TERMINAL[key] = (float(row["dH_kcal_per_mol"]), float(row["dS_cal_per_molK"]))

_load_tables()

def _internal_loop_entropy_cal(loop_nt_count: int) -> float:
    """Calcule l'entropie ΔS (cal/mol·K, négative) d'une boucle interne symétrique de taille `loop_nt_count` >= 4.
    Computes entropy ΔS (cal/mol·K, negative) for a symmetric internal loop of size `loop_nt_count` >= 4.

    Référence / Reference:
    SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct. 33:415-440, Eq. 5 & Table 5
    (extrapolation de Jacobson-Stockmayer au-delà de 30 nt).
    """
    if loop_nt_count < 4:
        return 0.0
    if loop_nt_count <= 30:
        dg37 = _INTERNAL_LOOP_DG37_KCAL[loop_nt_count]
    else:
        dg37 = _INTERNAL_LOOP_DG37_KCAL[30] + (
            _JACOBSON_STOCKMAYER_EXP * R_GAS_KCAL_MOL_K * _T37_KELVIN * math.log(loop_nt_count / 30.0)
        )
    # Boucle interne modélisée comme purement entropique (ΔH = 0, ΔG°_37 = -T_37 * ΔS)
    # Internal loop initiation modeled as purely entropic (ΔH = 0, ΔG°_37 = -T_37 * ΔS)
    return -(dg37 * 1000.0) / _T37_KELVIN


def classify_mismatch_pair(primer_base: str, template_base: str) -> str:
    """Qualifie chimiquement une paire (amorce, matrice) selon Kwok et al. (1990).
    Chemically classifies a (primer, template) base pair following Kwok et al. (1990).

    Returns:
        'match' | 'wobble' | 'purine_purine' | 'pyrimidine_pyrimidine' | 'transversion_other'
    """
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    a = primer_base.upper()
    b = template_base.upper()
    if comp.get(a, '') == b:
        return "match"
    if (a, b) in (('G', 'T'), ('T', 'G')):
        return "wobble"
    purines = {'A', 'G'}
    pyrimidines = {'C', 'T'}
    if a in purines and b in purines:
        return "purine_purine"
    if a in pyrimidines and b in pyrimidines:
        return "pyrimidine_pyrimidine"
    return "transversion_other"


def nn_duplex_energy(
    top_5to3: str, 
    bottom_under_top: str,
    temp_celsius: float = 65.0,
    include_init: bool = True
) -> Tuple[float, float, float]:
    """Calcule l'énergie NN d'un duplexe (mésappariements internes, terminaux et consécutifs).
    Computes NN duplex energy handling single internal, terminal, and consecutive mismatches.

    Modèle / Model:
    - Pas Watson-Crick : SantaLucia (1998)
    - Mésappariement interne isolé (1×1) : Allawi & SantaLucia (1997-1998), Peyret (1999)
    - Mésappariements terminaux (5' ou 3') : SantaLucia & Hicks (2004), Table 4 (`_MISMATCH_TERMINAL`)
    - Mésappariements internes consécutifs (k×k, k >= 2) : SantaLucia & Hicks (2004), Eq. 5 & Table 5
      (empilements terminaux fermants + entropie de boucle interne symétrique).
    
    Args:
        top_5to3: Séquence du brin haut (5'->3').
        bottom_under_top: Séquence du brin bas, antiparallèle, alignée sous le brin haut.
        temp_celsius: Température.
        include_init: Si True, inclut les pénalités d'initiation terminales sur les paires WC fermantes.
        
    Returns:
        (dH_kcal, dS_cal, dG_kcal)
    """
    if len(top_5to3) != len(bottom_under_top):
        raise ValueError("Les deux brins doivent avoir la même longueur.")

    n = len(top_5to3)
    if n == 0:
        return 0.0, 0.0, 0.0

    top = top_5to3.upper()
    bottom = bottom_under_top.upper()
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    is_wc = [complement.get(top[i], '') == bottom[i] for i in range(n)]
    wc_indices = [i for i, ok in enumerate(is_wc) if ok]

    temp_k = temp_celsius + 273.15

    # Aucun appariement Watson-Crick dans tout l'alignement : aucun duplexe stable formé
    # No Watson-Crick base pair anywhere in the alignment: no stable duplex formed
    if not wc_indices:
        if include_init:
            dh_unb = 2.0 * _INIT_AT[0]
            ds_unb = 2.0 * _INIT_AT[1]
            return dh_unb, ds_unb, dh_unb - temp_k * (ds_unb / 1000.0)
        return 0.0, 0.0, 0.0

    first_wc = wc_indices[0]
    last_wc = wc_indices[-1]

    dh_total = 0.0
    ds_total = 0.0

    if include_init:
        # Pénalités d'initiation sur les vraies paires Watson-Crick fermantes du cœur apparié
        # Initiation penalties on the actual closing Watson-Crick base pairs of the paired core
        for idx in (first_wc, last_wc):
            if (top[idx], bottom[idx]) in (('G', 'C'), ('C', 'G')):
                dh_total += _INIT_GC[0]
                ds_total += _INIT_GC[1]
            else:
                dh_total += _INIT_AT[0]
                ds_total += _INIT_AT[1]

    # Mésappariement terminal côté 5' du brin haut (adjacent à first_wc)
    # 5'-end terminal mismatch on top strand (immediately adjacent to first_wc)
    if first_wc >= 1:
        # Rotation antiparallèle à 180° : la paire WC (first_wc) passe à gauche, le mismatch (first_wc-1) à droite
        # 180° antiparallel rotation: WC pair (first_wc) on the left, terminal mismatch (first_wc-1) on the right
        tmm_left = f"{bottom[first_wc]}{bottom[first_wc - 1]}/{top[first_wc]}{top[first_wc - 1]}"
        if tmm_left in _MISMATCH_TERMINAL:
            dh, ds = _MISMATCH_TERMINAL[tmm_left]
            dh_total += dh
            ds_total += ds

    # Mésappariement terminal côté 3' du brin haut (adjacent à last_wc)
    # 3'-end terminal mismatch on top strand (immediately adjacent to last_wc)
    if last_wc <= n - 2:
        tmm_right = f"{top[last_wc:last_wc + 2]}/{bottom[last_wc:last_wc + 2]}"
        if tmm_right in _MISMATCH_TERMINAL:
            dh, ds = _MISMATCH_TERMINAL[tmm_right]
            dh_total += dh
            ds_total += ds

    # Parcours du cœur apparié [first_wc .. last_wc]
    # Traverse paired core [first_wc .. last_wc]
    i = first_wc
    while i < last_wc:
        if is_wc[i] and is_wc[i + 1]:
            key = _SEQ_TO_NN_KEY.get(top[i:i + 2])
            if key is not None:
                dh, ds = _NN_PARAMS[key]
                dh_total += dh
                ds_total += ds
            i += 1
        else:
            # Bloc de k_mm >= 1 mésappariements internes consécutifs entre deux paires WC (i et j)
            # Run of k_mm >= 1 consecutive internal mismatches bounded by WC pairs at i and j
            j = i + 1
            while j < last_wc and not is_wc[j]:
                j += 1
            k_mm = j - (i + 1)

            if k_mm == 1:
                # Mésappariement interne simple (1×1) : somme des deux pas `_MISMATCH_INTERNAL`
                # Single internal mismatch (1×1): sum of the two flanking `_MISMATCH_INTERNAL` steps
                for pos in (i, i + 1):
                    dinuc_top = top[pos:pos + 2]
                    dinuc_bot = bottom[pos:pos + 2]
                    key_fwd = f"{dinuc_top}/{dinuc_bot}"
                    key_rev = f"{dinuc_bot[::-1]}/{dinuc_top[::-1]}"
                    if key_fwd in _MISMATCH_INTERNAL:
                        dh, ds = _MISMATCH_INTERNAL[key_fwd]
                        dh_total += dh
                        ds_total += ds
                    elif key_rev in _MISMATCH_INTERNAL:
                        dh, ds = _MISMATCH_INTERNAL[key_rev]
                        dh_total += dh
                        ds_total += ds
                    else:
                        # Repli défensif (ex. caractère ambigu résiduel) : demi-boucle entropique
                        # Defensive fallback (e.g. residual ambiguous base): half internal loop penalty
                        ds_total += 0.5 * _internal_loop_entropy_cal(4)
            else:
                # Boucle interne symétrique k_mm × k_mm (k_mm >= 2) :
                # Empilement fermant gauche + empilement fermant droit + initiation entropique Table 5
                # Symmetric k_mm × k_mm internal loop (k_mm >= 2):
                # Left closing stack + right closing stack + Table 5 entropic loop initiation
                left_stack = f"{top[i:i + 2]}/{bottom[i:i + 2]}"
                right_stack = f"{bottom[j - 1:j + 1][::-1]}/{top[j - 1:j + 1][::-1]}"
                if left_stack in _MISMATCH_TERMINAL:
                    dh, ds = _MISMATCH_TERMINAL[left_stack]
                    dh_total += dh
                    ds_total += ds
                if right_stack in _MISMATCH_TERMINAL:
                    dh, ds = _MISMATCH_TERMINAL[right_stack]
                    dh_total += dh
                    ds_total += ds
                ds_total += _internal_loop_entropy_cal(2 * k_mm)

            i = j

    dg_kcal = dh_total - temp_k * (ds_total / 1000.0)
    return dh_total, ds_total, dg_kcal


def three_prime_extensible(
    primer_5to3: str, 
    template_under_primer: str, 
    enzyme: PolymeraseProfile,
    window: Optional[int] = None,
    allow_penultimate_wobble: Optional[bool] = None,
) -> Tuple[bool, Optional[int], Optional[str]]:
    """Évalue l'extensibilité du 3' selon la présence et la chimie des mésappariements.
    Evaluates 3' extensibility based on mismatch position and base-pair chemistry.
    
    Cette fonction implémente le veto cinétique basé sur le modèle ARMS
    (Newton et al. 1989, NAR 17(7):2503-2516) et la qualification wobble G·T
    (Kwok et al. 1990, NAR 18(4):999-1005).
    
    Returns:
        (extensible, first_bad_pos, severity)
        - first_bad_pos est la position (1-indexé depuis le 3') du premier mismatch fautif,
          ou None si tout est OK.
        - severity peut valoir "block" (zone absolue, ext bloquée), "wobble_weak" (wobble G·T
          toléré en position -2 si activé), "weak" (zone pénalisée), ou None.
    """
    if window is None:
        window = enzyme.three_prime_window
    if allow_penultimate_wobble is None:
        allow_penultimate_wobble = getattr(enzyme, "allow_penultimate_gt_wobble", False)
        
    p = primer_5to3.upper()
    t = template_under_primer.upper()
    
    # Parcours depuis le 3' (position terminale = 1)
    # Scan from 3' end (terminal nucleotide = position 1)
    for pos_from_3p in range(1, window + 1):
        idx = len(p) - pos_from_3p
        if idx < 0:
            break
            
        pair_class = classify_mismatch_pair(p[idx], t[idx])
        if pair_class != "match":
            if pos_from_3p <= enzyme.three_prime_absolute_window:
                # Si l'option wobble G·T pénultième (-2) est activée, un wobble G·T non-terminal fuit
                # If penultimate (-2) G·T wobble tolerance is enabled, a non-terminal G·T wobble leaks
                if allow_penultimate_wobble and pos_from_3p == 2 and pair_class == "wobble":
                    return True, pos_from_3p, "wobble_weak"
                return False, pos_from_3p, "block"
            else:
                # Pénalité forte mais extension possible (fuite ARMS)
                # Strong penalty but leaky extension possible (ARMS)
                return True, pos_from_3p, "weak"
            
    return True, None, None


def evaluate_primer_target_binding(
    primer_5to3: str, 
    target_region_5to3: str, 
    enzyme: PolymeraseProfile,
    max_mismatches: int, 
    dg_cutoff: float,
    temp_celsius: float = 65.0
) -> dict:
    """Évalue la liaison d'une amorce avec tolérance aux mésappariements."""
    # Le "bottom_under_top" est le brin cible inversé de gauche à droite.
    # The "bottom_under_top" is the target strand reversed from left to right.
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    bottom_under = target_region_5to3[::-1]
    
    n_mismatches = 0
    for i in range(len(primer_5to3)):
        if complement.get(primer_5to3[i], '') != bottom_under[i]:
            n_mismatches += 1
            
    # Énergie
    _, _, dg_binding = nn_duplex_energy(primer_5to3, bottom_under, temp_celsius)
    
    # Veto 3'
    _, bad_pos, severity = three_prime_extensible(primer_5to3, bottom_under, enzyme)
    
    # Pour la couverture cible, tout mismatch dans la fenêtre pose problème (sévérité "weak" ou "block")
    # On considère que la variante n'est pas bien couverte.
    covered = (n_mismatches <= max_mismatches) and (severity is None) and (dg_binding <= dg_cutoff)
    
    return {
        "n_mismatches": n_mismatches,
        "dg_binding": dg_binding,
        "three_prime_ok": severity is None,
        "first_bad_pos": bad_pos,
        "covered": covered
    }


def calculate_hybridization_dg(
    primer_seq: str,
    bottom_under_top: str,
    temp_celsius: float,
    backend,
    bd_lna=(),
    **backend_kwargs
) -> tuple[float, float]:
    """
    Calcule le dG d'hybridation en utilisant le backend par défaut pour un
    duplex parfait, et en appliquant un ddG (shift) calculé via nn_duplex_energy
    s'il y a des mésappariements.
    
    Args:
        primer_seq: Séquence de l'amorce 5'->3'
        bottom_under_top: Séquence cible antiparallèle sous l'amorce
        temp_celsius: Température d'hybridation
        backend: Le backend thermodynamique (ex: ViennaSaltShiftBackend)
        
    Returns:
        (dg_hyb_final, ddg_mismatch)
    """
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    n_mismatches = sum(1 for a, b in zip(primer_seq, bottom_under_top, strict=True) if comp.get(a, '') != b)

    def _revcomp(s):
        return "".join(comp.get(c, 'N') for c in reversed(s))

    if n_mismatches == 0:
        res_hyb = backend.calc_duplex(
            primer_seq, _revcomp(primer_seq),
            temp_celsius=temp_celsius,
            lna_positions_a=bd_lna,
            lna_positions_b=(),
            **backend_kwargs
        )
        return res_hyb.dg_kcal, 0.0

    res_hyb_perfect = backend.calc_duplex(
        primer_seq, _revcomp(primer_seq),
        temp_celsius=temp_celsius,
        lna_positions_a=bd_lna,
        lna_positions_b=(),
        **backend_kwargs
    )
    dg_hyb_base = res_hyb_perfect.dg_kcal

    perfect_bottom = "".join(comp.get(c, c) for c in primer_seq)
    _, _, dg_perfect_nn = nn_duplex_energy(primer_seq, perfect_bottom, temp_celsius)
    _, _, dg_mismatched_nn = nn_duplex_energy(primer_seq, bottom_under_top, temp_celsius)

    ddg_mismatch = dg_mismatched_nn - dg_perfect_nn
    return dg_hyb_base + ddg_mismatch, ddg_mismatch
