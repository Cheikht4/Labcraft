"""Constraints and filters for sequence optimization.

Filtres et contraintes pour l'optimisation des séquences.
"""
from typing import List, Dict, Tuple, Any, Optional
import RNA

from labcraft.diagnostics.amplifiable_dimer import is_amplifiable_dimer
from labcraft.diagnostics.enzyme import PolymeraseProfile
from labcraft.thermo.backends.base import DuplexEnergyBackend
from labcraft.thermo.backends.native import _NN_PARAMS, _SEQ_TO_NN_KEY

def calc_intrinsic_3p_dg(seq: str, window_3p: int = 6, temp_celsius: float = 63.0) -> float:
    """
    Calcule la stabilité intrinsèque (dG) des `window_3p` dernières bases
    appariées à leur complément parfait (modèle NN SantaLucia sans pénalité d'initiation).
    Plus c'est négatif, plus c'est "collant".
    """
    if len(seq) < window_3p:
        window_3p = len(seq)
        
    segment = seq[-window_3p:]
    
    dh_total = 0.0
    ds_total = 0.0
    
    for i in range(len(segment) - 1):
        dinuc = segment[i:i+2]
        key = _SEQ_TO_NN_KEY.get(dinuc)
        if key:
            dh, ds = _NN_PARAMS[key]
            dh_total += dh
            ds_total += ds
            
    temp_k = temp_celsius + 273.15
    return dh_total - temp_k * (ds_total / 1000.0)

def _check_pair_both_orders(
    seq_1: str,
    seq_2: str,
    enzyme: PolymeraseProfile,
    temp_celsius: float,
    backend: Optional[DuplexEnergyBackend] = None,
    **backend_kwargs: Any,
) -> Tuple[bool, float]:
    """Évalue une paire d'amorces dans les deux ordres de concaténation (1&2 et 2&1).
    Evaluates a primer pair in both concatenation orders (1&2 and 2&1).
    """
    salt_kwargs = {
        "na_mm": backend_kwargs.get("na_mm"),
        "k_mm": backend_kwargs.get("k_mm", 0.0),
        "tris_mm": backend_kwargs.get("tris_mm", 0.0),
        "mg_mm": backend_kwargs.get("mg_mm", 0.0),
        "dntp_mm": backend_kwargs.get("dntp_mm", 0.0),
    }

    if backend is not None and seq_1 == seq_2:
        res_aa = backend.calc_homodimer(seq_1, temp_celsius=temp_celsius, **backend_kwargs)
        struct_ab, mfe_ab = res_aa.structure, res_aa.dg_kcal
    elif backend is not None:
        res_ab = backend.calc_duplex(seq_1, seq_2, temp_celsius=temp_celsius, **backend_kwargs)
        struct_ab, mfe_ab = res_ab.structure, res_ab.dg_kcal
    else:
        struct_ab, mfe_ab = RNA.cofold(f"{seq_1}&{seq_2}")

    is_amp_ab, dg_ab, _, _ = is_amplifiable_dimer(
        seq_1, seq_2, struct_ab.replace('&', ''), mfe_ab, enzyme,
        temp_celsius=temp_celsius, **salt_kwargs,
    )
    if seq_1 == seq_2:
        return is_amp_ab, dg_ab

    if backend is not None:
        res_ba = backend.calc_duplex(seq_2, seq_1, temp_celsius=temp_celsius, **backend_kwargs)
        struct_ba, mfe_ba = res_ba.structure, res_ba.dg_kcal
    else:
        struct_ba, mfe_ba = RNA.cofold(f"{seq_2}&{seq_1}")

    is_amp_ba, dg_ba, _, _ = is_amplifiable_dimer(
        seq_2, seq_1, struct_ba.replace('&', ''), mfe_ba, enzyme,
        temp_celsius=temp_celsius, **salt_kwargs,
    )
    is_amp = is_amp_ab or is_amp_ba
    if is_amp_ab and is_amp_ba:
        worst_dg = min(dg_ab, dg_ba)
    elif is_amp_ab:
        worst_dg = dg_ab
    elif is_amp_ba:
        worst_dg = dg_ba
    else:
        worst_dg = min(dg_ab, dg_ba)
    return is_amp, worst_dg


def evaluate_variant_dimers(
    variant_seq: str,
    original_primers: List[Dict[str, Any]],
    modified_idx: int,
    enzyme: PolymeraseProfile,
    temp_celsius: float,
    backend: Optional[DuplexEnergyBackend] = None,
    **backend_kwargs: Any,
) -> List[Tuple[int, float]]:
    """
    Évalue tous les dimères impliquant le variant (hétérodimères et homodimère).
    Evaluates all dimers involving the variant (heterodimers and homodimer).
    
    Args:
        variant_seq: La séquence mutée / Mutated sequence.
        original_primers: Liste des amorces du panel (dict avec 'sequence') / Panel primers list.
        modified_idx: Index de l'amorce modifiée dans original_primers / Modified primer index.
        enzyme: Profil de polymérase / Polymerase profile.
        temp_celsius: Température d'évaluation / Evaluation temperature.
        backend: Backend thermodynamique optionnel (avec sel) / Optional salt-aware backend.
        
    Returns:
        Liste de (partner_idx, dg_3p) pour chaque dimère amplifiable trouvé.
    """
    amplifiable: List[Tuple[int, float]] = []
    
    # Doit être exécuté dans un contexte dna_session / dna_params(temp_celsius)
    # Must be executed inside a dna_session / dna_params(temp_celsius) context
    for j, p in enumerate(original_primers):
        seq_partner = variant_seq if j == modified_idx else p['sequence']
        is_amp, min_dg_3p = _check_pair_both_orders(
            variant_seq, seq_partner, enzyme, temp_celsius, backend=backend, **backend_kwargs
        )
        if is_amp:
            amplifiable.append((j, min_dg_3p))
            
    return amplifiable

def filter_no_new_artefacts(
    variant_seq: str,
    original_primers: List[Dict[str, Any]],
    modified_idx: int,
    enzyme: PolymeraseProfile,
    temp_celsius: float,
    backend: Optional[DuplexEnergyBackend] = None,
    **backend_kwargs: Any,
) -> Tuple[bool, float]:
    """
    Filtre un variant s'il possède un quelconque dimère amplifiable.
    Filters out a variant if it forms any amplifiable dimer.
    Retourne / Returns (is_valid, worst_dg_3p).
    """
    worst_dg_3p = 0.0
    
    for j, p in enumerate(original_primers):
        seq_partner = variant_seq if j == modified_idx else p['sequence']
        is_amp, min_dg_3p = _check_pair_both_orders(
            variant_seq, seq_partner, enzyme, temp_celsius, backend=backend, **backend_kwargs
        )
        
        if is_amp:
            return False, min_dg_3p
            
        if min_dg_3p < worst_dg_3p:
            worst_dg_3p = min_dg_3p
            
    return True, worst_dg_3p
