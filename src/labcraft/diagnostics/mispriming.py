"""Detection of inter-target off-target mispriming.
Détection du mésamorçage non spécifique inter-cibles (off-target mispriming).
"""
from __future__ import annotations

import re
from typing import List, Dict, Tuple, Optional, Any
from pydantic import BaseModel
from labcraft.lamp.domains import PhysicalPrimer
from labcraft.thermo.backends.base import DuplexEnergyBackend
from labcraft.diagnostics.enzyme import PolymeraseProfile
from labcraft.diagnostics.amplifiable_dimer import is_amplifiable_dimer


def _revcomp(seq: str) -> str:
    """Calcule le complément inverse d'une séquence nucléotidique.
    Computes reverse complement of a nucleotide sequence.
    """
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'U': 'T',
        'M': 'K', 'R': 'Y', 'W': 'W', 'S': 'S', 'Y': 'R', 'K': 'M',
        'V': 'B', 'H': 'D', 'D': 'H', 'B': 'V',
        'N': 'N'
    }
    clean_seq = re.sub(r'[^A-Za-z]', '', seq.upper())
    return "".join(complement.get(c, 'N') for c in reversed(clean_seq))


class MisprimingRisk(BaseModel):
    primer_name: str
    target_id: str
    delta_g: float
    structure: str
    site_sequence: str


def detect_inter_target_mispriming(
    primers: List[PhysicalPrimer],
    primer_to_panel: Dict[str, str],
    targets: Dict[str, str],
    backend: DuplexEnergyBackend,
    enzyme: PolymeraseProfile,
    temp_celsius: float,
    k_len: int = 5,
    **backend_kwargs
) -> List[MisprimingRisk]:
    """Détecte les risques d'amorçage croisé illégitime entre panels et cibles hétérologues.
    Detects off-target mispriming risks between primers and heterologous targets.

    Optimisation :
    1. Ancrage 3' de k_len=5 bases (une hybridation < 5 nt à 65°C est instable et incapable d'amorcer).
    2. Pré-filtrage rapide par taux d'appariement minimal (> 40%) sur la fenêtre d'alignement.
    3. Mise en cache des calculs de duplex pour variantes ou motifs identiques.
    """
    risks: List[MisprimingRisk] = []
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    duplex_cache: Dict[Tuple[str, str], Any] = {}

    cleaned_targets = {}
    for t_id, t_seq in targets.items():
        clean_seq = re.sub(r'[^A-Za-z]', '', t_seq.upper())
        if clean_seq:
            cleaned_targets[t_id] = (clean_seq, _revcomp(clean_seq))

    # Regrouper les amorces par séquence unique pour éviter les scans redondants
    # Group primers by unique sequence to avoid redundant scanning
    unique_primer_groups: Dict[Tuple[str, Tuple[int, ...], bool], List[PhysicalPrimer]] = {}
    for p in primers:
        p_seq = re.sub(r'[^A-Za-z]', '', p.sequence.upper())
        if len(p_seq) < k_len:
            continue
        key = (p_seq, p.lna_positions, p.blocked_3prime)
        unique_primer_groups.setdefault(key, []).append(p)

    for (p_seq, lna_pos, is_blocked), p_list in unique_primer_groups.items():
        p_3p = p_seq[-k_len:]
        k_mer = _revcomp(p_3p)
        p_len = len(p_seq)

        for t_id, (t_seq_sense, t_seq_rc) in cleaned_targets.items():
            # Vérifier si au moins une amorce du groupe cible une autre cible que t_id
            eligible_primers = [
                p for p in p_list 
                if (primer_to_panel.get(p.name) or primer_to_panel.get(p.name.split('#')[0])) != t_id
            ]
            if not eligible_primers:
                continue

            for strand_dir, s_seq in [('+', t_seq_sense), ('-', t_seq_rc)]:
                idx = 0
                while True:
                    idx = s_seq.find(k_mer, idx)
                    if idx == -1:
                        break

                    start = max(0, idx - 5)
                    end = min(len(s_seq), idx + p_len + 5)
                    site = s_seq[start:end]

                    # Pré-filtrage rapide
                    region = s_seq[idx:idx + p_len]
                    match_count = sum(1 for a, b in zip(reversed(p_seq), region) if comp.get(a, '') == b)
                    if match_count < 0.40 * p_len:
                        idx += 1
                        continue

                    cache_key = (p_seq, site)
                    if cache_key in duplex_cache:
                        res = duplex_cache[cache_key]
                    else:
                        try:
                            res = backend.calc_heterodimer(
                                p_seq, site,
                                temp_celsius=temp_celsius,
                                lna_positions_a=lna_pos,
                                lna_positions_b=(),
                                **backend_kwargs
                            )
                            duplex_cache[cache_key] = res
                        except Exception:
                            res = None

                    if res is not None:
                        struct = res.structure.replace('&', '')
                        mfe = res.dg_kcal
                        is_amp, dg_3p, ext_strand, blocked = is_amplifiable_dimer(
                            p_seq, site, struct, mfe,
                            enzyme, temp_celsius,
                            blocked_a=is_blocked, blocked_b=False,
                            lna_positions=lna_pos
                        )

                        if is_amp and ext_strand == 'a' and mfe <= enzyme.dimer_dg_threshold:
                            for p in eligible_primers:
                                risks.append(MisprimingRisk(
                                    primer_name=p.name,
                                    target_id=t_id,
                                    delta_g=mfe,
                                    structure=res.structure,
                                    site_sequence=site
                                ))

                    idx += 1

    return risks
