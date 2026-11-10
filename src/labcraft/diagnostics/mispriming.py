"""Detection of inter-target off-target mispriming.
Détection du mésamorçage non spécifique inter-cibles (off-target mispriming).
"""
from __future__ import annotations

import math
import re
from typing import List, Dict, Tuple, Optional, Any
from pydantic import BaseModel
from labcraft.lamp.domains import PhysicalPrimer
from labcraft.thermo.backends.base import DuplexEnergyBackend
from labcraft.diagnostics.enzyme import PolymeraseProfile
from labcraft.thermo.mismatch import calculate_hybridization_dg, three_prime_extensible


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
    delta_delta_g: Optional[float] = None
    affinity_ratio: Optional[float] = None
    position: Optional[int] = None
    strand: Optional[str] = None


def detect_inter_target_mispriming(
    primers: List[PhysicalPrimer],
    primer_to_panel: Dict[str, str],
    targets: Dict[str, str],
    backend: DuplexEnergyBackend,
    enzyme: PolymeraseProfile,
    temp_celsius: float,
    k_len: int = 5,
    ddg_max: float = 4.0,
    primer_to_target: Optional[Dict[str, str]] = None,
    legitimate_loci: Optional[Dict[str, Any]] = None,
    **backend_kwargs
) -> List[MisprimingRisk]:
    """Détecte les risques réels d'amorçage croisé illégitime entre panels et cibles hétérologues.
    Detects off-target mispriming risks between primers and heterologous targets.

    Méthode et critères :
    1. Grandeur mesurée : Énergie d'hybridation INTERMOLÉCULAIRE du duplexe amorce-matrice
       (identique au calcul pour un site légitime, excluant le repliement intramoléculaire de la cible).
    2. Extensibilité 3' : Ancrage 3' minimal de k_len=5 nt vérifié avec la règle d'élongation de la polymérase.
    3. Critère relatif : Signalement uniquement si ΔG s'approche à moins de ddg_max (4.0 kcal/mol)
       du duplexe parfait de l'amorce sur son site légitime (fraction d'affinité Ka_més / Ka_légit > 0.26%).
    """
    risks: List[MisprimingRisk] = []
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    r_gas_kcal = 1.9872e-3
    temp_k = temp_celsius + 273.15
    rt = r_gas_kcal * temp_k

    cleaned_targets = {}
    for t_id, t_seq in targets.items():
        clean_seq = re.sub(r'[^A-Za-z]', '', t_seq.upper())
        if clean_seq:
            cleaned_targets[t_id] = (clean_seq, _revcomp(clean_seq))

    # Regrouper les amorces par séquence unique
    # Group primers by unique sequence
    unique_primer_groups: Dict[Tuple[str, Tuple[int, ...], bool], List[PhysicalPrimer]] = {}
    for p in primers:
        p_seq = re.sub(r'[^A-Za-z]', '', p.sequence.upper())
        if len(p_seq) < k_len:
            continue
        key = (p_seq, p.lna_positions, p.blocked_3prime)
        unique_primer_groups.setdefault(key, []).append(p)

    for (p_seq, lna_pos, is_blocked), p_list in unique_primer_groups.items():
        if is_blocked:
            # Une amorce bloquée en 3' ne peut pas être étendue par la polymérase
            continue

        p_3p = p_seq[-k_len:]
        k_mer = _revcomp(p_3p)
        p_len = len(p_seq)

        # Calcul du duplex parfait de référence pour cette amorce
        # Perfect duplex reference for this primer
        bottom_perfect = "".join(comp.get(c, 'N') for c in p_seq)
        try:
            dg_perfect, _ = calculate_hybridization_dg(
                p_seq, bottom_perfect, temp_celsius, backend, bd_lna=lna_pos, **backend_kwargs
            )
        except Exception:
            dg_perfect = -15.0

        for t_id, (t_seq_sense, t_seq_rc) in cleaned_targets.items():
            eligible_primers = list(p_list)
            if not eligible_primers:
                continue

            t_len = len(t_seq_sense)
            loci_for_target = (
                legitimate_loci.get(t_id, {})
                if isinstance(legitimate_loci, dict) and isinstance(legitimate_loci.get(t_id), dict)
                else (legitimate_loci if isinstance(legitimate_loci, dict) else {})
            )

            for strand_dir, s_seq in [('+', t_seq_sense), ('-', t_seq_rc)]:
                # Attention à la convention de brin entre seeding/enumerate_complexes et k_mer :
                # Dans seeding/enumerate_complexes :
                # - brin '+' = l'amorce correspond au brin sens (+) de la cible -> elle s'hybride au brin antisens (-),
                #   donc k_mer (_revcomp(p_3p)) apparaît sur t_seq_rc (strand_dir == '-').
                # - brin '-' = le complément inverse de l'amorce apparaît sur le brin sens (+),
                #   donc k_mer apparaît sur t_seq_sense (strand_dir == '+').
                # Candidate genomic interval on clean_target [cand_start, cand_end) and seeding strand:
                seeding_strand = '-' if strand_dir == '+' else '+'

                idx = 0
                while True:
                    idx = s_seq.find(k_mer, idx)
                    if idx == -1:
                        break

                    if idx + p_len > len(s_seq):
                        idx += 1
                        continue

                    if strand_dir == '+':
                        cand_start = idx
                        cand_end = idx + p_len
                    else:
                        cand_start = t_len - (idx + p_len)
                        cand_end = t_len - idx

                    # Filtrer les amorces pour lesquelles ce site chevauche le locus légitime retenu en Partie B
                    # Filter out primers for which this site overlaps the legitimate locus retained in Part B
                    non_legit_primers = []
                    for p in eligible_primers:
                        par_name = p.parent_name or (p.name.split('#')[0] if '#' in p.name else p.name)
                        legit_info = (
                            loci_for_target.get(f"{par_name}_site")
                            or loci_for_target.get(f"{par_name}@{t_id}_site")
                            or loci_for_target.get(par_name)
                        )
                        is_legit_overlap = False
                        if isinstance(legit_info, dict) and "start" in legit_info and "end" in legit_info:
                            l_start = int(legit_info["start"])
                            l_end = int(legit_info["end"])
                            l_strand = legit_info.get("strand")
                            # Même brin (ou brin non précisé) et chevauchement des coordonnées
                            if (l_strand is None or l_strand == seeding_strand):
                                overlap = min(cand_end, l_end) - max(cand_start, l_start)
                                if overlap > 0:
                                    is_legit_overlap = True
                        if not is_legit_overlap:
                            non_legit_primers.append(p)

                    if not non_legit_primers:
                        idx += 1
                        continue

                    # Séquence sous l'amorce sur la matrice (5'->3' sur le brin cible)
                    region = s_seq[idx:idx + p_len]
                    # Brin sous l'amorce orienté 3'->5' (antiparallèle sous p_seq)
                    bottom_under_top = region[::-1]

                    # 1. Vérification de l'extensibilité 3'
                    extensible, _first_bad_pos, severity = three_prime_extensible(p_seq, bottom_under_top, enzyme)
                    if not extensible or severity == "block":
                        idx += 1
                        continue

                    # 2. Pré-filtrage rapide : taux d'appariement minimal (> 50% pour un risque réel)
                    match_count = sum(1 for a, b in zip(reversed(p_seq), region, strict=True) if comp.get(a, '') == b)
                    if match_count < 0.50 * p_len:
                        idx += 1
                        continue

                    # 3. Calcul de l'énergie d'hybridation intermoléculaire
                    try:
                        dg_hyb, _ddg = calculate_hybridization_dg(
                            p_seq, bottom_under_top, temp_celsius, backend, bd_lna=lna_pos, **backend_kwargs
                        )
                    except Exception:
                        dg_hyb = 0.0

                    # Un mésamorçage ne peut pas être plus stable qu'un match parfait
                    if dg_hyb < dg_perfect:
                        dg_hyb = dg_perfect

                    ddg_rel = max(0.0, dg_hyb - dg_perfect)

                    # 4. Critère relatif : ΔΔG <= ddg_max (affinité relative significative)
                    if dg_hyb < 0 and ddg_rel <= ddg_max:
                        aff_ratio = math.exp(-ddg_rel / rt) if rt > 0 else 0.0
                        struct = "(" * match_count + "." * (p_len - match_count)

                        for p in non_legit_primers:
                            risks.append(MisprimingRisk(
                                primer_name=p.name,
                                target_id=t_id,
                                delta_g=round(dg_hyb, 2),
                                structure=struct,
                                site_sequence=region,
                                delta_delta_g=round(ddg_rel, 2),
                                affinity_ratio=round(aff_ratio, 4),
                                position=cand_start,
                                strand=seeding_strand
                            ))

                    idx += 1

    # Trier par gravité (duplex le plus stable / affinité relative la plus forte)
    risks.sort(key=lambda r: (r.delta_g, -(r.affinity_ratio or 0.0)))
    return risks
