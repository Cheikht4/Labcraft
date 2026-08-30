"""Exhaustive complex enumeration / Énumération exhaustive des complexes.

Génère toutes les interactions bimoléculaires pour un panel de primers (monomères,
dimères, et complexes amorce-cible résolus conjointement ou séparément).
Generates all bimolecular interactions for a primer panel (monomers, dimers, and
primer-target complexes resolved jointly or separately).
"""
from __future__ import annotations

import re
import warnings
import numpy as np

from dataclasses import dataclass
from typing import Sequence, List, Tuple, Dict, Optional, Union
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole, _find_iupac_substring
from labcraft.lamp.stoichiometry import ConcentrationProfile, LAMP_DEFAULT_PROFILE
from labcraft.target.unfolding import calc_unfolding_penalty
from labcraft.thermo.backends.base import DuplexEnergyBackend
from labcraft.solver.types import EquilibriumProblem

@dataclass
class ComplexInfo:
    name: str
    stoichiometry: list[int]
    delta_g: float


def enumerate_complexes(
    primers: List[PhysicalPrimer],
    target_seq: Optional[Union[str, Dict[str, str]]] = None,
    backend: DuplexEnergyBackend = None,
    profile: Union[ConcentrationProfile, Dict[str, ConcentrationProfile]] = LAMP_DEFAULT_PROFILE,
    temp_celsius: float = 65.0,
    mon_molar: float | None = None,
    buffer: dict | None = None,
    unfolding_window: int = 150
) -> Tuple[EquilibriumProblem, List[str], List[str], Dict[str, Any]]:
    """Énumère toutes les espèces et génère le problème d'équilibre.

    Supporte une cible unique (chaîne) ou un dictionnaire multi-cibles ({t_id: t_seq})
    résolu conjointement dans le même système d'équilibre couplé.
    Supports a single target (str) or a multi-target dict ({t_id: t_seq})
    resolved jointly within the same coupled equilibrium system.

    Returns:
        (EquilibriumProblem, noms_des_especes, noms_des_complexes, unfolding_penalties)
    """
    backend_kwargs = {}
    if buffer:
        backend_kwargs = {
            'na_mm': buffer.get('na_mM', 50.0),
            'k_mm': buffer.get('k_mM', 0.0),
            'tris_mm': buffer.get('tris_mM', 0.0),
            'mg_mm': buffer.get('mg_mM', 0.0),
            'dntp_mm': buffer.get('dntp_mM', 0.0)
        }

    # Normalisation du dictionnaire des cibles
    # Normalization of targets dict
    targets_dict: Dict[str, str] = {}
    is_multi_target = False
    if isinstance(target_seq, dict):
        targets_dict = {k: v for k, v in target_seq.items() if v}
        is_multi_target = len(targets_dict) > 1
    elif isinstance(target_seq, str) and target_seq.strip():
        targets_dict = {"": target_seq.strip()}

    # 1. Identifier les sites cibles sur chaque cible
    # 1. Identify target sites across all targets
    target_sites: List[Dict[str, Any]] = []

    for t_id, raw_seq in targets_dict.items():
        clean_target = re.sub(r'[^A-Za-z]', '', raw_seq.upper())
        target_rc = _revcomp(clean_target)

        for p in primers:
            match_start = _find_iupac_substring(p.binding_domain, clean_target)
            strand = "+"
            if match_start == -1:
                match_start = _find_iupac_substring(p.binding_domain, target_rc)
                strand = "-"

            if match_start != -1:
                match_len = len(p.binding_domain)
                if strand == "+":
                    start, end = match_start, match_start + match_len
                else:
                    start = len(clean_target) - (match_start + match_len)
                    end = len(clean_target) - match_start

                if is_multi_target and t_id:
                    site_name = f"{p.name}@{t_id}_site"
                else:
                    site_name = f"{p.name}_site"

                target_sites.append({
                    "name": site_name,
                    "target_id": t_id,
                    "target_seq": clean_target,
                    "primer_name": p.name,
                    "start": start,
                    "end": end,
                    "strand": strand
                })
            else:
                if not is_multi_target:
                    warnings.warn(f"Le domaine de liaison de {p.name} n'est pas trouvé sur la cible.")

        # Détection de chevauchement stérique par cible
        sites_for_t = [s for s in target_sites if s.get("target_id") == t_id]
        for i, s1 in enumerate(sites_for_t):
            for j, s2 in enumerate(sites_for_t):
                if i < j and s1["strand"] == s2["strand"]:
                    overlap = max(0, min(s1["end"], s2["end"]) - max(s1["start"], s2["start"]))
                    if overlap > 0:
                        warnings.warn(f"Compétition stérique sur cible '{t_id}' : {s1['name']} et {s2['name']} se chevauchent de {overlap} bases.")

    # Espèces de base = amorces + sites cibles
    n_primers = len(primers)
    n_sites = len(target_sites)
    n_strands = n_primers + n_sites

    strand_names = [p.name for p in primers] + [s["name"] for s in target_sites]

    concentrations = np.zeros(n_strands)
    for i, p in enumerate(primers):
        if isinstance(profile, dict):
            # Prendre le premier profil disponible pour l'amorce
            p_prof = next(iter(profile.values())) if profile else LAMP_DEFAULT_PROFILE
        else:
            p_prof = profile
        concentrations[i] = p.nominal_concentration if p.nominal_concentration is not None else p_prof.get_concentration(p.role)

    for i, s in enumerate(target_sites):
        t_id = s.get("target_id", "")
        if isinstance(profile, dict) and t_id in profile:
            t_conc = profile[t_id].target
        elif isinstance(profile, ConcentrationProfile):
            t_conc = profile.target
        else:
            t_conc = LAMP_DEFAULT_PROFILE.target
        concentrations[n_primers + i] = t_conc

    complexes: List[ComplexInfo] = []

    # --- 2. Monomères libres ---
    for i, p in enumerate(primers):
        stoich = [0] * n_strands
        stoich[i] = 1
        complexes.append(ComplexInfo(f"{p.name}_free", stoich, 0.0))

    for i, s in enumerate(target_sites):
        stoich = [0] * n_strands
        stoich[n_primers + i] = 1
        complexes.append(ComplexInfo(f"{s['name']}_free", stoich, 0.0))

    # --- 3. Dimères d'amorces (Oligo entier vs Oligo entier) ---
    for i, p1 in enumerate(primers):
        for j, p2 in enumerate(primers):
            if j < i:
                continue

            try:
                if i == j:
                    res = backend.calc_homodimer(
                        p1.sequence,
                        temp_celsius=temp_celsius,
                        lna_positions=p1.lna_positions,
                        **backend_kwargs
                    )
                else:
                    res = backend.calc_heterodimer(
                        p1.sequence, p2.sequence,
                        temp_celsius=temp_celsius,
                        lna_positions_a=p1.lna_positions,
                        lna_positions_b=p2.lna_positions,
                        **backend_kwargs
                    )
                dg = res.dg_kcal
            except Exception:
                dg = 1.0

            if dg < 0:
                stoich = [0] * n_strands
                stoich[i] += 1
                stoich[j] += 1
                cname = f"{p1.name}_{p2.name}" if i != j else f"{p1.name}_homo"
                complexes.append(ComplexInfo(cname, stoich, dg))

    # --- 4. Complexes amorce-cible ---
    unfolding_penalties: Dict[str, Any] = {}

    for site_idx_rel, site_info in enumerate(target_sites):
        p_name = site_info["primer_name"]
        idx_p = next(i for i, p in enumerate(primers) if p.name == p_name)
        p = primers[idx_p]

        site_idx = n_primers + site_idx_rel
        site_name = site_info["name"]
        t_seq = site_info["target_seq"]
        t_id = site_info.get("target_id", "")

        offset = p.sequence.find(p.binding_domain)
        bd_lna = tuple(pos - offset for pos in p.lna_positions if offset <= pos < offset + len(p.binding_domain)) if offset != -1 else ()

        s0 = site_info["start"]
        e0 = site_info["end"]
        extracted_target = t_seq[s0:e0].upper()

        if len(extracted_target) != len(p.binding_domain):
            warnings.warn(f"Longueur inattendue pour le site {site_name} (site={len(extracted_target)}, amorce={len(p.binding_domain)}). Fallback sur match parfait.")
            res_hyb = backend.calc_duplex(
                p.binding_domain, _revcomp(p.binding_domain),
                temp_celsius=temp_celsius,
                lna_positions_a=bd_lna,
                lna_positions_b=(),
                **backend_kwargs
            )
            dg_hyb = res_hyb.dg_kcal
            extensible = True
            n_mismatches = 0
        else:
            if site_info["strand"] == "+":
                bottom_under_top = "".join({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}.get(c, c) for c in extracted_target)
            else:
                bottom_under_top = extracted_target[::-1]

            from labcraft.lamp.domains import IUPAC_MATCHABLE
            resolved_bottom = []
            for b_prim, b_targ in zip(p.binding_domain, bottom_under_top):
                comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
                set_p = IUPAC_MATCHABLE.get(b_prim, set())
                set_t_comp = set(comp.get(base, base) for base in IUPAC_MATCHABLE.get(b_targ, set([b_targ])))
                if set_p and set_t_comp and set_p.intersection(set_t_comp):
                    shared = list(set_p.intersection(set_t_comp))[0]
                    resolved_bottom.append(comp[shared])
                else:
                    resolved_bottom.append(b_targ)
            bottom_under_top = "".join(resolved_bottom)

            from labcraft.thermo.mismatch import calculate_hybridization_dg, three_prime_extensible
            dg_hyb, ddg_mismatch = calculate_hybridization_dg(
                p.binding_domain, bottom_under_top, temp_celsius, backend, bd_lna=bd_lna, **backend_kwargs
            )

            comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            n_mismatches = sum(1 for a, b in zip(p.binding_domain, bottom_under_top) if comp.get(a, '') != b)

            from labcraft.diagnostics.enzyme import get_enzyme
            enzyme = backend_kwargs.get("enzyme", get_enzyme("Bst2.0"))
            extensible, first_bad_pos, severity = three_prime_extensible(p.binding_domain, bottom_under_top, enzyme)

            site_info["mismatches"] = n_mismatches
            site_info["extensible"] = extensible

        # Accessibilité
        if p.role in (PrimerRole.LF, PrimerRole.LB):
            dg_unfold = 0.0
        else:
            W = unfolding_window
            win_start = max(0, s0 - W)
            win_end = min(len(t_seq), e0 + W)
            target_window = t_seq[win_start:win_end]
            local_start = s0 - win_start
            local_end = e0 - win_start

            dg_unfold = calc_unfolding_penalty(
                target_window, local_start, local_end,
                temp_celsius=temp_celsius, mon_molar=mon_molar
            )

        pen_data = {
            "dg_unfold": dg_unfold,
            "mismatches": site_info.get("mismatches", 0),
            "extensible": site_info.get("extensible", True)
        }

        if is_multi_target and t_id:
            if t_id not in unfolding_penalties:
                unfolding_penalties[t_id] = {}
            unfolding_penalties[t_id][site_name] = pen_data
            # Alias sans @t_id pour compatibilité
            unfolding_penalties[t_id][f"{p.name}_site"] = pen_data
        else:
            unfolding_penalties[site_name] = pen_data

        dg_eff = dg_hyb + dg_unfold

        if extensible and dg_eff < 0:
            stoich = [0] * n_strands
            stoich[idx_p] = 1
            stoich[site_idx] = 1
            complexes.append(ComplexInfo(f"{p.name}_on_{site_name}", stoich, dg_eff))

    stoich_matrix = np.array([c.stoichiometry for c in complexes], dtype=np.float64)
    dg_vector = np.array([c.delta_g for c in complexes], dtype=np.float64)
    complex_names = [c.name for c in complexes]

    prob = EquilibriumProblem(
        n_strands=n_strands,
        n_complexes=len(complexes),
        stoichiometry=stoich_matrix,
        delta_g=dg_vector,
        total_concentrations=concentrations,
        temperature_kelvin=273.15 + temp_celsius
    )
    return prob, strand_names, complex_names, unfolding_penalties


def _revcomp(seq: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complement.get(c, 'N') for c in reversed(seq))
