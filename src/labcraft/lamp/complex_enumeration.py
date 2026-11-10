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
from typing import List, Tuple, Dict, Optional, Union, Any
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole, _find_iupac_substring, IUPAC_MATCHABLE
from labcraft.lamp.stoichiometry import ConcentrationProfile, LAMP_DEFAULT_PROFILE
from labcraft.target.unfolding import calc_unfolding_penalty
from labcraft.target.seeding import find_candidate_sites
from labcraft.thermo.backends.base import DuplexEnergyBackend
from labcraft.thermo.mismatch import calculate_hybridization_dg, three_prime_extensible
from labcraft.diagnostics.enzyme import PolymeraseProfile, get_enzyme
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
    unfolding_window: int = 150,
    enzyme: Optional[PolymeraseProfile] = None,
    max_errors: int = 2,
    strict_3prime_len: int = 3,
    strict_3prime_tolerate: int = 2,
    primer_to_target: Optional[Dict[str, str]] = None
) -> Tuple[EquilibriumProblem, List[str], List[str], Dict[str, Any]]:
    """Énumère toutes les espèces et génère le problème d'équilibre.

    Supporte une cible unique (chaîne) ou un dictionnaire multi-cibles ({t_id: t_seq})
    résolu conjointement dans le même système d'équilibre couplé.
    Supports a single target (str) or a multi-target dict ({t_id: t_seq})
    resolved jointly within the same coupled equilibrium system.

    Returns:
        (EquilibriumProblem, noms_des_especes, noms_des_complexes, unfolding_penalties)
    """
    if enzyme is None:
        enzyme = get_enzyme("Bst2.0")

    backend_kwargs = {}
    if buffer:
        backend_kwargs = {
            'na_mm': buffer.get('na_mM', buffer.get('na_mm', 50.0)),
            'k_mm': buffer.get('k_mM', buffer.get('k_mm', 0.0)),
            'tris_mm': buffer.get('tris_mM', buffer.get('tris_mm', 0.0)),
            'mg_mm': buffer.get('mg_mM', buffer.get('mg_mm', 0.0)),
            'dntp_mm': buffer.get('dntp_mM', buffer.get('dntp_mm', 0.0))
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

    # Regrouper les amorces par amorce parente
    # Group primers by parent primer
    parent_to_primers: Dict[str, List[PhysicalPrimer]] = {}
    for p in primers:
        par = p.parent_name or (p.name.split('#')[0] if '#' in p.name else p.name)
        parent_to_primers.setdefault(par, []).append(p)

    # 1. Identifier les sites cibles sur chaque cible par LOCUS physique unique
    # 1. Identify target sites across all targets by unique physical LOCUS
    target_sites_map: Dict[Tuple[str, int, str], Dict[str, Any]] = {}
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    for t_id, raw_seq in targets_dict.items():
        clean_target = re.sub(r'[^A-Za-z]', '', raw_seq.upper())
        target_rc = _revcomp(clean_target)

        # Filtrer les parents pertinents pour cette cible si primer_to_target est fourni en multi-cible
        # Filter relevant parents for this target if primer_to_target is provided in multi-target mode
        eligible_parents = list(parent_to_primers.keys())
        if is_multi_target and primer_to_target:
            filtered = [
                par for par in eligible_parents
                if (primer_to_target.get(par) or primer_to_target.get(parent_to_primers[par][0].name)) == t_id
            ]
            if filtered:
                eligible_parents = filtered

        parent_matched: Dict[str, bool] = {par: False for par in eligible_parents}

        # Criblage permissif via find_candidate_sites pour cette cible
        # Permissive screening via find_candidate_sites for this target
        primers_for_seeding = [
            p for par in eligible_parents for p in parent_to_primers[par]
        ]
        candidate_records = find_candidate_sites(
            {t_id: clean_target},
            primers_for_seeding,
            max_errors=max_errors,
            strict_3prime_len=strict_3prime_len,
            strict_3prime_tolerate=strict_3prime_tolerate
        )
        cands_by_parent: Dict[str, List[Dict[str, Any]]] = {}
        for rec in candidate_records:
            cands_by_parent.setdefault(rec["primer_name"], []).append(rec)

        for par_name in eligible_parents:
            p_variants = parent_to_primers[par_name]
            # Collecter les loci candidats (start, end, strand)
            # Collect candidate loci (start, end, strand)
            candidate_loci_set: set[Tuple[int, int, str]] = set()

            for rec in cands_by_parent.get(par_name, []):
                s_pos = int(rec["position"])
                s_strand = rec["strand"]
                s_len = len(rec["site_seq"])
                candidate_loci_set.add((s_pos, s_pos + s_len, s_strand))

            # Fallback sur _find_iupac_substring pour les rôles non couverts par seeding (ex: PROBE)
            # ou séquences courtes <= strict_3prime_len
            # Fallback on _find_iupac_substring for roles not covered by seeding (e.g. PROBE)
            # or short sequences <= strict_3prime_len
            if not candidate_loci_set:
                for p in p_variants:
                    m_start = _find_iupac_substring(p.binding_domain, clean_target)
                    m_strand = "+"
                    if m_start == -1:
                        m_start = _find_iupac_substring(p.binding_domain, target_rc)
                        m_strand = "-"
                    if m_start != -1:
                        m_len = len(p.binding_domain)
                        if m_strand == "+":
                            s0, e0 = m_start, m_start + m_len
                        else:
                            s0 = len(clean_target) - (m_start + m_len)
                            e0 = len(clean_target) - m_start
                        candidate_loci_set.add((s0, e0, m_strand))

            if not candidate_loci_set:
                continue

            # Évaluer chaque locus candidat pour retenir le meilleur site (extensible prioritaire, dG le plus favorable)
            # Evaluate each candidate locus to select the best site (extensible first, most favorable dG)
            best_locus_eval = None
            best_locus_key = None

            for (s0, e0, strand) in sorted(candidate_loci_set):
                extracted_target = clean_target[s0:e0].upper()
                has_n = 'N' in extracted_target
                any_extensible = False
                best_dg_hyb = float('inf')
                min_mm = 999
                best_bad_pos = None
                best_severity = None
                variant_evals = []

                for p in p_variants:
                    offset = p.sequence.find(p.binding_domain)
                    bd_lna = tuple(pos - offset for pos in p.lna_positions if offset <= pos < offset + len(p.binding_domain)) if offset != -1 else ()

                    if len(extracted_target) != len(p.binding_domain):
                        res_hyb = backend.calc_duplex(
                            p.binding_domain, _revcomp(p.binding_domain),
                            temp_celsius=temp_celsius,
                            lna_positions_a=bd_lna,
                            lna_positions_b=(),
                            **backend_kwargs
                        )
                        dg_hyb = res_hyb.dg_kcal
                        ext = True
                        bad_pos = None
                        sev = None
                        n_mm = 0
                    else:
                        if strand == "+":
                            bot_raw = "".join(comp.get(c, c) for c in extracted_target)
                        else:
                            bot_raw = extracted_target[::-1]

                        resolved_bottom = []
                        for b_prim, b_targ in zip(p.binding_domain, bot_raw, strict=True):
                            set_p = IUPAC_MATCHABLE.get(b_prim, set())
                            set_t_comp = set(comp.get(base, base) for base in IUPAC_MATCHABLE.get(b_targ, set([b_targ])))
                            if set_p and set_t_comp and set_p.intersection(set_t_comp):
                                shared = sorted(set_p.intersection(set_t_comp))[0]
                                resolved_bottom.append(comp[shared])
                            else:
                                resolved_bottom.append(b_targ)
                        bottom_under_top = "".join(resolved_bottom)

                        dg_hyb, _ = calculate_hybridization_dg(
                            p.binding_domain, bottom_under_top, temp_celsius, backend, bd_lna=bd_lna, **backend_kwargs
                        )
                        ext, bad_pos, sev = three_prime_extensible(p.binding_domain, bottom_under_top, enzyme)
                        n_mm = sum(1 for a, b in zip(p.binding_domain, bottom_under_top, strict=True) if comp.get(a, '') != b)

                    variant_evals.append((p, dg_hyb, ext, n_mm, bad_pos, sev))
                    if ext:
                        any_extensible = True
                        if dg_hyb < best_dg_hyb:
                            best_dg_hyb = dg_hyb
                        if n_mm < min_mm:
                            min_mm = n_mm
                    else:
                        if not any_extensible:
                            if dg_hyb < best_dg_hyb:
                                best_dg_hyb = dg_hyb
                            if n_mm < min_mm:
                                min_mm = n_mm
                            best_bad_pos = bad_pos
                            best_severity = sev

                # Clé de tri : extensible d'abord, puis moins de mésappariements, puis dG le plus négatif
                # Sort key: extensible first, then fewer mismatches, then most negative dG
                rank_key = (1 if any_extensible else 0, -min_mm, -best_dg_hyb)
                if best_locus_key is None or rank_key > best_locus_key:
                    best_locus_key = rank_key
                    best_locus_eval = {
                        "start": s0,
                        "end": e0,
                        "strand": strand,
                        "any_extensible": any_extensible,
                        "best_dg_hyb": best_dg_hyb,
                        "min_mm": min_mm,
                        "first_bad_pos": best_bad_pos,
                        "severity": best_severity,
                        "has_n": has_n,
                        "variant_evals": variant_evals,
                        "role": p_variants[0].role
                    }

            if best_locus_eval is not None:
                parent_matched[par_name] = True
                s0 = best_locus_eval["start"]
                e0 = best_locus_eval["end"]
                strand = best_locus_eval["strand"]
                locus_key = (t_id, s0, strand)

                if is_multi_target and t_id:
                    site_name = f"{par_name}@{t_id}_site"
                else:
                    site_name = f"{par_name}_site"

                if locus_key not in target_sites_map:
                    target_sites_map[locus_key] = {
                        "name": site_name,
                        "target_id": t_id,
                        "target_seq": clean_target,
                        "parent_name": par_name,
                        "role": best_locus_eval["role"],
                        "start": s0,
                        "end": e0,
                        "strand": strand,
                        "binding_primers": list(p_variants),
                        "variant_evals": list(best_locus_eval["variant_evals"]),
                        "min_mm": best_locus_eval["min_mm"],
                        "any_extensible": best_locus_eval["any_extensible"],
                        "best_dg_hyb": best_locus_eval["best_dg_hyb"],
                        "first_bad_pos": best_locus_eval["first_bad_pos"],
                        "severity": best_locus_eval["severity"],
                        "has_n": best_locus_eval["has_n"]
                    }
                else:
                    target_sites_map[locus_key]["binding_primers"].extend(p_variants)
                    target_sites_map[locus_key]["variant_evals"].extend(best_locus_eval["variant_evals"])

        if not is_multi_target:
            for par_name, matched in parent_matched.items():
                if not matched:
                    warnings.warn(f"Le domaine de liaison de {par_name} n'est pas trouvé sur la cible.")

        # Détection de chevauchement stérique par cible entre loci distincts
        sites_for_t = [s for s in target_sites_map.values() if s.get("target_id") == t_id]
        for i, s1 in enumerate(sites_for_t):
            for j, s2 in enumerate(sites_for_t):
                if i < j and s1["strand"] == s2["strand"]:
                    overlap = max(0, min(s1["end"], s2["end"]) - max(s1["start"], s2["start"]))
                    if overlap > 0:
                        warnings.warn(f"Compétition stérique sur cible '{t_id}' : {s1['name']} et {s2['name']} se chevauchent de {overlap} bases.")

    target_sites = list(target_sites_map.values())

    # Espèces de base = amorces + sites cibles
    n_primers = len(primers)
    n_sites = len(target_sites)
    n_strands = n_primers + n_sites

    strand_names = [p.name for p in primers] + [s["name"] for s in target_sites]

    concentrations = np.zeros(n_strands)
    for i, p in enumerate(primers):
        if isinstance(profile, dict):
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

    # --- 4. Complexes amorce-cible (par locus) ---
    unfolding_penalties: Dict[str, Any] = {}

    for site_idx_rel, site_info in enumerate(target_sites):
        site_idx = n_primers + site_idx_rel
        site_name = site_info["name"]
        t_seq = site_info["target_seq"]
        t_id = site_info.get("target_id", "")
        s0 = site_info["start"]
        e0 = site_info["end"]
        extracted_target = t_seq[s0:e0].upper()
        role = site_info.get("role")
        parent_name = site_info.get("parent_name")

        # Accessibilité (unfolding penalty calculée une seule fois par locus physique)
        if role in (PrimerRole.LF, PrimerRole.LB):
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

        productive_variants = [
            p.name for (p, dg_hyb, ext, n_mm, bad_pos, sev) in site_info.get("variant_evals", [])
            if ext and (dg_hyb + dg_unfold) < 0 and not site_info.get("has_n", False)
        ]

        pen_data = {
            "dg_unfold": dg_unfold,
            "dg_hyb": site_info.get("best_dg_hyb"),
            "mismatches": site_info.get("min_mm", 0),
            "extensible": site_info.get("any_extensible", True) and not site_info.get("has_n", False),
            "first_bad_pos": site_info.get("first_bad_pos"),
            "severity": site_info.get("severity"),
            "position": s0,
            "start": s0,
            "end": e0,
            "strand": site_info.get("strand"),
            "has_n": site_info.get("has_n", False),
            "verdict": "indetermine" if site_info.get("has_n", False) else ("veto_3p" if not site_info.get("any_extensible", True) else ("parfait" if site_info.get("min_mm", 0) == 0 else "tolerable")),
            "n_variants": len(productive_variants) if productive_variants else len(site_info["binding_primers"]),
            "total_variants": len(site_info["binding_primers"]),
            "variants": productive_variants if productive_variants else [p.name for p in site_info["binding_primers"]]
        }

        if is_multi_target and t_id:
            if t_id not in unfolding_penalties:
                unfolding_penalties[t_id] = {}
            unfolding_penalties[t_id][site_name] = pen_data
            unfolding_penalties[t_id][f"{parent_name}_site"] = pen_data
        else:
            unfolding_penalties[site_name] = pen_data

        for (p, dg_hyb, extensible, _n_mismatches, _bad_pos, _sev) in site_info.get("variant_evals", []):
            idx_p = next(i for i, prim in enumerate(primers) if prim.name == p.name)
            dg_eff = dg_hyb + dg_unfold

            if extensible and not site_info.get("has_n", False) and dg_eff < 0:
                stoich = [0] * n_strands
                stoich[idx_p] = 1
                stoich[site_idx] = 1
                cname = f"{p.name}_on_{site_name}"
                complexes.append(ComplexInfo(cname, stoich, dg_eff))

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
