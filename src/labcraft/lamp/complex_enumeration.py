"""Exhaustive complex enumeration / Énumération exhaustive des complexes.

Génère toutes les interactions bimoléculaires pour un panel de primers.
"""
from __future__ import annotations

import re
import warnings
import numpy as np

from dataclasses import dataclass
from typing import Sequence, List, Tuple, Dict
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
    target_seq: str,
    backend: DuplexEnergyBackend,
    profile: ConcentrationProfile = LAMP_DEFAULT_PROFILE,
    temp_celsius: float = 65.0,
    mon_molar: float | None = None,
    buffer: dict | None = None,
    unfolding_window: int = 150
) -> Tuple[EquilibriumProblem, List[str], List[str], Dict[str, float]]:
    """Énumère toutes les espèces et génère le problème d'équilibre.
    
    1. Espèces de base = chaque amorce + sites cibles identifiés.
    2. Monomères libres (structurés en épingle).
    3. Dimères d'amorces (oligo entier vs oligo entier).
    4. Complexes amorce-cible (domaine de liaison vs cible).
    
    Returns:
        (EquilibriumProblem, noms_des_especes, noms_des_complexes, unfolding_penalties)
    """
    # 1. Identifier les sites cibles (et gérer les chevauchements)
    target_sites = []
    primer_to_site = {}
    
    backend_kwargs = {}
    if buffer:
        backend_kwargs = {
            'na_mm': buffer.get('na_mM', 50.0),
            'k_mm': buffer.get('k_mM', 0.0),
            'tris_mm': buffer.get('tris_mM', 0.0),
            'mg_mm': buffer.get('mg_mM', 0.0),
            'dntp_mm': buffer.get('dntp_mM', 0.0)
        }

    if target_seq:
        target_seq_upper = target_seq.upper()
        target_rc = _revcomp(target_seq_upper)
        
        for p in primers:
            match_start = _find_iupac_substring(p.binding_domain, target_seq_upper)
            strand = "+"
            if match_start == -1:
                # Chercher sur le brin - (donc dans target_rc, mais il faut remaper les indices)
                match_start = _find_iupac_substring(p.binding_domain, target_rc)
                strand = "-"
                
            if match_start != -1:
                match_len = len(p.binding_domain)
                if strand == "+":
                    start, end = match_start, match_start + match_len
                else:
                    # Si match sur le RC, l'indice 0 du RC est len - 1 du +.
                    # rc_start .. rc_end (exclusif) correspond à len - rc_end .. len - rc_start
                    start = len(target_seq) - (match_start + match_len)
                    end = len(target_seq) - match_start
                    
                site_name = f"{p.name}_site"
                target_sites.append({
                    "name": site_name,
                    "start": start,
                    "end": end,
                    "strand": strand
                })
                primer_to_site[p.name] = site_name
            else:
                warnings.warn(f"Le domaine de liaison de {p.name} n'est pas trouvé sur la cible.")

        # Détection de chevauchement stérique (Avertissement)
        for i, s1 in enumerate(target_sites):
            for j, s2 in enumerate(target_sites):
                if i < j and s1["strand"] == s2["strand"]:
                    overlap = max(0, min(s1["end"], s2["end"]) - max(s1["start"], s2["start"]))
                    if overlap > 0:
                        warnings.warn(f"Compétition stérique : {s1['name']} et {s2['name']} se chevauchent de {overlap} bases.")
                    
    # Espèces de base = amorces + sites cibles
    n_primers = len(primers)
    n_sites = len(target_sites)
    n_strands = n_primers + n_sites
    
    strand_names = [p.name for p in primers] + [s["name"] for s in target_sites]
    
    concentrations = np.zeros(n_strands)
    for i, p in enumerate(primers):
        concentrations[i] = p.nominal_concentration if p.nominal_concentration is not None else profile.get_concentration(p.role)
    for i in range(n_sites):
        concentrations[n_primers + i] = profile.target
        
    complexes = []
    
    # --- 2. Monomères libres ---
    # Pour l'instant, on suppose que le backend peut calculer l'épingle d'un monomère,
    # mais ViennaRNA est meilleur pour ça. 
    # Pour respecter la matrice : l'énergie du monomère libre est souvent mise à 0 (état de référence).
    # On va donc déclarer les monomères à 0.0 kcal/mol.
    # Si on calcule une épingle forte, on met un delta_G d'épingle < 0.
    # Mais le solveur Jalon 1 suppose que les composants de base sont les monomères déroulés (0.0).
    for i, p in enumerate(primers):
        stoich = [0] * n_strands
        stoich[i] = 1
        # L'énergie de la forme libre est prise comme référence 0
        # (les dimères seront calculés en relatif)
        complexes.append(ComplexInfo(f"{p.name}_free", stoich, 0.0))
        
    for i, s in enumerate(target_sites):
        stoich = [0] * n_strands
        stoich[n_primers + i] = 1
        complexes.append(ComplexInfo(f"{s['name']}_free", stoich, 0.0))
        
    # --- 3. Dimères d'amorces (Oligo entier vs Oligo entier) ---
    for i, p1 in enumerate(primers):
        for j, p2 in enumerate(primers):
            if j < i: continue # On ne compte qu'une fois la paire
            
            try:
                if i != j:
                    res = backend.calc_heterodimer(
                        p1.sequence, p2.sequence, 
                        temp_celsius=temp_celsius, 
                        lna_positions_a=p1.lna_positions,
                        lna_positions_b=p2.lna_positions,
                        **backend_kwargs
                    )
                else:
                    res = backend.calc_homodimer(
                        p1.sequence, 
                        temp_celsius=temp_celsius, 
                        lna_positions=p1.lna_positions,
                        **backend_kwargs
                    )
                dg = res.dg_kcal
            except ValueError as e:
                import logging
                logging.warning(f"Error calculating dimer {p1.name} - {p2.name}: {e}")
                dg = 1.0 # Ignorer ce complexe
            if dg < 0: # Ne retenir que les interactions stabilisantes
                stoich = [0] * n_strands
                stoich[i] += 1
                stoich[j] += 1
                cname = f"{p1.name}_{p2.name}" if i != j else f"{p1.name}_homo"
                complexes.append(ComplexInfo(cname, stoich, dg))

    # --- 4. Complexes amorce-cible ---
    unfolding_penalties = {}
    
    if target_seq:
        for idx_p, p in enumerate(primers):
            site_name = primer_to_site.get(p.name)
            if not site_name:
                continue
                
            site_idx = n_primers + next(k for k, s in enumerate(target_sites) if s["name"] == site_name)
            site_info = next(s for s in target_sites if s["name"] == site_name)
            
            # Extract lna_positions for the binding domain only
            offset = p.sequence.find(p.binding_domain)
            bd_lna = tuple(pos - offset for pos in p.lna_positions if offset <= pos < offset + len(p.binding_domain)) if offset != -1 else ()
            
            # Extraire la séquence RÉELLE du génome
            s0 = site_info["start"]
            e0 = site_info["end"]
            extracted_target = target_seq[s0:e0].upper()
            
            # Gestion des longueurs (fallback sur comportement précédent si longueur inattendue)
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
                # Orientation sous l'amorce (antiparallèle, gauche à droite = 5'->3' pour l'amorce, 3'->5' pour le template)
                if site_info["strand"] == "+":
                    # L'amorce s'hybride au brin -. La séquence du brin - sous l'amorce (lue 3'->5') est le complément exact de la cible extraite.
                    bottom_under_top = "".join({'A':'T','T':'A','C':'G','G':'C'}.get(c, c) for c in extracted_target)
                else:
                    # L'amorce s'hybride au brin +. Le brin + sous l'amorce (lu 3'->5') est l'inverse exact de la cible extraite.
                    bottom_under_top = extracted_target[::-1]
                
                # Résolution des IUPAC
                from labcraft.lamp.domains import IUPAC_MATCHABLE
                resolved_bottom = []
                for b_prim, b_targ in zip(p.binding_domain, bottom_under_top):
                    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
                    set_p = IUPAC_MATCHABLE.get(b_prim, set())
                    # b_targ is on the template strand
                    set_t_comp = set(comp.get(base, base) for base in IUPAC_MATCHABLE.get(b_targ, set([b_targ])))
                    if set_p and set_t_comp and set_p.intersection(set_t_comp):
                        shared = list(set_p.intersection(set_t_comp))[0]
                        resolved_bottom.append(comp[shared])
                    else:
                        resolved_bottom.append(b_targ)
                bottom_under_top = "".join(resolved_bottom)
                
                # Comptage exact des mésappariements et énergie
                from labcraft.thermo.mismatch import calculate_hybridization_dg
                dg_hyb, ddg_mismatch = calculate_hybridization_dg(
                    p.binding_domain, bottom_under_top, temp_celsius, backend, bd_lna=bd_lna, **backend_kwargs
                )
                
                comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
                n_mismatches = sum(1 for a, b in zip(p.binding_domain, bottom_under_top) if comp.get(a, '') != b)

                # Évaluation de l'extensibilité 3'
                from labcraft.diagnostics.enzyme import get_enzyme
                from labcraft.thermo.mismatch import three_prime_extensible
                # Fallback to Bst2.0 if not provided
                enzyme = backend_kwargs.get("enzyme", get_enzyme("Bst2.0"))
                extensible, first_bad_pos, severity = three_prime_extensible(p.binding_domain, bottom_under_top, enzyme)
                
                site_info["mismatches"] = n_mismatches
                site_info["extensible"] = extensible
            
            # Le calcul de l'accessibilité
            if p.role in (PrimerRole.LF, PrimerRole.LB):
                dg_unfold = 0.0
            else:
                W = unfolding_window
                s0 = site_info["start"]
                e0 = site_info["end"]
                win_start = max(0, s0 - W)
                win_end = min(len(target_seq), e0 + W)
                target_window = target_seq[win_start:win_end]
                local_start = s0 - win_start
                local_end = e0 - win_start
    
                dg_unfold = calc_unfolding_penalty(
                    target_window, local_start, local_end, 
                    temp_celsius=temp_celsius, mon_molar=mon_molar
                )
                
            unfolding_penalties[site_name] = {
                "dg_unfold": dg_unfold,
                "mismatches": site_info.get("mismatches", 0),
                "extensible": site_info.get("extensible", True)
            }
            
            # Couplage
            dg_eff = dg_hyb + dg_unfold
            
            # Application du veto ARMS 3'
            # Si le 3' n'est pas extensible, le complexe ne peut pas initier l'amplification.
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
