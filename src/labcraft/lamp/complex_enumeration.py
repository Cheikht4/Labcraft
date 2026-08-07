"""Exhaustive complex enumeration / Énumération exhaustive des complexes.

Génère toutes les interactions bimoléculaires pour un panel de primers.
"""
from __future__ import annotations

import re
import warnings
import numpy as np

from dataclasses import dataclass
from typing import Sequence
from labcraft.lamp.domains import PhysicalPrimer, _iupac_to_regex
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
    buffer: dict | None = None
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
            regex = _iupac_to_regex(p.binding_domain)
            # Chercher d'abord sur le brin +
            match = re.search(regex, target_seq_upper)
            strand = "+"
            if not match:
                # Chercher sur le brin - (donc dans target_rc, mais il faut remaper les indices)
                match = re.search(regex, target_rc)
                strand = "-"
                
            if match:
                if strand == "+":
                    start, end = match.start(), match.end()
                else:
                    # Si match sur le RC, l'indice 0 du RC est len - 1 du +.
                    # rc_start .. rc_end (exclusif) correspond à len - rc_end .. len - rc_start
                    start = len(target_seq) - match.end()
                    end = len(target_seq) - match.start()
                    
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
        concentrations[i] = profile.get_concentration(p.role)
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
                    res = backend.calc_heterodimer(p1.sequence, p2.sequence, temp_celsius=temp_celsius, **backend_kwargs)
                else:
                    res = backend.calc_homodimer(p1.sequence, temp_celsius=temp_celsius, **backend_kwargs)
                dg = res.dg_kcal
            except ValueError:
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
        for p in primers:
            site_name = primer_to_site.get(p.name)
            if not site_name:
                continue
                
            site_idx = n_primers + next(k for k, s in enumerate(target_sites) if s["name"] == site_name)
            site_info = next(s for s in target_sites if s["name"] == site_name)
            
            # Le backend calcule l'hybridation sur le domaine de liaison SEUL
            # On passe p.binding_domain et son reverse complement exact
            res_hyb = backend.calc_duplex(p.binding_domain, _revcomp(p.binding_domain), temp_celsius=temp_celsius, **backend_kwargs)
            dg_hyb = res_hyb.dg_kcal
            
            # Le calcul de l'accessibilité
            dg_unfold = calc_unfolding_penalty(
                target_seq, site_info["start"], site_info["end"], 
                temp_celsius=temp_celsius, mon_molar=mon_molar
            )
            unfolding_penalties[site_name] = dg_unfold
            
            # Couplage
            dg_eff = dg_hyb + dg_unfold
            
            if dg_eff < 0:
                stoich = [0] * n_strands
                stoich[i] = 1
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
