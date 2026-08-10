from __future__ import annotations
from labcraft.thermo.constants import R_GAS_KCAL_MOL_K

import numpy as np
from typing import Sequence
from dataclasses import dataclass

from labcraft.lamp.domains import PhysicalPrimer
from labcraft.lamp.complex_enumeration import enumerate_complexes
from labcraft.thermo.backends.base import DuplexEnergyBackend
from labcraft.solver.dual import solve_dual

@dataclass
class TitrationResult:
    dilution_factors: list[float]
    primer_names: list[str]
    free_fractions: np.ndarray # shape (n_dilutions, n_primers)
    dominant_complexes: dict[str, str] # primer_name -> dominant complex name at 1x
    complex_concentrations_1x: dict[str, float]

def simulate_titration(
    primers: Sequence[PhysicalPrimer],
    backend: DuplexEnergyBackend,
    base_concentrations: dict[str, float],
    dilutions: list[float] = [1.0, 0.5, 0.25, 0.125],
    temp_celsius: float = 65.0
) -> TitrationResult:
    """
    Calcule la compétition inter-amorces (sans cible) en fonction de la dilution.
    
    Args:
        primers: Liste des amorces (sans cible)
        backend: Moteur thermo
        base_concentrations: Dict {primer_name: conc_molar} à 1x
        dilutions: Facteurs multiplicatifs
        
    Returns:
        TitrationResult contenant les fractions libres
    """
    # 1. Énumération une seule fois (calcul de tous les delta G)
    # On passe None comme cible pour forcer le mode "sans cible"
    from labcraft.lamp.stoichiometry import ConcentrationProfile
    
    # On utilise un profil temporaire, on écrasera les concentrations totales
    # On met tout à 1e-6 pour ne pas planter
    dummy_profile = ConcentrationProfile(target=0, fip_bip=1e-6, f3_b3=1e-6, lf_lb=1e-6)
    
    prob_template, strands, complexes, _ = enumerate_complexes(
        primers, "", backend, profile=dummy_profile, temp_celsius=temp_celsius
    )
    
    n_primers = len(primers)
    n_dilutions = len(dilutions)
    
    free_fractions = np.zeros((n_dilutions, n_primers))
    
    # Vecteur des concentrations de base aligné sur strands
    c_tot_base = np.zeros(n_primers)
    for i, p in enumerate(primers):
        c_tot_base[i] = base_concentrations.get(p.name, 0.0)
        
    dom_complexes = {}
    c_complexes_1x = {}
        
    for d_idx, dil in enumerate(dilutions):
        # Mise à l'échelle
        c_tot_dil = c_tot_base * dil
        
        # Copie et mise à jour du problème
        import dataclasses
        prob = dataclasses.replace(prob_template, total_concentrations=c_tot_dil)
        
        res = solve_dual(prob)
        
        for i, p in enumerate(primers):
            free = res.free_concentrations[i]
            tot = c_tot_dil[i]
            free_fractions[d_idx, i] = (free / tot) * 100.0 if tot > 0 else 0.0
            
        # Si c'est la concentration standard (1.0), on extrait les complexes dominants
        if abs(dil - 1.0) < 1e-5:
            # Calculer la concentration de chaque complexe
            # C_j = exp(-dg_j/RT) * prod(u_i ^ A_ji)
            u = res.free_concentrations
            for j, c_name in enumerate(complexes):
                # c_name se termine par _free, on l'ignore pour la dominance des séquestrations
                if c_name.endswith("_free"):
                    continue
                # Calcul de la concentration
                RT = R_GAS_KCAL_MOL_K * prob.temperature_kelvin
                conc = np.exp(-prob.delta_g[j] / RT)
                for i in range(prob.n_strands):
                    if prob.stoichiometry[j, i] > 0:
                        conc *= (u[i] ** prob.stoichiometry[j, i])
                c_complexes_1x[c_name] = conc
                
            # Pour chaque primer avec fraction libre < 20%, trouver son pire séquestrateur
            for i, p in enumerate(primers):
                if free_fractions[d_idx, i] < 20.0:
                    worst_c = None
                    worst_val = -1
                    for j, c_name in enumerate(complexes):
                        if c_name.endswith("_free"):
                            continue
                        stoich_i = prob.stoichiometry[j, i]
                        if stoich_i > 0:
                            # Attention: un homodimère consomme 2 copies de l'amorce
                            # La quantité d'amorce séquestrée dans ce complexe est stoich_i * conc_complexe
                            c_val = stoich_i * c_complexes_1x.get(c_name, 0.0)
                            if c_val > worst_val:
                                worst_val = c_val
                                worst_c = c_name
                    if worst_c:
                        # Calculer le % de séquestration
                        pct = (worst_val / c_tot_dil[i]) * 100.0 if c_tot_dil[i] > 0 else 0
                        dom_complexes[p.name] = f"{worst_c} ({pct:.1f}%)"
                        
    return TitrationResult(
        dilution_factors=dilutions,
        primer_names=[p.name for p in primers],
        free_fractions=free_fractions,
        dominant_complexes=dom_complexes,
        complex_concentrations_1x=c_complexes_1x
    )
