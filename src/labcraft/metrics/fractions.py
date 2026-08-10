from dataclasses import dataclass
from typing import Dict, List
import numpy as np

@dataclass
class PrimerFractions:
    free: float
    hairpin: float
    homodimer: float
    heterodimer_intra: float # avec des amorces du même panel
    heterodimer_inter: float # hybridation croisée avec d'autres panels
    target_bound: float
    
    dominant_complex: str = ""
    dominant_fraction: float = 0.0
    
    @property
    def heterodimer(self) -> float:
        return self.heterodimer_intra + self.heterodimer_inter
        
    @property
    def sum(self) -> float:
        return self.free + self.hairpin + self.homodimer + self.heterodimer + self.target_bound

def compute_fractions(
    primer_names: List[str],
    complex_names: List[str],
    stoichiometry: np.ndarray,
    free_concentrations: np.ndarray,
    delta_g: np.ndarray,
    temp_celsius: float,
    primer_to_panel: Dict[str, str] = None
) -> Dict[str, PrimerFractions]:
    """
    Calcule la décomposition en pourcentages de chaque amorce.
    
    Args:
        primer_names: Noms des espèces de base (les N premières colonnes de la stoechiométrie)
        complex_names: Noms de tous les complexes (lignes de la stoechiométrie)
        stoichiometry: Matrice stoechiométrique (complexes x espèces)
        free_concentrations: Concentrations libres des espèces de base
        delta_g: Energies libres des complexes
        temp_celsius: Température
    """
    from labcraft.thermo.constants import R_GAS_KCAL_MOL_K
    R = R_GAS_KCAL_MOL_K
    RT = R * (273.15 + temp_celsius)
    
    u = np.log(free_concentrations)
    
    fractions = {}
    
    # Dictionnaires pour accumuler les contributions absolues par parent
    parent_totals = {}
    parent_contributions = {}
    parent_dominant = {}
    
    # 1. Obtenir le nom parent
    def get_parent(name: str) -> str:
        return name.split('#')[0] if '#' in name else name

    for j, p_name in enumerate(primer_names):
        if p_name.endswith("_site"):
            continue # On ne calcule pas les fractions de la cible
            
        parent = get_parent(p_name)
        if parent not in parent_totals:
            parent_totals[parent] = 0.0
            parent_contributions[parent] = {
                'free': 0.0, 'hairpin': 0.0, 'homo': 0.0,
                'hetero_intra': 0.0, 'hetero_inter': 0.0, 'target': 0.0
            }
            parent_dominant[parent] = []
            
        total = 0.0
        
        for i, c_name in enumerate(complex_names):
            coeff = stoichiometry[i, j]
            if coeff > 0:
                conc_complex = np.exp(-delta_g[i] / RT + np.dot(stoichiometry[i], u))
                contribution = coeff * conc_complex
                total += contribution
                
                # Classifier la contribution
                if c_name == f"{p_name}_free":
                    parent_contributions[parent]['free'] += contribution
                elif c_name == f"{p_name}_hairpin":
                    parent_contributions[parent]['hairpin'] += contribution
                    parent_dominant[parent].append((c_name, contribution))
                elif c_name == f"{p_name}_homo":
                    parent_contributions[parent]['homo'] += contribution
                    parent_dominant[parent].append((c_name, contribution))
                elif "_on_" in c_name:
                    parent_contributions[parent]['target'] += contribution
                else:
                    # C'est un complexe avec au moins une autre amorce
                    # Vérifions si TOUTES les amorces de ce complexe ont le MÊME parent
                    other_parents = set()
                    is_inter = False
                    
                    my_panel = primer_to_panel.get(p_name) if primer_to_panel else None
                    
                    for other_j, other_p in enumerate(primer_names):
                        if other_j != j and stoichiometry[i, other_j] > 0:
                            other_parent = get_parent(other_p)
                            other_parents.add(other_parent)
                            
                            if primer_to_panel:
                                other_panel = primer_to_panel.get(other_p)
                                if my_panel and other_panel and my_panel != other_panel:
                                    is_inter = True
                                    
                    # Règle : Si l'autre amorce a le même parent, c'est un homodimère pour le parent
                    if len(other_parents) == 1 and parent in other_parents:
                        parent_contributions[parent]['homo'] += contribution
                    else:
                        if is_inter:
                            parent_contributions[parent]['hetero_inter'] += contribution
                        else:
                            parent_contributions[parent]['hetero_intra'] += contribution
                            
                    parent_dominant[parent].append((c_name, contribution))
                    
        parent_totals[parent] += total
        
    for parent, total in parent_totals.items():
        if total > 0:
            dom_c = ""
            dom_f = 0.0
            if parent_dominant[parent]:
                # On peut regrouper les dominants de même nom ou garder le pire complexe individuel
                # Gardons le pire complexe absolu
                best_c, best_val = max(parent_dominant[parent], key=lambda x: x[1])
                dom_c = best_c
                dom_f = best_val / total
                
            fractions[parent] = PrimerFractions(
                free=parent_contributions[parent]['free'] / total,
                hairpin=parent_contributions[parent]['hairpin'] / total,
                homodimer=parent_contributions[parent]['homo'] / total,
                heterodimer_intra=parent_contributions[parent]['hetero_intra'] / total,
                heterodimer_inter=parent_contributions[parent]['hetero_inter'] / total,
                target_bound=parent_contributions[parent]['target'] / total,
                dominant_complex=dom_c,
                dominant_fraction=dom_f
            )
        else:
            fractions[parent] = PrimerFractions(0, 0, 0, 0, 0, 0)
            
    return fractions
