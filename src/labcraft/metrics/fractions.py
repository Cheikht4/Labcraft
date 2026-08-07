from dataclasses import dataclass
from typing import Dict, List
import numpy as np

@dataclass
class PrimerFractions:
    free: float
    hairpin: float
    homodimer: float
    heterodimer: float # avec des amorces du même set (intrapanel) ou cross-hybridation
    target_bound: float
    
    @property
    def sum(self) -> float:
        return self.free + self.hairpin + self.homodimer + self.heterodimer + self.target_bound

def compute_fractions(
    primer_names: List[str],
    complex_names: List[str],
    stoichiometry: np.ndarray,
    free_concentrations: np.ndarray,
    delta_g: np.ndarray,
    temp_celsius: float
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
    R = 0.00198720425864083
    RT = R * (273.15 + temp_celsius)
    
    u = np.log(free_concentrations)
    
    fractions = {}
    
    for j, p_name in enumerate(primer_names):
        if p_name.endswith("_site"):
            continue # On ne calcule pas les fractions de la cible
            
        f_free = 0.0
        f_hairpin = 0.0
        f_homo = 0.0
        f_hetero = 0.0
        f_target = 0.0
        
        # Total initial de l'amorce
        total = 0.0
        
        for i, c_name in enumerate(complex_names):
            coeff = stoichiometry[i, j]
            if coeff > 0:
                # Concentration du complexe = exp(-dg/RT + stoich . u)
                conc_complex = np.exp(-delta_g[i] / RT + np.dot(stoichiometry[i], u))
                
                # La contribution de ce complexe au total de l'amorce est coeff * conc_complex
                contribution = coeff * conc_complex
                total += contribution
                
                if c_name == f"{p_name}_free":
                    f_free += contribution
                elif c_name == f"{p_name}_hairpin":
                    f_hairpin += contribution
                elif c_name == f"{p_name}_homo":
                    f_homo += contribution
                elif "_on_" in c_name:
                    f_target += contribution
                else:
                    # Hétérodimère ou hybridation croisée
                    f_hetero += contribution
                    
        if total > 0:
            fractions[p_name] = PrimerFractions(
                free=f_free / total,
                hairpin=f_hairpin / total,
                homodimer=f_homo / total,
                heterodimer=f_hetero / total,
                target_bound=f_target / total
            )
        else:
            fractions[p_name] = PrimerFractions(0, 0, 0, 0, 0)
            
    return fractions
