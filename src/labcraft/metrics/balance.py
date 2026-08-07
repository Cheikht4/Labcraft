from typing import Dict, Tuple
from labcraft.metrics.fractions import PrimerFractions

def calculate_balance_index(fractions: Dict[str, PrimerFractions]) -> float:
    """
    Calcule l'indice de déséquilibre du panel.
    C'est le ratio de la variance sur la moyenne de l'occupation des cibles,
    ou plus simplement le ratio Max / Min de la fraction target_bound.
    """
    bounds = [f.target_bound for f in fractions.values() if f.target_bound > 0]
    if not bounds:
        return 0.0
    return max(bounds) / min(bounds)

def find_limiting_primer(fractions: Dict[str, PrimerFractions]) -> Tuple[str, float]:
    """
    Identifie l'amorce limitante (celle ayant la plus faible fraction liée à la cible).
    """
    # Filtre sur les amorces qui ont vocation à se lier à la cible (FIP, BIP, F3, B3, LF, LB)
    # Dans ce démonstrateur, on regarde toutes les amorces.
    limiting_name = min(fractions, key=lambda k: fractions[k].target_bound)
    return limiting_name, fractions[limiting_name].target_bound

