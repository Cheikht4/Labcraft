from dataclasses import dataclass
from typing import List, Dict
from labcraft.metrics.fractions import PrimerFractions
from labcraft.metrics.risk import RiskItem
from labcraft.metrics.balance import calculate_balance_index, find_limiting_primer

@dataclass
class PanelVerdict:
    status: str # "OK", "WARNING", "FAILURE"
    dominant_cause: str
    balance_index: float
    limiting_primer: str

def generate_verdict(
    fractions: Dict[str, PrimerFractions], 
    risks: List[RiskItem],
    target_concentration_molar: float = 1e-15
) -> PanelVerdict:
    """
    Génère le verdict global d'un panel sur une cible donnée.
    """
    balance = calculate_balance_index(fractions)
    lim_primer, lim_bound = find_limiting_primer(fractions)
    
    status = "OK"
    cause = "Equilibre thermodynamique favorable."
    
    # 1. Vérification de l'accessibilité de la cible (Si un primer lié à la cible a une occupation < 1%)
    # lim_bound est la fraction de l'amorce limitante liée à la cible.
    # Mais attention : l'occupation de la cible est T_bound / T_total.
    # T_bound / T_total = lim_bound * P_total / T_total.
    # Pour faire simple : on va utiliser les fractions.
    # Si la fraction liée productive d'une amorce clé est excessivement faible (ex: < 1e-6)
    if lim_bound < 1e-6:
        status = "FAILURE"
        cause = f"Échec d'accessibilité ou hybridation impossible : {lim_primer} ne peut pas se lier à la cible."
        return PanelVerdict(status, cause, balance, lim_primer)
        
    # 2. Vérification des risques d'artefacts
    if risks:
        top_risk = risks[0]
        if top_risk.severity == 10.0:
            status = "FAILURE"
            cause = f"Échec probable : Dimère amplifiable massif détecté ({top_risk.complex_name})."
        elif top_risk.severity >= 1.0 and top_risk.concentration > 1e-7:
            if status != "FAILURE":
                status = "WARNING"
                cause = f"Risque de perte d'efficacité : compétition par hybridation croisée ou homodimère ({top_risk.complex_name})."
                
    # 3. Déséquilibre du panel
    if status == "OK" and balance > 100.0:
        status = "WARNING"
        cause = f"Déséquilibre thermodynamique sévère (Indice {balance:.1f})."
        
    return PanelVerdict(status, cause, balance, lim_primer)
