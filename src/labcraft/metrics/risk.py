from dataclasses import dataclass
from typing import List

@dataclass
class RiskItem:
    complex_name: str
    severity: float
    concentration: float
    description: str

def evaluate_risks(
    complex_names: List[str], 
    concentrations: List[float], 
    amplifiable_flags: List[bool],
    is_warm_start: bool = False
) -> List[RiskItem]:
    """
    Évalue le risque d'artefacts basé sur la concentration et l'extensibilité.
    """
    risks = []
    
    for c_name, conc, is_amp in zip(complex_names, concentrations, amplifiable_flags):
        if "_free" in c_name or "_on_" in c_name:
            continue # Pas de risque pour les espèces libres ou liées à la cible
            
        severity = 0.0
        desc = ""
        
        if is_amp:
            severity = 10.0 # Dimère amplifiable = risque majeur
            desc = "Dimère amplifiable (3' extensible)"
            if is_warm_start:
                # La WarmStart ne change pas la thermodynamique à 65°C, 
                # mais elle évite l'amplification à froid. Le risque subsiste pendant l'amplification.
                desc += " - Note: La WarmStart mitige le risque de formation pré-incubation."
        else:
            # Dimère bloquant (appauvrit l'amorce mais ne s'amplifie pas)
            severity = 1.0 
            desc = "Dimère bloquant (non extensible)"
            if "homo" in c_name:
                desc = "Homodimère bloquant"
                
        if conc > 1e-9: # Ne reporter que les artefacts ayant une concentration > 1 nM
            risks.append(RiskItem(c_name, severity, conc, desc))
            
    # Trier par risque global (sévérité * concentration)
    risks.sort(key=lambda r: r.severity * r.concentration, reverse=True)
    return risks

