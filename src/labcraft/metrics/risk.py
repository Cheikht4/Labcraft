from dataclasses import dataclass, field
from typing import List, Optional

@dataclass
class RiskItem:
    complex_name: str
    concentration: float
    severity: float 
    description: str
    
    # Nouvelles métriques pour les fiches (Jalon 6.1)
    seq_a: str = ""
    seq_b: str = ""
    structure: str = ""
    delta_g: float = 0.0
    delta_g_3p: float = 0.0
    alignment_ascii: str = ""
    alignment_columns: list = field(default_factory=list)
    arrow_metrics: dict = field(default_factory=dict)
    is_blocked_veto: bool = False

def evaluate_risks(
    complex_names: List[str],  
    concentrations: List[float], 
    amplifiable_flags: List[bool],
    is_warm_start: bool = True,
    dimer_details: Optional[List[dict]] = None
) -> List[RiskItem]:
    """
    Évalue le risque d'artefacts basé sur la concentration et l'extensibilité.
    """
    risks = []
    
    for i, (c_name, conc, is_amp) in enumerate(zip(complex_names, concentrations, amplifiable_flags)):
        if "_free" in c_name or "_on_" in c_name:
            continue # Pas de risque pour les espèces libres ou liées à la cible
            
        severity = 0.0
        desc = ""
        
        details = dimer_details[i] if dimer_details and i < len(dimer_details) else {}

        if is_amp:
            # Sévérité pondérée par la marge d'amplifiabilité (ΔG 3' end) et la concentration
            margin = max(0.0, - (details.get("delta_g_3p", 0.0) + 3.0)) # ex: dg_3p=-5 -> margin = 2
            severity = 10.0 + margin
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
        
        if is_amp or conc > 1e-9:
            risks.append(RiskItem(
                complex_name=c_name,
                concentration=conc,
                severity=severity,
                description=desc,
                seq_a=details.get("seq_a", ""),
                seq_b=details.get("seq_b", ""),
                structure=details.get("structure", ""),
                delta_g=details.get("delta_g", 0.0),
                delta_g_3p=details.get("delta_g_3p", 0.0),
                alignment_ascii=details.get("alignment", ""),
                alignment_columns=details.get("alignment_columns", []),
                arrow_metrics=details.get("arrow_metrics", {}),
                is_blocked_veto=details.get("is_blocked_veto", False)
            ))
            
    # Trier par risque global (sévérité * concentration)
    risks.sort(key=lambda r: r.severity * r.concentration, reverse=True)
    return risks
