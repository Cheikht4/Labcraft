"""Non-sequential recommendations.

Recommandations qualitatives (non-séquentielles) pour l'optimisation.
"""
from typing import List, Dict, Any, Optional
from labcraft.metrics.verdict import PanelVerdict
from labcraft.metrics.risk import RiskItem

def generate_recommendations(verdict: PanelVerdict, risks: Optional[List[RiskItem]] = None) -> List[str]:
    """
    Génère des recommandations qualitatives basées sur le verdict ET les risques réels.
    
    La détection des dimères amplifiables se fonde sur les RiskItem (severity >= 10.0),
    pas sur une recherche de sous-chaîne dans les causes du verdict.
    Un panel portant des dimères amplifiables ne doit jamais être qualifié de sain.
    
    Args:
        verdict: Le verdict global du panel.
        risks: La liste des risques évalués (RiskItem).
        
    Returns:
        Une liste de recommandations en clair.
    """
    recs = []
    
    # Détection basée sur les risques réels (source de vérité)
    # Detection based on actual risks (source of truth)
    has_amplifiable = False
    has_blocking = False
    
    if risks:
        has_amplifiable = any(r.severity >= 10.0 and r.concentration > 1e-9 for r in risks)
        has_blocking = any(r.severity < 10.0 and r.concentration > 1e-9 for r in risks)
    
    # Fallback sur le verdict si les risques ne sont pas passés (rétrocompatibilité)
    # Fallback on verdict if risks are not passed (backward compatibility)
    if risks is None:
        has_amplifiable = any("amplifiable" in issue.cause.lower() for issue in verdict.issues)
        has_blocking = any("bloquant" in issue.cause.lower() for issue in verdict.issues)
    
    if has_amplifiable:
        recs.append(
            "Enzyme : Envisagez l'utilisation d'une polymérase WarmStart. "
            "Cela empêchera l'extension des dimères amplifiables à basse température "
            "lors de la préparation du mix."
        )
        recs.append(
            "Technologie de détection : Si les dimères amplifiables ne peuvent être éliminés, "
            "envisagez de remplacer les colorants intercalants (ex: SYBR Green, EvaGreen) "
            "par des sondes spécifiques (ex: TaqMan, LAMP-OSD) pour masquer le signal parasite."
        )
        
    if has_blocking:
        recs.append(
            "Asymétrie de concentration : Certains dimères bloquants séquestrent massivement "
            "des amorces. Envisagez de baisser la concentration de l'amorce la plus abondante "
            "impliquée dans ces dimères, ou d'augmenter celle de l'amorce séquestrée."
        )
        
    if not recs:
        recs.append("Aucune recommandation non-séquentielle critique. Le panel semble sain.")
        
    return recs
