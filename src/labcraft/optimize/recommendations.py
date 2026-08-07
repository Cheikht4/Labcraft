"""Non-sequential recommendations.

Recommandations qualitatives (non-séquentielles) pour l'optimisation.
"""
from typing import List, Dict, Any
from labcraft.metrics.verdict import PanelVerdict

def generate_recommendations(verdict: PanelVerdict) -> List[str]:
    """
    Génère des recommandations qualitatives basées sur le verdict.
    
    Args:
        verdict: Le verdict global du panel.
        
    Returns:
        Une liste de recommandations en clair.
    """
    recs = []
    
    # Analyse des problèmes
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
