import jinja2
import os
from typing import Dict, List
from labcraft.metrics.verdict import PanelVerdict
from labcraft.metrics.fractions import PrimerFractions
from labcraft.metrics.risk import RiskItem

def render_report(
    verdict: PanelVerdict,
    fractions: Dict[str, PrimerFractions],
    risks: List[RiskItem],
    metadata: dict
) -> str:
    """
    Rend le rapport HTML à partir du template Jinja2.
    """
    # Chemin vers le template (relatif à ce fichier)
    template_dir = os.path.join(os.path.dirname(__file__), "templates")
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(template_dir),
        undefined=jinja2.StrictUndefined,
        autoescape=True
    )
    
    # On ajoute un filtre pour formatter les pourcentages
    def format_pct(value):
        return f"{value * 100:.1f}%"
    env.filters['pct'] = format_pct
    
    template = env.get_template("report.html.j2")
    
    # Groupement des issues par cible
    issues_by_target = {}
    for issue in verdict.issues:
        if issue.target_site not in issues_by_target:
            issues_by_target[issue.target_site] = []
        issues_by_target[issue.target_site].append(issue)
        
    # Séparation des dimères amplifiables et bloquants pour un affichage clair
    amplifiable_dimers = [r for r in risks if r.severity >= 10.0]
    blocking_dimers = [r for r in risks if r.severity < 10.0]
    
    return template.render(
        verdict=verdict,
        issues_by_target=issues_by_target,
        fractions=fractions,
        amplifiable_dimers=amplifiable_dimers,
        blocking_dimers=blocking_dimers,
        metadata=metadata
    )
