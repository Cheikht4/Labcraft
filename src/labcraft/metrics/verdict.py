from dataclasses import dataclass
from typing import List, Dict
from labcraft.metrics.fractions import PrimerFractions
from labcraft.metrics.risk import RiskItem

@dataclass
class PrimerIssue:
    primer_name: str
    target_site: str
    occupation: float
    cause: str
    is_critical: bool

@dataclass
class PanelVerdict:
    status: str # "OK", "WARNING", "FAILURE"
    issues: List[PrimerIssue]
    global_cause: str

def generate_verdict(
    fractions: Dict[str, PrimerFractions], 
    target_occupations: Dict[str, Dict[str, float]],
    risks: List[RiskItem]
) -> PanelVerdict:
    """
    Évalue le panel et identifie toutes les amorces en difficulté
    en nommant le mécanisme dominant réel pour chacune.
    """
    issues = []
    
    if target_occupations:
        # Mode avec cible(s)
        for target_id, occupations in target_occupations.items():
            for primer_name, f in fractions.items():
                site_name = f"{primer_name}_site"
                occ = occupations.get(site_name, 0.0)
                
                # Si le primer est censé cibler cette cible (ex: se termine par _A pour SynthA)
                # ou s'il est physiquement en difficulté
                if occ < 0.1 and not ("LF" in primer_name or "LB" in primer_name): # Moins de 10% d'occupation = amorce d'initiation en difficulté
                    # Filtrons pour ne signaler que les amorces qui SONT du panel de la cible
                    # (On suppose ici que F3_A cible SynthA, etc.)
                    target_suffix = target_id.replace("Synth", "")
                    if not primer_name.endswith(f"_{target_suffix}"):
                        continue
                    
                    is_crit = occ < 0.01
                    
                    # Déterminer la cause dominante
                    if f.free < 0.1:
                        dom_frac_pct = f.dominant_fraction * 100
                        if f.heterodimer_inter > f.homodimer and f.heterodimer_inter > f.heterodimer_intra:
                            cause = f"Hybridation croisée inter-jeux séquestrant {primer_name} à {dom_frac_pct:.1f}% dans {f.dominant_complex}."
                        else:
                            cause = f"Séquestration de {primer_name} à {dom_frac_pct:.1f}% dans {f.dominant_complex}."
                    else:
                        if site_name not in occupations:
                            cause = "Site absent de la cible ou amorce non appariée."
                        else:
                            cause = f"Inaccessibilité du site cible pour {primer_name} (barrière d'ouverture structurale ou thermodynamique défavorable)."
                        
                    issues.append(PrimerIssue(
                        primer_name=primer_name,
                        target_site=target_id,
                        occupation=occ,
                        cause=cause,
                        is_critical=is_crit
                    ))
    else:
        # Mode compétition pure sans cible
        for primer_name, f in fractions.items():
            if f.free < 0.20:
                dom_frac_pct = f.dominant_fraction * 100
                if f.heterodimer_inter > f.homodimer and f.heterodimer_inter > f.heterodimer_intra:
                    cause = f"Hybridation croisée inter-jeux séquestrant {primer_name} à {dom_frac_pct:.1f}% dans {f.dominant_complex}."
                else:
                    cause = f"Séquestration de {primer_name} à {dom_frac_pct:.1f}% dans {f.dominant_complex}."
                
                issues.append(PrimerIssue(
                    primer_name=primer_name,
                    target_site="Compétition Sans Cible",
                    occupation=f.free,
                    cause=cause,
                    is_critical=True
                ))
            
    # Détermination du statut global
    status = "OK"
    if target_occupations:
        global_cause = "Tous les sites d'initiation présentent une accessibilité favorable."
    else:
        global_cause = "Equilibre thermodynamique favorable."
    
    if any(issue.is_critical for issue in issues):
        status = "FAILURE"
        global_cause = "Échec critique sur une ou plusieurs cibles (voir détails)."
    elif issues:
        status = "WARNING"
        global_cause = "Accessibilité d'initiation réduite sur certains sites."
        
    # Surcharge du verdict si on a des dimères amplifiables massifs
    has_amplifiable = False
    for r in risks:
        if r.severity == 10.0 and r.concentration > 1e-9:
            has_amplifiable = True
            break
            
    if has_amplifiable:
        status = "FAILURE"
        global_cause = "Échec probable : Présence de dimères amplifiables massifs."
        
    return PanelVerdict(status, issues, global_cause)
