"""Panel verdict generation / Génération du verdict de panel.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Dict, Optional
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
    risks: List[RiskItem],
    loop_primer_parents: set[str] = None,
    primer_to_panel: Optional[Dict[str, str]] = None
) -> PanelVerdict:
    """
    Évalue le panel et identifie toutes les amorces en difficulté
    en nommant le mécanisme dominant réel pour chacune.
    """
    if loop_primer_parents is None:
        loop_primer_parents = set()
        
    issues = []
    
    if target_occupations:
        # Mode avec cible(s)
        for target_id, occupations in target_occupations.items():
            for primer_name, f in fractions.items():
                # Résoudre le parent pour le lookup du site cible
                # Resolve parent for target site lookup
                parent_name = primer_name.split('#')[0] if '#' in primer_name else primer_name
                
                # Filtrer selon le panel d'appartenance si disponible
                if primer_to_panel is not None:
                    p_panel = primer_to_panel.get(primer_name) or primer_to_panel.get(parent_name)
                    if p_panel:
                        target_clean = target_id.replace("Synth", "")
                        panel_clean = p_panel.replace("Synth", "")
                        if target_id != p_panel and target_clean != panel_clean:
                            continue
                else:
                    target_suffix = target_id.replace("Synth", "")
                    if not (parent_name.endswith(f"_{target_suffix}") or f"_{target_suffix}_" in parent_name or parent_name.startswith(f"{target_suffix}_")):
                        if len(target_occupations) > 1:
                            continue

                site_name = f"{parent_name}_site"
                occ = occupations.get(site_name)
                if occ is None:
                    for k, v in occupations.items():
                        if k == site_name or k.startswith(f"{parent_name}@") or k.startswith(f"{parent_name}_"):
                            occ = v
                            break
                if occ is None:
                    occ = 0.0
                
                # Moins de 10% d'occupation = amorce d'initiation en difficulté
                if occ < 0.1 and parent_name not in loop_primer_parents:
                    is_crit = occ < 0.01
                    
                    # Déterminer la cause dominante
                    if f.free < 0.1:
                        dom_frac_pct = f.dominant_fraction * 100
                        if f.heterodimer_inter > f.homodimer and f.heterodimer_inter > f.heterodimer_intra:
                            cause = f"Hybridation croisée inter-jeux séquestrant {primer_name} à {dom_frac_pct:.1f}% dans {f.dominant_complex}."
                        else:
                            cause = f"Séquestration de {primer_name} à {dom_frac_pct:.1f}% dans {f.dominant_complex}."
                    else:
                        if site_name not in occupations and not any(k.startswith(f"{parent_name}@") or k.startswith(f"{parent_name}_") for k in occupations):
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
                is_crit = f.free < 0.05
                
                if f.heterodimer_inter > f.homodimer and f.heterodimer_inter > f.heterodimer_intra:
                    cause = f"Hybridation croisée inter-jeux séquestrant {primer_name} à {dom_frac_pct:.1f}% dans {f.dominant_complex}."
                else:
                    cause = f"Séquestration de {primer_name} à {dom_frac_pct:.1f}% dans {f.dominant_complex}."
                    
                issues.append(PrimerIssue(
                    primer_name=primer_name,
                    target_site="Compétition Sans Cible",
                    occupation=f.free,
                    cause=cause,
                    is_critical=is_crit
                ))
                
    # Évaluation globale
    has_critical_dimers = any(r.is_amplifiable and r.severity >= 10.0 for r in risks)
    has_critical_issues = any(i.is_critical for i in issues)
    
    if has_critical_dimers or has_critical_issues:
        status = "FAILURE"
        global_cause = "Présence de dimères d'amorces amplifiables majeurs ou séquestration critique interdisant l'amplification."
    elif len(issues) > 0 or any(r.severity >= 5.0 for r in risks):
        status = "WARNING"
        global_cause = "Compétitions ou pertes d'accessibilité modérées identifiées."
    else:
        status = "OK"
        global_cause = "Aucun risque majeur d'interaction parasite identifié."
        
    return PanelVerdict(status=status, issues=issues, global_cause=global_cause)
