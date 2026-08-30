"""Panel verdict generation / Génération du verdict de panel.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Dict, Optional
import numpy as np
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
        # Mode avec cible(s) : évaluation par amorce parente / locus
        # Mode with target(s): evaluation per parent primer / locus
        for target_id, occupations in target_occupations.items():
            # Regrouper les fractions par amorce parente
            parent_to_variants: Dict[str, List[str]] = {}
            for primer_name in fractions.keys():
                parent_name = primer_name.split('#')[0] if '#' in primer_name else primer_name
                parent_to_variants.setdefault(parent_name, []).append(primer_name)

            for parent_name, var_names in parent_to_variants.items():
                # Filtrer selon le panel d'appartenance si disponible
                if primer_to_panel is not None:
                    p_panel = primer_to_panel.get(parent_name) or primer_to_panel.get(var_names[0])
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

                # Représentant des fractions pour le parent (moyenne ou premier)
                f_list = [fractions[v] for v in var_names if v in fractions]
                mean_free = float(np.mean([f.free for f in f_list])) if f_list else 1.0
                dom_complex = f_list[0].dominant_complex if f_list else ""
                dom_frac_pct = float(np.mean([f.dominant_fraction for f in f_list])) * 100 if f_list else 0.0
                mean_inter = float(np.mean([f.heterodimer_inter for f in f_list])) if f_list else 0.0
                mean_intra = float(np.mean([f.heterodimer_intra for f in f_list])) if f_list else 0.0
                mean_homo = float(np.mean([f.homodimer for f in f_list])) if f_list else 0.0

                # Moins de 10% d'occupation = amorce d'initiation en difficulté
                if occ < 0.1 and parent_name not in loop_primer_parents:
                    is_crit = occ < 0.01

                    # Déterminer la cause dominante
                    if mean_free < 0.1:
                        if mean_inter > mean_homo and mean_inter > mean_intra:
                            cause = f"Hybridation croisée inter-jeux séquestrant {parent_name} à {dom_frac_pct:.1f}% dans {dom_complex}."
                        else:
                            cause = f"Séquestration de {parent_name} à {dom_frac_pct:.1f}% dans {dom_complex}."
                    else:
                        has_site = (site_name in occupations or any(k.startswith(f"{parent_name}@") or k.startswith(f"{parent_name}_") for k in occupations))
                        if not has_site or occ == 0.0:
                            cause = "Site absent de la cible ou amorce non appariée."
                        else:
                            if len(var_names) > 1:
                                cause = f"Inaccessibilité du site cible ou dilution dégénérée ({len(var_names)} variantes au panel, conc. utile fractionnée)."
                            else:
                                cause = f"Inaccessibilité du site cible pour {parent_name} (barrière d'ouverture structurale ou thermodynamique défavorable)."

                    issues.append(PrimerIssue(
                        primer_name=parent_name,
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
