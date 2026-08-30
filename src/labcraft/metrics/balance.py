"""Multiplex balance metrics / Métriques d'équilibre et de balance multiplexe.
"""
from typing import Dict, Tuple, List, Optional
import numpy as np

def calculate_multiplex_balance(
    primer_to_panel: Dict[str, str],
    target_occupations: Dict[str, Dict[str, float]],
    free_fractions: Dict[str, float],
    loop_primer_parents: set[str] = None
) -> Tuple[Dict[str, dict], Optional[float]]:
    """
    Calcule la balance des panels en multiplexe.
    Utilise l'occupation côté SITE des amorces d'INITIATION (F3, B3, FIP, BIP).
    Les amorces de boucle (LF, LB) sont ignorées pour la balance d'initiation.
    
    Returns:
        panel_summaries: Dict avec pour chaque panel:
            - mean_occupation: occupation moyenne d'initiation (ou None si non analysé)
            - min_occupation: occupation minimale (ou None si non analysé)
            - limiting_primer: nom de l'amorce limitante
            - mean_free: fraction libre moyenne des amorces du panel
        cv: Coefficient de variation inter-panels (std/mean), ou None si < 2 panels ou cibles non analysées.
    """
    if loop_primer_parents is None:
        loop_primer_parents = set()
        
    panels = set(primer_to_panel.values())
    
    panel_summaries = {}
    panel_means = []
    
    for panel in sorted(panels):
        # Vérifier si la cible/panel a été analysée dans target_occupations
        target_occs_for_panel = target_occupations.get(panel)
        # Si panel porte 'SynthA' ou 'A'
        if target_occs_for_panel is None:
            for t_k, t_v in target_occupations.items():
                if t_k == panel or t_k.replace("Synth", "") == panel.replace("Synth", ""):
                    target_occs_for_panel = t_v
                    break

        initiation_occs = {}
        panel_free_fracs = []
        
        for p_name, p_panel in primer_to_panel.items():
            if p_panel == panel:
                free_val = free_fractions.get(p_name, 0.0)
                if free_val == 0.0 and '#' in p_name:
                    free_val = free_fractions.get(p_name.split('#')[0], 0.0)
                panel_free_fracs.append(free_val)
                
                parent_name = p_name.split('#')[0] if '#' in p_name else p_name
                
                if parent_name not in loop_primer_parents:
                    site_name = f"{parent_name}_site"
                    occ = None
                    if target_occs_for_panel is not None:
                        occ = target_occs_for_panel.get(site_name)
                        if occ is None:
                            for k, v in target_occs_for_panel.items():
                                if k == site_name or k.startswith(f"{parent_name}@") or k.startswith(f"{parent_name}_"):
                                    occ = v
                                    break

                    if occ is not None and parent_name not in initiation_occs:
                        initiation_occs[parent_name] = occ
                    
        mean_free = float(np.mean(panel_free_fracs)) if panel_free_fracs else 0.0

        if target_occs_for_panel is None or not initiation_occs:
            # Panel ou cible non analysé
            panel_summaries[panel] = {
                "mean_occupation": None,
                "min_occupation": None,
                "limiting_primer": "Non analysé",
                "mean_free": mean_free
            }
        else:
            mean_occ = float(np.mean(list(initiation_occs.values())))
            min_occ = float(min(initiation_occs.values()))
            limiting_primer = min(initiation_occs, key=initiation_occs.get)
            
            panel_summaries[panel] = {
                "mean_occupation": mean_occ,
                "min_occupation": min_occ,
                "limiting_primer": limiting_primer,
                "mean_free": mean_free
            }
            panel_means.append(mean_occ)
        
    cv = None
    if len(panel_means) == len(panels) and len(panel_means) > 1:
        mean_all = float(np.mean(panel_means))
        if mean_all > 0:
            std_all = float(np.std(panel_means))
            cv = std_all / mean_all
        else:
            cv = 0.0
            
    return panel_summaries, cv
