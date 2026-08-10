from typing import Dict, Tuple, List, Optional
import numpy as np

def calculate_multiplex_balance(
    primer_to_panel: Dict[str, str],
    target_occupations: Dict[str, Dict[str, float]],
    free_fractions: Dict[str, float]
) -> Tuple[Dict[str, dict], Optional[float]]:
    """
    Calcule la balance des panels en multiplexe.
    Utilise l'occupation côté SITE des amorces d'INITIATION (F3, B3, FIP, BIP).
    Les amorces de boucle (LF, LB) sont ignorées pour la balance d'initiation.
    
    Returns:
        panel_summaries: Dict avec pour chaque panel:
            - mean_occupation: occupation moyenne d'initiation
            - min_occupation: occupation minimale
            - limiting_primer: nom de l'amorce limitante
            - mean_free: fraction libre moyenne des amorces du panel
        cv: Coefficient de variation inter-panels (std/mean), ou None si < 2 panels.
            Proche de 0 = panels équilibrés ; croissant = déséquilibre.
    """
    panels = set(primer_to_panel.values())
    
    panel_summaries = {}
    panel_means = []
    
    for panel in panels:
        # Extraire les amorces d'initiation de ce panel et leurs occupations
        # On peut trouver l'occupation via target_occupations[panel][site_name]
        initiation_occs = {}
        panel_free_fracs = []
        
        for p_name, p_panel in primer_to_panel.items():
            if p_panel == panel:
                # Résoudre le parent pour le lookup des fractions libres
                # Resolve parent for free fractions lookup
                free_val = free_fractions.get(p_name, 0.0)
                if free_val == 0.0 and '#' in p_name:
                    free_val = free_fractions.get(p_name.split('#')[0], 0.0)
                panel_free_fracs.append(free_val)
                # Vérifier si c'est une amorce d'initiation
                # (Les boucles ont souvent LF ou LB dans le nom)
                if "LF" not in p_name and "LB" not in p_name:
                    # Résoudre le parent pour le lookup du site cible
                    # Resolve parent for target site lookup
                    parent_name = p_name.split('#')[0] if '#' in p_name else p_name
                    site_name = f"{parent_name}_site"
                    occ = target_occupations.get(panel, {}).get(site_name, 0.0)
                    # Éviter les doublons : on ne garde qu'une entrée par parent
                    # Avoid duplicates: keep only one entry per parent
                    if parent_name not in initiation_occs:
                        initiation_occs[parent_name] = occ
                    
        if initiation_occs:
            mean_occ = np.mean(list(initiation_occs.values()))
            min_occ = min(initiation_occs.values())
            limiting_primer = min(initiation_occs, key=initiation_occs.get)
        else:
            mean_occ = 0.0
            min_occ = 0.0
            limiting_primer = "N/A"
            
        mean_free = np.mean(panel_free_fracs) if panel_free_fracs else 0.0
        
        panel_summaries[panel] = {
            "mean_occupation": mean_occ,
            "min_occupation": min_occ,
            "limiting_primer": limiting_primer,
            "mean_free": mean_free
        }
        
        panel_means.append(mean_occ)
        
    cv = None
    if len(panel_means) > 1:
        mean_all = np.mean(panel_means)
        if mean_all > 0:
            std_all = np.std(panel_means)
            cv = std_all / mean_all
        else:
            cv = 0.0
            
    return panel_summaries, cv
