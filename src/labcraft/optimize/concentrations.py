from typing import List, Dict, Any, Tuple
import copy
import dataclasses
import numpy as np

from labcraft.solver.types import EquilibriumProblem, ConvergenceError
from labcraft.solver.dual import solve_dual
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.metrics.balance import calculate_multiplex_balance
from labcraft.diagnostics.amplifiable_dimer import evaluate_pair_amplifiable

def get_dangerous_dimers(
    prob_template: EquilibriumProblem,
    complex_names: List[str],
    species_names: List[str],
    primers: List[PhysicalPrimer],
    temp_celsius: float, 
    backend: Any, 
    enzyme: Any
) -> set[str]:
    dangerous = set()
    primer_map = {p.name: p for p in primers}
    
    pair_cache = {}
    
    for c_idx, c_name in enumerate(complex_names):
        row = prob_template.stoichiometry[c_idx]
        is_dimer = True
        total_primers = 0
        consumed_primers = []
        
        for s_idx, stoich in enumerate(row):
            if stoich > 0:
                s_name = species_names[s_idx]
                if s_name.endswith("_site"):
                    is_dimer = False
                    break
                if s_name in primer_map:
                    total_primers += int(stoich)
                    consumed_primers.extend([primer_map[s_name]] * int(stoich))
                else:
                    is_dimer = False
                    break
                    
        if is_dimer and total_primers == 2 and len(consumed_primers) == 2:
            p1, p2 = consumed_primers
            key = tuple(sorted([p1.name, p2.name]))
            if key not in pair_cache:
                is_amp = evaluate_pair_amplifiable(p1, p2, backend, enzyme, temp_celsius)[0]
                pair_cache[key] = is_amp
                
            if pair_cache[key]:
                dangerous.add(c_name)
                
    return dangerous

def optimize_concentrations(
    prob_template: EquilibriumProblem,
    species_names: List[str],
    primers: List[PhysicalPrimer],
    target_dict: Dict[str, str],
    primer_to_panel: Dict[str, str],
    original_free_fractions: Dict[str, float],
    original_target_occupations: Dict[str, float],
    complex_names: List[str],
    temp_celsius: float,
    backend: Any,
    enzyme: Any = None,
    max_iter: int = 3,
    min_initiation_occupation: float = 0.01,
    weight_dangerous: float = 1e6, # Normalized assuming primers sum to ~1e-6 M
    weight_floor: float = 1e4,
    weight_cv: float = 1.0,
    weight_regularization: float = 1.0,
    min_improvement_ratio: float = 0.05,
    min_absolute_improvement: float = 1e-3, # Min fraction of total primer conc
    max_occ_drop_ratio: float = 0.20,
    min_tradeoff_factor: float = 2.0, # Gain relatif en dimères doit être >= N * Perte relative max en occupation
    min_dimer_reduction_ratio: float = 0.20 # Doit réduire les dimères d'au moins 20%
) -> List[Dict[str, Any]]:
    
    # Identifier les dimères dangereux
    dangerous_complexes = get_dangerous_dimers(
        prob_template, complex_names, species_names, primers, temp_celsius, backend, enzyme
    )
    
    # Indexer
    primer_idx = {p.name: species_names.index(p.name) for p in primers if p.name in species_names}
    
    # Trouver les complexes d'initiation
    initiation_complexes = {}
    for i, c in enumerate(complex_names):
        for p in primers:
            if f"{p.name}_on_{p.name}_site" == c:
                initiation_complexes[p.name] = i
                
    dangerous_indices = [i for i, c in enumerate(complex_names) if c in dangerous_complexes]
    
    # Prépare mapping panel inverse
    panels = set(primer_to_panel.values())
    
    # Initialize iteration counter
    rejected_count = 0
    
    def evaluate(c_candidats: np.ndarray) -> Tuple[float, dict]:
        nonlocal rejected_count
        prob = dataclasses.replace(prob_template, total_concentrations=c_candidats)
        try:
            res = solve_dual(prob)
        except ConvergenceError:
            rejected_count += 1
            return float('inf'), {}
            
        score = 0.0
        
        # 1. Dimères dangereux (normalisé)
        # sum of total concentrations is around ~10 uM -> 1e-5 M. We divide by 1e-6 to get a number around 1.0
        total_primer_conc = np.sum(prob_template.total_concentrations)
        dang_frac = np.sum(res.concentrations[dangerous_indices])
        normalized_dang = dang_frac / total_primer_conc if total_primer_conc > 0 else 0.0
        score += normalized_dang * weight_dangerous
        
        # 2. Plancher d'amorçage
        initiations = {}
        target_occs_for_balance = {panel: {} for panel in panels}
        free_fracs_for_balance = {}
        
        floor_penalty = 0.0
        reg_penalty = 0.0
        
        for p in primers:
            p_name = p.name
            if p_name in primer_idx:
                idx = primer_idx[p_name]
                free_c = res.free_concentrations[idx]
                tot_c = c_candidats[idx]
                free_fracs_for_balance[p_name] = free_c / tot_c if tot_c > 0 else 0.0
            
            if p_name in initiation_complexes:
                c_idx = initiation_complexes[p_name]
                site_name = f"{p_name}_site"
                try:
                    site_idx = species_names.index(site_name)
                    site_conc = c_candidats[site_idx]
                except ValueError:
                    site_conc = 0.0
                    
                occ = res.concentrations[c_idx] / site_conc if site_conc > 0 else 0.0
                initiations[p_name] = occ
                
                panel = primer_to_panel.get(p_name)
                if panel:
                    target_occs_for_balance[panel][site_name] = occ
                    
                if occ < min_initiation_occupation:
                    deficit = (min_initiation_occupation - occ) / min_initiation_occupation
                    floor_penalty += (deficit ** 2)
                    
            if p_name in primer_idx:
                c_origine = prob_template.total_concentrations[primer_idx[p_name]]
                c_candidat = c_candidats[primer_idx[p_name]]
                if c_origine > 0:
                    delta_rel = abs(c_candidat - c_origine) / c_origine
                    reg_penalty += delta_rel
                    
        score += floor_penalty * weight_floor
        score += reg_penalty * weight_regularization
                    
        # 3. CV Balance inter-panels
        cv = 0.0
        panel_summaries, calc_cv = calculate_multiplex_balance(
            primer_to_panel, target_occs_for_balance, free_fracs_for_balance
        )
        if calc_cv is not None:
            cv = calc_cv
            score += cv * weight_cv
            
        return score, {
            'dangerous': dang_frac, 
            'normalized_dangerous': normalized_dang,
            'initiations': initiations, 
            'cv': cv,
            'score_terms': {
                'dangerous': normalized_dang * weight_dangerous,
                'floor': floor_penalty * weight_floor,
                'regularization': reg_penalty * weight_regularization,
                'cv': cv * weight_cv
            }
        }
        
    c_current = np.copy(prob_template.total_concentrations)
    
    # Plages de concentrations selon le rôle
    ranges = {
        PrimerRole.FIP: np.arange(0.8e-6, 2.1e-6, 0.2e-6),
        PrimerRole.BIP: np.arange(0.8e-6, 2.1e-6, 0.2e-6),
        PrimerRole.LF: np.arange(0.2e-6, 1.1e-6, 0.2e-6),
        PrimerRole.LB: np.arange(0.2e-6, 1.1e-6, 0.2e-6),
        PrimerRole.F3: np.arange(0.1e-6, 0.5e-6, 0.1e-6),
        PrimerRole.B3: np.arange(0.1e-6, 0.5e-6, 0.1e-6),
    }
    
    best_score, best_metrics = evaluate(c_current)
    stopped_on_max = True
    
    # Calculer occupations initiales pour le garde-fou 20%
    occ_initial = best_metrics['initiations']
    
    for iteration in range(max_iter):
        changed = False
        for p in primers:
            if p.name not in primer_idx or p.role not in ranges:
                continue
                
            idx = primer_idx[p.name]
            best_c = c_current[idx]
            
            panel = primer_to_panel.get(p.name)
            
            for c_test in ranges[p.role]:
                # Contrainte de hiérarchie pour le même panel
                valid_hierarchy = True
                
                # Check bounds based on role against current concentrations
                if p.role in (PrimerRole.FIP, PrimerRole.BIP):
                    for other_p in primers:
                        if primer_to_panel.get(other_p.name) == panel and other_p.name in primer_idx:
                            other_c = c_current[primer_idx[other_p.name]]
                            if other_p.role in (PrimerRole.LF, PrimerRole.LB, PrimerRole.F3, PrimerRole.B3) and c_test < other_c:
                                valid_hierarchy = False
                                break
                                
                if p.role in (PrimerRole.LF, PrimerRole.LB):
                    for other_p in primers:
                        if primer_to_panel.get(other_p.name) == panel and other_p.name in primer_idx:
                            other_c = c_current[primer_idx[other_p.name]]
                            if other_p.role in (PrimerRole.FIP, PrimerRole.BIP) and c_test > other_c:
                                valid_hierarchy = False
                                break
                            if other_p.role in (PrimerRole.F3, PrimerRole.B3) and c_test < other_c:
                                valid_hierarchy = False
                                break
                                
                if p.role in (PrimerRole.F3, PrimerRole.B3):
                    for other_p in primers:
                        if primer_to_panel.get(other_p.name) == panel and other_p.name in primer_idx:
                            other_c = c_current[primer_idx[other_p.name]]
                            if other_p.role in (PrimerRole.FIP, PrimerRole.BIP, PrimerRole.LF, PrimerRole.LB) and c_test > other_c:
                                valid_hierarchy = False
                                break
                                
                if not valid_hierarchy:
                    continue
                    
                c_candidats = np.copy(c_current)
                c_candidats[idx] = c_test
                score, metrics = evaluate(c_candidats)
                
                # 1. Calcul du gain relatif sur les dimères
                dang_initial = best_metrics['normalized_dangerous']
                dang_candidat = metrics['normalized_dangerous']
                rel_gain_dimers = (dang_initial - dang_candidat) / dang_initial if dang_initial > 0 else 0.0
                
                # 2. Calcul du coût (perte relative maximale d'occupation)
                # Ce calcul itère sur toutes les amorces (F3, B3, FIP, BIP, ET LF, LB)
                # Cela assure une protection explicite et proportionnée des amorces de boucle,
                # dont la baisse de concentration impacterait fortement leur occupation.
                max_rel_occ_drop = 0.0
                for p_name, initial_occ in occ_initial.items():
                    candidat_occ = metrics['initiations'].get(p_name, 0.0)
                    if initial_occ > 0 and candidat_occ < initial_occ:
                        drop = (initial_occ - candidat_occ) / initial_occ
                        if drop > max_rel_occ_drop:
                            max_rel_occ_drop = drop
                            
                # Garde-fou absolu d'occupation (ex: 20%)
                if max_rel_occ_drop > max_occ_drop_ratio:
                    # print(f"REJECTED {p.name} {c_test}: max_rel_occ_drop {max_rel_occ_drop} > {max_occ_drop_ratio}")
                    continue
                    
                # Critère de compromis : le gain doit valoir le coût
                if rel_gain_dimers > 0 and max_rel_occ_drop > 0:
                    if rel_gain_dimers < min_tradeoff_factor * max_rel_occ_drop:
                        # print(f"REJECTED {p.name} {c_test}: gain {rel_gain_dimers} < {min_tradeoff_factor} * {max_rel_occ_drop}")
                        continue # Rejeté : compromis défavorable
                        
                # Exiger une baisse significative des dimères (ex: 20%)
                if dang_initial > 0 and rel_gain_dimers < min_dimer_reduction_ratio:
                    # print(f"REJECTED {p.name} {c_test}: gain {rel_gain_dimers} < {min_dimer_reduction_ratio}")
                    continue # Rejeté : baisse trop faible
                
                # Seuil d'amélioration absolue et relative du score global
                if score < best_score:
                    rel_improvement = (best_score - score) / best_score if best_score > 0 else 0
                    abs_improvement = best_score - score
                    
                    if rel_improvement > min_improvement_ratio and abs_improvement > min_absolute_improvement:
                        best_score = score
                        best_metrics = metrics
                        best_c = c_test
                        changed = True
                    
            c_current[idx] = best_c
            
        if not changed:
            stopped_on_max = False
            break
            
    # Résultat
    results = []
    final_score, final_metrics = evaluate(c_current)
    
    if stopped_on_max:
        # Journaliser (on pourrait utiliser logging, ici on stocke juste l'info)
        print(f"Optimize concentrations stopped on max_iter={max_iter}")
        
    for p in primers:
        if p.name in primer_idx:
            idx = primer_idx[p.name]
            orig_c = prob_template.total_concentrations[idx]
            new_c = c_current[idx]
            
            if abs(orig_c - new_c) > 1e-9:
                results.append({
                    "primer_name": p.name,
                    "original_conc": orig_c,
                    "suggested_conc": new_c,
                    "reason": "Réduction des dimères dangereux ou amélioration de l'équilibre multiplexe."
                })
                
    return results
