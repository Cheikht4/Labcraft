from typing import List, Dict, Any, Tuple, Union
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
    enzyme: Any,
    **backend_kwargs
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
                is_amp = evaluate_pair_amplifiable(p1, p2, backend, enzyme, temp_celsius, **backend_kwargs)[0]
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
    min_dimer_reduction_ratio: float = 0.20, # Doit réduire les dimères d'au moins 20%
    **backend_kwargs
) -> Union[Dict[str, Any], List[Any]]:
    
    # Identifier les dimères dangereux
    dangerous_complexes = get_dangerous_dimers(
        prob_template, complex_names, species_names, primers, temp_celsius, backend, enzyme, **backend_kwargs
    )
    
    # Indexer les amorces et regrouper par amorce parente
    # Index primers and group by parent primer
    primer_idx = {p.name: species_names.index(p.name) for p in primers if p.name in species_names}
    parent_to_primers: Dict[str, List[PhysicalPrimer]] = {}
    for p in primers:
        if p.name in primer_idx:
            par = p.parent_name or (p.name.split('#')[0] if '#' in p.name else p.name)
            parent_to_primers.setdefault(par, []).append(p)

    # Trouver les complexes d'initiation : uniquement pour les amorces non dégénérées (sans '#')
    # ayant leur propre site exact `f"{p.name}_on_{p.name}_site"`, conformément au modèle d'origine.
    # Find initiation complexes: only for non-degenerate primers (without '#')
    # having their own exact site `f"{p.name}_on_{p.name}_site"`, matching the original model.
    parent_initiation_info: Dict[str, Tuple[int, List[int]]] = {}
    for par_name, p_list in parent_to_primers.items():
        if len(p_list) == 1:
            p = p_list[0]
            c_name = f"{p.name}_on_{p.name}_site"
            site_name = f"{p.name}_site"
            if c_name in complex_names and site_name in species_names:
                parent_initiation_info[par_name] = (
                    species_names.index(site_name),
                    [complex_names.index(c_name)]
                )
                
    dangerous_indices = [i for i, c in enumerate(complex_names) if c in dangerous_complexes]
    
    # Prépare mapping panel inverse
    panels = set(primer_to_panel.values())
    
    def evaluate(c_candidats: np.ndarray) -> Tuple[float, dict]:
        prob = dataclasses.replace(prob_template, total_concentrations=c_candidats)
        try:
            res = solve_dual(prob)
        except ConvergenceError:
            return float('inf'), {}
            
        score = 0.0
        
        # 1. Dimères dangereux (normalisé)
        total_primer_conc = np.sum(prob_template.total_concentrations)
        dang_frac = float(np.sum(res.concentrations[dangerous_indices]))
        normalized_dang = dang_frac / total_primer_conc if total_primer_conc > 0 else 0.0
        score += normalized_dang * weight_dangerous
        
        # 2. Plancher d'amorçage (évalué par amorce parente / locus)
        # 2. Initiation floor (evaluated per parent primer / locus)
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

        for par_name, p_list in parent_to_primers.items():
            if par_name in parent_initiation_info:
                site_idx, c_indices = parent_initiation_info[par_name]
                site_conc = c_candidats[site_idx]
                bound_conc = float(np.sum(res.concentrations[c_indices]))
                occ = bound_conc / site_conc if site_conc > 0 else 0.0
                initiations[par_name] = occ
                
                panel = primer_to_panel.get(par_name) or primer_to_panel.get(p_list[0].name)
                if panel and panel in target_occs_for_balance:
                    target_occs_for_balance[panel][f"{par_name}_site"] = occ
                    
                if occ < min_initiation_occupation:
                    deficit = (min_initiation_occupation - occ) / min_initiation_occupation
                    floor_penalty += (deficit ** 2)

            c_origine_par = sum(prob_template.total_concentrations[primer_idx[p.name]] for p in p_list)
            c_candidat_par = sum(c_candidats[primer_idx[p.name]] for p in p_list)
            if c_origine_par > 0:
                delta_rel = abs(c_candidat_par - c_origine_par) / c_origine_par
                reg_penalty += delta_rel
                    
        score += floor_penalty * weight_floor
        score += reg_penalty * weight_regularization
                    
        # 3. CV Balance inter-panels
        cv = 0.0
        _, calc_cv = calculate_multiplex_balance(
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
    
    # Plages de concentrations totales selon le rôle (au niveau de l'amorce parente)
    # Total concentration ranges by role (at the parent primer level)
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
    
    # État initial de référence pour le calcul cumulé du gain et du compromis
    # Initial reference state for cumulative gain and tradeoff calculation
    occ_initial = best_metrics['initiations']
    initial_metrics = best_metrics
    
    for _iteration in range(max_iter):
        changed = False
        for par_name, p_list in parent_to_primers.items():
            role = p_list[0].role
            if role not in ranges:
                continue
                
            n_vars = len(p_list)
            par_indices = [primer_idx[p.name] for p in p_list]
            best_c_total = float(sum(c_current[idx] for idx in par_indices))
            
            panel = primer_to_panel.get(par_name) or primer_to_panel.get(p_list[0].name)
            
            for c_test in ranges[role]:
                # Contrainte de hiérarchie pour le même panel (comparée aux concentrations totales par parent)
                valid_hierarchy = True
                
                for other_par, other_list in parent_to_primers.items():
                    if other_par == par_name:
                        continue
                    other_panel = primer_to_panel.get(other_par) or primer_to_panel.get(other_list[0].name)
                    if other_panel != panel:
                        continue
                    other_role = other_list[0].role
                    other_c_total = float(sum(c_current[primer_idx[op.name]] for op in other_list))
                    
                    if role in (PrimerRole.FIP, PrimerRole.BIP):
                        if other_role in (PrimerRole.LF, PrimerRole.LB, PrimerRole.F3, PrimerRole.B3) and c_test < other_c_total - 1e-12:
                            valid_hierarchy = False
                            break
                    elif role in (PrimerRole.LF, PrimerRole.LB):
                        if other_role in (PrimerRole.FIP, PrimerRole.BIP) and c_test > other_c_total + 1e-12:
                            valid_hierarchy = False
                            break
                        if other_role in (PrimerRole.F3, PrimerRole.B3) and c_test < other_c_total - 1e-12:
                            valid_hierarchy = False
                            break
                    elif role in (PrimerRole.F3, PrimerRole.B3):
                        if other_role in (PrimerRole.FIP, PrimerRole.BIP, PrimerRole.LF, PrimerRole.LB) and c_test > other_c_total + 1e-12:
                            valid_hierarchy = False
                            break
                                
                if not valid_hierarchy:
                    continue
                    
                # Tous les variants d'une même amorce parente bougent ensemble (équirépartition)
                # All variants of the same parent primer move together (equal split)
                c_candidats = np.copy(c_current)
                for idx in par_indices:
                    c_candidats[idx] = c_test / n_vars
                score, metrics = evaluate(c_candidats)
                if not metrics:
                    continue
                
                # 1. Calcul du gain relatif cumulé sur les dimères par rapport à l'état INITIAL
                # 1. Cumulative relative dimer gain wrt INITIAL state
                dang_initial = initial_metrics['dangerous']
                dang_candidat = metrics['dangerous']
                rel_gain_dimers = (dang_initial - dang_candidat) / dang_initial if dang_initial > 0 else 0.0
                
                # 2. Calcul du coût cumulé (perte relative maximale d'occupation par rapport à l'état INITIAL)
                # 2. Cumulative cost (max relative occupation drop wrt INITIAL state)
                max_rel_occ_drop = 0.0
                for p_key, initial_occ in occ_initial.items():
                    candidat_occ = metrics['initiations'].get(p_key, 0.0)
                    if initial_occ > 0 and candidat_occ < initial_occ:
                        drop = (initial_occ - candidat_occ) / initial_occ
                        if drop > max_rel_occ_drop:
                            max_rel_occ_drop = drop
                            
                # Garde-fou absolu d'occupation (ex: 20%)
                if max_rel_occ_drop > max_occ_drop_ratio:
                    continue
                    
                # Critère de compromis cumulé : le gain cumulé doit valoir le coût cumulé
                if rel_gain_dimers > 0 and max_rel_occ_drop > 0:
                    if rel_gain_dimers < min_tradeoff_factor * max_rel_occ_drop:
                        continue
                        
                # Exiger une baisse cumulée significative des dimères (ex: 20%)
                if dang_initial > 0 and rel_gain_dimers < min_dimer_reduction_ratio:
                    continue
                
                # Seuil d'amélioration absolue et relative du score global par rapport à l'étape courante
                if score < best_score:
                    rel_improvement = (best_score - score) / best_score if best_score > 0 else 0
                    abs_improvement = best_score - score
                    
                    if rel_improvement > min_improvement_ratio and abs_improvement > min_absolute_improvement:
                        best_score = score
                        best_metrics = metrics
                        best_c_total = float(c_test)
                        changed = True
                    
            for idx in par_indices:
                c_current[idx] = best_c_total / n_vars
            
        if not changed:
            stopped_on_max = False
            break
            
    # Résultat final calculé par rapport à l'état initial
    # Final result computed wrt initial state
    results = []
    _, final_metrics = evaluate(c_current)
    
    if stopped_on_max:
        print(f"Optimize concentrations stopped on max_iter={max_iter}")
        
    for par_name, p_list in parent_to_primers.items():
        par_indices = [primer_idx[p.name] for p in p_list]
        orig_c = float(sum(prob_template.total_concentrations[idx] for idx in par_indices))
        new_c = float(sum(c_current[idx] for idx in par_indices))
        
        if abs(orig_c - new_c) > 1e-9:
            results.append({
                "primer_name": par_name,
                "original_conc": orig_c,
                "suggested_conc": new_c,
                "reason": "Réduction des dimères dangereux ou amélioration de l'équilibre multiplexe."
            })
                
    if results:
        dang_before = initial_metrics['dangerous']
        dang_after = final_metrics['dangerous']
        final_rel_gain_dimers = (dang_before - dang_after) / dang_before if dang_before > 0 else 0.0
        
        final_max_rel_occ_drop = 0.0
        for p_key, initial_occ in occ_initial.items():
            final_occ = final_metrics['initiations'].get(p_key, 0.0)
            if initial_occ > 0 and final_occ < initial_occ:
                drop = (initial_occ - final_occ) / initial_occ
                if drop > final_max_rel_occ_drop:
                    final_max_rel_occ_drop = drop
                    
        final_tradeoff_factor = (
            final_rel_gain_dimers / final_max_rel_occ_drop
            if final_max_rel_occ_drop > 0 else float('inf')
        )
        return {
            "suggestions": results,
            "metrics": {
                "dangerous_before": dang_before,
                "dangerous_after": dang_after,
                "dangerous_gain_rel": final_rel_gain_dimers,
                "initiations_before": initial_metrics['initiations'],
                "initiations_after": final_metrics['initiations'],
                "max_occ_drop_rel": final_max_rel_occ_drop,
                "cv_before": initial_metrics['cv'],
                "cv_after": final_metrics['cv'],
                "tradeoff_factor": final_tradeoff_factor
            }
        }
    return []
