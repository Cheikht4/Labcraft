from typing import List, Dict, Any, Tuple
import copy
import dataclasses
import numpy as np

from labcraft.solver.types import EquilibriumProblem
from labcraft.solver.dual import solve_dual
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.metrics.balance import calculate_initiation_balance_cv
from labcraft.diagnostics.amplifiable_dimer import evaluate_pair_amplifiable

def get_dangerous_dimers(primers: List[PhysicalPrimer], temp_celsius: float, backend: Any, enzyme: Any) -> set[str]:
    dangerous = set()
    for p1 in primers:
        for p2 in primers:
            if evaluate_pair_amplifiable(p1, p2, backend, enzyme, temp_celsius)[0]:
                dangerous.add(f"{p1.name}_on_{p2.name}")
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
    max_iter: int = 3
) -> List[Dict[str, Any]]:
    
    # Identifier les dimères dangereux
    dangerous_complexes = get_dangerous_dimers(primers, temp_celsius, backend, enzyme)
    
    # Indexer
    primer_idx = {p.name: species_names.index(p.name) for p in primers if p.name in species_names}
    target_indices = {s for s in species_names if s.endswith('_site')}
    
    # Trouver les complexes d'initiation
    initiation_complexes = {}
    for i, c in enumerate(complex_names):
        for p in primers:
            if f"{p.name}_on_{p.name}_site" == c:
                initiation_complexes[p.name] = i
                
    dangerous_indices = [i for i, c in enumerate(complex_names) if c in dangerous_complexes]
    
    def evaluate(c_candidats: np.ndarray) -> Tuple[float, dict]:
        prob = dataclasses.replace(prob_template, total_concentrations=c_candidats)
        try:
            res = solve_dual(prob)
        except Exception:
            return float('inf'), {}
            
        score = 0.0
        
        # 1. Dimères dangereux
        dang_frac = np.sum(res.concentrations[dangerous_indices])
        score += dang_frac * 1e9 # Forte pénalité
        
        # 2. Plancher d'amorçage
        initiations = {}
        for p_name, c_idx in initiation_complexes.items():
            occ = res.concentrations[c_idx] / c_candidats[primer_idx[p_name]] if c_candidats[primer_idx[p_name]] > 0 else 0
            initiations[p_name] = occ
            if occ < 0.001:
                score += 1e6 # Plancher
                
        # 3. CV Balance
        if initiations:
            cv = np.std(list(initiations.values())) / (np.mean(list(initiations.values())) + 1e-9)
            score += cv * 100
            
        return score, {'dangerous': dang_frac, 'initiations': initiations, 'cv': cv if initiations else 0}
        
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
    
    best_score, _ = evaluate(c_current)
    
    for _ in range(max_iter):
        changed = False
        for p in primers:
            if p.name not in primer_idx or p.role not in ranges:
                continue
                
            idx = primer_idx[p.name]
            best_c = c_current[idx]
            
            for c_test in ranges[p.role]:
                c_candidats = np.copy(c_current)
                c_candidats[idx] = c_test
                score, metrics = evaluate(c_candidats)
                
                if score < best_score:
                    best_score = score
                    best_c = c_test
                    changed = True
                    
            c_current[idx] = best_c
            
        if not changed:
            break
            
    # Résultat
    results = []
    final_score, final_metrics = evaluate(c_current)
    
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
                    "reason": "Réduction des dimères dangereux ou amélioration de l'équilibre."
                })
                
    return results
