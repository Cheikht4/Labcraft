"""Target accessibility (ΔG_unfolding) / Accessibilité de la cible.

Calcul de la pénalité d'ouverture de la cible via la fonction de partition
de ViennaRNA (ensemble free energy).
"""
import RNA

from labcraft.thermo.vienna import dna_params
from labcraft.target.cache import unfolding_cache


def calc_unfolding_penalty(
    target_seq: str, bind_start: int, bind_end: int,
    temp_celsius: float = 65.0, mon_molar: float | None = None
) -> float:
    """Calcule le coût d'ouverture ΔG_unfolding d'un site sur la cible.
    
    La séquence doit être orientée 5'->3'. bind_start et bind_end sont 0-based.
    La cible sera analysée dans son intégralité telle que fournie (fenêtrage
    doit être fait par l'appelant).
    
    Args:
        target_seq: Séquence 5'->3' de la fenêtre cible.
        bind_start: Index (0-based) de début du site de liaison.
        bind_end: Index (0-based) de fin du site de liaison (exclusif).
        temp_celsius: Température de repliement.
        mon_molar: Concentration totale en cations monovalents.
        
    Returns:
        Pénalité d'ouverture en kcal/mol (toujours >= 0).
    """
    if bind_end <= bind_start:
        raise ValueError("bind_end doit être strictement supérieur à bind_start.")
    if bind_start < 0 or bind_end > len(target_seq):
        raise ValueError(f"Site ({bind_start}, {bind_end}) hors limites de la séquence ({len(target_seq)}).")
        
    # Clé de cache : inclut le mon_molar si présent
    cache_key = (target_seq, bind_start, bind_end, temp_celsius, mon_molar)
    cached_val = unfolding_cache.get(cache_key)
    if cached_val is not None:
        return cached_val

    with dna_params(temp_celsius=temp_celsius, mon_molar=mon_molar):
        # 1. Ensemble libre
        fc_free = RNA.fold_compound(target_seq)
        # cvar.pf_scale est parfois requis pour stabiliser le partitionnement sur de longues séquences
        # fc_free.pf() rescale automatiquement si on l'appelle tel quel.
        _, g_free = fc_free.pf()
        
        # 2. Ensemble contraint (le site est forcé d'être non apparié)
        fc_constrained = RNA.fold_compound(target_seq)
        
        # Dans ViennaRNA, les indices sont 1-based.
        # Donc base 0-based `i` correspond à la position `i+1`.
        # hc_add_up(idx) force la base à la position 1-based `idx` à être non appariée.
        for i in range(bind_start, bind_end):
            fc_constrained.hc_add_up(i + 1)
            
        _, g_constrained = fc_constrained.pf()
        
    # ΔG_unfolding = G_constrained - G_free
    dg_unfolding = g_constrained - g_free
    
    # Plancher à 0 avec une légère tolérance numérique (bruit de float)
    if dg_unfolding < 0.0 and dg_unfolding > -1e-2:
        dg_unfolding = 0.0
    elif dg_unfolding < 0.0:
        # Théoriquement impossible en physique exacte, G_constrained >= G_free
        # puisque l'ensemble contraint est un sous-ensemble de l'ensemble libre.
        import logging
        logging.warning(f"Negative unfolding penalty calculated: {dg_unfolding:.2f} kcal/mol. This is theoretically impossible.")
        dg_unfolding = 0.0
        
    unfolding_cache.put(cache_key, dg_unfolding)
    return dg_unfolding
