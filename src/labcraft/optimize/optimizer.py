"""Main optimizer logic.

Logique principale de l'optimiseur.
"""
from typing import List, Dict, Any, Optional, Tuple

from labcraft.optimize.mutations import enumerate_variants
from labcraft.optimize.constraints import filter_no_new_artefacts, calc_intrinsic_3p_dg
from labcraft.optimize.ranking import score_variants
from labcraft.diagnostics.enzyme import BST, PolymeraseProfile
from labcraft.thermo.backends.base import DuplexEnergyBackend
from labcraft.thermo.vienna import dna_session


def optimize_primer(
    primer_name: str,
    primers: List[Dict[str, Any]],
    temp_celsius: float = 63.0,
    max_mutations: int = 2,
    window_3p: int = 6,
    enzyme: Optional[PolymeraseProfile] = None,
    backend: Optional[DuplexEnergyBackend] = None,
    **backend_kwargs: Any,
) -> List[Dict[str, Any]]:
    """
    Optimise une amorce ciblée dans un panel pour supprimer les dimères amplifiables.
    Optimizes a targeted primer in a panel to eliminate amplifiable dimers.
    
    Args:
        primer_name: Nom de l'amorce à optimiser / Name of the primer to optimize.
        primers: Liste des amorces du panel (dict avec 'name', 'sequence', etc.) / Panel primers list.
        temp_celsius: Température de la réaction / Reaction temperature in Celsius.
        max_mutations: Nombre max de mutations (1 ou 2) / Maximum substitutions (1 or 2).
        window_3p: Taille de la fenêtre 3' à explorer / 3' window size to explore.
        enzyme: Profil de polymérase configuré (défaut : BST) / Configured polymerase profile (default: BST).
        backend: Backend thermodynamique optionnel (ex. ViennaSaltShiftBackend) / Optional salt-aware backend.
        **backend_kwargs: Paramètres ioniques optionnels (na_mm, mg_mm, dntp_mm, etc.) / Optional ionic kwargs.
        
    Returns:
        Une liste des meilleurs variants (top 3) avec leurs métadonnées / Top 3 ranked variants.
    """
    enzyme_profile = enzyme if enzyme is not None else BST

    # Trouver l'amorce / Locate target primer
    modified_idx = -1
    for i, p in enumerate(primers):
        if p['name'] == primer_name:
            modified_idx = i
            break
            
    if modified_idx == -1:
        raise ValueError(f"Amorce '{primer_name}' introuvable dans le panel.")
        
    original_seq = str(primers[modified_idx]['sequence'])
    
    # 1. Énumération des variants / Enumerate candidate variants
    variants = enumerate_variants(original_seq, max_mutations, window_3p)
    
    # Calcul de la stabilité intrinsèque de l'original / Compute original 3' intrinsic stability
    orig_intrinsic_dg = calc_intrinsic_3p_dg(original_seq, window_3p, temp_celsius)
    
    # 2. Application des filtres durs (évaluation thermodynamique) / Apply strict thermodynamic filters
    valid_variants: List[Tuple[str, List[Tuple[int, str, str]], float, float]] = []
    
    # L'évaluation s'exécute dans une session ADN unique (Mathews 2004)
    # Evaluation runs inside a single DNA parameter session (Mathews 2004)
    with dna_session(temp_celsius=temp_celsius):
        for seq, muts in variants:
            # Filtre : rejeter tout variant qui rend l'extrémité 3' plus collante que l'original
            # Filter: reject any variant making the 3' end stickier than the original
            intrinsic_dg = calc_intrinsic_3p_dg(seq, window_3p, temp_celsius)
            if intrinsic_dg < orig_intrinsic_dg:
                continue
                
            # Filtre strict : aucun dimère amplifiable impliquant le variant (hétéro et homo)
            # Strict filter: zero amplifiable dimers involving the variant (hetero and homo)
            is_valid, worst_dg_3p = filter_no_new_artefacts(
                seq,
                primers,
                modified_idx,
                enzyme_profile,
                temp_celsius,
                backend=backend,
                **backend_kwargs,
            )
            
            if is_valid:
                valid_variants.append((seq, muts, worst_dg_3p, intrinsic_dg))
                
    # 3. Classement des survivants / Rank surviving variants
    ranked = score_variants(valid_variants)
    
    # Retourne les 3 meilleurs / Return top 3
    return ranked[:3]
