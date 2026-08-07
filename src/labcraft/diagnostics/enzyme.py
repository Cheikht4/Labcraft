from dataclasses import dataclass

@dataclass
class PolymeraseProfile:
    """Profil enzymatique pour le diagnostic des dimères.
    
    Attributes:
        name: Nom de la polymérase (ex: Bst 2.0).
        strand_displacement: Vrai si l'enzyme a une forte activité de déplacement de brin.
        warm_start: Vrai si l'enzyme est inhibée à température ambiante.
        dimer_dg_threshold: Seuil de ΔG pour qu'un dimère 3'-apparié soit classé amplifiable (kcal/mol).
    """
    name: str
    strand_displacement: bool
    warm_start: bool
    dimer_dg_threshold: float
    # Modèle du Veto 3' selon ARMS (Newton et al. 1989, NAR 17(7):2503-2516)
    # - Amorces 7b/7c (pos 3) : effet fort
    # - Amorces 8b/8c (pos 4) : effet marginal
    # - Amorces 7f/7g (pos 7) : effet fuyant
    three_prime_window: int = 3  # Fenêtre de forte pénalité (jusqu'à la position 3)
    three_prime_absolute_window: int = 2  # Fenêtre de blocage absolu de l'extension (pos 1 et 2)

# Profils standards de polymérases
BST = PolymeraseProfile(
    name="Bst",
    strand_displacement=True,
    warm_start=False,
    dimer_dg_threshold=-2.0
)

BST_2_0 = PolymeraseProfile(
    name="Bst 2.0",
    strand_displacement=True,
    warm_start=False,
    dimer_dg_threshold=-2.0
)

BST_2_0_WARM_START = PolymeraseProfile(
    name="Bst 2.0 WarmStart",
    strand_displacement=True,
    warm_start=True,
    dimer_dg_threshold=-2.0 # Le seuil reste le même, la physique à 65°C ne change pas
)

TAQ = PolymeraseProfile(
    name="Taq",
    strand_displacement=False,
    warm_start=False,
    dimer_dg_threshold=-3.0 # Moins agressive sur l'extension de structures complexes
)

# La Bst 3.0 est plus agressive, avec une activité transcriptase inverse et déplacement de brin
# renforcée. Elle tolère des appariements 3' plus faibles.
# Le seuil de -1.5 kcal/mol est une valeur PROVISOIRE à calibrer expérimentalement.
# Des amplifications non spécifiques avec Bst 3.0 sont rapportées dans la littérature 
# (e.g. Lopez-Jimena 2018), nécessitant un seuil moins strict que Bst 2.0.
BST_3_0 = PolymeraseProfile(
    name="Bst 3.0",
    strand_displacement=True,
    warm_start=False,
    dimer_dg_threshold=-1.5
)

ENZYME_REGISTRY = {
    "bst": BST,
    "bst2.0": BST_2_0,
    "bst2.0_ws": BST_2_0_WARM_START,
    "bst3.0": BST_3_0,
    "taq": TAQ
}

def get_enzyme(spec) -> PolymeraseProfile:
    """
    Récupère un profil de polymérase depuis une chaîne ou un dictionnaire.
    Permet de surcharger certains paramètres (dimer_dg_threshold, three_prime_window).
    """
    if isinstance(spec, str):
        name_key = spec.lower()
        if name_key not in ENZYME_REGISTRY:
            raise ValueError(f"Enzyme inconnue: {spec}. Disponibles: {list(ENZYME_REGISTRY.keys())}")
        return ENZYME_REGISTRY[name_key]
        
    elif isinstance(spec, dict):
        if "name" not in spec:
            raise ValueError("La spécification de l'enzyme doit contenir une clé 'name'.")
        
        name_key = spec["name"].lower()
        if name_key not in ENZYME_REGISTRY:
            raise ValueError(f"Enzyme inconnue: {spec['name']}. Disponibles: {list(ENZYME_REGISTRY.keys())}")
            
        base_profile = ENZYME_REGISTRY[name_key]
        
        # Surcharges optionnelles
        dg_thresh = spec.get("dimer_dg_threshold", base_profile.dimer_dg_threshold)
        tp_window = spec.get("three_prime_window", base_profile.three_prime_window)
        tp_abs_window = spec.get("three_prime_absolute_window", base_profile.three_prime_absolute_window)
        
        return PolymeraseProfile(
            name=base_profile.name,
            strand_displacement=base_profile.strand_displacement,
            warm_start=base_profile.warm_start,
            dimer_dg_threshold=float(dg_thresh),
            three_prime_window=int(tp_window),
            three_prime_absolute_window=int(tp_abs_window)
        )
        
    else:
        raise TypeError("La spécification de l'enzyme doit être une chaîne ou un dictionnaire.")
