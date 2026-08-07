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
