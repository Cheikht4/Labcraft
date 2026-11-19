from dataclasses import dataclass
from typing import Any, Dict, Union

@dataclass
class PolymeraseProfile:
    """Profil enzymatique pour le diagnostic des dimères et du veto 3'.
    Enzymatic profile for dimer diagnostics and 3'-end extension veto.
    
    Attributes:
        name: Nom de la polymérase (ex: Bst 2.0) / Polymerase name.
        strand_displacement: Vrai si forte activité de déplacement de brin / Strand-displacement activity.
        warm_start: Vrai si inhibée à température ambiante / Inhibited at ambient temperature.
        dimer_dg_threshold: Seuil heuristique interne de ΔG_3' local (kcal/mol, évalué à la température
            et au sel de réaction sans pénalité d'initiation sur <= 6 pb contiguës en 3') en dessous
            duquel un dimère 3'-apparié avec matrice 5' est classé amplifiable.
            Internal heuristic local ΔG_3' threshold (kcal/mol, evaluated at reaction temperature and
            salt without initiation penalty over <= 6 contiguous 3' bp) below which a 3'-paired dimer
            with 5' template overhang is classified as amplifiable.
    """
    name: str
    strand_displacement: bool
    warm_start: bool
    dimer_dg_threshold: float
    # Modèle du Veto 3' selon ARMS (Newton et al. 1989, NAR 17(7):2503-2516)
    # 3' Veto model following ARMS (Newton et al. 1989, NAR 17(7):2503-2516)
    # - Amorces 7b/7c (pos 3) : effet fort / strong effect
    # - Amorces 8b/8c (pos 4) : effet marginal / marginal effect
    # - Amorces 7f/7g (pos 7) : effet fuyant / leaky effect
    three_prime_window: int = 3  # Fenêtre de forte pénalité (jusqu'à la position 3) / Strong-penalty window
    three_prime_absolute_window: int = 2  # Fenêtre de blocage absolu de l'extension (pos 1 et 2) / Absolute block window
    # Option cinétique (Kwok et al. 1990) : tolérer un wobble G:T / T:G en position pénultième (-2)
    # Kinetic option (Kwok et al. 1990): tolerate a G:T / T:G wobble at penultimate position (-2)
    allow_penultimate_gt_wobble: bool = False


# NOTE SUR L'ORIGINE DES SEUILS `dimer_dg_threshold` (HEURISTIQUE NON CALIBRÉE) :
# Aucune publication biochimique ne fournit de seuil universel de ΔG_3' propre à Bst, Bst 2.0,
# Bst 3.0, GspSSD ou Taq à 63–65 °C. Dans LabCraft, `local_dg` est évalué à la température réelle
# de réaction et au sel du tampon sans pénalité d'initiation sur <= 6 pb contiguës en 3'
# (`amplifiable_dimer.py`) : entre (37 °C, 1 M Na+) et (63 °C, tampon LAMP standard 50 mM Na+,
# 10 mM K+, 20 mM Tris, 6 mM Mg2+, 1,4 mM dNTP), ΔG_3' remonte de +1,50 kcal/mol sur 3 pb
# (`5'-GGA-3'`, -3,14 -> -1,64), +2,37 sur 4 pb (`5'-GCAT-3'`, -4,56 -> -2,20), +3,16 sur 5 pb
# (`5'-GCATG-3'`, -6,02 -> -2,86) et +4,00 kcal/mol sur 6 pb (`5'-GCATGC-3'`, -8,26 -> -4,26).
# Les valeurs ci-dessous (-2.0 pour Bst / Bst 2.0 / GspSSD, -1.5 pour Bst 3.0, -3.0 pour Taq)
# constituent une heuristique non calibrée ; les données de paillasse actuelles ne permettent pas
# de discriminer entre -2.0 et -3.5 kcal/mol sur les 5 cas de référence.
#
# NOTE ON `dimer_dg_threshold` PROVENANCE (UNCALIBRATED HEURISTIC):
# No peer-reviewed biochemical study publishes a universal polymerase-specific ΔG_3' cutoff for
# Bst, Bst 2.0, Bst 3.0, GspSSD, or Taq at 63–65 °C. In LabCraft, `local_dg` is evaluated at the
# actual reaction temperature and buffer salt without initiation penalty over <= 6 contiguous 3' bp
# (`amplifiable_dimer.py`): going from (37 °C, 1 M Na+) to (63 °C, standard LAMP buffer 50 mM Na+,
# 10 mM K+, 20 mM Tris, 6 mM Mg2+, 1.4 mM dNTP) shifts ΔG_3' upward by +1.50 kcal/mol on 3 bp
# (`5'-GGA-3'`, -3.14 -> -1.64), +2.37 on 4 bp (`5'-GCAT-3'`, -4.56 -> -2.20), +3.16 on 5 bp
# (`5'-GCATG-3'`, -6.02 -> -2.86), and +4.00 kcal/mol on 6 bp (`5'-GCATGC-3'`, -8.26 -> -4.26).
# The values below (-2.0 for Bst / Bst 2.0 / GspSSD, -1.5 for Bst 3.0, -3.0 for Taq) are an
# uncalibrated heuristic; current wet-lab data do not discriminate between -2.0 and -3.5 kcal/mol
# across the 5 reference cases.

# Profils standards de polymérases / Standard polymerase profiles
BST = PolymeraseProfile(
    name="Bst",
    strand_displacement=True,
    warm_start=False,
    dimer_dg_threshold=-2.0,  # Heuristique interne non publiée / Unpublished internal heuristic
)

BST_2_0 = PolymeraseProfile(
    name="Bst 2.0",
    strand_displacement=True,
    warm_start=False,
    dimer_dg_threshold=-2.0,  # Heuristique interne non publiée / Unpublished internal heuristic
)

BST_2_0_WARM_START = PolymeraseProfile(
    name="Bst 2.0 WarmStart",
    strand_displacement=True,
    warm_start=True,
    # Même seuil heuristique (-2.0) : l'aptamère WarmStart n'agit qu'à froid, pas à 63–65 °C
    # Same heuristic threshold (-2.0): WarmStart aptamer only inhibits below ~45 °C, not at 63–65 °C
    dimer_dg_threshold=-2.0,
)

TAQ = PolymeraseProfile(
    name="Taq",
    strand_displacement=False,
    warm_start=False,
    # Heuristique interne non publiée (absence de déplacement de brin isotherme)
    # Unpublished internal heuristic (no isothermal strand displacement)
    dimer_dg_threshold=-3.0,
)

# Bst 3.0 : activité transcriptase inverse et déplacement de brin renforcés ; des amplifications
# non spécifiques accrues sont rapportées (p. ex. Lopez-Jimena et al. 2018), mais la valeur -1.5
# kcal/mol est une heuristique interne provisoire non issue d'une mesure enzymatique publiée.
# Bst 3.0: stronger reverse-transcriptase and strand-displacement activity; higher non-specific
# background is reported (e.g. Lopez-Jimena et al. 2018), but -1.5 kcal/mol is a provisional
# internal heuristic not derived from a published enzymatic measurement.
BST_3_0 = PolymeraseProfile(
    name="Bst 3.0",
    strand_displacement=True,
    warm_start=False,
    dimer_dg_threshold=-1.5,
)

# GspSSD (polymérase des kits Optigene, ex: ISO-004) :
# Non calibré expérimentalement ; aucune valeur propre de seuil 3' ou de fenêtre ARMS n'est publiée dans la littérature.
# Reprend à l'identique les paramètres de Bst 2.0 sans inventer de seuil propre.
# GspSSD (polymerase used in Optigene kits, e.g. ISO-004):
# Not experimentally calibrated; no published specific 3' threshold or ARMS window values exist in the literature.
# Reuses Bst 2.0 parameters identically without inventing a custom threshold.
GSPSSD = PolymeraseProfile(
    name="GspSSD",
    strand_displacement=True,
    warm_start=False,
    dimer_dg_threshold=-2.0
)

ENZYME_REGISTRY = {
    "bst": BST,
    "bst2.0": BST_2_0,
    "bst 2.0": BST_2_0,
    "bst2.0ws": BST_2_0_WARM_START,
    "bst2.0warmstart": BST_2_0_WARM_START,
    "bst2.0_ws": BST_2_0_WARM_START,
    "bst2.0_warmstart": BST_2_0_WARM_START,
    "warmstart": BST_2_0_WARM_START,
    "bst3.0": BST_3_0,
    "bst 3.0": BST_3_0,
    "gspssd": GSPSSD,
    "gspssd2.0": GSPSSD,
    "iso004": GSPSSD,
    "taq": TAQ
}

def normalize_enzyme_name(name: str) -> str:
    # Retire les espaces, tirets et majuscules
    return name.lower().replace(" ", "").replace("-", "")

def get_enzyme(spec: Union[str, Dict[str, Any]]) -> PolymeraseProfile:
    """
    Récupère un profil de polymérase depuis une chaîne ou un dictionnaire.
    Permet de surcharger certains paramètres (dimer_dg_threshold, three_prime_window).
    """
    if isinstance(spec, str):
        name_key = normalize_enzyme_name(spec)
        # Fallback pour certains alias avec underscores (si on n'a pas tout remplacé)
        if name_key not in ENZYME_REGISTRY and spec.lower() in ENZYME_REGISTRY:
            name_key = spec.lower()
            
        if name_key not in ENZYME_REGISTRY:
            raise ValueError(f"Enzyme inconnue: {spec}. Disponibles: {list(ENZYME_REGISTRY.keys())}")
        return ENZYME_REGISTRY[name_key]
        
    elif isinstance(spec, dict):
        if "name" not in spec:
            raise ValueError("La spécification de l'enzyme doit contenir une clé 'name'.")
        
        name_key = normalize_enzyme_name(spec["name"])
        if name_key not in ENZYME_REGISTRY and spec["name"].lower() in ENZYME_REGISTRY:
            name_key = spec["name"].lower()
            
        if name_key not in ENZYME_REGISTRY:
            raise ValueError(f"Enzyme inconnue: {spec['name']}. Disponibles: {list(ENZYME_REGISTRY.keys())}")
            
        base_profile = ENZYME_REGISTRY[name_key]
        
        # Surcharges optionnelles / Optional overrides
        dg_thresh = spec.get("dimer_dg_threshold", base_profile.dimer_dg_threshold)
        tp_window = spec.get("three_prime_window", base_profile.three_prime_window)
        tp_abs_window = spec.get("three_prime_absolute_window", base_profile.three_prime_absolute_window)
        allow_wobble = spec.get("allow_penultimate_gt_wobble", base_profile.allow_penultimate_gt_wobble)
        
        return PolymeraseProfile(
            name=base_profile.name,
            strand_displacement=base_profile.strand_displacement,
            warm_start=base_profile.warm_start,
            dimer_dg_threshold=float(dg_thresh),
            three_prime_window=int(tp_window),
            three_prime_absolute_window=int(tp_abs_window),
            allow_penultimate_gt_wobble=bool(allow_wobble),
        )
        
    else:
        raise TypeError("La spécification de l'enzyme doit être une chaîne ou un dictionnaire.")
