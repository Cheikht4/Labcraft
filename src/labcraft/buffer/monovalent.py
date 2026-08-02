"""Monovalent cation calculation / Calcul des cations monovalents."""

def get_total_monovalent(na_molar: float, k_molar: float, tris_molar: float) -> float:
    """Calcule la concentration totale effective en cations monovalents.
    
    Formule (Owczarzy 2008) : [Mon+] = [Na+] + [K+] + 0.5 * [Tris+]
    Note : Le pH du Tris chute de 8.3 (25°C) à 6.9 (95°C), mais l'article
    indique que cela n'affecte pas la stabilité du duplex dans cette plage.
    
    Args:
        na_molar: Concentration en sodium (Molar).
        k_molar: Concentration en potassium (Molar).
        tris_molar: Concentration en Tris (Molar).
        
    Returns:
        Concentration effective en monovalents (Molar).
    """
    return na_molar + k_molar + 0.5 * tris_molar
