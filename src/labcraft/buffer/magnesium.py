"""Free magnesium calculation / Calcul du magnésium libre."""

def get_free_magnesium(mg_total_molar: float, dntp_total_molar: float) -> float:
    """Calcule la concentration de magnésium libre après chélation par les dNTPs.
    
    Niveau 1: soustraction directe (Owczarzy 2008).
    
    Args:
        mg_total_molar: Concentration totale en magnésium (Molar).
        dntp_total_molar: Concentration totale en dNTPs (Molar).
        
    Returns:
        Concentration en magnésium libre (Molar), >= 0.0.
    """
    return max(0.0, mg_total_molar - dntp_total_molar)
