"""Folding energy cache / Cache des énergies de repliement."""

class UnfoldingCache:
    """Cache simple pour éviter de recalculer pf() sur les mêmes fenêtres.
    
    Les calculs O(N^3) de ViennaRNA sont très coûteux, on met en cache
    la valeur de ΔG_unfolding.
    """
    
    def __init__(self, max_size: int = 10000):
        self.cache: dict[tuple, float] = {}
        self.max_size = max_size
        
    def get(self, key: tuple) -> float | None:
        return self.cache.get(key)
        
    def put(self, key: tuple, value: float) -> None:
        if len(self.cache) >= self.max_size:
            # Nettoyage arbitraire brutal si débordement (pour ne pas plomber la RAM)
            self.cache.clear()
        self.cache[key] = value
        
    def clear(self) -> None:
        self.cache.clear()


# Instance globale pour le module target
unfolding_cache = UnfoldingCache()
