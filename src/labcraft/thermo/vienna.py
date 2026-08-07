"""ViennaRNA interface with parameter context managers.

Interface ViennaRNA avec gestionnaires de contexte pour les paramètres.
ATTENTION: params_load_* modifie l'état global de la bibliothèque C.
Utiliser dna_params() / rna_params().
"""
import contextlib
import threading

try:
    import RNA
    if not hasattr(RNA, "params_load_DNA_Mathews2004"):
        raise ImportError("ViennaRNA version is too old. RNA.params_load_DNA_Mathews2004 not found.")
    HAS_VIENNA = True
except ImportError as err:
    HAS_VIENNA = False
    _VIENNA_IMPORT_ERROR = err

_VIENNA_LOCK = threading.Lock()
_vienna_local = threading.local()

def _get_depth() -> int:
    if not hasattr(_vienna_local, 'depth'):
        _vienna_local.depth = 0
    return _vienna_local.depth

def _inc_depth():
    _vienna_local.depth = _get_depth() + 1

def _dec_depth():
    _vienna_local.depth = _get_depth() - 1



def _ensure_vienna():
    """Lève une exception si ViennaRNA n'est pas disponible."""
    if not HAS_VIENNA:
        raise ImportError(f"ViennaRNA is required for this module: {_VIENNA_IMPORT_ERROR}")


@contextlib.contextmanager
def dna_params(temp_celsius: float = 65.0, mon_molar: float | None = None):
    """Context manager pour charger les paramètres ADN et l'environnement.
    
    Verrouille l'état global de la bibliothèque C pour éviter toute course,
    charge les paramètres ADN, applique la température de réaction, et
    (si supporté) la concentration en sel monovalent, puis
    restaure l'état global à la sortie.
    
    Note: Le magnésium n'est pas modélisé ici, car la fonction pf()
    de ViennaRNA ne prend en compte (officiellement) que les sels monovalents.
    
    Args:
        temp_celsius: Température de repliement en degrés Celsius.
        mon_molar: Concentration totale en cations monovalents (Na+ + K+ + Tris/2).
    """
    _ensure_vienna()
    
    if _get_depth() > 0:
        _inc_depth()
        try:
            yield
        finally:
            _dec_depth()
        return
        
    with _VIENNA_LOCK:
        _inc_depth()
        # Snapshot
        saved_temp = RNA.cvar.temperature
        saved_salt = RNA.cvar.salt if hasattr(RNA.cvar, 'salt') else None
        
        try:
            RNA.cvar.temperature = temp_celsius
            if mon_molar is not None and hasattr(RNA.cvar, 'salt'):
                RNA.cvar.salt = mon_molar
                
            RNA.params_load_DNA_Mathews2004()
            yield
        finally:
            _dec_depth()
            # Restauration
            RNA.cvar.temperature = saved_temp
            if saved_salt is not None:
                RNA.cvar.salt = saved_salt
            RNA.params_load_RNA_Turner2004()


@contextlib.contextmanager
def rna_params(temp_celsius: float = 37.0, mon_molar: float | None = None):
    """Context manager pour charger les paramètres ARN par défaut."""
    _ensure_vienna()
    
    if _get_depth() > 0:
        _inc_depth()
        try:
            yield
        finally:
            _dec_depth()
        return
        
    with _VIENNA_LOCK:
        _inc_depth()
        saved_temp = RNA.cvar.temperature
        saved_salt = RNA.cvar.salt if hasattr(RNA.cvar, 'salt') else None
        
        try:
            RNA.cvar.temperature = temp_celsius
            if mon_molar is not None and hasattr(RNA.cvar, 'salt'):
                RNA.cvar.salt = mon_molar
                
            RNA.params_load_RNA_Turner2004()
            yield
        finally:
            _dec_depth()
            RNA.cvar.temperature = saved_temp
            if saved_salt is not None:
                RNA.cvar.salt = saved_salt
            RNA.params_load_RNA_Turner2004()
