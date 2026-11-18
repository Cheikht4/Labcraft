"""ViennaRNA interface with parameter context managers.

Interface ViennaRNA avec gestionnaires de contexte pour les paramètres.
ATTENTION: params_load_* modifie l'état global de la bibliothèque C.
Utiliser dna_params() / rna_params() / dna_session().
"""
import contextlib
import threading
from typing import Optional

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
    """Retourne la profondeur d'imbrication courante du thread / Current thread nesting depth."""
    if not hasattr(_vienna_local, "depth"):
        _vienna_local.depth = 0
    return _vienna_local.depth


def _inc_depth() -> None:
    _vienna_local.depth = _get_depth() + 1


def _dec_depth() -> None:
    _vienna_local.depth = _get_depth() - 1


def _get_mode() -> Optional[str]:
    """Retourne le jeu de paramètres actif ('DNA', 'RNA' ou None) / Active parameter set."""
    return getattr(_vienna_local, "mode", None)


def _set_mode(mode: Optional[str]) -> None:
    _vienna_local.mode = mode


def _ensure_vienna() -> None:
    """Lève une exception si ViennaRNA n'est pas disponible / Raise if ViennaRNA is missing."""
    if not HAS_VIENNA:
        raise ImportError(f"ViennaRNA is required for this module: {_VIENNA_IMPORT_ERROR}")


@contextlib.contextmanager
def dna_params(temp_celsius: float = 65.0, mon_molar: float | None = None):
    """Context manager pour charger les paramètres ADN et l'environnement.
    Context manager loading DNA parameters and reaction environment.

    Verrouille l'état global de la bibliothèque C au premier niveau d'entrée,
    charge Mathews 2004 une seule fois pour toute la portée englobante, ajuste
    température et sel aux niveaux imbriqués sans recharger la table C, puis
    restaure l'état global (Turner 2004 ARN, température, sel) à la sortie finale.
    Locks global C state at outermost entry, loads Mathews 2004 DNA parameters once
    for the enclosing scope, updates temperature/salt on nested calls without
    reloading parameter tables, and restores RNA Turner 2004 state on final exit.
    """
    _ensure_vienna()
    has_salt = hasattr(RNA.cvar, "salt")

    if _get_depth() > 0:
        _inc_depth()
        prev_mode = _get_mode()
        prev_temp = RNA.cvar.temperature
        prev_salt = RNA.cvar.salt if has_salt else None
        temp_changed = prev_temp != temp_celsius
        salt_changed = mon_molar is not None and has_salt and prev_salt != mon_molar
        try:
            if temp_changed:
                RNA.cvar.temperature = temp_celsius
            if salt_changed:
                RNA.cvar.salt = mon_molar
            if prev_mode != "DNA":
                RNA.params_load_DNA_Mathews2004()
                _set_mode("DNA")
            yield
        finally:
            _dec_depth()
            if temp_changed:
                RNA.cvar.temperature = prev_temp
            if salt_changed and prev_salt is not None:
                RNA.cvar.salt = prev_salt
            if prev_mode == "RNA":
                RNA.params_load_RNA_Turner2004()
                _set_mode("RNA")
        return

    with _VIENNA_LOCK:
        _inc_depth()
        saved_temp = RNA.cvar.temperature
        saved_salt = RNA.cvar.salt if has_salt else None
        try:
            RNA.cvar.temperature = temp_celsius
            if mon_molar is not None and has_salt:
                RNA.cvar.salt = mon_molar
            RNA.params_load_DNA_Mathews2004()
            _set_mode("DNA")
            yield
        finally:
            _dec_depth()
            _set_mode(None)
            RNA.cvar.temperature = saved_temp
            if saved_salt is not None and has_salt:
                RNA.cvar.salt = saved_salt
            RNA.params_load_RNA_Turner2004()


@contextlib.contextmanager
def rna_params(temp_celsius: float = 37.0, mon_molar: float | None = None):
    """Context manager pour charger les paramètres ARN par défaut.
    Context manager loading default RNA parameters.
    """
    _ensure_vienna()
    has_salt = hasattr(RNA.cvar, "salt")

    if _get_depth() > 0:
        _inc_depth()
        prev_mode = _get_mode()
        prev_temp = RNA.cvar.temperature
        prev_salt = RNA.cvar.salt if has_salt else None
        temp_changed = prev_temp != temp_celsius
        salt_changed = mon_molar is not None and has_salt and prev_salt != mon_molar
        try:
            if temp_changed:
                RNA.cvar.temperature = temp_celsius
            if salt_changed:
                RNA.cvar.salt = mon_molar
            if prev_mode != "RNA":
                RNA.params_load_RNA_Turner2004()
                _set_mode("RNA")
            yield
        finally:
            _dec_depth()
            if temp_changed:
                RNA.cvar.temperature = prev_temp
            if salt_changed and prev_salt is not None:
                RNA.cvar.salt = prev_salt
            if prev_mode == "DNA":
                RNA.params_load_DNA_Mathews2004()
                _set_mode("DNA")
        return

    with _VIENNA_LOCK:
        _inc_depth()
        saved_temp = RNA.cvar.temperature
        saved_salt = RNA.cvar.salt if has_salt else None
        try:
            RNA.cvar.temperature = temp_celsius
            if mon_molar is not None and has_salt:
                RNA.cvar.salt = mon_molar
            RNA.params_load_RNA_Turner2004()
            _set_mode("RNA")
            yield
        finally:
            _dec_depth()
            _set_mode(None)
            RNA.cvar.temperature = saved_temp
            if saved_salt is not None and has_salt:
                RNA.cvar.salt = saved_salt
            RNA.params_load_RNA_Turner2004()


@contextlib.contextmanager
def dna_session(temp_celsius: float = 65.0, mon_molar: float | None = None):
    """Ouvre une session ADN unique si ViennaRNA est présent (no-op sinon).
    Opens a single-load DNA parameter session when ViennaRNA is available (no-op otherwise).
    """
    if not HAS_VIENNA:
        yield
        return
    with dna_params(temp_celsius=temp_celsius, mon_molar=mon_molar):
        yield
