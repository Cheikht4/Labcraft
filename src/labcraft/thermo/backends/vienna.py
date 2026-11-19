"""ViennaRNA backend (duplexfold/cofold, no GPL) / Backend ViennaRNA"""
from typing import Any, Optional
import warnings
import RNA
from labcraft.thermo.vienna import dna_params, _ensure_vienna
from labcraft.thermo.backends.base import DuplexEnergyBackend, DuplexResult

_MAX_CACHE_SIZE = 65536
_SHARED_VIENNA_CACHE: dict[tuple[Any, ...], DuplexResult] = {}


def _cache_put(cache: dict[tuple[Any, ...], DuplexResult], key: tuple[Any, ...], value: DuplexResult) -> None:
    """Insère dans le cache borné partagé / Inserts into the shared bounded cache."""
    if len(cache) >= _MAX_CACHE_SIZE:
        cache.clear()
    cache[key] = value


def _has_intermolecular_pair(structure: str, len_a: int) -> bool:
    """Vérifie qu'au moins une paire de bases traverse la jonction inter-brins `len_a`.
    Checks that at least one base pair spans across the intermolecular junction `len_a`.
    """
    st = structure.replace("&", "")
    stack: list[int] = []
    for idx, ch in enumerate(st):
        if ch == "(":
            stack.append(idx)
        elif ch == ")" and stack:
            open_idx = stack.pop()
            if open_idx < len_a <= idx:
                return True
    return False


def _duplexfold_connected(seq1: str, seq2: str, temp_celsius: float) -> tuple[str, float]:
    """Repliement strictement intermoléculaire (`RNA.duplexfold`) ramené à la longueur totale `len(seq1)+len(seq2)`.
    Strictly intermolecular duplex fold (`RNA.duplexfold`) mapped onto full `len(seq1)+len(seq2)` length.
    """
    d = RNA.duplexfold(seq1, seq2)
    raw_struct = str(d.structure) if d.structure else ""
    if "(" not in raw_struct or "&" not in raw_struct:
        return "." * (len(seq1) + len(seq2)), 0.0
    left_st, right_st = raw_struct.split("&", 1)
    end1 = int(d.i)
    start1 = max(0, end1 - len(left_st))
    start2 = max(0, int(d.j) - 1)
    full_st = (
        "." * start1
        + left_st
        + "." * max(0, len(seq1) - end1)
        + "." * start2
        + right_st
        + "." * max(0, len(seq2) - (start2 + len(right_st)))
    )
    return full_st, float(d.energy)


class ViennaRNABackend(DuplexEnergyBackend):
    """Backend ViennaRNA avec cache MFE partagé entre instances (pur en séquence et température).
    ViennaRNA backend with cross-instance shared MFE cache (pure in sequence and temperature).
    """

    def __init__(self) -> None:
        _ensure_vienna()
        self._cache = _SHARED_VIENNA_CACHE

    def calc_heterodimer(self, seq1: str, seq2: str, *, temp_celsius: float = 65.0, **kwargs: Any) -> DuplexResult:
        key = ("hetero", seq1, seq2, temp_celsius)
        if key in self._cache:
            return self._cache[key]
        with dna_params(temp_celsius):
            seq = f"{seq1}&{seq2}"
            structure, mfe = RNA.cofold(seq)
            if "(" in structure and not _has_intermolecular_pair(structure, len(seq1)):
                # Repliement cofold non connexe (épingles intramoléculaires séparées) : exiger un duplex inter-brins
                # Disconnected cofold (separate intramolecular hairpins): enforce an intermolecular duplex
                structure, mfe = _duplexfold_connected(seq1, seq2, temp_celsius)
            res = DuplexResult(0.0, 0.0, float(mfe), 0.0, structure, temp_celsius)
            _cache_put(self._cache, key, res)
            return res

    def calc_homodimer(self, seq: str, *, temp_celsius: float = 65.0, **kwargs: Any) -> DuplexResult:
        key = ("homo", seq, temp_celsius)
        if key in self._cache:
            return self._cache[key]
        with dna_params(temp_celsius):
            seq2 = f"{seq}&{seq}"
            structure, mfe = RNA.cofold(seq2)
            if "(" in structure and not _has_intermolecular_pair(structure, len(seq)):
                structure, mfe = _duplexfold_connected(seq, seq, temp_celsius)
            res = DuplexResult(0.0, 0.0, float(mfe), 0.0, structure, temp_celsius)
            _cache_put(self._cache, key, res)
            return res

    def calc_hairpin(self, seq: str, *, temp_celsius: float = 65.0, **kwargs: Any) -> DuplexResult:
        """Calcule l'énergie libre MFE d'une épingle intramoléculaire via `RNA.fold`.

        Note : `ViennaRNABackend` simple n'applique pas de paramètres LNA (`lna_positions` est ignoré
        et déclenche un `UserWarning` s'il est non vide ; utiliser `ViennaSaltShiftBackend` pour
        inclure la stabilisation thermodynamique des bases LNA sur les épingles).

        Computes monomolecular hairpin MFE free energy via `RNA.fold`.
        Note: bare `ViennaRNABackend` does not model LNA bases (`lna_positions` is ignored and emits a
        `UserWarning` when non-empty; use `ViennaSaltShiftBackend` for LNA hairpin stabilization).
        """
        lna_positions = kwargs.get("lna_positions")
        if lna_positions:
            warnings.warn(
                "ViennaRNABackend.calc_hairpin ignore lna_positions ; utilisez ViennaSaltShiftBackend "
                "pour appliquer la stabilisation LNA sur les épingles. / "
                "ViennaRNABackend.calc_hairpin ignores lna_positions; use ViennaSaltShiftBackend for LNA.",
                UserWarning,
                stacklevel=2,
            )
        mon_molar: Optional[float] = kwargs.get("mon_molar")
        salt_key = round(float(mon_molar), 6) if mon_molar is not None else None
        key = ("hairpin", seq, temp_celsius, salt_key)
        if key in self._cache:
            return self._cache[key]
        with dna_params(temp_celsius, mon_molar=mon_molar):
            structure, mfe = RNA.fold(seq)
            res = DuplexResult(0.0, 0.0, float(mfe), 0.0, structure, temp_celsius)
            _cache_put(self._cache, key, res)
            return res

    def calc_duplex(self, seq1: str, seq2: str, *, temp_celsius: float = 65.0, **kwargs: Any) -> DuplexResult:
        return self.calc_heterodimer(seq1, seq2, temp_celsius=temp_celsius, **kwargs)
