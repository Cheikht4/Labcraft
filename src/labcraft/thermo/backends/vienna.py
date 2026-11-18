"""ViennaRNA backend (duplexfold/cofold, no GPL) / Backend ViennaRNA"""
from typing import Any
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
            res = DuplexResult(0.0, 0.0, float(mfe), 0.0, structure, temp_celsius)
            _cache_put(self._cache, key, res)
            return res

    def calc_hairpin(self, seq: str, *, temp_celsius: float = 65.0, **kwargs: Any) -> DuplexResult:
        key = ("hairpin", seq, temp_celsius)
        if key in self._cache:
            return self._cache[key]
        with dna_params(temp_celsius):
            structure, mfe = RNA.fold(seq)
            res = DuplexResult(0.0, 0.0, float(mfe), 0.0, structure, temp_celsius)
            _cache_put(self._cache, key, res)
            return res

    def calc_duplex(self, seq1: str, seq2: str, *, temp_celsius: float = 65.0, **kwargs: Any) -> DuplexResult:
        return self.calc_heterodimer(seq1, seq2, temp_celsius=temp_celsius, **kwargs)
