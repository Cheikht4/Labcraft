"""Abstract base class for duplex energy backends.

Classe abstraite pour les backends de calcul d'énergie de duplexes.
Trois implémentations : NativeBackend (NN pur, sans GPL),
Primer3Backend (référence, GPLv2), ViennaRNABackend (cofold, sans GPL).
Sélectionnable via --duplex-backend.
Aucun module métier n'importe primer3 directement.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass


@dataclass(frozen=True)
class DuplexResult:
    """Result of a duplex energy calculation.
    
    Résultat d'un calcul d'énergie de duplexe.
    
    Attributes:
        dh_kcal: Enthalpy in kcal/mol / Enthalpie en kcal/mol.
        ds_cal_per_k: Entropy in cal/(mol·K) / Entropie en cal/(mol·K).
        dg_kcal: Free energy in kcal/mol at the specified temperature /
            Énergie libre en kcal/mol à la température spécifiée.
        tm_celsius: Melting temperature in °C / Température de fusion en °C.
        structure: Aligned structure string (dot-bracket or ASCII art) /
            Structure alignée (notation point-parenthèse ou ASCII art).
        temperature_celsius: Temperature at which dG was computed /
            Température à laquelle ΔG a été calculé.
    """
    dh_kcal: float
    ds_cal_per_k: float
    dg_kcal: float
    tm_celsius: float
    structure: str
    temperature_celsius: float


class DuplexEnergyBackend(ABC):
    """Abstract backend for computing duplex thermodynamics.
    
    Backend abstrait pour le calcul de la thermodynamique des duplexes.
    """
    
    @abstractmethod
    def calc_heterodimer(
        self, seq1: str, seq2: str, *, temp_celsius: float = 65.0,
        na_mm: float = 50.0, mg_mm: float = 0.0,
    ) -> DuplexResult:
        """Compute heterodimer thermodynamics / Calcul thermodynamique d'un hétérodimère."""
        ...
    
    @abstractmethod
    def calc_homodimer(
        self, seq: str, *, temp_celsius: float = 65.0,
        na_mm: float = 50.0, mg_mm: float = 0.0,
    ) -> DuplexResult:
        """Compute homodimer thermodynamics / Calcul thermodynamique d'un homodimère."""
        ...
    
    @abstractmethod
    def calc_hairpin(
        self, seq: str, *, temp_celsius: float = 65.0,
        na_mm: float = 50.0, mg_mm: float = 0.0,
    ) -> DuplexResult:
        """Compute hairpin thermodynamics / Calcul thermodynamique d'une épingle à cheveux."""
        ...
    
    @abstractmethod
    def calc_duplex(
        self, seq1: str, seq2: str, *, temp_celsius: float = 65.0,
        na_mm: float = 50.0, mg_mm: float = 0.0,
    ) -> DuplexResult:
        """Compute perfect-match duplex thermodynamics.
        
        Calcul thermodynamique d'un duplexe parfaitement apparié.
        """
        ...
