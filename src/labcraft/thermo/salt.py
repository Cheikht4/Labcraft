"""Ionic corrections (Owczarzy 2004, 2008) / Corrections ioniques.

Modèles orthogonaux de correction saline pour la prédiction thermodynamique.
Ces modèles sont injectés dans le pipeline via SaltCorrectedBackend.
"""
from __future__ import annotations

import math
from abc import ABC, abstractmethod
import warnings

from .backends.base import DuplexEnergyBackend, DuplexResult


class SaltModel(ABC):
    """Interface abstraite pour les modèles de correction saline."""

    @abstractmethod
    def correct_tm(
        self, tm_ref_celsius: float, f_gc: float, na_target_molar: float,
        na_ref_molar: float = 1.0
    ) -> float:
        """Chemin empirique : calcule le Tm corrigé.
        
        Args:
            tm_ref_celsius: Tm de référence en degrés Celsius à la concentration `na_ref_molar`.
            f_gc: Fraction de paires G-C dans le duplexe (0.0 à 1.0).
            na_target_molar: Concentration cible en cations (Molar).
            na_ref_molar: Concentration de référence en cations (Molar), défaut 1.0 M.
            
        Returns:
            Le Tm corrigé en degrés Celsius.
        """
        ...

    @abstractmethod
    def corrected_thermodynamics(
        self, dh_kcal: float, ds_cal: float, na_target_molar: float, n_bp: int
    ) -> tuple[float, float]:
        """Chemin entropique (Solveur) : calcule les dH et dS corrigés.
        
        Args:
            dh_kcal: Enthalpie de référence à 1 M (kcal/mol).
            ds_cal: Entropie de référence à 1 M (cal/mol·K).
            na_target_molar: Concentration cible en cations (Molar).
            n_bp: Longueur du duplexe en paires de bases.
            
        Returns:
            Tuple (dH_corrige_kcal, dS_corrige_cal).
        """
        ...


class Owczarzy2004SaltModel(SaltModel):
    """Modèle Owczarzy 2004 pour les cations monovalents (Na+).
    
    Ref: Biochemistry 2004, 43, 3537-3554.
    """

    def correct_tm(
        self, tm_ref_celsius: float, f_gc: float, na_target_molar: float,
        na_ref_molar: float = 1.0
    ) -> float:
        if na_target_molar <= 0 or na_ref_molar <= 0:
            raise ValueError("Sodium concentrations must be > 0.")
            
        if not (0.05 <= na_target_molar <= 1.1) or not (0.05 <= na_ref_molar <= 1.1):
            warnings.warn("Owczarzy 2004 eq 22 is calibrated for [Na+] between 0.05 M and 1.1 M.", RuntimeWarning)

        tm_ref_k = tm_ref_celsius + 273.15
        inv_tm_ref = 1.0 / tm_ref_k
        
        ln_na2 = math.log(na_target_molar)
        ln_na1 = math.log(na_ref_molar)
        
        term1 = (4.29 * f_gc - 3.95) * 1e-5 * (ln_na2 - ln_na1)
        term2 = 9.40 * 1e-6 * (ln_na2**2 - ln_na1**2)
        
        inv_tm_target = inv_tm_ref + term1 + term2
        return (1.0 / inv_tm_target) - 273.15

    def corrected_thermodynamics(
        self, dh_kcal: float, ds_cal: float, na_target_molar: float, n_bp: int
    ) -> tuple[float, float]:
        if na_target_molar <= 0:
            raise ValueError("Sodium concentration must be > 0.")
            
        # Delta H remains unchanged
        # Delta S corrected: dS(Na) = dS(1M) + 0.368 * (N - 1) * ln[Na+]
        ds_corr = ds_cal + 0.368 * (n_bp - 1) * math.log(na_target_molar)
        return dh_kcal, ds_corr


class SaltCorrectedBackend(DuplexEnergyBackend):
    """Décorateur (Wrapper) injectant une correction saline sur un backend."""

    def __init__(self, backend: DuplexEnergyBackend, salt_model: SaltModel):
        self._backend = backend
        self._salt_model = salt_model

    def _apply_correction(self, result: DuplexResult, temp_celsius: float, na_mm: float, seq: str, ct_molar: float) -> DuplexResult:
        """Applique la correction entropique et recalcule dg et tm."""
        if na_mm == 1000.0:
            return result # Pas de correction nécessaire si on est déjà à 1 M (la réf)

        # On suppose que le backend de base calcule par défaut à 1 M
        # Le salt model recalcule dH, dS
        dh_kcal, ds_cal = self._salt_model.corrected_thermodynamics(
            dh_kcal=result.dh_kcal,
            ds_cal=result.ds_cal_per_k,
            na_target_molar=na_mm / 1000.0,
            n_bp=len(seq)
        )
        
        temp_k = temp_celsius + 273.15
        dg_kcal = dh_kcal - temp_k * (ds_cal / 1000.0)
        
        # Recalcul du Tm avec la même convention de concentration
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        try:
            rev_comp = "".join(complement[b] for b in reversed(seq))
            is_sym = (seq == rev_comp)
        except KeyError:
            is_sym = False
            
        x = 1.0 if is_sym else 4.0
        r_gas = 1.9872
        
        if ct_molar > 0:
            tm_k = (dh_kcal * 1000.0) / (ds_cal + r_gas * math.log(ct_molar / x))
            tm_celsius = tm_k - 273.15
        else:
            tm_celsius = float('nan')

        return DuplexResult(
            dh_kcal=dh_kcal,
            ds_cal_per_k=ds_cal,
            dg_kcal=dg_kcal,
            tm_celsius=tm_celsius,
            structure=result.structure,
            temperature_celsius=temp_celsius
        )

    def calc_heterodimer(
        self, seq1: str, seq2: str, *, temp_celsius: float = 65.0,
        na_mm: float = 50.0, mg_mm: float = 0.0,
        ct_molar: float | None = None
    ) -> DuplexResult:
        # Appel du backend interne à 1000 mM (1 M Na+) pour obtenir la référence
        res_1m = self._backend.calc_heterodimer(
            seq1, seq2, temp_celsius=temp_celsius, 
            na_mm=1000.0, mg_mm=mg_mm, ct_molar=ct_molar
        )
        # On utilise le ct_molar par défaut du backend si non spécifié (NativeBackend fournit default_ct_molar)
        eff_ct = ct_molar if ct_molar is not None else getattr(self._backend, 'default_ct_molar', 2e-6)
        return self._apply_correction(res_1m, temp_celsius, na_mm, seq1, eff_ct)

    def calc_homodimer(
        self, seq: str, *, temp_celsius: float = 65.0,
        na_mm: float = 50.0, mg_mm: float = 0.0,
        ct_molar: float | None = None
    ) -> DuplexResult:
        res_1m = self._backend.calc_homodimer(
            seq, temp_celsius=temp_celsius,
            na_mm=1000.0, mg_mm=mg_mm, ct_molar=ct_molar
        )
        eff_ct = ct_molar if ct_molar is not None else getattr(self._backend, 'default_ct_molar', 2e-6)
        return self._apply_correction(res_1m, temp_celsius, na_mm, seq, eff_ct)

    def calc_hairpin(
        self, seq: str, *, temp_celsius: float = 65.0,
        na_mm: float = 50.0, mg_mm: float = 0.0,
        ct_molar: float | None = None
    ) -> DuplexResult:
        res_1m = self._backend.calc_hairpin(
            seq, temp_celsius=temp_celsius,
            na_mm=1000.0, mg_mm=mg_mm, ct_molar=ct_molar
        )
        eff_ct = ct_molar if ct_molar is not None else getattr(self._backend, 'default_ct_molar', 2e-6)
        return self._apply_correction(res_1m, temp_celsius, na_mm, seq, eff_ct)

    def calc_duplex(
        self, seq1: str, seq2: str, *, temp_celsius: float = 65.0,
        na_mm: float = 50.0, mg_mm: float = 0.0,
        ct_molar: float | None = None
    ) -> DuplexResult:
        return self.calc_heterodimer(
            seq1, seq2, temp_celsius=temp_celsius,
            na_mm=na_mm, mg_mm=mg_mm, ct_molar=ct_molar
        )
