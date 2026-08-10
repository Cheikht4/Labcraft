"""Ionic corrections (Owczarzy 2004, 2008) / Corrections ioniques.

Modèles orthogonaux de correction saline pour la prédiction thermodynamique.
Ces modèles sont injectés dans le pipeline via SaltCorrectedBackend.
"""
from __future__ import annotations

import math
import warnings
from abc import ABC, abstractmethod

from .backends.base import DuplexEnergyBackend, DuplexResult
from labcraft.buffer.magnesium import get_free_magnesium
from labcraft.buffer.monovalent import get_total_monovalent


class SaltModel(ABC):
    """Interface abstraite pour les modèles de correction saline."""

    @abstractmethod
    def correct_tm(
        self, tm_ref_celsius: float, f_gc: float, n_bp: int,
        na_molar: float = 0.0, k_molar: float = 0.0, tris_molar: float = 0.0,
        mg_molar: float = 0.0, dntp_molar: float = 0.0,
        na_ref_molar: float = 1.0
    ) -> float:
        """Chemin empirique : calcule le Tm corrigé.
        
        Args:
            tm_ref_celsius: Tm de référence en degrés Celsius à la concentration `na_ref_molar`.
            f_gc: Fraction de paires G-C dans le duplexe (0.0 à 1.0).
            n_bp: Nombre de paires de bases.
            na_molar, k_molar, tris_molar, mg_molar, dntp_molar: Concentrations du tampon.
            na_ref_molar: Concentration de référence en cations (Molar), défaut 1.0 M.
            
        Returns:
            Le Tm corrigé en degrés Celsius.
        """
        ...

    @abstractmethod
    def corrected_thermodynamics(
        self, dh_kcal: float, ds_cal: float, f_gc: float, n_bp: int,
        na_molar: float = 0.0, k_molar: float = 0.0, tris_molar: float = 0.0,
        mg_molar: float = 0.0, dntp_molar: float = 0.0,
    ) -> tuple[float, float]:
        """Chemin entropique (Solveur) : calcule les dH et dS corrigés.
        
        Repose sur l'identité ΔS°_corrige = ΔS°(ref) + ΔH° * Δ(1/Tm).
        """
        ...


class UnifiedSaltModel(SaltModel):
    """Modèle unifié d'Owczarzy (2004, 2008) pour Na+, Mg2+, et régimes mixtes.
    
    Ref: 
    - Owczarzy et al. 2004, Biochemistry 43:3537 (Na+)
    - Owczarzy et al. 2008, Biochemistry 47:5336 (Mg2+ et Mixte)
    """

    def __init__(self, use_legacy_entropy: bool = False):
        self.use_legacy_entropy = use_legacy_entropy

    def _calc_delta_inv_tm(
        self, f_gc: float, n_bp: int, mon_total: float, mg_free: float, na_ref_molar: float
    ) -> float:
        """Calcule le terme de correction empirique Δ(1/Tm) en K^-1."""
        if n_bp <= 1:
            return 0.0
        if mon_total == 0.0 and mg_free == 0.0:
            return 0.0  # Aucune correction si l'utilisateur donne 0 absolu partout

        # Eq 16 params (Table 2)
        a = 3.92e-5
        b = -9.11e-6
        c = 6.26e-5
        d = 1.42e-5
        e = -4.82e-4
        f = 5.25e-4
        g = 8.31e-5

        if mon_total == 0.0:
            # Magnésium pur
            if not (0.5e-3 <= mg_free <= 600e-3):
                warnings.warn("Mg2+ concentration outside calibrated range (0.5 - 600 mM).", RuntimeWarning)
            
            ln_mg = math.log(mg_free)
            delta_inv_tm = (
                a + b * ln_mg
                + f_gc * (c + d * ln_mg)
                + (1.0 / (2.0 * (n_bp - 1))) * (e + f * ln_mg + g * (ln_mg ** 2))
            )
            return delta_inv_tm

        r_ratio = math.sqrt(mg_free) / mon_total

        if r_ratio < 0.22:
            # Monovalent dominant (Owczarzy 2004, Eq 22)
            if not (0.05 <= mon_total <= 1.1) or not (0.05 <= na_ref_molar <= 1.1):
                warnings.warn("Monovalent concentration outside calibrated range (0.05 - 1.1 M).", RuntimeWarning)
                
            ln_mon = math.log(mon_total)
            ln_ref = math.log(na_ref_molar)
            term1 = (4.29 * f_gc - 3.95) * 1e-5 * (ln_mon - ln_ref)
            term2 = 9.40 * 1e-6 * (ln_mon**2 - ln_ref**2)
            return term1 + term2
            
        elif r_ratio < 6.0:
            # Régime Mixte (Owczarzy 2008, Eqs 18-20)
            ln_mon = math.log(mon_total)
            sqrt_mon = math.sqrt(mon_total)
            
            a = 3.92e-5 * (0.843 - 0.352 * sqrt_mon * ln_mon)
            d = 1.42e-5 * (1.279 - 4.03e-3 * ln_mon - 8.03e-3 * (ln_mon ** 2))
            g = 8.31e-5 * (0.486 - 0.258 * ln_mon + 5.25e-3 * (ln_mon ** 3))
            
            ln_mg = math.log(mg_free)
            delta_inv_tm = (
                a + b * ln_mg
                + f_gc * (c + d * ln_mg)
                + (1.0 / (2.0 * (n_bp - 1))) * (e + f * ln_mg + g * (ln_mg ** 2))
            )
            return delta_inv_tm
            
        else:
            # Magnésium dominant
            if not (0.5e-3 <= mg_free <= 600e-3):
                warnings.warn("Mg2+ concentration outside calibrated range (0.5 - 600 mM).", RuntimeWarning)
            
            ln_mg = math.log(mg_free)
            delta_inv_tm = (
                a + b * ln_mg
                + f_gc * (c + d * ln_mg)
                + (1.0 / (2.0 * (n_bp - 1))) * (e + f * ln_mg + g * (ln_mg ** 2))
            )
            return delta_inv_tm

    def correct_tm(
        self, tm_ref_celsius: float, f_gc: float, n_bp: int,
        na_molar: float = 0.0, k_molar: float = 0.0, tris_molar: float = 0.0,
        mg_molar: float = 0.0, dntp_molar: float = 0.0,
        na_ref_molar: float = 1.0
    ) -> float:
        mon_total = get_total_monovalent(na_molar, k_molar, tris_molar)
        mg_free = get_free_magnesium(mg_molar, dntp_molar)
        
        delta_inv = self._calc_delta_inv_tm(f_gc, n_bp, mon_total, mg_free, na_ref_molar)
        
        tm_ref_k = tm_ref_celsius + 273.15
        inv_tm_target = (1.0 / tm_ref_k) + delta_inv
        
        return (1.0 / inv_tm_target) - 273.15

    def corrected_thermodynamics(
        self, dh_kcal: float, ds_cal: float, f_gc: float, n_bp: int,
        na_molar: float = 0.0, k_molar: float = 0.0, tris_molar: float = 0.0,
        mg_molar: float = 0.0, dntp_molar: float = 0.0,
    ) -> tuple[float, float]:
        mon_total = get_total_monovalent(na_molar, k_molar, tris_molar)
        mg_free = get_free_magnesium(mg_molar, dntp_molar)
        
        if self.use_legacy_entropy:
            # Ancien modèle grossier (ne gère que Na+)
            if mon_total > 0:
                ds_corr = ds_cal + 0.368 * (n_bp - 1) * math.log(mon_total)
            else:
                ds_corr = ds_cal
            return dh_kcal, ds_corr
            
        delta_inv = self._calc_delta_inv_tm(f_gc, n_bp, mon_total, mg_free, 1.0)
        
        # Identité: dS_corr = dS(1M) + dH(1M) * delta(1/Tm)
        # Note: dH est en kcal/mol, delta(1/Tm) est en K^-1. On doit multiplier dh par 1000 pour avoir des cal.
        # ds_cal est en cal/mol.K.
        dh_cal = dh_kcal * 1000.0
        ds_corr = ds_cal + dh_cal * delta_inv
        
        return dh_kcal, ds_corr


# Rétro-compatibilité pour les tests existants et nom explicite
Owczarzy2004SaltModel = UnifiedSaltModel


class SaltCorrectedBackend(DuplexEnergyBackend):
    """Décorateur (Wrapper) injectant une correction saline sur un backend."""

    def __init__(self, backend: DuplexEnergyBackend, salt_model: SaltModel):
        self._backend = backend
        self._salt_model = salt_model

    def _apply_correction(
        self, result: DuplexResult, temp_celsius: float, seq: str, 
        na_mm: float, mg_mm: float, k_mm: float, tris_mm: float, dntp_mm: float, 
        ct_molar: float
    ) -> DuplexResult:
        """Applique la correction unifiée entropique et recalcule dg et tm."""
        # On suppose que le backend de base calcule par défaut à 1 M Na+ (pas de Mg, pas de K)
        # Si on est déjà aux conditions de réf, on pourrait shunter, mais on laisse le UnifiedSaltModel gérer 
        # (si mon=1.0 et mg=0, delta=0).
        
        seq = seq.upper()
        gc_count = seq.count('G') + seq.count('C')
        f_gc = gc_count / len(seq) if seq else 0.0
        
        dh_kcal, ds_cal = self._salt_model.corrected_thermodynamics(
            dh_kcal=result.dh_kcal,
            ds_cal=result.ds_cal_per_k,
            f_gc=f_gc,
            n_bp=len(seq),
            na_molar=na_mm / 1000.0,
            k_molar=k_mm / 1000.0,
            tris_molar=tris_mm / 1000.0,
            mg_molar=mg_mm / 1000.0,
            dntp_molar=dntp_mm / 1000.0
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
        k_mm: float = 0.0, tris_mm: float = 0.0, dntp_mm: float = 0.0,
        ct_molar: float | None = None,
        **kwargs
    ) -> DuplexResult:
        res_1m = self._backend.calc_heterodimer(
            seq1, seq2, temp_celsius=temp_celsius, 
            na_mm=1000.0, mg_mm=0.0, k_mm=0.0, tris_mm=0.0, dntp_mm=0.0,
            ct_molar=ct_molar, **kwargs
        )
        eff_ct = ct_molar if ct_molar is not None else getattr(self._backend, 'default_ct_molar', 2e-6)
        return self._apply_correction(
            res_1m, temp_celsius, seq1, na_mm, mg_mm, k_mm, tris_mm, dntp_mm, eff_ct
        )

    def calc_homodimer(
        self, seq: str, *, temp_celsius: float = 65.0,
        na_mm: float = 50.0, mg_mm: float = 0.0,
        k_mm: float = 0.0, tris_mm: float = 0.0, dntp_mm: float = 0.0,
        ct_molar: float | None = None,
        **kwargs
    ) -> DuplexResult:
        res_1m = self._backend.calc_homodimer(
            seq, temp_celsius=temp_celsius,
            na_mm=1000.0, mg_mm=0.0, k_mm=0.0, tris_mm=0.0, dntp_mm=0.0,
            ct_molar=ct_molar, **kwargs
        )
        eff_ct = ct_molar if ct_molar is not None else getattr(self._backend, 'default_ct_molar', 2e-6)
        return self._apply_correction(
            res_1m, temp_celsius, seq, na_mm, mg_mm, k_mm, tris_mm, dntp_mm, eff_ct
        )

    def calc_hairpin(
        self, seq: str, *, temp_celsius: float = 65.0,
        na_mm: float = 50.0, mg_mm: float = 0.0,
        k_mm: float = 0.0, tris_mm: float = 0.0, dntp_mm: float = 0.0,
        ct_molar: float | None = None,
        **kwargs
    ) -> DuplexResult:
        res_1m = self._backend.calc_hairpin(
            seq, temp_celsius=temp_celsius,
            na_mm=1000.0, mg_mm=0.0, k_mm=0.0, tris_mm=0.0, dntp_mm=0.0,
            ct_molar=ct_molar, **kwargs
        )
        eff_ct = ct_molar if ct_molar is not None else getattr(self._backend, 'default_ct_molar', 2e-6)
        return self._apply_correction(
            res_1m, temp_celsius, seq, na_mm, mg_mm, k_mm, tris_mm, dntp_mm, eff_ct
        )

    def calc_duplex(
        self, seq1: str, seq2: str, *, temp_celsius: float = 65.0,
        na_mm: float = 50.0, mg_mm: float = 0.0,
        k_mm: float = 0.0, tris_mm: float = 0.0, dntp_mm: float = 0.0,
        ct_molar: float | None = None,
        **kwargs
    ) -> DuplexResult:
        return self.calc_heterodimer(
            seq1, seq2, temp_celsius=temp_celsius,
            na_mm=na_mm, mg_mm=mg_mm, k_mm=k_mm, tris_mm=tris_mm, dntp_mm=dntp_mm,
            ct_molar=ct_molar, **kwargs
        )

def sodium_equivalent_for_folding(
    mon_molar: float, mg_molar: float, dntp_molar: float = 0.0,
    f_gc: float = 0.5, n_bp: int = 20, cap_molar: float = 1.0
) -> float:
    """Monovalent effectif (mol/L) reproduisant la correction saline Owczarzy
    Na+Mg du tampon réel, pour un repliement ViennaRNA qui ne connaît que le sel
    monovalent. Dérivé du modèle interne, pas de von Ahsen."""
    mg_free = get_free_magnesium(mg_molar, dntp_molar)
    if mg_free <= 0:
        return mon_molar
    
    if mon_molar <= 0:
        mon_molar = 1e-3
        
    model = UnifiedSaltModel()
    
    # Target delta inv Tm using actual buffer
    tgt = model._calc_delta_inv_tm(f_gc, n_bp, mon_molar, mg_free, na_ref_molar=1.0)
    
    lo, hi = mon_molar, cap_molar
    for _ in range(60):
        mid = (lo * hi) ** 0.5
        d = model._calc_delta_inv_tm(f_gc, n_bp, mid, 0.0, 1.0)
        if d > tgt:   # d décroît quand mid augmente
            lo = mid
        else:
            hi = mid
            
    return (lo * hi) ** 0.5
