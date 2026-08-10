"""Native NN backend (pure Python, no GPL).

Backend NN natif (Python pur, sans GPL).
Implémente le modèle Nearest-Neighbor de SantaLucia 1998 pour les duplexes
parfaitement appariés, avec correction de symétrie.
"""
from __future__ import annotations
from labcraft.thermo.constants import R_GAS_KCAL_MOL_K

import math
from typing import Dict, Tuple

from .base import DuplexEnergyBackend, DuplexResult

# Constantes SantaLucia 1998 (Table 2 PNAS 95:1460)
# Format : (dH en kcal/mol, dS en cal/mol·K)
_NN_PARAMS: Dict[str, Tuple[float, float]] = {
    'AA/TT': (-7.9, -22.2),
    'AT/TA': (-7.2, -20.4),
    'TA/AT': (-7.2, -21.3),
    'CA/GT': (-8.5, -22.7),
    'GT/CA': (-8.4, -22.4),
    'CT/GA': (-7.8, -21.0),
    'GA/CT': (-8.2, -22.2),
    'CG/GC': (-10.6, -27.2),
    'GC/CG': (-9.8, -24.4),
    'GG/CC': (-8.0, -19.9),
}

_INIT_GC = (0.1, -2.8)
_INIT_AT = (2.3, 4.1)
_SYM = (0.0, -1.4)

_SEQ_TO_NN_KEY = {
    'AA': 'AA/TT', 'TT': 'AA/TT',
    'AT': 'AT/TA',
    'TA': 'TA/AT',
    'CA': 'CA/GT', 'TG': 'CA/GT',
    'GT': 'GT/CA', 'AC': 'GT/CA',
    'CT': 'CT/GA', 'AG': 'CT/GA',
    'GA': 'GA/CT', 'TC': 'GA/CT',
    'CG': 'CG/GC',
    'GC': 'GC/CG',
    'GG': 'GG/CC', 'CC': 'GG/CC',
}


def _is_self_complementary(seq: str) -> bool:
    """Détermine si la séquence est auto-complémentaire."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    try:
        rev_comp = "".join(complement[b] for b in reversed(seq))
        return seq == rev_comp
    except KeyError:
        return False


class NativeBackend(DuplexEnergyBackend):
    """Native backend implementing SantaLucia 1998 thermodynamics.
    
    Ne traite que les duplexes parfaitement appariés.
    """
    
    def __init__(self, default_ct_molar: float = 2e-6):
        """Initialise le backend avec une concentration totale par défaut.
        
        Args:
            default_ct_molar: Concentration totale des brins C_T en M.
                              Vaut 2 µM par défaut (pour la validation Owczarzy).
        """
        self.default_ct_molar = default_ct_molar

    def _calc_perfect_duplex(
        self, seq: str, temp_celsius: float, ct_molar: float, saltcorr: float = 0.0,
        lna_positions: tuple[int, ...] = ()
    ) -> DuplexResult:
        """Helper for perfect duplexes (homodimer or matching heterodimer)"""
        from labcraft.thermo.lna import load_lna_params, _LNA_PARAMS_MXL, _LNA_PARAMS_XLN
        load_lna_params()
        
        seq = seq.upper()
        if not seq or len(seq) < 2:
            raise ValueError("Sequence too short for NN computation.")
            
        dh_total = 0.0
        ds_total = 0.0
        
        # Initiation penalties
        for end in (seq[0], seq[-1]):
            if end in ('G', 'C'):
                dh_total += _INIT_GC[0]
                ds_total += _INIT_GC[1]
            elif end in ('A', 'T'):
                dh_total += _INIT_AT[0]
                ds_total += _INIT_AT[1]
            else:
                raise ValueError(f"Invalid base '{end}' in sequence.")
                
        # Symmetry correction
        self_comp = _is_self_complementary(seq)
        if self_comp:
            dh_total += _SYM[0]
            ds_total += _SYM[1]
            x = 1.0
        else:
            x = 4.0
            
        # Nearest-Neighbor terms
        for i in range(len(seq) - 1):
            dinuc = seq[i:i+2]
            try:
                key = _SEQ_TO_NN_KEY[dinuc]
                dh, ds = _NN_PARAMS[key]
                dh_total += dh
                ds_total += ds
            except KeyError:
                import warnings
                warnings.warn(f"Ignored unknown dinucleotide '{dinuc}' in sequence.")
        # McTigue 2004 LNA Corrections
        for idx in lna_positions:
            if idx > 0:
                base5 = seq[idx-1]
                base_lna = seq[idx]
                key_mxl = (base5, base_lna)
                if key_mxl in _LNA_PARAMS_MXL:
                    dh_mxl, ds_mxl = _LNA_PARAMS_MXL[key_mxl]
                    dh_total += dh_mxl
                    ds_total += ds_mxl
            if idx < len(seq) - 1:
                base_lna = seq[idx]
                base3 = seq[idx+1]
                key_xln = (base_lna, base3)
                if key_xln in _LNA_PARAMS_XLN:
                    dh_xln, ds_xln = _LNA_PARAMS_XLN[key_xln]
                    dh_total += dh_xln
                    ds_total += ds_xln
                    
        # Salt correction (if provided by an external wrapper, else 0.0)
        ds_total += saltcorr
        
        # Free Energy
        temp_k = temp_celsius + 273.15
        dg_kcal = dh_total - temp_k * (ds_total / 1000.0)
        
        # Melting Temperature
        r_gas = R_GAS_KCAL_MOL_K * 1000.0
        # Tm convention SantaLucia : Tm = dH / (dS + R*ln(Ct/x))
        if ct_molar > 0:
            tm_k = (dh_total * 1000.0) / (ds_total + r_gas * math.log(ct_molar / x))
            tm_celsius = tm_k - 273.15
        else:
            tm_celsius = float('nan')
            
        # Generate dot-bracket structure for perfect duplex
        structure = "(" * len(seq) + "+" + ")" * len(seq)
        
        return DuplexResult(
            dh_kcal=round(dh_total, 2),
            ds_cal_per_k=round(ds_total, 2),
            dg_kcal=dg_kcal,
            tm_celsius=tm_celsius,
            structure=structure,
            temperature_celsius=temp_celsius
        )

    def calc_heterodimer(
        self, seq1: str, seq2: str, *, temp_celsius: float = 65.0,
        na_mm: float = 50.0, mg_mm: float = 0.0,
        k_mm: float = 0.0, tris_mm: float = 0.0, dntp_mm: float = 0.0,
        ct_molar: float | None = None,
        **kwargs
    ) -> DuplexResult:
        """Compute heterodimer thermodynamics (must be perfect match)."""
        lna_positions = kwargs.get('lna_positions', ())
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        try:
            comp_seq = "".join(complement[c] for c in reversed(seq1.upper()))
        except KeyError:
            raise ValueError("Invalid character in sequence for complementation.")
            
        if comp_seq != seq2.upper():
            raise NotImplementedError("NativeBackend only supports perfectly matched heterodimers.")
            
        ct = ct_molar if ct_molar is not None else self.default_ct_molar
        return self._calc_perfect_duplex(seq1, temp_celsius, ct, saltcorr=0.0, lna_positions=lna_positions)
        
    def calc_homodimer(
        self, seq: str, *, temp_celsius: float = 65.0,
        na_mm: float = 50.0, mg_mm: float = 0.0,
        k_mm: float = 0.0, tris_mm: float = 0.0, dntp_mm: float = 0.0,
        ct_molar: float | None = None,
        lna_positions: tuple[int, ...] = ()
    ) -> DuplexResult:
        """Compute homodimer thermodynamics (must be self-complementary)."""
        if not _is_self_complementary(seq):
            raise ValueError("Sequence is not self-complementary, cannot form perfect homodimer.")
        ct = ct_molar if ct_molar is not None else self.default_ct_molar
        return self._calc_perfect_duplex(seq, temp_celsius, ct, saltcorr=0.0, lna_positions=lna_positions)

    def calc_hairpin(
        self, seq: str, *, temp_celsius: float = 65.0,
        na_mm: float = 50.0, mg_mm: float = 0.0,
        k_mm: float = 0.0, tris_mm: float = 0.0, dntp_mm: float = 0.0,
        ct_molar: float | None = None,
    ) -> DuplexResult:
        raise NotImplementedError("NativeBackend does not support hairpins.")

    def calc_duplex(
        self, seq1: str, seq2: str, *, temp_celsius: float = 65.0,
        na_mm: float = 50.0, mg_mm: float = 0.0,
        k_mm: float = 0.0, tris_mm: float = 0.0, dntp_mm: float = 0.0,
        ct_molar: float | None = None,
        **kwargs
    ) -> DuplexResult:
        return self.calc_heterodimer(
            seq1, seq2, temp_celsius=temp_celsius, 
            na_mm=na_mm, mg_mm=mg_mm, 
            k_mm=k_mm, tris_mm=tris_mm, dntp_mm=dntp_mm, 
            ct_molar=ct_molar, **kwargs
        )
