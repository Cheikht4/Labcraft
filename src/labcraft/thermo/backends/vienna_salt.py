"""ViennaRNA backend with salt shift correction.

Applique un décalage salin (ΔΔG) au ΔG natif de ViennaRNA en estimant
l'enthalpie et l'entropie de la structure via le modèle du plus proche voisin.
"""
from typing import Dict, Tuple, Optional
import math
from .base import DuplexEnergyBackend, DuplexResult
from .vienna import ViennaRNABackend
from .native import _NN_PARAMS, _SEQ_TO_NN_KEY, _INIT_GC, _INIT_AT
from labcraft.thermo.salt import UnifiedSaltModel

def _extract_pairs(structure: str) -> dict[int, int]:
    pairs = {}
    stack = []
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                j = stack.pop()
                pairs[j] = i
                pairs[i] = j
    return pairs

def estimate_helix_thermo(seq: str, structure: str, lna_positions: tuple[int, ...] = ()) -> Tuple[float, float, int, float, float, float]:
    """Estime dH_total, dS_total, n_bp, f_gc, ddH_lna, ddS_lna from sequence and dot-bracket structure."""
    from labcraft.thermo.lna import load_lna_params, _LNA_PARAMS_MXL, _LNA_PARAMS_XLN
    load_lna_params()
    seq = seq.replace('&', '').upper()
    pairs = _extract_pairs(structure)
    
    dh_total = 0.0
    ds_total = 0.0
    n_bp = 0
    gc_count = 0
    ddH_lna = 0.0
    ddS_lna = 0.0
    
    for i in range(len(seq) - 1):
        if i in pairs and (i+1) in pairs:
            j = pairs[i]
            j_minus_1 = pairs[i+1]
            if j_minus_1 == j - 1 and i < j and (i+1) < (j-1):
                dinuc_top = seq[i:i+2]
                try:
                    key = _SEQ_TO_NN_KEY[dinuc_top]
                    dh, ds = _NN_PARAMS[key]
                    dh_total += dh
                    ds_total += ds
                except KeyError:
                    pass
                
                # LNA corrections for this step
                if i in lna_positions:
                    # i is 5' LNA, i+1 is 3' neighbor (TOP STRAND)
                    key_xln = (seq[i], seq[i+1])
                    if key_xln in _LNA_PARAMS_XLN:
                        dh_xln, ds_xln = _LNA_PARAMS_XLN[key_xln]
                        dh_total += dh_xln
                        ds_total += ds_xln
                        ddH_lna += dh_xln
                        ddS_lna += ds_xln
                if (i+1) in lna_positions:
                    # i is 5' neighbor, i+1 is 3' LNA (TOP STRAND)
                    key_mxl = (seq[i], seq[i+1])
                    if key_mxl in _LNA_PARAMS_MXL:
                        dh_mxl, ds_mxl = _LNA_PARAMS_MXL[key_mxl]
                        dh_total += dh_mxl
                        ds_total += ds_mxl
                        ddH_lna += dh_mxl
                        ddS_lna += ds_mxl
                        
                # BOTTOM STRAND
                # top strand is 5'-seq[i] seq[i+1]-3'
                # bottom strand is 3'-seq[j] seq[j-1]-5'
                # So bottom strand 5' to 3' is seq[j-1] to seq[j]
                if (j-1) in lna_positions:
                    # j-1 is 5' LNA, j is 3' neighbor (BOTTOM STRAND)
                    key_xln = (seq[j-1], seq[j])
                    if key_xln in _LNA_PARAMS_XLN:
                        dh_xln, ds_xln = _LNA_PARAMS_XLN[key_xln]
                        dh_total += dh_xln
                        ds_total += ds_xln
                        ddH_lna += dh_xln
                        ddS_lna += ds_xln
                if j in lna_positions:
                    # j-1 is 5' neighbor, j is 3' LNA (BOTTOM STRAND)
                    key_mxl = (seq[j-1], seq[j])
                    if key_mxl in _LNA_PARAMS_MXL:
                        dh_mxl, ds_mxl = _LNA_PARAMS_MXL[key_mxl]
                        dh_total += dh_mxl
                        ds_total += ds_mxl
                        ddH_lna += dh_mxl
                        ddS_lna += ds_mxl
                
    for i in pairs:
        if i < pairs[i]:
            j = pairs[i]
            is_left_end = True
            if (i-1) in pairs and pairs[i-1] == j+1:
                is_left_end = False
                
            is_right_end = True
            if (i+1) in pairs and pairs[i+1] == j-1:
                is_right_end = False
                
            pair_base_1, pair_base_2 = seq[i], seq[j]
            
            if is_left_end:
                if (pair_base_1, pair_base_2) in [('G', 'C'), ('C', 'G')]:
                    dh_total += _INIT_GC[0]
                    ds_total += _INIT_GC[1]
                else:
                    dh_total += _INIT_AT[0]
                    ds_total += _INIT_AT[1]
            if is_right_end:
                if (pair_base_1, pair_base_2) in [('G', 'C'), ('C', 'G')]:
                    dh_total += _INIT_GC[0]
                    ds_total += _INIT_GC[1]
                else:
                    dh_total += _INIT_AT[0]
                    ds_total += _INIT_AT[1]
                    
            n_bp += 1
            if seq[i] in ('G', 'C') or seq[j] in ('G', 'C'):
                gc_count += 1
                
    f_gc = gc_count / n_bp if n_bp > 0 else 0.0
    return dh_total, ds_total, n_bp, f_gc, ddH_lna, ddS_lna

class ViennaSaltShiftBackend(DuplexEnergyBackend):
    def __init__(self, salt_model: UnifiedSaltModel, default_ct_molar: float = 2e-6):
        self.vienna = ViennaRNABackend()
        self.salt_model = salt_model
        self.default_ct_molar = default_ct_molar
        
    def _apply_shift(self, res: DuplexResult, seq: str, ct_molar: float, kwargs: dict, lna_positions: tuple[int, ...] = ()) -> DuplexResult:
        if not res.structure or '(' not in res.structure:
            return res
            
        dh_helix, ds_helix, n_bp, f_gc, ddH_lna, ddS_lna = estimate_helix_thermo(seq, res.structure, lna_positions)
        
        na_molar = kwargs.get('na_mm', 50.0) / 1000.0
        k_molar = kwargs.get('k_mm', 0.0) / 1000.0
        tris_molar = kwargs.get('tris_mm', 0.0) / 1000.0
        mg_molar = kwargs.get('mg_mm', 0.0) / 1000.0
        dntp_molar = kwargs.get('dntp_mm', 0.0) / 1000.0
        
        dh_corr, ds_corr = self.salt_model.corrected_thermodynamics(
            dh_helix, ds_helix, f_gc, n_bp, 
            na_molar, k_molar, tris_molar, mg_molar, dntp_molar
        )
        
        temp_k = res.temperature_celsius + 273.15
        
        # Calculate ddG_LNA
        ddg_lna = ddH_lna - temp_k * (ddS_lna / 1000.0)
        
        # Salt applies on top of the already LNA-corrected dh_helix / ds_helix
        ddg_salt = -temp_k * (ds_corr - ds_helix) / 1000.0
        
        dg_final = res.dg_kcal + ddg_salt + ddg_lna
        
        # Recalculate Tm coherently.
        r_gas = 1.9872
        
        x = 4.0
        if '&' in seq:
            p1, p2 = seq.split('&', 1)
            if p1 == p2:
                x = 1.0
        else:
            x = 1.0
            
        if res.tm_celsius is not None and not math.isnan(res.tm_celsius):
            if ct_molar > 0:
                new_tm_k = (dh_corr * 1000.0) / (ds_corr + r_gas * math.log(ct_molar / x))
                new_tm_celsius = new_tm_k - 273.15
            else:
                new_tm_celsius = res.tm_celsius
        else:
            new_tm_celsius = float('nan')
            
        return DuplexResult(
            dh_kcal=round(dh_corr, 2),
            ds_cal_per_k=round(ds_corr, 2),
            dg_kcal=dg_final,
            tm_celsius=new_tm_celsius,
            structure=res.structure,
            temperature_celsius=res.temperature_celsius
        )

    def calc_heterodimer(self, seq1: str, seq2: str, *, temp_celsius: float = 65.0, **kwargs) -> DuplexResult:
        res = self.vienna.calc_heterodimer(seq1, seq2, temp_celsius=temp_celsius, **kwargs)
        ct = kwargs.get('ct_molar', self.default_ct_molar)
        if ct is None: ct = self.default_ct_molar
        lna_pos_a = kwargs.get('lna_positions_a', ())
        lna_pos_b = kwargs.get('lna_positions_b', ())
        # Remap lna_positions on concatenated sequence 'seq1seq2' (without &)
        mapped_lna = list(lna_pos_a) + [pos + len(seq1) for pos in lna_pos_b]
        return self._apply_shift(res, f"{seq1}&{seq2}", ct, kwargs, tuple(mapped_lna))
        
    def calc_homodimer(self, seq: str, *, temp_celsius: float = 65.0, **kwargs) -> DuplexResult:
        res = self.vienna.calc_homodimer(seq, temp_celsius=temp_celsius, **kwargs)
        ct = kwargs.get('ct_molar', self.default_ct_molar)
        if ct is None: ct = self.default_ct_molar
        lna_pos = kwargs.get('lna_positions', ())
        mapped_lna = list(lna_pos) + [pos + len(seq) for pos in lna_pos]
        return self._apply_shift(res, f"{seq}&{seq}", ct, kwargs, tuple(mapped_lna))
        
    def calc_hairpin(self, seq: str, *, temp_celsius: float = 65.0, **kwargs) -> DuplexResult:
        res = self.vienna.calc_hairpin(seq, temp_celsius=temp_celsius, **kwargs)
        ct = kwargs.get('ct_molar', self.default_ct_molar)
        if ct is None: ct = self.default_ct_molar
        lna_pos = kwargs.get('lna_positions', ())
        return self._apply_shift(res, seq, ct, kwargs, tuple(lna_pos))
        
    def calc_duplex(self, seq1: str, seq2: str, *, temp_celsius: float = 65.0, **kwargs) -> DuplexResult:
        return self.calc_heterodimer(seq1, seq2, temp_celsius=temp_celsius, **kwargs)
