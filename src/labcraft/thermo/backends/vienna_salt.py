from labcraft.thermo.constants import R_GAS_KCAL_MOL_K
"""ViennaRNA backend with salt shift correction.

Applique un décalage salin (ΔΔG) au ΔG natif de ViennaRNA en estimant
l'enthalpie et l'entropie de la structure via le modèle du plus proche voisin.
"""
from typing import Any, Dict, Tuple
import math
from .base import DuplexEnergyBackend, DuplexResult
from .vienna import ViennaRNABackend
from .native import _NN_PARAMS, _SEQ_TO_NN_KEY, _INIT_GC, _INIT_AT
from labcraft.thermo.salt import UnifiedSaltModel

def _extract_pairs(structure: str) -> dict[int, int]:
    """Extrait les paires d'indices (0-based) d'une structure dot-bracket sans '&'.
    Extracts 0-based base-pair index mapping from a dot-bracket structure without '&'.
    """
    pairs: dict[int, int] = {}
    stack: list[int] = []
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                j = stack.pop()
                pairs[j] = i
                pairs[i] = j
    return pairs


def _estimate_structure_thermo_impl(
    seq: str,
    structure: str,
    lna_positions: tuple[int, ...] = (),
    *,
    include_loops: bool = False,
) -> Tuple[float, float, int, float, float, float, bool]:
    """Estime (dH, dS, n_bp, f_gc, ddH_lna, ddS_lna, has_loops) sur une structure secondaire.
    Estimates (dH, dS, n_bp, f_gc, ddH_lna, ddS_lna, has_loops) on a secondary structure.
    """
    from labcraft.thermo.lna import load_lna_params, _LNA_PARAMS_MXL, _LNA_PARAMS_XLN
    load_lna_params()

    cut_idx = seq.find('&')
    seq_clean = seq.replace('&', '').upper()
    struct_clean = structure.replace('&', '')
    pairs = _extract_pairs(struct_clean)
    
    dh_total = 0.0
    ds_total = 0.0
    n_bp = 0
    gc_count = 0
    ddH_lna = 0.0
    ddS_lna = 0.0
    
    for i in range(len(seq_clean) - 1):
        if include_loops and cut_idx != -1 and (i + 1) == cut_idx:
            # Ne pas empiler à travers la coupure inter-brins '&' en mode structure complète
            # Do not stack across the intermolecular '&' strand break in full-structure mode
            continue
        if i in pairs and (i + 1) in pairs:
            j = pairs[i]
            j_minus_1 = pairs[i + 1]
            if include_loops and cut_idx != -1 and j == cut_idx:
                continue
            if j_minus_1 == j - 1 and i < j and (i + 1) < (j - 1):
                dinuc_top = seq_clean[i:i + 2]
                try:
                    key = _SEQ_TO_NN_KEY[dinuc_top]
                    dh, ds = _NN_PARAMS[key]
                    dh_total += dh
                    ds_total += ds
                except KeyError:
                    import warnings
                    warnings.warn(f"Ignored unknown dinucleotide '{dinuc_top}' in sequence.", stacklevel=2)
                
                # Corrections LNA sur le brin supérieur (5'->3') / Top-strand LNA corrections
                if i in lna_positions:
                    key_xln = (seq_clean[i], seq_clean[i + 1])
                    if key_xln in _LNA_PARAMS_XLN:
                        dh_xln, ds_xln = _LNA_PARAMS_XLN[key_xln]
                        dh_total += dh_xln
                        ds_total += ds_xln
                        ddH_lna += dh_xln
                        ddS_lna += ds_xln
                if (i + 1) in lna_positions:
                    key_mxl = (seq_clean[i], seq_clean[i + 1])
                    if key_mxl in _LNA_PARAMS_MXL:
                        dh_mxl, ds_mxl = _LNA_PARAMS_MXL[key_mxl]
                        dh_total += dh_mxl
                        ds_total += ds_mxl
                        ddH_lna += dh_mxl
                        ddS_lna += ds_mxl
                        
                # Corrections LNA sur le brin inférieur (5'->3' : j-1 vers j) / Bottom-strand LNA corrections
                if (j - 1) in lna_positions:
                    key_xln = (seq_clean[j - 1], seq_clean[j])
                    if key_xln in _LNA_PARAMS_XLN:
                        dh_xln, ds_xln = _LNA_PARAMS_XLN[key_xln]
                        dh_total += dh_xln
                        ds_total += ds_xln
                        ddH_lna += dh_xln
                        ddS_lna += ds_xln
                if j in lna_positions:
                    key_mxl = (seq_clean[j - 1], seq_clean[j])
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
            if (i - 1) in pairs and pairs[i - 1] == j + 1 and not (include_loops and cut_idx != -1 and i == cut_idx):
                is_left_end = False
                
            is_right_end = True
            if (i + 1) in pairs and pairs[i + 1] == j - 1 and not (include_loops and cut_idx != -1 and (i + 1) == cut_idx):
                is_right_end = False
                
            pair_base_1, pair_base_2 = seq_clean[i], seq_clean[j]
            
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
            if seq_clean[i] in ('G', 'C') or seq_clean[j] in ('G', 'C'):
                gc_count += 1

    # Détection des renflements (bulges) et boucles internes entre paires emboîtées (i, j) -> (p, q)
    # Detect bulges and internal loops between directly nested base pairs (i, j) -> (p, q)
    ordered_pairs = sorted((i, j) for i, j in pairs.items() if i < j)
    has_loops = False
    for idx_p in range(len(ordered_pairs) - 1):
        i, j = ordered_pairs[idx_p]
        p, q = ordered_pairs[idx_p + 1]
        if not (i < p < q < j):
            continue
        if p == i + 1 and q == j - 1 and not (cut_idx != -1 and p == cut_idx):
            continue
        if cut_idx != -1 and ((i < cut_idx <= p) or (q < cut_idx <= j)):
            continue
        if any(k in pairs for k in range(i + 1, p)) or any(k in pairs for k in range(q + 1, j)):
            continue

        has_loops = True
        if not include_loops:
            continue

        l1 = p - i - 1
        l2 = j - q - 1
        # Retirer les deux initiations bimoléculaires fictives ajoutées aux bornes internes de la boucle
        # Remove the two spurious bimolecular initiations added at the internal loop boundaries
        for b1, b2 in ((seq_clean[i], seq_clean[j]), (seq_clean[p], seq_clean[q])):
            init_term = _INIT_GC if (b1, b2) in (('G', 'C'), ('C', 'G')) else _INIT_AT
            dh_total -= init_term[0]
            ds_total -= init_term[1]

        t37_k = 310.15
        if (l1 == 1 and l2 == 0) or (l1 == 0 and l2 == 1):
            # Renflement d'1 nt (SantaLucia & Hicks 2004, Table 4) : empilement NN conservé + pénalité +4,0 kcal/mol à 37 °C
            # Single-nt bulge (SantaLucia & Hicks 2004, Table 4): intervening NN stack retained + +4.0 kcal/mol at 37 °C
            stack_key = _SEQ_TO_NN_KEY.get(seq_clean[i] + seq_clean[p])
            if stack_key:
                dh_s, ds_s = _NN_PARAMS[stack_key]
                dh_total += dh_s
                ds_total += ds_s
            dh_total += 2.6
            ds_total += -4.5
        elif l1 == 0 or l2 == 0:
            # Renflement de k >= 2 nt (SantaLucia & Hicks 2004, Table 4 & extrapolation Jacobson-Stockmayer)
            # Multi-nt bulge k >= 2 (SantaLucia & Hicks 2004, Table 4 & Jacobson-Stockmayer extrapolation)
            k_bulge = l1 + l2
            if k_bulge == 2:
                dg37_bulge = 2.9
            elif k_bulge == 3:
                dg37_bulge = 3.1
            else:
                dg37_bulge = 3.2 + 2.44 * R_GAS_KCAL_MOL_K * t37_k * math.log(k_bulge / 3.0)
            # Pénalité de fermeture A·T (+0,5 kcal/mol par paire A·T adjacente, SantaLucia & Hicks 2004)
            # A·T closing penalty (+0.5 kcal/mol per closing A·T pair, SantaLucia & Hicks 2004)
            for b1, b2 in ((seq_clean[i], seq_clean[j]), (seq_clean[p], seq_clean[q])):
                if (b1, b2) not in (('G', 'C'), ('C', 'G')):
                    dg37_bulge += 0.5
            ds_total += -(dg37_bulge / t37_k) * 1000.0
        else:
            # Boucle interne l1 x l2 (Allawi & SantaLucia 1997 + SantaLucia & Hicks 2004 Table 5)
            # Internal loop l1 x l2 (Allawi & SantaLucia 1997 + SantaLucia & Hicks 2004 Table 5)
            from labcraft.thermo.mismatch import (
                _MISMATCH_INTERNAL,
                _MISMATCH_TERMINAL,
                _internal_loop_entropy_cal,
            )

            if l1 == 1 and l2 == 1:
                for dinuc_t, dinuc_b in (
                    (seq_clean[i] + seq_clean[i + 1], seq_clean[j] + seq_clean[j - 1]),
                    (seq_clean[i + 1] + seq_clean[p], seq_clean[j - 1] + seq_clean[q]),
                ):
                    key_fwd = f"{dinuc_t}/{dinuc_b}"
                    key_rev = f"{dinuc_b[::-1]}/{dinuc_t[::-1]}"
                    if key_fwd in _MISMATCH_INTERNAL:
                        dh_m, ds_m = _MISMATCH_INTERNAL[key_fwd]
                        dh_total += dh_m
                        ds_total += ds_m
                    elif key_rev in _MISMATCH_INTERNAL:
                        dh_m, ds_m = _MISMATCH_INTERNAL[key_rev]
                        dh_total += dh_m
                        ds_total += ds_m
            else:
                left_stack = f"{seq_clean[i]}{seq_clean[i + 1]}/{seq_clean[j]}{seq_clean[j - 1]}"
                right_stack = f"{seq_clean[q]}{seq_clean[q + 1]}/{seq_clean[p]}{seq_clean[p - 1]}"
                if left_stack in _MISMATCH_TERMINAL:
                    dh_l, ds_l = _MISMATCH_TERMINAL[left_stack]
                    dh_total += dh_l
                    ds_total += ds_l
                if right_stack in _MISMATCH_TERMINAL:
                    dh_r, ds_r = _MISMATCH_TERMINAL[right_stack]
                    dh_total += dh_r
                    ds_total += ds_r
                ds_total += _internal_loop_entropy_cal(l1 + l2)
                if l1 != l2:
                    ds_total += -(0.3 * abs(l1 - l2) / t37_k) * 1000.0

    # Compter aussi les épingles intramoléculaires comme structures à boucles pour l'auto-cohérence Gibbs-Helmholtz
    # Also flag intramolecular hairpins as loop-bearing structures for Gibbs-Helmholtz self-consistency
    if not has_loops:
        for i_hp, j_hp in ordered_pairs:
            if not ((i_hp + 1) in pairs and pairs[i_hp + 1] == j_hp - 1):
                unpaired_inside = j_hp - i_hp - 1
                if unpaired_inside > 0 and not (cut_idx != -1 and i_hp < cut_idx <= j_hp):
                    has_loops = True
                    break
                
    f_gc = gc_count / n_bp if n_bp > 0 else 0.0
    return dh_total, ds_total, n_bp, f_gc, ddH_lna, ddS_lna, has_loops


def estimate_helix_thermo(
    seq: str,
    structure: str,
    lna_positions: tuple[int, ...] = (),
    *,
    include_loops: bool = False,
) -> Tuple[float, float, int, float, float, float]:
    """Estime (dH_total, dS_total, n_bp, f_gc, ddH_lna, ddS_lna) à partir de la séquence et de la structure.
    Estimates (dH_total, dS_total, n_bp, f_gc, ddH_lna, ddS_lna) from sequence and dot-bracket structure.
    """
    dh, ds, n_bp, f_gc, ddh_lna, dds_lna, _ = _estimate_structure_thermo_impl(
        seq, structure, lna_positions, include_loops=include_loops
    )
    return dh, ds, n_bp, f_gc, ddh_lna, dds_lna


def estimate_structure_thermo(
    seq: str,
    structure: str,
    lna_positions: tuple[int, ...] = (),
) -> Tuple[float, float, int, float, float, float]:
    """Estime la thermodynamique complète incluant renflements (bulges) et boucles internes (SantaLucia 2004).
    Estimates full secondary-structure thermodynamics including bulges and internal loops (SantaLucia 2004).
    """
    return estimate_helix_thermo(seq, structure, lna_positions, include_loops=True)

_MAX_SALT_CACHE_SIZE = 65536
_SHARED_SALT_CACHE: Dict[Tuple[Any, ...], DuplexResult] = {}


class ViennaSaltShiftBackend(DuplexEnergyBackend):
    def __init__(
        self,
        salt_model: UnifiedSaltModel,
        default_ct_molar: float = 2e-6,
        *,
        couple_loops_in_salt_dg: bool = False,
    ) -> None:
        self.vienna = ViennaRNABackend()
        self.salt_model = salt_model
        self.default_ct_molar = default_ct_molar
        self.couple_loops_in_salt_dg = couple_loops_in_salt_dg
        if (
            type(salt_model) is UnifiedSaltModel
            and default_ct_molar == 2e-6
            and not couple_loops_in_salt_dg
        ):
            self._cache: Dict[Tuple[Any, ...], DuplexResult] = _SHARED_SALT_CACHE
        else:
            self._cache = {}
        
    def _apply_shift(self, res: DuplexResult, seq: str, ct_molar: float, kwargs: dict[str, Any], lna_positions: tuple[int, ...] = ()) -> DuplexResult:
        if not res.structure or '(' not in res.structure:
            return res
            
        dh_helix, ds_helix, n_bp, f_gc, ddH_lna, ddS_lna, has_loops = _estimate_structure_thermo_impl(
            seq, res.structure, lna_positions, include_loops=self.couple_loops_in_salt_dg
        )
        
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
        
        # Calcul du décalage LNA / Calculate ddG_LNA
        ddg_lna = ddH_lna - temp_k * (ddS_lna / 1000.0)
        
        # La condensation ionique (Owczarzy) s'applique sur les tiges empilées (ou structure complète si activé)
        # Ionic screening (Owczarzy) applies to stacked stems (or full structure when enabled)
        ddg_salt = -temp_k * (ds_corr - ds_helix) / 1000.0
        
        dg_final = res.dg_kcal + ddg_salt + ddg_lna

        # Pour les structures comportant des renflements (bulges) ou boucles internes, recalculer dH/dS/Tm
        # de manière auto-cohérente (Gibbs-Helmholtz : dG_final = dH_out - T * dS_out) afin que le Tm reflète
        # fidèlement la déstabilisation des boucles tout en préservant l'échelle dG_final calibrée.
        # For structures containing bulges or internal loops, compute Gibbs-Helmholtz self-consistent
        # dH/dS/Tm (dG_final = dH_out - T * dS_out) so Tm faithfully reflects loop destabilization
        # while preserving the calibrated dG_final equilibrium scale.
        if has_loops:
            if self.couple_loops_in_salt_dg:
                dh_out = dh_corr
            else:
                dh_struct, _, _, _, _, _, _ = _estimate_structure_thermo_impl(
                    seq, res.structure, lna_positions, include_loops=True
                )
                dh_out = dh_struct
            ds_out = ((dh_out - dg_final) / temp_k) * 1000.0 if temp_k > 0 else ds_corr
        else:
            dh_out = dh_corr
            ds_out = ds_corr
        
        # Recalcul cohérent du Tm / Coherent Tm recalculation
        r_gas = R_GAS_KCAL_MOL_K * 1000.0
        
        x = 4.0
        if '&' in seq:
            p1, p2 = seq.split('&', 1)
            if p1 == p2:
                x = 1.0
            
            if res.tm_celsius is not None and not math.isnan(res.tm_celsius):
                if ct_molar > 0:
                    denom = ds_out + r_gas * math.log(ct_molar / x)
                    new_tm_k = (dh_out * 1000.0) / denom if denom != 0.0 else float('nan')
                    new_tm_celsius = new_tm_k - 273.15
                else:
                    new_tm_celsius = res.tm_celsius
            else:
                new_tm_celsius = float('nan')
        else:
            if res.tm_celsius is not None and not math.isnan(res.tm_celsius) and ds_out != 0.0:
                new_tm_k = (dh_out * 1000.0) / ds_out
                new_tm_celsius = new_tm_k - 273.15
            else:
                new_tm_celsius = float('nan')
            
        return DuplexResult(
            dh_kcal=round(dh_out, 2),
            ds_cal_per_k=round(ds_out, 2),
            dg_kcal=dg_final,
            tm_celsius=new_tm_celsius,
            structure=res.structure,
            temperature_celsius=res.temperature_celsius
        )

    def calc_heterodimer(self, seq1: str, seq2: str, *, temp_celsius: float = 65.0, **kwargs: Any) -> DuplexResult:
        cache_key = ("hetero", seq1, seq2, temp_celsius, tuple(sorted((k, tuple(v) if isinstance(v, (list, tuple)) else v) for k, v in kwargs.items())))
        if cache_key in self._cache:
            return self._cache[cache_key]
        res = self.vienna.calc_heterodimer(seq1, seq2, temp_celsius=temp_celsius, **kwargs)
        ct = kwargs.get('ct_molar', self.default_ct_molar)
        if ct is None: ct = self.default_ct_molar
        lna_pos_a = kwargs.get('lna_positions_a', ())
        lna_pos_b = kwargs.get('lna_positions_b', ())
        # Remap lna_positions on concatenated sequence 'seq1seq2' (without &)
        mapped_lna = list(lna_pos_a) + [pos + len(seq1) for pos in lna_pos_b]
        shifted = self._apply_shift(res, f"{seq1}&{seq2}", ct, kwargs, tuple(mapped_lna))
        if len(self._cache) >= _MAX_SALT_CACHE_SIZE:
            self._cache.clear()
        self._cache[cache_key] = shifted
        return shifted
        
    def calc_homodimer(self, seq: str, *, temp_celsius: float = 65.0, **kwargs: Any) -> DuplexResult:
        cache_key = ("homo", seq, temp_celsius, tuple(sorted((k, tuple(v) if isinstance(v, (list, tuple)) else v) for k, v in kwargs.items())))
        if cache_key in self._cache:
            return self._cache[cache_key]
        res = self.vienna.calc_homodimer(seq, temp_celsius=temp_celsius, **kwargs)
        ct = kwargs.get('ct_molar', self.default_ct_molar)
        if ct is None: ct = self.default_ct_molar
        lna_pos = kwargs.get('lna_positions', ())
        mapped_lna = list(lna_pos) + [pos + len(seq) for pos in lna_pos]
        shifted = self._apply_shift(res, f"{seq}&{seq}", ct, kwargs, tuple(mapped_lna))
        if len(self._cache) >= _MAX_SALT_CACHE_SIZE:
            self._cache.clear()
        self._cache[cache_key] = shifted
        return shifted
        
    def calc_hairpin(self, seq: str, *, temp_celsius: float = 65.0, **kwargs: Any) -> DuplexResult:
        """Évalue une épingle intramoléculaire à la force ionique équivalente sodium de ViennaRNA + LNA.
        Evaluates a monomolecular hairpin at ViennaRNA's sodium-equivalent ionic strength + LNA shift.
        """
        cache_key = ("hairpin", seq, temp_celsius, tuple(sorted((k, tuple(v) if isinstance(v, (list, tuple)) else v) for k, v in kwargs.items())))
        if cache_key in self._cache:
            return self._cache[cache_key]

        from labcraft.buffer.monovalent import get_total_monovalent
        from labcraft.thermo.salt import sodium_equivalent_for_folding

        eff_mon: Any = kwargs.get("mon_molar")
        if eff_mon is None:
            na_mm = float(kwargs.get("na_mm", 50.0))
            k_mm = float(kwargs.get("k_mm", 0.0))
            tris_mm = float(kwargs.get("tris_mm", 0.0))
            mg_mm = float(kwargs.get("mg_mm", 0.0))
            dntp_mm = float(kwargs.get("dntp_mm", 0.0) or 0.0)
            if abs(na_mm - 1000.0) < 1e-6 and k_mm == 0.0 and tris_mm == 0.0 and mg_mm == 0.0:
                eff_mon = None
            else:
                mon_true = get_total_monovalent(na_mm, k_mm, tris_mm) / 1000.0
                eff_mon = sodium_equivalent_for_folding(mon_true, mg_mm / 1000.0, dntp_mm / 1000.0)

        res = self.vienna.calc_hairpin(seq, temp_celsius=temp_celsius, mon_molar=eff_mon)
        if not res.structure or "(" not in res.structure:
            if len(self._cache) >= _MAX_SALT_CACHE_SIZE:
                self._cache.clear()
            self._cache[cache_key] = res
            return res

        lna_pos = tuple(kwargs.get("lna_positions", ()))
        dh_st, _, _, _, ddH_lna, ddS_lna, _ = _estimate_structure_thermo_impl(
            seq, res.structure, lna_pos, include_loops=True
        )
        temp_k = temp_celsius + 273.15
        ddg_lna = ddH_lna - temp_k * (ddS_lna / 1000.0)
        dg_final = res.dg_kcal + ddg_lna
        ds_st = ((dh_st - dg_final) / temp_k) * 1000.0 if temp_k > 0 else 0.0
        tm_hp = (dh_st * 1000.0) / ds_st - 273.15 if ds_st != 0.0 else float("nan")
        shifted = DuplexResult(
            dh_kcal=round(dh_st, 2),
            ds_cal_per_k=round(ds_st, 2),
            dg_kcal=dg_final,
            tm_celsius=tm_hp,
            structure=res.structure,
            temperature_celsius=temp_celsius,
        )
        if len(self._cache) >= _MAX_SALT_CACHE_SIZE:
            self._cache.clear()
        self._cache[cache_key] = shifted
        return shifted
        
    def calc_duplex(self, seq1: str, seq2: str, *, temp_celsius: float = 65.0, **kwargs: Any) -> DuplexResult:
        return self.calc_heterodimer(seq1, seq2, temp_celsius=temp_celsius, **kwargs)
