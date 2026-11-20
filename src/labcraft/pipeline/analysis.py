"""Orchestration scientifique de l'analyse thermodynamique de panels d'amorces.

Scientific orchestration pipeline for multiplex primer panel thermodynamic analysis.
"""
from __future__ import annotations

import copy
import importlib.metadata
import time
import warnings
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

from labcraft.cli.config import PanelConfig, PrimerSetConfig, build_engine_from_config
from labcraft.cli.parsers import read_multi_fasta, validate_target_sequence
from labcraft.diagnostics.amplifiable_dimer import (
    enrich_dimer_boltzmann_3p,
    evaluate_pair_amplifiable,
)
from labcraft.diagnostics.enzyme import PolymeraseProfile
from labcraft.diagnostics.mispriming import MisprimingRisk, detect_inter_target_mispriming
from labcraft.diagnostics.probe_tm import check_probes_tm
from labcraft.lamp.complex_enumeration import enumerate_complexes
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.lamp.stoichiometry import LAMP_DEFAULT_PROFILE
from labcraft.metrics.balance import calculate_multiplex_balance
from labcraft.metrics.fractions import PrimerFractions, compute_fractions
from labcraft.metrics.risk import (
    RiskItem,
    evaluate_risks,
    summarize_amplifiable_dimers,
)
from labcraft.metrics.verdict import (
    DEFAULT_DIMER_WARN_RATIO,
    DEFAULT_ROLE_IMBALANCE_RATIO,
    PanelVerdict,
    generate_verdict,
)
from labcraft.optimize.concentrations import optimize_concentrations
from labcraft.optimize.recommendations import generate_recommendations
from labcraft.report.alignment import dotbracket_to_alignment, get_alignment_columns
from labcraft.report.json_export import build_analyze_json_payload
from labcraft.solver.dual import solve_dual
from labcraft.solver.types import EquilibriumProblem, EquilibriumResult
from labcraft.thermo.backends.base import DuplexEnergyBackend
from labcraft.thermo.vienna import dna_session


class EmptyPanelError(ValueError):
    """Erreur levée lorsqu'un panel d'amorces est vide.
    Raised when a primer panel contains zero primers.
    """


class TargetFileError(RuntimeError):
    """Erreur levée lors de la lecture ou validation d'un fichier cible.
    Raised when reading or validating a target sequence file fails.
    """

    def __init__(self, filepath: str, cause: Exception) -> None:
        super().__init__(f"{filepath}: {cause}")
        self.filepath = filepath
        self.cause = cause


@dataclass
class PanelAnalysisResult:
    """Résultat complet de l'analyse thermodynamique d'un panel ou multiplexe.
    Complete thermodynamic analysis result for a single panel or joint multiplex.
    """

    config: PanelConfig
    panel_label: Optional[str]
    primers: List[PhysicalPrimer]
    primer_to_target: Dict[str, str]
    primer_to_panel_map: Dict[str, str]
    effective_concs_uM: Dict[str, float]
    temp_celsius: float
    enzyme: PolymeraseProfile
    has_true_target: bool
    fractions: Dict[str, PrimerFractions]
    target_occupations: Dict[str, Dict[str, float]]
    unfolding_penalties: Dict[str, Any]
    risks: List[RiskItem]
    dimer_decomposition: Dict[str, Any]
    interaction_matrix: Dict[str, Dict[str, float]]
    verdict: PanelVerdict
    recommendations: List[str]
    probe_tm_results: List[Dict[str, Any]]
    optimization_results: Union[Dict[str, Any], List[Any]]
    panel_summaries: Dict[str, Dict[str, Any]]
    balance_cv: Optional[float]
    mispriming_risks: List[MisprimingRisk]
    max_residual: float
    warnings: List[str]
    metadata: Dict[str, Any]
    json_payload: Dict[str, Any]


def _get_parent_name(name: str) -> str:
    """Retourne le nom d'amorce parente sans suffixe de variant `#k`.
    Returns the parent primer name stripped of any `#k` variant suffix.
    """
    return name.split("#")[0] if "#" in name else name


def _format_nm(value: float) -> str:
    """Formate une concentration en nM avec précision adaptée aux traces.
    Formats a concentration in nM with extra precision below 0.01 nM.
    """
    return f"{value:.4f}" if 0.0 < value < 0.01 else f"{value:.2f}"


DEFAULT_BOLTZMANN_MAX_DIMERS: Optional[int] = None
DEFAULT_BOLTZMANN_DEDUP_3P_WINDOW: int = 0
_BoltzmannDedupKey = Tuple[Tuple[str, str], Tuple[str, str]]


def _make_boltzmann_dedup_key(
    pa_par: str,
    pa_name: str,
    seq_a: str,
    pb_par: str,
    pb_name: str,
    seq_b: str,
    window_3p: int,
) -> _BoltzmannDedupKey:
    """Construit la clé canonique de regroupement Boltzmann :
    - `window_3p <= 0` (défaut `0`, désactivé) : clé exacte par variant physique `((name_1, seq_1), (name_2, seq_2))`
      sans aucune fusion (évaluation exhaustive) ;
    - `window_3p > 0` : APPROXIMATION `((parent_1, 3'_1), (parent_2, 3'_2))` sur `window_3p` nt (peut masquer un
      variant à risque lorsque l'extrémité 3' s'apparie sur une base dégénérée interne du brin partenaire,
      cf. `LB_3#2_homo` à `0,00 %` face à `LB_3#1_homo` à `42,76 %` et `LB_3#1_LB_3#2` à `43,20 %`).

    Builds the canonical Boltzmann grouping key:
    - `window_3p <= 0` (default `0`, disabled): exact per-variant key `((name_1, seq_1), (name_2, seq_2))` with
      zero merging (exhaustive evaluation);
    - `window_3p > 0`: APPROXIMATION `((parent_1, 3'_1), (parent_2, 3'_2))` over `window_3p` nt (may mask risky
      variants when the 3' end base-pairs against an internal degenerate base on the partner strand, e.g.
      `LB_3#2_homo` at `0.00%` vs `LB_3#1_homo` at `42.76%` and `LB_3#1_LB_3#2` at `43.20%`).
    """
    w = int(window_3p)
    if w <= 0:
        end_a = (pa_name, seq_a)
        end_b = (pb_name, seq_b)
    else:
        end_a = (pa_par, seq_a[-w:])
        end_b = (pb_par, seq_b[-w:])
    return (end_a, end_b) if end_a <= end_b else (end_b, end_a)


def evaluate_dimer_complexes(
    prob: EquilibriumProblem,
    res: EquilibriumResult,
    complexes: List[str],
    primers: List[PhysicalPrimer],
    backend: DuplexEnergyBackend,
    enzyme: PolymeraseProfile,
    temp_celsius: float,
    backend_kwargs: Dict[str, Any],
    primer_to_panel_map: Dict[str, str],
    *,
    dimer_conc_max_nm: float = 0.23,
    boltzmann_max_dimers: Optional[int] = DEFAULT_BOLTZMANN_MAX_DIMERS,
    boltzmann_dedup_3p_window: int = DEFAULT_BOLTZMANN_DEDUP_3P_WINDOW,
    boltzmann_stats_out: Optional[Dict[str, Any]] = None,
) -> Tuple[List[bool], List[float], List[Dict[str, Any]]]:
    """Évalue chaque complexe dimérique sans cible dans une session ADN unique.
    Evaluates each target-free dimer complex within a single DNA parameter session.
    """
    with dna_session(temp_celsius=temp_celsius):
        return _evaluate_dimer_complexes_impl(
            prob,
            res,
            complexes,
            primers,
            backend,
            enzyme,
            temp_celsius,
            backend_kwargs,
            primer_to_panel_map,
            dimer_conc_max_nm=dimer_conc_max_nm,
            boltzmann_max_dimers=boltzmann_max_dimers,
            boltzmann_dedup_3p_window=boltzmann_dedup_3p_window,
            boltzmann_stats_out=boltzmann_stats_out,
        )


def _evaluate_dimer_complexes_impl(
    prob: EquilibriumProblem,
    res: EquilibriumResult,
    complexes: List[str],
    primers: List[PhysicalPrimer],
    backend: DuplexEnergyBackend,
    enzyme: PolymeraseProfile,
    temp_celsius: float,
    backend_kwargs: Dict[str, Any],
    primer_to_panel_map: Dict[str, str],
    *,
    dimer_conc_max_nm: float = 0.23,
    boltzmann_max_dimers: Optional[int] = DEFAULT_BOLTZMANN_MAX_DIMERS,
    boltzmann_dedup_3p_window: int = DEFAULT_BOLTZMANN_DEDUP_3P_WINDOW,
    boltzmann_stats_out: Optional[Dict[str, Any]] = None,
) -> Tuple[List[bool], List[float], List[Dict[str, Any]]]:
    """Évalue chaque complexe dimérique sans cible pour l'extensibilité 3' et la décomposition.
    Evaluates each target-free dimer complex for 3' extensibility and panel decomposition.
    """
    amplifiable_flags: List[bool] = []
    concs: List[float] = []
    dimer_details: List[Dict[str, Any]] = []
    dimer_records: List[
        Tuple[
            int,
            str,
            PhysicalPrimer,
            PhysicalPrimer,
            bool,
            float,
            Tuple[str, str],
            _BoltzmannDedupKey,
        ]
    ] = []

    for i, c_name in enumerate(complexes):
        conc = float(res.concentrations[i])
        concs.append(conc)

        is_amp = False
        details: Dict[str, Any] = {}
        has_target = any(
            k >= len(primers) and val > 0 for k, val in enumerate(prob.stoichiometry[i])
        )

        if not has_target:
            p_counts = {
                k: int(prob.stoichiometry[i][k])
                for k, val in enumerate(prob.stoichiometry[i])
                if val > 0 and k < len(primers)
            }
            if sum(p_counts.values()) == 2:
                ks = list(p_counts.keys())
                p_a = primers[ks[0]]
                p_b = primers[ks[0]] if len(ks) == 1 else primers[ks[1]]
                is_amp, _dg_3p, details = evaluate_pair_amplifiable(
                    p_a, p_b, backend, enzyme, temp_celsius, **backend_kwargs
                )
                pa_par = str(details.get("parent_a") or _get_parent_name(p_a.name))
                pb_par = str(details.get("parent_b") or _get_parent_name(p_b.name))
                parent_pair = (pa_par, pb_par) if pa_par <= pb_par else (pb_par, pa_par)
                dedup_key = _make_boltzmann_dedup_key(
                    pa_par,
                    p_a.name,
                    p_a.sequence,
                    pb_par,
                    p_b.name,
                    p_b.sequence,
                    boltzmann_dedup_3p_window,
                )
                dimer_records.append((i, c_name, p_a, p_b, is_amp, conc, parent_pair, dedup_key))
                details["panel_a"] = primer_to_panel_map.get(
                    pa_par, primer_to_panel_map.get(p_a.name, "")
                )
                details["panel_b"] = primer_to_panel_map.get(
                    pb_par, primer_to_panel_map.get(p_b.name, "")
                )

                if details and "alignment" not in details:
                    details["alignment"] = dotbracket_to_alignment(
                        details["seq_a"], details["seq_b"], details["structure"]
                    )
                if details and "alignment_columns" not in details:
                    cols = get_alignment_columns(
                        details["seq_a"],
                        details["seq_b"],
                        details["structure"],
                        details.get("extensible_strand"),
                    )
                    details["alignment_columns"] = cols
                    three_prime_idx = next(
                        (idx for idx, c in enumerate(cols) if c.get("role") == "three_prime"),
                        -1,
                    )
                    template_count = sum(1 for c in cols if c.get("role") == "template")
                    before_arrow = (
                        sum(1 for c in cols[:three_prime_idx] if not c.get("is_truncation"))
                        if three_prime_idx != -1
                        else 0
                    )
                    details["arrow_metrics"] = {
                        "show": three_prime_idx != -1 and template_count > 0,
                        "margin_cols": before_arrow,
                        "width_cols": template_count,
                    }
        amplifiable_flags.append(is_amp)
        dimer_details.append(details)

    # `conc` du solveur est en mol/L (M) ; `dimer_conc_max_nm` est en nmol/L (1 nM = 1e-9 M).
    # Enrichissement Boltzmann 3' (informatif uniquement, sans effet sur le verdict) :
    # 1. Par défaut (`boltzmann_dedup_3p_window = 0`, `boltzmann_max_dimers = None`), l'évaluation est
    #    EXHAUSTIVE : 100 % des dimères amplifiables + 100 % des dimères non amplifiables vérifiant
    #    `conc >= 0.1 * dimer_conc_max` (0,023 nM = 2,3e-11 M par défaut), sans dédoublonnage ni plafond ;
    # 2. Accélérateurs optionnels (APPROXIMATIONS pour très gros panels, pouvant masquer un variant à risque
    #    tel que `LB_3#2_homo` à 0,00 % face à `LB_3#1_homo` à 42,76 % lorsque la région matrice interne du
    #    brin partenaire porte une base dégénérée) : `boltzmann_dedup_3p_window > 0` fusionne par clé
    #    `((parent_A, 3'_A), (parent_B, 3'_B))` et `boltzmann_max_dimers > 0` plafonne à K entrées ;
    # 3. Pour chaque couple parent, la `P_Boltzmann` maximale sur ses représentants évalués est rapportée.
    # Solver `conc` is in mol/L (M); `dimer_conc_max_nm` is in nmol/L (1 nM = 1e-9 M).
    # Boltzmann 3' enrichment (informational only, never alters the verdict):
    # 1. By default (`boltzmann_dedup_3p_window = 0`, `boltzmann_max_dimers = None`), evaluation is
    #    EXHAUSTIVE: 100% of amplifiable dimers + 100% of non-amplifiable dimers with `conc >= 0.1 * dimer_conc_max`;
    # 2. Optional accelerators (APPROXIMATIONS for very large panels that can mask risky variants when the
    #    partner strand's internal template carries a degenerate base, e.g. `LB_3#2_homo` at 0.00% vs
    #    `LB_3#1_homo` at 42.76%): `boltzmann_dedup_3p_window > 0` and `boltzmann_max_dimers > 0`;
    # 3. For each parent pair, the maximum `P_Boltzmann` across its evaluated representatives is reported.
    boltzmann_min_conc_molar = 0.1 * dimer_conc_max_nm * 1e-9
    amp_records = [rec for rec in dimer_records if rec[4]]
    eligible_non_amp = [
        rec for rec in dimer_records if (not rec[4]) and rec[5] >= boltzmann_min_conc_molar
    ]

    best_by_dedup_key: Dict[
        _BoltzmannDedupKey,
        Tuple[
            int,
            str,
            PhysicalPrimer,
            PhysicalPrimer,
            bool,
            float,
            Tuple[str, str],
            _BoltzmannDedupKey,
        ],
    ] = {}
    for rec in eligible_non_amp:
        dkey = rec[7]
        prev = best_by_dedup_key.get(dkey)
        if prev is None or (rec[5], -rec[0]) > (prev[5], -prev[0]):
            best_by_dedup_key[dkey] = rec

    dedup_reps = list(best_by_dedup_key.values())
    dedup_reps.sort(key=lambda rec: (-rec[5], rec[0]))
    if boltzmann_max_dimers is None or int(boltzmann_max_dimers) <= 0:
        selected_non_amp = dedup_reps
        skipped_dedup_reps = []
    else:
        cap = int(boltzmann_max_dimers)
        selected_non_amp = dedup_reps[:cap]
        skipped_dedup_reps = dedup_reps[cap:]

    parent_max_prob: Dict[Tuple[str, str], float] = {}
    for idx_rec, _cname, p_a, p_b, _is_amp, _conc_rec, ppair, _dkey in amp_records:
        enrich_dimer_boltzmann_3p(
            dimer_details[idx_rec], p_a, p_b, True, temp_celsius, enzyme=enzyme, **backend_kwargs
        )
        p_val = float(dimer_details[idx_rec].get("boltzmann_3p_prob") or 0.0)
        if p_val > parent_max_prob.get(ppair, -1.0):
            parent_max_prob[ppair] = p_val

    for idx_rec, _cname, p_a, p_b, _is_amp, _conc_rec, ppair, _dkey in selected_non_amp:
        enrich_dimer_boltzmann_3p(
            dimer_details[idx_rec], p_a, p_b, False, temp_celsius, enzyme=enzyme, **backend_kwargs
        )
        p_val = float(dimer_details[idx_rec].get("boltzmann_3p_prob") or 0.0)
        if p_val > parent_max_prob.get(ppair, -1.0):
            parent_max_prob[ppair] = p_val

    for idx_rec, _cname, _p_a, _p_b, _is_amp, _conc_rec, ppair, _dkey in (
        amp_records + selected_non_amp
    ):
        if ppair in parent_max_prob:
            dimer_details[idx_rec]["parent_pair_max_boltzmann_3p_prob"] = parent_max_prob[ppair]

    max_prob_by_parent_str: Dict[str, float] = {
        f"{pp[0]} × {pp[1]}": round(val, 6)
        for pp, val in sorted(parent_max_prob.items(), key=lambda kv: (-kv[1], kv[0]))
    }

    boltzmann_mode = "dedup_approx" if int(boltzmann_dedup_3p_window) > 0 else "exhaustive"
    if boltzmann_stats_out is not None:
        selected_names = {rec[1] for rec in selected_non_amp}
        elig_ppairs = {rec[6] for rec in dedup_reps}
        eval_ppairs = {rec[6] for rec in selected_non_amp}
        boltzmann_stats_out["boltzmann_mode"] = boltzmann_mode
        boltzmann_stats_out["evaluated_amp_dimers"] = len(amp_records)
        boltzmann_stats_out["eligible_non_amp_dimers"] = len(eligible_non_amp)
        boltzmann_stats_out["eligible_dedup_groups"] = len(dedup_reps)
        boltzmann_stats_out["skipped_dedup_groups"] = len(skipped_dedup_reps)
        boltzmann_stats_out["eligible_parent_pairs"] = len(elig_ppairs)
        boltzmann_stats_out["evaluated_non_amp_dimers"] = len(selected_non_amp)
        boltzmann_stats_out["skipped_parent_pairs"] = len(elig_ppairs - eval_ppairs)
        boltzmann_stats_out["skipped_eligible_dimers"] = len(eligible_non_amp) - len(
            selected_non_amp
        )
        boltzmann_stats_out["eligible_dedup_keys"] = {rec[7] for rec in dedup_reps}
        boltzmann_stats_out["evaluated_dedup_keys"] = {rec[7] for rec in selected_non_amp}
        boltzmann_stats_out["eligible_parent_pair_keys"] = elig_ppairs
        boltzmann_stats_out["evaluated_parent_pair_keys"] = eval_ppairs
        boltzmann_stats_out["max_prob_by_parent_pair"] = max_prob_by_parent_str
        boltzmann_stats_out["evaluated_names"] = selected_names
        boltzmann_stats_out["skipped_names"] = {
            rec[1] for rec in eligible_non_amp if rec[1] not in selected_names
        }

    return amplifiable_flags, concs, dimer_details


def load_targets_from_config(
    config_obj: PanelConfig,
    explicit_config: bool = True,
    progress_callback: Optional[Callable[[str], None]] = None,
) -> Tuple[Dict[str, str], List[str]]:
    """Charge et valide les séquences cibles référencées par un `PanelConfig`.
    Loads and validates target sequences referenced by a `PanelConfig`.
    """
    targets_dict: Dict[str, str] = {}
    target_warnings: List[str] = []

    for t in config_obj.targets:
        try:
            records = read_multi_fasta(t.sequence_file)
            for h, seq in records:
                t_id = h.split()[0]
                if t_id == t.id or not explicit_config:
                    if t_id in targets_dict:
                        continue
                    clean_seq, w_msg = validate_target_sequence(t_id, seq, warn_dashes=True)
                    targets_dict[t_id] = clean_seq
                    if w_msg:
                        if progress_callback:
                            progress_callback(f"Avertissement : {w_msg}")
                        target_warnings.append(w_msg)
        except Exception as exc:
            raise TargetFileError(t.sequence_file, exc) from exc

    if not targets_dict:
        targets_dict = {"Compétition Sans Cible": ""}

    return targets_dict, target_warnings


def split_colliding_panel_configs(
    config_obj: PanelConfig,
    force_joint: bool = False,
) -> Tuple[List[PanelConfig], bool]:
    """Détecte les panels alternatifs partageant une même cible et les sépare sauf `force_joint`.
    Detects alternative panels sharing a target and splits them unless `force_joint` is set.
    """
    target_to_panels: Dict[str, List[PrimerSetConfig]] = {}
    for pset in config_obj.primer_sets:
        target_to_panels.setdefault(pset.target, []).append(pset)

    has_collisions = any(len(panels) > 1 for panels in target_to_panels.values())
    if has_collisions and not force_joint:
        config_list: List[PanelConfig] = []
        for pset in config_obj.primer_sets:
            new_config = copy.deepcopy(config_obj)
            new_config.primer_sets = [pset]
            config_list.append(new_config)
        return config_list, True

    return [config_obj], False


def run_single_panel_analysis(
    current_config: PanelConfig,
    targets_dict: Dict[str, str],
    *,
    panel_label: Optional[str] = None,
    target_warnings: Optional[List[str]] = None,
    inputs_fp: Optional[Dict[str, Any]] = None,
    file_hash: str = "in-memory",
    max_variants: int = 64,
    errors: int = 2,
    strict_3prime: int = 3,
    strict_3prime_tolerate: int = 2,
    dimer_conc_max: float = 0.23,
    dimer_warn_ratio: float = DEFAULT_DIMER_WARN_RATIO,
    role_imbalance_ratio: float = DEFAULT_ROLE_IMBALANCE_RATIO,
    boltzmann_max_dimers: Optional[int] = DEFAULT_BOLTZMANN_MAX_DIMERS,
    boltzmann_dedup_3p_window: int = DEFAULT_BOLTZMANN_DEDUP_3P_WINDOW,
    separate_targets: bool = False,
    progress_callback: Optional[Callable[[str], None]] = None,
) -> PanelAnalysisResult:
    """Exécute l'analyse thermodynamique complète d'un `PanelConfig` unitaire ou conjoint.
    Runs the complete thermodynamic analysis of a single or joint `PanelConfig`.
    """
    emit = progress_callback or (lambda _msg: None)
    has_true_target = any(k != "Compétition Sans Cible" for k in targets_dict.keys())

    (
        primers,
        primer_to_target,
        backend,
        backend_kwargs,
        mon_molar,
        enzyme,
        temp_celsius,
        profiles,
    ) = build_engine_from_config(current_config, targets_dict, max_variants=max_variants)

    primer_to_panel = primer_to_target
    primer_to_panel_map: Dict[str, str] = {}
    for ps in current_config.primer_sets:
        panel_id = ps.panel_name if ps.panel_name else ps.target
        for role_name in ps.primers:
            primer_to_panel_map[f"{role_name}_{panel_id}"] = panel_id
    for p in primers:
        p_par = p.parent_name or _get_parent_name(p.name)
        if p_par in primer_to_panel_map:
            primer_to_panel_map[p.name] = primer_to_panel_map[p_par]
        elif p.name not in primer_to_panel_map:
            primer_to_panel_map[p.name] = primer_to_target.get(p.name, "")

    buf_cfg = current_config.experiment.buffer
    if buf_cfg:
        emit(
            f"Conditions utilisées : T={temp_celsius:.1f} °C, Enzyme={enzyme.name}, "
            f"Na={buf_cfg.na_mM} mM, Mg={buf_cfg.mg_mM} mM, dNTP={buf_cfg.dntp_mM} mM"
        )
    else:
        emit(
            f"Conditions utilisées : T={temp_celsius:.1f} °C, Enzyme={enzyme.name}, "
            "Référence (1.0 M Na+, 0 mM Mg2+)"
        )

    if not primers:
        raise EmptyPanelError(
            "Erreur : Le panel d'amorces est vide. Veuillez définir au moins un jeu d'amorces (primer_sets)."
        )

    # Concentration effective de chaque amorce (somme des variants dégénérés), en µM
    # Effective concentration of each primer (sum over degenerate variants), in uM
    effective_concs_uM: Dict[str, float] = {}
    for p in primers:
        p_par = p.parent_name or _get_parent_name(p.name)
        effective_concs_uM[p_par] = (
            effective_concs_uM.get(p_par, 0.0) + float(p.nominal_concentration or 0.0) * 1e6
        )
    emit(
        "Concentrations effectives par amorce (µM) : "
        + ", ".join(f"{n}={v:.3g}" for n, v in effective_concs_uM.items())
    )

    # Simulation d'équilibre couplé ou séparé dans une session ADN active
    # Coupled or per-target equilibrium simulation within an active DNA session
    target_occupations: Dict[str, Dict[str, float]] = {}
    all_fractions: Dict[str, PrimerFractions] = {}
    all_risks: List[RiskItem] = []
    max_residual_global = 0.0
    all_unfolding_penalties: Dict[str, Any] = {}
    interaction_matrix: Dict[str, Dict[str, float]] = {}
    boltzmann_eligible_dedup_groups: int = 0
    boltzmann_skipped_dedup_groups: int = 0
    boltzmann_eligible_parent_pairs: int = 0
    boltzmann_skipped_parent_pairs: int = 0
    boltzmann_skipped_eligible_dimers: int = 0
    boltzmann_max_prob_by_parent_pair: Dict[str, float] = {}

    prob: Optional[EquilibriumProblem] = None
    strands: Optional[List[str]] = None
    complexes: Optional[List[str]] = None

    buf_dump = buf_cfg.model_dump() if buf_cfg else None

    with dna_session(temp_celsius=temp_celsius):
        with warnings.catch_warnings(record=True) as captured_warnings:
            warnings.simplefilter("always")

            if not separate_targets and len(targets_dict) > 1:
                # Mode CONJOINT (par défaut en multiplexe multi-cibles)
                # JOINT mode (default in multi-target multiplex)
                emit(
                    f"Building joint complex network and solving coupled thermodynamic equilibrium "
                    f"for {len(targets_dict)} targets ({', '.join(targets_dict.keys())})..."
                )

                prob, strands, complexes, unfolding_penalties = enumerate_complexes(
                    primers,
                    targets_dict if has_true_target else None,
                    backend,
                    profile=profiles,
                    temp_celsius=temp_celsius,
                    mon_molar=mon_molar,
                    buffer=buf_dump,
                    enzyme=enzyme,
                    max_errors=errors,
                    strict_3prime_len=strict_3prime,
                    strict_3prime_tolerate=strict_3prime_tolerate,
                    primer_to_target=primer_to_target,
                )

                res = solve_dual(prob)
                max_residual_global = res.max_residual

                boltz_stats_joint: Dict[str, Any] = {}
                amplifiable_flags, concs, dimer_details = evaluate_dimer_complexes(
                    prob,
                    res,
                    complexes,
                    primers,
                    backend,
                    enzyme,
                    temp_celsius,
                    backend_kwargs,
                    primer_to_panel_map,
                    dimer_conc_max_nm=dimer_conc_max,
                    boltzmann_max_dimers=boltzmann_max_dimers,
                    boltzmann_dedup_3p_window=boltzmann_dedup_3p_window,
                    boltzmann_stats_out=boltz_stats_joint,
                )
                boltzmann_eligible_dedup_groups = int(
                    boltz_stats_joint.get("eligible_dedup_groups", 0)
                )
                boltzmann_skipped_dedup_groups = int(
                    boltz_stats_joint.get("skipped_dedup_groups", 0)
                )
                boltzmann_eligible_parent_pairs = int(
                    boltz_stats_joint.get("eligible_parent_pairs", 0)
                )
                boltzmann_skipped_parent_pairs = int(
                    boltz_stats_joint.get("skipped_parent_pairs", 0)
                )
                boltzmann_skipped_eligible_dimers = int(
                    boltz_stats_joint.get("skipped_eligible_dimers", 0)
                )
                boltzmann_max_prob_by_parent_pair = dict(
                    boltz_stats_joint.get("max_prob_by_parent_pair", {})
                )

                risks = evaluate_risks(
                    complexes,
                    concs,
                    amplifiable_flags,
                    is_warm_start=enzyme.warm_start,
                    dimer_details=dimer_details,
                )
                fractions = compute_fractions(
                    strands,
                    complexes,
                    prob.stoichiometry,
                    res.free_concentrations,
                    prob.delta_g,
                    temp_celsius,
                    primer_to_panel,
                )

                all_risks = risks
                all_fractions = fractions

                for t_id in targets_dict.keys():
                    target_occupations[t_id] = {}

                for s in strands:
                    if s.endswith("_site"):
                        idx = strands.index(s)
                        tot_s = float(prob.total_concentrations[idx])
                        free_s = float(res.free_concentrations[idx])
                        occ = (tot_s - free_s) / tot_s if tot_s > 0 else 0.0
                        if "@" in s:
                            p_part, rest_part = s.split("@", 1)
                            t_part = rest_part.replace("_site", "")
                            if t_part in target_occupations:
                                target_occupations[t_part][s] = occ
                                target_occupations[t_part][f"{p_part}_site"] = occ
                        else:
                            for t_id in targets_dict.keys():
                                target_occupations[t_id][s] = occ

                if isinstance(unfolding_penalties, dict) and any(
                    isinstance(v, dict) for v in unfolding_penalties.values()
                ):
                    all_unfolding_penalties = unfolding_penalties
                else:
                    for t_id in targets_dict.keys():
                        all_unfolding_penalties[t_id] = unfolding_penalties
            else:
                # Mode Séparé (ou cible unique)
                # Separate mode (or single target)
                sep_evaluated_names: set[str] = set()
                sep_skipped_names: set[str] = set()
                sep_eligible_dkeys: set[_BoltzmannDedupKey] = set()
                sep_evaluated_dkeys: set[_BoltzmannDedupKey] = set()
                sep_eligible_ppairs: set[Tuple[str, str]] = set()
                sep_evaluated_ppairs: set[Tuple[str, str]] = set()
                sep_max_prob_by_parent: Dict[str, float] = {}
                for t_id, t_seq in targets_dict.items():
                    emit(f"Solving for target {t_id}...")
                    profile = profiles.get(t_id)
                    if profile is None and profiles:
                        profile = next(iter(profiles.values()))
                    eff_profile = profile if profile is not None else LAMP_DEFAULT_PROFILE

                    prob, strands, complexes, unfolding_penalties = enumerate_complexes(
                        primers,
                        t_seq if t_seq else None,
                        backend,
                        profile=eff_profile,
                        temp_celsius=temp_celsius,
                        mon_molar=mon_molar,
                        buffer=buf_dump,
                        enzyme=enzyme,
                        max_errors=errors,
                        strict_3prime_len=strict_3prime,
                        strict_3prime_tolerate=strict_3prime_tolerate,
                        primer_to_target=primer_to_target,
                    )

                    res = solve_dual(prob)
                    max_residual_global = max(max_residual_global, res.max_residual)

                    boltz_stats_sep: Dict[str, Any] = {}
                    amplifiable_flags, concs, dimer_details = evaluate_dimer_complexes(
                        prob,
                        res,
                        complexes,
                        primers,
                        backend,
                        enzyme,
                        temp_celsius,
                        backend_kwargs,
                        primer_to_panel_map,
                        dimer_conc_max_nm=dimer_conc_max,
                        boltzmann_max_dimers=boltzmann_max_dimers,
                        boltzmann_dedup_3p_window=boltzmann_dedup_3p_window,
                        boltzmann_stats_out=boltz_stats_sep,
                    )
                    sep_evaluated_names.update(boltz_stats_sep.get("evaluated_names", set()))
                    sep_skipped_names.update(boltz_stats_sep.get("skipped_names", set()))
                    sep_eligible_dkeys.update(
                        boltz_stats_sep.get("eligible_dedup_keys", set())
                    )
                    sep_evaluated_dkeys.update(
                        boltz_stats_sep.get("evaluated_dedup_keys", set())
                    )
                    sep_eligible_ppairs.update(
                        boltz_stats_sep.get("eligible_parent_pair_keys", set())
                    )
                    sep_evaluated_ppairs.update(
                        boltz_stats_sep.get("evaluated_parent_pair_keys", set())
                    )
                    for pp_lbl, pp_prob in boltz_stats_sep.get(
                        "max_prob_by_parent_pair", {}
                    ).items():
                        if float(pp_prob) > sep_max_prob_by_parent.get(pp_lbl, -1.0):
                            sep_max_prob_by_parent[pp_lbl] = float(pp_prob)

                    risks = evaluate_risks(
                        complexes,
                        concs,
                        amplifiable_flags,
                        is_warm_start=enzyme.warm_start,
                        dimer_details=dimer_details,
                    )
                    fractions = compute_fractions(
                        strands,
                        complexes,
                        prob.stoichiometry,
                        res.free_concentrations,
                        prob.delta_g,
                        temp_celsius,
                        primer_to_panel,
                    )

                    for r in risks:
                        if not any(ar.complex_name == r.complex_name for ar in all_risks):
                            all_risks.append(r)

                    target_occupations[t_id] = {}
                    for s in strands:
                        if s.endswith("_site"):
                            idx = strands.index(s)
                            tot_s = float(prob.total_concentrations[idx])
                            free_s = float(res.free_concentrations[idx])
                            occ = (tot_s - free_s) / tot_s if tot_s > 0 else 0.0
                            target_occupations[t_id][s] = occ

                    all_fractions.update(fractions)
                    if t_id not in all_unfolding_penalties:
                        all_unfolding_penalties[t_id] = unfolding_penalties
                boltzmann_eligible_dedup_groups = len(sep_eligible_dkeys)
                boltzmann_skipped_dedup_groups = len(sep_eligible_dkeys - sep_evaluated_dkeys)
                boltzmann_eligible_parent_pairs = len(sep_eligible_ppairs)
                boltzmann_skipped_parent_pairs = len(sep_eligible_ppairs - sep_evaluated_ppairs)
                boltzmann_skipped_eligible_dimers = len(sep_skipped_names - sep_evaluated_names)
                boltzmann_max_prob_by_parent_pair = dict(
                    sorted(sep_max_prob_by_parent.items(), key=lambda kv: (-kv[1], kv[0]))
                )

            # Matrice d'interaction globale pour la carte thermique / Global interaction matrix for heatmap
            unique_parents = list(dict.fromkeys(_get_parent_name(p.name) for p in primers))
            for parent1 in unique_parents:
                interaction_matrix[parent1] = {}
                for parent2 in unique_parents:
                    interaction_matrix[parent1][parent2] = 0.0

            if complexes and prob is not None:
                for i_c, c_name in enumerate(complexes):
                    if "_on_" in c_name or "_free" in c_name or c_name.endswith("_hairpin"):
                        continue
                    dg = float(prob.delta_g[i_c])
                    stoich = prob.stoichiometry[i_c]
                    p_indices = [
                        k for k, val in enumerate(stoich) if val > 0 and k < len(primers)
                    ]
                    if len(p_indices) == 1:
                        k = p_indices[0]
                        par = _get_parent_name(primers[k].name)
                        if interaction_matrix[par][par] == 0.0 or dg < interaction_matrix[par][par]:
                            interaction_matrix[par][par] = dg
                    elif len(p_indices) == 2:
                        k1, k2 = p_indices
                        par1 = _get_parent_name(primers[k1].name)
                        par2 = _get_parent_name(primers[k2].name)
                        if (
                            interaction_matrix[par1][par2] == 0.0
                            or dg < interaction_matrix[par1][par2]
                        ):
                            interaction_matrix[par1][par2] = dg
                            interaction_matrix[par2][par1] = dg

        # Verdict et décomposition des dimères / Verdict and dimer decomposition
        loop_primer_parents = {
            _get_parent_name(p.name) for p in primers if p.role in (PrimerRole.LF, PrimerRole.LB)
        }
        parent_nominal_concs: Dict[str, float] = {}
        parent_variant_counts: Dict[str, int] = {}
        for p in primers:
            par_key = p.parent_name or _get_parent_name(p.name)
            parent_nominal_concs[par_key] = parent_nominal_concs.get(par_key, 0.0) + float(
                p.nominal_concentration or 0.0
            )
            parent_variant_counts[par_key] = parent_variant_counts.get(par_key, 0) + 1
        verdict = generate_verdict(
            all_fractions,
            target_occupations,
            all_risks,
            loop_primer_parents=loop_primer_parents,
            primer_to_panel=primer_to_panel_map,
            primer_to_target=primer_to_target,
            unfolding_penalties=all_unfolding_penalties,
            dimer_conc_max_nm=dimer_conc_max,
            role_imbalance_ratio=role_imbalance_ratio,
            dimer_warn_ratio=dimer_warn_ratio,
            primer_nominal_concs=parent_nominal_concs,
            primer_variant_counts=parent_variant_counts,
        )

        dimer_decomposition = summarize_amplifiable_dimers(
            all_risks, primer_to_panel=primer_to_panel_map
        )
        if dimer_decomposition["total_nm"] > 0.0:
            dom_f = dimer_decomposition["dominant_family"]
            dom_txt = (
                f" | Famille dominante : {dom_f['label']} ({_format_nm(dom_f['conc_nm'])} nM)"
                if dom_f
                else ""
            )
            if dimer_decomposition["is_multiplex"]:
                intra_txt = ", ".join(
                    f"{pan}={_format_nm(val)} nM"
                    for pan, val in dimer_decomposition["intra_by_panel_nm"].items()
                )
                emit(
                    f"Dimères amplifiables : {_format_nm(dimer_decomposition['total_nm'])} nM "
                    f"(intra : {intra_txt} ; inter-panels : {_format_nm(dimer_decomposition['inter_nm'])} nM{dom_txt})"
                )
            else:
                emit(f"Dimères amplifiables : {_format_nm(dimer_decomposition['total_nm'])} nM{dom_txt}")

        # Recommandations, sondes TaqMan, concentrations et mésamorçage inter-cibles
        # Recommendations, TaqMan probes, concentrations, and inter-target mispriming
        recommendations = generate_recommendations(
            verdict, risks=all_risks, primer_to_panel=primer_to_panel_map
        )

        probe_tm_results: List[Dict[str, Any]] = []
        if any(p.role == PrimerRole.PROBE for p in primers):
            probe_tm_results = check_probes_tm(primers, backend, temp_celsius, **backend_kwargs)

        optimization_results: Union[Dict[str, Any], List[Any]] = []
        free_fractions = {p: f.free for p, f in all_fractions.items()}
        if prob is not None and complexes is not None and strands is not None:
            optimization_results = optimize_concentrations(
                prob_template=prob,
                species_names=strands,
                primers=primers,
                target_dict=targets_dict,
                primer_to_panel=primer_to_panel,
                original_free_fractions=free_fractions,
                original_target_occupations={},
                complex_names=complexes,
                temp_celsius=temp_celsius,
                backend=backend,
                enzyme=enzyme,
                **backend_kwargs,
            )

        panel_summaries: Dict[str, Dict[str, Any]] = {}
        balance_cv: Optional[float] = None
        if target_occupations:
            panel_summaries, balance_cv = calculate_multiplex_balance(
                primer_to_target,
                target_occupations,
                free_fractions,
                loop_primer_parents=loop_primer_parents,
                primer_to_target=primer_to_target,
            )

        mispriming_risks: List[MisprimingRisk] = []
        if has_true_target and targets_dict:
            mispriming_risks = detect_inter_target_mispriming(
                primers,
                primer_to_target,
                targets_dict,
                backend,
                enzyme,
                temp_celsius,
                primer_to_target=primer_to_target,
                legitimate_loci=all_unfolding_penalties,
                **backend_kwargs,
            )

    try:
        pkg_version = importlib.metadata.version("labcraft")
    except importlib.metadata.PackageNotFoundError:
        pkg_version = "0.0.1"

    timestamp_utc = time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())
    multiplex_mode_str = (
        "Conjoint (joint multi-target equilibrium)"
        if not separate_targets
        else "Séparé (cible par cible)"
    )
    all_warning_messages = list(target_warnings or []) + [
        str(w.message) for w in captured_warnings
    ]

    metadata: Dict[str, Any] = {
        "version": pkg_version,
        "timestamp": timestamp_utc,
        "file_hash": file_hash,
        "max_residual": max_residual_global,
        "temperature": temp_celsius,
        "buffer": (
            buf_dump if buf_dump is not None else "Reference conditions (1.0 M Na+, 0 mM Mg2+)"
        ),
        "enzyme": enzyme.name,
        "dimer_dg_threshold": enzyme.dimer_dg_threshold,
        "dimer_conc_max_nm": dimer_conc_max,
        "dimer_warn_floor_nm": dimer_warn_ratio * dimer_conc_max,
        "boltzmann_max_dimers": boltzmann_max_dimers,
        "boltzmann_dedup_3p_window": boltzmann_dedup_3p_window,
        "boltzmann_eligible_dedup_groups": boltzmann_eligible_dedup_groups,
        "boltzmann_skipped_dedup_groups": boltzmann_skipped_dedup_groups,
        "boltzmann_eligible_parent_pairs": boltzmann_eligible_parent_pairs,
        "boltzmann_skipped_parent_pairs": boltzmann_skipped_parent_pairs,
        "boltzmann_skipped_eligible_dimers": boltzmann_skipped_eligible_dimers,
        "boltzmann_max_prob_by_parent_pair": boltzmann_max_prob_by_parent_pair,
        "concentrations_fip_bip": "Per target (check config)",
        "concentrations_target": "Per target (check config)",
        "interaction_matrix": interaction_matrix,
        "primer_names": list(interaction_matrix.keys()),
        "unfolding_penalties": all_unfolding_penalties,
        "target_occupations": target_occupations,
        "has_true_target": has_true_target,
        "primer_to_panel": primer_to_panel,
        "panel_summaries": panel_summaries,
        "balance_cv": balance_cv,
        "mispriming_risks": mispriming_risks,
        "chemistry": current_config.experiment.chemistry,
        "optimization_results": optimization_results,
        "probe_tm_results": probe_tm_results,
        "recommendations": recommendations,
        "loop_primer_parents": list(loop_primer_parents),
        "primers_parent": primers,
        "has_lna": any(len(p.lna_positions) > 0 for p in primers),
        "dimer_decomposition": dimer_decomposition,
        "primer_concentrations_uM": effective_concs_uM,
        "multiplex_mode": multiplex_mode_str,
        "warnings": all_warning_messages,
    }

    boltzmann_mode_str = (
        "dedup_approx" if int(boltzmann_dedup_3p_window) > 0 else "exhaustive"
    )
    conditions_json: Dict[str, Any] = {
        "temperature_C": temp_celsius,
        "enzyme": enzyme.name,
        "buffer": buf_dump,
        "dimer_conc_max_nm": dimer_conc_max,
        "dimer_warn_ratio": dimer_warn_ratio,
        "role_imbalance_ratio": role_imbalance_ratio,
        "boltzmann_mode": boltzmann_mode_str,
        "boltzmann_max_dimers": boltzmann_max_dimers,
        "boltzmann_dedup_3p_window": boltzmann_dedup_3p_window,
        "dimer_dg_threshold": enzyme.dimer_dg_threshold,
        "max_variants": max_variants,
        "seeding": {
            "errors": errors,
            "strict_3prime": strict_3prime,
            "strict_3prime_tolerate": strict_3prime_tolerate,
        },
        "multiplex_mode": multiplex_mode_str,
        "target_copies_per_uL": {t.id: t.copies_per_uL for t in current_config.targets},
    }

    json_payload = build_analyze_json_payload(
        version=pkg_version,
        timestamp=timestamp_utc,
        inputs=inputs_fp or {"config": None, "primers": None, "targets": None},
        conditions=conditions_json,
        primer_concentrations_uM=effective_concs_uM,
        fractions=all_fractions,
        unfolding_penalties=all_unfolding_penalties,
        target_occupations=target_occupations,
        risks=all_risks,
        dimer_decomposition=dimer_decomposition,
        mispriming_risks=mispriming_risks,
        verdict=verdict,
        panel_summaries=panel_summaries,
        boltzmann_mode=boltzmann_mode_str,
        boltzmann_eligible_dedup_groups=boltzmann_eligible_dedup_groups,
        boltzmann_skipped_dedup_groups=boltzmann_skipped_dedup_groups,
        boltzmann_eligible_parent_pairs=boltzmann_eligible_parent_pairs,
        boltzmann_skipped_parent_pairs=boltzmann_skipped_parent_pairs,
        boltzmann_skipped_eligible_dimers=boltzmann_skipped_eligible_dimers,
        boltzmann_max_prob_by_parent_pair=boltzmann_max_prob_by_parent_pair,
    )

    return PanelAnalysisResult(
        config=current_config,
        panel_label=panel_label,
        primers=primers,
        primer_to_target=primer_to_target,
        primer_to_panel_map=primer_to_panel_map,
        effective_concs_uM=effective_concs_uM,
        temp_celsius=temp_celsius,
        enzyme=enzyme,
        has_true_target=has_true_target,
        fractions=all_fractions,
        target_occupations=target_occupations,
        unfolding_penalties=all_unfolding_penalties,
        risks=all_risks,
        dimer_decomposition=dimer_decomposition,
        interaction_matrix=interaction_matrix,
        verdict=verdict,
        recommendations=recommendations,
        probe_tm_results=probe_tm_results,
        optimization_results=optimization_results,
        panel_summaries=panel_summaries,
        balance_cv=balance_cv,
        mispriming_risks=mispriming_risks,
        max_residual=max_residual_global,
        warnings=all_warning_messages,
        metadata=metadata,
        json_payload=json_payload,
    )


def run_analysis_pipeline(
    config_obj: PanelConfig,
    *,
    explicit_config: bool = True,
    inputs_fp: Optional[Dict[str, Any]] = None,
    file_hash: str = "in-memory",
    max_variants: int = 64,
    errors: int = 2,
    strict_3prime: int = 3,
    strict_3prime_tolerate: int = 2,
    dimer_conc_max: float = 0.23,
    dimer_warn_ratio: float = DEFAULT_DIMER_WARN_RATIO,
    role_imbalance_ratio: float = DEFAULT_ROLE_IMBALANCE_RATIO,
    boltzmann_max_dimers: Optional[int] = DEFAULT_BOLTZMANN_MAX_DIMERS,
    boltzmann_dedup_3p_window: int = DEFAULT_BOLTZMANN_DEDUP_3P_WINDOW,
    separate_targets: bool = False,
    force_joint: bool = False,
    progress_callback: Optional[Callable[[str], None]] = None,
) -> List[PanelAnalysisResult]:
    """Pipeline complet d'analyse : chargement des cibles, découpage éventuel et résolution.
    Full analysis pipeline: target loading, collision splitting if needed, and equilibrium solving.
    """
    emit = progress_callback or (lambda _msg: None)
    with dna_session(temp_celsius=config_obj.experiment.temperature_C):
        targets_dict, target_warnings = load_targets_from_config(
            config_obj,
            explicit_config=explicit_config,
            progress_callback=progress_callback,
        )
        config_list, had_collisions = split_colliding_panel_configs(
            config_obj, force_joint=force_joint
        )
        if had_collisions:
            emit(
                "Attention: Plusieurs panels alternatifs ciblant la même séquence ont été détectés."
            )
            emit("Analyse séparée des panels...")

        results: List[PanelAnalysisResult] = []
        for i_cfg, current_config in enumerate(config_list):
            panel_label: Optional[str] = None
            if len(config_list) > 1:
                panel_label = current_config.primer_sets[0].panel_name or str(i_cfg)
                emit(f"--- Analyse du panel {panel_label} ---")

            res = run_single_panel_analysis(
                current_config,
                targets_dict,
                panel_label=panel_label,
                target_warnings=target_warnings,
                inputs_fp=inputs_fp,
                file_hash=file_hash,
                max_variants=max_variants,
                errors=errors,
                strict_3prime=strict_3prime,
                strict_3prime_tolerate=strict_3prime_tolerate,
                dimer_conc_max=dimer_conc_max,
                dimer_warn_ratio=dimer_warn_ratio,
                role_imbalance_ratio=role_imbalance_ratio,
                boltzmann_max_dimers=boltzmann_max_dimers,
                boltzmann_dedup_3p_window=boltzmann_dedup_3p_window,
                separate_targets=separate_targets,
                progress_callback=progress_callback,
            )
            results.append(res)
        return results
