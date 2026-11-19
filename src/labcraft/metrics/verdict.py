"""Panel verdict generation / Génération du verdict de panel.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple
import numpy as np
from labcraft.metrics.fractions import PrimerFractions
from labcraft.metrics.risk import RiskItem

# Seuils heuristiques non calibrés gouvernant la frontière INFO / WARNING :
# - `DEFAULT_ROLE_IMBALANCE_RATIO = 20.0` : un rapport d'occupation brute >= 20× entre amorces
#   d'initiation jumelles (F3 vs B3, FIP vs BIP) correspond à un déséquilibre thermodynamique
#   ou de dilution ΔΔG = RT ln(20) ≈ 2,0 kcal/mol à 63 °C (RT ≈ 0,668 kcal/mol).
# - `DEFAULT_DIMER_WARN_RATIO = 0.5` : WARNING lorsque les dimères amplifiables atteignent 50 %
#   du plafond critique `dimer_conc_max_nm` (et INFO en deçà).
# Uncalibrated heuristic thresholds governing the INFO / WARNING boundary:
# - `DEFAULT_ROLE_IMBALANCE_RATIO = 20.0`: a raw occupation ratio >= 20× between twin initiation
#   primers (F3 vs B3, FIP vs BIP) corresponds to ΔΔG = RT ln(20) ≈ 2.0 kcal/mol at 63 °C.
# - `DEFAULT_DIMER_WARN_RATIO = 0.5`: WARNING when amplifiable dimers reach 50% of `dimer_conc_max_nm`.
DEFAULT_ROLE_IMBALANCE_RATIO: float = 20.0
DEFAULT_DIMER_WARN_RATIO: float = 0.5

_TWIN_ROLE_MAP: Dict[str, str] = {
    "F3": "B3",
    "B3": "F3",
    "FIP": "BIP",
    "BIP": "FIP",
}


def _extract_role_and_group(parent_name: str, panel_map: Optional[Dict[str, str]], target_id: str) -> Tuple[str, str]:
    """Infère (rôle_normalisé, groupe_de_panel) pour apparier F3<->B3 et FIP<->BIP d'un même panel.
    Infers (normalized_role, panel_group) to pair F3<->B3 and FIP<->BIP within the same panel.
    """
    panel_id = (panel_map.get(parent_name) if panel_map else None) or ""
    parts = parent_name.split("_")
    first_up = parts[0].upper()
    last_up = parts[-1].upper()
    if first_up in _TWIN_ROLE_MAP:
        grp = panel_id or ("_".join(parts[1:]) if len(parts) > 1 else target_id)
        return first_up, grp
    if last_up in _TWIN_ROLE_MAP:
        grp = panel_id or ("_".join(parts[:-1]) if len(parts) > 1 else target_id)
        return last_up, grp
    return "", panel_id or target_id


def _fmt_fr_1dec(val: float) -> str:
    """Formate un flottant à 1 décimale avec virgule française / Formats float to 1 decimal with French comma."""
    return f"{val:.1f}".replace(".", ",")


def _fmt_nm(val: float) -> str:
    """Formate une concentration en nM / Formats a concentration in nM."""
    return f"{val:.4f}" if 0.0 < val < 0.01 else f"{val:.2f}"


@dataclass
class PrimerIssue:
    primer_name: str
    target_site: str
    occupation: float
    cause: str
    is_critical: bool
    level: str = "WARNING"  # "INFO", "WARNING", "FAILURE"

    def __post_init__(self) -> None:
        if self.is_critical:
            self.level = "FAILURE"


@dataclass
class PanelVerdict:
    status: str  # "OK", "INFO", "WARNING", "FAILURE"
    issues: List[PrimerIssue]
    global_cause: str


def generate_verdict(
    fractions: Dict[str, PrimerFractions],
    target_occupations: Dict[str, Dict[str, float]],
    risks: List[RiskItem],
    loop_primer_parents: Optional[set[str]] = None,
    primer_to_panel: Optional[Dict[str, str]] = None,
    primer_to_target: Optional[Dict[str, str]] = None,
    unfolding_penalties: Optional[Dict[str, Any]] = None,
    dimer_conc_max_nm: Optional[float] = None,
    *,
    role_imbalance_ratio: float = DEFAULT_ROLE_IMBALANCE_RATIO,
    dimer_warn_ratio: float = DEFAULT_DIMER_WARN_RATIO,
    primer_nominal_concs: Optional[Dict[str, float]] = None,
    primer_variant_counts: Optional[Dict[str, int]] = None,
) -> PanelVerdict:
    """Évalue le panel selon quatre niveaux (`OK`, `INFO`, `WARNING`, `FAILURE`) :
    - Une occupation d'initiation `< 10 %` non nulle est classée `INFO` par défaut (comportement
      attendu sur génome complet), et ne passe en `WARNING` que si l'amorce subit une séquestration
      notable (`mean_free < 0.20`) ou un déséquilibre brut `>= role_imbalance_ratio` (20× par défaut,
      sans normalisation par `N_var`) face à l'amorce jumelle de même rôle (`F3` vs `B3`, `FIP` vs `BIP`),
      en distinguant explicitement la dilution par dégénérescence d'un blocage thermodynamique du site.
    - Les dimères amplifiables sous `dimer_conc_max_nm` sont classés `WARNING` si
      `total_amp_dimer_conc_nm >= dimer_warn_ratio * dimer_conc_max_nm`, et `INFO` sinon.

    Evaluates the panel across four levels (`OK`, `INFO`, `WARNING`, `FAILURE`):
    - Non-zero initiation occupation `< 10%` defaults to `INFO` unless free fraction `< 0.20` or raw
      occupation lags its twin role (`F3` vs `B3`, `FIP` vs `BIP`) by `>= role_imbalance_ratio` (default 20×,
      without `N_var` normalization), explicitly distinguishing degeneracy dilution from site blockage.
    - Sub-threshold amplifiable dimers trigger `WARNING` at `>= dimer_warn_ratio * dimer_conc_max_nm`, else `INFO`.
    """
    if loop_primer_parents is None:
        loop_primer_parents = set()

    assoc_map = primer_to_target if primer_to_target is not None else primer_to_panel
    issues: List[PrimerIssue] = []

    if target_occupations:
        for target_id, occupations in target_occupations.items():
            parent_to_variants: Dict[str, List[str]] = {}
            for primer_name in fractions.keys():
                parent_name = primer_name.split("#")[0] if "#" in primer_name else primer_name
                parent_to_variants.setdefault(parent_name, []).append(primer_name)

            # Première passe : collecter l'occupation réelle θ de chaque amorce parente sur cette cible
            # First pass: collect actual site occupation θ of each parent primer on this target
            target_parents: List[Tuple[str, List[str], float, bool]] = []
            role_occ_lookup: Dict[Tuple[str, str], Tuple[str, float]] = {}

            for parent_name, var_names in parent_to_variants.items():
                if assoc_map is not None:
                    p_target = assoc_map.get(parent_name) or assoc_map.get(var_names[0])
                    if p_target and p_target != target_id:
                        continue
                else:
                    site_name_check = f"{parent_name}_site"
                    has_site_in_this_target = (
                        site_name_check in occupations
                        or any(
                            k.startswith(f"{parent_name}@") or k.startswith(f"{parent_name}_")
                            for k in occupations
                        )
                    )
                    if not has_site_in_this_target and len(target_occupations) > 1:
                        continue

                site_name = f"{parent_name}_site"
                occ = occupations.get(site_name)
                if occ is None:
                    for k, v in occupations.items():
                        if k == site_name or k.startswith(f"{parent_name}@") or k.startswith(f"{parent_name}_"):
                            occ = v
                            break
                has_site = (
                    site_name in occupations
                    or any(
                        k.startswith(f"{parent_name}@") or k.startswith(f"{parent_name}_")
                        for k in occupations
                    )
                )
                occ_val = float(occ) if occ is not None else 0.0
                target_parents.append((parent_name, var_names, occ_val, has_site))

                if parent_name not in loop_primer_parents:
                    role_key, grp_key = _extract_role_and_group(parent_name, primer_to_panel, target_id)
                    if role_key:
                        prev = role_occ_lookup.get((grp_key, role_key))
                        if prev is None or occ_val > prev[1]:
                            role_occ_lookup[(grp_key, role_key)] = (parent_name, occ_val)

            # Seconde passe : diagnostic causal et niveau (`INFO` / `WARNING` / `FAILURE`)
            # Second pass: causal attribution and severity level (`INFO` / `WARNING` / `FAILURE`)
            for parent_name, var_names, occ_val, has_site in target_parents:
                if occ_val >= 0.1 or parent_name in loop_primer_parents:
                    continue

                f_list = [fractions[v] for v in var_names if v in fractions]
                mean_free = float(np.mean([f.free for f in f_list])) if f_list else 1.0
                dom_complex = f_list[0].dominant_complex if f_list else ""
                dom_frac_pct = float(np.mean([f.dominant_fraction for f in f_list])) * 100.0 if f_list else 0.0
                mean_inter = float(np.mean([f.heterodimer_inter for f in f_list])) if f_list else 0.0
                mean_intra = float(np.mean([f.heterodimer_intra for f in f_list])) if f_list else 0.0
                mean_homo = float(np.mean([f.homodimer for f in f_list])) if f_list else 0.0

                is_crit = (occ_val == 0.0) or (mean_free < 0.05)

                site_name = f"{parent_name}_site"
                pen_for_target = (
                    unfolding_penalties.get(target_id, {}) if isinstance(unfolding_penalties, dict) else {}
                )
                site_pen = pen_for_target.get(site_name) or pen_for_target.get(f"{parent_name}@{target_id}_site")

                role_key, grp_key = _extract_role_and_group(parent_name, primer_to_panel, target_id)
                twin_role = _TWIN_ROLE_MAP.get(role_key, "")
                twin_entry = role_occ_lookup.get((grp_key, twin_role)) if twin_role else None
                twin_name = twin_entry[0] if twin_entry is not None else ""
                twin_occ = twin_entry[1] if twin_entry is not None else 0.0
                imbalance = (twin_occ / occ_val) if (occ_val > 0.0 and twin_occ > 0.0) else 0.0
                has_role_imbalance = (
                    role_imbalance_ratio > 0.0 and imbalance >= role_imbalance_ratio
                )

                n_vars = (
                    int(primer_variant_counts[parent_name])
                    if (primer_variant_counts and parent_name in primer_variant_counts)
                    else len(var_names)
                )
                if is_crit:
                    issue_level = "FAILURE"
                elif mean_free < 0.20 or has_role_imbalance:
                    issue_level = "WARNING"
                else:
                    issue_level = "INFO"

                if isinstance(site_pen, dict) and site_pen.get("verdict") == "veto_3p":
                    bad_pos = site_pen.get("first_bad_pos", 1)
                    cause = f"Veto d'extension en 3' pour {parent_name} (mésappariement bloquant à pos -{bad_pos} du 3')."
                elif isinstance(site_pen, dict) and site_pen.get("verdict") == "indetermine":
                    cause = f"Site indéterminé pour {parent_name} (présence de base(s) N dans le domaine de liaison)."
                elif mean_free < 0.20:
                    if mean_inter > mean_homo and mean_inter > mean_intra:
                        cause = (
                            f"Hybridation croisée inter-jeux séquestrant {parent_name} "
                            f"à {dom_frac_pct:.1f}% dans {dom_complex}."
                        )
                    else:
                        cause = f"Séquestration de {parent_name} à {dom_frac_pct:.1f}% dans {dom_complex}."
                elif not has_site or occ_val == 0.0:
                    cause = "Site absent de la cible ou amorce non appariée."
                elif has_role_imbalance:
                    if n_vars > 1:
                        nom_molar = primer_nominal_concs.get(parent_name) if primer_nominal_concs else None
                        if nom_molar is not None and nom_molar > 0.0:
                            per_var_nm = (nom_molar / n_vars) * 1e9
                            var_detail = f"{parent_name} : {n_vars} variants, ~{_fmt_fr_1dec(per_var_nm)} nM par variant"
                        else:
                            var_detail = f"{parent_name} : {n_vars} variants"
                        cause = (
                            f"Occupation réduite par dégénérescence ({var_detail} ; "
                            f"déséquilibre {_fmt_fr_1dec(imbalance)}× face à {twin_name})."
                        )
                    else:
                        cause = (
                            f"Blocage thermodynamique relatif du site cible pour {parent_name} "
                            f"(déséquilibre {_fmt_fr_1dec(imbalance)}× face à {twin_name} : "
                            "barrière d'ouverture structurale ou affinité défavorable)."
                        )
                else:
                    if n_vars > 1:
                        cause = (
                            f"Inaccessibilité du site cible ou dilution dégénérée "
                            f"({n_vars} variantes au panel, conc. utile fractionnée)."
                        )
                    else:
                        cause = (
                            f"Inaccessibilité du site cible pour {parent_name} "
                            "(barrière d'ouverture structurale ou thermodynamique défavorable)."
                        )

                issues.append(
                    PrimerIssue(
                        primer_name=parent_name,
                        target_site=target_id,
                        occupation=occ_val,
                        cause=cause,
                        is_critical=is_crit,
                        level=issue_level,
                    )
                )
    else:
        # Mode compétition pure sans cible / Target-free competition mode
        for primer_name, f in fractions.items():
            if f.free < 0.20:
                dom_frac_pct = f.dominant_fraction * 100.0
                is_crit = f.free < 0.05
                if f.heterodimer_inter > f.homodimer and f.heterodimer_inter > f.heterodimer_intra:
                    cause = (
                        f"Hybridation croisée inter-jeux séquestrant {primer_name} "
                        f"à {dom_frac_pct:.1f}% dans {f.dominant_complex}."
                    )
                else:
                    cause = f"Séquestration de {primer_name} à {dom_frac_pct:.1f}% dans {f.dominant_complex}."

                issues.append(
                    PrimerIssue(
                        primer_name=primer_name,
                        target_site="Compétition Sans Cible",
                        occupation=f.free,
                        cause=cause,
                        is_critical=is_crit,
                        level="FAILURE" if is_crit else "WARNING",
                    )
                )

    # Évaluation globale à 4 niveaux : OK / INFO / WARNING / FAILURE
    # 4-level global evaluation: OK / INFO / WARNING / FAILURE
    total_amp_dimer_conc_nm = (
        sum(r.concentration for r in risks if r.is_amplifiable and r.severity >= 10.0) * 1e9
    )
    has_any_amp_dimers = any(r.is_amplifiable and r.severity >= 10.0 and r.concentration > 0.0 for r in risks)

    if dimer_conc_max_nm is not None:
        has_critical_dimers = total_amp_dimer_conc_nm >= dimer_conc_max_nm
        warn_floor_nm = dimer_warn_ratio * dimer_conc_max_nm
        has_warn_dimers = (not has_critical_dimers) and (total_amp_dimer_conc_nm >= warn_floor_nm)
        has_info_dimers = (not has_critical_dimers) and (not has_warn_dimers) and has_any_amp_dimers
    else:
        warn_floor_nm = 0.0
        has_critical_dimers = any(r.is_amplifiable and r.severity >= 10.0 for r in risks)
        has_warn_dimers = False
        has_info_dimers = False

    has_critical_issues = any(i.is_critical or i.level == "FAILURE" for i in issues)
    has_warn_issues = any((not i.is_critical) and i.level == "WARNING" for i in issues)
    has_info_issues = any((not i.is_critical) and i.level == "INFO" for i in issues)

    if has_critical_dimers or has_critical_issues:
        status = "FAILURE"
        global_cause = (
            "Présence de dimères d'amorces amplifiables majeurs ou séquestration critique "
            "interdisant l'amplification."
        )
    elif has_warn_dimers or has_warn_issues or any(5.0 <= r.severity < 10.0 for r in risks):
        status = "WARNING"
        deg_warn_causes = [
            i.cause for i in issues if i.level == "WARNING" and "dégénérescence" in i.cause
        ]
        if has_warn_dimers and dimer_conc_max_nm is not None:
            global_cause = (
                f"Dimères amplifiables modérés ({_fmt_nm(total_amp_dimer_conc_nm)} nM >= plancher WARNING "
                f"{_fmt_nm(warn_floor_nm)} nM, < seuil FAILURE {dimer_conc_max_nm:.2f} nM)."
            )
        elif deg_warn_causes:
            global_cause = " ".join(deg_warn_causes)
        else:
            global_cause = "Déséquilibre d'occupation entre amorces de même rôle ou compétitions modérées identifiées."
    elif has_info_dimers or has_info_issues:
        status = "INFO"
        if has_info_dimers and dimer_conc_max_nm is not None:
            global_cause = (
                f"Traces mineures de dimères amplifiables ({_fmt_nm(total_amp_dimer_conc_nm)} nM < plancher WARNING "
                f"{_fmt_nm(warn_floor_nm)} nM) et occupations d'initiation équilibrées."
            )
        else:
            global_cause = (
                "Occupations d'initiation basses mais équilibrées entre amorces de même rôle "
                "(comportement attendu sur génome complet), sans dimère amplifiable significatif."
            )
    else:
        status = "OK"
        global_cause = "Aucun risque majeur d'interaction parasite identifié."

    return PanelVerdict(status=status, issues=issues, global_cause=global_cause)
