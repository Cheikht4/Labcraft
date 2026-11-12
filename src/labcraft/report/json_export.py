"""Export JSON pour la commande `labcraft analyze`.
JSON export for the `labcraft analyze` command.
"""
from __future__ import annotations

import json
import os
from typing import Any, Dict, List, Optional

from labcraft.metrics.fractions import PrimerFractions
from labcraft.metrics.risk import RiskItem
from labcraft.metrics.verdict import PanelVerdict


def json_path_for_output(output_path: str) -> str:
    """Calcule le chemin JSON associé à un fichier de sortie HTML.
    Computes the JSON sidecar path for an HTML output path.
    """
    if output_path.lower().endswith(".html"):
        return output_path[:-5] + ".json"
    base, _ = os.path.splitext(output_path)
    return base + ".json"


def build_analyze_json_payload(
    *,
    version: str,
    timestamp: str,
    inputs: Dict[str, Any],
    conditions: Dict[str, Any],
    primer_concentrations_uM: Dict[str, float],
    fractions: Dict[str, PrimerFractions],
    unfolding_penalties: Dict[str, Any],
    target_occupations: Dict[str, Dict[str, float]],
    risks: List[RiskItem],
    dimer_decomposition: Dict[str, Any],
    mispriming_risks: List[Any],
    verdict: PanelVerdict,
    panel_summaries: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Construit le dictionnaire JSON structuré de `labcraft analyze`.
    Builds the structured JSON payload for `labcraft analyze`.
    """
    free_fractions = {p: round(float(f.free), 6) for p, f in fractions.items()}
    fractions_detail = {
        p: {
            "free": float(f.free),
            "homodimer": float(f.homodimer),
            "heterodimer_intra": float(f.heterodimer_intra),
            "heterodimer_inter": float(f.heterodimer_inter),
        }
        for p, f in fractions.items()
    }

    # Sites par cible / Sites by target
    sites_by_target: Dict[str, Dict[str, Any]] = {}
    for t_id, sites_map in (unfolding_penalties or {}).items():
        if not isinstance(sites_map, dict):
            continue
        occ_map = target_occupations.get(t_id, {})
        t_sites: Dict[str, Any] = {}
        for site_key, pen in sites_map.items():
            if not isinstance(pen, dict):
                continue
            occ_val = occ_map.get(site_key)
            if occ_val is None and "@" in site_key:
                occ_val = occ_map.get(f"{site_key.split('@', 1)[0]}_site")
            t_sites[site_key] = {
                "position": pen.get("position"),
                "start": pen.get("start"),
                "end": pen.get("end"),
                "strand": pen.get("strand"),
                "mismatches": pen.get("mismatches"),
                "dg_hyb": pen.get("dg_hyb"),
                "dg_open": pen.get("dg_unfold"),
                "occupation": occ_val,
                "extensible": pen.get("extensible"),
            }
        sites_by_target[t_id] = t_sites

    # Liste complète des dimères / Full list of dimers
    dimers_list: List[Dict[str, Any]] = []
    blocking_total_nm = 0.0
    for r in risks:
        conc_nm = float(r.concentration) * 1e9
        is_amp = bool(r.is_amplifiable)
        if not is_amp:
            blocking_total_nm += conc_nm
        is_inter = bool(r.panel_a and r.panel_b and r.panel_a != r.panel_b)
        dimers_list.append({
            "name": r.complex_name,
            "conc_nm": conc_nm,
            "dg_global": float(r.delta_g),
            "dg_3prime": float(r.delta_g_3p),
            "amplifiable": is_amp,
            "scope": "inter" if is_inter else "intra",
            "panel_a": r.panel_a,
            "panel_b": r.panel_b,
            "parent_a": r.parent_a,
            "parent_b": r.parent_b,
            "extensible_primer": r.extensible_primer,
            "template_primer": r.template_primer,
        })

    # Alertes de mésamorçage / Mispriming alerts
    mispriming_list: List[Dict[str, Any]] = []
    for mr in mispriming_risks or []:
        if hasattr(mr, "model_dump"):
            mispriming_list.append(mr.model_dump())
        elif isinstance(mr, dict):
            mispriming_list.append(mr)

    # Totaux alignés sur le rapport HTML / Totals aligned with HTML report
    totals = {
        "amplifiable_dimers_total_nm": float(dimer_decomposition.get("total_nm", 0.0)),
        "amplifiable_dimers_intra_by_panel_nm": dict(dimer_decomposition.get("intra_by_panel_nm", {})),
        "amplifiable_dimers_inter_nm": float(dimer_decomposition.get("inter_nm", 0.0)),
        "blocking_dimers_total_nm": blocking_total_nm,
        "is_multiplex": bool(dimer_decomposition.get("is_multiplex", False)),
        "dominant_family": dimer_decomposition.get("dominant_family"),
        "families": dimer_decomposition.get("families", []),
    }

    return {
        "version": version,
        "timestamp": timestamp,
        "inputs": inputs,
        "conditions": conditions,
        "primer_concentrations_uM": primer_concentrations_uM,
        "free_fractions": free_fractions,
        "fractions": fractions_detail,
        "sites_by_target": sites_by_target,
        "dimers": dimers_list,
        "totals": totals,
        "mispriming_risks": mispriming_list,
        "initiation_access_by_panel": panel_summaries or {},
        "verdict": verdict.status,
        "verdict_cause": verdict.global_cause,
        "issues": [
            {
                "primer_name": iss.primer_name,
                "target_site": iss.target_site,
                "occupation": float(iss.occupation),
                "cause": iss.cause,
                "is_critical": bool(iss.is_critical),
            }
            for iss in verdict.issues
        ],
    }


def _json_default(obj: Any) -> Any:
    """Convertit les types non natifs (numpy, enums, ensembles) pour json.dump.
    Converts non-native types (numpy, enums, sets) for json.dump.
    """
    if hasattr(obj, "item"):
        return obj.item()
    if hasattr(obj, "value"):
        return obj.value
    if isinstance(obj, (set, tuple)):
        return list(obj)
    return str(obj)


def write_analyze_json(output_html_path: str, payload: Dict[str, Any]) -> str:
    """Écrit le fichier JSON à côté du rapport HTML et retourne son chemin.
    Writes the JSON file alongside the HTML report and returns its path.
    """
    jpath = json_path_for_output(output_html_path)
    with open(jpath, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False, default=_json_default)
    return jpath
