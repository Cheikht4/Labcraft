from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple


@dataclass
class RiskItem:
    complex_name: str
    concentration: float
    severity: float 
    description: str
    
    # Nouvelles métriques pour les fiches (Jalon 6.1)
    seq_a: str = ""
    seq_b: str = ""
    structure: str = ""
    delta_g: float = 0.0
    delta_g_3p: float = 0.0
    alignment_ascii: str = ""
    alignment_columns: list = field(default_factory=list)
    arrow_metrics: dict = field(default_factory=dict)
    is_blocked_veto: bool = False
    parent_a: str = ""
    parent_b: str = ""
    panel_a: str = ""
    panel_b: str = ""
    extensible_primer: str = ""
    template_primer: str = ""

    @property
    def is_amplifiable(self) -> bool:
        return self.severity >= 10.0

def evaluate_risks(
    complex_names: List[str],  
    concentrations: List[float], 
    amplifiable_flags: List[bool],
    is_warm_start: bool = False,
    dimer_details: Optional[List[dict]] = None
) -> List[RiskItem]:
    """
    Évalue le risque d'artefacts basé sur la concentration et l'extensibilité.
    """
    risks = []
    
    for i, (c_name, conc, is_amp) in enumerate(zip(complex_names, concentrations, amplifiable_flags, strict=True)):
        if "_free" in c_name or "_on_" in c_name:
            continue # Pas de risque pour les espèces libres ou liées à la cible
            
        severity = 0.0
        desc = ""
        
        details = dimer_details[i] if dimer_details and i < len(dimer_details) else {}

        if is_amp:
            # Sévérité pondérée par la marge d'amplifiabilité (ΔG 3' end) et la concentration
            margin = max(0.0, - (details.get("delta_g_3p", 0.0) + 3.0)) # ex: dg_3p=-5 -> margin = 2
            severity = 10.0 + margin
            desc = "Dimère amplifiable (3' extensible)"
            if is_warm_start:
                # La WarmStart ne change pas la thermodynamique à 65°C, 
                # mais elle évite l'amplification à froid. Le risque subsiste pendant l'amplification.
                desc += " - Note: La WarmStart mitige le risque de formation pré-incubation."
        else:
            # Dimère bloquant (appauvrit l'amorce mais ne s'amplifie pas)
            severity = 1.0 
            desc = "Dimère bloquant (non extensible)"
            if "homo" in c_name:
                desc = "Homodimère bloquant"
        
        if is_amp or conc > 1e-9:
            risks.append(RiskItem(
                complex_name=c_name,
                concentration=conc,
                severity=severity,
                description=desc,
                seq_a=details.get("seq_a", ""),
                seq_b=details.get("seq_b", ""),
                structure=details.get("structure", ""),
                delta_g=details.get("delta_g", 0.0),
                delta_g_3p=details.get("delta_g_3p", 0.0),
                alignment_ascii=details.get("alignment", ""),
                alignment_columns=details.get("alignment_columns", []),
                arrow_metrics=details.get("arrow_metrics", {}),
                is_blocked_veto=details.get("is_blocked_veto", False),
                parent_a=details.get("parent_a", ""),
                parent_b=details.get("parent_b", ""),
                panel_a=details.get("panel_a", ""),
                panel_b=details.get("panel_b", ""),
                extensible_primer=details.get("extensible_primer_parent", ""),
                template_primer=details.get("template_primer_parent", ""),
            ))
            
    # Trier par risque global (sévérité * concentration)
    risks.sort(key=lambda r: r.severity * r.concentration, reverse=True)
    return risks


def _infer_parent_pair(risk: RiskItem) -> Tuple[str, str]:
    """Extrait la paire d'amorces parentes d'un RiskItem (métadonnées ou nom du complexe).
    Extracts the parent primer pair from a RiskItem (explicit metadata or complex name fallback).
    """
    if risk.parent_a and risk.parent_b:
        return risk.parent_a, risk.parent_b
    cname = risk.complex_name
    if ":" in cname:
        left, right = cname.split(":", 1)
        return left.split("#")[0], right.split("#")[0]
    if cname.startswith("homo_"):
        p = cname[len("homo_"):].split("#")[0]
        return p, p
    return cname.split("#")[0], ""


def summarize_amplifiable_dimers(
    risks: List[RiskItem],
    primer_to_panel: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """Décompose les dimères amplifiables par panel (intra vs inter-panels) et par famille parente.
    Decomposes amplifiable dimers by panel (intra vs inter-panel) and by parent primer family.
    """
    amp_risks = [r for r in risks if r.is_amplifiable and r.concentration > 0.0]

    known_panels: List[str] = []
    if primer_to_panel:
        for pan in primer_to_panel.values():
            if pan and pan not in known_panels:
                known_panels.append(pan)

    intra_by_panel_nm: Dict[str, float] = {pan: 0.0 for pan in known_panels}
    inter_nm = 0.0
    total_nm = 0.0

    fam_buckets: Dict[Tuple[str, str], Dict[str, Any]] = {}

    for r in amp_risks:
        c_nm = float(r.concentration) * 1e9
        total_nm += c_nm
        pa, pb = _infer_parent_pair(r)

        pan_a = r.panel_a or (primer_to_panel.get(pa, "") if primer_to_panel else "")
        pan_b = r.panel_b or (primer_to_panel.get(pb, "") if primer_to_panel else "")
        for p in (pan_a, pan_b):
            if p and p not in known_panels:
                known_panels.append(p)
                intra_by_panel_nm.setdefault(p, 0.0)

        is_inter = bool(pan_a and pan_b and pan_a != pan_b)
        if is_inter:
            inter_nm += c_nm
        elif pan_a:
            intra_by_panel_nm[pan_a] = intra_by_panel_nm.get(pan_a, 0.0) + c_nm

        key = tuple(sorted((pa, pb)))
        if key not in fam_buckets:
            label = pa if not pb else f"{pa} × {pb}"
            fam_buckets[key] = {
                "label": label,
                "parent_a": pa,
                "parent_b": pb,
                "conc_nm": 0.0,
                "is_inter_panel": is_inter,
                "panel_a": pan_a,
                "panel_b": pan_b,
                "_ext_weights": {},
            }
        bucket = fam_buckets[key]
        bucket["conc_nm"] += c_nm
        if r.extensible_primer:
            ext_key = (r.extensible_primer, r.template_primer)
            bucket["_ext_weights"][ext_key] = bucket["_ext_weights"].get(ext_key, 0.0) + c_nm

    families: List[Dict[str, Any]] = []
    for bucket in sorted(fam_buckets.values(), key=lambda x: x["conc_nm"], reverse=True):
        ext_weights = bucket.pop("_ext_weights")
        ext_primer, tpl_primer = ("", "")
        if ext_weights:
            ext_primer, tpl_primer = max(ext_weights.items(), key=lambda kv: kv[1])[0]
        bucket["extensible_primer"] = ext_primer
        bucket["template_primer"] = tpl_primer
        bucket["share_pct"] = (bucket["conc_nm"] / total_nm * 100.0) if total_nm > 0.0 else 0.0
        families.append(bucket)

    return {
        "total_nm": total_nm,
        "is_multiplex": len(known_panels) > 1,
        "panels": known_panels,
        "intra_by_panel_nm": intra_by_panel_nm,
        "inter_nm": inter_nm,
        "families": families,
        "dominant_family": families[0] if families else None,
    }
