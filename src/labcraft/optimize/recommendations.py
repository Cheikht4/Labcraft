"""Non-sequential recommendations.

Recommandations qualitatives (non-séquentielles) pour l'optimisation.
"""
from typing import Dict, List, Optional
from labcraft.metrics.verdict import PanelVerdict
from labcraft.metrics.risk import RiskItem, summarize_amplifiable_dimers


def _fmt_nm(val: float) -> str:
    """Formate une concentration en nM (4 décimales sous 0,01 nM).
    Formats a concentration in nM (4 decimals below 0.01 nM).
    """
    return f"{val:.4f}" if 0.0 < val < 0.01 else f"{val:.2f}"


def generate_recommendations(
    verdict: PanelVerdict,
    risks: Optional[List[RiskItem]] = None,
    primer_to_panel: Optional[Dict[str, str]] = None,
) -> List[str]:
    """
    Génère des recommandations qualitatives basées sur le verdict ET les risques réels.
    Generates qualitative recommendations based on verdict AND evaluated dimer risks.

    La détection des dimères amplifiables se fonde sur les RiskItem (severity >= 10.0)
    et leur décomposition par famille parente ; aucune recommandation sur les dimères
    n'est émise si le total est nul.
    """
    recs: List[str] = []

    decomp = summarize_amplifiable_dimers(risks or [], primer_to_panel=primer_to_panel) if risks else None
    has_amplifiable = bool(decomp and decomp["total_nm"] > 0.0)
    has_blocking = bool(risks and any(r.severity < 10.0 and r.concentration > 1e-9 for r in risks))

    # Fallback sur le verdict si les risques ne sont pas passés (rétrocompatibilité)
    # Fallback on verdict if risks are not passed (backward compatibility)
    if risks is None:
        has_amplifiable = any("amplifiable" in issue.cause.lower() for issue in verdict.issues)
        has_blocking = any("bloquant" in issue.cause.lower() for issue in verdict.issues)

    if has_amplifiable and decomp is not None:
        total_nm = decomp["total_nm"]
        dom = decomp["dominant_family"]
        if dom is not None:
            context_part = ""
            if decomp["is_multiplex"]:
                context_part = f" (dont {_fmt_nm(decomp['inter_nm'])} nM inter-panels)"
            recs.append(
                f"Famille de dimères amplifiables dominante : {dom['label']} "
                f"({_fmt_nm(dom['conc_nm'])} nM sur {_fmt_nm(total_nm)} nM au total{context_part}, "
                f"soit {dom['share_pct']:.0f} % du total)."
            )
            ext_p = dom.get("extensible_primer") or ""
            tpl_p = dom.get("template_primer") or ""
            if ext_p and tpl_p and ext_p != tpl_p:
                recs.append(
                    f"Extrémité 3' étendue : l'extrémité 3' de {ext_p} s'apparie et est étendue "
                    f"sur {tpl_p} servant de matrice. Déplacer ou modifier les bases 3' de {ext_p}, "
                    f"ou réduire sa concentration."
                )
            elif ext_p:
                recs.append(
                    f"Extrémité 3' étendue : l'extrémité 3' de {ext_p} est auto-complémentaire et étendue. "
                    f"Modifier les bases 3' de {ext_p} ou réduire sa concentration."
                )
            else:
                recs.append(
                    f"Extrémité 3' étendue au sein de la famille {dom['label']}. "
                    "Modifier l'extrémité 3' de l'amorce responsable ou ajuster les concentrations."
                )
    elif has_amplifiable:
        recs.append(
            "Des dimères amplifiables sont signalés : modifier l'extrémité 3' des amorces impliquées "
            "ou réduire leur concentration."
        )

    if has_blocking:
        recs.append(
            "Asymétrie de concentration : Certains dimères bloquants séquestrent massivement "
            "des amorces. Envisagez de baisser la concentration de l'amorce la plus abondante "
            "impliquée dans ces dimères, ou d'augmenter celle de l'amorce séquestrée."
        )

    if not recs:
        if verdict.status in ("WARNING", "FAILURE"):
            # Un panel en WARNING ou FAILURE ne doit jamais être qualifié de sain
            # A panel in WARNING or FAILURE must never be described as healthy
            recs.append(
                f"Diagnostic ({verdict.status}) : Aucune recommandation non-séquentielle sur les dimères, "
                f"mais le panel présente des limitations sur cible ({verdict.global_cause}). "
                "Consultez le diagnostic par cible et par amorce ci-dessus."
            )
        else:
            recs.append("Aucune recommandation non-séquentielle critique. Le panel semble sain.")

    return recs
