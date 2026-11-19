"""Banc d'essai reproductible de référence LabCraft (`validation/reference_benchmarks.md`).
Reproducible LabCraft reference benchmark suite (`validation/reference_benchmarks.md`).

Produit à partir du code actuel :
(a) le tableau des 5 cas de paillasse de référence (dimères amplifiables à -2.0 et -3.5 kcal/mol,
    verdict, cause globale, indicateur `suboptimal_3p_masked` et `boltzmann_skipped_eligible_dimers`) ;
(b) le comparatif exact `labcraft coverage` entre les commits `58ae5d9` et `8cf5bd8` sur les 3 génomes
    GenBank de référence et les 2 cohortes synthétiques déterministes (120 DENV3 + 120 DENV4, graine 20261119) ;
(c) le benchmark de temps de calcul de `enrich_dimer_boltzmann_3p` sur 81 espèces (3 321 paires)
    sans plafond vs avec le plafond `--boltzmann-max-dimers 25`.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from typing import Any, Dict, List, Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))
if str(REPO_ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "scripts"))

from generate_coverage_cohorts import DEFAULT_SEED, generate_default_cohorts  # noqa: E402

from labcraft.cli.parsers import build_config_from_cli, read_multi_fasta, reverse_complement  # noqa: E402
from labcraft.diagnostics import amplifiable_dimer as amp_mod  # noqa: E402
from labcraft.diagnostics import enzyme as enz_mod  # noqa: E402
from labcraft.diagnostics.amplifiable_dimer import is_amplifiable_dimer  # noqa: E402
from labcraft.diagnostics.enzyme import (  # noqa: E402
    BST,
    BST_2_0,
    BST_2_0_WARM_START,
    BST_3_0,
    GSPSSD,
    TAQ,
)
from labcraft.lamp.domains import _find_iupac_substring  # noqa: E402
from labcraft.metrics.verdict import (  # noqa: E402
    DEFAULT_DIMER_WARN_RATIO,
    DEFAULT_ROLE_IMBALANCE_RATIO,
)
from labcraft.pipeline.analysis import (  # noqa: E402
    DEFAULT_BOLTZMANN_MAX_DIMERS,
    run_analysis_pipeline,
)

GT_DIR = REPO_ROOT / "tests" / "fixtures" / "ground_truth"
REAL_DIR = REPO_ROOT / "tests" / "fixtures" / "real_primer_files"


def _git_output(args: List[str]) -> str:
    """Exécute une commande git en lecture seule / Runs a read-only git command."""
    res = subprocess.run(
        ["git", "-C", str(REPO_ROOT), *args],
        check=True,
        capture_output=True,
        text=True,
    )
    return res.stdout.strip()


def _extract_commit_tree(commit_ref: str, dest_dir: Path) -> None:
    """Extrait l'arbre d'un commit via `git archive` / Extracts a commit tree via `git archive`."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    proc_arc = subprocess.Popen(
        ["git", "-C", str(REPO_ROOT), "archive", commit_ref],
        stdout=subprocess.PIPE,
    )
    assert proc_arc.stdout is not None
    subprocess.run(["tar", "-x", "-C", str(dest_dir)], stdin=proc_arc.stdout, check=True)
    proc_arc.wait()
    if proc_arc.returncode != 0:
        raise RuntimeError(f"git archive {commit_ref} failed with code {proc_arc.returncode}")


def _slice_around_amplicon(target_seq: str, primer_file: Path, pad: int = 100) -> str:
    """Découpe une fenêtre locale autour de l'amplicon / Slices a local window around the amplicon."""
    recs = dict(read_multi_fasta(str(primer_file)))
    spans: List[Tuple[int, int]] = []
    for name, seq in recs.items():
        if name.endswith(("FIP", "BIP")):
            continue
        for cand in (seq, reverse_complement(seq)):
            pos = _find_iupac_substring(cand, target_seq)
            if pos >= 0:
                spans.append((pos, pos + len(cand)))
    lo = max(0, min(a for a, _ in spans) - pad)
    hi = min(len(target_seq), max(b for _, b in spans) + pad)
    return target_seq[lo:hi]


def run_section_dg3p_shifts() -> List[Tuple[str, float, float, float, float]]:
    """Calcule le décalage exact de `local_dg` entre (37 °C, 1 M Na+) et (63 °C, tampon LAMP par défaut).
    Computes exact `local_dg` shift from (37 °C, 1 M Na+) to (63 °C, default LAMP buffer).
    """
    motifs = [
        ("3 pb (`5'-GGA-3'`)", "TTTTTTTGGA", "AAAATCCTAA", ".......(((...)))...."),
        ("4 pb (`5'-GCAT-3'`)", "TTTTTTGCAT", "AAAAATGCAA", "......((((..))))...."),
        ("5 pb (`5'-GCATG-3'`)", "TTTTTGCATG", "AAAACATGCA", ".....(((((.)))))...."),
        ("6 pb (`5'-GCATGC-3'`)", "TTTTGCATGC", "AAAAGCATGC", "....((((((....))))))"),
    ]
    rows: List[Tuple[str, float, float, float, float]] = []
    for label, sa, sb, st in motifs:
        _, dg_37_1m, _, _ = is_amplifiable_dimer(
            sa, sb, st, -5.0, BST_2_0, temp_celsius=37.0, na_mm=1000.0
        )
        _, dg_63_1m, _, _ = is_amplifiable_dimer(
            sa, sb, st, -5.0, BST_2_0, temp_celsius=63.0, na_mm=1000.0
        )
        _, dg_63_lamp, _, _ = is_amplifiable_dimer(
            sa,
            sb,
            st,
            -5.0,
            BST_2_0,
            temp_celsius=63.0,
            na_mm=50.0,
            k_mm=10.0,
            tris_mm=20.0,
            mg_mm=6.0,
            dntp_mm=1.4,
        )
        rows.append((label, dg_37_1m, dg_63_1m, dg_63_lamp, dg_63_lamp - dg_37_1m))
    return rows


def run_section_a_five_cases() -> Tuple[List[Dict[str, Any]], Dict[str, List[Dict[str, Any]]]]:
    """Exécute `run_analysis_pipeline` sur les 5 cas de référence à -2.0 et -3.5 kcal/mol.
    Runs `run_analysis_pipeline` on the 5 reference cases at -2.0 and -3.5 kcal/mol.
    """
    cases = [
        (
            "singleplex_denv3",
            "SUCCÈS paillasse",
            GT_DIR / "primer_dengue_3_panel3.txt",
            GT_DIR / "DEN3_M93130.fasta",
            0.8,
            0.4,
            0.2,
        ),
        (
            "singleplex_denv4",
            "SUCCÈS paillasse",
            GT_DIR / "parida_DENV4.txt",
            GT_DIR / "Dengue_4_S1028.fasta",
            0.8,
            0.4,
            0.2,
        ),
        (
            "singleplex_yfv",
            "SUCCÈS paillasse",
            GT_DIR / "nunes_YFV.txt",
            GT_DIR / "YFV_17D_17DD.fasta",
            0.8,
            0.4,
            0.2,
        ),
        (
            "duplex_half",
            "ÉCHEC paillasse (faux +)",
            GT_DIR / "duplex_YFV_DENV4_primers.txt",
            GT_DIR / "duplex_YFV_DENV4_targets.fasta",
            0.4,
            0.2,
            0.1,
        ),
        (
            "duplex_full",
            "ÉCHEC paillasse (faux +)",
            GT_DIR / "duplex_YFV_DENV4_primers.txt",
            GT_DIR / "duplex_YFV_DENV4_targets.fasta",
            0.8,
            0.4,
            0.2,
        ),
    ]

    rows: List[Dict[str, Any]] = []
    masked_by_case: Dict[str, List[Dict[str, Any]]] = {}
    orig_bst2 = enz_mod.BST_2_0
    try:
        for case_id, bench_label, p_path, t_path, c_fb, c_lp, c_f3 in cases:
            res_by_thresh: Dict[float, Any] = {}
            for thresh in (-2.0, -3.5):
                mod_p = replace(orig_bst2, dimer_dg_threshold=thresh)
                enz_mod.ENZYME_REGISTRY["bst2.0"] = mod_p
                enz_mod.ENZYME_REGISTRY["bst 2.0"] = mod_p

                cfg = build_config_from_cli(
                    primers_path=str(p_path),
                    targets_path=str(t_path),
                    temperature=63.0,
                    conc_fip_bip=c_fb,
                    conc_loop=c_lp,
                    conc_f3_b3=c_f3,
                )
                res_by_thresh[thresh] = run_analysis_pipeline(
                    cfg,
                    explicit_config=False,
                    boltzmann_max_dimers=DEFAULT_BOLTZMANN_MAX_DIMERS,
                )[0]

            r20 = res_by_thresh[-2.0]
            r35 = res_by_thresh[-3.5]
            dom20 = r20.dimer_decomposition.get("dominant_family")
            dom_str = (
                f"{dom20['label']} ({dom20['conc_nm']:.4f} nM)"
                if dom20 is not None
                else "—"
            )
            masked_items = [
                {
                    "name": r.complex_name,
                    "conc_nm": float(r.concentration) * 1e9,
                    "prob": float(r.boltzmann_3p_prob or 0.0),
                    "dg": float(r.delta_g),
                }
                for r in r20.risks
                if r.suboptimal_3p_masked
            ]
            masked_by_case[case_id] = masked_items
            rows.append(
                {
                    "case_id": case_id,
                    "bench": bench_label,
                    "tot_20": float(r20.dimer_decomposition["total_nm"]),
                    "inter_20": float(r20.dimer_decomposition["inter_nm"]),
                    "verdict_20": r20.verdict.status,
                    "cause_20": r20.verdict.global_cause,
                    "tot_35": float(r35.dimer_decomposition["total_nm"]),
                    "inter_35": float(r35.dimer_decomposition["inter_nm"]),
                    "verdict_35": r35.verdict.status,
                    "dom_20": dom_str,
                    "masked_count_cap25": len(masked_items),
                    "skipped_cap25": int(
                        r20.json_payload["totals"]["boltzmann_skipped_eligible_dimers"]
                    ),
                }
            )
    finally:
        enz_mod.ENZYME_REGISTRY["bst2.0"] = orig_bst2
        enz_mod.ENZYME_REGISTRY["bst 2.0"] = orig_bst2

    return rows, masked_by_case


def _run_cov_json_for_tree(
    src_dir: Path,
    primers: Path,
    fasta: Path,
    out_html: Path,
    panel: Optional[str],
    errors: int,
) -> Dict[str, Any]:
    """Exécute `labcraft coverage` depuis un arbre source donné (`src_dir`).
    Runs `labcraft coverage` using a specific extracted `src_dir`.
    """
    cmd = [
        sys.executable,
        "-c",
        f"import sys; sys.path.insert(0, {str(src_dir)!r}); from labcraft.cli.main import app; app()",
        "coverage",
        "-p",
        str(primers),
        "-f",
        str(fasta),
        "-o",
        str(out_html),
        "--temperature",
        "63.0",
        "--na",
        "50.0",
        "--mg",
        "8.0",
        "--dntp",
        "1.4",
        "--errors",
        str(errors),
    ]
    if panel:
        cmd.extend(["--panel", panel])
    env = dict(os.environ)
    env["PYTHONPATH"] = str(src_dir)
    subprocess.run(cmd, check=True, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    data: Dict[str, Any] = json.loads(out_html.with_suffix(".json").read_text(encoding="utf-8"))
    return data


def _classify_strain_3tier(strain: Dict[str, Any]) -> str:
    """Classe une souche en COUVERT / A_RISQUE / ECHEC / INDETERMINE."""
    if strain.get("indeterminate"):
        return "INDETERMINE"
    if not strain.get("thermo"):
        return "ECHEC"
    init_roles = ("F3", "B3", "FIP", "BIP")
    evals = strain.get("evals", {})
    has_init_tolerable = any(evals.get(r, {}).get("verdict") == "tolerable" for r in init_roles)
    has_loop_issue = any(
        evals.get(r, {}).get("verdict") in ("tolerable", "non_viable", "veto_3p", "absent")
        for r in ("LF", "LB")
        if r in evals
    )
    if has_init_tolerable or has_loop_issue:
        return "A_RISQUE"
    return "COUVERT"


def _compare_coverage_pair(label: str, j58: Dict[str, Any], j8c: Dict[str, Any]) -> Dict[str, Any]:
    """Compare deux sorties JSON `coverage` souche par souche et site par site."""
    s58 = {s["strain_id"]: s for s in j58["strains"]}
    s8c = {s["strain_id"]: s for s in j8c["strains"]}
    tier3_changes: List[Tuple[str, str, str]] = []
    site_changes: List[Tuple[str, str, Any, Any]] = []
    ddg_diffs: List[Tuple[float, float, str, str, float, float, str, str, int]] = []

    c58 = {"COUVERT": 0, "A_RISQUE": 0, "ECHEC": 0, "INDETERMINE": 0}
    c8c = {"COUVERT": 0, "A_RISQUE": 0, "ECHEC": 0, "INDETERMINE": 0}

    for sid, a in s58.items():
        b = s8c[sid]
        ca = _classify_strain_3tier(a)
        cb = _classify_strain_3tier(b)
        c58[ca] += 1
        c8c[cb] += 1
        if ca != cb:
            tier3_changes.append((sid, ca, cb))
        for role, ea in a["evals"].items():
            eb = b["evals"].get(role, {})
            if ea.get("verdict") != eb.get("verdict"):
                site_changes.append((sid, role, ea.get("verdict"), eb.get("verdict")))
            if ea.get("ddg") is not None and eb.get("ddg") is not None:
                diff = float(eb["ddg"]) - float(ea["ddg"])
                ddg_diffs.append(
                    (
                        abs(diff),
                        diff,
                        sid,
                        role,
                        float(ea["ddg"]),
                        float(eb["ddg"]),
                        str(ea["verdict"]),
                        str(eb["verdict"]),
                        int(ea["mismatches"]),
                    )
                )

    ddg_diffs.sort(key=lambda item: (-item[0], item[2], item[3]))
    nonzero = [d for d in ddg_diffs if d[0] > 1e-9]
    return {
        "label": label,
        "total_strains": int(j58["total_strains"]),
        "counts_58": c58,
        "counts_8c": c8c,
        "thermo_rate_58": float(j58["thermo_rate"]),
        "thermo_rate_8c": float(j8c["thermo_rate"]),
        "tier3_changes": len(tier3_changes),
        "site_changes": len(site_changes),
        "total_evals": len(ddg_diffs),
        "nonzero_evals": len(nonzero),
        "top5": ddg_diffs[:5],
    }


def run_section_b_coverage() -> Tuple[Dict[str, Tuple[int, str]], List[Dict[str, Any]]]:
    """Compare `58ae5d9` vs `8cf5bd8` sur les 3 génomes de référence et les 2 cohortes synthétiques."""
    with tempfile.TemporaryDirectory(prefix="labcraft_ref_bench_") as tmp_s:
        tmp = Path(tmp_s)
        wt_58 = tmp / "tree_58ae5d9"
        wt_8c = tmp / "tree_8cf5bd8"
        _extract_commit_tree("58ae5d9", wt_58)
        _extract_commit_tree("8cf5bd8", wt_8c)

        cohorts_dir = tmp / "cohorts"
        cohorts = generate_default_cohorts(cohorts_dir, seed=DEFAULT_SEED)
        d3_path, d3_n, d3_sha = cohorts["denv3"]
        d4_path, d4_n, d4_sha = cohorts["denv4"]

        runs_spec: List[Tuple[str, Path, Path, Optional[str], int]] = [
            (
                "Génome GenBank réf. DENV3 (`DEN3_M93130.fasta`, panel 3, `--errors 2`)",
                REAL_DIR / "primer_dengue_3.txt",
                GT_DIR / "DEN3_M93130.fasta",
                "3",
                2,
            ),
            (
                "Génome GenBank réf. DENV4 (`Dengue_4_S1028.fasta`, Parida, `--errors 2`)",
                GT_DIR / "parida_DENV4.txt",
                GT_DIR / "Dengue_4_S1028.fasta",
                None,
                2,
            ),
            (
                "Génome GenBank réf. YFV (`YFV_17D_17DD.fasta`, Nunes, `--errors 2`)",
                GT_DIR / "nunes_YFV.txt",
                GT_DIR / "YFV_17D_17DD.fasta",
                None,
                2,
            ),
            (
                "Cohorte synthétique DENV3 (`120` mutants, `--errors 2`, défaut CLI)",
                REAL_DIR / "primer_dengue_3.txt",
                d3_path,
                "3",
                2,
            ),
            (
                "Cohorte synthétique DENV4 (`120` mutants, `--errors 2`, défaut CLI)",
                GT_DIR / "parida_DENV4.txt",
                d4_path,
                None,
                2,
            ),
            (
                "Cohorte synthétique DENV3 (`120` mutants, `--errors 3`, boucles `3×3` incluses)",
                REAL_DIR / "primer_dengue_3.txt",
                d3_path,
                "3",
                3,
            ),
            (
                "Cohorte synthétique DENV4 (`120` mutants, `--errors 3`, boucles `3×3` incluses)",
                GT_DIR / "parida_DENV4.txt",
                d4_path,
                None,
                3,
            ),
        ]

        results: List[Dict[str, Any]] = []
        for idx, (label, pfile, ffile, pan, err) in enumerate(runs_spec):
            j58 = _run_cov_json_for_tree(
                wt_58 / "src", pfile, ffile, tmp / f"out58_{idx}.html", pan, err
            )
            j8c = _run_cov_json_for_tree(
                wt_8c / "src", pfile, ffile, tmp / f"out8c_{idx}.html", pan, err
            )
            results.append(_compare_coverage_pair(label, j58, j8c))

        cohort_meta = {
            "denv3": (d3_n, d3_sha),
            "denv4": (d4_n, d4_sha),
        }
        return cohort_meta, results


def _clear_thermo_caches() -> None:
    """Vide les caches MFE/sel/Boltzmann entre deux mesures de temps / Clears thermo caches."""
    from labcraft.thermo.backends import vienna as v_mod
    from labcraft.thermo.backends import vienna_salt as vs_mod
    from labcraft.target import unfolding as unf_mod

    v_mod._SHARED_VIENNA_CACHE.clear()
    vs_mod._SHARED_SALT_CACHE.clear()
    amp_mod._SHARED_PF_3P_CACHE.clear()
    if hasattr(unf_mod, "_UNFOLDING_CACHE"):
        unf_mod._UNFOLDING_CACHE.clear()


def run_section_c_boltzmann_perf() -> List[Dict[str, Any]]:
    """Mesure le temps et le nombre d'appels de `enrich_dimer_boltzmann_3p` sur 81 espèces (3 321 paires)."""
    from labcraft.pipeline import analysis as an_mod

    rows: List[Dict[str, Any]] = []
    with tempfile.TemporaryDirectory(prefix="labcraft_boltz_bench_") as tmp_s:
        tmp = Path(tmp_s)
        d4_id, d4_seq = read_multi_fasta(str(GT_DIR / "Dengue_4_S1028.fasta"))[0]
        yfv_id, yfv_seq = read_multi_fasta(str(GT_DIR / "YFV_17D_17DD.fasta"))[0]
        d4_win = _slice_around_amplicon(d4_seq, GT_DIR / "parida_DENV4.txt")
        yfv_win = _slice_around_amplicon(yfv_seq, GT_DIR / "nunes_YFV.txt")
        local_dup = tmp / "dup_local.fasta"
        local_dup.write_text(f">{yfv_id}\n{yfv_win}\n>{d4_id}\n{d4_win}\n", encoding="utf-8")

        scenarios = [
            (
                "81 espèces — fenêtres locales (`test_step3_redesign_v2`)",
                local_dup,
                DEFAULT_BOLTZMANN_MAX_DIMERS,
                f"défaut (`--boltzmann-max-dimers {DEFAULT_BOLTZMANN_MAX_DIMERS}`)",
            ),
            (
                "81 espèces — fenêtres locales (`test_step3_redesign_v2`)",
                local_dup,
                999999,
                "sans plafond (`999999`)",
            ),
            (
                "81 espèces — `duplex_full` génomes viraux complets",
                GT_DIR / "duplex_YFV_DENV4_targets.fasta",
                DEFAULT_BOLTZMANN_MAX_DIMERS,
                f"défaut (`--boltzmann-max-dimers {DEFAULT_BOLTZMANN_MAX_DIMERS}`)",
            ),
            (
                "81 espèces — `duplex_full` génomes viraux complets",
                GT_DIR / "duplex_YFV_DENV4_targets.fasta",
                999999,
                "sans plafond (`999999`)",
            ),
        ]

        orig_fn = an_mod.enrich_dimer_boltzmann_3p
        for scen_label, target_fasta, cap_val, mode_label in scenarios:
            cfg = build_config_from_cli(
                primers_path=str(GT_DIR / "duplex_YFV_DENV4_primers.txt"),
                targets_path=str(target_fasta),
                temperature=63.0,
                conc_fip_bip=0.8,
                conc_loop=0.4,
                conc_f3_b3=0.2,
            )

            # 1. Référence hors Boltzmann (no-op)
            _clear_thermo_caches()
            an_mod.enrich_dimer_boltzmann_3p = lambda details, *a, **kw: details  # type: ignore[assignment]
            try:
                t0 = time.perf_counter()
                run_analysis_pipeline(cfg, explicit_config=False, boltzmann_max_dimers=cap_val)
                t_without = time.perf_counter() - t0
            finally:
                an_mod.enrich_dimer_boltzmann_3p = orig_fn

            # 2. Mesure instrumentée avec Boltzmann
            calls = 0
            pf_time = 0.0

            def _wrapped(details: Dict[str, Any], *a: Any, **kw: Any) -> Dict[str, Any]:
                nonlocal calls, pf_time
                calls += 1
                ts = time.perf_counter()
                out = orig_fn(details, *a, **kw)
                pf_time += time.perf_counter() - ts
                return out

            _clear_thermo_caches()
            an_mod.enrich_dimer_boltzmann_3p = _wrapped  # type: ignore[assignment]
            try:
                t0 = time.perf_counter()
                res_with = run_analysis_pipeline(
                    cfg, explicit_config=False, boltzmann_max_dimers=cap_val
                )[0]
                t_with = time.perf_counter() - t0
            finally:
                an_mod.enrich_dimer_boltzmann_3p = orig_fn

            skipped = int(res_with.json_payload["totals"]["boltzmann_skipped_eligible_dimers"])
            rows.append(
                {
                    "scenario": scen_label,
                    "mode": mode_label,
                    "calls": calls,
                    "skipped": skipped,
                    "pf_time": pf_time,
                    "t_without": t_without,
                    "t_with": t_with,
                    "share_pct": (pf_time / t_with) * 100.0 if t_with > 0 else 0.0,
                    "overhead_pct": (pf_time / t_without) * 100.0 if t_without > 0 else 0.0,
                }
            )
    return rows


def render_markdown() -> str:
    """Exécute l'ensemble des bancs d'essai et formate `validation/reference_benchmarks.md`."""
    code_commit_short = _git_output(["log", "-1", "--format=%h", "--", "src"])
    code_commit_full = _git_output(["log", "-1", "--format=%H", "--", "src"])
    head_short = _git_output(["rev-parse", "--short", "HEAD"])

    dg3p_rows = run_section_dg3p_shifts()
    five_rows, masked_by_case = run_section_a_five_cases()
    cohort_meta, cov_rows = run_section_b_coverage()
    perf_rows = run_section_c_boltzmann_perf()

    lines: List[str] = [
        "# Bancs d'essai reproductibles de référence (`reference_benchmarks`)",
        "",
        "Document généré automatiquement par `python scripts/reference_benchmarks.py`.",
        "Toute valeur numérique citée dans `CHANGELOG.md` doit provenir exclusivement de ce fichier.",
        "",
        "## Métadonnées et conditions expérimentales",
        "",
        f"- **Commit de code évalué (`src/`)** : `{code_commit_short}` (`{code_commit_full}`)",
        f"- **HEAD lors de l'exécution** : `{head_short}`",
        "- **Température de réaction (`analyze` & `coverage`)** : `63.0 °C`",
        (
            "- **Tampon LAMP par défaut (`analyze`)** : `50.0 mM Na+`, `10.0 mM K+`, "
            "`20.0 mM Tris`, `6.0 mM Mg2+`, `1.4 mM dNTP`"
        ),
        (
            "- **Tampon explicite du comparatif `coverage` (Section B)** : `50.0 mM Na+`, "
            "`0.0 mM K+`, `0.0 mM Tris`, `8.0 mM Mg2+`, `1.4 mM dNTP`, `--ddg-max 3.0 kcal/mol`"
        ),
        (
            f"- **Profils enzymatiques (`dimer_dg_threshold`)** : "
            f"`Bst = {BST.dimer_dg_threshold:.1f}`, "
            f"`Bst 2.0 = {BST_2_0.dimer_dg_threshold:.1f}`, "
            f"`Bst 2.0 WarmStart = {BST_2_0_WARM_START.dimer_dg_threshold:.1f}`, "
            f"`GspSSD = {GSPSSD.dimer_dg_threshold:.1f}`, "
            f"`Bst 3.0 = {BST_3_0.dimer_dg_threshold:.1f}`, "
            f"`Taq = {TAQ.dimer_dg_threshold:.1f}` kcal/mol"
        ),
        (
            f"- **Seuils de verdict et Boltzmann (`analyze`)** : "
            f"`--dimer-conc-max = 0.230 nM`, "
            f"`--dimer-warn-ratio = {DEFAULT_DIMER_WARN_RATIO:.3f}` "
            f"(plancher `WARNING = {DEFAULT_DIMER_WARN_RATIO * 0.23:.3f} nM`), "
            f"`--role-imbalance-ratio = {DEFAULT_ROLE_IMBALANCE_RATIO:.1f}`, "
            f"`--boltzmann-max-dimers = {DEFAULT_BOLTZMANN_MAX_DIMERS}` "
            f"(plancher d'éligibilité `0.1 × 0.230 = 0.023 nM`)"
        ),
        "",
        "### Décalage de `local_dg` (`37.0 °C, 1000 mM Na+` → `63.0 °C, tampon LAMP par défaut`)",
        "",
        "| Motif 3' contigu (sans initiation) | `37.0 °C, 1 M Na+` (`kcal/mol`) | `63.0 °C, 1 M Na+` (`kcal/mol`) | `63.0 °C, tampon LAMP` (`kcal/mol`) | Décalage total (`kcal/mol`) |",
        "|---|---:|---:|---:|---:|",
    ]
    for label, d37, d63_1m, d63_lamp, shift in dg3p_rows:
        lines.append(
            f"| {label} | `{d37:+.2f}` | `{d63_1m:+.2f}` | `{d63_lamp:+.2f}` | `{shift:+.2f}` |"
        )

    lines.extend(
        [
            "",
            "---",
            "",
            "## A. Cinq cas de référence de paillasse (`-2.0` vs `-3.5 kcal/mol`, génomes complets)",
            "",
            "| Cas | Paillasse | Dimères amp. `-2.0` (`nM`) | dont inter `-2.0` (`nM`) | Dimères amp. `-3.5` (`nM`) | dont inter `-3.5` (`nM`) | Famille dominante (`-2.0`) | Verdict (`-2.0` / `-3.5`) | `suboptimal_3p_masked` (`cap=25`) | Éligibles non évalués (`cap=25`) |",
            "|---|---|---:|---:|---:|---:|---|---|---:|---:|",
        ]
    )
    for r in five_rows:
        lines.append(
            f"| `{r['case_id']}` | {r['bench']} | `{r['tot_20']:.4f}` | `{r['inter_20']:.4f}` | "
            f"`{r['tot_35']:.4f}` | `{r['inter_35']:.4f}` | `{r['dom_20']}` | "
            f"`{r['verdict_20']}` / `{r['verdict_35']}` | `{r['masked_count_cap25']}` | "
            f"`{r['skipped_cap25']}` |"
        )

    lines.extend(
        [
            "",
            "### Causes globales (`verdict.global_cause` à `-2.0 kcal/mol`)",
            "",
        ]
    )
    for r in five_rows:
        lines.append(f"- `{r['case_id']}` (`{r['verdict_20']}`) : {r['cause_20']}")

    lines.extend(
        [
            "",
            f"### Détail des dimères `suboptimal_3p_masked = True` (`--boltzmann-max-dimers {DEFAULT_BOLTZMANN_MAX_DIMERS}`)",
            "",
        ]
    )
    for case_id, items in masked_by_case.items():
        if not items:
            lines.append(f"- `{case_id}` : `0` dimère `suboptimal_3p_masked`.")
        else:
            lines.append(f"- `{case_id}` : `{len(items)}` dimère(s) `suboptimal_3p_masked` :")
            for it in items:
                lines.append(
                    f"  - `{it['name']}` : `conc = {it['conc_nm']:.4f} nM`, "
                    f"`P_Boltzmann = {it['prob'] * 100.0:.2f} %`, `ΔG = {it['dg']:.2f} kcal/mol`"
                )

    d3_n, d3_sha = cohort_meta["denv3"]
    d4_n, d4_sha = cohort_meta["denv4"]
    lines.extend(
        [
            "",
            "---",
            "",
            "## B. Comparatif `labcraft coverage` : commit `58ae5d9` vs commit `8cf5bd8`",
            "",
            "Cohortes synthétiques de non-régression générées par `scripts/generate_coverage_cohorts.py` "
            f"(graine fixe `SEED = {DEFAULT_SEED}` — mutants synthétiques et non souches cliniques de terrain) :",
            f"- `synthetic_cohort_denv3_120.fasta` : `{d3_n}` séquences, `sha256 = {d3_sha}`",
            f"- `synthetic_cohort_denv4_120.fasta` : `{d4_n}` séquences, `sha256 = {d4_sha}`",
            "",
            "| Jeu évalué | `58ae5d9` (`COUVERT / A_RISQUE / ECHEC`) | `8cf5bd8` (`COUVERT / A_RISQUE / ECHEC`) | Changements de classe | Changements `SiteVerdict` | Évaluations `|Δ(ΔΔG)| > 0` |",
            "|---|---|---|---|---:|---:|",
        ]
    )
    for c in cov_rows:
        a = c["counts_58"]
        b = c["counts_8c"]
        lines.append(
            f"| {c['label']} | `{a['COUVERT']} / {a['A_RISQUE']} / {a['ECHEC']}` (`{c['thermo_rate_58']:.1f} %`) | "
            f"`{b['COUVERT']} / {b['A_RISQUE']} / {b['ECHEC']}` (`{c['thermo_rate_8c']:.1f} %`) | "
            f"`{c['tier3_changes']} / {c['total_strains']}` | `{c['site_changes']}` | "
            f"`{c['nonzero_evals']} / {c['total_evals']}` |"
        )

    lines.extend(
        [
            "",
            "### Top 5 des plus grands `|Δ(ΔΔG)|` (`8cf5bd8 - 58ae5d9`, à `63.0 °C`) sur les cohortes synthétiques",
            "",
        ]
    )
    for c in cov_rows:
        if c["total_strains"] < 100:
            continue
        lines.append(f"#### {c['label']}")
        for _abs_d, diff, sid, role, d58, d8c, v58, v8c, nmm in c["top5"]:
            lines.append(
                f"- `{sid}` [`{role}`, `{nmm}` mm] : `ddG(58ae5d9) = {d58:+.4f}` → "
                f"`ddG(8cf5bd8) = {d8c:+.4f} kcal/mol` (`Δ(ΔΔG) = {diff:+.4f} kcal/mol`, "
                f"`{v58}` → `{v8c}`)"
            )
        lines.append("")

    lines.extend(
        [
            "---",
            "",
            "## C. Performance `enrich_dimer_boltzmann_3p` sur 81 espèces (`3 321` paires)",
            "",
            "| Scénario | Mode | Appels Boltzmann | Éligibles non évalués | `pf_time` (`s`) | Total hors Boltzmann (`s`) | Total avec Boltzmann (`s`) | Part `pf / avec` | Surcoût `pf / hors` |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for p in perf_rows:
        lines.append(
            f"| {p['scenario']} | {p['mode']} | `{p['calls']}` | `{p['skipped']}` | "
            f"`{p['pf_time']:.3f}` | `{p['t_without']:.2f}` | `{p['t_with']:.2f}` | "
            f"`{p['share_pct']:.2f} %` | `{p['overhead_pct']:+.2f} %` |"
        )
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run reproducible LabCraft reference benchmarks and write markdown report."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "validation" / "reference_benchmarks.md",
        help="Output Markdown path.",
    )
    args = parser.parse_args()

    md = render_markdown()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(md, encoding="utf-8")
    print(md, end="")


if __name__ == "__main__":
    main()
