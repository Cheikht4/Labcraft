"""Tests du Point 1 : cohérence du seuil ΔG3' (`dimer_dg_threshold`) et de `--dimer-conc-max`.

Tests for Point 1: coherence of local 3' dimer ΔG threshold (`dimer_dg_threshold`) and
`--dimer-conc-max` documentation across `enzyme.py`, `README.md`, and CLI `--help`, plus
non-normative regression snapshot at `-2.0` vs `-3.5 kcal/mol`.
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest
from typer.testing import CliRunner

from labcraft.cli.main import app
from labcraft.cli.parsers import build_config_from_cli, read_multi_fasta, reverse_complement
from labcraft.diagnostics import enzyme as enz_mod
from labcraft.diagnostics.enzyme import (
    BST,
    BST_2_0,
    BST_2_0_WARM_START,
    BST_3_0,
    GSPSSD,
    TAQ,
)
from labcraft.lamp.domains import _find_iupac_substring
from labcraft.pipeline.analysis import run_analysis_pipeline

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
GT_DIR = REPO_ROOT / "tests" / "fixtures" / "ground_truth"


def _slice_around_amplicon(target_seq: str, primer_file: Path, pad: int = 100) -> str:
    """Découpe une fenêtre locale autour de l'amplicon LAMP / Slices local window around LAMP amplicon."""
    recs = dict(read_multi_fasta(str(primer_file)))
    spans = []
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


def test_enzyme_dg3p_thresholds_code_doc_cli_coherence() -> None:
    """Vérifie l'alignement exact des seuils ΔG3' (-2.0, -1.5, -3.0) et de la mention d'heuristique non calibrée.
    Verifies exact alignment of ΔG3' thresholds (-2.0, -1.5, -3.0) and uncalibrated heuristic wording.
    """
    assert BST.dimer_dg_threshold == -2.0
    assert BST_2_0.dimer_dg_threshold == -2.0
    assert BST_2_0_WARM_START.dimer_dg_threshold == -2.0
    assert GSPSSD.dimer_dg_threshold == -2.0
    assert BST_3_0.dimer_dg_threshold == -1.5
    assert TAQ.dimer_dg_threshold == -3.0

    # Le code doit documenter l'heuristique non calibrée, sans mention de co-calibration ni PrimerExplorer
    # Code comments must state uncalibrated heuristic, without co-calibration or PrimerExplorer comparison
    enzyme_src = (REPO_ROOT / "src" / "labcraft" / "diagnostics" / "enzyme.py").read_text(encoding="utf-8")
    assert "heuristique non calibrée" in enzyme_src
    assert "ne permettent pas" in enzyme_src and "discriminer entre -2.0 et -3.5" in enzyme_src
    assert "co-calibr" not in enzyme_src
    assert "PrimerExplorer" not in enzyme_src

    readme = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    assert "-2.0 kcal/mol" in readme
    assert "-1.5 kcal/mol" in readme
    assert "-3.0 kcal/mol" in readme
    assert "heuristique non calibrée ; les données de paillasse actuelles ne permettent pas de discriminer entre -2.0 et -3.5" in readme
    assert "0.154 – 0.336 nM" in readme
    assert "FIP_Nunes_YFV × F3_ParidaD4" in readme
    assert "co-calibr" not in readme
    assert "facteur 2,3" not in readme

    runner = CliRunner()
    for subcmd in ("analyze", "coverage", "redesign"):
        res = runner.invoke(app, [subcmd, "--help"])
        assert res.exit_code == 0
        clean = " ".join(res.stdout.replace("│", " ").split())
        assert "-2.0 kcal/mol" in clean
        assert "-1.5 kcal/mol" in clean
        assert "-3.0 kcal/mol" in clean
        assert "heuristique non calibrée" in clean
        assert "discriminer entre -2.0 et -3.5" in clean
        assert "co-calibr" not in clean

    for subcmd in ("analyze", "redesign"):
        res = runner.invoke(app, [subcmd, "--help"])
        clean = " ".join(res.stdout.replace("│", " ").split())
        assert "0.154 – 0.336 nM" in clean
        assert "FIP_Nunes_YFV × F3_ParidaD4" in clean
        assert "facteur 2,3" not in clean


@pytest.mark.regression
def test_regression_separation_at_minus_2_0_vs_minus_3_5_on_reference_cases(tmp_path: Path) -> None:
    """TEST DE RÉGRESSION (ET NON DE VALIDATION BIOLOGIQUE) : photographie le comportement actuel
    du seuil provisoire `0.23 nM` à `-2.0` et `-3.5 kcal/mol` sur DENV4 et le duplex YFV+DENV4.
    Ce test montre précisément que ces cas ne contraignent pas `dimer_dg_threshold` entre `-2.0` et
    `-3.5` ; il doit être mis à jour librement lors d'une future recalibration sur cohorte élargie
    (ou désélectionné via `pytest -m "not regression"`).

    REGRESSION SNAPSHOT (NOT BIOLOGICAL VALIDATION): records the current behavior of the provisional
    `0.23 nM` cutoff at both `-2.0` and `-3.5 kcal/mol` on DENV4 and the YFV+DENV4 duplex.
    This test illustrates that current reference cases do not constrain `dimer_dg_threshold` between
    `-2.0` and `-3.5`; update freely when recalibrating on an expanded cohort (or skip via
    `pytest -m "not regression"`).
    """
    d4_id, d4_seq = read_multi_fasta(str(GT_DIR / "Dengue_4_S1028.fasta"))[0]
    yfv_id, yfv_seq = read_multi_fasta(str(GT_DIR / "YFV_17D_17DD.fasta"))[0]
    d4_win = _slice_around_amplicon(d4_seq, GT_DIR / "parida_DENV4.txt")
    yfv_win = _slice_around_amplicon(yfv_seq, GT_DIR / "nunes_YFV.txt")

    p_d4 = tmp_path / "d4.fasta"
    p_dup = tmp_path / "dup.fasta"
    p_d4.write_text(f">{d4_id}\n{d4_win}\n", encoding="utf-8")
    p_dup.write_text(f">{yfv_id}\n{yfv_win}\n>{d4_id}\n{d4_win}\n", encoding="utf-8")

    orig_bst2 = enz_mod.BST_2_0
    try:
        for thresh in (-2.0, -3.5):
            mod_p = replace(orig_bst2, dimer_dg_threshold=thresh)
            enz_mod.ENZYME_REGISTRY["bst2.0"] = mod_p
            enz_mod.ENZYME_REGISTRY["bst 2.0"] = mod_p

            cfg_ok = build_config_from_cli(
                primers_path=str(GT_DIR / "parida_DENV4.txt"),
                targets_path=str(p_d4),
                temperature=63.0,
                conc_fip_bip=0.8,
                conc_loop=0.4,
                conc_f3_b3=0.2,
            )
            res_ok = run_analysis_pipeline(cfg_ok, explicit_config=False)[0]
            assert res_ok.dimer_decomposition["total_nm"] < 0.23
            assert res_ok.verdict.status != "FAILURE"

            cfg_fail = build_config_from_cli(
                primers_path=str(GT_DIR / "duplex_YFV_DENV4_primers.txt"),
                targets_path=str(p_dup),
                temperature=63.0,
                conc_fip_bip=0.4,
                conc_loop=0.2,
                conc_f3_b3=0.1,
            )
            res_fail = run_analysis_pipeline(cfg_fail, explicit_config=False)[0]
            assert res_fail.dimer_decomposition["total_nm"] > 0.23
            assert res_fail.verdict.status == "FAILURE"
    finally:
        enz_mod.ENZYME_REGISTRY["bst2.0"] = orig_bst2
        enz_mod.ENZYME_REGISTRY["bst 2.0"] = orig_bst2
