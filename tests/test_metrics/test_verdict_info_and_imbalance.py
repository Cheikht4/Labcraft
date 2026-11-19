"""Tests du niveau `INFO` et des critères relatifs de `verdict.py` (Option B2).

Tests for the `INFO` level, twin-role occupation imbalance (`--role-imbalance-ratio`),
degeneracy dilution attribution vs thermodynamic site blockage, and relative dimer
warning floor (`--dimer-warn-ratio`).
"""
from __future__ import annotations

import json
from pathlib import Path

from typer.testing import CliRunner

from labcraft.cli.main import app
from labcraft.cli.parsers import build_config_from_cli
from labcraft.metrics.fractions import PrimerFractions
from labcraft.metrics.verdict import generate_verdict
from labcraft.pipeline.analysis import run_analysis_pipeline

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
GT_DIR = REPO_ROOT / "tests" / "fixtures" / "ground_truth"


def _free_frac(occ: float) -> PrimerFractions:
    return PrimerFractions(
        free=0.95,
        hairpin=0.0,
        homodimer=0.0,
        heterodimer_intra=0.0,
        heterodimer_inter=0.0,
        target_bound=occ,
        dominant_complex="free",
        dominant_fraction=0.95,
    )


def test_verdict_info_vs_degeneracy_warning_vs_thermodynamic_blockage() -> None:
    """Vérifie :
    1) Occupations < 10 % équilibrées entre rôles jumeaux -> `INFO` ;
    2) Déséquilibre >= 20× sur amorce dégénérée (`N_var > 1`) -> `WARNING` explicite de dégénérescence
       avec nombre de variants et concentration par variant (~3,1 nM) ;
    3) Déséquilibre >= 20× sur amorce non dégénérée (`N_var == 1`) -> `WARNING` de blocage thermodynamique.
    """
    # 1) Occupations basses mais équilibrées (ratio 4× < 20×) -> INFO
    fracs_bal = {
        "F3_P": _free_frac(0.005),
        "B3_P": _free_frac(0.020),
        "FIP_P": _free_frac(0.030),
        "BIP_P": _free_frac(0.045),
    }
    occs_bal = {
        "T1": {"F3_P_site": 0.005, "B3_P_site": 0.020, "FIP_P_site": 0.030, "BIP_P_site": 0.045}
    }
    v_info = generate_verdict(fracs_bal, occs_bal, risks=[], dimer_conc_max_nm=0.23)
    assert v_info.status == "INFO"
    assert all(i.level == "INFO" for i in v_info.issues)

    # 2) Déséquilibre >= 20× dû à 64 variants dégénérés à 0,2 µM (~3,1 nM/variant) -> WARNING dégénérescence
    fracs_deg = {"F3_Nunes_YFV": _free_frac(0.0734)}
    for k in range(1, 65):
        fracs_deg[f"B3_Nunes_YFV#{k}"] = _free_frac(0.0013 / 64.0)
    occs_deg = {"U17066.1": {"F3_Nunes_YFV_site": 0.0734, "B3_Nunes_YFV_site": 0.001315}}
    v_deg = generate_verdict(
        fracs_deg,
        occs_deg,
        risks=[],
        dimer_conc_max_nm=0.23,
        primer_nominal_concs={"F3_Nunes_YFV": 0.2e-6, "B3_Nunes_YFV": 0.2e-6},
    )
    assert v_deg.status == "WARNING"
    b3_issue = next(i for i in v_deg.issues if i.primer_name == "B3_Nunes_YFV")
    assert b3_issue.level == "WARNING"
    assert "Occupation réduite par dégénérescence" in b3_issue.cause
    assert "B3_Nunes_YFV : 64 variants, ~3,1 nM par variant" in b3_issue.cause
    assert "Blocage thermodynamique" not in b3_issue.cause

    # 3) Déséquilibre >= 20× sur amorce non dégénérée (N_var == 1) -> WARNING blocage thermodynamique
    fracs_blk = {"F3_S": _free_frac(0.080), "B3_S": _free_frac(0.002)}
    occs_blk = {"T1": {"F3_S_site": 0.080, "B3_S_site": 0.002}}
    v_blk = generate_verdict(fracs_blk, occs_blk, risks=[], dimer_conc_max_nm=0.23)
    assert v_blk.status == "WARNING"
    b3_blk = next(i for i in v_blk.issues if i.primer_name == "B3_S")
    assert b3_blk.level == "WARNING"
    assert "Blocage thermodynamique relatif du site cible pour B3_S" in b3_blk.cause
    assert "dégénérescence" not in b3_blk.cause


def test_five_reference_cases_verdict_and_cli_info_exit_zero(tmp_path: Path) -> None:
    """Vérifie le verdict des 5 cas de référence et la sortie CLI/HTML/JSON du niveau `INFO` (exit code 0) :
    - `singleplex_denv3` -> `INFO` (code retour CLI = 0, HTML + JSON `INFO`)
    - `singleplex_denv4` -> `WARNING` (dimères amplifiables >= 0.5 * 0.23 nM)
    - `singleplex_yfv`   -> `WARNING` (dégénérescence B3_Nunes_YFV : 64 variants, ~3,1 nM par variant)
    - `duplex_half`      -> `FAILURE`
    - `duplex_full`      -> `FAILURE`
    """
    runner = CliRunner()
    out_d3 = tmp_path / "d3.html"
    res_cli = runner.invoke(
        app,
        [
            "analyze",
            "-p", str(GT_DIR / "primer_dengue_3_panel3.txt"),
            "-t", str(GT_DIR / "DEN3_M93130.fasta"),
            "--temperature", "63.0",
            "--conc-fip-bip", "0.8",
            "--conc-loop", "0.4",
            "--conc-f3-b3", "0.2",
            "-o", str(out_d3),
        ],
    )
    assert res_cli.exit_code == 0, res_cli.stdout
    assert "Verdict Global : INFO" in out_d3.read_text(encoding="utf-8")
    d3_json = json.loads((tmp_path / "d3.json").read_text(encoding="utf-8"))
    assert d3_json["verdict"] == "INFO"
    assert all(iss["level"] == "INFO" for iss in d3_json["issues"])

    # Vérifier les options CLI --role-imbalance-ratio et --dimer-warn-ratio dans --help
    res_help = runner.invoke(app, ["analyze", "--help"])
    clean_help = " ".join(res_help.stdout.replace("│", " ").split())
    assert "--role-imbalance-ratio" in clean_help
    assert "--dimer-warn-ratio" in clean_help
    assert "2,0 kcal/mol" in clean_help

    cases = [
        ("singleplex_denv4", GT_DIR / "parida_DENV4.txt", GT_DIR / "Dengue_4_S1028.fasta", 0.8, 0.4, 0.2, "WARNING", "Dimères amplifiables modérés"),
        ("singleplex_yfv", GT_DIR / "nunes_YFV.txt", GT_DIR / "YFV_17D_17DD.fasta", 0.8, 0.4, 0.2, "WARNING", "B3_Nunes_YFV : 16 variants, ~12,5 nM par variant"),
        ("duplex_half", GT_DIR / "duplex_YFV_DENV4_primers.txt", GT_DIR / "duplex_YFV_DENV4_targets.fasta", 0.4, 0.2, 0.1, "FAILURE", "dimères d'amorces amplifiables majeurs"),
        ("duplex_full", GT_DIR / "duplex_YFV_DENV4_primers.txt", GT_DIR / "duplex_YFV_DENV4_targets.fasta", 0.8, 0.4, 0.2, "FAILURE", "dimères d'amorces amplifiables majeurs"),
    ]
    for name, pf, tf, c_fb, c_lp, c_f3, exp_status, exp_substr in cases:
        cfg = build_config_from_cli(
            primers_path=str(pf),
            targets_path=str(tf),
            temperature=63.0,
            conc_fip_bip=c_fb,
            conc_loop=c_lp,
            conc_f3_b3=c_f3,
        )
        res = run_analysis_pipeline(cfg, explicit_config=False)[0]
        assert res.verdict.status == exp_status, f"{name}: got {res.verdict.status} ({res.verdict.global_cause})"
        assert exp_substr in res.verdict.global_cause, f"{name}: got cause {res.verdict.global_cause}"
