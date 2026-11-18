"""Tests unitaires et d'intégration de l'Étape 4 : pipeline modulaire et optimiseur salin.

Unit and integration tests for Step 4: modular `labcraft.pipeline.analysis` Python API,
salt/enzyme-aware `optimize_primer`, and strict CI typing enforcement.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from labcraft.cli.config import BufferConfig, ExperimentConfig, PanelConfig, PrimerConfig, PrimerSetConfig
from labcraft.cli.main import app
from labcraft.cli.parsers import build_config_from_cli, read_multi_fasta, reverse_complement
from labcraft.diagnostics.enzyme import BST_2_0, PolymeraseProfile
from labcraft.lamp.domains import _find_iupac_substring
from labcraft.optimize.constraints import evaluate_variant_dimers
from labcraft.optimize.optimizer import optimize_primer
from labcraft.pipeline.analysis import (
    EmptyPanelError,
    run_analysis_pipeline,
    run_single_panel_analysis,
    split_colliding_panel_configs,
)
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.thermo.salt import UnifiedSaltModel
from labcraft.thermo.vienna import dna_session

REPO_ROOT = Path(__file__).resolve().parent.parent
GT_DIR = REPO_ROOT / "tests" / "fixtures" / "ground_truth"


def _slice_around_amplicon(target_seq: str, primer_file: Path, pad: int = 100) -> str:
    """Extrait une fenêtre locale autour de l'amplicon LAMP / Slices local window around LAMP amplicon."""
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


def test_pipeline_analysis_api_parity_with_cli(tmp_path: Path) -> None:
    """Vérifie la parité numérique exacte entre l'API Python `run_analysis_pipeline` et `labcraft analyze`.
    Verifies exact numerical parity between `run_analysis_pipeline` Python API and `labcraft analyze` CLI.
    """
    d4_id, d4_seq = read_multi_fasta(str(GT_DIR / "Dengue_4_S1028.fasta"))[0]
    d4_win = _slice_around_amplicon(d4_seq, GT_DIR / "parida_DENV4.txt")
    target_path = tmp_path / "DENV4.fasta"
    target_path.write_text(f">{d4_id}\n{d4_win}\n", encoding="utf-8")
    primers_path = GT_DIR / "parida_DENV4.txt"

    # 1. Appel direct via l'API Python modulaire / Direct call via modular Python API
    cfg = build_config_from_cli(
        primers_path=str(primers_path),
        targets_path=str(target_path),
        temperature=63.0,
        conc_fip_bip=0.8,
        conc_loop=0.4,
        conc_f3_b3=0.2,
    )
    logs: list[str] = []
    results = run_analysis_pipeline(
        cfg,
        explicit_config=False,
        progress_callback=logs.append,
    )
    assert len(results) == 1
    api_res = results[0]
    assert any("Conditions utilisées" in line for line in logs)

    # 2. Appel via le CLI Typer / Call via Typer CLI
    out_html = tmp_path / "cli_d4.html"
    runner = CliRunner()
    cli_res = runner.invoke(
        app,
        [
            "analyze",
            "-p", str(primers_path),
            "-t", str(target_path),
            "--temperature", "63.0",
            "--conc-fip-bip", "0.8",
            "--conc-loop", "0.4",
            "--conc-f3-b3", "0.2",
            "-o", str(out_html),
        ],
    )
    assert cli_res.exit_code == 0, cli_res.stdout
    cli_json = json.loads((tmp_path / "cli_d4.json").read_text(encoding="utf-8"))

    # Parité machine (1e-12) sur le bilan d'équilibre et le verdict
    # Machine parity (1e-12) on equilibrium totals and verdict
    assert api_res.verdict.status == cli_json["verdict"]
    assert abs(
        api_res.dimer_decomposition["total_nm"]
        - cli_json["totals"]["amplifiable_dimers_total_nm"]
    ) < 1e-12
    for p_name, p_fracs in cli_json["fractions"].items():
        assert abs(api_res.fractions[p_name].free - p_fracs["free"]) < 1e-12


def test_pipeline_collision_splitting_and_empty_panel() -> None:
    """Vérifie le découpage automatique des panels entrant en collision et l'erreur `EmptyPanelError`.
    Verifies automatic splitting of colliding panels and `EmptyPanelError` on empty panels.
    """
    with pytest.raises(EmptyPanelError):
        run_single_panel_analysis(PanelConfig(), {"Compétition Sans Cible": ""})

    cfg = PanelConfig(
        experiment=ExperimentConfig(temperature_C=65.0, buffer=BufferConfig()),
        targets=[],
        primer_sets=[],
    )
    # Deux jeux visant la même cible virtuelle / Two primer sets targeting the same virtual target
    from labcraft.cli.config import TargetConfig

    cfg_coll = PanelConfig(
        experiment=ExperimentConfig(temperature_C=65.0),
        targets=[TargetConfig(id="T1", sequence_file="dummy.fasta")],
        primer_sets=[
            PrimerSetConfig(target="T1", panel_name="P1", primers={"F3": PrimerConfig(seq="ACGTACGTACGTACGT")}),
            PrimerSetConfig(target="T1", panel_name="P2", primers={"F3": PrimerConfig(seq="TGCATGCATGCATGCA")}),
        ],
    )
    split_cfgs, had_coll = split_colliding_panel_configs(cfg_coll, force_joint=False)
    assert had_coll is True
    assert len(split_cfgs) == 2
    joint_cfgs, had_coll_joint = split_colliding_panel_configs(cfg_coll, force_joint=True)
    assert had_coll_joint is False
    assert len(joint_cfgs) == 1
    assert len(cfg.primer_sets) == 0


def test_optimize_primer_accepts_salt_backend_and_custom_enzyme() -> None:
    """Vérifie que `optimize_primer` prend en compte le `DuplexEnergyBackend` salin et le `PolymeraseProfile`.
    Verifies `optimize_primer` honors the configured salt `DuplexEnergyBackend` and `PolymeraseProfile`.
    """
    panel = [
        {"name": "P_A", "sequence": "TTTTTTTGCATGCAT"},
        {"name": "P_B", "sequence": "AAAAAAAATGCATGC"},
    ]
    salt_backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    salt_kwargs = {"na_mm": 50.0, "k_mm": 0.0, "tris_mm": 0.0, "mg_mm": 8.0, "dntp_mm": 1.4}

    with dna_session(temp_celsius=65.0):
        # Les paramètres salins modifient l'énergie 3' évaluée par rapport aux conditions 1 M Na+
        # Salt parameters shift the evaluated 3' dimer dG compared to 1 M Na+ reference conditions
        amp_ref = evaluate_variant_dimers(panel[0]["sequence"], panel, 0, BST_2_0, 65.0)
        amp_salt = evaluate_variant_dimers(
            panel[0]["sequence"], panel, 0, BST_2_0, 65.0, backend=salt_backend, **salt_kwargs
        )
        assert len(amp_ref) == 1 and len(amp_salt) == 1
        assert amp_ref[0][1] != pytest.approx(amp_salt[0][1], abs=1e-3)

    strict_enzyme = PolymeraseProfile(
        name="StrictCustomBst",
        strand_displacement=True,
        three_prime_window=5,
        three_prime_absolute_window=2,
        dimer_dg_threshold=-2.5,
        warm_start=True,
    )
    ranked = optimize_primer(
        primer_name="P_A",
        primers=panel,
        temp_celsius=65.0,
        max_mutations=2,
        window_3p=6,
        enzyme=strict_enzyme,
        backend=salt_backend,
        **salt_kwargs,
    )
    assert len(ranked) >= 1
    best = ranked[0]
    assert 1 <= best["num_mutations"] <= 2

    # Vérifier que le meilleur variant n'a plus aucun dimère amplifiable dans le tampon salin configuré
    # Verify the best variant has zero amplifiable dimers in the configured salt buffer
    with dna_session(temp_celsius=65.0):
        rem = evaluate_variant_dimers(
            best["sequence"], panel, 0, strict_enzyme, 65.0, backend=salt_backend, **salt_kwargs
        )
        assert rem == []


def test_ci_workflow_enforces_blocking_mypy() -> None:
    """Vérifie que `.github/workflows/ci.yml` exécute `mypy` sans `continue-on-error: true`.
    Verifies `.github/workflows/ci.yml` runs `mypy` without `continue-on-error: true`.
    """
    ci_text = (REPO_ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
    assert "mypy src/labcraft" in ci_text
    assert "continue-on-error" not in ci_text
