"""Tests CLI de bout en bout pour les correctifs de la Phase 1.

End-to-end CLI tests via CliRunner for Phase 1 fixes (default dimer threshold,
version/metadata, sequence alphabet validation, and multiplex report decomposition).
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

from labcraft import __version__
from labcraft.cli.main import app
from labcraft.cli.parsers import read_multi_fasta, reverse_complement
from labcraft.lamp.domains import _find_iupac_substring

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
GT_DIR = REPO_ROOT / "tests" / "fixtures" / "ground_truth"


def _slice_around_amplicon(target_seq: str, primer_file: Path, pad: int = 100) -> str:
    """Extrait la fenêtre cible autour de l'amplicon pour accélérer les tests CLI.
    Slices target sequence tightly around the LAMP amplicon to speed up CLI tests.
    """
    recs = dict(read_multi_fasta(str(primer_file)))
    spans = []
    for name, seq in recs.items():
        if name.endswith(("FIP", "BIP")):
            continue
        for cand in (seq, reverse_complement(seq)):
            pos = _find_iupac_substring(cand, target_seq)
            if pos >= 0:
                spans.append((pos, pos + len(cand)))
    assert spans, f"Aucun domaine trouvé pour {primer_file.name}"
    lo = max(0, min(a for a, _ in spans) - pad)
    hi = min(len(target_seq), max(b for _, b in spans) + pad)
    return target_seq[lo:hi]


@pytest.fixture(scope="module")
def truncated_gt_targets(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """Prépare les cibles de vérité terrain tronquées autour de chaque amplicon.
    Prepares truncated ground-truth target FASTA files around each amplicon.
    """
    tdir = tmp_path_factory.mktemp("gt_trunc")

    d3_id, d3_seq = read_multi_fasta(str(GT_DIR / "DEN3_M93130.fasta"))[0]
    d4_id, d4_seq = read_multi_fasta(str(GT_DIR / "Dengue_4_S1028.fasta"))[0]
    yfv_id, yfv_seq = read_multi_fasta(str(GT_DIR / "YFV_17D_17DD.fasta"))[0]

    d3_win = _slice_around_amplicon(d3_seq, GT_DIR / "primer_dengue_3_panel3.txt")
    d4_win = _slice_around_amplicon(d4_seq, GT_DIR / "parida_DENV4.txt")
    yfv_win = _slice_around_amplicon(yfv_seq, GT_DIR / "nunes_YFV.txt")

    p_d3 = tdir / "DEN3.fasta"
    p_d4 = tdir / "DENV4.fasta"
    p_yfv = tdir / "YFV.fasta"
    p_dup = tdir / "DUPLEX.fasta"

    p_d3.write_text(f">{d3_id}\n{d3_win}\n", encoding="utf-8")
    p_d4.write_text(f">{d4_id}\n{d4_win}\n", encoding="utf-8")
    p_yfv.write_text(f">{yfv_id}\n{yfv_win}\n", encoding="utf-8")
    p_dup.write_text(f">{yfv_id}\n{yfv_win}\n>{d4_id}\n{d4_win}\n", encoding="utf-8")
    return {"d3": p_d3, "d4": p_d4, "yfv": p_yfv, "dup": p_dup}


@pytest.fixture(scope="module")
def duplex_half_cli_run(
    truncated_gt_targets: dict[str, Path],
    tmp_path_factory: pytest.TempPathFactory,
) -> tuple[int, str, str]:
    """Exécute une fois le duplex demi-concentration via CliRunner pour partage entre tests.
    Runs the half-concentration duplex CLI case once for reuse across assertions.
    """
    out_dir = tmp_path_factory.mktemp("dup_half")
    out_html = out_dir / "duplex_half.html"
    runner = CliRunner()
    res = runner.invoke(
        app,
        [
            "analyze",
            "-p", str(GT_DIR / "duplex_YFV_DENV4_primers.txt"),
            "-t", str(truncated_gt_targets["dup"]),
            "--temperature", "63.0",
            "--conc-fip-bip", "0.4",
            "--conc-loop", "0.2",
            "--conc-f3-b3", "0.1",
            "-o", str(out_html),
        ],
    )
    return res.exit_code, res.stdout, out_html.read_text(encoding="utf-8")


def _sum_amplifiable_families_nm(html: str) -> float:
    """Additionne les concentrations (en nM) des familles amplifiables affichées dans le HTML.
    Sums amplifiable family concentrations (in nM) from the HTML decomposition table.
    """
    vals = re.findall(r"<td>([0-9]+\.[0-9]{4})\s+nM</td>", html)
    return sum(float(v) for v in vals)


@pytest.mark.parametrize(
    ("case_id", "primers_file", "target_key", "c_fip", "c_loop", "c_f3", "ref_nm", "expect_failure"),
    [
        ("singleplex_denv3", "primer_dengue_3_panel3.txt", "d3", "0.8", "0.4", "0.2", 0.0, False),
        ("singleplex_denv4", "parida_DENV4.txt", "d4", "0.8", "0.4", "0.2", 0.1540, False),
        ("singleplex_yfv", "nunes_YFV.txt", "yfv", "0.8", "0.4", "0.2", 0.0011, False),
        ("duplex_half", "duplex_YFV_DENV4_primers.txt", "dup", "0.4", "0.2", "0.1", 0.3514, True),
        ("duplex_full", "duplex_YFV_DENV4_primers.txt", "dup", "0.8", "0.4", "0.2", 1.3498, True),
    ],
)
def test_phase1_ground_truth_default_dimer_threshold(
    case_id: str,
    primers_file: str,
    target_key: str,
    c_fip: str,
    c_loop: str,
    c_f3: str,
    ref_nm: float,
    expect_failure: bool,
    truncated_gt_targets: dict[str, Path],
    duplex_half_cli_run: tuple[int, str, str],
    tmp_path: Path,
) -> None:
    """Les 3 singleplex réussis ne sont pas en FAILURE, les 2 duplex échoués sont en FAILURE.
    The 3 working singleplex assays are not FAILURE; the 2 failed duplex assays are FAILURE.
    """
    if case_id == "duplex_half":
        exit_code, stdout, html = duplex_half_cli_run
    else:
        out_html = tmp_path / f"{case_id}.html"
        runner = CliRunner()
        res = runner.invoke(
            app,
            [
                "analyze",
                "-p", str(GT_DIR / primers_file),
                "-t", str(truncated_gt_targets[target_key]),
                "--temperature", "63.0",
                "--conc-fip-bip", c_fip,
                "--conc-loop", c_loop,
                "--conc-f3-b3", c_f3,
                "-o", str(out_html),
            ],
        )
        exit_code, stdout, html = res.exit_code, res.stdout, out_html.read_text(encoding="utf-8")

    assert exit_code == 0, stdout
    if expect_failure:
        assert "Verdict Global : FAILURE" in html
    else:
        assert "Verdict Global : FAILURE" not in html

    total_nm = _sum_amplifiable_families_nm(html)
    if ref_nm == 0.0:
        assert total_nm == 0.0
        # Aucune recommandation sur les dimères quand le total est nul
        # No dimer recommendation when total amplifiable dimers is zero
        rec_section = html.split("<h2>Recommandations</h2>")[1].split("<h2>Concentrations Suggérées</h2>")[0]
        assert "Famille de dimères amplifiables dominante" not in html
        assert "WarmStart" not in rec_section
        assert "TaqMan" not in rec_section
    else:
        assert abs(total_nm - ref_nm) / ref_nm < 0.01, f"{case_id}: {total_nm} vs ref {ref_nm}"


def test_phase1_duplex_report_decomposition_and_recommendations_cli(
    duplex_half_cli_run: tuple[int, str, str],
) -> None:
    """Vérifie le rapport multiplexe : 1 ligne par site, famille dominante, sans LNA ni textes génériques.
    Verifies duplex HTML report: one row per site, dominant family named, no LNA banner or generic advice.
    """
    exit_code, stdout, html = duplex_half_cli_run
    assert exit_code == 0, stdout

    # 1.4.1 Une seule ligne par site dans le tableau d'accessibilité
    # 1.4.1 Single row per site in the accessibility table
    acc_section = html.split("Accessibilité à l'initiation de la cible")[1].split("Matrice d'Interactions")[0]
    for expected_prefix in (
        "F3_Nunes_YFV", "B3_Nunes_YFV", "FIP_Nunes_YFV", "BIP_Nunes_YFV", "LF_Nunes_YFV", "LB_Nunes_YFV",
        "F3_ParidaD4", "B3_ParidaD4", "FIP_ParidaD4", "BIP_ParidaD4", "LF_ParidaD4", "LB_ParidaD4",
    ):
        matches = re.findall(rf"{re.escape(expected_prefix)}(?:@[^\s<]+)?_site\b", acc_section)
        assert len(matches) == 1, f"Site '{expected_prefix}' trouvé {len(matches)} fois dans l'accessibilité"

    # 1.4.2 Absence du bandeau LNA sur un panel sans LNA
    # 1.4.2 No LNA banner on a panel without LNA bases
    assert "Support LNA (Phase 2)" not in html

    # 1.4.3 Décomposition multiplexe attendue (inter-panels 0.31 / 0.35, FIP_Nunes_YFV × F3_ParidaD4 à 0.29 nM)
    # 1.4.3 Expected multiplex decomposition (inter-panel 0.31 / 0.35, FIP_Nunes_YFV x F3_ParidaD4 at 0.29 nM)
    assert "Décomposition des Dimères Amplifiables" in html
    assert "0.35 nM" in html
    assert "0.31 nM" in html
    assert "FIP_Nunes_YFV × F3_ParidaD4" in html
    assert "0.29 nM" in html

    # 1.4.4 Recommandations ciblées nommant la famille dominante et l'amorce étendue en 3' (sans textes génériques)
    # 1.4.4 Targeted recommendations naming dominant family and 3'-extended primer (no generic texts)
    rec_section = html.split("<h2>Recommandations</h2>")[1].split("<h2>Concentrations Suggérées</h2>")[0]
    assert "FIP_Nunes_YFV × F3_ParidaD4" in rec_section
    assert "F3_ParidaD4" in rec_section
    assert "WarmStart" not in rec_section
    assert "TaqMan" not in rec_section


def test_phase1_lna_banner_shown_only_when_lna_present(tmp_path: Path) -> None:
    """Le bandeau 'Support LNA (Phase 2)' apparaît uniquement lorsqu'une amorce contient une base LNA (+).
    The 'Support LNA (Phase 2)' banner appears only when at least one primer contains an LNA base (+).
    """
    primers_file = tmp_path / "lna_primers.txt"
    target_file = tmp_path / "lna_target.fasta"
    out_html = tmp_path / "lna_report.html"

    target_seq = (
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "CATGTCGTAGCTAGCATGCA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "GGTCAACTGCTAGCTAGCTA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    )
    target_file.write_text(f">T1\n{target_seq}\n", encoding="utf-8")
    primers_file.write_text(
        ">P_F3\nCATG+TCGTAGCTAGCATGCA\n"
        ">P_B3\nTAGCTAGCTAGCAGTTGACC\n",
        encoding="utf-8",
    )

    runner = CliRunner()
    res = runner.invoke(
        app,
        ["analyze", "-p", str(primers_file), "-t", str(target_file), "--temperature", "63.0", "-o", str(out_html)],
    )
    assert res.exit_code == 0, res.stdout
    assert "Support LNA (Phase 2)" in out_html.read_text(encoding="utf-8")


def test_phase1_version_metadata_and_help() -> None:
    """Vérifie labcraft --version, l'identité d'auteur, le README sans badge PyPI et l'aide --dimer-conc-max.
    Checks labcraft --version, author metadata, README without PyPI badge, and --dimer-conc-max help text.
    """
    runner = CliRunner()
    res_ver = runner.invoke(app, ["--version"])
    assert res_ver.exit_code == 0
    assert __version__ in res_ver.stdout

    res_help = runner.invoke(app, ["analyze", "--help"])
    assert res_help.exit_code == 0
    help_clean = " ".join(res_help.stdout.replace("│", " ").split())
    assert "0.23" in help_clean
    assert "calibré sur 5 cas de paillasse" in help_clean
    assert "facteur 2,3" in help_clean

    citation = (REPO_ROOT / "CITATION.cff").read_text(encoding="utf-8")
    pyproject = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    readme = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    assert "Touré" in citation
    assert "Touré" in pyproject
    assert "Diagne" not in citation
    assert "Diagne" not in pyproject
    assert "![PyPI]" not in readme
    assert "calibré sur 5 cas de paillasse" in readme


def test_phase1_alphabet_validation_cli(tmp_path: Path) -> None:
    """Une amorce contenant 'X' renvoie exit code 1 avec le nom de l'amorce ; '-' dans la cible avertit.
    A primer containing 'X' fails with exit code 1 naming the primer; '-' in target strips and warns.
    """
    bad_primers = tmp_path / "bad_primers.txt"
    bad_primers.write_text(">BadPanel_F3\nACGTXGCTAGCTAGCT\n", encoding="utf-8")

    ok_target = tmp_path / "target.fasta"
    ok_target.write_text(">T1\nAAAAACGTAGCTAGCTAGCTAAAAA\n", encoding="utf-8")

    runner = CliRunner()
    res_bad = runner.invoke(
        app,
        ["analyze", "-p", str(bad_primers), "-t", str(ok_target), "-o", str(tmp_path / "out.html")],
    )
    assert res_bad.exit_code == 1
    assert "BadPanel_F3" in res_bad.stdout
    assert "'X'" in res_bad.stdout
    assert "position 5" in res_bad.stdout

    # Cible avec tirets d'alignement et base U -> nettoyés avec avertissement donnant le nombre exact
    # Target with alignment dashes and U base -> cleaned with warning stating exact dash count
    good_primers = tmp_path / "u_primers.txt"
    good_primers.write_text(
        ">P_F3\nAUGCUCGUAGCUAGCAUGCA\n"
        ">P_B3\nUAGCUAGCUAGCAGUUGACC\n",
        encoding="utf-8",
    )
    gapped_target = tmp_path / "gapped.fasta"
    gapped_target.write_text(
        ">GappedT\nAAAAA---AUGCUCGUAGCUAGCAUGCA--AAAAAAGGTCAACTGCTAGCTAGCTAAAAAA\n",
        encoding="utf-8",
    )
    out_gap = tmp_path / "gap_report.html"
    res_gap = runner.invoke(
        app,
        ["analyze", "-p", str(good_primers), "-t", str(gapped_target), "--temperature", "63.0", "-o", str(out_gap)],
    )
    assert res_gap.exit_code == 0, res_gap.stdout
    assert "5 tiret(s)" in res_gap.stdout
    assert "5 tiret(s)" in out_gap.read_text(encoding="utf-8")

    # Vérifie aussi labcraft coverage avec amorce contenant 'X' -> exit code 1
    # Also verify labcraft coverage rejects a primer containing 'X' with exit code 1
    res_cov_bad = runner.invoke(
        app,
        ["coverage", "-p", str(bad_primers), "-f", str(ok_target), "-o", str(tmp_path / "cov_bad.html")],
    )
    assert res_cov_bad.exit_code == 1
    assert "BadPanel_F3" in res_cov_bad.stdout
    assert "'X'" in res_cov_bad.stdout
