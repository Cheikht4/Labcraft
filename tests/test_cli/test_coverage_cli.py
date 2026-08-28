import pytest
from typer.testing import CliRunner
import tempfile
import csv
import json
import subprocess
import sys
import time
from pathlib import Path

from labcraft.cli.main import app
from labcraft.cli.config import PanelConfig, ExperimentConfig, BufferConfig


def test_coverage_cli_with_csv(tmp_path: Path):
    """Test du mode couverture avec fichier CSV de sites pré-calculé."""
    fasta_file = tmp_path / "strains.fasta"
    dummy_seq = "GAAGAAGCTGTGCAGCCTG" + "N"*20 + "TGGCTTTTGGGCCTGACTTC" + "N"*20 + "GCACGGTGTAGCAGACTAG" + "N"*20 + "TCCCCACGACGGAGCTACAG"
    
    fasta_text = ">S1_parfaite\n" + dummy_seq + "\n>S2_mism_interne\n" + dummy_seq + "\n>S3_veto_3p\n" + dummy_seq + "\n>S4_loop_fail\n" + dummy_seq + "\n>S5_count_lost\n" + dummy_seq + "\n"
    fasta_file.write_text(fasta_text)
    
    csv_file = tmp_path / "sites.csv"
    with open(csv_file, "w") as f:
        writer = csv.DictWriter(f, fieldnames=["strain_id", "primer_role", "position", "strand", "site_seq", "n_mismatches", "panel"])
        writer.writeheader()
        
        strains = ["S1_parfaite", "S2_mism_interne", "S3_veto_3p", "S4_loop_fail", "S5_count_lost"]
        
        for s in strains:
            f3_seq = "GCCACCTTAAGCCACAGTA"
            f3_mm = 0
            if s == "S2_mism_interne":
                f3_seq = "GCCACCTTAAGCGACAGTA"
                f3_mm = 1
            elif s == "S3_veto_3p":
                f3_seq = "GCCACCTTAAGCCACAGTT"
                f3_mm = 1
            elif s == "S5_count_lost":
                f3_mm = 3
                
            writer.writerow({"strain_id": s, "primer_role": "F3", "position": 0, "strand": "+", "site_seq": f3_seq, "n_mismatches": f3_mm, "panel": "Test"})
            writer.writerow({"strain_id": s, "primer_role": "B3", "position": 0, "strand": "-", "site_seq": "CCCTCCCATGACACAAC", "n_mismatches": 0, "panel": "Test"})
            writer.writerow({"strain_id": s, "primer_role": "FIP", "position": 0, "strand": "+", "site_seq": "GAAGAAGCTGTGCAGCCTG", "n_mismatches": 0, "panel": "Test"})
            writer.writerow({"strain_id": s, "primer_role": "BIP", "position": 0, "strand": "+", "site_seq": "CTAGTCTGCTACACCGTGC", "n_mismatches": 0, "panel": "Test"})
            
            lf_seq = "CCTTGGACGGGGAA" if s == "S4_loop_fail" else "CCTTGGACGGGGCT"
            lf_mm = 2 if s == "S4_loop_fail" else 0
            writer.writerow({"strain_id": s, "primer_role": "LF", "position": 0, "strand": "+", "site_seq": lf_seq, "n_mismatches": lf_mm, "panel": "Test"})
            
    primers_file = tmp_path / "primers.fasta"
    primers_file.write_text(">Test_F3\nGCCACCTTAAGCCACAGTA\n>Test_B3\nGTTGTGTCATGGGAGGG\n>Test_F1\nGAAGTCAGGCCCAAAAGCCA\n>Test_F2\nGAAGAAGCTGTGCAGCCTG\n>Test_B1\nTCCCCACGACGGAGCTACAG\n>Test_B2\nCTAGTCTGCTACACCGTGC\n>Test_LF\nCCTTGGACGGGGCT\n")
    
    out_html = tmp_path / "report.html"
    
    runner = CliRunner()
    result = runner.invoke(app, [
        "coverage",
        "-p", str(primers_file),
        "-f", str(fasta_file),
        "-s", str(csv_file),
        "-o", str(out_html),
        "--temperature", "63.0",
        "--dg-threshold", "-3.0",
        "--max-mismatches", "2"
    ])
    
    assert result.exit_code == 0, f"Error: {result.stdout}"
    
    html_content = out_html.read_text()
    assert "S1_parfaite" in html_content
    assert "S2_mism_interne" in html_content
    assert "S3_veto_3p" in html_content
    assert "S5_count_lost" in html_content
    assert "S4_loop_fail" in html_content
    
    assert "Divergences Comptage vs Thermodynamique" in html_content
    assert "S3_veto_3p" in html_content
    assert "S5_count_lost" in html_content


def test_coverage_cli_auto_seeding_and_export(tmp_path: Path):
    """Test du criblage intégré (sans -s) et vérification de la ré-exécutabilité via --export-sites."""
    fasta_file = tmp_path / "strains.fasta"
    # Séquence contenant F3, B3(RC), F2, B2(RC), LF
    seq_s1 = (
        "GCCACCTTAAGCCACAGTA" + "A"*30 +
        "GAAGAAGCTGTGCAGCCTG" + "A"*30 +
        "CCTTGGACGGGGCT" + "A"*30 +
        "GCACGGTGTAGCAGACTAG" + "A"*30 +
        "CCCTCCCATGACACAAC"
    )
    fasta_file.write_text(f">Strain_1\n{seq_s1}\n")
    
    primers_file = tmp_path / "primers.fasta"
    primers_file.write_text(
        ">Test_F3\nGCCACCTTAAGCCACAGTA\n"
        ">Test_B3\nGTTGTGTCATGGGAGGG\n"
        ">Test_F1\nGAAGTCAGGCCCAAAAGCCA\n"
        ">Test_F2\nGAAGAAGCTGTGCAGCCTG\n"
        ">Test_B1\nTCCCCACGACGGAGCTACAG\n"
        ">Test_B2\nCTAGTCTGCTACACCGTGC\n"
        ">Test_LF\nCCTTGGACGGGGCT\n"
    )
    
    out_html1 = tmp_path / "report1.html"
    exported_csv = tmp_path / "exported_sites.csv"
    
    runner = CliRunner()
    # 1. Exécution sans -s avec --export-sites
    res1 = runner.invoke(app, [
        "coverage",
        "-p", str(primers_file),
        "-f", str(fasta_file),
        "-o", str(out_html1),
        "--export-sites", str(exported_csv),
        "--temperature", "63.0",
        "--dg-threshold", "-3.0"
    ])
    assert res1.exit_code == 0, f"Error res1: {res1.stdout}"
    assert exported_csv.exists()
    
    # 2. Exécution avec le CSV exporté
    out_html2 = tmp_path / "report2.html"
    res2 = runner.invoke(app, [
        "coverage",
        "-p", str(primers_file),
        "-f", str(fasta_file),
        "-s", str(exported_csv),
        "-o", str(out_html2),
        "--temperature", "63.0",
        "--dg-threshold", "-3.0"
    ])
    assert res2.exit_code == 0, f"Error res2: {res2.stdout}"
    
    # Comparaison des résultats JSON
    json1 = json.loads((tmp_path / "report1.json").read_text())
    json2 = json.loads((tmp_path / "report2.json").read_text())
    assert json1["covered_thermo"] == json2["covered_thermo"]
    assert json1["covered_count"] == json2["covered_count"]
    assert json1["strains"] == json2["strains"]


def test_coverage_real_dengue_single_command(tmp_path: Path):
    """Test sur le jeu de données réel DEN-3 (amorces primer_dengue_3.txt et génome DEN3_M93130.fasta)."""
    fixture_dir = Path(__file__).resolve().parent.parent / "fixtures" / "real_primer_files"
    primers_path = fixture_dir / "primer_dengue_3.txt"
    genome_path = fixture_dir / "DEN3_M93130.fasta"
    
    out_html = tmp_path / "coverage_den3.html"
    
    runner = CliRunner()
    result = runner.invoke(app, [
        "coverage",
        "-p", str(primers_path),
        "-f", str(genome_path),
        "--panel", "3",
        "-o", str(out_html),
        "--temperature", "63.0",
        "--dg-threshold", "-5.0"
    ])
    
    assert result.exit_code == 0, f"Error: {result.stdout}"
    
    json_data = json.loads((tmp_path / "coverage_den3.json").read_text())
    assert json_data["total_strains"] == 1
    assert json_data["covered_thermo"] == 1
    assert json_data["covered_count"] == 1
    
    strain_evals = json_data["strains"][0]["evals"]
    assert strain_evals["F3"]["verdict"] == "parfait"
    assert strain_evals["F3"]["mismatches"] == 0


def test_main_module_executable():
    """Vérifie que le module labcraft.cli.main est exécutable directement avec python -m."""
    res = subprocess.run(
        [sys.executable, "-m", "labcraft.cli.main", "--help"],
        capture_output=True,
        text=True
    )
    assert res.returncode == 0
    assert "Usage: " in res.stdout or "Options" in res.stdout


def test_mg_dntp_warning_ratio():
    """Vérifie qu'un avertissement est émis dès que le rapport Mg/dNTP < 3.0."""
    with pytest.warns(UserWarning, match="Rapport Mg/dNTP faible"):
        PanelConfig(
            experiment=ExperimentConfig(
                buffer=BufferConfig(mg_mM=2.8, dntp_mM=1.4)  # ratio = 2.0 < 3.0
            )
        )


def test_coverage_seeding_performance(tmp_path: Path):
    """Test de performance : 200 souches x 6 amorces criblées et analysées en ~1 seconde."""
    # Créer un multi-FASTA de 200 souches
    strains_lines = []
    base_seq = (
        "GCCACCTTAAGCCACAGTA" + "A"*50 +
        "GAAGAAGCTGTGCAGCCTG" + "A"*50 +
        "CCTTGGACGGGGCT" + "A"*50 +
        "GCACGGTGTAGCAGACTAG" + "A"*50 +
        "CCCTCCCATGACACAAC" + "A"*100
    )
    for i in range(200):
        strains_lines.append(f">Strain_{i}\n{base_seq}\n")
        
    fasta_file = tmp_path / "200_strains.fasta"
    fasta_file.write_text("".join(strains_lines))
    
    primers_file = tmp_path / "primers.fasta"
    primers_file.write_text(
        ">Test_F3\nGCCACCTTAAGCCACAGTA\n"
        ">Test_B3\nGTTGTGTCATGGGAGGG\n"
        ">Test_F1\nGAAGTCAGGCCCAAAAGCCA\n"
        ">Test_F2\nGAAGAAGCTGTGCAGCCTG\n"
        ">Test_B1\nTCCCCACGACGGAGCTACAG\n"
        ">Test_B2\nCTAGTCTGCTACACCGTGC\n"
        ">Test_LF\nCCTTGGACGGGGCT\n"
    )
    
    out_html = tmp_path / "perf_report.html"
    
    runner = CliRunner()
    t0 = time.time()
    res = runner.invoke(app, [
        "coverage",
        "-p", str(primers_file),
        "-f", str(fasta_file),
        "-o", str(out_html),
        "--temperature", "63.0"
    ])
    t1 = time.time()
    
    assert res.exit_code == 0, f"Error: {res.stdout}"
    elapsed = t1 - t0
    print(f"\n[PERFORMANCE] 200 souches analysées en {elapsed:.3f} s.")
    assert elapsed < 5.0  # Large marge, typiquement < 1.0 s
