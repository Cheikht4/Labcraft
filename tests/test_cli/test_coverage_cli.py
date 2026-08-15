import pytest
from typer.testing import CliRunner
import tempfile
import csv
from pathlib import Path

from labcraft.cli.main import app

def test_coverage_cli(tmp_path: Path):
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
                f3_mm = 3  # Count lost, but sequence is perfect so Thermo is covered
                
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
    
    if result.exit_code != 0:
        import traceback
        traceback.print_exception(type(result.exception), result.exception, result.exception.__traceback__)
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
