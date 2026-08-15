import pytest
import yaml
import os

def test_multiplex_integration(tmp_path):
    target1_fasta = tmp_path / "t1.fasta"
    target1_fasta.write_text(">T1\nATCGATCGATCGATCGATCGATCGATCGATCG\n")
    
    target2_fasta = tmp_path / "t2.fasta"
    target2_fasta.write_text(">T2\nGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n")
    
    config_yaml = tmp_path / "multiplex.yaml"
    config_yaml.write_text(f"""
experiment:
  name: "Multiplex Test"
  chemistry: "LAMP"
  temperature_C: 65.0
  enzyme: "bst2.0"
  buffer:
    na_mM: 50.0
    mg_mM: 8.0
    dntp_mM: 1.4

targets:
  - id: "T1"
    sequence_file: "{target1_fasta}"
    copies_per_uL: 100.0
  - id: "T2"
    sequence_file: "{target2_fasta}"
    copies_per_uL: 100.0

primer_sets:
  - target: "T1"
    primers:
      F3:
        seq: "ATCGATCGATCGATCG"
        conc_uM: 0.2
        domains:
          F3: "ATCGATCGATCGATCG"
  - target: "T2"
    primers:
      F3:
        seq: "GCTAGCTAGCTAGCTA"
        conc_uM: 0.2
        domains:
          F3: "GCTAGCTAGCTAGCTA"
""")

    # Try to load and run the whole main analysis path logic for mispriming
    from labcraft.cli.main import app
    from typer.testing import CliRunner
    # Just running analyze should work and not crash
    report_file = tmp_path / "report.html"
    runner = CliRunner()
    result = runner.invoke(app, ["analyze", "-c", str(config_yaml), "-o", str(report_file)])
    assert result.exit_code == 0
    
    assert report_file.exists()
