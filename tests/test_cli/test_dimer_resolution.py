import pytest
import os
import tempfile
from typer.testing import CliRunner
from labcraft.cli.main import app
import yaml

def test_dimer_resolution_stoichiometry():
    """Test que les dimères sont bien résolus par stoechiométrie même avec beaucoup d'underscores."""
    runner = CliRunner()
    
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as f_fasta:
        f_fasta.write(">dummy\nATGC\n")
        dummy_fasta = f_fasta.name
    
    # Créer un panel temporaire
    panel_data = {
        "experiment": {
            "name": "Test",
            "chemistry": "LAMP",
            "temperature_C": 65.0,
            "enzyme": "bst2.0",
        },
        "targets": [
            {
                "id": "Target_With_Many_Underscores",
                "sequence_file": dummy_fasta
            }
        ],
        "primer_sets": [
            {
                "target": "Target_With_Many_Underscores",
                "primers": {
                    "FIP": {"seq": "ATGCATGCATGCATGC", "conc_uM": 1.6},
                    "BIP": {"seq": "GCATGCATGCATGCAT", "conc_uM": 1.6}
                }
            }
        ]
    }
    
    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as f:
        yaml.dump(panel_data, f)
        config_path = f.name
        
    out_path = config_path.replace(".yaml", "_out.html")
    
    try:
        result = runner.invoke(app, [config_path, "-o", out_path])
        assert result.exit_code == 0, f"Erreur CLI: {result.stdout}"
        
        # Verify the report was created and has no errors
        assert os.path.exists(out_path)
        with open(out_path, "r") as f:
            html = f.read()
        
        assert "StopIteration" not in result.stdout
    finally:
        if os.path.exists(config_path): os.remove(config_path)
        if os.path.exists(out_path): os.remove(out_path)
        if os.path.exists(dummy_fasta): os.remove(dummy_fasta)

def test_no_buffer_config():
    """Test qu'une config sans bloc buffer ne plante pas (UnboundLocalError)."""
    runner = CliRunner()
    
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as f_fasta:
        f_fasta.write(">dummy\nATGC\n")
        dummy_fasta = f_fasta.name
        
    panel_data = {
        "experiment": {
            "name": "Test No Buffer",
            "chemistry": "LAMP",
            "temperature_C": 65.0,
            "enzyme": "bst2.0"
        },
        "targets": [
            {
                "id": "Target",
                "sequence_file": dummy_fasta
            }
        ],
        "primer_sets": [
            {
                "target": "Target",
                "primers": {
                    "FIP": {"seq": "ATGCATGCATGCATGC", "conc_uM": 1.6},
                    "BIP": {"seq": "GCATGCATGCATGCAT", "conc_uM": 1.6}
                }
            }
        ]
    }
    import yaml
    from labcraft.cli.main import app
    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as f2:
        yaml.dump(panel_data, f2)
        config_path = f2.name
    
    out_path = config_path.replace(".yaml", "_out.html")
    try:
        result = runner.invoke(app, [config_path, "-o", out_path])
        assert result.exit_code == 0, f"Erreur CLI: {result.stdout}"
        assert os.path.exists(out_path)
    finally:
        if os.path.exists(config_path): os.remove(config_path)
        if os.path.exists(out_path): os.remove(out_path)
        if os.path.exists(dummy_fasta): os.remove(dummy_fasta)
