import pytest
import os
import tempfile
from typer.testing import CliRunner
from labcraft.cli.main import app
import yaml

panel_data = {
    "experiment": {
        "name": "Test",
        "chemistry": "LAMP",
        "temperature_C": 65.0,
        "enzyme": "bst2.0",
        "buffer": {
            "na_mM": 50.0,
            "k_mM": 10.0,
            "tris_mM": 20.0,
            "mg_mM": 6.0,
            "dntp_mM": 1.4
        }
    },
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

runner = CliRunner()
result = runner.invoke(app, ["analyze", config_path, "-o", out_path])
print(result.stdout)
if result.exception:
    import traceback
    traceback.print_exception(type(result.exception), result.exception, result.exception.__traceback__)
