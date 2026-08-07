import pytest
from labcraft.cli.config import PanelConfig, build_engine_from_config
from labcraft.diagnostics.enzyme import get_enzyme
import yaml

def test_enzyme_tolerance():
    assert get_enzyme("bst 2.0").name == "Bst 2.0"
    assert get_enzyme("BST2.0").name == "Bst 2.0"
    assert get_enzyme("Bst2.0_WarmStart").name == "Bst 2.0 WarmStart"
    assert get_enzyme("bst 2.0 warmstart").name == "Bst 2.0 WarmStart"
    assert get_enzyme("taq").name == "Taq"

    with pytest.raises(ValueError):
        get_enzyme("unknown_enzyme")

def test_config_buffer_parsing():
    yaml_str = """
experiment:
  name: "Test"
  chemistry: LAMP
  temperature_C: 65.0
  enzyme: bst2.0
  buffer:
    na_mM: 50.0
    mg_mM: 8.0
targets: []
primer_sets: []
    """
    data = yaml.safe_load(yaml_str)
    config = PanelConfig.model_validate(data)
    
    assert config.experiment.buffer.na_mM == 50.0
    assert config.experiment.buffer.mg_mM == 8.0
    
    _, _, backend, backend_kwargs, _, _, _, _ = build_engine_from_config(config, {})
    assert backend_kwargs['mg_mm'] == 8.0

