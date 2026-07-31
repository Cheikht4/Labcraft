"""Tests for experiment schemas.

Tests pour les schémas d'expérience.
"""
from __future__ import annotations

import warnings
from pathlib import Path

import pytest
from pydantic import ValidationError

from labcraft.schemas.experiment import LabCraftInput


VALID_LAMP_YAML = """
experiment:
  name: "Panel Dengue / Fievre Jaune"
  chemistry: LAMP
  temperature_C: 65.0
  enzyme: Bst2.0_WarmStart

buffer:
  Na_mM: 50
  K_mM: 10
  Tris_mM: 20
  Mg_total_mM: 8.0
  dNTP_total_mM: 1.4
  EDTA_mM: 0.0
  chelation_model: simple

targets:
  - id: DENV
    sequence_file: denv.fasta
    type: RNA
    copies_per_uL: 1000
  - id: YFV
    sequence_file: yfv.fasta
    type: RNA
    copies_per_uL: 1000

primer_sets:
  - target: DENV
    primers:
      F3:  {seq: "ACCTGG", conc_uM: 0.2}
      B3:  {seq: "TTGCAA", conc_uM: 0.2}
      FIP: {seq: "GCTTTTTTCAGG", conc_uM: 1.6, domains: {F1c: "GCTT", linker: "TTTT", F2: "CAGG"}}
      BIP: {seq: "ACG", conc_uM: 1.6, domains: auto}
      LF:  {seq: "TCA", conc_uM: 0.8}
      LB:  {seq: "GAC", conc_uM: 0.8}
  - target: YFV
    primers:
      F3:  {seq: "A", conc_uM: 0.2}
      B3:  {seq: "T", conc_uM: 0.2}
      FIP: {seq: "G", conc_uM: 1.6, domains: auto}
      BIP: {seq: "C", conc_uM: 1.6, domains: auto}

analysis:
  temperature_scan: [50, 70, 0.5]
  specificity_db: human_transcriptome.fasta
  amplifiable_dimer_threshold_kcal: -2.0
  duplex_backend: primer3
  max_complex_size: 2
"""

VALID_PCR_YAML = """
experiment:
  name: "Simple PCR"
  chemistry: PCR
  temperature_C: 95.0
  enzyme: Taq

buffer:
  Mg_total_mM: 2.0
  dNTP_total_mM: 0.2

targets:
  - id: target1
    sequence_file: test.fasta

primer_sets:
  - target: target1
    primers:
      Forward: {seq: "ACGT", conc_uM: 0.5}
      Reverse: {seq: "TGCA", conc_uM: 0.5}
"""

def test_load_valid_pcr():
    """Test loading minimal valid PCR yaml.
    Test du chargement d'un YAML PCR valide minimal.
    """
    import yaml
    data = yaml.safe_load(VALID_PCR_YAML)
    model = LabCraftInput.model_validate(data)
    
    assert model.experiment.name == "Simple PCR"
    assert model.experiment.chemistry == "PCR"
    assert model.experiment.enzyme == "Taq"
    assert len(model.primer_sets[0].primers) == 2


def test_load_valid_lamp():
    """Test loading valid full LAMP yaml.
    Test du chargement d'un YAML LAMP complet valide.
    """
    import yaml
    data = yaml.safe_load(VALID_LAMP_YAML)
    model = LabCraftInput.model_validate(data)
    
    assert model.experiment.name == "Panel Dengue / Fievre Jaune"
    assert len(model.targets) == 2
    assert model.analysis.max_complex_size == 2
    assert model.buffer.Mg_total_mM == 8.0


def test_target_not_found():
    """Test error when target referenced in primer_set is not found.
    Test d'erreur si target non trouvé dans primer_set.target.
    """
    import yaml
    data = yaml.safe_load(VALID_PCR_YAML)
    data["primer_sets"][0]["target"] = "unknown_target"
    
    with pytest.raises(ValidationError, match="Target 'unknown_target' referenced by primer set not found"):
        LabCraftInput.model_validate(data)


def test_invalid_dna_sequence():
    """Test error when sequence contains invalid characters.
    Test d'erreur si séquence contient des caractères invalides.
    """
    import yaml
    data = yaml.safe_load(VALID_PCR_YAML)
    data["primer_sets"][0]["primers"]["Forward"]["seq"] = "ACGTZ"
    
    with pytest.raises(ValidationError, match="contains invalid DNA IUPAC characters"):
        LabCraftInput.model_validate(data)


def test_lamp_without_fip_bip():
    """Test error if LAMP is missing FIP/BIP.
    Test d'erreur si LAMP sans FIP/BIP.
    """
    import yaml
    data = yaml.safe_load(VALID_LAMP_YAML)
    del data["primer_sets"][0]["primers"]["FIP"]
    
    with pytest.raises(ValidationError, match="For LAMP, primer set must contain at least"):
        LabCraftInput.model_validate(data)


def test_negative_temperature():
    """Test error on negative temperature.
    Test d'erreur si température négative.
    """
    import yaml
    data = yaml.safe_load(VALID_PCR_YAML)
    data["experiment"]["temperature_C"] = -5.0
    
    with pytest.raises(ValidationError, match="greater than 0"):
        LabCraftInput.model_validate(data)


def test_default_values():
    """Test default values are populated.
    Test que les valeurs par défaut sont correctes.
    """
    import yaml
    data = yaml.safe_load(VALID_PCR_YAML)
    model = LabCraftInput.model_validate(data)
    
    assert model.buffer.Na_mM == 50.0
    assert model.targets[0].type == "DNA"
    assert model.targets[0].copies_per_uL == 1000.0
    assert model.analysis.amplifiable_dimer_threshold_kcal == -2.0
    assert model.analysis.duplex_backend == "primer3"
    assert model.analysis.max_complex_size == 2


def test_from_yaml(tmp_path: Path):
    """Test loading from a temporary YAML file.
    Test du from_yaml avec un fichier temporaire.
    """
    yaml_file = tmp_path / "test.yaml"
    yaml_file.write_text(VALID_PCR_YAML, encoding="utf-8")
    
    model = LabCraftInput.from_yaml(yaml_file)
    assert model.experiment.chemistry == "PCR"


def test_duplex_backend_and_complex_size():
    """Test validation of duplex_backend and max_complex_size.
    Test de la validation du duplex_backend et max_complex_size.
    """
    import yaml
    data = yaml.safe_load(VALID_PCR_YAML)
    data["analysis"] = {
        "duplex_backend": "invalid_backend",
        "max_complex_size": 5
    }
    
    with pytest.raises(ValidationError) as exc_info:
        LabCraftInput.model_validate(data)
        
    errors = str(exc_info.value)
    assert "Input should be 'native', 'primer3' or 'vienna'" in errors
    assert "Input should be less than or equal to 4" in errors
