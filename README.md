# LabCraft
**LabCraft: building reliable multiplex amplification assays, one thermodynamic brick at a time.**

LabCraft is an open-source thermodynamic simulation engine designed for multiplex primer panel design, supporting both LAMP (Loop-Mediated Isothermal Amplification) and PCR (Polymerase Chain Reaction) assays.

![CI](https://img.shields.io/badge/CI-passing-success)
![PyPI](https://img.shields.io/badge/PyPI-v0.0.1-blue)
![License](https://img.shields.io/badge/License-MIT-green)

## Features
- Thermodynamic calculations for various duplexes
- Advanced diagnostics and risk evaluation
- Scalable pipeline for multiplex assay design

## Installation
Currently in pre-alpha. To install in development mode:
```bash
git clone https://github.com/Cheikht4/Labcraft.git
cd Labcraft
pip install -e .
```

## Usage

You can launch an analysis using the command line interface:
```bash
labcraft examples/denv2_lopez.yaml -o report.html
```

You can also run validation scripts:
```bash
python validation/validate_parida.py
python validation/validate_dengue.py
```

### Configuration Schema (YAML)

LabCraft uses a canonical configuration schema:
```yaml
experiment:
  name: "DENV2 LAMP 6-plex"
  chemistry: LAMP
  temperature_C: 65.0
  enzyme: bst2.0
  buffer: # optional
    na_mM: 50
    k_mM: 10
    tris_mM: 20
    mg_mM: 6
    dntp_mM: 1.4

targets:
  - id: DENV2
    sequence_file: path/to/target.fasta
    copies_per_uL: 1000

primer_sets:
  - target: DENV2
    primers:
      F3:  { seq: "...", conc_uM: 0.05 }
      B3:  { seq: "...", conc_uM: 0.05 }
      FIP: { seq: "...", conc_uM: 0.4, domains: { F2: "...", F1c: "...", linker: "" } }
      BIP: { seq: "...", conc_uM: 0.4 }
      LF:  { seq: "...", conc_uM: 0.2 }
      LB:  { seq: "...", conc_uM: 0.2 }
```

## Quick Start
```python
import labcraft

print(f"LabCraft version: {labcraft.__version__}")
```

## Documentation
Documentation will be available at [https://github.com/Cheikht4/Labcraft](https://github.com/Cheikht4/Labcraft).

## Citation
If you use LabCraft, please cite our forthcoming publication.

## Availability & Requirements
- **Python:** ≥ 3.10
- **License:** MIT License
- **Note:** `primer3-py` (GPLv2) is a runtime dependency.
