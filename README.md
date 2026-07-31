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
git clone https://github.com/labcraft-dev/labcraft.git
cd labcraft
pip install -e .
```

## Quick Start
```python
import labcraft

print(f"LabCraft version: {labcraft.__version__}")
```

## Documentation
Documentation will be available at [https://labcraft-dev.github.io/labcraft](https://labcraft-dev.github.io/labcraft).

## Citation
If you use LabCraft, please cite our forthcoming publication.

## Availability & Requirements
- **Python:** ≥ 3.10
- **License:** MIT License
- **Note:** `primer3-py` (GPLv2) is a runtime dependency.
