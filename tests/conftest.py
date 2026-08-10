"""Shared test fixtures / Fixtures de test partagées."""
import pytest

# NUPACK availability check / Vérification de la disponibilité de NUPACK
try:
    import nupack  # noqa: F401
    HAS_NUPACK = True
except ImportError:
    HAS_NUPACK = False

requires_nupack = pytest.mark.skipif(
    not HAS_NUPACK,
    reason="NUPACK not installed (optional validation dependency)"
)

import sys
import pathlib
project_root = pathlib.Path(__file__).parents[1]
sys.path.insert(0, str(project_root))
