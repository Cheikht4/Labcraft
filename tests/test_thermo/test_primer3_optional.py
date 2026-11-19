"""Tests vérifiant que `primer3-py` (GPL-2.0) est purement optionnel (`dev`) et absent du runtime MIT.
Tests verifying `primer3-py` (GPL-2.0) is purely optional (`dev`) and absent from the MIT runtime.
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

try:
    import tomllib  # type: ignore[import-not-found]
except ModuleNotFoundError:  # Python 3.10
    import yaml as tomllib  # pragma: no cover

REPO_ROOT = Path(__file__).resolve().parent.parent.parent


def test_primer3_moved_to_dev_optional_and_absent_from_src(monkeypatch: pytest.MonkeyPatch) -> None:
    """Vérifie que `pyproject.toml` place `primer3-py` dans `[project.optional-dependencies].dev`,
    qu'aucun fichier de `src/labcraft/` n'importe `primer3`, et que les tests dédiés sont ignorés
    proprement (`pytest.skip.Exception`) lorsque `primer3` est absent.

    Verifies `pyproject.toml` puts `primer3-py` under `[project.optional-dependencies].dev`,
    no file in `src/labcraft/` imports `primer3`, and Primer3 tests skip cleanly when `primer3` is absent.
    """
    pyproject_text = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    if hasattr(tomllib, "loads") and tomllib.__name__ == "tomllib":
        data = tomllib.loads(pyproject_text)
        runtime_deps = data["project"]["dependencies"]
        dev_deps = data["project"]["optional-dependencies"]["dev"]
        assert not any("primer3" in d for d in runtime_deps)
        assert any(d.startswith("primer3-py") for d in dev_deps)
    else:
        deps_block = pyproject_text.split("[project.optional-dependencies]")[0]
        dev_block = pyproject_text.split("[project.optional-dependencies]")[1]
        assert "primer3-py" not in deps_block
        assert "primer3-py" in dev_block

    # Aucun import de `primer3` dans `src/labcraft` / Zero `primer3` imports in `src/labcraft`
    for py_path in (REPO_ROOT / "src" / "labcraft").rglob("*.py"):
        tree = ast.parse(py_path.read_text(encoding="utf-8"), filename=str(py_path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                assert all(alias.name.split(".")[0] != "primer3" for alias in node.names), str(py_path)
            elif isinstance(node, ast.ImportFrom) and node.module:
                assert node.module.split(".")[0] != "primer3", str(py_path)

    # Simulation de l'absence de `primer3` : `importorskip` déclenche proprement `Skipped`
    # Simulate missing `primer3`: `importorskip` cleanly raises `Skipped`
    from tests.test_thermo.test_step2_mismatch_thermo import (
        test_consecutive_tandem_mismatches_match_primer3_santalucia_2004,
    )
    from tests.test_thermo.test_step5_bulges_and_boltzmann_3p import (
        test_point2_santalucia_hicks_2004_matches_primer3_reference,
    )

    monkeypatch.setitem(sys.modules, "primer3", None)
    with pytest.raises(pytest.skip.Exception):
        test_consecutive_tandem_mismatches_match_primer3_santalucia_2004()
    with pytest.raises(pytest.skip.Exception):
        test_point2_santalucia_hicks_2004_matches_primer3_reference()
