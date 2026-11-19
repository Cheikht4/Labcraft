"""Vérifie que tous les chemins, plages de lignes et valeurs numériques de `CHANGELOG.md`
proviennent réellement du dépôt et de `validation/reference_benchmarks.md`.

Verifies that all file paths, line ranges, and numerical values in `CHANGELOG.md`
exist in the repository and come verbatim from `validation/reference_benchmarks.md`.
"""
from __future__ import annotations

from pathlib import Path
import re
from typing import List, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent

_BACKTICK_SPAN_RE = re.compile(r"`([^`\n]+)`")
_PATH_TOKEN_RE = re.compile(
    r"^(?P<path>"
    r"(?:(?:src|tests|scripts|validation|examples|\.github)/[A-Za-z0-9_./-]+)"
    r"|(?:pyproject\.toml|README\.md|CHANGELOG\.md|CITATION\.cff)"
    r"|(?:[A-Za-z0-9_./-]+\.(?:py|j2|toml|ya?ml|md|cff|csv|fasta|txt))"
    r")"
    r"(?::(?P<lines>[0-9]+(?:\s*[-–,]\s*[0-9]+)*))?$"
)


def _extract_changelog_citations(changelog_text: str) -> List[Tuple[str, str]]:
    """Extrait les couples `(chemin, spec_lignes)` cités entre backticks dans `CHANGELOG.md`.
    Extracts `(path, line_spec)` pairs cited inside backticks in `CHANGELOG.md`.
    """
    citations: List[Tuple[str, str]] = []
    for match in _BACKTICK_SPAN_RE.finditer(changelog_text):
        span = match.group(1).strip()
        candidates = [span]
        if " " in span:
            candidates.extend(tok.strip(",;()") for tok in span.split())
        for cand in candidates:
            m = _PATH_TOKEN_RE.match(cand)
            if m:
                citations.append((m.group("path"), m.group("lines") or ""))
    return citations


def test_all_changelog_paths_and_line_ranges_exist() -> None:
    """Vérifie que chaque chemin cité dans `CHANGELOG.md` existe sur le disque et que ses numéros
    de lignes sont compris dans les bornes réelles du fichier.

    Verifies that every path cited in `CHANGELOG.md` exists on disk and that any cited line numbers
    fall within the actual line count of the referenced file.
    """
    changelog_path = REPO_ROOT / "CHANGELOG.md"
    assert changelog_path.is_file(), "CHANGELOG.md est introuvable à la racine du dépôt."
    bench_path = REPO_ROOT / "validation" / "reference_benchmarks.md"
    assert bench_path.is_file(), "validation/reference_benchmarks.md est introuvable."

    text = changelog_path.read_text(encoding="utf-8")
    citations = _extract_changelog_citations(text)
    assert len(citations) >= 35, f"Trop peu de chemins extraits de CHANGELOG.md ({len(citations)})."

    missing_files: List[str] = []
    invalid_lines: List[str] = []

    for rel_path, line_spec in citations:
        target = REPO_ROOT / rel_path
        if not target.exists():
            missing_files.append(rel_path)
            continue
        if line_spec and target.is_file():
            n_lines = len(target.read_text(encoding="utf-8").splitlines())
            for part in [p.strip() for p in line_spec.split(",") if p.strip()]:
                norm = part.replace("–", "-")
                if "-" in norm:
                    lo_s, hi_s = [x.strip() for x in norm.split("-", 1)]
                    lo, hi = int(lo_s), int(hi_s)
                    if not (1 <= lo <= hi <= n_lines):
                        invalid_lines.append(f"{rel_path}:{line_spec} (fichier de {n_lines} lignes)")
                else:
                    ln = int(norm)
                    if not (1 <= ln <= n_lines):
                        invalid_lines.append(f"{rel_path}:{line_spec} (fichier de {n_lines} lignes)")

    assert not missing_files, f"Chemins introuvables cités dans CHANGELOG.md : {missing_files}"
    assert not invalid_lines, f"Plages de lignes hors bornes dans CHANGELOG.md : {invalid_lines}"


def test_all_changelog_decimal_figures_come_from_reference_benchmarks() -> None:
    """Vérifie qu'aucun nombre décimal (`[+-]?\\d+[.,]\\d+`) présent dans `CHANGELOG.md`
    n'est absent de `validation/reference_benchmarks.md`.

    Verifies that every decimal number (`[+-]?\\d+[.,]\\d+`) in `CHANGELOG.md` appears verbatim
    in `validation/reference_benchmarks.md`.
    """
    changelog_text = (REPO_ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    bench_text = (REPO_ROOT / "validation" / "reference_benchmarks.md").read_text(encoding="utf-8")

    # Retirer les hash de commits hexadécimaux de 7 caractères entre backticks avant extraction
    # Strip 7-char hex commit hashes in backticks before scanning decimal numbers
    cleaned = re.sub(r"`[0-9a-f]{7}`", "", changelog_text)
    decimals = sorted(set(re.findall(r"[+-]?\d+[.,]\d+", cleaned)))
    assert len(decimals) >= 20
    missing_decimals = [val for val in decimals if val not in bench_text]
    assert not missing_decimals, (
        f"Valeurs décimales de CHANGELOG.md absentes de validation/reference_benchmarks.md : {missing_decimals}"
    )
