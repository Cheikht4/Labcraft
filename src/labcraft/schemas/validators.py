"""Validation utilities for LabCraft schemas.

Utilitaires de validation pour les schémas LabCraft.
"""
from __future__ import annotations

import re

# Regex for DNA IUPAC
# Regex pour l'ADN IUPAC
DNA_IUPAC_REGEX = re.compile(r"^[ACGTRYSWKMBDHVN]+$", re.IGNORECASE)

def is_valid_dna_iupac(seq: str) -> bool:
    """Check if sequence contains only valid DNA IUPAC characters.
    Vérifie si la séquence ne contient que des caractères ADN IUPAC valides.
    """
    return bool(DNA_IUPAC_REGEX.match(seq))
