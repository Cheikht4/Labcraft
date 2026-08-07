"""FIP/BIP domain decomposition / Décomposition en domaines fonctionnels.

Modélisation des amorces physiques et extraction de leurs domaines de liaison.
"""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass
from enum import Enum


class PrimerRole(Enum):
    F3 = "F3"
    B3 = "B3"
    FIP = "FIP"
    BIP = "BIP"
    LF = "LF"
    LB = "LB"
    # PCR
    FWD = "FWD"
    REV = "REV"
    PROBE = "PROBE"


IUPAC_DNA = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
}

def _iupac_to_regex(seq: str) -> str:
    return "".join(IUPAC_DNA.get(c.upper(), c) for c in seq)

def _revcomp(seq: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                  'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
                  'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                  'D': 'H', 'H': 'D', 'N': 'N'}
    return "".join(complement.get(c.upper(), 'N') for c in reversed(seq))


@dataclass(frozen=True)
class PhysicalPrimer:
    """Représente un oligonucléotide physique dans le tube.
    
    Pour FIP/BIP, l'amorce est composite (5'-c1-linker-2-3').
    Le domaine 3' (binding_domain) s'hybride à la cible initiale.
    Le domaine 5' (tail_domain) sert à la formation des boucles.
    """
    name: str
    sequence: str
    role: PrimerRole
    binding_domain: str
    tail_domain: str | None = None
    linker: str | None = None
    
    @classmethod
    def from_simple(cls, name: str, sequence: str, role: PrimerRole) -> PhysicalPrimer:
        """Crée une amorce simple (ex: F3, B3, LF, LB)."""
        return cls(
            name=name, sequence=sequence, role=role,
            binding_domain=sequence, tail_domain=None, linker=None
        )

    @classmethod
    def from_domains(
        cls, name: str, role: PrimerRole, tail: str, f2_b2: str, linker: str = ""
    ) -> PhysicalPrimer:
        """Déclaration manuelle explicite des domaines."""
        seq = tail + linker + f2_b2
        return cls(
            name=name, sequence=seq, role=role,
            binding_domain=f2_b2, tail_domain=tail, linker=linker if linker else None
        )

    @classmethod
    def from_alignment(
        cls, name: str, sequence: str, role: PrimerRole, target_seq: str,
        allow_mismatches: bool = False
    ) -> PhysicalPrimer:
        """Autodétection des domaines par alignement sur la cible.
        
        Tente de trouver la coupure (tail / binding) en cherchant le domaine 3'
        (F2/B2) et le domaine 5' (F1c/B1c) sur le brin cible ou son reverse-complement.
        
        Args:
            name: Nom de l'amorce.
            sequence: Séquence complète de l'amorce.
            role: FIP ou BIP.
            target_seq: Séquence de la cible (un seul brin de référence).
            allow_mismatches: Si True, autorise l'alignement imparfait (non implémenté par défaut).
        """
        if allow_mismatches:
            warnings.warn("L'alignement avec mismatches n'est pas supporté par défaut.", UserWarning)
            
        sequence = sequence.upper()
        target_seq = target_seq.upper()
        target_rc = _revcomp(target_seq)
        
        # Pour une FIP/BIP, on sait que le domaine de liaison est à l'extrémité 3'.
        # On va chercher le plus long suffixe de l'amorce qui s'aligne (sens ou RC) sur la cible.
        # Longueur minimale d'un domaine : typiquement 12-15 bp.
        min_domain_len = 12
        best_binding_len = 0
        
        # Chercher le suffixe (binding domain)
        for i in range(len(sequence) - min_domain_len, 0, -1):
            suffix = sequence[i:]
            regex = _iupac_to_regex(suffix)
            # F2 s'hybride sur le brin (+), donc F2 est le RC d'un bout de (+).
            # Ou B2 s'hybride sur le brin (-), donc B2 est de même sens que (+).
            # On cherche donc le suffixe tel quel ou son RC dans la cible.
            if re.search(regex, target_seq) or re.search(regex, target_rc):
                # On veut le suffixe maximum
                if len(suffix) > best_binding_len:
                    best_binding_len = len(suffix)
                    
        if best_binding_len == 0:
            if target_seq: # Only warn if we genuinely couldn't find it when target is provided
                warnings.warn(f"Le domaine de liaison de {name} n'est pas trouvé sur la cible.", UserWarning)
            return cls.from_simple(name, sequence, role)
            
        # Chercher le préfixe (tail domain)
        binding_idx = len(sequence) - best_binding_len
        best_tail_len = 0
        
        for i in range(min_domain_len, binding_idx + 1):
            prefix = sequence[:i]
            regex = _iupac_to_regex(prefix)
            if re.search(regex, target_seq) or re.search(regex, target_rc):
                best_tail_len = i
                
        if best_tail_len == 0:
            if target_seq:
                warnings.warn(f"Impossible d'autodétecter le domaine 5' pour {name}.", UserWarning)
            # On fallback sur une séparation brutale (tout ce qui n'est pas binding est tail)
            tail = sequence[:binding_idx]
            binding = sequence[binding_idx:]
            return cls(name=name, sequence=sequence, role=role, binding_domain=binding, tail_domain=tail, linker=None)
            
        tail = sequence[:best_tail_len]
        binding = sequence[binding_idx:]
        linker = sequence[best_tail_len:binding_idx]
        
        # Validation d'ambiguïté (si l'amorce est quasi homopolymérique etc.)
        # Normalement l'expression régulière fait un match simple.
        
        return cls(
            name=name, sequence=sequence, role=role,
            binding_domain=binding, tail_domain=tail, linker=linker if linker else None
        )
