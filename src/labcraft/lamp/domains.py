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

IUPAC_EXPANSION = {
    'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
    'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'C', 'G', 'T']
}

def expand_degenerate(seq: str, max_variants: int = 16) -> list[str]:
    """
    Développe une séquence avec des codes IUPAC dégénérés en toutes ses séquences concrètes.
    Lève une ValueError si le nombre de variants dépasse max_variants.
    """
    import itertools
    import math
    
    seq = seq.upper().replace('U', 'T')
    options = [IUPAC_EXPANSION.get(c, [c]) for c in seq]
    
    total_variants = math.prod(len(opts) for opts in options)
    if total_variants > max_variants:
        raise ValueError(f"La séquence {seq} générerait {total_variants} variants (limite: {max_variants}).")
        
    return ["".join(p) for p in itertools.product(*options)]


def _match_iupac_substring(query: str, target: str) -> bool:
    """
    Vérifie si la séquence query est présente dans target avec tolérance IUPAC bilatérale.
    Une position correspond si l'intersection des bases possibles des deux côtés n'est pas vide.
    """
    return _find_iupac_substring(query, target) != -1

def _find_iupac_substring(query: str, target: str) -> int:
    """
    Retourne l'indice de la première occurrence de query dans target avec tolérance IUPAC bilatérale,
    ou -1 si non trouvé.
    """
    q_len = len(query)
    t_len = len(target)
    if q_len == 0 or q_len > t_len:
        return -1
        
    q_sets = [set(IUPAC_EXPANSION.get(c.upper(), [c.upper()])) for c in query]
    t_sets = [set(IUPAC_EXPANSION.get(c.upper(), [c.upper()])) for c in target]
    
    for i in range(t_len - q_len + 1):
        match = True
        for j in range(q_len):
            if not q_sets[j].intersection(t_sets[i+j]):
                match = False
                break
        if match:
            return i
            
    return -1
def _revcomp(seq: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                  'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
                  'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                  'D': 'H', 'H': 'D', 'N': 'N'}
    return "".join(complement.get(c.upper(), 'N') for c in reversed(seq.upper().replace('U', 'T')))


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
    nominal_concentration: float | None = None
    parent_name: str | None = None
    blocked_3prime: bool = False
    label_5prime: str | None = None
    lna_positions: tuple[int, ...] = ()
    
    @classmethod
    def from_simple(
        cls, name: str, sequence: str, role: PrimerRole, 
        nominal_concentration: float | None = None, parent_name: str | None = None,
        blocked_3prime: bool = False, label_5prime: str | None = None
    ) -> PhysicalPrimer:
        """Crée une amorce simple (ex: F3, B3, LF, LB, PROBE)."""
        from labcraft.thermo.lna import parse_lna_sequence
        bare_seq, lna_pos = parse_lna_sequence(sequence)
        return cls(
            name=name, sequence=bare_seq, role=role,
            binding_domain=bare_seq, tail_domain=None, linker=None,
            nominal_concentration=nominal_concentration, parent_name=parent_name,
            blocked_3prime=blocked_3prime, label_5prime=label_5prime,
            lna_positions=tuple(lna_pos)
        )

    @classmethod
    def from_domains(
        cls, name: str, role: PrimerRole, tail: str, f2_b2: str, linker: str = "", 
        nominal_concentration: float | None = None, parent_name: str | None = None,
        blocked_3prime: bool = False, label_5prime: str | None = None
    ) -> PhysicalPrimer:
        """Déclaration manuelle explicite des domaines."""
        from labcraft.thermo.lna import parse_lna_sequence
        
        # Le parsing global permet de conserver l'information LNA
        seq = tail + linker + f2_b2
        bare_seq, lna_pos = parse_lna_sequence(seq)
        
        # On parse aussi les domaines individuellement pour les stocker nus
        bare_tail, _ = parse_lna_sequence(tail)
        bare_linker, _ = parse_lna_sequence(linker)
        bare_f2_b2, _ = parse_lna_sequence(f2_b2)
        
        return cls(
            name=name, sequence=bare_seq, role=role,
            binding_domain=bare_f2_b2, tail_domain=bare_tail, linker=bare_linker if bare_linker else None,
            nominal_concentration=nominal_concentration, parent_name=parent_name,
            blocked_3prime=blocked_3prime, label_5prime=label_5prime,
            lna_positions=tuple(lna_pos)
        )

    @classmethod
    def from_alignment(
        cls, name: str, sequence: str, role: PrimerRole, target_seq: str,
        allow_mismatches: bool = False, nominal_concentration: float | None = None, parent_name: str | None = None,
        blocked_3prime: bool = False, label_5prime: str | None = None
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
            
        from labcraft.thermo.lna import parse_lna_sequence
        bare_seq, lna_pos = parse_lna_sequence(sequence.upper())
        sequence = bare_seq
        
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
            if _match_iupac_substring(suffix, target_seq) or _match_iupac_substring(suffix, target_rc):
                # On veut le suffixe maximum
                if len(suffix) > best_binding_len:
                    best_binding_len = len(suffix)
                    
        if best_binding_len == 0:
            if target_seq: # Only warn if we genuinely couldn't find it when target is provided
                warnings.warn(f"Le domaine de liaison de {name} n'est pas trouvé sur la cible.", UserWarning)
            # Reconstruct original sequence to pass to from_simple
            # Wait, from_simple will parse again, so we should just construct it directly here
            return cls(
                name=name, sequence=sequence, role=role,
                binding_domain=sequence, tail_domain=None, linker=None,
                nominal_concentration=nominal_concentration, parent_name=parent_name,
                blocked_3prime=blocked_3prime, label_5prime=label_5prime,
                lna_positions=tuple(lna_pos)
            )
            
        # Chercher le préfixe (tail domain)
        binding_idx = len(sequence) - best_binding_len
        best_tail_len = 0
        
        for i in range(min_domain_len, binding_idx + 1):
            prefix = sequence[:i]
            if _match_iupac_substring(prefix, target_seq) or _match_iupac_substring(prefix, target_rc):
                best_tail_len = i
                
        if best_tail_len == 0:
            if target_seq:
                warnings.warn(f"Impossible d'autodétecter le domaine 5' pour {name}.", UserWarning)
            # On fallback sur une séparation brutale (tout ce qui n'est pas binding est tail)
            tail = sequence[:binding_idx]
            binding = sequence[binding_idx:]
            return cls(
                name=name, sequence=sequence, role=role, 
                binding_domain=binding, tail_domain=tail, linker=None, 
                nominal_concentration=nominal_concentration, parent_name=parent_name,
                blocked_3prime=blocked_3prime, label_5prime=label_5prime,
                lna_positions=tuple(lna_pos)
            )
            
        tail = sequence[:best_tail_len]
        binding = sequence[binding_idx:]
        linker = sequence[best_tail_len:binding_idx]
        
        # Validation d'ambiguïté (si l'amorce est quasi homopolymérique etc.)
        # Normalement l'expression régulière fait un match simple.
        
        return cls(
            name=name, sequence=sequence, role=role,
            binding_domain=binding, tail_domain=tail, linker=linker if linker else None,
            nominal_concentration=nominal_concentration, parent_name=parent_name,
            blocked_3prime=blocked_3prime, label_5prime=label_5prime,
            lna_positions=tuple(lna_pos)
        )
