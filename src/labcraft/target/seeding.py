"""Module de criblage permissif de sites candidats pour le mode couverture.
Permissive candidate site seeding module for coverage analysis.

Origine / Origin:
Adapté de 'lamp_coverage.py' (dépôt primer-analysis-suite, même auteur).
Rôle : criblage permissif produisant des sites CANDIDATS, re-notés ensuite
par la couche thermodynamique de LabCraft.
Role: permissive screening generating CANDIDATE sites, which are then re-scored
by the thermodynamic layer of LabCraft.

Modèle de criblage / Screening model:
LabCraft crible sur SUBSTITUTIONS SEULES ({s<=n}), par cohérence avec le modèle
thermodynamique NN qui suppose un alignement sans brèche. La recherche tolérante
aux indels reste dans l'outil d'origine (lamp_coverage.py).
LabCraft screens on SUBSTITUTIONS ONLY ({s<=n}), in consistency with the NN
thermodynamic model which assumes gapless alignment. Indel-tolerant search
remains in the original tool (lamp_coverage.py).
"""

from __future__ import annotations
from typing import Dict, List, Optional, Tuple, Any
import regex

from labcraft.lamp.domains import PhysicalPrimer, PrimerRole, _revcomp, IUPAC_MATCHABLE

# Dictionnaire IUPAC vers classes d'expression régulière
# IUPAC dictionary to regular expression classes
IUPAC_DICT = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
}


def seq_to_regex(seq: str) -> str:
    """Convertit une séquence avec codes IUPAC en motif d'expression régulière.
    Converts an IUPAC sequence into a regular expression pattern.
    """
    pattern = ""
    for char in seq.upper():
        if char.isalpha():
            pattern += IUPAC_DICT.get(char, char)
    return pattern


def build_primer_regex(
    primer_seq: str,
    max_errors: int,
    strict_3prime_len: int = 3,
    strict_3prime_tolerate: int = 0,
    is_rc: bool = False
) -> str:
    """Construit une expression régulière tolérante en 5' et stricte en 3' (substitutions seules).
    Builds a regex pattern tolerant in 5' and strict in 3' (substitutions only).
    """
    tolerate_positions = set()
    if strict_3prime_tolerate == 1:
        tolerate_positions = {2}
    elif strict_3prime_tolerate == 2:
        tolerate_positions = {1, 2}

    if strict_3prime_len > 0 and len(primer_seq) > strict_3prime_len:
        if not is_rc:
            # Brin sens / Sense strand
            seq_5 = primer_seq[:-strict_3prime_len]
            pattern_5 = seq_to_regex(seq_5)
            pattern_3 = ""
            for k in range(strict_3prime_len, 0, -1):
                base = primer_seq[-k]
                base_regex = seq_to_regex(base)
                if k in tolerate_positions:
                    pattern_3 += f"(?:{base_regex}){{s<=1}}"
                else:
                    pattern_3 += base_regex
            return f"(?e)(?:{pattern_5}){{s<={max_errors}}}{pattern_3}"
        else:
            # Brin anti-sens / Antisense strand (primer_seq is already primer_rc)
            primer_rc = primer_seq
            seq_5_rc = primer_rc[strict_3prime_len:]
            pattern_5_rc = seq_to_regex(seq_5_rc)
            pattern_3_rc = ""
            for k in range(1, strict_3prime_len + 1):
                base = primer_rc[k - 1]
                base_regex = seq_to_regex(base)
                if k in tolerate_positions:
                    pattern_3_rc += f"(?:{base_regex}){{s<=1}}"
                else:
                    pattern_3_rc += base_regex
            return f"(?e){pattern_3_rc}(?:{pattern_5_rc}){{s<={max_errors}}}"
    else:
        pattern = seq_to_regex(primer_seq)
        return f"(?e)({pattern}){{s<={max_errors}}}"


def primer_matches_sequence(
    target_seq: str,
    primer_seq: str,
    max_errors: int,
    strict_3prime_len: int = 3,
    strict_3prime_tolerate: int = 0
) -> Optional[Tuple[int, int, str]]:
    """Vérifie si l'amorce s'hybride à la séquence cible (brin sens ou anti-sens).
    Checks whether the primer matches the target sequence (sense or antisense strand).

    Retourne (start, end, strand) sur le brin sens (+), ou None si aucun match.
    Returns (start, end, strand) on the sense (+) strand, or None if no match.
    """
    # Brin sens / Sense strand
    regex_pattern = build_primer_regex(primer_seq, max_errors, strict_3prime_len, strict_3prime_tolerate, is_rc=False)
    match_sense = regex.search(regex_pattern, target_seq, regex.BESTMATCH)
    if match_sense:
        return match_sense.start(), match_sense.end(), '+'

    # Brin anti-sens / Antisense strand
    primer_rc = _revcomp(primer_seq)
    regex_rc_pattern = build_primer_regex(primer_rc, max_errors, strict_3prime_len, strict_3prime_tolerate, is_rc=True)
    match_antisense = regex.search(regex_rc_pattern, target_seq, regex.BESTMATCH)
    if match_antisense:
        return match_antisense.start(), match_antisense.end(), '-'

    return None


def count_iupac_mismatches(primer_seq: str, site_seq: str, strand: str) -> int:
    """Compte les mésappariements entre l'amorce et le site extrait par intersection IUPAC.
    Counts mismatches between primer and extracted site using IUPAC set intersection.
    """
    if len(primer_seq) != len(site_seq):
        raise ValueError(
            f"Longueur de séquence incohérente pour le comptage de mésappariements: "
            f"amorce ({len(primer_seq)} nt) vs site ({len(site_seq)} nt)."
        )

    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    if strand == '+':
        bottom_under_top = "".join(comp.get(c, 'N') for c in site_seq)
    else:
        bottom_under_top = site_seq[::-1]

    mismatches = 0
    for a, b in zip(primer_seq, bottom_under_top):
        a_set = IUPAC_MATCHABLE.get(a.upper())
        b_comp = comp.get(b.upper(), '')
        b_comp_set = IUPAC_MATCHABLE.get(b_comp)
        if a_set is None or b_comp_set is None or not a_set.intersection(b_comp_set):
            mismatches += 1
    return mismatches


def find_candidate_sites(
    strains: Dict[str, str],
    primers: List[PhysicalPrimer],
    max_errors: int = 2,
    strict_3prime_len: int = 3,
    strict_3prime_tolerate: int = 0,
    panel_name: str = "DefaultPanel"
) -> List[Dict[str, Any]]:
    """Crible permissivement les génomes de souches pour extraire tous les sites candidats.
    Permissively screens strain genomes to extract all candidate binding sites.

    Optimisation : criblage effectué UNE SEULE FOIS par rôle/amorce parente non développée.
    Optimization: screening performed ONCE per unexpanded parent role/primer.
    """
    # Rôles d'amorces participant à l'hybridation initiale
    # Primer roles involved in initial hybridization
    valid_roles = {
        PrimerRole.F3, PrimerRole.B3, PrimerRole.FIP, PrimerRole.BIP,
        PrimerRole.LF, PrimerRole.LB, PrimerRole.FWD, PrimerRole.REV
    }

    # Filtrer et pré-compiler les expressions régulières uniques
    # Filter and pre-compile unique regex patterns
    seen_seeds = set()
    compiled_primers = []
    for p in primers:
        if p.role not in valid_roles:
            continue
        p_seq = (p.parent_binding_domain or p.binding_domain).upper()
        p_key = (p.role, p.parent_name or p.name, p_seq)
        if p_key in seen_seeds:
            continue
        seen_seeds.add(p_key)

        pat_sense = build_primer_regex(p_seq, max_errors, strict_3prime_len, strict_3prime_tolerate, is_rc=False)
        pat_rc = build_primer_regex(_revcomp(p_seq), max_errors, strict_3prime_len, strict_3prime_tolerate, is_rc=True)
        compiled_sense = regex.compile(pat_sense, regex.BESTMATCH)
        compiled_rc = regex.compile(pat_rc, regex.BESTMATCH)
        p_display_name = p.parent_name or p.name
        compiled_primers.append((p.role, p_display_name, p_seq, compiled_sense, compiled_rc))

    candidate_records: List[Dict[str, Any]] = []

    for strain_id, genome_seq in strains.items():
        genome_upper = genome_seq.upper()
        for role, p_name, p_seq, reg_sense, reg_rc in compiled_primers:
            match = reg_sense.search(genome_upper)
            strand = '+'
            if match is None:
                match = reg_rc.search(genome_upper)
                strand = '-'

            if match is not None:
                start, end = match.start(), match.end()
                site_seq = genome_upper[start:end]
                n_mm = count_iupac_mismatches(p_seq, site_seq, strand)
                candidate_records.append({
                    "strain_id": strain_id,
                    "primer_role": role.value,
                    "primer_name": p_name,
                    "position": start,
                    "strand": strand,
                    "site_seq": site_seq,
                    "n_mismatches": n_mm,
                    "panel": panel_name
                })

    return candidate_records
