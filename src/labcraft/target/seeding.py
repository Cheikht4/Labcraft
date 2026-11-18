"""Module de criblage permissif de sites candidats pour le mode couverture.
Permissive candidate site seeding module for coverage analysis.

Origine / Origin:
Adapté de 'lamp_coverage.py' (dépôt primer-analysis-suite, même auteur).
Rôle : criblage permissif produisant des sites CANDIDATS, re-notés ensuite
par la couche thermodynamique de LabCraft.
Role: permissive screening generating CANDIDATE sites, which are then re-scored
by the thermodynamic layer of LabCraft.

Modèle de criblage / Screening model:
Le criblage est PERMISSIF par construction : il ne doit filtrer que ce que la
couche thermodynamique ne saurait pas noter (les indels, via substitutions
seules {s<=n}), jamais ce qu'elle sait noter (les mésappariements, y compris
en 3'). Par conséquent, strict_3prime_tolerate=2 par défaut afin de transmettre
les sites mutés en 3' à three_prime_extensible qui portera le diagnostic veto_3p.
The screening is PERMISSIVE by design: it filters only what thermodynamics
cannot score (indels, via substitutions only {s<=n}), never what it can score
(mismatches, including 3' terminal ones). Therefore, strict_3prime_tolerate=2
by default so that 3' mutated sites are passed to three_prime_extensible for
proper veto_3p diagnosis.
"""

from __future__ import annotations
from typing import Dict, List, Optional, Tuple, Any
import regex

from labcraft.lamp.domains import PhysicalPrimer, PrimerRole, _revcomp, IUPAC_MATCHABLE

# Dictionnaire IUPAC vers classes d'expression régulière (inclut N sur la cible pour ne pas rejeter les sites contenant des N)
# IUPAC dictionary to regular expression classes (includes N on target so sites with N are not rejected as substitutions)
IUPAC_DICT = {
    'A': '[AN]', 'C': '[CN]', 'G': '[GN]', 'T': '[TN]', 'U': '[TN]',
    'R': '[AGN]', 'Y': '[CTN]', 'S': '[GCN]', 'W': '[ATN]',
    'K': '[GTN]', 'M': '[ACN]', 'B': '[CGTN]', 'D': '[AGTN]',
    'H': '[ACTN]', 'V': '[ACGN]', 'N': '[ACGTN]'
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
    strict_3prime_tolerate: int = 2,
    is_rc: bool = False
) -> str:
    """Construit une expression régulière tolérante en 5' et permissive en 3' (substitutions seules).
    Builds a regex pattern tolerant in 5' and permissive in 3' (substitutions only).
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


def _search_best_non_n(
    compiled_reg: Any,
    genome_upper: str,
) -> Optional[Tuple[int, int]]:
    """Recherche le meilleur site non dominé par des bases 'N' (franchit les paddings de N en O(1)).
    Searches for the best site not dominated by 'N' bases (skips 'N' padding runs in O(1)).
    """
    pos = 0
    g_len = len(genome_upper)
    while pos < g_len:
        m = compiled_reg.search(genome_upper, pos)
        if m is None:
            return None
        start, end = m.start(), m.end()
        win_len = end - start
        if win_len <= 0:
            pos = start + 1
            continue
        site_seq = genome_upper[start:end]
        max_n = win_len // 2
        if site_seq.count('N') <= max_n:
            return start, end
        # Sauter directement au-delà d'un bloc contigu de 'N' en tête sans rater de fenêtre valide
        # Jump directly past a leading contiguous 'N' run without skipping any valid window
        n_end = start
        while n_end < g_len and genome_upper[n_end] == 'N':
            n_end += 1
        pos = max(start + 1, n_end - max_n)
    return None


def _has_3prime_mismatch(primer_seq: str, site_seq: str, strand: str, window_3p: int = 2) -> bool:
    """Indique si l'une des `window_3p` dernières bases 3' présente un mésappariement IUPAC.
    Returns True if any of the last `window_3p` 3' bases has an IUPAC mismatch.
    """
    if window_3p <= 0 or not primer_seq or len(primer_seq) != len(site_seq):
        return False
    w = min(window_3p, len(primer_seq))
    p_tail = primer_seq[-w:]
    s_tail = site_seq[-w:] if strand == '+' else site_seq[:w]
    return count_iupac_mismatches(p_tail, s_tail, strand) > 0


def primer_matches_sequence(
    target_seq: str,
    primer_seq: str,
    max_errors: int,
    strict_3prime_len: int = 3,
    strict_3prime_tolerate: int = 2
) -> Optional[Tuple[int, int, str]]:
    """Vérifie si l'amorce s'hybride à la séquence cible en évaluant les DEUX brins (+ et -).
    Checks whether the primer matches the target sequence by evaluating BOTH strands (+ and -).

    Retourne (start, end, strand) sur le brin sens (+), ou None si aucun match.
    Returns (start, end, strand) on the sense (+) strand, or None if no match.
    """
    target_upper = target_seq.upper()
    p_upper = primer_seq.upper()

    # Brin sens / Sense strand (+)
    reg_sense = regex.compile(
        build_primer_regex(p_upper, max_errors, strict_3prime_len, strict_3prime_tolerate, is_rc=False),
        regex.BESTMATCH,
    )
    hit_sense = _search_best_non_n(reg_sense, target_upper)

    # Brin anti-sens / Antisense strand (-)
    primer_rc = _revcomp(p_upper)
    reg_rc = regex.compile(
        build_primer_regex(primer_rc, max_errors, strict_3prime_len, strict_3prime_tolerate, is_rc=True),
        regex.BESTMATCH,
    )
    hit_rc = _search_best_non_n(reg_rc, target_upper)

    if hit_sense is not None and hit_rc is None:
        return hit_sense[0], hit_sense[1], '+'
    if hit_rc is not None and hit_sense is None:
        return hit_rc[0], hit_rc[1], '-'
    if hit_sense is None and hit_rc is None:
        return None

    # Les deux brins matchent : retenir le site physiquement supérieur (3' intact d'abord, puis moins de mismatches)
    # Both strands match: retain the physically superior locus (intact 3' end first, then fewer mismatches)
    assert hit_sense is not None and hit_rc is not None
    seq_s = target_upper[hit_sense[0]:hit_sense[1]]
    seq_r = target_upper[hit_rc[0]:hit_rc[1]]
    mm_s = count_iupac_mismatches(p_upper, seq_s, '+')
    mm_r = count_iupac_mismatches(p_upper, seq_r, '-')
    bad3_s = _has_3prime_mismatch(p_upper, seq_s, '+')
    bad3_r = _has_3prime_mismatch(p_upper, seq_r, '-')

    if (bad3_r, mm_r) < (bad3_s, mm_s):
        return hit_rc[0], hit_rc[1], '-'
    return hit_sense[0], hit_sense[1], '+'


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
    for a, b in zip(primer_seq, bottom_under_top, strict=True):
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
    strict_3prime_tolerate: int = 2,
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
            # Évaluer systématiquement les DEUX brins (+ et -) sans court-circuit par un match parasite sur '+'
            # Systematically evaluate BOTH strands (+ and -) without short-circuiting on a noisy '+' hit
            hit_sense = _search_best_non_n(reg_sense, genome_upper)
            hit_rc = _search_best_non_n(reg_rc, genome_upper)

            strand_hits: List[Tuple[int, int, str, str, int, bool]] = []
            if hit_sense is not None:
                s_seq = genome_upper[hit_sense[0]:hit_sense[1]]
                mm_s = count_iupac_mismatches(p_seq, s_seq, '+')
                b3_s = _has_3prime_mismatch(p_seq, s_seq, '+')
                strand_hits.append((hit_sense[0], hit_sense[1], '+', s_seq, mm_s, b3_s))
            if hit_rc is not None:
                r_seq = genome_upper[hit_rc[0]:hit_rc[1]]
                mm_r = count_iupac_mismatches(p_seq, r_seq, '-')
                b3_r = _has_3prime_mismatch(p_seq, r_seq, '-')
                strand_hits.append((hit_rc[0], hit_rc[1], '-', r_seq, mm_r, b3_r))

            if len(strand_hits) == 2:
                h0, h1 = strand_hits[0], strand_hits[1]
                # Si les deux brins ont <= max_errors (vraie répétition inversée), garder les deux ;
                # sinon éliminer le match parasite issu du budget supplémentaire toléré en 3'.
                # If both strands have <= max_errors (true inverted repeat), keep both;
                # otherwise drop the noisy hit admitted only via the extra 3' tolerance budget.
                if not (h0[4] <= max_errors and h1[4] <= max_errors):
                    k0 = (1 if h0[5] else 0, h0[4])
                    k1 = (1 if h1[5] else 0, h1[4])
                    if k0 < k1:
                        strand_hits = [h0]
                    elif k1 < k0:
                        strand_hits = [h1]

            for start, _end, strand, site_seq, n_mm, _b3 in strand_hits:
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
