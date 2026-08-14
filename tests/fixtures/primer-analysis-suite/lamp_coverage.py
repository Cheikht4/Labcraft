#!/usr/bin/env python3
import argparse
import sys
import os
import re
import itertools
import regex
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.Seq import Seq
import concurrent.futures
import math

# Import optionnel de tqdm pour les barres de progression / Optional tqdm import for progress bars
try:
    from tqdm import tqdm
except ImportError:
    # Fallback si tqdm n'est pas installé / Fallback if tqdm is not installed
    def tqdm(iterable=None, **kwargs):
        desc = kwargs.get('desc', '')
        if desc and iterable is None:
            class _DummyBar:
                def update(self, n=1): pass
                def set_postfix_str(self, s): pass
                def close(self): pass
            return _DummyBar()
        if desc:
            print(f"{desc}...")
        return iterable if iterable is not None else []

# Dictionnaire IUPAC vers expression régulière
IUPAC_DICT = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
}

# =============================================================================
# Types canoniques d'amorces (utilisés pour la détection des multi-versions)
# Canonical primer types (used for multi-version detection)
# =============================================================================

# LAMP : noms de positions canoniques / LAMP: canonical position names
LAMP_CANONICAL_TYPES = {
    'F3', 'B3', 'F2', 'F1', 'B1', 'B2',
    'FLOOP', 'BLOOP', 'STEMF', 'STEMB', 'FIP', 'BIP'
}
# PCR : types canoniques (après résolution des alias) / PCR: canonical types (after alias resolution)
PCR_CANONICAL_TYPES = {'F', 'R', 'P'}


def get_base_type(primer_id, is_pcr=False):
    """
    Retourne le type de base d'une amorce en supprimant le suffixe numérique de version.
    Ex : F31 → F3, P2 → P, FLOOP3 → FLOOP.
    Si le suffixe ne correspond pas à un type canonique, retourne le nom original.

    Returns the base type of a primer by stripping the numeric version suffix.
    E.g.: F31 → F3, P2 → P, FLOOP3 → FLOOP.
    If the stripped name is not a canonical type, returns the original name.
    """
    canonical = PCR_CANONICAL_TYPES if is_pcr else LAMP_CANONICAL_TYPES
    stripped = primer_id.rstrip('0123456789')
    if stripped != primer_id and stripped in canonical:
        return stripped
    return primer_id


def group_primer_versions(primers_dict, is_pcr=False):
    """
    Groupe les amorces d'un set par leur type de base.
    Ex : {'F31': seq, 'F32': seq, 'B3': seq} → {'F3': ['F31','F32'], 'B3': ['B3']}

    Groups primers of a set by their base type.
    E.g.: {'F31': seq, 'F32': seq, 'B3': seq} → {'F3': ['F31','F32'], 'B3': ['B3']}
    """
    groups = defaultdict(list)
    for primer_id in primers_dict.keys():
        base = get_base_type(primer_id, is_pcr)
        groups[base].append(primer_id)
    return dict(groups)


def get_all_nonempty_subsets(lst):
    """Retourne tous les sous-ensembles non-vides d'une liste sous forme de tuples triés."""
    subsets = []
    for r in range(1, len(lst) + 1):
        for combo in itertools.combinations(lst, r):
            subsets.append(tuple(sorted(combo)))
    return subsets


def find_best_subsets_by_pool(type_groups, primer_matches_set, max_exhaustive=50_000):
    """
    Calcule la meilleure combinaison d'amorces pour chaque taille de pool M,
    de M_min (1 version par type) à M_max (toutes les versions de tous les types).
    Si le nombre de combinaisons de sous-ensembles <= max_exhaustive: recherche exhaustive.
    Sinon: algorithme glouton progressif.

    Calculates the best primer combination for each pool size M,
    from M_min (1 version per type) to M_max (all versions of all types).
    If subset combinations count <= max_exhaustive: exhaustive search.
    Otherwise: progressive greedy algorithm.
    """
    types = list(type_groups.keys())
    
    # Générer tous les sous-ensembles non-vides pour chaque type
    subsets_by_type = {}
    for t in types:
        subsets_by_type[t] = get_all_nonempty_subsets(type_groups[t])
        
    # Calculer le nombre total de combinaisons de sous-ensembles
    n_combos = 1
    for t in types:
        n_combos *= len(subsets_by_type[t])
        
    # Pré-calculer les unions de matchs pour chaque sous-ensemble pour accélérer
    union_match = {}
    for t in types:
        union_match[t] = {}
        for s in subsets_by_type[t]:
            union_match[t][s] = set.union(*[primer_matches_set.get(v, set()) for v in s]) if s else set()

    best_by_pool = {}  # M -> (best_combo: {type: tuple_of_versions}, best_coverage: set, algo: str)

    if n_combos <= max_exhaustive:
        # Recherche exhaustive parmi toutes les combinaisons de sous-ensembles
        for combo_subsets in itertools.product(*[subsets_by_type[t] for t in types]):
            m_size = sum(len(s) for s in combo_subsets)
            match_sets = [union_match[types[i]][combo_subsets[i]] for i in range(len(types))]
            coverage = set.intersection(*match_sets) if match_sets else set()
            
            if m_size not in best_by_pool or len(coverage) > len(best_by_pool[m_size][1]):
                best_combo_dict = {types[i]: combo_subsets[i] for i in range(len(types))}
                best_by_pool[m_size] = (best_combo_dict, coverage, 'exhaustive')
    else:
        # Algorithme glouton progressif (greedy forward selection)
        # Étape 1 : M = n (1 version par type). Recherche exhaustive sur les tailles = 1 par type.
        n_combos_1 = 1
        for t in types:
            n_combos_1 *= len(type_groups[t])
            
        best_combo_1 = {}
        best_cov_1 = set()
        
        if n_combos_1 <= max_exhaustive:
            v_lists = [type_groups[t] for t in types]
            for combo_versions in itertools.product(*v_lists):
                match_sets = [primer_matches_set.get(v, set()) for v in combo_versions]
                coverage = set.intersection(*match_sets) if match_sets else set()
                if len(coverage) > len(best_cov_1) or not best_combo_1:
                    best_cov_1 = coverage
                    best_combo_1 = {types[i]: (combo_versions[i],) for i in range(len(types))}
        else:
            chosen = []
            for t in types:
                best_v = max(type_groups[t], key=lambda v: len(primer_matches_set.get(v, set())))
                chosen.append(best_v)
            match_sets = [primer_matches_set.get(v, set()) for v in chosen]
            best_cov_1 = set.intersection(*match_sets) if match_sets else set()
            best_combo_1 = {types[i]: (chosen[i],) for i in range(len(types))}
            
        best_by_pool[len(types)] = (best_combo_1, best_cov_1, 'greedy')
        
        # Étapes suivantes : M de n+1 à sum(|V_t|)
        curr_combo = {t: list(best_combo_1[t]) for t in types}
        m_start = len(types)
        m_max = sum(len(type_groups[t]) for t in types)
        
        for m in range(m_start + 1, m_max + 1):
            best_candidate_combo = None
            best_candidate_cov = None
            best_candidate_len = -1
            
            for t in types:
                if len(curr_combo[t]) < len(type_groups[t]):
                    for v in type_groups[t]:
                        if v not in curr_combo[t]:
                            cand_combo = {type_name: list(versions) for type_name, versions in curr_combo.items()}
                            cand_combo[t].append(v)
                            
                            match_sets = []
                            for tn, vers in cand_combo.items():
                                union_set = set.union(*[primer_matches_set.get(vi, set()) for vi in vers])
                                match_sets.append(union_set)
                            cov = set.intersection(*match_sets) if match_sets else set()
                            
                            if len(cov) > best_candidate_len:
                                best_candidate_len = len(cov)
                                best_candidate_cov = cov
                                best_candidate_combo = {tn: tuple(sorted(vers)) for tn, vers in cand_combo.items()}
                                
            if best_candidate_combo:
                best_by_pool[m] = (best_candidate_combo, best_candidate_cov, 'greedy')
                curr_combo = {t: list(best_candidate_combo[t]) for t in types}
            else:
                break
                
    return best_by_pool, n_combos


def analyze_marginal_value(type_groups, primer_matches_set):
    """
    Calcule pour chaque type d'amorce multi-versions la couverture optimale pour k versions,
    lorsque tous les autres types d'amorces utilisent toutes leurs versions.

    Calculates for each multi-version primer type the optimal coverage for k versions,
    when all other primer types use all of their versions.
    """
    types = list(type_groups.keys())
    
    # U_{-t} = intersection_{t' != t} (union_{v' in V_{t'}} matches(v'))
    unions_other = {}
    for t in types:
        other_unions = []
        for t_other in types:
            if t_other != t:
                u = set.union(*[primer_matches_set.get(v, set()) for v in type_groups[t_other] if primer_matches_set.get(v) is not None])
                other_unions.append(u)
        unions_other[t] = set.intersection(*other_unions) if other_unions else None

    marginal_by_type = {}
    
    for t in types:
        versions = type_groups[t]
        if len(versions) <= 1:
            continue
            
        marginal_by_type[t] = []
        u_other = unions_other[t]
        
        for k in range(1, len(versions) + 1):
            best_subset = None
            best_cov = set()
            best_len = -1
            
            for subset in itertools.combinations(versions, k):
                u_subset = set.union(*[primer_matches_set.get(v, set()) for v in subset if primer_matches_set.get(v) is not None])
                if u_other is not None:
                    cov = u_subset.intersection(u_other)
                else:
                    cov = u_subset
                    
                if len(cov) > best_len or best_subset is None:
                    best_len = len(cov)
                    best_cov = cov
                    best_subset = tuple(sorted(subset))
                    
            marginal_by_type[t].append((k, best_subset, best_cov))
            
    return marginal_by_type


def pigeonhole_search_simple(seq_upper, kmer, kmer_offset):
    kmer_len = 18
    max_kmer_err = 2
    n_segments = 3
    seg_len = 6
    max_cands = 15
    n = len(seq_upper)
    for seg_i in range(n_segments):
        seg_start_in_kmer = seg_i * seg_len
        segment = kmer[seg_start_in_kmer : seg_start_in_kmer + seg_len]
        cands = 0
        hit = seq_upper.find(segment)
        while hit >= 0 and cands < max_cands:
            est_kmer_start = hit - seg_start_in_kmer
            if 0 <= est_kmer_start <= n - kmer_len:
                window = seq_upper[est_kmer_start : est_kmer_start + kmer_len]
                mismatches = sum(a != b for a, b in zip(kmer, window))
                if mismatches <= max_kmer_err:
                    return max(0, est_kmer_start - kmer_offset)
            cands += 1
            hit = seq_upper.find(segment, hit + 1)
    return None


def extract_amplicon_region(seq_id, targets, gapped_targets, ref_seq_id, ref_seq, ref_gapped, amp_start, amp_end, amp_len, is_msa, mapping, blast_results, kmers):
    # 1. Chemin MSA
    if is_msa and mapping:
        gapped_start = mapping[amp_start]
        gapped_end   = mapping[amp_end - 1] + 1
        seq_gapped = gapped_targets.get(seq_id, '')
        if seq_gapped:
            region_gapped = seq_gapped[gapped_start:gapped_end]
            region = region_gapped.replace('-', '')
            if region:
                return region
                
    # 2. Chemin BLAST
    if blast_results and seq_id in blast_results:
        hit = blast_results[seq_id]
        seq = targets.get(seq_id, '')
        if seq:
            if hit['sstart'] < hit['send']:
                est_sstart = hit['sstart'] - (hit['qstart'] - 1)
                est_send = hit['send'] + (amp_len - hit['qend'])
            else:
                est_sstart = hit['send'] - (hit['qstart'] - 1)
                est_send = hit['sstart'] + (amp_len - hit['qend'])
            s_start = max(0, est_sstart - 1)
            s_end = min(len(seq), est_send)
            region = seq[s_start:s_end]
            if region:
                return region

    # 3. Chemin Pigeonhole
    seq = targets.get(seq_id, '')
    if seq and kmers:
        seq_upper = seq.upper()
        # Essai brin sens
        for kmer_offset, kmer, _ in kmers:
            amp_est = pigeonhole_search_simple(seq_upper, kmer, kmer_offset)
            if amp_est is not None:
                r_end = min(len(seq), amp_est + amp_len)
                return seq[amp_est:r_end]
        # Essai brin anti-sens
        for kmer_offset, _, kmer_rc in kmers:
            amp_est = pigeonhole_search_simple(seq_upper, kmer_rc, kmer_offset)
            if amp_est is not None:
                r_end = min(len(seq), amp_est + amp_len)
                return seq[amp_est:r_end]
                
    return None


def check_relaxed_primer_match(seq_amp, primer_seq):
    # Convertit l'amorce en regex
    pattern = seq_to_regex(primer_seq)
    # Recherche avec max 6 erreurs
    regex_pattern = f"(?e)({pattern}){{e<=6}}"
    match = regex.search(regex_pattern, seq_amp, regex.BESTMATCH)
    if match:
        return True
    
    # Brin anti-sens
    primer_rc = str(Seq(primer_seq).reverse_complement())
    pattern_rc = seq_to_regex(primer_rc)
    regex_rc_pattern = f"(?e)({pattern_rc}){{e<=6}}"
    match_rc = regex.search(regex_rc_pattern, seq_amp, regex.BESTMATCH)
    if match_rc:
        return True
        
    return False


def diagnose_nonmatching_sequences(
    non_matched_ids, targets, gapped_targets, primers_dict, primer_positions, set_id,
    valid_order_matches, is_pcr, max_n_pct_diag, primer_matches_set
):
    """
    Diagnostique les séquences non-matchées à partir de la couverture des amorces et de l'alignement/BLAST.
    
    Classifie en:
    - too_variable             : au moins 1 amorce se fixe (divergence réelle)
    - truncated_poor_divergent : 0 amorce se fixe, région tronquée ou mauvaise qualité (>max_n_pct_diag % de N dans l'amplicon)
    """
    # ── Séquence de référence (premier match valide) pour l'amplicon ──
    ref_seq_id = next(iter(valid_order_matches), None)
    if ref_seq_id is None:
        return {
            'too_variable': [],
            'truncated_poor_divergent': non_matched_ids,
            'details': {sid: {'status': 'truncated_poor_divergent', 'matched': [], 'unmatched': {p: 'not_found' for p in primers_dict}} for sid in non_matched_ids},
            'ref_amp_len': 0,
            'method': 'Primer Binding Analysis'
        }

    ref_seq = targets.get(ref_seq_id, '')
    ref_pos = primer_positions.get(set_id, {}).get(ref_seq_id, {})
    if not ref_seq or not ref_pos:
        return {
            'too_variable': [],
            'truncated_poor_divergent': non_matched_ids,
            'details': {sid: {'status': 'truncated_poor_divergent', 'matched': [], 'unmatched': {p: 'not_found' for p in primers_dict}} for sid in non_matched_ids},
            'ref_amp_len': 0,
            'method': 'Primer Binding Analysis'
        }

    start_type = 'F' if is_pcr else 'F3'
    end_type   = 'R' if is_pcr else 'B3'
    amp_start, amp_end = None, None
    for pid, pos in ref_pos.items():
        bt = get_base_type(pid, is_pcr)
        if bt == start_type and (amp_start is None or pos[0] < amp_start):
            amp_start = pos[0]
        if bt == end_type and (amp_end is None or pos[1] > amp_end):
            amp_end = pos[1]

    if amp_start is None:
        starts = [p[0] for p in ref_pos.values()]
        amp_start = min(starts) if starts else None
    if amp_end is None:
        ends = [p[1] for p in ref_pos.values()]
        amp_end = max(ends) if ends else None

    if amp_start is None or amp_end is None or amp_end <= amp_start:
        return {
            'too_variable': [],
            'truncated_poor_divergent': non_matched_ids,
            'details': {sid: {'status': 'truncated_poor_divergent', 'matched': [], 'unmatched': {p: 'not_found' for p in primers_dict}} for sid in non_matched_ids},
            'ref_amp_len': 0,
            'method': 'Primer Binding Analysis'
        }

    amp_len = amp_end - amp_start
    ref_amp = ref_seq[amp_start:amp_end]

    # Détection MSA
    ref_gapped = gapped_targets.get(ref_seq_id, '')
    ref_gapped_len = len(ref_gapped)
    sample = list(non_matched_ids)[:30]
    is_msa = ref_gapped_len > 0 and all(len(gapped_targets.get(sid, '')) == ref_gapped_len for sid in sample)

    mapping = []
    if is_msa:
        for i, char in enumerate(ref_gapped):
            if char != '-':
                mapping.append(i)

    # BLAST (si pas MSA)
    blast_results = {}
    if not is_msa:
        import shutil
        import tempfile
        import subprocess
        blastn_path = shutil.which("blastn") or "/opt/homebrew/bin/blastn"
        if not os.path.exists(blastn_path):
            blastn_path = "/usr/local/bin/blastn" if os.path.exists("/usr/local/bin/blastn") else None
            
        if blastn_path:
            try:
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as q_file:
                    q_file.write(f">ref\n{ref_amp}\n")
                    q_path = q_file.name
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as s_file:
                    for seq_id in non_matched_ids:
                        seq = targets.get(seq_id, '')
                        if seq:
                            s_file.write(f">{seq_id}\n{seq}\n")
                    s_path = s_file.name
                cmd = [
                    blastn_path, "-query", q_path, "-subject", s_path,
                    "-outfmt", "6 sseqid pident length qstart qend sstart send gaps",
                    "-perc_identity", "50"
                ]
                proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                if proc.returncode == 0:
                    for line in proc.stdout.strip().split("\n"):
                        if not line: continue
                        parts = line.split("\t")
                        if len(parts) >= 8:
                            sseqid, pident, aln_len, qstart, qend, sstart, send = parts[0], float(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]), int(parts[5]), int(parts[6])
                            qcov = (qend - qstart + 1) / amp_len * 100
                            if qcov >= 50.0:
                                if sseqid not in blast_results or qcov > blast_results[sseqid]['qcov']:
                                    blast_results[sseqid] = {
                                        'sstart': sstart, 'send': send, 'qstart': qstart, 'qend': qend, 'qcov': qcov
                                    }
                os.remove(q_path)
                os.remove(s_path)
            except:
                pass

    # Kmers Pigeonhole
    kmers = []
    for i in range(8):
        offset = (amp_len * i) // 9
        pos = amp_start + offset
        if pos + 18 <= len(ref_seq):
            kmer = ref_seq[pos:pos+18].upper()
            if 'N' not in kmer and len(kmer) == 18:
                kmers.append((offset, kmer, str(Seq(kmer).reverse_complement())))

    too_variable = []
    truncated_poor_divergent = []
    details = {}

    for seq_id in non_matched_ids:
        matched_primers = []
        unmatched_primers = []
        for p_id in primers_dict:
            if seq_id in primer_matches_set.get(p_id, set()):
                matched_primers.append(p_id)
            else:
                unmatched_primers.append(p_id)
                
        # 1. Règle absolue : 0 amorce fixée -> tronqué/pauvre/très divergent
        if len(matched_primers) == 0:
            truncated_poor_divergent.append(seq_id)
            details[seq_id] = {
                'status': 'truncated_poor_divergent',
                'matched': [],
                'unmatched': {p: 'not_found' for p in unmatched_primers}
            }
            continue
            
        # 2. Si au moins 1 amorce fixée, extraire la région de l'amplicon
        seq_amp = extract_amplicon_region(
            seq_id, targets, gapped_targets, ref_seq_id, ref_seq, ref_gapped,
            amp_start, amp_end, amp_len, is_msa, mapping, blast_results, kmers
        )
        
        # 3. Vérifier la qualité de l'amplicon extrait
        if seq_amp:
            n_pct = (seq_amp.count('N') / len(seq_amp)) * 100
            if n_pct > max_n_pct_diag:
                # Mauvaise qualité -> classé en truncated/poor/divergent
                truncated_poor_divergent.append(seq_id)
                details[seq_id] = {
                    'status': 'truncated_poor_divergent',
                    'matched': matched_primers,
                    'unmatched': {p: 'poor_quality' for p in unmatched_primers}
                }
                continue
                
        # 4. Si qualité OK (ou région non extraite mais au moins 1 amorce fixée) -> trop variable
        too_variable.append(seq_id)
        unmatched_status = {}
        for p_id in unmatched_primers:
            primer_seq = primers_dict[p_id]
            if seq_amp and check_relaxed_primer_match(seq_amp, primer_seq):
                unmatched_status[p_id] = 'too_many_errors'
            else:
                unmatched_status[p_id] = 'not_found'
                
        details[seq_id] = {
            'status': 'too_variable',
            'matched': matched_primers,
            'unmatched': unmatched_status
        }

    return {
        'too_variable': too_variable,
        'truncated_poor_divergent': truncated_poor_divergent,
        'details': details,
        'ref_amp_len': amp_len,
        'method': 'Primer Binding + Align Analysis'
    }



def seq_to_regex(seq):
    """Convertit une séquence avec codes IUPAC en expression régulière."""

    pattern = ""
    for char in seq.upper():
        if char.isalpha():
            pattern += IUPAC_DICT.get(char, char)
    return pattern

def primer_matches_sequence(target_seq, primer_seq, max_errors, strict_3prime_len=3, strict_3prime_tolerate=0):
    """
    Vérifie si l'amorce match la séquence cible (brin sens ou anti-sens).
    Retourne un tuple (start, end) de la position du match sur le brin sens,
    ou None s'il n'y a pas de match.
    """
    # Définition des positions tolérées / Definition of tolerated positions
    # Position 1 = base -1, Position 2 = base -2, Position 3 = base -3 (compté depuis l'extrémité 3' / counted from 3' end)
    tolerate_positions = set()
    if strict_3prime_tolerate == 1:
        tolerate_positions = {2}
    elif strict_3prime_tolerate == 2:
        tolerate_positions = {1, 2}

    def build_regex(seq, max_e, strict_len, is_rc=False):
        if strict_len > 0 and len(seq) > strict_len:
            # Choix du type d'erreur pour la partie 5' : substitutions uniquement (s) en mode ARMS pour éviter les décalages d'indels
            # Choice of error type for the 5' part: substitutions only (s) in ARMS mode to avoid indel shifts
            type_e = 's' if strict_3prime_tolerate > 0 else 'e'
            
            if not is_rc:
                # Brin sens / Sense strand
                seq_5 = seq[:-strict_len]
                pattern_5 = seq_to_regex(seq_5)
                
                # Construction de la zone 3' base par base / Building the 3' region base by base
                pattern_3 = ""
                # De la base la plus éloignée de 3' vers la base terminale (5' vers 3') / From the base furthest from 3' to the terminal base (5' to 3')
                for k in range(strict_len, 0, -1):
                    base = seq[-k]
                    base_regex = seq_to_regex(base)
                    if k in tolerate_positions:
                        pattern_3 += f"(?:{base_regex}){{s<=1}}"
                    else:
                        pattern_3 += base_regex
                
                return f"(?e)(?:{pattern_5}){{{type_e}<={max_e}}}{pattern_3}"
            else:
                # Brin anti-sens / Antisense strand
                primer_rc = seq # Ici seq est déjà primer_rc / Here seq is already primer_rc
                seq_5_rc = primer_rc[strict_len:]
                pattern_5_rc = seq_to_regex(seq_5_rc)
                
                # Le 3' de l'amorce originale correspond au début (5') de primer_rc / The 3' of the original primer corresponds to the beginning (5') of primer_rc
                pattern_3_rc = ""
                for k in range(1, strict_len + 1):
                    base = primer_rc[k-1]
                    base_regex = seq_to_regex(base)
                    if k in tolerate_positions:
                        pattern_3_rc += f"(?:{base_regex}){{s<=1}}"
                    else:
                        pattern_3_rc += base_regex
                
                return f"(?e){pattern_3_rc}(?:{pattern_5_rc}){{{type_e}<={max_e}}}"
        else:
            pattern = seq_to_regex(seq)
            return f"(?e)({pattern}){{e<={max_e}}}"

    # Brin sens
    regex_pattern = build_regex(primer_seq, max_errors, strict_3prime_len, is_rc=False)
    match_sense = regex.search(regex_pattern, target_seq, regex.BESTMATCH)
    if match_sense:
        return match_sense.span()
        
    # Brin anti-sens
    primer_rc = str(Seq(primer_seq).reverse_complement())
    regex_rc_pattern = build_regex(primer_rc, max_errors, strict_3prime_len, is_rc=True)
    match_antisense = regex.search(regex_rc_pattern, target_seq, regex.BESTMATCH)
    
    if match_antisense:
        return match_antisense.span()
        
    return None

def load_primers(filepath, is_pcr=False):
    """
    Charge les amorces et les regroupe par set.
    Attend un format de nom du type >SetID_PrimerName (ex: >2_F3).
    """
    primer_sets = defaultdict(dict)
    
    try:
        with open(filepath, 'r') as f:
            first_line = f.readline()
            if not first_line: return {}
            
            f.seek(0)
            first_char = f.read(1)
            f.seek(0)
            
            records = []
            if first_char == '>':
                records = list(SeqIO.parse(f, "fasta"))
            else:
                # Fichier texte basique
                primer_count = 1
                for line in f:
                    line = line.strip()
                    if not line: continue
                    parts = line.split()
                    if len(parts) >= 2:
                        name = parts[0]
                        seq = "".join(parts[1:])
                    else:
                        name = f"DefaultSet_Primer{primer_count}"
                        seq = parts[0]
                        primer_count += 1
                    
                    seq = "".join(c for c in seq if c.isalpha())
                    records.append({'id': name, 'seq': seq})
                    
            # Initialisation de la liste de collecte pour gérer les doublons d'ID d'amorce
            # Initialize collection list to handle duplicate primer IDs within the same set
            collected = []

            for record in records:
                if hasattr(record, 'id'):
                    name = record.id
                    seq = str(record.seq).upper()
                else:
                    name = record['id']
                    seq = record['seq'].upper()
                    
                if '_' in name:
                    parts = name.rsplit('_', 1)
                    set_id   = parts[0]
                    primer_id = parts[1]

                    # NOUVEAU : si le dernier segment est un entier pur → numéro de version
                    # NEW: if the last segment is a pure integer → version number
                    # Ex: "SetA_F3_2" → set_id="SetA", primer_id="F3", version=2
                    # La résolution des doublons (post-traitement ci-dessous) gérera le renommage
                    # Duplicate resolution (post-processing below) will handle renaming
                    if primer_id.isdigit():
                        inner = set_id.rsplit('_', 1)
                        if len(inner) == 2:
                            set_id    = inner[0]
                            primer_id = inner[1]  # vrai type d'amorce / actual primer type
                        # Sinon : set_id court, garder les parties telles quelles
                        # Otherwise: short set_id, keep parts as-is
                else:
                    set_id = "Default"
                    primer_id = name
                    
                # Gestion des alias / Alias resolution
                if is_pcr:
                    alias_map = {
                        'FWD': 'F', 'FORWARD': 'F', 'F3': 'F', 'FP': 'F',
                        'REV': 'R', 'REVERSE': 'R', 'B3': 'R', 'RP': 'R',
                        'PROBE': 'P', 'SONDE': 'P'
                    }
                else:
                    alias_map = {
                        'BLP': 'BLOOP',
                        'FLP': 'FLOOP',
                        'LOOPF': 'FLOOP',
                        'LOOPB': 'BLOOP',
                        'LF': 'FLOOP',
                        'LB': 'BLOOP',
                        'F1C': 'F1',
                        'B1C': 'B1'
                    }
                # Résolution flexible des alias (tolère suffixes/préfixes numériques ou alphabétiques)
                # Flexible alias resolution (tolerates numeric/alpha suffixes and prefixes)
                primer_id_upper = primer_id.upper()
                resolved = alias_map.get(primer_id_upper)  # 1. Correspondance exacte / Exact match
                if resolved is None:
                    # 2. Recherche par préfixe (ex: FP1 → FP → F, RP2 → RP → R)
                    # Prefix search (e.g. FP1 → FP → F, RP2 → RP → R)
                    for key in sorted(alias_map.keys(), key=len, reverse=True):
                        if primer_id_upper.startswith(key):
                            resolved = alias_map[key]
                            break
                if resolved is None:
                    # 3. Suppression du suffixe non-alpha final (ex: Probe1 → Probe, RP2 → RP)
                    # Strip trailing non-alpha suffix (e.g. Probe1 → Probe, RP2 → RP)
                    stripped = re.sub(r'[^A-Z]+$', '', primer_id_upper)
                    resolved = alias_map.get(stripped)
                if resolved is None:
                    # 4. Recherche par suffixe (ex: SFP → FP → F, SRP → RP → R)
                    # Suffix search (e.g. SFP → FP → F, SRP → RP → R)
                    for key in sorted(alias_map.keys(), key=len, reverse=True):
                        if primer_id_upper.endswith(key):
                            resolved = alias_map[key]
                            break
                if resolved is not None:
                    primer_id = resolved
                # Sinon, on conserve le nom original / Otherwise keep original name

                # Collecter pour post-traitement des doublons / Collect for duplicate post-processing
                collected.append((set_id, primer_id, seq))

            # Post-traitement : renommer les doublons d'ID au sein du même set
            # Ex : P, P → P1, P2 ; R, R → R1, R2 (sondes/amorces multiples)
            # Post-processing: rename duplicate IDs within the same set
            # E.g.: P, P → P1, P2 ; R, R → R1, R2 (multiple probes/primers)
            id_count = Counter((s, p) for s, p, _ in collected)
            id_seen  = defaultdict(int)
            for s_id, p_id, p_seq in collected:
                key = (s_id, p_id)
                if id_count[key] > 1:
                    id_seen[key] += 1
                    final_id = f"{p_id}{id_seen[key]}"
                else:
                    final_id = p_id
                primer_sets[s_id][final_id] = p_seq

    except Exception as e:
        print(f"Error reading primers file / Erreur lecture fichier amorces : {e}")
        sys.exit(1)

    return primer_sets

def auto_split_fip_bip(primer_sets, targets, txt):
    """
    Tente de séparer automatiquement les amorces nommées FIP ou BIP
    en trouvant l'alignement de leurs sous-parties sur les premières séquences cibles.
    """
    # Prendre jusqu'à 200 cibles pour tester (au cas où les premières auraient des mutations divergentes)
    test_targets = list(targets.values())[:200]
    if not test_targets: return
    
    for set_id, primers in list(primer_sets.items()):
        for combo_name in ['FIP', 'BIP']:
            if combo_name in primers:
                seq_combo = primers[combo_name]
                found_split = False
                
                # On cherche d'abord avec 0 erreur, puis 1, puis 2
                for allowed_err in [0, 1, 2]:
                    if found_split: break
                    # Ensuite on teste l'absence de linker (0), puis les linkers biologiquement valides (4 à 9 bases)
                    # Next we test the absence of linker (0), then biologically valid linkers (4 to 9 bases)
                    for linker_len in [0, 4, 5, 6, 7, 8, 9]:
                        if found_split: break
                        for t_seq in test_targets:
                            if found_split: break
                            
                            valid_splits = []
                            L = len(seq_combo)
                            for i in range(14, L - 13):
                                if i + linker_len > L - 14:
                                    continue
                                
                                if linker_len > 0:
                                    linker_seq = seq_combo[i:i+linker_len]
                                    # Contrainte biologique: un linker doit être une répétition de base (ex: TTTT, AAAA)
                                    if len(set(linker_seq)) != 1:
                                        continue
                                
                                part1 = seq_combo[:i]
                                part2 = seq_combo[i+linker_len:]
                                
                                pos1 = primer_matches_sequence(t_seq, part1, max_errors=allowed_err, strict_3prime_len=0)
                                pos2 = primer_matches_sequence(t_seq, part2, max_errors=allowed_err, strict_3prime_len=0)
                                
                                if pos1 and pos2:
                                    dist = abs(pos2[0] - pos1[0])
                                    if dist < 250:
                                        if combo_name == 'FIP' and pos2[0] < pos1[0]:
                                            valid_splits.append((part1, part2))
                                        elif combo_name == 'BIP' and pos1[0] < pos2[0]:
                                            valid_splits.append((part1, part2))
                        
                            if valid_splits:
                                # En cas d'ambiguïté (plusieurs coupes valides avec 0 erreur), 
                                # on choisit la coupe qui équilibre le mieux la taille des deux amorces.
                                best_split = min(valid_splits, key=lambda x: abs(len(x[0]) - len(x[1])))
                                part1, part2 = best_split
                                
                                if combo_name == 'FIP':
                                    primers['F1'] = str(Seq(part1).reverse_complement())
                                    primers['F2'] = part2
                                else:
                                    primers['B1'] = str(Seq(part1).reverse_complement())
                                    primers['B2'] = part2
                                    
                                found_split = True
                                print(txt['split_success'].format(set_id, combo_name, combo_name[0]+'1', combo_name[0]+'2', linker_len))
                
                if not found_split:
                    import re
                    # Fallback heuristique : si la souche compatible n'est pas dans les 200 premières,
                    # on tente une coupe "à l'aveugle" en cherchant au moins 4 "T" au milieu (méthode manuelle classique)
                    match = re.search(r'T{4,}', seq_combo[14:-14])
                    if match:
                        linker_start = match.start() + 14
                        part1 = seq_combo[:linker_start]
                        part2 = seq_combo[linker_start+4:]  # On ignore exactement 4 "T"
                        
                        if combo_name == 'FIP':
                            primers['F1'] = str(Seq(part1).reverse_complement())
                            primers['F2'] = part2
                        else:
                            primers['B1'] = str(Seq(part1).reverse_complement())
                            primers['B2'] = part2
                            
                        found_split = True
                        print(txt['split_success'].format(set_id, combo_name, combo_name[0]+'1', combo_name[0]+'2', '4 (fallback heuristique)'))
                
                if found_split:
                    del primers[combo_name]
                else:
                    print(txt['split_fail'].format(set_id, combo_name))

def process_primer_matches_chunk(args):
    """
    Fonction worker pour traiter un sous-ensemble (chunk) de séquences cibles
    contre tous les sets d'amorces en parallèle.
    """
    targets_chunk, primer_sets, errors, strict_3prime, strict_3prime_tolerate = args
    
    local_primer_matches = {}
    local_primer_positions = {}
    
    for seq_id, seq in targets_chunk:
        for set_id, primers in primer_sets.items():
            if set_id not in local_primer_matches:
                local_primer_matches[set_id] = {}
                local_primer_positions[set_id] = {}
            if seq_id not in local_primer_positions[set_id]:
                local_primer_positions[set_id][seq_id] = {}
                
            for primer_id, primer_seq in primers.items():
                pos = primer_matches_sequence(seq, primer_seq, errors, strict_3prime, strict_3prime_tolerate)
                if pos:
                    if primer_id not in local_primer_matches[set_id]:
                        local_primer_matches[set_id][primer_id] = set()
                    local_primer_matches[set_id][primer_id].add(seq_id)
                    local_primer_positions[set_id][seq_id][primer_id] = pos
                    
    return local_primer_matches, local_primer_positions, len(targets_chunk) * sum(len(p) for p in primer_sets.values())

def main():
    parser = argparse.ArgumentParser(description="Évalue la couverture et l'ordre des amorces LAMP / Evaluate LAMP primer coverage and structural order.")
    parser.add_argument("-t", "--target", required=True, help="Fichier FASTA cible / Target FASTA file.")
    parser.add_argument("-p", "--primers", required=True, help="Fichier FASTA ou TXT des amorces / Primers FASTA or TXT file.")
    parser.add_argument("-o", "--output", required=True, help="Fichier de rapport en sortie / Output report file.")
    parser.add_argument("-e", "--errors", type=int, default=0, help="Nombre max d'erreurs hors zone 3' / Max errors outside 3' region. Def: 0")
    parser.add_argument("-s", "--strict-3prime", type=int, default=3, dest="strict_3prime", help="Taille zone 3' stricte / Strict 3' region size. Def: 3")
    parser.add_argument("--strict-3prime-tolerate", type=int, choices=[0, 1, 2], default=0, help="Niveau de tolérance en zone 3' (0: tout strict, 1: pos 2 tolérée, 2: pos 1 et 2 tolérées). / Tolerance level in the 3' region (0: all strict, 1: pos 2 tolerated, 2: pos 1 and 2 tolerated).")
    parser.add_argument("--strict-intersection", action="store_true", help="Exige que toutes les amorces du fichier matchent la cible. / Requires all primers in the file to match the target.")
    parser.add_argument("--max-n-pct", type=float, default=0, dest="max_n_pct",
        help="Exclure les séquences dont le pourcentage de N dépasse ce seuil / Exclude sequences with N percentage above this threshold. 0=désactivé/disabled (défaut/default). Def: 0")
    
    # Options de sortie
    parser.add_argument("--summary-only", action="store_true", help="N'affiche que les statistiques / Output only summary statistics.")
    parser.add_argument("--combine", action="store_true", help="Couverture combinatoire 2 à 2 / Calculate 2-by-2 multiplexing coverage.")
    parser.add_argument("--export-seqs", action="store_true", help="Exporte les séquences validées / Export validated sequences per set.")
    parser.add_argument("--pcr", action="store_true", help="Mode PCR : Gère les amorces Forward, Reverse et Sonde / PCR mode: handles Fwd, Rev and Probe primers.")
    parser.add_argument("--diagnose-nonmatch", action="store_true", dest="diagnose_nonmatch",
        help="Diagnostique les séquences non-matchées : qualité insuffisante vs vrai non-match. / Diagnose non-matching sequences: poor quality vs true non-match.")
    parser.add_argument("--diag-n-pct", type=float, default=5.0, dest="diag_n_pct",
        help="Seuil de %% de N dans l'amplicon pour qualifier une séquence de mauvaise qualité (utilisé avec --diagnose-nonmatch). / N%% threshold in the amplicon to flag a sequence as poor quality (used with --diagnose-nonmatch). Def: 5.0")
    
    # Langue
    parser.add_argument("--lng", type=str, default="en", choices=["en", "fr"], help="Langue du rapport / Report language (en, fr). Def: en")
    
    args = parser.parse_args()
    
    # Dictionnaire de traduction
    lang = args.lng.lower()
    T = {
        'fr': {
            'target_load': "Chargement des séquences cibles...",
            'target_err': "Erreur : Aucune séquence cible trouvée.",
            'target_loaded': "séquences cibles chargées.",
            'primer_load': "Chargement des amorces...",
            'primer_err': "Erreur : Aucune amorce trouvée.",
            'primer_loaded': "set(s) d'amorces identifié(s).",
            'analysis': "Analyse en cours avec tolérance = {} erreur(s) (et {} nt stricts en 3')...",
            'progression': "Progression : {}%",
            'report_gen': "Génération du rapport...",
            'report_title': "Rapport de Couverture des Amorces LAMP",
            'target_file': "Fichier cible",
            'primer_file': "Fichier amorces",
            'tolerance': "Tolérance d'erreurs par amorce",
            'strict_3': "avec {} nt stricts en 3'",
            'active_options': "Options actives",
            'set_title': "--- SET D'AMORCES : {} ---",
            'indiv_match': "Matchs individuels par amorce :",
            'global_raw': "Match Global du Set (Intersection Brute, toutes amorces présentes)",
            'global_base': "Match de Base du Set (Intersection Validation : amorces essentielles uniquement)",
            'global_valid': "Match Global Valide (Intersection Base + Ordre Correct structurel LAMP)",
            'excluded_label': "Séquences de mauvaise qualité exclues (>{}% de N)",
            'total_analysed': "Total de séquences analysées",
            'amplified_seqs': "Séquences amplifiées théoriquement par le Set {} :",
            'table_header': "Séquence_ID\tTaille_Amplicon\tStatut_Ordre\tOrdre_Observe",
            'order_correct': "Ordre Correct",
            'order_incorrect': "Ordre Incorrect",
            'no_seq': "(Aucune séquence ne remplit les critères d'amplification)",
            'export_success': "-> Fichier d'export généré avec les séquences valides : {}",
            'export_err': "-> Erreur lors de l'export des séquences : {}",
            'combine_title': "Calcul Combinatoire (Multiplexage 2 à 2)",
            'combine_subtitle': "Couverture théorique si les sets sont combinés dans la même réaction (Union des séquences valides) :",
            'done': "Terminé. Rapport principal sauvegardé dans : {}",
            'write_err': "Erreur lors de l'écriture du rapport : {}",
            'seqs_word': "séquences",
            'split_success': "Auto-split réussi pour Set {} ({}) -> {} et {} (linker de {} nt ignoré).",
            'split_fail': "Attention : Impossible de déterminer la coupe automatique pour Set {} ({}).",
            'type_union_label':  "Union",
            'marginal_title':     "Valeur ajoutée par version supplémentaire (par type) :",
            'pool_title':         "Meilleures combinaisons selon le nombre total d'amorces (pool size) :",
            'pool_tested':        "combinaison(s) de sous-ensembles testée(s)",
            'pool_exhaustive':    "recherche exhaustive",
            'pool_greedy':        "algorithme glouton progressif",
            'diag_title':         "Diagnostic des séquences non-matchées (méthode : {})",
            'diag_ref_amp':       "Amplicon de référence : {} nt (séquence de réf : {})",
            'diag_too_variable':  "🚫 Séquences trop variables (au moins 1 amorce se fixe, divergence réelle)",
            'diag_truncated_poor_divergent': "❓ Séquences tronquées, de mauvaise qualité ou très divergentes (0 amorce se fixe)",
            'of_non_matched':     "des non-matchés",
            'too_variable_list_title': "Liste des séquences trop variables :",
            'truncated_list_title': "Liste des séquences tronquées, de mauvaise qualité ou très divergentes :",
            'too_many_errors':    "Trop d'erreurs",
            'not_found_word':     "Introuvable",
            'poor_quality_word':  "Mauvaise qualité",
            'matched_word':       "Amorce(s) fixée(s)",
            'unmatched_word':     "Non-fixée(s)",
            'diag_potential_cov': "→ Couverture potentielle si séquences tronquées, de mauvaise qualité ou très divergentes exclues",
        },
        'en': {
            'target_load': "Loading target sequences...",
            'target_err': "Error: No target sequences found.",
            'target_loaded': "target sequences loaded.",
            'primer_load': "Loading primers...",
            'primer_err': "Error: No primers found.",
            'primer_loaded': "primer set(s) identified.",
            'analysis': "Analysis in progress with tolerance = {} error(s) (and {} strict 3' nt)...",
            'progression': "Progress : {}%",
            'report_gen': "Generating report...",
            'report_title': "LAMP Primer Coverage Report",
            'target_file': "Target file",
            'primer_file': "Primer file",
            'tolerance': "Error tolerance per primer",
            'strict_3': "with {} strict 3' nt",
            'active_options': "Active options",
            'set_title': "--- PRIMER SET : {} ---",
            'indiv_match': "Individual matches per primer:",
            'global_raw': "Set Global Match (Raw Intersection, all primers present)",
            'global_base': "Set Base Match (Validation Intersection: essential primers only)",
            'global_valid': "Set Valid Global Match (Base Intersection + Structurally Correct LAMP Order)",
            'excluded_label': "Bad quality sequences excluded (>{}% N)",
            'total_analysed': "Total sequences analysed",
            'amplified_seqs': "Theoretically amplified sequences by Set {} :",
            'table_header': "Sequence_ID\tAmplicon_Size\tOrder_Status\tObserved_Order",
            'order_correct': "Correct Order",
            'order_incorrect': "Incorrect Order",
            'no_seq': "(No sequences meet the amplification criteria)",
            'export_success': "-> Export file generated with valid sequences: {}",
            'export_err': "-> Error exporting sequences: {}",
            'combine_title': "Combinatorial Calculation (2-by-2 Multiplexing)",
            'combine_subtitle': "Theoretical coverage if sets are combined in the same reaction (Union of valid sequences):",
            'done': "Done. Main report saved in: {}",
            'write_err': "Error writing report: {}",
            'seqs_word': "sequences",
            'split_success': "Auto-split successful for Set {} ({}) -> {} and {} ({} nt linker ignored).",
            'split_fail': "Warning: Could not automatically determine split for Set {} ({}).",
            'type_union_label':  "Union",
            'marginal_title':     "Marginal value added per additional version:",
            'pool_title':         "Best combinations by total number of primers (pool size):",
            'pool_tested':        "subset combination(s) tested",
            'pool_exhaustive':    "exhaustive search",
            'pool_greedy':        "progressive greedy algorithm",
            'diag_title':         "Non-matching sequences diagnosis (method: {})",
            'diag_ref_amp':       "Reference amplicon: {} nt (reference seq: {})",
            'diag_too_variable':  "🚫 Too variable sequences (at least 1 primer binds, real divergence)",
            'diag_truncated_poor_divergent': "❓ Truncated, poor quality or highly divergent (0 primers bind)",
            'of_non_matched':     "of non-matched",
            'too_variable_list_title': "List of too variable sequences:",
            'truncated_list_title': "List of truncated, poor or highly divergent sequences:",
            'too_many_errors':    "Too many errors",
            'not_found_word':     "Not found",
            'poor_quality_word':  "Poor quality",
            'matched_word':       "Matched primer(s)",
            'unmatched_word':     "Unmatched",
            'diag_potential_cov': "→ Potential coverage if truncated, poor quality or highly divergent sequences excluded",
        }
    }
    
    txt = T[lang]
    if args.pcr:
        if lang == 'fr':
            txt['report_title'] = "Rapport de Couverture des Amorces PCR"
            txt['global_base']  = "Couverture Validation PCR (F∪ ∩ R∪ ∩ Sondes∪)"
            txt['global_valid'] = "Match Global Valide (Validation + Ordre PCR Correct)"
        else:
            txt['report_title'] = "PCR Primer Coverage Report"
            txt['global_base']  = "PCR Validation Coverage (F∪ ∩ R∪ ∩ Probes∪)"
            txt['global_valid'] = "Set Valid Global Match (Validation + Correct PCR Order)"
    
    print(txt['target_load'])
    targets = {}
    try:
        # Chargement avec barre de progression / Loading with progress bar
        all_records = list(tqdm(
            SeqIO.parse(args.target, "fasta"),
            desc="📂 Chargement / Loading" if args.lng == 'fr' else "📂 Loading sequences",
            unit=" seq",
            colour="cyan"
        ))
        for record in all_records:
            raw_seq   = str(record.seq).upper()
            clean_seq = raw_seq.replace('-', '')
            if len(clean_seq) > 100:
                targets[record.description] = clean_seq   # Sans gaps / Without gaps (for matching)
        
        gapped_targets = {record.description: str(record.seq).upper() for record in all_records}
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    total_targets = len(targets)
    if total_targets == 0:
        print(txt['target_err'])
        sys.exit(1)
    print(f"{total_targets} {txt['target_loaded']}")

    # Pré-filtre qualité : exclut les séquences dont le % de N dépasse le seuil
    # Quality pre-filter: exclude sequences whose N percentage exceeds the threshold
    n_excluded = 0
    if args.max_n_pct > 0:
        bad_qual = {sid for sid, seq in targets.items()
                    if seq.count('N') / len(seq) * 100 > args.max_n_pct}
        n_excluded = len(bad_qual)
        if n_excluded > 0:
            targets = {sid: seq for sid, seq in targets.items() if sid not in bad_qual}
            if lang == 'fr':
                print(f"  ⚠️  {n_excluded} séquence(s) exclues (>  {args.max_n_pct}% de N dans la séquence).")
            else:
                print(f"  ⚠️  {n_excluded} sequence(s) excluded (> {args.max_n_pct}% N in sequence).")
        total_targets = len(targets)
        if total_targets == 0:
            print(txt['target_err'])
            sys.exit(1)
    
    print(txt['primer_load'])
    primer_sets = load_primers(args.primers, is_pcr=args.pcr)
    if not primer_sets:
        print(txt['primer_err'])
        sys.exit(1)
    
    suspected_pcr = False
    if not args.pcr:
        # Tentative de séparation auto de FIP et BIP
        auto_split_fip_bip(primer_sets, targets, txt)
        
        # Vérification si un essai PCR est suspecté en mode LAMP / Checking if a PCR assay is suspected in LAMP mode
        for s_id, primers in primer_sets.items():
            set_keys = {k.upper() for k in primers.keys()}
            if 'FP' in set_keys or 'RP' in set_keys or 'PROBE' in set_keys or (set_keys.issubset({'F', 'R', 'P'}) and not set_keys.intersection({'F3', 'B3', 'F2', 'F1', 'B1', 'B2', 'FIP', 'BIP'})):
                suspected_pcr = True
                break
                
        if suspected_pcr:
            print("\n" + "!" * 70)
            if args.lng == 'fr':
                print("⚠️  ATTENTION : Les amorces ressemblent à un essai PCR (FP, RP, Probe).")
                print("   Le script s'exécute actuellement en mode LAMP.")
                print("   Veuillez relancer avec l'option '--pcr' pour une analyse correcte.")
            else:
                print("⚠️  WARNING: The primers look like a PCR assay (FP, RP, Probe).")
                print("   The script is currently running in LAMP mode.")
                print("   Please rerun with the '--pcr' option for a correct analysis.")
            print("!" * 70 + "\n")
    
    print(f"{len(primer_sets)} {txt['primer_loaded']}")
    print(txt['analysis'].format(args.errors, args.strict_3prime))
    
    primer_matches = defaultdict(lambda: defaultdict(set))
    primer_positions = defaultdict(lambda: defaultdict(dict))
    bad_seqs_per_set = defaultdict(set)  # inutilisé ici mais conservé pour compatibilité / kept for compatibility

    # Boucle principale d'analyse / Main analysis loop
    total_steps = len(targets) * sum(len(p) for p in primer_sets.values())
    bar_label = "🔬 Analyse" if args.lng == 'fr' else "🔬 Analysing"

    workers = min(8, max(1, (os.cpu_count() or 4) - 2))
    print(f"\n🚀 Utilisation de {workers} processus pour l'analyse / Using {workers} processes for analysis\n")

    targets_items = list(targets.items())
    
    total_primers = sum(len(p) for p in primer_sets.values())
    
    # Calcul dynamique de la taille de chunk pour une barre de progression fluide
    # On veut suffisamment de chunks pour alimenter les workers, mais pas trop gros 
    # pour que l'affichage se mette à jour fréquemment (idéalement ~1000 comparaisons par chunk max).
    base_chunk = math.ceil(len(targets_items) / (workers * 4))
    max_chunk = max(1, 1000 // (total_primers if total_primers > 0 else 1))
    chunk_size = min(base_chunk, max_chunk)
    chunk_size = max(1, chunk_size)
    
    chunks = []
    for i in range(0, len(targets_items), chunk_size):
        chunk = targets_items[i:i+chunk_size]
        chunks.append((chunk, primer_sets, args.errors, args.strict_3prime, args.strict_3prime_tolerate))
        
    with tqdm(total=total_steps, desc=bar_label, unit=" seq", colour="green") as pbar:
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
            futures = [executor.submit(process_primer_matches_chunk, chunk) for chunk in chunks]
            
            for future in concurrent.futures.as_completed(futures):
                local_matches, local_positions, steps_done = future.result()
                
                # Fusionner les résultats / Merge results
                for set_id, primers_dict in local_matches.items():
                    for primer_id, seq_ids in primers_dict.items():
                        primer_matches[set_id][primer_id].update(seq_ids)
                        
                for set_id, seq_dict in local_positions.items():
                    for seq_id, primer_dict in seq_dict.items():
                        primer_positions[set_id][seq_id].update(primer_dict)
                        
                pbar.update(steps_done)

    print(txt['report_gen'])
    
    if args.pcr:
        MASTER_ORDER = ['F', 'P', 'R']
    else:
        # Ordre général LAMP incluant à la fois les LOOP et les STEM (s'ils sont là)
        MASTER_ORDER = ['F3', 'F2', 'FLOOP', 'F1', 'STEMF', 'STEMB', 'B1', 'BLOOP', 'B2', 'B3']
    MASTER_ORDER_RC = list(reversed(MASTER_ORDER))
    
    # Définition des amorces essentielles (utilisées pour l'intersection relaxée)
    # Definition of essential primers (used for relaxed intersection)
    ESSENTIAL_LAMP = {'F3', 'B3', 'F2', 'F1', 'B1', 'B2', 'FIP', 'BIP'}
    ESSENTIAL_PCR = {'F', 'R'}
    
    # Stockage des sets de séquences valides pour le calcul combinatoire
    valid_sequences_per_set = {}
    
    try:
        with open(args.output, 'w') as out:
            if not args.pcr and suspected_pcr:
                if args.lng == 'fr':
                    out.write("⚠️  ATTENTION : Détecté comme essai PCR mais exécuté en mode LAMP. Veuillez utiliser --pcr.\n")
                else:
                    out.write("⚠️  WARNING: Detected as PCR assay but executed in LAMP mode. Please use --pcr.\n")
                out.write("=" * 80 + "\n\n")
                
            out.write(f"{txt['report_title']}\n")
            out.write("=" * max(len(txt['report_title']), 38) + "\n")
            out.write(f"{txt['target_file']} : {args.target} ({total_targets} {txt['seqs_word']})\n")
            out.write(f"{txt['primer_file']} : {args.primers}\n")
            out.write(f"{txt['tolerance']} : {args.errors} ({txt['strict_3'].format(args.strict_3prime)})\n")
            
            # Afficher les options actives
            active_options = []
            if args.summary_only: active_options.append("--summary-only")
            if args.combine: active_options.append("--combine")
            if args.export_seqs: active_options.append("--export-seqs")
            if args.strict_intersection: active_options.append("--strict-intersection")
            if active_options:
                out.write(f"{txt['active_options']} : {', '.join(active_options)}\n")
            out.write("\n")
            
            for set_id, primers in primer_sets.items():
                out.write(txt['set_title'].format(set_id) + "\n")

                # Nombre de séquences exclues (global) / Number of excluded sequences (global)
                effective_targets = total_targets  # ici, exclusion globale / global exclusion here

                # Affichage des exclusions de mauvaise qualité / Display of quality exclusions
                if args.max_n_pct > 0:
                    out.write(f"{txt['excluded_label'].format(args.max_n_pct)} : {n_excluded}\n")
                    out.write(f"{txt['total_analysed']} : {effective_targets}\n")

                set_matches_list = []
                out.write(f"{txt['indiv_match']}\n")
                denom = effective_targets if effective_targets > 0 else 1
                for primer_id, primer_seq in primers.items():
                    matches = primer_matches[set_id][primer_id]
                    set_matches_list.append(matches)
                    match_pct = (len(matches) / denom) * 100
                    out.write(f"  - {primer_id} : {match_pct:.2f}% ({len(matches)}/{effective_targets})\n")
                
                # ── Groupement par type de base (multi-versions) ───────────────────────
                # Group by base type (multi-version support)
                # Ex: {F3: [F31, F32], B3: [B3], FIP: [FIP1, FIP2, FIP3]}
                type_groups = group_primer_versions(primers, args.pcr)

                # ── Union de chaque type d'amorce ─────────────────────────────────────
                # Union per primer type: a sequence is covered if ≥1 version of that type detects it
                type_unions = {}
                for base_type, versions in type_groups.items():
                    union_seqs = set.union(*[primer_matches[set_id][v] for v in versions if primer_matches[set_id].get(v) is not None])
                    type_unions[base_type] = union_seqs

                    # Afficher l'union si plusieurs versions existent pour ce type
                    # Display the union if multiple versions exist for this type
                    if len(versions) > 1:
                        u_pct = (len(union_seqs) / denom) * 100
                        names = '+'.join(versions)
                        out.write(f"  \u2192 {txt['type_union_label']} {base_type} ({names}) : {u_pct:.2f}% ({len(union_seqs)}/{effective_targets})\n")

                # ── Intersection brute (tous les types d'amorces présents matchent) ────
                # Raw intersection (all present primer types must match, taking the union of variants per type)
                all_unions = list(type_unions.values())
                intersection_matches = set.intersection(*all_unions) if all_unions else set()

                # ── Couverture poolée (union par type, intersection des types) ─────────
                # Pooled coverage (union per type, then intersect across types)
                if args.pcr:
                    # PCR : types essentiels = F, R (P optionnel mais inclus s'il existe)
                    # PCR: essential types = F, R (P optional but included if present)
                    essential_types = [t for t in type_groups if t in PCR_CANONICAL_TYPES]
                    if not essential_types:
                        essential_types = list(type_groups.keys())
                else:
                    # LAMP : types essentiels = {F3,B3,F2,F1,B1,B2,FIP,BIP} ∩ types présents
                    # LAMP: essential types = {F3,B3,F2,F1,B1,B2,FIP,BIP} ∩ present types
                    essential_types = [t for t in type_groups if t in ESSENTIAL_LAMP]
                    if not essential_types:
                        essential_types = list(type_groups.keys())

                essential_union_list = [type_unions[t] for t in essential_types]

                if args.strict_intersection:
                    validation_matches = intersection_matches
                else:
                    validation_matches = set.intersection(*essential_union_list) if essential_union_list else set()

                # ── Vérification de l'ordre et calcul de la taille de l'amplicon ──────
                # Order verification and amplicon size calculation
                valid_order_matches = []
                seq_details = []

                for seq_id in validation_matches:
                    positions = primer_positions[set_id][seq_id]

                    # Tri des noms d'amorces selon la coordonnée de départ
                    # Sort primers by their start position on the target
                    sorted_primers = sorted(positions.keys(), key=lambda p: positions[p][0])

                    if args.pcr:
                        # Remapper les versions d'amorces PCR vers leur type de base (F, R, P)
                        # Remap versioned PCR primers to their base type (F, R, P)
                        base_type_pos = defaultdict(list)
                        for p, pos in positions.items():
                            bt = get_base_type(p, is_pcr=True)
                            base_type_pos[bt].append(pos)
                        
                        if 'F' not in base_type_pos or 'R' not in base_type_pos:
                            is_correct_order = False
                        else:
                            f_starts = [pos[0] for pos in base_type_pos['F']]
                            f_ends = [pos[1] for pos in base_type_pos['F']]
                            r_starts = [pos[0] for pos in base_type_pos['R']]
                            r_ends = [pos[1] for pos in base_type_pos['R']]
                            
                            min_f_start, max_f_end = min(f_starts), max(f_ends)
                            min_r_start, max_r_end = min(r_starts), max(r_ends)
                            
                            if 'P' in base_type_pos:
                                p_starts = [pos[0] for pos in base_type_pos['P']]
                                p_ends = [pos[1] for pos in base_type_pos['P']]
                                
                                # Sens : F < P < R
                                is_sense = (min_f_start < min(p_starts)) and (max(p_ends) < max_r_end)
                                # Anti-sens : R < P < F
                                is_antisense = (min_r_start < min(p_starts)) and (max(p_ends) < max_f_end)
                                is_correct_order = is_sense or is_antisense
                            else:
                                is_correct_order = (min_f_start < max_r_end) or (min_r_start < max_f_end)
                    else:
                        # LAMP : remapper les versions vers leur type de base pour la vérification d'ordre
                        # LAMP: remap versioned primers to their base type for order check
                        base_type_pos = {}
                        for p, pos in positions.items():
                            bt = get_base_type(p, is_pcr=False)
                            # Si plusieurs versions matchent, garder la position la plus en 5'
                            # If multiple versions match, keep the most 5' position
                            if bt not in base_type_pos or pos[0] < base_type_pos[bt][0]:
                                base_type_pos[bt] = pos
                        sorted_base_types = sorted(base_type_pos.keys(), key=lambda t: base_type_pos[t][0])
                        expected_sense = [p for p in MASTER_ORDER if p in base_type_pos]
                        expected_anti  = [p for p in MASTER_ORDER_RC if p in base_type_pos]
                        is_correct_order = (sorted_base_types == expected_sense) or (sorted_base_types == expected_anti)

                    if is_correct_order:
                        valid_order_matches.append(seq_id)

                    starts = [pos[0] for pos in positions.values()]
                    ends   = [pos[1] for pos in positions.values()]
                    amplicon_size = max(ends) - min(starts)

                    status = txt['order_correct'] if is_correct_order else txt['order_incorrect']
                    observed_order = "-".join(sorted_primers)
                    seq_details.append((seq_id, amplicon_size, status, observed_order))

                # Sauvegarde pour combine
                valid_sequences_per_set[set_id] = set(valid_order_matches)

                # ── Calcul des pourcentages ────────────────────────────────────────────
                # Percentage calculation
                raw_match_pct   = (len(intersection_matches) / denom) * 100
                val_match_pct   = (len(validation_matches)   / denom) * 100
                valid_match_pct = (len(valid_order_matches)  / denom) * 100 if denom > 0 else 0

                # ── Affichage statistiques ────────────────────────────────────────────
                # Statistics display
                has_multi_versions = any(len(vs) > 1 for vs in type_groups.values())

                out.write(f"\n{txt['global_raw']} : {raw_match_pct:.2f}% ({len(intersection_matches)}/{effective_targets})\n")
                # Couverture poolée (union par type) — affichée si différente du brut
                # Pooled coverage (union per type) — shown if different from raw
                if validation_matches != intersection_matches:
                    out.write(f"{txt['global_base']} : {val_match_pct:.2f}% ({len(validation_matches)}/{effective_targets})\n")
                out.write(f"{txt['global_valid']} : {valid_match_pct:.2f}% ({len(valid_order_matches)}/{effective_targets})\n")

                # ── Diagnostic des séquences non-matchées (si --diagnose-nonmatch) ──────
                # Non-matching sequences diagnosis (if --diagnose-nonmatch)
                if args.diagnose_nonmatch:
                    all_seq_ids   = set(targets.keys())
                    non_matched   = list(all_seq_ids - validation_matches)
                    n_non_matched = len(non_matched)
                    
                    if n_non_matched > 0:
                        # Trouver la ref : séquence avec ordre correct / Find ref: sequence with correct order
                        ref_id = valid_order_matches[0] if valid_order_matches else None
                        
                        diag = diagnose_nonmatching_sequences(
                            non_matched_ids    = non_matched,
                            targets            = targets,
                            gapped_targets     = gapped_targets,
                            primers_dict       = primers,
                            primer_positions   = primer_positions,
                            set_id             = set_id,
                            valid_order_matches= valid_order_matches,
                            is_pcr             = args.pcr,
                            max_n_pct_diag     = args.diag_n_pct,
                            primer_matches_set = primer_matches[set_id]
                        )
                        
                        out.write(f"\n{'─'*40}\n")
                        out.write(f"{txt['diag_title'].format(diag.get('method', 'Primer Binding Analysis'))} ({n_non_matched} {txt['seqs_word']}) :\n")
                        
                        # Ligne amplicon de référence si disponible
                        # Reference amplicon line if available
                        if diag['ref_amp_len'] > 0 and ref_id:
                            out.write(f"  {txt['diag_ref_amp'].format(diag['ref_amp_len'], ref_id[:60])}\n")
                        
                        # Afficher les catégories de diagnostic basées sur le nombre de liaisons d'amorces
                        # Display diagnostic categories based on primer binding counts
                        tv  = diag['too_variable']
                        tpd = diag['truncated_poor_divergent']
                        
                        tv_pct  = (len(tv) / n_non_matched) * 100 if n_non_matched > 0 else 0
                        tpd_pct = (len(tpd) / n_non_matched) * 100 if n_non_matched > 0 else 0
                        
                        out.write(f"  {txt['diag_too_variable']} : {len(tv)} ({tv_pct:.1f}% {txt['of_non_matched']})\n")
                        out.write(f"  {txt['diag_truncated_poor_divergent']} : {len(tpd)} ({tpd_pct:.1f}% {txt['of_non_matched']})\n")
                        
                        if len(tpd) > 0:
                            corrected_total = effective_targets - len(tpd)
                            if corrected_total > 0:
                                potential_pct = (len(valid_order_matches) / corrected_total) * 100
                                out.write(f"  {txt['diag_potential_cov']} : {potential_pct:.2f}% ({len(valid_order_matches)}/{corrected_total})\n")
                        
                        # Si l'option --summary-only n'est pas active, lister les IDs détaillés pour chaque catégorie
                        # If --summary-only is not active, list detailed IDs for each category
                        if not args.summary_only:
                            details = diag['details']
                            if tv:
                                out.write(f"\n  {txt['too_variable_list_title']}\n")
                                for seq_id in sorted(tv):
                                    det = details[seq_id]
                                    matched_str = ", ".join(det['matched']) if det['matched'] else "-"
                                    unmatched_parts = []
                                    for p_id, status in sorted(det['unmatched'].items()):
                                        lbl = txt.get(status + '_word', status) if status in ['not_found', 'poor_quality'] else txt.get('too_many_errors', 'Too many errors')
                                        unmatched_parts.append(f"{p_id} ({lbl})")
                                    unmatched_str = ", ".join(unmatched_parts) if unmatched_parts else "-"
                                    out.write(f"    - {seq_id}\n")
                                    out.write(f"      * {txt['matched_word']} : {matched_str}\n")
                                    out.write(f"      * {txt['unmatched_word']} : {unmatched_str}\n")
                            if tpd:
                                out.write(f"\n  {txt['truncated_list_title']}\n")
                                for seq_id in sorted(tpd):
                                    det = details[seq_id]
                                    matched_str = ", ".join(det['matched']) if det['matched'] else "-"
                                    unmatched_parts = []
                                    for p_id, status in sorted(det['unmatched'].items()):
                                        lbl = txt.get(status + '_word', status) if status in ['not_found', 'poor_quality'] else txt.get('too_many_errors', 'Too many errors')
                                        unmatched_parts.append(f"{p_id} ({lbl})")
                                    unmatched_str = ", ".join(unmatched_parts) if unmatched_parts else "-"
                                    out.write(f"    - {seq_id}\n")
                                    if det['matched']:
                                        out.write(f"      * {txt['matched_word']} : {matched_str}\n")
                                    out.write(f"      * {txt['unmatched_word']} : {unmatched_str}\n")
                                    
                        out.write(f"{'─'*40}\n")


                # ── Analyse combinatoire (uniquement si multi-versions) ───────────────
                # Combinatorial analysis (only if multi-version types exist)
                if has_multi_versions:
                    # 1. Progression marginale par type
                    marginal_data = analyze_marginal_value(type_groups, dict(primer_matches[set_id]))
                    if marginal_data:
                        out.write(f"\n{txt['marginal_title']}\n")
                        for base_type, steps in sorted(marginal_data.items()):
                            out.write(f"  Type {base_type} :\n")
                            prev_cov_len = 0
                            for k, subset, cov_set in steps:
                                pct = (len(cov_set) / denom) * 100
                                gain_pct = ((len(cov_set) - prev_cov_len) / denom) * 100 if k > 1 else 0
                                gain_str = f" (+{gain_pct:.2f}%)" if k > 1 else ""
                                versions_str = "+".join(subset)
                                out.write(f"    {k} version(s) : {versions_str} \u2192 {pct:.2f}% ({len(cov_set)}/{effective_targets}){gain_str}\n")
                                prev_cov_len = len(cov_set)

                    # 2. Combinaisons optimales par taille de pool
                    best_by_pool, n_combos = find_best_subsets_by_pool(
                        type_groups, dict(primer_matches[set_id])
                    )
                    out.write(f"\n{txt['pool_title']}\n")
                    algo_label = txt['pool_exhaustive'] if n_combos <= 50_000 else txt['pool_greedy']
                    out.write(f"  ({n_combos} {txt['pool_tested']}, {algo_label}) :\n")
                    
                    prev_cov_len = 0
                    for m in sorted(best_by_pool.keys()):
                        combo_dict, cov_set, step_algo = best_by_pool[m]
                        pct = (len(cov_set) / denom) * 100
                        gain_pct = ((len(cov_set) - prev_cov_len) / denom) * 100 if m > len(type_groups) else 0
                        gain_str = f" (+{gain_pct:.2f}%)" if m > len(type_groups) else ""
                        
                        combo_parts = []
                        for bt in sorted(combo_dict.keys()):
                            versions = combo_dict[bt]
                            if len(versions) > 1:
                                combo_parts.append(f"{bt}=[{'+'.join(versions)}]")
                            else:
                                combo_parts.append(f"{bt}={versions[0]}")
                        combo_str = " + ".join(combo_parts)
                        
                        out.write(f"  - Pool de {m} amorce(s) : {combo_str} \u2192 {pct:.2f}% ({len(cov_set)}/{effective_targets}){gain_str}\n")
                        prev_cov_len = len(cov_set)



                
                # 1. Option : Ne pas afficher les séquences (si --summary-only)
                if not args.summary_only:
                    out.write(f"\n{txt['amplified_seqs'].format(set_id)}\n")
                    out.write(f"{txt['table_header']}\n")
                    
                    if seq_details:
                        seq_details.sort(key=lambda x: x[0])
                        for detail in seq_details:
                            out.write(f"{detail[0]}\t{detail[1]}\t{detail[2]}\t{detail[3]}\n")
                    else:
                        out.write(f"{txt['no_seq']}\n")
                
                # 2. Option : Exporter les séquences dans un fichier propre (--export-seqs)
                if args.export_seqs:
                    base_name, ext = os.path.splitext(args.output)
                    export_filename = f"{base_name}_Set{set_id}_seqs.txt"
                    try:
                        with open(export_filename, 'w') as f_export:
                            sorted_valid = sorted(valid_order_matches)
                            for s in sorted_valid:
                                f_export.write(f"{s}\n")
                        out.write("\n" + txt['export_success'].format(export_filename) + "\n")
                    except Exception as e:
                        out.write("\n" + txt['export_err'].format(e) + "\n")
                
                out.write("\n" + "="*40 + "\n\n")

            # 3. Option : Calcul Combinatoire (--combine)
            if args.combine and len(primer_sets) >= 2:
                out.write(f"{txt['combine_title']}\n")
                out.write("-" * len(txt['combine_title']) + "\n")
                out.write(f"{txt['combine_subtitle']}\n\n")
                
                set_keys = list(primer_sets.keys())
                combine_results = []
                for s1, s2 in itertools.combinations(set_keys, 2):
                    union_set = valid_sequences_per_set[s1] | valid_sequences_per_set[s2]
                    combine_pct = (len(union_set) / total_targets) * 100 if total_targets > 0 else 0
                    combine_results.append((len(union_set), combine_pct, s1, s2))
                
                # Tri décroissant selon le nombre de séquences trouvées
                combine_results.sort(key=lambda x: x[0], reverse=True)
                
                for union_len, combine_pct, s1, s2 in combine_results:
                    out.write(f"  - Set {s1} + Set {s2} : {combine_pct:.2f}% ({union_len}/{total_targets})\n")
                out.write("\n" + "="*40 + "\n\n")

        print(txt['done'].format(args.output))
        
    except Exception as e:
        print(txt['write_err'].format(e))

if __name__ == "__main__":
    main()
