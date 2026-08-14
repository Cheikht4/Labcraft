#!/usr/bin/env python3
import argparse
import sys
import re
import regex
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Import optionnel de tqdm pour les barres de progression / Optional tqdm import for progress bars
try:
    from tqdm import tqdm
except ImportError:
    # Fallback si tqdm n'est pas installé / Fallback if tqdm is not installed
    def tqdm(iterable, **kwargs):
        desc = kwargs.get('desc', '')
        if desc:
            print(f"{desc}...")
        return iterable
    tqdm.write = print

import concurrent.futures
import os

# Dictionnaire IUPAC vers expression régulière
# IUPAC dictionary to regular expression / Dictionnaire IUPAC vers expression régulière
IUPAC_DICT = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
}

def seq_to_regex(seq):
    """Convertit une séquence IUPAC en expression régulière.
    Convert an IUPAC sequence to a regular expression.
    """
    pattern = ""
    for char in seq.upper():
        if char.isalpha():
            pattern += IUPAC_DICT.get(char, char)
    return pattern

def find_best_match(target_seq, primer_seq, max_errors):
    """
    Cherche la meilleure correspondance de l'amorce dans la séquence cible.
    Retourne (start, end, match_string, errors) ou None.
    Finds the best match of the primer in the target sequence.
    Returns (start, end, match_string, errors) or None.
    """
    pattern_str = seq_to_regex(primer_seq)
    if not pattern_str:
        return None
    # (?e) active le mode flou (fuzzy) dans la bibliothèque regex
    # (?e) enables fuzzy matching mode in the regex library
    regex_pattern = f"(?e)({pattern_str}){{e<={max_errors}}}"
    
    matches = regex.finditer(regex_pattern, target_seq, regex.BESTMATCH)
    
    best_match = None
    min_errors = float('inf')
    
    for match in matches:
        counts = match.fuzzy_counts
        total_errors = sum(counts)
        
        if total_errors < min_errors:
            min_errors = total_errors
            best_match = match
            
    if best_match:
        return best_match.start(), best_match.end(), best_match.group(), min_errors
    
    return None

def primer_matches_sequence_simple(target_seq, primer_seq, max_errors):
    """
    Version simplifiée pour auto_split: cherche l'amorce sur le brin sens uniquement.
    Simplified version for auto_split: searches only on the sense strand.
    Retourne (start, end) ou None.
    """
    pattern_str = seq_to_regex(primer_seq)
    if not pattern_str:
        return None
    regex_pattern = f"(?e)({pattern_str}){{e<={max_errors}}}"
    m = regex.search(regex_pattern, target_seq, regex.BESTMATCH)
    if m:
        return m.span()
    return None

def get_ungapped_mapping(gapped_seq):
    """
    Renvoie la séquence sans trous et un tableau mappant l'index sans trous à l'index avec trous.
    Returns the ungapped sequence and an array mapping ungapped index to gapped index.
    """
    ungapped = []
    mapping = []
    for i, char in enumerate(gapped_seq):
        if char != '-':
            ungapped.append(char)
            mapping.append(i)
    return "".join(ungapped), mapping

def auto_split_fip_bip(primers_dict, ref_ungapped_str, max_errors=2):
    """
    Détecte et sépare automatiquement les amorces FIP et BIP en leurs sous-parties.
    Imite la logique de lamp_coverage.py : cherche la meilleure coupure en testant
    différentes longueurs de linker (0 à 9 nt) et différents niveaux d'erreur.
    
    Automatically detects and splits FIP and BIP primers into their sub-parts.
    Mirrors the logic from lamp_coverage.py: tries different linker lengths (0-9 nt)
    and different error levels to find the best split.
    
    Paramètres / Parameters:
        primers_dict  : dictionnaire {nom: sequence} des amorces / dict {name: seq}
        ref_ungapped_str : séquence de référence sans gaps (string)
        max_errors    : nombre max d'erreurs tolérées / max allowed errors
    
    Retourne un nouveau dictionnaire d'amorces avec FIP/BIP remplacés par leurs sous-parties.
    Returns a new primer dict with FIP/BIP replaced by their sub-parts.
    """
    result = dict(primers_dict)

    for combo_name in ['FIP', 'BIP']:
        if combo_name not in result:
            continue

        seq_combo = result[combo_name]
        found_split = False

        # Essaie 0, 1, puis 2 erreurs. / Try 0, then 1, then 2 errors.
        for allowed_err in [0, 1, 2]:
            if found_split:
                break
            # Essaie plusieurs longueurs de linker (0 = pas de linker, puis 4, 1, 2…9 nt)
            # Try different linker lengths (0 = no linker, then 4, 1, 2…9 nt)
            for linker_len in [0, 4, 1, 2, 3, 5, 6, 7, 8, 9]:
                if found_split:
                    break

                L = len(seq_combo)
                valid_splits = []

                for i in range(14, L - 13):
                    if i + linker_len > L - 14:
                        continue

                    # Contrainte biologique : le linker doit être une base répétée (ex : TTTT)
                    # Biological constraint: linker must be a repeated single base (e.g. TTTT)
                    if linker_len > 0:
                        linker_seq = seq_combo[i:i + linker_len]
                        if len(set(linker_seq)) != 1:
                            continue

                    part1 = seq_combo[:i]
                    part2 = seq_combo[i + linker_len:]

                    pos1 = primer_matches_sequence_simple(ref_ungapped_str, part1, max_errors=allowed_err)
                    pos2 = primer_matches_sequence_simple(ref_ungapped_str, part2, max_errors=allowed_err)

                    if pos1 and pos2:
                        dist = abs(pos2[0] - pos1[0])
                        if dist < 250:
                            # FIP : F1c (part1) est en aval de F2 (part2) → pos2 < pos1
                            # BIP : B1c (part1) est en amont de B2 (part2) → pos1 < pos2
                            if combo_name == 'FIP' and pos2[0] < pos1[0]:
                                valid_splits.append((part1, part2, linker_len))
                            elif combo_name == 'BIP' and pos1[0] < pos2[0]:
                                valid_splits.append((part1, part2, linker_len))

                if valid_splits:
                    # En cas d'ambiguïté, on choisit la coupure la plus équilibrée
                    # If ambiguous, pick the most balanced split
                    best = min(valid_splits, key=lambda x: abs(len(x[0]) - len(x[1])))
                    part1, part2, lk = best

                    if combo_name == 'FIP':
                        # F1c = reverse complement de part1, F2 = part2
                        result['F1c'] = str(Seq(part1).reverse_complement())
                        result['F2']  = part2
                        print(f"  [FIP split] F1c={result['F1c']} | F2={result['F2']} (linker {lk} nt)")
                    else:
                        # B1c = reverse complement de part1, B2 = part2
                        result['B1c'] = str(Seq(part1).reverse_complement())
                        result['B2']  = part2
                        print(f"  [BIP split] B1c={result['B1c']} | B2={result['B2']} (linker {lk} nt)")

                    del result[combo_name]
                    found_split = True
                    break  # Linker trouvé / linker found

        # Fallback heuristique : cherche TTTT dans la partie centrale de la séquence
        # Heuristic fallback: look for TTTT in the middle part of the sequence
        if not found_split:
            m = re.search(r'T{4,}', seq_combo[14:-14])
            if m:
                linker_start = m.start() + 14
                part1 = seq_combo[:linker_start]
                part2 = seq_combo[linker_start + 4:]

                if combo_name == 'FIP':
                    result['F1c'] = str(Seq(part1).reverse_complement())
                    result['F2']  = part2
                    print(f"  [FIP split heuristique] F1c={result['F1c']} | F2={result['F2']}")
                else:
                    result['B1c'] = str(Seq(part1).reverse_complement())
                    result['B2']  = part2
                    print(f"  [BIP split heuristique] B1c={result['B1c']} | B2={result['B2']}")

                del result[combo_name]
                found_split = True

        if not found_split:
            print(f"  [!] Impossible de déterminer la coupure automatique pour {combo_name}. L'amorce sera alignée telle quelle.")

    return result

def load_primers(filepath):
    """
    Charge les amorces depuis un fichier (FASTA ou Texte simple).
    Supports FASTA files and plain text files with 'Name Sequence' format.
    """
    primers = []
    try:
        with open(filepath, 'r') as f:
            first_line = f.readline()
            if not first_line:
                return []

            # Détection de fichier RTF (TextEdit Mac) / RTF file detection (TextEdit Mac)
            if first_line.startswith('{\\rtf'):
                print(f"Erreur : Le fichier d'amorces '{filepath}' est au format RTF.")
                print("Veuillez convertir ce fichier en Texte Brut (Format > Convertir au format Texte Brut).")
                sys.exit(1)

            f.seek(0)
            first_char = f.read(1)
            f.seek(0)

            if first_char == '>':
                # Fichier FASTA / FASTA file
                primers = list(SeqIO.parse(f, "fasta"))
            else:
                # Fichier texte simple / Plain text file
                primer_count = 1
                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    parts = line.split()
                    if len(parts) >= 2:
                        name = parts[0]
                        seq = "".join(parts[1:])
                    else:
                        name = f"Primer_{primer_count}"
                        seq = parts[0]
                        primer_count += 1

                    # Nettoyage des caractères parasites / Clean up stray characters
                    seq = "".join(c for c in seq if c.isalpha())
                    primers.append(SeqRecord(Seq(seq), id=name, description=""))
    except Exception as e:
        print(f"Erreur lors de la lecture du fichier d'amorces : {e}")
        sys.exit(1)

    return primers

def inject_gaps(primer_seq, gapped_ref_sub):
    """
    Injecte les tirets de la séquence de référence dans la séquence de l'amorce
    pour conserver l'alignement en colonnes.
    Injects gaps from the reference into the primer sequence to preserve column alignment.
    """
    result = ""
    p_idx = 0
    for char in gapped_ref_sub:
        if char == '-':
            result += '-'
        else:
            if p_idx < len(primer_seq):
                result += primer_seq[p_idx]
                p_idx += 1
            else:
                result += '-'
    return result

def get_expected_strand(primer_id):
    pid = primer_id.upper()
    for tag in ['_RP', '_R1', '_R2', '_R3', '_BIP', '_B3', '_R_']:
        if tag in pid or pid.endswith('RP') or pid.endswith('R1') or pid.endswith('R2') or pid.endswith('_R'):
            return 'REV'
    for tag in ['_FP', '_F1', '_F2', '_F3', '_FIP', '_F_']:
        if tag in pid or pid.endswith('FP') or pid.endswith('F1') or pid.endswith('F2') or pid.endswith('_F'):
            return 'FWD'
    return 'ANY'

def align_one_primer(primer_id, primer_seq, ref_ungapped_str, ref_gapped_str,
                     ungapped_to_gapped, msa_len, max_errors, silent=False):
    """
    Aligne une seule amorce sur la séquence de référence.
    Aligns a single primer against the reference sequence.
    Retourne / Returns : SeqRecord prêt à être écrit, ou None si non trouvé.
    """
    expected = get_expected_strand(primer_id)

    # Recherche brin sens / Forward strand search
    match_fwd = find_best_match(ref_ungapped_str, primer_seq, max_errors) if expected in ['FWD', 'ANY'] else None

    # Recherche brin antisens (reverse complement) / Reverse strand search
    primer_rc_str = str(Seq(primer_seq).reverse_complement())
    match_rev = find_best_match(ref_ungapped_str, primer_rc_str, max_errors) if expected in ['REV', 'ANY'] else None

    best_is_fwd = True
    best_match = None

    if match_fwd and match_rev:
        best_is_fwd = match_fwd[3] <= match_rev[3]
        best_match = match_fwd if best_is_fwd else match_rev
    elif match_fwd:
        best_match = match_fwd
        best_is_fwd = True
    elif match_rev:
        best_match = match_rev
        best_is_fwd = False

    if not best_match:
        return None

    start_ungapped, end_ungapped, _, errors = best_match

    # Conversion en coordonnées gappées / Convert to gapped coordinates
    start_gapped = ungapped_to_gapped[start_ungapped]
    end_gapped   = ungapped_to_gapped[end_ungapped - 1] + 1

    # Région de référence avec gaps pour injecter les tirets dans l'amorce
    # Gapped reference region to inject gaps into the primer
    gapped_ref_sub = ref_gapped_str[start_gapped:end_gapped]

    actual_seq = primer_seq if best_is_fwd else primer_rc_str
    aligned_sub = inject_gaps(actual_seq, gapped_ref_sub)

    padded_seq = "-" * start_gapped + aligned_sub + "-" * (msa_len - end_gapped)

    out_id  = primer_id if best_is_fwd else primer_id + "+c"
    strand  = "Forward" if best_is_fwd else "Reverse"
    desc    = f"Errors={errors} Pos={start_gapped+1}-{end_gapped} Strand={strand}"

    if not silent:
        tqdm.write(f"  [+] Trouvé : {out_id} (Erreurs: {errors}, Position MSA: {start_gapped+1}-{end_gapped})")
    return SeqRecord(Seq(padded_seq), id=out_id, description=desc)

global_targets_mappings = []

def init_worker(mappings):
    global global_targets_mappings
    global_targets_mappings = mappings

def process_primer_task(args):
    """
    Fonction exécutée par chaque processus pour aligner une amorce.
    """
    primer_id, primer_seq, ref_ungapped_str, ref_gapped_str, ungapped_to_gapped, msa_len, max_errors = args
    
    record = align_one_primer(
        primer_id, primer_seq,
        ref_ungapped_str, ref_gapped_str,
        ungapped_to_gapped, msa_len,
        max_errors, silent=True
    )
    
    if record:
        pos_info = [p for p in record.description.split() if p.startswith("Pos=")][0].split("=")[1]
        err_info = [p for p in record.description.split() if p.startswith("Errors=")][0].split("=")[1]
        msg = f"  [+] Trouvé : {record.id} (Erreurs: {err_info}, Position MSA: {pos_info})"
        return record, msg

    # Fallback : recherche sur toutes les autres séquences cibles du MSA
    for alt_id, alt_ungapped, alt_gapped, alt_map in global_targets_mappings:
        record = align_one_primer(
            primer_id, primer_seq,
            alt_ungapped, alt_gapped,
            alt_map, msa_len,
            max_errors,
            silent=True
        )
        if record:
            record.description += f" AltRef={alt_id}"
            pos_info = [p for p in record.description.split() if p.startswith("Pos=")][0].split("=")[1]
            err_info = [p for p in record.description.split() if p.startswith("Errors=")][0].split("=")[1]
            msg = f"  [+] Trouvé (via alt ref: {alt_id[:30]}) : {record.id} (Erreurs: {err_info}, Position MSA: {pos_info})"
            return record, msg
            
    return None, f"  [-] Non trouvé dans toute la base / Not found in entire database : {primer_id}"

def main():
    parser = argparse.ArgumentParser(
        description="Aligne des amorces (dont FIP/BIP) sur un MSA ou un génome simple. / "
                    "Aligns primers (including FIP/BIP) onto a MSA or a single genome.")
    parser.add_argument("-t", "--target",  required=True,
                        help="Fichier FASTA cible (simple ou MSA aligné) / Target FASTA (single or aligned MSA)")
    parser.add_argument("-p", "--primers", required=True,
                        help="Fichier FASTA ou TXT des amorces / Primers FASTA or TXT file")
    parser.add_argument("-o", "--output",  required=True,
                        help="Fichier FASTA de sortie / Output FASTA file")
    parser.add_argument("-e", "--errors",  type=int, default=2,
                        help="Nombre max d'erreurs tolérées / Max errors tolerated. Défaut/Default: 2")

    args = parser.parse_args()

    # Chargement de la séquence cible avec barre de progression / Load target with progress bar
    try:
        targets = list(tqdm(
            SeqIO.parse(args.target, "fasta"),
            desc="📂 Chargement séquences / Loading sequences",
            unit=" seq",
            colour="cyan"
        ))
    except Exception as e:
        print(f"Erreur lors de la lecture du fichier cible : {e}")
        sys.exit(1)

    if not targets:
        print("Erreur: Le fichier cible est vide. / Error: target file is empty.")
        sys.exit(1)

    # Chargement des amorces / Load primers
    primers_list = load_primers(args.primers)

    # Construction du dictionnaire {nom: sequence} / Build {name: sequence} dict
    primers_dict = {p.id: str(p.seq).upper() for p in primers_list}

    # Recherche de la meilleure séquence de référence (celle qui matche le plus de primers)
    # Find the best reference sequence (the one matching the highest number of primers)
    best_ref_idx = 0
    best_match_count = -1
    
    print("\n🔍 Recherche de la meilleure séquence de référence pour l'alignement...")
    for idx, rec in enumerate(targets):
        ungapped_str, _ = get_ungapped_mapping(str(rec.seq).upper())
        matches = 0
        for p_seq in primers_dict.values():
            if primer_matches_sequence_simple(ungapped_str, p_seq, max_errors=args.errors):
                matches += 1
        
        if matches > best_match_count:
            best_match_count = matches
            best_ref_idx = idx
            # Arrêt rapide si on trouve une séquence qui matche toutes les amorces
            # Short-circuit if we find a sequence matching all primers
            if matches == len(primers_dict):
                break
                
    ref_record      = targets[best_ref_idx]
    ref_gapped_str  = str(ref_record.seq).upper()
    msa_len         = len(ref_gapped_str)
    print(f"Séquence de référence sélectionnée : {ref_record.id} (Matches: {best_match_count}/{len(primers_dict)}, Longueur MSA: {msa_len})")
    
    # Mappage ungapped -> gapped pour la référence / ungapped -> gapped mapping
    ref_ungapped_str, ungapped_to_gapped = get_ungapped_mapping(ref_gapped_str)

    # ─────────────────────────────────────────────────────
    # Séparation automatique des amorces FIP et BIP
    # Automatic splitting of FIP and BIP primers
    # ─────────────────────────────────────────────────────
    # Recherche les amorces dont le nom se termine par FIP ou BIP (ex: '1_FIP', 'Set2_BIP', 'FIP')
    # Find primers whose name ends with FIP or BIP (e.g. '1_FIP', 'Set2_BIP', 'FIP')
    fip_bip_keys = [k for k in primers_dict if k.upper().endswith('FIP') or k.upper().endswith('BIP')]
    if fip_bip_keys:
        print("\nDétection d'amorces FIP/BIP – tentative de séparation automatique...")
        print("FIP/BIP primers detected – attempting automatic split...\n")

        new_primers_dict = {}
        for key, seq in primers_dict.items():
            key_upper = key.upper()
            if key_upper.endswith('FIP'):
                combo_type = 'FIP'
            elif key_upper.endswith('BIP'):
                combo_type = 'BIP'
            else:
                new_primers_dict[key] = seq
                continue

            # Détermine le préfixe (ex: '1_' pour '1_FIP')
            # Determine the prefix (e.g. '1_' for '1_FIP')
            prefix = key[: len(key) - len(combo_type)]

            # Crée un dict temporaire avec juste cette amorce pour la passer à auto_split_fip_bip
            temp = {combo_type: seq}
            split_result = auto_split_fip_bip(temp, ref_ungapped_str, max_errors=args.errors)

            if combo_type not in split_result:  # Split réussi / split succeeded
                for part_name, part_seq in split_result.items():
                    # Réinjecte le préfixe original / Re-inject original prefix
                    new_primers_dict[prefix + part_name] = part_seq
            else:  # Pas de split, on garde l'amorce telle quelle / No split, keep as-is
                new_primers_dict[key] = seq

        primers_dict = new_primers_dict
        print()

    # ─────────────────────────────────────────────────────
    # Pré-calcul des mappages des autres séquences cibles pour le fallback
    # Precompute mappings of other target sequences for fallback search
    # ─────────────────────────────────────────────────────
    targets_mappings = []
    for rec in targets:
        if rec.id == ref_record.id:
            continue
        rec_gapped = str(rec.seq).upper()
        rec_ungapped, rec_map = get_ungapped_mapping(rec_gapped)
        targets_mappings.append((rec.id, rec_ungapped, rec_gapped, rec_map))

    # ─────────────────────────────────────────────────────
    # Alignement de chaque amorce avec ProcessPoolExecutor
    # Align each primer using ProcessPoolExecutor
    # ─────────────────────────────────────────────────────
    out_records = list(targets)  # Toutes les séquences cibles en premier / All target seqs first

    workers = min(8, max(1, (os.cpu_count() or 4) - 2))
    print(f"\n🚀 Utilisation de {workers} processus pour l'alignement / Using {workers} processes for alignment")

    tasks = []
    for primer_id, primer_seq in primers_dict.items():
        tasks.append((primer_id, primer_seq, ref_ungapped_str, ref_gapped_str, ungapped_to_gapped, msa_len, args.errors))

    with concurrent.futures.ProcessPoolExecutor(max_workers=workers, initializer=init_worker, initargs=(targets_mappings,)) as executor:
        futures = {executor.submit(process_primer_task, task): task[0] for task in tasks}
        for future in tqdm(
            concurrent.futures.as_completed(futures),
            total=len(futures),
            desc="🧬 Alignement amorces / Aligning primers",
            unit=" primer",
            colour="magenta"
        ):
            record, msg = future.result()
            if record:
                out_records.append(record)
            tqdm.write(msg)

    # Écriture du fichier de sortie / Write output file
    try:
        SeqIO.write(out_records, args.output, "fasta")
        print(f"\nAlignement terminé. Résultat sauvegardé dans {args.output}")
        print(f"Alignment done. Result saved to {args.output}")
    except Exception as e:
        print(f"Erreur lors de l'écriture du fichier : {e}")

if __name__ == "__main__":
    main()
