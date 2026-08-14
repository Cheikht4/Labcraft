from Bio import SeqIO
from collections import defaultdict
from lamp_coverage import load_primers, auto_split_fip_bip

targets = {}
for record in SeqIO.parse("sequences_DENGUE_1_NCBI_03_04_2026.fasta", "fasta"):
    clean_seq = str(record.seq).upper().replace('-', '')
    if len(clean_seq) > 100:
        targets[record.description] = clean_seq
        if len(targets) >= 5:
            break

primer_sets = load_primers('primer_dengue_1_FIPBIP.txt')
txt = {
    'split_success': "Auto-split réussi pour Set {} ({}) -> {} et {} (linker de {} nt ignoré).",
    'split_fail': "Attention : Impossible de déterminer la coupe automatique pour Set {} ({})."
}
auto_split_fip_bip(primer_sets, targets, txt)

print("Comparing F1:")
print("1 F1 :", primer_sets['1']['F1'])
print("1S F1:", primer_sets['1S']['F1'])
print("Are they equal?", primer_sets['1']['F1'] == primer_sets['1S']['F1'])

print("\nComparing B2:")
print("1 B2 :", primer_sets['1']['B2'])
print("1S B2:", primer_sets['1S']['B2'])
print("Are they equal?", primer_sets['1']['B2'] == primer_sets['1S']['B2'])

