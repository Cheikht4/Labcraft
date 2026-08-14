from Bio import SeqIO
from lamp_coverage import load_primers, primer_matches_sequence

targets = []
for record in SeqIO.parse("sequences_DENGUE_1_NCBI_03_04_2026.fasta", "fasta"):
    clean_seq = str(record.seq).upper().replace('-', '')
    if len(clean_seq) > 100:
        targets.append(clean_seq)
        if len(targets) >= 5:
            break

primer_sets = load_primers('Nawar_dengue1.txt')

fip = primer_sets['1']['FIP']
print("Testing 1_FIP:", fip)
# Test manual split based on TTTT
part1 = "GCTGCGTTGTGTCTTGGGAGG"
part2 = "CTGTACGCATGGGGTAGC"
print("part1 (len %d): %s" % (len(part1), part1))
print("part2 (len %d): %s" % (len(part2), part2))
for i, t in enumerate(targets):
    p1 = primer_matches_sequence(t, part1, max_errors=2, strict_3prime_len=0)
    p2 = primer_matches_sequence(t, part2, max_errors=2, strict_3prime_len=0)
    print("Target %d -> pos1: %s, pos2: %s" % (i, p1, p2))

d_fip = primer_sets['DENV1']['FIP']
print("\nTesting DENV1_FIP:", d_fip)
found = False
for allowed_err in [0, 1, 2, 3, 4]:
    if found: break
    for linker_len in [0, 1, 2, 3, 4, 5, 6]:
        if found: break
        for i in range(14, len(d_fip)-13):
            p1 = d_fip[:i]
            p2 = d_fip[i+linker_len:]
            pos1 = primer_matches_sequence(targets[0], p1, max_errors=allowed_err, strict_3prime_len=0)
            pos2 = primer_matches_sequence(targets[0], p2, max_errors=allowed_err, strict_3prime_len=0)
            if pos1 and pos2:
                print("Found match with err=%d, linker=%d, i=%d" % (allowed_err, linker_len, i))
                print("p1:", p1)
                print("p2:", p2)
                found = True
                break

if not found:
    print("No split found for DENV1_FIP on target 0 even with 4 errors!")
