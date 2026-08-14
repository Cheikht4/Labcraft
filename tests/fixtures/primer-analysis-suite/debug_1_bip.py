from Bio import SeqIO
from lamp_coverage import primer_matches_sequence

bip = "CCCAACACCAGGGGAAGCTGTTTTTTTGTTGTTGTGCGGGGG"
targets = []
for record in SeqIO.parse("sequences_DENGUE_1_NCBI_03_04_2026.fasta", "fasta"):
    t = str(record.seq).upper().replace('-', '')
    if len(t) > 100:
        targets.append(t)
        if len(targets) >= 50:
            break

found = False
for count, t in enumerate(targets):
    if found: break
    for allowed_err in [0, 1, 2, 3, 4]:
        if found: break
        for linker_len in [4, 7]:
            if found: break
            for i in range(14, len(bip)-13):
                if i + linker_len > len(bip) - 14:
                    continue
                p1 = bip[:i]
                p2 = bip[i+linker_len:]
                pos1 = primer_matches_sequence(t, p1, max_errors=allowed_err, strict_3prime_len=0)
                pos2 = primer_matches_sequence(t, p2, max_errors=allowed_err, strict_3prime_len=0)
                if pos1 and pos2:
                    print(f"Match on target {count} | err={allowed_err} | linker={linker_len} | i={i}")
                    print(f"p1: {p1}")
                    print(f"p2: {p2}")
                    found = True
                    break

if not found:
    print("NO MATCH AT ALL on first 50 targets, even with 4 errors.")
