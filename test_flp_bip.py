import RNA
bip = "CTGTAGCTCCGTCGTGGGGATTTTCTAGTCTGCTACACCGTGC"
flp = "CCTTGGACGGGGCT"
seq = f"{flp}&{bip}"
RNA.cvar.temperature = 63.0
struct, mfe = RNA.cofold(seq)
print(f"FLP+BIP structure:\n{flp}\n{struct}\n{bip}\nmfe: {mfe}")

seq2 = f"{bip}&{flp}"
struct2, mfe2 = RNA.cofold(seq2)
print(f"BIP+FLP structure:\n{bip}\n{struct2}\n{flp}\nmfe: {mfe2}")
