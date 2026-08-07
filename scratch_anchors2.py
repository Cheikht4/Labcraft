import RNA
RNA.cvar.temperature = 65.0
RNA.params_load_DNA_Mathews2004()

# Let's find a stem around 8.5
seq_dengue = "CGCGATAT" + "AAAA" + "ATATCGCG" # 8 bp stem
fc3 = RNA.fold_compound(seq_dengue)
s3, g_free3 = fc3.pf()
fc3_cons = RNA.fold_compound(seq_dengue)
for i in range(1, 9):
    fc3_cons.hc_add_up(i)
_, g_cons3 = fc3_cons.pf()
print(f"Stem 8bp dG={g_cons3-g_free3:.3f}")

seq_dengue2 = "GCGCATATGC" + "AAAA" + "GCATATGCGC" # 10 bp
fc3 = RNA.fold_compound(seq_dengue2)
s3, g_free3 = fc3.pf()
fc3_cons = RNA.fold_compound(seq_dengue2)
for i in range(1, 11):
    fc3_cons.hc_add_up(i)
_, g_cons3 = fc3_cons.pf()
print(f"Stem 10bp dG={g_cons3-g_free3:.3f}")

seq_dengue3 = "CGCGCGCGCG" + "AAAA" + "CGCGCGCGCG" # 10 bp all GC
fc3 = RNA.fold_compound(seq_dengue3)
s3, g_free3 = fc3.pf()
fc3_cons = RNA.fold_compound(seq_dengue3)
for i in range(1, 11):
    fc3_cons.hc_add_up(i)
_, g_cons3 = fc3_cons.pf()
print(f"Stem 10bp GC dG={g_cons3-g_free3:.3f}")
