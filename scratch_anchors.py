import RNA

# Setting temperature to 65.0
RNA.cvar.temperature = 65.0
RNA.params_load_DNA_Mathews2004()

# 1. Negative control
seq_unstructured = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
fc1 = RNA.fold_compound(seq_unstructured)
_, g_free1 = fc1.pf()
fc1.hc_add_up(10)
_, g_cons1 = fc1.pf()
print(f"Negative control: G_free={g_free1:.3f}, G_cons={g_cons1:.3f}, dG={g_cons1-g_free1:.3f}")

# 2. Strong hairpin match
seq_hp = "GGGGGGGGGGGGCCAAAAAAGGCCCCCCCCCCCC"
fc2 = RNA.fold_compound(seq_hp)
s, g_free2 = fc2.pf()
print(f"Hairpin 65C: {s} {g_free2:.3f}")
# site binding one arm
fc2_cons = RNA.fold_compound(seq_hp)
for i in range(1, 13): # constrain GGGGGGGGGGGG
    fc2_cons.hc_add_up(i)
_, g_cons2 = fc2_cons.pf()
print(f"Hairpin constrained: G_cons={g_cons2:.3f}, dG={g_cons2-g_free2:.3f}")

# 3. Dengue equivalent (synthetic) giving ~8.5 kcal/mol
# Let's try some stem
seq_dengue = "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG" # just a very strong GC stem
seq_dengue = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC" + "AAAA" + "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"
fc3 = RNA.fold_compound(seq_dengue)
s3, g_free3 = fc3.pf()
print(f"Dengue stem 65C: {s3} {g_free3:.3f}")
fc3_cons = RNA.fold_compound(seq_dengue)
for i in range(1, 17):
    fc3_cons.hc_add_up(i)
_, g_cons3 = fc3_cons.pf()
print(f"Dengue dG={g_cons3-g_free3:.3f}")
