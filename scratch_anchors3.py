import RNA
RNA.cvar.temperature = 65.0
RNA.params_load_DNA_Mathews2004()

seq_dengue4 = "GCGCGCATCG" + "AAAA" + "CGATGCGCGC" 
fc3 = RNA.fold_compound(seq_dengue4)
_, g_free3 = fc3.pf()
fc3_cons = RNA.fold_compound(seq_dengue4)
for i in range(1, 11):
    fc3_cons.hc_add_up(i)
_, g_cons3 = fc3_cons.pf()
print(f"Stem GCGCGCATCG dG={g_cons3-g_free3:.3f}")
