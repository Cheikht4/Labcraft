import RNA

print(RNA.__version__)
print(hasattr(RNA, 'params_load_DNA_Mathews2004'))

# Check how to do partition function with constraints
fc = RNA.fold_compound("ACGTACGTACGTACGT")
# Load DNA params
RNA.params_load_DNA_Mathews2004()
fc = RNA.fold_compound("ACGTACGTACGTACGT")
fc.pf()
print(fc.pf()) # Returns (ensemble_structure, free_energy)

# Try constraints
fc.hc_add_up(1, 4) # Forces bases 1 to 4 to be unpaired
res = fc.pf()
print("Constrained:", res)
