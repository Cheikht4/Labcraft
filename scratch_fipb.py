import RNA
from labcraft.diagnostics.enzyme import BST_2_0
from labcraft.diagnostics.amplifiable_dimer import is_amplifiable_dimer, get_pairs
p_a = "ATCGATCGATCGATCGATCGTTTTCAGTACGCATAGCTAGCTAG" # FIP_A
p_b = "GATCGATCGATCGATCGATCTTTTCTAGCTGCATGCATGATCG" # BIP_B
struct, mfe = RNA.cofold(f"{p_a}&{p_b}")
print("struct:", struct, "mfe:", mfe)
pairs = get_pairs(struct)
l_a = len(p_a)
if l_a - 1 in pairs: print("A 3' pairs with", pairs[l_a - 1])
if l_a + len(p_b) - 1 in pairs: print("B 3' pairs with", pairs[l_a + len(p_b) - 1])
is_amp, dg = is_amplifiable_dimer(p_a, p_b, struct, mfe, BST_2_0)
print("Amplifiable:", is_amp, "dg:", dg)
