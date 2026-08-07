import RNA
from labcraft.diagnostics.enzyme import BST_2_0
from labcraft.diagnostics.amplifiable_dimer import is_amplifiable_dimer, get_pairs
p_a = "ATGCATGCATGC"
p_b = "TTTTGCATGCATGCAT"
struct, mfe = RNA.cofold(f"{p_a}&{p_b}")
print("struct:", struct, "mfe:", mfe)
pairs = get_pairs(struct)
print("pairs:", pairs)
l_a = len(p_a)
print("L_A:", l_a)
print("3' of A is index", l_a - 1)
if l_a - 1 in pairs:
    print("A 3' pairs with", pairs[l_a - 1])
else:
    print("A 3' is NOT PAIRED!")
