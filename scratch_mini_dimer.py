import RNA
from labcraft.diagnostics.enzyme import BST_2_0
from labcraft.diagnostics.amplifiable_dimer import is_amplifiable_dimer, get_pairs

p_a = "GATCGATCGGCCTTAAAAA" # F3
p_b = "TTTTTTTTAAGGCC"      # B3
# F3: 19 nt. B3: 14 nt.
# We expect B3's 3' end (AAGGCC) to pair with F3's GGCCTT.
# GGCCTT is at index 8..13 in F3.
# Let's check RNA.cofold
struct, mfe = RNA.cofold(f"{p_a}&{p_b}")
print(f"F3: {p_a}")
print(f"B3: {p_b}")
print(f"Structure: {struct}  MFE: {mfe}")
is_amp, dg = is_amplifiable_dimer(p_a, p_b, struct, mfe, BST_2_0, 65.0)
print(f"Amplifiable: {is_amp}  dG_3p: {dg}")

# Let's also check B3 as first arg and F3 as second
struct2, mfe2 = RNA.cofold(f"{p_b}&{p_a}")
is_amp2, dg2 = is_amplifiable_dimer(p_b, p_a, struct2, mfe2, BST_2_0, 65.0)
print(f"Amplifiable (rev): {is_amp2}  dG_3p: {dg2}")
