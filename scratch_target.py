import RNA
from labcraft.target.unfolding import calc_unfolding_penalty
import random

random.seed(123)
def random_seq(n):
    return "".join(random.choices(["A", "C", "G", "T"], k=n))

def test_target():
    for i in range(10000):
        random_1 = random_seq(20)
        random_2 = random_seq(20)
        
        f3_a = "GCTAGCAATCGTACGCATAG"
        stem_left = "GGCATGCCTA"
        revcomp_stem = "TAGGCATGCC"
        
        seq = random_1 + stem_left + f3_a + revcomp_stem + random_2
        
        dg_c1 = calc_unfolding_penalty(seq, 0, 20, temp_celsius=65.0)
        dg_c2 = calc_unfolding_penalty(seq, 60, 80, temp_celsius=65.0)
        
        if dg_c1 < 1.0 and dg_c2 < 1.0:
            dg_f3 = calc_unfolding_penalty(seq, 30, 50, temp_celsius=65.0)
            if dg_f3 > 6.0:
                print(f">SynthA\n{seq}")
                print(f"F3_A unfold: {dg_f3}")
                print(f"Control 1 unfold: {dg_c1}")
                print(f"Control 2 unfold: {dg_c2}")
                break

test_target()
