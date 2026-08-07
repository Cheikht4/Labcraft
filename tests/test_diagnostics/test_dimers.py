import pytest
from labcraft.diagnostics.enzyme import BST_2_0
from labcraft.diagnostics.amplifiable_dimer import is_amplifiable_dimer
import RNA

def test_amplifiable_dimer_cases():
    """Valide la détection d'amplifiabilité selon les règles métier (3' apparié + matrice)."""
    
    # Cas 1 : Dimère stable avec 3' apparié et matrice (Extensible)
    # A: 5' ATGCATGCATGC 3'
    # B: 3' TACGTACGTACG 5'
    # B's 3' end doesn't matter, A's 3' binds B's 5' region (with matrix on B)
    # Let A = "ATGCATGCATGC", B = "GCATGCATGCAT"
    # A (5'->3') ATGCATGCATGC
    # B (3'->5') TACGTACGTACG
    # They are perfectly complementary! 
    # The 3' of A is C. It pairs with G at the 5' end of B. 
    # Wait, if A pairs exactly at the 5' end of B, there is NO matrix left to copy!
    # Let's verify: A pairs with B.
    
    # Let's make a case WHERE there IS matrix.
    # A = "ATGCATGCGCG" (11 nt)
    # B = "TTTTCGCGCATGCAT" (15 nt)
    # A (5'->3'):          A T G C A T G C G C G
    # B (3'->5'):  T A C G T A C G C G C T T T T
    # B in 5'->3' is TTTTCGCGCATGCAT.
    
    p_a = "ATGCATGCGCG"
    p_b = "TTTTCGCGCATGCAT"
    struct, mfe = RNA.cofold(f"{p_a}&{p_b}")
    
    is_amp, dg_3p = is_amplifiable_dimer(p_a, p_b, struct, mfe, BST_2_0)
    assert is_amp, "Le dimère a un 3' apparié et une matrice 5' sur le partenaire, il DOIT être amplifiable."
    
    # Cas 2 : Dimère avec 3' apparié mais SANS matrice (Blunt end / non extensible)
    # p_a = "ATGC"
    # p_b = "GCAT"
    p_a = "ATGCATGCGCG"
    p_b = "CGCGCATGCAT"
    struct, mfe = RNA.cofold(f"{p_a}&{p_b}")
    is_amp, dg_3p = is_amplifiable_dimer(p_a, p_b, struct, mfe, BST_2_0)
    assert not is_amp, "Le 3' est apparié au bord 5' exact, pas de matrice au-delà, NE DOIT PAS être amplifiable."
    
    # Cas 3 : Dimère instable au 3' ou non apparié au 3'
    # p_a = "ATGCATGCATGC" + "TTTT" (3' mismatch)
    # p_b = "GCATGCATGCATTTT"
    p_a = "ATGCATGCGCGTTTT"
    p_b = "CGCGCATGCATTTT"
    struct, mfe = RNA.cofold(f"{p_a}&{p_b}")
    is_amp, dg_3p = is_amplifiable_dimer(p_a, p_b, struct, mfe, BST_2_0)
    assert not is_amp, "Le 3' de A est composé de TTTT non appariés, NE DOIT PAS être amplifiable."

