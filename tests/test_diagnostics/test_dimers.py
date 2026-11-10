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
    
    is_amp, dg_3p, ext, blocked_veto = is_amplifiable_dimer(p_a, p_b, struct, mfe, BST_2_0)
    assert is_amp, "Le dimère a un 3' apparié et une matrice 5' sur le partenaire, il DOIT être amplifiable."
    
    # Cas 2 : Dimère avec 3' apparié mais SANS matrice (Blunt end / non extensible)
    # p_a = "ATGC"
    # p_b = "GCAT"
    p_a = "ATGCATGCGCG"
    p_b = "CGCGCATGCAT"
    struct, mfe = RNA.cofold(f"{p_a}&{p_b}")
    is_amp, dg_3p, ext, blocked_veto = is_amplifiable_dimer(p_a, p_b, struct, mfe, BST_2_0)
    assert not is_amp, "Le 3' est apparié au bord 5' exact, pas de matrice au-delà, NE DOIT PAS être amplifiable."
    
    # Cas 3 : Dimère instable au 3' ou non apparié au 3'
    # p_a = "ATGCATGCATGC" + "TTTT" (3' mismatch)
    # p_b = "GCATGCATGCATTTT"
    p_a = "ATGCATGCGCGTTTT"
    p_b = "CGCGCATGCATTTT"
    struct, mfe = RNA.cofold(f"{p_a}&{p_b}")
    is_amp, dg_3p, ext, blocked_veto = is_amplifiable_dimer(p_a, p_b, struct, mfe, BST_2_0)
    assert not is_amp, "Le 3' de A est composé de TTTT non appariés, NE DOIT PAS être amplifiable."


def test_part_c_local_dg_3p_salt_correction_and_monotonic_mg():
    """Partie C : Vérifie que :
    1. Le ΔG local en 3' d'un même segment diffère entre 1 M Na+ et 50 mM Na + 8 mM Mg + 1,4 mM dNTP
       dans le sens physique donné par le modèle Owczarzy 2004/2008.
    2. Le ΔG local varie de façon monotone (devient plus négatif / plus stable) lorsque le Mg2+ libre augmente.
    """
    from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
    from labcraft.diagnostics.amplifiable_dimer import evaluate_pair_amplifiable
    from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
    from labcraft.thermo.salt import UnifiedSaltModel

    # Amorce LB_3#1 du panel 3 qui forme un homodimère avec segment 3' de 3 paires (GGA / CCT)
    # à 1 M Na+ : local_dg = -2.05 kcal/mol
    p_lb1 = PhysicalPrimer("LB_3#1", "AAGGTCATCCAACTTAGTAGGA", PrimerRole.LB, "AAGGTCATCCAACTTAGTAGGA")
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())

    # 1. Comparaison 1 M Na+ (1000 mM Na, 0 mM Mg) et 50 mM Na + 8 mM Mg + 1.4 mM dNTP à 63 °C
    is_amp_1m, dg_1m, _ = evaluate_pair_amplifiable(
        p_lb1, p_lb1, backend, BST_2_0, 63.0, na_mm=1000.0, mg_mm=0.0, dntp_mm=0.0
    )
    is_amp_lamp, dg_lamp, _ = evaluate_pair_amplifiable(
        p_lb1, p_lb1, backend, BST_2_0, 63.0, na_mm=50.0, mg_mm=8.0, dntp_mm=1.4
    )

    assert dg_1m == pytest.approx(-2.05, abs=0.02)
    assert is_amp_1m is True
    # Avec sodium_equivalent_for_folding (~0.2 M) + correction entropique SantaLucia 1998,
    # le ΔG local en 3' est moins stable (plus positif) qu'à 1 M Na+ (-1.66 > -2.05 kcal/mol)
    # et passe au-dessus du seuil -2.0 kcal/mol de Bst 2.0 (n'est plus amplifiable).
    # With sodium_equivalent_for_folding (~0.2 M) + SantaLucia 1998 entropic correction,
    # the local 3' dG is less stable (more positive) than at 1 M Na+ (-1.66 > -2.05 kcal/mol)
    # and crosses above the -2.0 kcal/mol threshold of Bst 2.0 (no longer amplifiable).
    assert dg_lamp > dg_1m, f"Attendu dg_lamp ({dg_lamp:.3f}) > dg_1m ({dg_1m:.3f})"
    assert dg_lamp > BST_2_0.dimer_dg_threshold
    assert is_amp_lamp is False

    # 2. Monotonie avec le Mg2+ libre (à 50 mM Na, 1.4 mM dNTP, Mg total croissant de 2 à 20 mM)
    mg_values = [2.0, 4.0, 6.0, 8.0, 12.0, 20.0]
    dgs = []
    for mg in mg_values:
        _, dg_val, _ = evaluate_pair_amplifiable(
            p_lb1, p_lb1, backend, BST_2_0, 63.0, na_mm=50.0, mg_mm=mg, dntp_mm=1.4
        )
        dgs.append(dg_val)

    for i in range(len(dgs) - 1):
        assert dgs[i + 1] < dgs[i], (
            f"Non-monotonie du ΔG local 3' avec le Mg libre : Mg={mg_values[i]} -> {dgs[i]:.4f}, "
            f"Mg={mg_values[i+1]} -> {dgs[i+1]:.4f}"
        )


def test_part_c_cli_panel3_dimer_reduction(tmp_path):
    """Partie C (CLI) : Vérifie via CliRunner sur `labcraft analyze` que la correction saline du ΔG local 3'
    est bien appliquée lors d'une exécution réelle sur le panel 3."""
    from typer.testing import CliRunner
    from pathlib import Path
    from labcraft.cli.main import app
    from labcraft.cli.parsers import read_multi_fasta

    fixture_dir = Path(__file__).resolve().parent.parent / "fixtures" / "real_primer_files"
    p3_file = fixture_dir / "primer_dengue_3.txt"
    d3_fasta = fixture_dir / "DEN3_M93130.fasta"

    p3_records = [
        f">{h}\n{s}" for h, s in read_multi_fasta(str(p3_file)) if h.startswith("3_")
    ]
    primers_only_3 = tmp_path / "panel_3.txt"
    primers_only_3.write_text("\n".join(p3_records) + "\n")
    out_html = tmp_path / "report_p3.html"

    runner = CliRunner()
    res = runner.invoke(app, [
        "analyze", "-p", str(primers_only_3), "-t", str(d3_fasta), "-o", str(out_html),
        "--temperature", "63.0", "--na", "50.0", "--mg", "8.0", "--dntp", "1.4"
    ])
    assert res.exit_code == 0, f"CLI failed: {res.stdout}"
    html = out_html.read_text()
    # LB_3#1_homo avait dg_3p = -2.05 à 1 M Na+, mais avec la correction saline (-1.72 > -2.0),
    # il ne doit plus apparaître à -2.05 kcal/mol dans les dimères amplifiables.
    assert "-2.05 kcal/mol" not in html


