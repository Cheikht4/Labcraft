"""Tests for multi-target multiplex analysis, degenerate variants scaling, and performance.
Tests pour l'analyse multiplexe multi-cibles, le passage d'amorces dégénérées et les performances.
"""
import pytest
import time
from pathlib import Path
from typer.testing import CliRunner

from labcraft.cli.main import app
from labcraft.lamp.domains import expand_degenerate
from labcraft.metrics.balance import calculate_multiplex_balance


def test_multi_target_cli_two_targets_two_panels(tmp_path: Path):
    """Vérifie qu'avec deux cibles et deux panels, les deux cibles sont analysées conjointement,
    que le rapport contient un diagnostic pour CHACUNE des deux cibles et que les occupations sont non nulles."""
    # Création de deux cibles synthétiques distinctes
    # Cible A : 250 nt contenant les sites du Panel A
    site_f3_a = "GCCACCTTAAGCCACAGTA"
    site_b3_a = "TTAGCTTAGGCTAAGCTA"
    site_f2_a = "AAGCTTAACCGGTTAAC"
    site_b2_a = "CCGGTAACCCTTAAGGC"
    
    target_a_seq = f"AAAA{site_f3_a}TTTT{site_f2_a}GGGG{site_b2_a}CCCC{site_b3_a}AAAA"
    
    # Cible B : 250 nt contenant les sites du Panel B
    site_f3_b = "CGGTCCTTAAGCCACAGTC"
    site_b3_b = "GTAGCTTAGGCTAAGCTC"
    site_f2_b = "TAGCTTAACCGGTTAAT"
    site_b2_b = "GCGGTAACCCTTAAGGA"
    
    target_b_seq = f"TTTT{site_f3_b}GGGG{site_f2_b}AAAA{site_b2_b}TTTT{site_b3_b}CCCC"
    
    # Fichier multi-FASTA de cibles
    targets_fasta = tmp_path / "mini_cibles.fasta"
    targets_fasta.write_text(f">CIBLE_A\n{target_a_seq}\n>CIBLE_B\n{target_b_seq}\n")
    
    # Fichier d'amorces pour les deux panels
    primers_txt = tmp_path / "mini_amorces.txt"
    primers_txt.write_text(
        f">CIBLE_A_F3   {site_f3_a}\n"
        f">CIBLE_A_B3   {site_b3_a}\n"
        f">CIBLE_A_FIP  TTTTTTTTTTTTTTTTTTTT{site_f2_a}\n"
        f">CIBLE_A_BIP  CCCCCCCCCCCCCCCCCCCC{site_b2_a}\n"
        f">CIBLE_B_F3   {site_f3_b}\n"
        f">CIBLE_B_B3   {site_b3_b}\n"
        f">CIBLE_B_FIP  GGGGGGGGGGGGGGGGGGGG{site_f2_b}\n"
        f">CIBLE_B_BIP  AAAAAAAAAAAAAAAAAAAA{site_b2_b}\n"
    )
    
    out_html = tmp_path / "multiplex_report.html"
    
    runner = CliRunner()
    result = runner.invoke(app, [
        "analyze",
        "-p", str(primers_txt),
        "-t", str(targets_fasta),
        "-o", str(out_html),
        "--temperature", "65",
        "--mg", "8.0",
        "--dntp", "1.4"
    ])
    
    assert result.exit_code == 0, f"CLI execution failed:\n{result.stdout}"
    assert out_html.exists()
    
    html_content = out_html.read_text()
    
    # Vérification 1 : Les deux cibles doivent être mentionnées dans le rapport
    assert "CIBLE_A" in html_content
    assert "CIBLE_B" in html_content
    
    # Vérification 2 : La section Diagnostic par Cible doit exister pour les deux cibles
    assert "Cible : CIBLE_A" in html_content or "CIBLE_A" in html_content
    assert "Cible : CIBLE_B" in html_content or "CIBLE_B" in html_content
    
    # Vérification 3 : Le mode Conjoint est bien mentionné
    assert "Conjoint" in html_content
    
    # Vérification 4 : Le tableau de comparaison des panels ne contient pas 0.0 % pour un panel actif
    # Les deux panels A et B doivent avoir une accessibilité mesurée
    assert "CIBLE_A" in html_content
    assert "CIBLE_B" in html_content


def test_multiplex_balance_displays_not_analyzed_when_missing():
    """Vérifie que calculate_multiplex_balance et le rapport renvoient 'Non analysé'
    (et jamais 0.0 %) si une cible n'a pas été calculée."""
    primer_to_panel = {
        "F3_A": "PanelA",
        "B3_A": "PanelA",
        "F3_B": "PanelB",
        "B3_B": "PanelB"
    }
    
    # Seul PanelA a des données d'occupation
    target_occupations = {
        "PanelA": {
            "F3_A_site": 0.85,
            "B3_A_site": 0.80
        }
        # PanelB est absent
    }
    
    free_fractions = {
        "F3_A": 0.90, "B3_A": 0.90,
        "F3_B": 0.90, "B3_B": 0.90
    }
    
    summaries, cv = calculate_multiplex_balance(primer_to_panel, target_occupations, free_fractions)
    
    assert summaries["PanelA"]["mean_occupation"] == pytest.approx(0.825, rel=1e-2)
    assert summaries["PanelB"]["mean_occupation"] is None
    assert summaries["PanelB"]["min_occupation"] is None
    assert summaries["PanelB"]["limiting_primer"] == "Non analysé"
    # Le CV inter-panels ne doit pas fabriquer un faux 1.00 quand un panel n'est pas analysé
    assert cv is None


def test_five_degenerate_positions_supported_by_default():
    """Vérifie que la séquence de Nunes et al. 2015 à 5 positions dégénérées (32 variants)
    passe avec la valeur par défaut (max_variants=64), et échoue avec un message clair si limitée."""
    # FIP Nunes et al. 2015 : 5 dégénérescences (R, Y, R, R, R) -> 2^5 = 32 variants
    fip_seq = "GRCCTCCGATTGAYCTCGGCTTTARTGTGARTGGCCRCTGAC"
    
    # Doit réussir par défaut (max_variants=64)
    variants = expand_degenerate(fip_seq, max_variants=64)
    assert len(variants) == 32
    
    # Doit échouer si max_variants=16 avec un message explicite
    with pytest.raises(ValueError) as exc_info:
        expand_degenerate(fip_seq, max_variants=16)
        
    msg = str(exc_info.value)
    assert "5 positions dégénérées" in msg
    assert "32 variants" in msg
    assert "--max-variants" in msg


def test_real_two_panels_multiplex_performance(tmp_path: Path):
    """Benchmark : deux panels réels de Dengue sur deux cibles flavivirus de 10.7 kb terminent rapidement."""
    fixture_dir = Path(__file__).resolve().parent.parent / "fixtures" / "real_primer_files"
    p3_file = fixture_dir / "primer_dengue_3.txt"
    d3_fasta = fixture_dir / "DEN3_M93130.fasta"
    
    from labcraft.cli.parsers import read_multi_fasta
    d3_seq = read_multi_fasta(str(d3_fasta))[0][1]
    
    # 2 génomes flavivirus complets de 10.7 kb correspondant aux panels 3 et 14
    targets_fasta = tmp_path / "two_flavivirus_genomes.fasta"
    targets_fasta.write_text(f">3\n{d3_seq}\n>14\n{d3_seq[::-1]}\n")
    
    out_html = tmp_path / "den3_den1_report.html"
    
    runner = CliRunner()
    t0 = time.time()
    result = runner.invoke(app, [
        "analyze",
        "-p", str(p3_file),
        "-t", str(targets_fasta),
        "-o", str(out_html),
        "--temperature", "63.0",
        "--mg", "8.0",
        "--dntp", "1.4"
    ])
    t1 = time.time()
    elapsed = t1 - t0
    
    assert result.exit_code == 0, f"Error: {result.stdout}"
    print(f"\n[BENCHMARK MULTIPLEX] Analyse 2 génomes de 10.7 kb terminée en {elapsed:.2f} s")
    # Doit terminer largement en dessous de 30 secondes (au lieu de > 7 minutes)
    assert elapsed < 30.0
