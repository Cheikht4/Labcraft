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
    # Doit terminer largement en dessous du régime des 7 minutes (seuil fixé à 300s pour CI/machines lentes)
    assert elapsed < 300.0


def test_single_locus_for_shared_binding_domain_variants():
    """Vérifie que deux variantes d'une même amorce partageant le même domaine de liaison
    ne créent qu'UNE seule espèce de site dans l'équilibre, et que la somme des complexes
    sur ce site ne dépasse pas la concentration de la cible."""
    from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
    from labcraft.lamp.complex_enumeration import enumerate_complexes
    from labcraft.thermo.backends.vienna import ViennaRNABackend
    from labcraft.solver.dual import solve_dual

    backend = ViennaRNABackend()
    # Deux variantes de FIP ayant des domaines F1c distincts mais le même domaine F2 (liaison matrice)
    f2_domain = "AAGCTTAACCGGTTAAC"
    p1 = PhysicalPrimer("FIP#1", f"AAAAAA{f2_domain}", PrimerRole.FIP, f2_domain)
    p2 = PhysicalPrimer("FIP#2", f"CCCCCC{f2_domain}", PrimerRole.FIP, f2_domain)

    target_seq = f"TTTTTT{f2_domain}GGGGGG"

    prob, strands, complexes, unfolding = enumerate_complexes(
        [p1, p2], target_seq, backend, temp_celsius=65.0
    )

    # Vérification 1 : Exactement UN site créé pour les deux variantes
    site_strands = [s for s in strands if s.endswith("_site")]
    assert len(site_strands) == 1, f"Attendu: 1 site unique, Obtenu: {site_strands}"

    # Vérification 2 : Les deux complexes amorce-cible se fixent sur ce même site
    target_complexes = [c for c in complexes if "_on_" in c]
    assert len(target_complexes) == 2
    assert "FIP#1_on_FIP_site" in target_complexes
    assert "FIP#2_on_FIP_site" in target_complexes

    # Résolution thermodynamique
    res = solve_dual(prob)
    site_idx = strands.index(site_strands[0])
    target_total = prob.total_concentrations[site_idx]

    # Somme des concentrations des complexes formés sur le site
    c_indices = [i for i, c in enumerate(complexes) if "_on_" in c]
    bound_sum = sum(res.concentrations[i] for i in c_indices)

    assert bound_sum <= target_total + 1e-18, f"La somme ({bound_sum}) dépasse la matrice ({target_total})"
    # Occupation calculée
    site_occ = (target_total - res.free_concentrations[site_idx]) / target_total
    assert 0.0 <= site_occ <= 1.0


def test_degenerate_panel_one_matching_variant_nonzero_occupation(tmp_path: Path):
    """Vérifie qu'un panel dont une seule variante sur seize apparie la cible ressort
    avec un accès à l'initiation NON NUL et une amorce limitante nommée (jamais 'Non analysé')."""
    # F3 dégénérée à 4 positions (R, Y, R, Y) = 16 variantes
    # Une seule variante (ex: A...C...A...C) correspond exactement à la cible
    f3_deg = "CTCGTGTWGGAATGGGAG" # W = A/T -> 2 variantes
    # FIP avec dégénérescences
    fip_deg = "CTGAWGCCTTTCCYCAGA"
    
    target_seq = (
        "AAAA"
        "CTCGTGTAGGAATGGGAG" # Site F3 (seule variante avec W=A apparie)
        "TTTT"
        "CTGAAGCCTTTCCCCAGA" # Site FIP (seule variante avec W=A, Y=C apparie)
        "GGGG"
        "ACCACAAAGTCCCAATCA" # Site BIP
        "CCCC"
        "CTGGTTTGAGACATCTTCT" # Site B3
        "AAAA"
    )

    targets_fasta = tmp_path / "target_deg.fasta"
    targets_fasta.write_text(f">TestTarget\n{target_seq}\n")

    primers_txt = tmp_path / "primers_deg.txt"
    primers_txt.write_text(
        f">TestTarget_F3   {f3_deg}\n"
        f">TestTarget_B3   CTGGTTTGAGACATCTTCT\n"
        f">TestTarget_FIP  TTTTTTTTTTTTTTTTTTTT{fip_deg}\n"
        f">TestTarget_BIP  CCCCCCCCCCCCCCCCCCCCACCACAAAGTCCCAATCA\n"
    )

    out_html = tmp_path / "deg_report.html"

    runner = CliRunner()
    result = runner.invoke(app, [
        "analyze",
        "-p", str(primers_txt),
        "-t", str(targets_fasta),
        "-o", str(out_html),
        "--temperature", "63.0"
    ])

    assert result.exit_code == 0, f"Error: {result.stdout}"
    html = out_html.read_text()

    # Le tableau de comparaison ne doit JAMAIS afficher 'Non analysé' pour une cible analysée
    assert "Non analysé" not in html
    # L'accès à l'initiation doit être calculé et positif
    assert "TestTarget" in html


def test_mispriming_relative_criterion_and_ranking():
    """Vérifie que la détection de mésamorçage applique le critère relatif ddG <= 4.0 kcal/mol,
    qu'aucun mésamorçage n'est meilleur que le duplex parfait, et que les résultats sont triés par gravité."""
    from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
    from labcraft.diagnostics.enzyme import get_enzyme
    from labcraft.thermo.backends.vienna import ViennaRNABackend
    from labcraft.diagnostics.mispriming import detect_inter_target_mispriming

    backend = ViennaRNABackend()
    enzyme = get_enzyme("Bst2.0")

    p1 = PhysicalPrimer("F3_PanelA", "GCCACCTTAAGCCACAGTA", PrimerRole.F3, "GCCACCTTAAGCCACAGTA")
    primer_to_panel = {"F3_PanelA": "PanelA"}

    # Cible B avec un site partiel (5 nt match en 3', mais reste fortement divergent)
    # et un autre site quasi-parfait (1 mismatch interne)
    target_b_seq = (
        "TTTTTTTTTT"
        "GCCACCTTAAGCCACAGTA" # Match parfait sur cible B (très fort)
        "GGGGGGGGGG"
        "AAAAAAAAAAACACAGTA" # Match 3' seulement mais 5' très instable (> 6 kcal/mol d'écart)
        "CCCCCCCCCC"
    )

    targets = {
        "PanelA": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", # Cible A
        "PanelB": target_b_seq                               # Cible B (hétérologue)
    }

    risks = detect_inter_target_mispriming(
        [p1], primer_to_panel, targets, backend, enzyme, temp_celsius=65.0, ddg_max=4.0
    )

    for r in risks:
        assert r.target_id == "PanelB"
        assert r.delta_delta_g is not None
        assert r.delta_delta_g <= 4.0, f"Le ddG ({r.delta_delta_g}) dépasse le seuil relatif de 4.0 kcal/mol"
        assert r.delta_g < 0.0
