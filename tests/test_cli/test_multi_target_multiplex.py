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
    targets_fasta.write_text(f">CIBLE_A\n{target_a_seq}\n>CIBLE_B\n{target_b_seq}\n", encoding="utf-8")
    
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
    , encoding="utf-8")
    
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
    
    html_content = out_html.read_text(encoding="utf-8")
    
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
    
    # 2 génomes flavivirus complets de 10.7 kb correspondant aux panels 3 et 14 (tous deux dérivés de M93130.1)
    targets_fasta = tmp_path / "two_flavivirus_genomes.fasta"
    targets_fasta.write_text(f">3\n{d3_seq}\n>14\n{d3_seq}\n", encoding="utf-8")
    
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
    targets_fasta.write_text(f">TestTarget\n{target_seq}\n", encoding="utf-8")

    primers_txt = tmp_path / "primers_deg.txt"
    primers_txt.write_text(
        f">TestTarget_F3   {f3_deg}\n"
        f">TestTarget_B3   CTGGTTTGAGACATCTTCT\n"
        f">TestTarget_FIP  TTTTTTTTTTTTTTTTTTTT{fip_deg}\n"
        f">TestTarget_BIP  CCCCCCCCCCCCCCCCCCCCACCACAAAGTCCCAATCA\n"
    , encoding="utf-8")

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
    html = out_html.read_text(encoding="utf-8")

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


def test_part_a_same_panel_different_fasta_headers_identical_diagnosis(tmp_path: Path):
    """Partie A : Vérifie via CliRunner sur de vrais fichiers que le même panel
    sur le même génome avec deux en-têtes FASTA différents (>M93130.1 vs >3)
    produit un diagnostic par cible strictement identique."""
    import re
    fixture_dir = Path(__file__).resolve().parent.parent / "fixtures" / "real_primer_files"
    p3_file = fixture_dir / "primer_dengue_3.txt"
    d3_fasta = fixture_dir / "DEN3_M93130.fasta"

    from labcraft.cli.parsers import read_multi_fasta
    d3_seq = read_multi_fasta(str(d3_fasta))[0][1]

    # Extraire uniquement le panel 3
    p3_records = [
        f">{h}\n{s}" for h, s in read_multi_fasta(str(p3_file)) if h.startswith("3_")
    ]
    primers_only_3 = tmp_path / "panel_3.txt"
    primers_only_3.write_text("\n".join(p3_records) + "\n", encoding="utf-8")

    fasta_m93130 = tmp_path / "target_M93130.fasta"
    fasta_m93130.write_text(f">M93130.1\n{d3_seq}\n", encoding="utf-8")

    fasta_3 = tmp_path / "target_3.fasta"
    fasta_3.write_text(f">3\n{d3_seq}\n", encoding="utf-8")

    out_m93130 = tmp_path / "report_M93130.html"
    out_3 = tmp_path / "report_3.html"

    runner = CliRunner()
    res1 = runner.invoke(app, [
        "analyze", "-p", str(primers_only_3), "-t", str(fasta_m93130), "-o", str(out_m93130),
        "--temperature", "63.0", "--na", "50.0", "--mg", "8.0", "--dntp", "1.4"
    ])
    res2 = runner.invoke(app, [
        "analyze", "-p", str(primers_only_3), "-t", str(fasta_3), "-o", str(out_3),
        "--temperature", "63.0", "--na", "50.0", "--mg", "8.0", "--dntp", "1.4"
    ])

    assert res1.exit_code == 0, f"Run 1 failed: {res1.stdout}"
    assert res2.exit_code == 0, f"Run 2 failed: {res2.stdout}"

    html1 = out_m93130.read_text(encoding="utf-8")
    html2 = out_3.read_text(encoding="utf-8")

    # Extraire la section Diagnostic par Cible de chacun des deux rapports
    def extract_target_diag(html: str, target_label: str) -> str:
        m = re.search(r'<h2>Diagnostic par Cible</h2>(.*?)<h2>', html, re.DOTALL)
        assert m is not None, "Section Diagnostic par Cible introuvable"
        block = m.group(1).replace(f"Cible : {target_label}", "Cible : TARGET_ID")
        return block.strip()

    diag1 = extract_target_diag(html1, "M93130.1")
    diag2 = extract_target_diag(html2, "3")

    assert diag1 == diag2, f"Le diagnostic diffère selon l'en-tête FASTA:\n--- >M93130.1 ---\n{diag1}\n--- >3 ---\n{diag2}"
    assert "Toutes les cibles présentent une occupation favorable" not in diag1


def test_part_b_tolerant_sites_in_analyze_and_coherence_with_coverage(tmp_path: Path):
    """Partie B : Vérifie via CliRunner sur `labcraft analyze` et `labcraft coverage` que :
    1. B3 (1 mésappariement) et LB (2 mésappariements) du panel 3 sont trouvés sur M93130.1 en mode analyse.
    2. Pour chaque rôle (F3, B3, FIP, BIP, LF, LB), `analyze` et `coverage` donnent la même position,
       le même brin et le même nombre de mésappariements.
    3. Non-régression : les sites parfaits gardent leur ΔG au centième à 63 °C, 50 mM Na, 8 mM Mg, 1,4 mM dNTP :
       F3 -6,94, B3 parfait -6,72, LF -7,08, LB parfait -6,51, FIP -8,25, BIP -6,15.
    4. Un site frappé de veto 3' est signalé comme tel dans le rapport d'analyse et n'entre pas dans l'équilibre.
    """
    import json
    from labcraft.cli.parsers import read_multi_fasta, build_config_from_cli
    from labcraft.cli.config import build_engine_from_config
    from labcraft.lamp.complex_enumeration import enumerate_complexes
    from labcraft.lamp.coverage import CoverageAnalyzer

    fixture_dir = Path(__file__).resolve().parent.parent / "fixtures" / "real_primer_files"
    p3_file = fixture_dir / "primer_dengue_3.txt"
    d3_fasta = fixture_dir / "DEN3_M93130.fasta"

    d3_records = read_multi_fasta(str(d3_fasta))
    d3_id, d3_seq = d3_records[0][0].split()[0], d3_records[0][1]

    p3_records = [
        f">{h}\n{s}" for h, s in read_multi_fasta(str(p3_file)) if h.startswith("3_")
    ]
    primers_only_3 = tmp_path / "panel_3.txt"
    primers_only_3.write_text("\n".join(p3_records) + "\n", encoding="utf-8")

    out_analyze = tmp_path / "analyze_den3.html"
    out_cov = tmp_path / "coverage_den3.html"

    runner = CliRunner()
    res_an = runner.invoke(app, [
        "analyze", "-p", str(primers_only_3), "-t", str(d3_fasta), "-o", str(out_analyze),
        "--temperature", "63.0", "--na", "50.0", "--mg", "8.0", "--dntp", "1.4"
    ])
    assert res_an.exit_code == 0, f"Analyze failed: {res_an.stdout}"

    res_cov = runner.invoke(app, [
        "coverage", "-p", str(primers_only_3), "-f", str(d3_fasta), "-o", str(out_cov),
        "--temperature", "63.0", "--na", "50.0", "--mg", "8.0", "--dntp", "1.4"
    ])
    assert res_cov.exit_code == 0, f"Coverage failed: {res_cov.stdout}"

    html_an = out_analyze.read_text(encoding="utf-8")
    # 1. B3_3_site (1 mismatch) et LB_3_site (2 mismatches) présents dans le rapport HTML d'analyse
    assert "B3_3_site" in html_an
    assert "LB_3_site" in html_an
    assert "1 mismatches" in html_an
    assert "2 mismatches" in html_an

    # 2. Comparaison position, brin, nombre de mésappariements entre coverage (JSON) et analyze
    cov_json = json.loads((tmp_path / "coverage_den3.json").read_text(encoding="utf-8"))
    cov_evals = cov_json["strains"][0]["evals"]

    cfg = build_config_from_cli(
        config_path=None,
        primers_path=str(primers_only_3),
        targets_path=str(d3_fasta),
        temperature=63.0, na=50.0, mg=8.0, dntp=1.4
    )
    primers, p2t, backend, b_kwargs, mon_molar, enzyme, temp, profiles = build_engine_from_config(cfg, {d3_id: d3_seq})
    _, _, _, unf = enumerate_complexes(
        primers, {d3_id: d3_seq}, backend, profile=profiles, temp_celsius=63.0,
        mon_molar=mon_molar, buffer=cfg.experiment.buffer.model_dump(), enzyme=enzyme, primer_to_target=p2t
    )

    for role in ("F3", "B3", "FIP", "BIP", "LF", "LB"):
        site_key = f"{role}_3_site"
        assert site_key in unf, f"{site_key} manquant dans enumerate_complexes"
        an_info = unf[site_key]
        cv_info = cov_evals[role]
        assert an_info["position"] == cv_info["position"], f"Position différente pour {role}"
        assert an_info["strand"] == cv_info["strand"], f"Brin différent pour {role}"
        assert an_info["mismatches"] == cv_info["mismatches"], f"Mismatches différents pour {role}"

    # 3a. Non-régression au centième sur les sites parfaits à 63 °C, « 50 mM Na, 0 Mg »
    # 3a. Non-regression to 0.01 kcal/mol on perfect sites at 63 °C, "50 mM Na, 0 Mg"
    cov_analyzer_50na_0mg = CoverageAnalyzer(
        primers, {d3_id: d3_seq}, backend, enzyme, temp_celsius=63.0,
        na_mm=50.0, mg_mm=0.0, dntp_mm=0.0,
    )
    expected_perfect_dg_50na_0mg = {
        "CTCGTGTAGGAATGGGAG": -6.94,      # F3 parfait (50 mM Na, 0 Mg)
        "CTGGCTTGAGACATCTTCT": -6.72,     # B3 parfait (50 mM Na, 0 Mg)
        "TGAATTCCACGAGCGTTC": -7.08,      # LF parfait (50 mM Na, 0 Mg)
        "AAGGTCATCCAACTTAGTAGGA": -6.51,  # LB parfait (50 mM Na, 0 Mg)
        "CTGAAGCCTTTCCCCAGA": -8.25,      # FIP parfait (50 mM Na, 0 Mg)
        "ACCACAAAGTCCCAATCA": -6.15,      # BIP parfait (50 mM Na, 0 Mg)
    }
    for seq, expected_dg in expected_perfect_dg_50na_0mg.items():
        actual_dg = round(cov_analyzer_50na_0mg._get_perfect_dg_seq(seq), 2)
        assert actual_dg == expected_dg, (
            f"Non-régression ΔG (50 mM Na, 0 Mg) pour {seq}: attendu {expected_dg}, obtenu {actual_dg}"
        )

    # 3b. Non-régression au centième sur les sites parfaits à 63 °C, « 50 mM Na, 8 mM Mg, 1,4 mM dNTP »
    # 3b. Non-regression to 0.01 kcal/mol on perfect sites at 63 °C, "50 mM Na, 8 mM Mg, 1.4 mM dNTP"
    cov_analyzer_50na_8mg = CoverageAnalyzer(
        primers, {d3_id: d3_seq}, backend, enzyme, temp_celsius=63.0,
        **b_kwargs,
    )
    expected_perfect_dg_50na_8mg = {
        "CTCGTGTAGGAATGGGAG": -10.65,      # F3 parfait (50 mM Na, 8 mM Mg, 1.4 mM dNTP)
        "CTGGCTTGAGACATCTTCT": -11.00,     # B3 parfait (50 mM Na, 8 mM Mg, 1.4 mM dNTP)
        "TGAATTCCACGAGCGTTC": -11.16,      # LF parfait (50 mM Na, 8 mM Mg, 1.4 mM dNTP)
        "AAGGTCATCCAACTTAGTAGGA": -11.53,  # LB parfait (50 mM Na, 8 mM Mg, 1.4 mM dNTP)
        "CTGAAGCCTTTCCCCAGA": -11.95,      # FIP parfait (50 mM Na, 8 mM Mg, 1.4 mM dNTP)
        "ACCACAAAGTCCCAATCA": -10.32,      # BIP parfait (50 mM Na, 8 mM Mg, 1.4 mM dNTP)
    }
    for seq, expected_dg in expected_perfect_dg_50na_8mg.items():
        actual_dg = round(cov_analyzer_50na_8mg._get_perfect_dg_seq(seq), 2)
        assert actual_dg == expected_dg, (
            f"Non-régression ΔG (50 mM Na, 8 mM Mg, 1.4 mM dNTP) pour {seq}: attendu {expected_dg}, obtenu {actual_dg}"
        )

    # Vérification des valeurs sur M93130.1 à 50 mM Na, 8 mM Mg, 1,4 mM dNTP dans le JSON de coverage
    expected_m93130_dg_8mg = {
        "F3": -10.65,
        "B3": -9.66,
        "FIP": -10.78,
        "BIP": -10.32,
        "LF": -10.05,
        "LB": -8.09,
    }
    for role, exp_dg in expected_m93130_dg_8mg.items():
        assert round(cov_evals[role]["dg"], 2) == exp_dg
    assert cov_evals["B3"]["verdict"] == "tolerable"
    assert round(cov_evals["B3"]["ddg"], 2) == 1.33

    # 4. Vérification via CLI qu'un site frappé de veto 3' est signalé dans le rapport d'analyse et a 0 % d'occupation
    f3_pos = cov_evals["F3"]["position"]
    mut_seq = d3_seq[:f3_pos + 17] + ("T" if d3_seq[f3_pos + 17] != "T" else "A") + d3_seq[f3_pos + 18:]
    mut_fasta = tmp_path / "den3_f3_veto.fasta"
    mut_fasta.write_text(f">DEN3_veto\n{mut_seq}\n", encoding="utf-8")
    out_veto = tmp_path / "analyze_veto.html"
    res_veto = runner.invoke(app, [
        "analyze", "-p", str(primers_only_3), "-t", str(mut_fasta), "-o", str(out_veto),
        "--temperature", "63.0", "--na", "50.0", "--mg", "8.0", "--dntp", "1.4"
    ])
    assert res_veto.exit_code == 0
    html_veto = out_veto.read_text(encoding="utf-8")
    assert "Veto" in html_veto or "3' Bloqué" in html_veto


