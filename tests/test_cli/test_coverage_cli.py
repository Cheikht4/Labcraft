import pytest
from typer.testing import CliRunner
import tempfile
import csv
import json
import subprocess
import sys
import time
import random
from pathlib import Path

from labcraft.cli.main import app
from labcraft.cli.config import PanelConfig, ExperimentConfig, BufferConfig
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.lamp.coverage import CoverageAnalyzer, SiteVerdict
from labcraft.target.seeding import find_candidate_sites, build_primer_regex
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.thermo.salt import UnifiedSaltModel
from labcraft.diagnostics.enzyme import get_enzyme


def test_coverage_cli_with_csv(tmp_path: Path):
    """Test du mode couverture avec fichier CSV de sites pré-calculé."""
    fasta_file = tmp_path / "strains.fasta"
    dummy_seq = "GAAGAAGCTGTGCAGCCTG" + "N"*20 + "TGGCTTTTGGGCCTGACTTC" + "N"*20 + "GCACGGTGTAGCAGACTAG" + "N"*20 + "TCCCCACGACGGAGCTACAG"
    
    fasta_text = ">S1_parfaite\n" + dummy_seq + "\n>S2_mism_interne\n" + dummy_seq + "\n>S3_veto_3p\n" + dummy_seq + "\n>S4_loop_fail\n" + dummy_seq + "\n>S5_count_lost\n" + dummy_seq + "\n"
    fasta_file.write_text(fasta_text, encoding="utf-8")
    
    csv_file = tmp_path / "sites.csv"
    with open(csv_file, "w", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["strain_id", "primer_role", "primer_name", "position", "strand", "site_seq", "n_mismatches", "panel"])
        writer.writeheader()
        
        strains = ["S1_parfaite", "S2_mism_interne", "S3_veto_3p", "S4_loop_fail", "S5_count_lost"]
        
        for s in strains:
            f3_seq = "GCCACCTTAAGCCACAGTA"
            f3_mm = 0
            if s == "S2_mism_interne":
                f3_seq = "GCCACCTTTGCCACAGTA"  # Benign mismatch (ddG <= 3.0)
                f3_mm = 1
            elif s == "S3_veto_3p":
                f3_seq = "GCCACCTTAAGCCACAGTT"  # 3' mismatch -> veto 3'
                f3_mm = 1
            elif s == "S5_count_lost":
                f3_mm = 3
                
            writer.writerow({"strain_id": s, "primer_role": "F3", "primer_name": "Test_F3", "position": 0, "strand": "+", "site_seq": f3_seq, "n_mismatches": f3_mm, "panel": "Test"})
            writer.writerow({"strain_id": s, "primer_role": "B3", "primer_name": "Test_B3", "position": 0, "strand": "-", "site_seq": "CCCTCCCATGACACAAC", "n_mismatches": 0, "panel": "Test"})
            writer.writerow({"strain_id": s, "primer_role": "FIP", "primer_name": "Test_FIP", "position": 0, "strand": "+", "site_seq": "GAAGAAGCTGTGCAGCCTG", "n_mismatches": 0, "panel": "Test"})
            writer.writerow({"strain_id": s, "primer_role": "BIP", "primer_name": "Test_BIP", "position": 0, "strand": "+", "site_seq": "CTAGTCTGCTACACCGTGC", "n_mismatches": 0, "panel": "Test"})
            
            lf_seq = "CCTTGGACGGGGAA" if s == "S4_loop_fail" else "CCTTGGACGGGGCT"
            lf_mm = 2 if s == "S4_loop_fail" else 0
            writer.writerow({"strain_id": s, "primer_role": "LF", "primer_name": "Test_LF", "position": 0, "strand": "+", "site_seq": lf_seq, "n_mismatches": lf_mm, "panel": "Test"})
            
    primers_file = tmp_path / "primers.fasta"
    primers_file.write_text(">Test_F3\nGCCACCTTAAGCCACAGTA\n>Test_B3\nGTTGTGTCATGGGAGGG\n>Test_F1\nGAAGTCAGGCCCAAAAGCCA\n>Test_F2\nGAAGAAGCTGTGCAGCCTG\n>Test_B1\nTCCCCACGACGGAGCTACAG\n>Test_B2\nCTAGTCTGCTACACCGTGC\n>Test_LF\nCCTTGGACGGGGCT\n", encoding="utf-8")
    
    out_html = tmp_path / "report.html"
    
    runner = CliRunner()
    result = runner.invoke(app, [
        "coverage",
        "-p", str(primers_file),
        "-f", str(fasta_file),
        "-s", str(csv_file),
        "-o", str(out_html),
        "--temperature", "63.0",
        "--ddg-max", "3.0",
        "--max-mismatches", "2"
    ])
    
    assert result.exit_code == 0, f"Error: {result.stdout}"
    
    html_content = out_html.read_text(encoding="utf-8")
    assert "S1_parfaite" in html_content
    assert "S2_mism_interne" in html_content
    assert "S3_veto_3p" in html_content
    assert "S5_count_lost" in html_content
    assert "S4_loop_fail" in html_content
    
    assert "Divergences Comptage vs Thermodynamique" in html_content
    assert "S3_veto_3p" in html_content
    assert "S5_count_lost" in html_content


def test_coverage_cli_auto_seeding_and_export(tmp_path: Path):
    """Test du criblage intégré (sans -s) et vérification de la ré-exécutabilité via --export-sites."""
    fasta_file = tmp_path / "strains.fasta"
    # Séquence contenant F3, F2, LF, F1c(-), B1c(+), B2(RC), B3(RC)
    seq_s1 = (
        "GCCACCTTAAGCCACAGTA" + "A"*30 +
        "GAAGAAGCTGTGCAGCCTG" + "A"*30 +
        "CCTTGGACGGGGCT" + "A"*30 +
        "TGGCTTTTGGGCCTGACTTC" + "A"*30 +
        "TCCCCACGACGGAGCTACAG" + "A"*30 +
        "GCACGGTGTAGCAGACTAG" + "A"*30 +
        "CCCTCCCATGACACAAC"
    )
    fasta_file.write_text(f">Strain_1\n{seq_s1}\n", encoding="utf-8")
    
    primers_file = tmp_path / "primers.fasta"
    primers_file.write_text(
        ">Test_F3\nGCCACCTTAAGCCACAGTA\n"
        ">Test_B3\nGTTGTGTCATGGGAGGG\n"
        ">Test_F1\nGAAGTCAGGCCCAAAAGCCA\n"
        ">Test_F2\nGAAGAAGCTGTGCAGCCTG\n"
        ">Test_B1\nTCCCCACGACGGAGCTACAG\n"
        ">Test_B2\nCTAGTCTGCTACACCGTGC\n"
        ">Test_LF\nCCTTGGACGGGGCT\n"
    , encoding="utf-8")
    
    out_html1 = tmp_path / "report1.html"
    exported_csv = tmp_path / "exported_sites.csv"
    
    runner = CliRunner()
    # 1. Exécution sans -s avec --export-sites
    res1 = runner.invoke(app, [
        "coverage",
        "-p", str(primers_file),
        "-f", str(fasta_file),
        "-o", str(out_html1),
        "--export-sites", str(exported_csv),
        "--temperature", "63.0",
        "--ddg-max", "3.0"
    ])
    assert res1.exit_code == 0, f"Error res1: {res1.stdout}"
    assert exported_csv.exists()
    
    # Vérifier que le CSV exporté contient la colonne primer_name
    with open(exported_csv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        assert "primer_name" in reader.fieldnames
        rows = list(reader)
        assert len(rows) > 0
        assert all(row["primer_name"] != "" for row in rows)
    
    # 2. Exécution avec le CSV exporté
    out_html2 = tmp_path / "report2.html"
    res2 = runner.invoke(app, [
        "coverage",
        "-p", str(primers_file),
        "-f", str(fasta_file),
        "-s", str(exported_csv),
        "-o", str(out_html2),
        "--temperature", "63.0",
        "--ddg-max", "3.0"
    ])
    assert res2.exit_code == 0, f"Error res2: {res2.stdout}"
    
    # Comparaison des résultats JSON
    json1 = json.loads((tmp_path / "report1.json").read_text(encoding="utf-8"))
    json2 = json.loads((tmp_path / "report2.json").read_text(encoding="utf-8"))
    assert json1["covered_thermo"] == json2["covered_thermo"]
    assert json1["covered_count"] == json2["covered_count"]
    assert json1["strains"] == json2["strains"]


def test_coverage_reference_amplifiable_default_options(tmp_path: Path):
    """Vérifie qu'avec les options PAR DÉFAUT (sans aucun seuil précisé),
    le génome de référence M93130.1 est déclaré amplifiable pour le panel DEN-3.
    """
    fixture_dir = Path(__file__).resolve().parent.parent / "fixtures" / "real_primer_files"
    primers_path = fixture_dir / "primer_dengue_3.txt"
    genome_path = fixture_dir / "DEN3_M93130.fasta"
    
    out_html = tmp_path / "coverage_den3_default.html"
    
    runner = CliRunner()
    result = runner.invoke(app, [
        "coverage",
        "-p", str(primers_path),
        "-f", str(genome_path),
        "--panel", "3",
        "-o", str(out_html)
    ])
    
    assert result.exit_code == 0, f"Error: {result.stdout}"
    
    json_data = json.loads((tmp_path / "coverage_den3_default.json").read_text(encoding="utf-8"))
    assert json_data["total_strains"] == 1
    assert json_data["covered_thermo"] == 1, "Le génome de référence doit être amplifiable par défaut."
    assert json_data["covered_count"] == 1
    
    strain = json_data["strains"][0]
    assert strain["thermo"] is True
    assert strain["count"] is True
    # F3 doit être parfait (0 mm)
    assert strain["evals"]["F3"]["verdict"] == "parfait"
    # B3 avec 1 mismatch interne bénin doit être tolérable
    assert strain["evals"]["B3"]["verdict"] == "tolerable"
    assert strain["evals"]["B3"]["ddg"] is not None
    assert strain["evals"]["B3"]["ddg"] <= 3.0


def test_seeding_indel_rejection():
    """Vérifie que le motif de criblage rejette les insertions et délétions (substitutions seules)."""
    # var004-like : insertion d'un A supplémentaire dans le site F3 (longueur 19 au lieu de 18)
    # F3 = CTCGTGTWGGAATGGGAG (18 nt)
    # Cible mutée avec insertion : CTCGTGTAGGAATGGAAGAG (20 nt)
    f3_primer = PhysicalPrimer('3_F3', 'CTCGTGTWGGAATGGGAG', PrimerRole.F3, 'CTCGTGTWGGAATGGGAG')
    target_with_insertion = "A"*50 + "CTCGTGTAGGAATGGAAGAG" + "A"*50
    
    strains = {"var004": target_with_insertion}
    candidates = find_candidate_sites(strains, [f3_primer], max_errors=2, strict_3prime_len=3)
    
    for c in candidates:
        # Aucun site de longueur différente de 18 nt ne doit être retourné
        assert len(c["site_seq"]) == 18, f"Le site trouvé a une longueur incorrecte: {len(c['site_seq'])}"


def test_independent_counting_rule():
    """Vérifie que la règle par comptage prend le minimum de mésappariements sur toutes les lignes,
    indépendamment du choix de la variante thermodynamique."""
    # Variante 1: 1 mismatch bloquant en 3' (VETO_3P, mm=1)
    # Variante 2: 1 mismatch interne bénin (TOLERABLE, mm=2 dans CSV pour tester l'indépendance)
    p1 = PhysicalPrimer('F3_1', 'GCCACCTTAAGCCACAGTT', PrimerRole.F3, 'GCCACCTTAAGCCACAGTT')
    p2 = PhysicalPrimer('F3_2', 'GCCACCTTTAGCCACAGTA', PrimerRole.F3, 'GCCACCTTTAGCCACAGTA')
    
    primers = [
        p1, p2,
        PhysicalPrimer('B3', 'GTTGTGTCATGGGAGGG', PrimerRole.B3, 'GTTGTGTCATGGGAGGG'),
        PhysicalPrimer('FIP', 'TGGCTTTTGGGCCTGACTTCTTTTTTGAAGAAGCTGTGCAGCCTG', PrimerRole.FIP, 'GAAGAAGCTGTGCAGCCTG', 'TGGCTTTTGGGCCTGACTTC'),
        PhysicalPrimer('BIP', 'CTGTAGCTCCGTCGTGGGGATTTTCTAGTCTGCTACACCGTGC', PrimerRole.BIP, 'CTAGTCTGCTACACCGTGC', 'CTGTAGCTCCGTCGTGGGGA'),
    ]
    
    fasta_dict = {"Strain_1": "GCCACCTTAAGCCACAGTA" + "A"*50}
    
    csv_records = [
        {"strain_id": "Strain_1", "primer_role": "F3", "primer_name": "F3_1", "position": 0, "strand": "+", "site_seq": "GCCACCTTAAGCCACAGTA", "n_mismatches": 1},
        {"strain_id": "Strain_1", "primer_role": "F3", "primer_name": "F3_2", "position": 0, "strand": "+", "site_seq": "GCCACCTTAAGCCACAGTA", "n_mismatches": 2},
        {"strain_id": "Strain_1", "primer_role": "B3", "position": 0, "strand": "-", "site_seq": "CCCTCCCATGACACAAC", "n_mismatches": 0},
        {"strain_id": "Strain_1", "primer_role": "FIP", "position": 0, "strand": "+", "site_seq": "GAAGAAGCTGTGCAGCCTG", "n_mismatches": 0},
        {"strain_id": "Strain_1", "primer_role": "BIP", "position": 0, "strand": "+", "site_seq": "CTAGTCTGCTACACCGTGC", "n_mismatches": 0},
    ]
    
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    enzyme = get_enzyme("Bst2.0")
    analyzer = CoverageAnalyzer(primers, fasta_dict, backend, enzyme, temp_celsius=63.0, ddg_max=3.0, max_mismatches_count=1)
    
    verdicts = analyzer.analyze_strains(csv_records)
    assert len(verdicts) == 1
    v = verdicts[0]
    
    # Thermodynamiquement: F3_2 est choisie (TOLERABLE l'emporte sur VETO_3P)
    # Comptage: min(1, 2) = 1 <= max_mismatches_count(1) -> amplifiable_count = True
    assert v.evaluations["F3"].verdict == SiteVerdict.TOLERABLE
    assert v.is_amplifiable_thermo is True
    assert v.is_amplifiable_count is True


def test_divergence_two_ways():
    """Vérifie que la table de divergence peut avoir des souches dans les deux sens :
    - Thermo=True, Count=False (ex: 3 mismatches très bénins ddG <= 3.0 avec seuil count=2)
    - Thermo=False, Count=True (ex: 1 mismatch bloquant en 3' veto_3p avec seuil count=2)
    """
    p_f3 = PhysicalPrimer('F3', 'GCCACCTTAAGCCACAGTA', PrimerRole.F3, 'GCCACCTTAAGCCACAGTA')
    primers = [
        p_f3,
        PhysicalPrimer('B3', 'GTTGTGTCATGGGAGGG', PrimerRole.B3, 'GTTGTGTCATGGGAGGG'),
        PhysicalPrimer('FIP', 'TGGCTTTTGGGCCTGACTTCTTTTTTGAAGAAGCTGTGCAGCCTG', PrimerRole.FIP, 'GAAGAAGCTGTGCAGCCTG', 'TGGCTTTTGGGCCTGACTTC'),
        PhysicalPrimer('BIP', 'CTGTAGCTCCGTCGTGGGGATTTTCTAGTCTGCTACACCGTGC', PrimerRole.BIP, 'CTAGTCTGCTACACCGTGC', 'CTGTAGCTCCGTCGTGGGGA'),
    ]
    
    fasta_dict = {
        "Strain_Veto": "GCCACCTTAAGCCACAGTT" + "A"*50,
        "Strain_ThermoWon": "GCCACCTTTAGCCACAGTA" + "A"*50,
    }
    
    csv_records = [
        # Strain_Veto: 1 mm en 3' -> Veto 3' (Thermo=False, Count=True)
        {"strain_id": "Strain_Veto", "primer_role": "F3", "position": 0, "strand": "+", "site_seq": "GCCACCTTAAGCCACAGTT", "n_mismatches": 1},
        {"strain_id": "Strain_Veto", "primer_role": "B3", "position": 0, "strand": "-", "site_seq": "CCCTCCCATGACACAAC", "n_mismatches": 0},
        {"strain_id": "Strain_Veto", "primer_role": "FIP", "position": 0, "strand": "+", "site_seq": "GAAGAAGCTGTGCAGCCTG", "n_mismatches": 0},
        {"strain_id": "Strain_Veto", "primer_role": "BIP", "position": 0, "strand": "+", "site_seq": "CTAGTCTGCTACACCGTGC", "n_mismatches": 0},
        
        # Strain_ThermoWon: 3 mm déclarés au comptage mais ddg <= 3.0 -> Thermo=True, Count=False (si max_mismatches=2)
        {"strain_id": "Strain_ThermoWon", "primer_role": "F3", "position": 0, "strand": "+", "site_seq": "GCCACCTTTAGCCACAGTA", "n_mismatches": 3},
        {"strain_id": "Strain_ThermoWon", "primer_role": "B3", "position": 0, "strand": "-", "site_seq": "CCCTCCCATGACACAAC", "n_mismatches": 0},
        {"strain_id": "Strain_ThermoWon", "primer_role": "FIP", "position": 0, "strand": "+", "site_seq": "GAAGAAGCTGTGCAGCCTG", "n_mismatches": 0},
        {"strain_id": "Strain_ThermoWon", "primer_role": "BIP", "position": 0, "strand": "+", "site_seq": "CTAGTCTGCTACACCGTGC", "n_mismatches": 0},
    ]
    
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    enzyme = get_enzyme("Bst2.0")
    analyzer = CoverageAnalyzer(primers, fasta_dict, backend, enzyme, temp_celsius=63.0, ddg_max=3.0, max_mismatches_count=2)
    verdicts = {v.strain_id: v for v in analyzer.analyze_strains(csv_records)}
    
    # Sens 1 : Count=True, Thermo=False
    assert verdicts["Strain_Veto"].is_amplifiable_count is True
    assert verdicts["Strain_Veto"].is_amplifiable_thermo is False
    
    # Sens 2 : Count=False, Thermo=True
    assert verdicts["Strain_ThermoWon"].is_amplifiable_count is False
    assert verdicts["Strain_ThermoWon"].is_amplifiable_thermo is True


def test_cli_mg_dntp_warning_via_runner(tmp_path: Path):
    """Vérifie que l'avertissement de chélation Mg/dNTP est émis via la CLI."""
    fixture_dir = Path(__file__).resolve().parent.parent / "fixtures" / "real_primer_files"
    primers_path = fixture_dir / "primer_dengue_3.txt"
    genome_path = fixture_dir / "DEN3_M93130.fasta"
    
    out_html = tmp_path / "coverage_warn.html"
    
    runner = CliRunner()
    with pytest.warns(UserWarning, match="Rapport Mg/dNTP faible"):
        result = runner.invoke(app, [
            "coverage",
            "-p", str(primers_path),
            "-f", str(genome_path),
            "--panel", "3",
            "-o", str(out_html),
            "--mg", "2.8",
            "--dntp", "1.4"
        ])
        assert result.exit_code == 0


def test_coverage_200_strains_real_performance(tmp_path: Path):
    """Test sur 201 souches mutées à 1 % dérivées de DEN-3 M93130.1 complet."""
    fixture_dir = Path(__file__).resolve().parent.parent / "fixtures" / "real_primer_files"
    primers_path = fixture_dir / "primer_dengue_3.txt"
    genome_path = fixture_dir / "DEN3_M93130.fasta"
    
    from labcraft.cli.parsers import read_multi_fasta
    ref_seq = read_multi_fasta(str(genome_path))[0][1]
    
    random.seed(42)
    strains_lines = [f">M93130_ref\n{ref_seq}"]
    
    bases = ['A', 'C', 'G', 'T']
    # Générer 200 souches mutées à 1%
    for i in range(1, 201):
        seq_chars = list(ref_seq)
        n_muts = int(len(seq_chars) * 0.01)
        mut_indices = random.sample(range(len(seq_chars)), n_muts)
        for idx in mut_indices:
            orig = seq_chars[idx]
            choices = [b for b in bases if b != orig]
            seq_chars[idx] = random.choice(choices)
        strains_lines.append(f">Strain_var_{i:03d}\n{''.join(seq_chars)}")
        
    multi_fasta_path = tmp_path / "201_strains_dengue.fasta"
    multi_fasta_path.write_text("\n".join(strains_lines) + "\n", encoding="utf-8")
    
    out_html = tmp_path / "coverage_201_strains.html"
    
    runner = CliRunner()
    t0 = time.time()
    result = runner.invoke(app, [
        "coverage",
        "-p", str(primers_path),
        "-f", str(multi_fasta_path),
        "--panel", "3",
        "-o", str(out_html),
        "--temperature", "63.0"
    ])
    t1 = time.time()
    
    assert result.exit_code == 0, f"Error: {result.stdout}"
    total_time = t1 - t0
    
    json_data = json.loads((tmp_path / "coverage_201_strains.json").read_text(encoding="utf-8"))
    assert json_data["total_strains"] == 201
    assert json_data["covered_thermo"] >= 1
    assert json_data["covered_count"] >= 1
    
    print(f"\n[BENCHMARK] 201 souches DEN-3 analysées en {total_time:.2f} s. "
          f"Couverture Thermo: {json_data['covered_thermo']}/201 ({json_data['covered_thermo']/201*100:.1f}%), "
          f"Comptage: {json_data['covered_count']}/201 ({json_data['covered_count']/201*100:.1f}%)")


def test_veto_3p_reachable_single_command(tmp_path: Path):
    """Vérifie qu'une souche mutée sur la base 3' terminale d'une amorce d'initiation (ex: F3),
    analysée par labcraft coverage sans option de seuil particulière, ressort VETO_3P (et non absent)."""
    fixture_dir = Path(__file__).resolve().parent.parent / "fixtures" / "real_primer_files"
    primers_path = fixture_dir / "primer_dengue_3.txt"
    genome_path = fixture_dir / "DEN3_M93130.fasta"
    
    from labcraft.cli.parsers import read_multi_fasta
    ref_seq = read_multi_fasta(str(genome_path))[0][1]
    
    # F3 se situe en pos 5420 sur M93130.1: CTCGTGTAGGAATGGGAG (18 nt)
    # Mutons la base terminale en 3' (position 5420 + 17 = 5437) de G vers T
    f3_pos = 5420
    f3_term_pos = f3_pos + 17
    assert ref_seq[f3_term_pos] == 'G'
    mutated_seq = ref_seq[:f3_term_pos] + 'T' + ref_seq[f3_term_pos + 1:]
    
    fasta_file = tmp_path / "mutated_f3_term.fasta"
    fasta_file.write_text(f">DEN3_F3_3p_mut\n{mutated_seq}\n", encoding="utf-8")
    
    out_html = tmp_path / "coverage_f3_term.html"
    
    runner = CliRunner()
    result = runner.invoke(app, [
        "coverage",
        "-p", str(primers_path),
        "-f", str(fasta_file),
        "--panel", "3",
        "-o", str(out_html),
        "--temperature", "63.0"
    ])
    
    assert result.exit_code == 0, f"Error: {result.stdout}"
    
    json_data = json.loads((tmp_path / "coverage_f3_term.json").read_text(encoding="utf-8"))
    strain_eval = json_data["strains"][0]["evals"]["F3"]
    assert strain_eval["verdict"] == "veto_3p", f"Attendu: veto_3p, Obtenu: {strain_eval['verdict']}"
    assert strain_eval["first_bad_pos"] == 1
    assert strain_eval["severity"] == "block"
    assert strain_eval["position"] == f3_pos
    assert strain_eval["strand"] == "+"
    assert json_data["strains"][0]["thermo"] is False


def test_absent_site_when_heavily_mutated(tmp_path: Path):
    """Vérifie qu'une souche avec plus de 2 substitutions dans la zone 5' ressort bien ABSENT."""
    fixture_dir = Path(__file__).resolve().parent.parent / "fixtures" / "real_primer_files"
    primers_path = fixture_dir / "primer_dengue_3.txt"
    genome_path = fixture_dir / "DEN3_M93130.fasta"
    
    from labcraft.cli.parsers import read_multi_fasta
    ref_seq = read_multi_fasta(str(genome_path))[0][1]
    
    # Mutons 4 bases dans le site F3 (pos 5420)
    f3_pos = 5420
    mutated_chars = list(ref_seq)
    mutated_chars[f3_pos] = 'A' if ref_seq[f3_pos] != 'A' else 'T'
    mutated_chars[f3_pos + 1] = 'A' if ref_seq[f3_pos + 1] != 'A' else 'T'
    mutated_chars[f3_pos + 2] = 'A' if ref_seq[f3_pos + 2] != 'A' else 'T'
    mutated_chars[f3_pos + 3] = 'A' if ref_seq[f3_pos + 3] != 'A' else 'T'
    
    fasta_file = tmp_path / "heavily_mutated.fasta"
    fasta_file.write_text(f">DEN3_heavy_mut\n{''.join(mutated_chars)}\n", encoding="utf-8")
    
    out_html = tmp_path / "coverage_heavy.html"
    
    runner = CliRunner()
    result = runner.invoke(app, [
        "coverage",
        "-p", str(primers_path),
        "-f", str(fasta_file),
        "--panel", "3",
        "-o", str(out_html),
        "--temperature", "63.0"
    ])
    
    assert result.exit_code == 0, f"Error: {result.stdout}"
    json_data = json.loads((tmp_path / "coverage_heavy.json").read_text(encoding="utf-8"))
    strain_eval = json_data["strains"][0]["evals"]["F3"]
    assert strain_eval["verdict"] == "absent"


def test_loop_primer_3prime_mutation_remains_amplifiable(tmp_path: Path):
    """Vérifie qu'une mutation 3' terminale sur une amorce de boucle (LF/LB)
    donne veto_3p sur l'amorce de boucle mais permet à la souche de rester amplifiable."""
    fixture_dir = Path(__file__).resolve().parent.parent / "fixtures" / "real_primer_files"
    primers_path = fixture_dir / "primer_dengue_3.txt"
    genome_path = fixture_dir / "DEN3_M93130.fasta"
    
    from labcraft.cli.parsers import read_multi_fasta
    ref_seq = read_multi_fasta(str(genome_path))[0][1]
    
    # LF est sur brin - en pos 5536 (18 nt). Mutons sa base terminale 3'
    # LF primer = TGAATTCCAYGAGCGTTC
    # Sur génome (brin +), le 3' terminal de LF correspond au complément inverse de la base 5' du site ou 3' du primer
    lf_pos = 5536
    mutated_chars = list(ref_seq)
    mutated_chars[lf_pos] = 'A' if ref_seq[lf_pos] != 'A' else 'T'
    
    fasta_file = tmp_path / "mutated_lf.fasta"
    fasta_file.write_text(f">DEN3_LF_mut\n{''.join(mutated_chars)}\n", encoding="utf-8")
    
    out_html = tmp_path / "coverage_lf.html"
    
    runner = CliRunner()
    result = runner.invoke(app, [
        "coverage",
        "-p", str(primers_path),
        "-f", str(fasta_file),
        "--panel", "3",
        "-o", str(out_html),
        "--temperature", "63.0"
    ])
    
    assert result.exit_code == 0, f"Error: {result.stdout}"
    json_data = json.loads((tmp_path / "coverage_lf.json").read_text(encoding="utf-8"))
    strain = json_data["strains"][0]
    # LF a un veto_3p mais l'amplification globale reste True
    assert strain["evals"]["LF"]["verdict"] == "veto_3p"
    assert strain["thermo"] is True
    assert strain["count"] is True


def test_coverage_analyzer_warning_when_no_buffer():
    """Vérifie qu'un avertissement est émis si CoverageAnalyzer est construit sans ViennaSaltShiftBackend."""
    from labcraft.thermo.backends.vienna import ViennaRNABackend
    p = PhysicalPrimer('F3', 'GCCACCTTAAGCCACAGTA', PrimerRole.F3, 'GCCACCTTAAGCCACAGTA')
    enzyme = get_enzyme("Bst2.0")
    
    with pytest.warns(UserWarning, match="Aucun modèle de sel/tampon"):
        CoverageAnalyzer([p], {"S1": "GCCACCTTAAGCCACAGTA"}, ViennaRNABackend(), enzyme)


def test_part_h_three_ns_in_f3_site_gives_indetermine_and_preserves_coverage_rate(tmp_path: Path):
    """
    Partie H : Vérifie qu'une souche contenant 3 N au milieu du site F3 donne le verdict
    'indetermine' (et non 'absent'), ne fait pas baisser le taux de couverture, et
    affiche le statut indéterminé en mode analyze.
    Part H: Verify that a strain with 3 N bases in the middle of the F3 site gets the
    'indetermine' verdict (not 'absent'), does not lower the coverage rate, and
    displays the indeterminate status in analyze mode.
    """
    fixture_dir = Path(__file__).resolve().parent.parent / "fixtures" / "real_primer_files"
    primers_path = fixture_dir / "primer_dengue_3.txt"
    genome_path = fixture_dir / "DEN3_M93130.fasta"

    from labcraft.cli.parsers import read_multi_fasta
    ref_seq = read_multi_fasta(str(genome_path))[0][1]

    # Site F3 sur M93130.1 en position 5420 (18 nt) ; on remplace 3 bases au milieu par NNN
    # F3 site on M93130.1 at position 5420 (18 nt); replace 3 bases in the middle with NNN
    f3_pos = 5420
    seq_with_3n = ref_seq[:f3_pos + 7] + "NNN" + ref_seq[f3_pos + 10:]

    fasta_2strains = tmp_path / "strains_with_3n.fasta"
    fasta_2strains.write_text(f">M93130_ref\n{ref_seq}\n>M93130_F3_3N\n{seq_with_3n}\n", encoding="utf-8")

    out_cov_html = tmp_path / "coverage_3n.html"
    runner = CliRunner()
    res_cov = runner.invoke(app, [
        "coverage",
        "-p", str(primers_path),
        "-f", str(fasta_2strains),
        "--panel", "3",
        "-o", str(out_cov_html),
    ])
    assert res_cov.exit_code == 0, f"Coverage CLI failed: {res_cov.stdout}"

    cov_json = json.loads((tmp_path / "coverage_3n.json").read_text(encoding="utf-8"))
    assert cov_json["total_strains"] == 2
    assert cov_json["indeterminate_strains"] == 1
    assert cov_json["evaluable_strains"] == 1
    assert cov_json["covered_thermo"] == 1
    assert cov_json["thermo_rate"] == 100.0

    strains_by_id = {s["strain_id"]: s for s in cov_json["strains"]}
    s_3n = strains_by_id["M93130_F3_3N"]
    assert s_3n["indeterminate"] is True
    assert s_3n["thermo"] is None
    assert s_3n["evals"]["F3"]["verdict"] == "indetermine"

    # Vérification également en mode analyze sur la souche contenant 3 N (panel 3)
    # Also verify in analyze mode on the strain containing 3 N (panel 3)
    p3_records = [
        f">{h}\n{s}" for h, s in read_multi_fasta(str(primers_path)) if h.startswith("3_")
    ]
    primers_only_3 = tmp_path / "panel_3.txt"
    primers_only_3.write_text("\n".join(p3_records) + "\n", encoding="utf-8")

    fasta_1strain = tmp_path / "strain_3n_only.fasta"
    fasta_1strain.write_text(f">M93130_F3_3N\n{seq_with_3n}\n", encoding="utf-8")
    out_ana_html = tmp_path / "analyze_3n.html"
    res_ana = runner.invoke(app, [
        "analyze",
        "-p", str(primers_only_3),
        "-t", str(fasta_1strain),
        "-o", str(out_ana_html),
    ])
    assert res_ana.exit_code == 0, f"Analyze CLI failed: {res_ana.stdout}"
    ana_html_text = out_ana_html.read_text(encoding="utf-8")
    assert "Indéterminé" in ana_html_text


def test_coverage_mg_dependence_and_m93130_verdicts(tmp_path: Path):
    """Vérifie que :
    1. La couverture à 0, 8 et 20 mM de Mg donne des ΔG absolus strictement différents.
    2. Les verdicts et ΔΔG de M93130.1 sont inchangés, avec B3 toujours tolérable à +1,33.
    Verify that:
    1. Coverage at 0, 8, and 20 mM Mg gives strictly different absolute dG values.
    2. M93130.1 verdicts and ddG are unchanged, with B3 still tolerable at +1.33.
    """
    fixture_dir = Path(__file__).resolve().parent.parent / "fixtures" / "real_primer_files"
    primers_path = fixture_dir / "primer_dengue_3.txt"
    genome_path = fixture_dir / "DEN3_M93130.fasta"

    runner = CliRunner()
    evals_by_mg = {}
    for mg_val in (0.0, 8.0, 20.0):
        out_html = tmp_path / f"cov_mg_{int(mg_val)}.html"
        dntp_val = "0.0" if mg_val == 0.0 else "1.4"
        res = runner.invoke(app, [
            "coverage",
            "-p", str(primers_path),
            "-f", str(genome_path),
            "--panel", "3",
            "-o", str(out_html),
            "--temperature", "63.0",
            "--na", "50.0",
            "--mg", str(mg_val),
            "--dntp", dntp_val,
        ])
        assert res.exit_code == 0, f"Coverage CLI failed for Mg={mg_val}: {res.stdout}"
        data = json.loads((tmp_path / f"cov_mg_{int(mg_val)}.json").read_text(encoding="utf-8"))
        evals_by_mg[mg_val] = data["strains"][0]["evals"]

    for role in ("F3", "B3", "FIP", "BIP", "LF", "LB"):
        dg_0 = round(evals_by_mg[0.0][role]["dg"], 2)
        dg_8 = round(evals_by_mg[8.0][role]["dg"], 2)
        dg_20 = round(evals_by_mg[20.0][role]["dg"], 2)
        # Les ΔG absolus doivent être strictement différents et de plus en plus négatifs avec le Mg
        # Absolute dG values must be strictly different and increasingly negative with Mg
        assert dg_20 < dg_8 < dg_0, (
            f"ΔG non strictement décroissants en Mg pour {role}: 0 mM={dg_0}, 8 mM={dg_8}, 20 mM={dg_20}"
        )
        # Les verdicts et ΔΔG de M93130.1 restent inchangés quelle que soit la concentration en Mg
        # M93130.1 verdicts and ddG remain unchanged across Mg concentrations
        assert evals_by_mg[0.0][role]["verdict"] == evals_by_mg[8.0][role]["verdict"] == evals_by_mg[20.0][role]["verdict"]
        assert round(evals_by_mg[0.0][role]["ddg"], 2) == round(evals_by_mg[8.0][role]["ddg"], 2) == round(evals_by_mg[20.0][role]["ddg"], 2)

    # B3 toujours tolérable à +1,33 kcal/mol
    # B3 always tolerable at +1.33 kcal/mol
    for mg_val in (0.0, 8.0, 20.0):
        assert evals_by_mg[mg_val]["B3"]["verdict"] == "tolerable"
        assert round(evals_by_mg[mg_val]["B3"]["ddg"], 2) == 1.33


