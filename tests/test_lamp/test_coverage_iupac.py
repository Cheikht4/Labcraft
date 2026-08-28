import pytest
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole, _revcomp
from labcraft.lamp.coverage import CoverageAnalyzer, SiteVerdict
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.thermo.salt import UnifiedSaltModel
from labcraft.diagnostics.enzyme import get_enzyme


def test_coverage_iupac_perfect_match():
    """Vérifie qu'une amorce contenant des codes IUPAC (W, Y, R) appariant parfaitement
    la cible produit 0 mésappariement et un verdict PARFAIT sans pénalité indue.
    """
    # 3_F3 de primer_dengue_3.txt : CTCGTGTWGGAATGGGAG (W = A ou T)
    primer_f3 = PhysicalPrimer('3_F3', 'CTCGTGTWGGAATGGGAG', PrimerRole.F3, 'CTCGTGTWGGAATGGGAG')
    # 3_B3 : CTGGYTTGAGACATCTTCT (Y = C ou T)
    primer_b3 = PhysicalPrimer('3_B3', 'CTGGYTTGAGACATCTTCT', PrimerRole.B3, 'CTGGYTTGAGACATCTTCT')
    # 3_F2 : CTGAWGCCTTTCCYCAGA (W=A/T, Y=C/T)
    primer_fip = PhysicalPrimer('3_FIP', 'GGGAARACGGTGTGGTTCTGAWGCCTTTCCYCAGA', PrimerRole.FIP, 'CTGAWGCCTTTCCYCAGA', 'GGGAARACGGTGTGGTT')
    # 3_B2 : ACCACRAAGTCCCAATCA (R=A/G)
    primer_bip = PhysicalPrimer('3_BIP', 'TCATTTCCRGCTTTGATGCACCACRAAGTCCCAATCA', PrimerRole.BIP, 'ACCACRAAGTCCCAATCA', 'TCATTTCCRGCTTTGATGC')

    primers = [primer_f3, primer_b3, primer_fip, primer_bip]

    # Génome avec les bases correspondantes sur DEN-3
    # F3 site : CTCGTGTAGGAATGGGAG (le W apparie A)
    site_f3 = "CTCGTGTAGGAATGGGAG"
    # B3 est sur brin -, donc site extrait sur brin + est le reverse complement de B3 (len 19)
    site_b3 = _revcomp("CTGGCTTGAGACATCTTCT")
    # FIP binding (F2) sur brin +
    site_f2 = "CTGATGCCTTTCCCCAGA"
    # BIP binding (B2) sur brin -
    site_b2 = _revcomp("ACCACAAAGTCCCAATCA")

    fasta_dict = {
        "DEN3_Ref": site_f3 + "A"*50 + site_f2 + "A"*50
    }

    csv_records = [
        {"strain_id": "DEN3_Ref", "primer_role": "F3", "position": 0, "strand": "+", "site_seq": site_f3, "n_mismatches": 0},
        {"strain_id": "DEN3_Ref", "primer_role": "B3", "position": 0, "strand": "-", "site_seq": site_b3, "n_mismatches": 0},
        {"strain_id": "DEN3_Ref", "primer_role": "FIP", "position": 0, "strand": "+", "site_seq": site_f2, "n_mismatches": 0},
        {"strain_id": "DEN3_Ref", "primer_role": "BIP", "position": 0, "strand": "-", "site_seq": site_b2, "n_mismatches": 0},
    ]

    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    enzyme = get_enzyme("Bst2.0")
    analyzer = CoverageAnalyzer(primers, fasta_dict, backend, enzyme, temp_celsius=63.0, dg_threshold=-6.0, max_mismatches_count=2)

    verdicts = analyzer.analyze_strains(csv_records)
    assert len(verdicts) == 1
    v = verdicts[0]

    assert v.evaluations["F3"].n_mismatches_count == 0
    assert v.evaluations["F3"].verdict == SiteVerdict.PARFAIT
    assert v.evaluations["B3"].n_mismatches_count == 0
    assert v.evaluations["B3"].verdict == SiteVerdict.PARFAIT
    assert v.evaluations["FIP"].n_mismatches_count == 0
    assert v.evaluations["FIP"].verdict == SiteVerdict.PARFAIT
    assert v.evaluations["BIP"].n_mismatches_count == 0
    assert v.evaluations["BIP"].verdict == SiteVerdict.PARFAIT
    assert v.is_amplifiable_thermo is True
    assert v.is_amplifiable_count is True


def test_coverage_variants_selection():
    """Vérifie que la meilleure variante est retenue pour un rôle dégénéré."""
    # Variante 1 : veto en 3'
    p1 = PhysicalPrimer('F3_1', 'GCCACCTTAAGCCACAGTT', PrimerRole.F3, 'GCCACCTTAAGCCACAGTT')
    # Variante 2 : parfaite
    p2 = PhysicalPrimer('F3_2', 'GCCACCTTAAGCCACAGTA', PrimerRole.F3, 'GCCACCTTAAGCCACAGTA')

    primers = [
        p1, p2,
        PhysicalPrimer('B3', 'GTTGTGTCATGGGAGGG', PrimerRole.B3, 'GTTGTGTCATGGGAGGG'),
        PhysicalPrimer('FIP', 'TGGCTTTTGGGCCTGACTTCTTTTTTGAAGAAGCTGTGCAGCCTG', PrimerRole.FIP, 'GAAGAAGCTGTGCAGCCTG', 'TGGCTTTTGGGCCTGACTTC'),
        PhysicalPrimer('BIP', 'CTGTAGCTCCGTCGTGGGGATTTTCTAGTCTGCTACACCGTGC', PrimerRole.BIP, 'CTAGTCTGCTACACCGTGC', 'CTGTAGCTCCGTCGTGGGGA'),
    ]

    fasta_dict = {"Strain_1": "GCCACCTTAAGCCACAGTA" + "A"*50}

    csv_records = [
        {"strain_id": "Strain_1", "primer_role": "F3", "primer_name": "F3_1", "position": 0, "strand": "+", "site_seq": "GCCACCTTAAGCCACAGTA", "n_mismatches": 1},
        {"strain_id": "Strain_1", "primer_role": "F3", "primer_name": "F3_2", "position": 0, "strand": "+", "site_seq": "GCCACCTTAAGCCACAGTA", "n_mismatches": 0},
        {"strain_id": "Strain_1", "primer_role": "B3", "position": 0, "strand": "-", "site_seq": "CCCTCCCATGACACAAC", "n_mismatches": 0},
        {"strain_id": "Strain_1", "primer_role": "FIP", "position": 0, "strand": "+", "site_seq": "GAAGAAGCTGTGCAGCCTG", "n_mismatches": 0},
        {"strain_id": "Strain_1", "primer_role": "BIP", "position": 0, "strand": "+", "site_seq": "CTAGTCTGCTACACCGTGC", "n_mismatches": 0},
    ]

    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    enzyme = get_enzyme("Bst2.0")
    analyzer = CoverageAnalyzer(primers, fasta_dict, backend, enzyme, temp_celsius=63.0, dg_threshold=-6.0)

    verdicts = analyzer.analyze_strains(csv_records)
    assert len(verdicts) == 1
    v = verdicts[0]

    # La variante 2 (parfaite) doit l'emporter sur la variante 1 (veto 3')
    assert v.evaluations["F3"].verdict == SiteVerdict.PARFAIT
    assert v.evaluations["F3"].primer_name == "F3_2"
    assert v.is_amplifiable_thermo is True
