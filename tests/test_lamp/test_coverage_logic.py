import pytest
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.lamp.coverage import CoverageAnalyzer, SiteVerdict
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.thermo.salt import UnifiedSaltModel
from labcraft.diagnostics.enzyme import get_enzyme

def test_coverage_analyzer():
    primers = [
        PhysicalPrimer('F3', 'GCCACCTTAAGCCACAGTA', PrimerRole.F3, 'GCCACCTTAAGCCACAGTA'),
        PhysicalPrimer('B3', 'GTTGTGTCATGGGAGGG', PrimerRole.B3, 'GTTGTGTCATGGGAGGG'),
        PhysicalPrimer('FIP', 'TGGCTTTTGGGCCTGACTTCTTTTTTGAAGAAGCTGTGCAGCCTG', PrimerRole.FIP, 'GAAGAAGCTGTGCAGCCTG', 'TGGCTTTTGGGCCTGACTTC'),
        PhysicalPrimer('BIP', 'CTGTAGCTCCGTCGTGGGGATTTTCTAGTCTGCTACACCGTGC', PrimerRole.BIP, 'CTAGTCTGCTACACCGTGC', 'CTGTAGCTCCGTCGTGGGGA'),
    ]
    
    fasta_dict = {
        "Strain_Perf": "GCCACCTTAAGCCACAGTA" + "A"*50 + "GTTGTGTCATGGGAGGG" + "C"*50,
        "Strain_Veto": "GCCACCTTAAGCCACAGTT" + "A"*50 + "GTTGTGTCATGGGAGGG" + "C"*50,
    }
    
    csv_records = [
        {"strain_id": "Strain_Perf", "primer_role": "F3", "n_mismatches": 0, "position": 0, "strand": "+", "site_seq": "GCCACCTTAAGCCACAGTA"},
        {"strain_id": "Strain_Veto", "primer_role": "F3", "n_mismatches": 0, "position": 0, "strand": "+", "site_seq": "GCCACCTTAAGCCACAGTT"},
    ]
    
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    enzyme = get_enzyme("Bst2.0")
    
    analyzer = CoverageAnalyzer(primers, fasta_dict, backend, enzyme, dg_threshold=-6.0)
    verdicts = analyzer.analyze_strains(csv_records)
    
    assert len(verdicts) == 2
    
    perf = [v for v in verdicts if v.strain_id == "Strain_Perf"][0]
    veto = [v for v in verdicts if v.strain_id == "Strain_Veto"][0]
    
    assert perf.evaluations["F3"].verdict == SiteVerdict.PARFAIT
    assert veto.evaluations["F3"].verdict == SiteVerdict.VETO_3P
    
    # Init primers missing => False
    assert perf.is_amplifiable_thermo == False # B3, FIP, BIP missing
