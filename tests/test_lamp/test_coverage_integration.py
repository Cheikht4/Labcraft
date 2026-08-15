import pytest
import os
import tempfile
import csv
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.lamp.coverage import CoverageAnalyzer, SiteVerdict
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.thermo.salt import UnifiedSaltModel
from labcraft.diagnostics.enzyme import get_enzyme

def test_full_coverage_suite():
    fasta_dict = {
        "Strain_Perf": "GCCACCTTAAGCCACAGTA" + "A"*50,
        "Strain_Tol":  "GCCACCTTAAGCGACAGTA" + "A"*50,
        "Strain_Veto": "GCCACCTTAAGCCACAGTT" + "A"*50,
        "Strain_Drop": "GCCACATTAAGACAAAGTA" + "A"*50,
        "Strain_Loop": "GCCACCTTAAGCCACAGTA" + "A"*50
    }
    
    primers = [
        PhysicalPrimer('F3', 'GCCACCTTAAGCCACAGTA', PrimerRole.F3, 'GCCACCTTAAGCCACAGTA'),
        PhysicalPrimer('B3', 'GTTGTGTCATGGGAGGG', PrimerRole.B3, 'GTTGTGTCATGGGAGGG'),
        PhysicalPrimer('FIP', 'TGGCTTTTGGGCCTGACTTCTTTTTTGAAGAAGCTGTGCAGCCTG', PrimerRole.FIP, 'GAAGAAGCTGTGCAGCCTG', 'TGGCTTTTGGGCCTGACTTC'),
        PhysicalPrimer('BIP', 'CTGTAGCTCCGTCGTGGGGATTTTCTAGTCTGCTACACCGTGC', PrimerRole.BIP, 'CTAGTCTGCTACACCGTGC', 'CTGTAGCTCCGTCGTGGGGA'),
        PhysicalPrimer('LF', 'CCTTGGACGGGGCT', PrimerRole.LF, 'CCTTGGACGGGGCT'),
    ]
    
    csv_records = []
    for s_id in fasta_dict.keys():
        csv_records.append({"strain_id": s_id, "primer_name": "F3", "primer_role": "F3", "position": 0, "strand": "+", "site_seq": ""})
        
        # Perfect matches for the others
        csv_records.append({"strain_id": s_id, "primer_name": "B3", "position": 0, "strand": "-", "site_seq": "CCCTCCCATGACACAAC"})
        csv_records.append({"strain_id": s_id, "primer_name": "FIP", "position": 0, "strand": "+", "site_seq": "GAAGAAGCTGTGCAGCCTG"})
        csv_records.append({"strain_id": s_id, "primer_name": "BIP", "position": 0, "strand": "+", "site_seq": "CTAGTCTGCTACACCGTGC"})
        
        if s_id == "Strain_Loop":
            csv_records.append({"strain_id": s_id, "primer_name": "LF", "position": 0, "strand": "+", "site_seq": "CCTTGGACGGGGAA"})
        else:
            csv_records.append({"strain_id": s_id, "primer_name": "LF", "position": 0, "strand": "+", "site_seq": "CCTTGGACGGGGCT"})

    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    enzyme = get_enzyme("Bst2.0")
    analyzer = CoverageAnalyzer(primers, fasta_dict, backend, enzyme, temp_celsius=63.0, dg_threshold=-3.0, max_mismatches_count=2)
    
    verdicts = {v.strain_id: v for v in analyzer.analyze_strains(csv_records)}
    
    assert verdicts["Strain_Perf"].is_amplifiable_thermo == True
    assert verdicts["Strain_Perf"].is_amplifiable_count == True
    assert verdicts["Strain_Perf"].evaluations["F3"].verdict == SiteVerdict.PARFAIT
    
    assert verdicts["Strain_Tol"].is_amplifiable_thermo == True
    assert verdicts["Strain_Tol"].is_amplifiable_count == True
    assert verdicts["Strain_Tol"].evaluations["F3"].verdict == SiteVerdict.TOLERABLE
    assert verdicts["Strain_Tol"].evaluations["F3"].n_mismatches_count == 1
    
    assert verdicts["Strain_Veto"].is_amplifiable_thermo == False
    assert verdicts["Strain_Veto"].is_amplifiable_count == True
    assert verdicts["Strain_Veto"].evaluations["F3"].verdict == SiteVerdict.VETO_3P
    
    assert verdicts["Strain_Drop"].is_amplifiable_thermo == False
    assert verdicts["Strain_Drop"].is_amplifiable_count == False
    assert verdicts["Strain_Drop"].evaluations["F3"].verdict == SiteVerdict.ABSENT
    
    assert verdicts["Strain_Loop"].is_amplifiable_thermo == True
    assert verdicts["Strain_Loop"].evaluations["LF"].verdict == SiteVerdict.VETO_3P

def test_performance():
    pass
