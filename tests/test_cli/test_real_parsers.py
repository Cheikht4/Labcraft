import pytest
from pathlib import Path
from labcraft.cli.parsers import parse_primer_file
from labcraft.cli.parsers import read_multi_fasta

FIXTURE_DIR = Path(__file__).parents[1] / "fixtures" / "real_primer_files"

def test_real_primer_dengue_3():
    # 6 oligos physiques, avec target
    primer_file = FIXTURE_DIR / "primer_dengue_3.txt"
    target_fasta = FIXTURE_DIR / "DEN3_M93130.fasta"
        
    target_seq = read_multi_fasta(str(target_fasta))[0][1]
    
    with pytest.warns(UserWarning, match="n'a pas pu être vérifiée"):
        sets = parse_primer_file(str(primer_file), [("DEN3_M93130", target_seq)])
        
    assert len(sets) == 2
    # Check that it produces 6 oligonucleotides
    assert len(sets[0].primers) == 6
    for role in ["F3", "B3", "FIP", "BIP", "LF", "LB"]:
        assert role in sets[0].primers
        assert role in sets[1].primers

def test_real_amorceparida():
    primer_file = FIXTURE_DIR / "amorceparida.txt"
    target_fasta = FIXTURE_DIR / "DEN3_M93130.fasta"
    
    target_seq = read_multi_fasta(str(target_fasta))[0][1]
    
    with pytest.warns(UserWarning, match="n'a pas pu être vérifiée"):
        sets = parse_primer_file(str(primer_file), [("DEN3_M93130", target_seq)])
    assert len(sets) == 4
    for i in range(4):
        assert len(sets[i].primers) == 6

def test_real_nawar_dengue1():
    primer_file = FIXTURE_DIR / "Nawar_dengue1.txt"
        
    # Nawar has two panels: '1' and 'DENV1'. 
    # Target name provided doesn't exactly match, but includes them
    sets = parse_primer_file(str(primer_file), [("DENV1_Genome", "ATGC")])
    assert len(sets) == 2
    assert sets[0].target == "DENV1_Genome"
    assert sets[1].target == "DENV1_Genome"
    assert len(sets[0].primers) == 6
    assert len(sets[1].primers) == 6
    for role in ["F3", "B3", "FIP", "BIP", "LF", "LB"]:
        assert role in sets[0].primers
        assert role in sets[1].primers

def test_type_inconnu():
    with pytest.raises(Exception, match="contains no valid role"):
        from labcraft.cli.parsers import parse_primer_name
        parse_primer_name("Panel1_STEMF")
