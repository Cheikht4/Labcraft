import pytest
import tempfile
from labcraft.cli.parsers import (
    read_multi_fasta, parse_primer_name, parse_primer_file, 
    reconstruct_fip_bip, ParseError
)
from labcraft.cli.config import PrimerSetConfig

def test_read_multi_fasta_standard():
    with tempfile.NamedTemporaryFile("w", delete=False) as f:
        f.write(">Seq1\nATGC\n>Seq2\nCGTA\n")
        path = f.name
        
    records = read_multi_fasta(path)
    assert len(records) == 2
    assert records[0] == ("Seq1", "ATGC")
    assert records[1] == ("Seq2", "CGTA")

def test_read_multi_fasta_oneline():
    with tempfile.NamedTemporaryFile("w", delete=False) as f:
        f.write(">Seq1    ATGC\n>Seq2\tCGTA\n")
        path = f.name
        
    records = read_multi_fasta(path)
    assert len(records) == 2
    assert records[0] == ("Seq1", "ATGC")
    assert records[1] == ("Seq2", "CGTA")

def test_read_multi_fasta_duplicates():
    with tempfile.NamedTemporaryFile("w", delete=False) as f:
        f.write(">Seq1\nATGC\n>Seq1\nCGTA\n")
        path = f.name
        
    with pytest.raises(ParseError, match="Duplicate header"):
        read_multi_fasta(path)

def test_parse_primer_name_aliases():
    assert parse_primer_name("Panel1_FLOOP")[1] == "LF"
    assert parse_primer_name("Panel1_FLP")[1] == "LF"
    assert parse_primer_name("Panel1_BLOOP")[1] == "LB"
    assert parse_primer_name("Panel1_BLP")[1] == "LB"

def test_parse_primer_name_versions():
    panel, role, version = parse_primer_name("DENV3_F3_1")
    assert panel == "DENV3"
    assert role == "F3"
    assert version == "1"
    
    panel, role, version = parse_primer_name("DENV3_F3")
    assert panel == "DENV3"
    assert role == "F3"
    assert version is None

def test_parse_primer_name_unknown():
    with pytest.raises(ParseError, match="contains no valid role"):
        parse_primer_name("Panel1_STEMF")

def test_reconstruct_fip_bip_no_target():
    primers = {"F1": "ATGC", "F2": "CGTA"}
    reconstruct_fip_bip("Panel1", primers, [])
    assert "FIP" in primers
    assert "F1" not in primers
    assert primers["FIP"] == "GCATCGTA"

def test_reconstruct_fip_bip_with_target_f1_rc():
    target = "AAAAACGTAAAAAAAGCATA"
    primers = {"F1": "GCAT", "F2": "CGTA"}
    reconstruct_fip_bip("Panel1", primers, [("Panel1", target)])
    assert primers["FIP"] == "ATGCCGTA"

def test_reconstruct_fip_bip_with_target_f1c():
    target = "AAAAACGTAAAAAAAGCATA"
    primers = {"F1": "ATGC", "F2": "CGTA"}
    reconstruct_fip_bip("Panel1", primers, [("Panel1", target)])
    assert primers["FIP"] == "ATGCCGTA"

def test_parse_primer_file_versions():
    with tempfile.NamedTemporaryFile("w", delete=False) as f:
        f.write(">Panel_F3_1\nAT\n>Panel_B3_1\nCG\n>Panel_F3_2\nTA\n>Panel_B3_2\nGC\n")
        path = f.name
        
    with pytest.warns(UserWarning, match="Multiple versions"):
        sets = parse_primer_file(path, [("Panel", "ATCGTA")])
        
    assert len(sets) == 2
    assert sets[0].target == "Panel"
    assert sets[1].target == "Panel"
    assert sets[0].primers["F3"].seq == "AT"
    assert sets[1].primers["F3"].seq == "TA"
