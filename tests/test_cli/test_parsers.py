import pytest
import tempfile
from labcraft.cli.parsers import (
    read_multi_fasta, parse_primer_name, parse_primer_file, 
    reconstruct_fip_bip, ParseError
)
from labcraft.cli.config import PrimerSetConfig

def test_read_multi_fasta_standard():
    with tempfile.NamedTemporaryFile("w", delete=False, encoding="utf-8") as f:
        f.write(">Seq1\nATGC\n>Seq2\nCGTA\n")
        path = f.name
        
    records = read_multi_fasta(path)
    assert len(records) == 2
    assert records[0] == ("Seq1", "ATGC")
    assert records[1] == ("Seq2", "CGTA")

def test_read_multi_fasta_oneline():
    with tempfile.NamedTemporaryFile("w", delete=False, encoding="utf-8") as f:
        f.write(">Seq1    ATGC\n>Seq2\tCGTA\n")
        path = f.name
        
    records = read_multi_fasta(path)
    assert len(records) == 2
    assert records[0] == ("Seq1", "ATGC")
    assert records[1] == ("Seq2", "CGTA")

def test_read_multi_fasta_duplicates():
    with tempfile.NamedTemporaryFile("w", delete=False, encoding="utf-8") as f:
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
    """Vérifie qu'une cible manquante lève une erreur explicite (pas de RC silencieux).
    Verify that a missing target raises an explicit error (no silent fallback)."""
    primers = {"F1": "ATGC", "F2": "CGTA"}
    with pytest.raises(ParseError, match=">Panel1_FIP"):
        reconstruct_fip_bip("Panel1", primers, [])

def test_reconstruct_fip_bip_unmatched_f1_or_b1_raises():
    """Vérifie qu'un F1/B1 introuvable sur la cible lève une ParseError explicite.
    Verify that an unlocatable F1/B1 on target raises an explicit ParseError."""
    target = "AAAAACGTAAAAAAAGCATA"
    with pytest.raises(ParseError, match=">Panel1_FIP"):
        reconstruct_fip_bip("Panel1", {"F1": "GGGGGGGGGG", "F2": "CGTA"}, [("Panel1", target)])
    with pytest.raises(ParseError, match=">Panel1_BIP"):
        reconstruct_fip_bip("Panel1", {"B1": "GGGGGGGGGG", "B2": "CGTA"}, [("Panel1", target)])

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

def test_reconstruct_bip_with_target_inverse_of_f():
    """Vérifie que B1c reste sur le brin sens (+) opposé à B2 (-), à l'inverse de F1c (-).
    Verify that B1c stays on sense (+) strand opposite to B2 (-), inverse of F1c (-)."""
    # CTGACTAG sur brin sens (+), son RC est CTAGTCAG (-)
    # ACCAGTAG sur brin sens (+) correspond au site cible de B2 = CTACTGGT (-)
    target = "AAAACTGACTAGTTTTACCAGTAGCCCC"

    # Cas 1 : B1 déjà en B1c (brin sens +) -> conservé tel quel
    primers_b1c = {"B1": "CTGACTAG", "B2": "CTACTGGT"}
    reconstruct_fip_bip("Panel1", primers_b1c, [("Panel1", target)], linker="TTTT")
    assert primers_b1c["BIP"] == "CTGACTAGTTTTCTACTGGT"
    assert "B1" not in primers_b1c

    # Cas 2 : B1 en brin antisens (-) -> complémenté inverse vers +
    primers_b1_minus = {"B1": "CTAGTCAG", "B2": "CTACTGGT"}
    reconstruct_fip_bip("Panel1", primers_b1_minus, [("Panel1", target)], linker="TTTT")
    assert primers_b1_minus["BIP"] == "CTGACTAGTTTTCTACTGGT"

def test_direct_full_fip_bip_precedence():
    """Vérifie que >Panel_FIP et >Panel_BIP complets priment sans reconstruction.
    Verify that full >Panel_FIP and >Panel_BIP take precedence without reconstruction."""
    full_fip = "ATGCCCCCACGT"
    full_bip = "CTGACTAGAAAACTACTGGT"
    primers = {
        "FIP": full_fip,
        "F1": "GGGGGGGGGG",
        "F2": "CGTA",
        "BIP": full_bip,
        "B1": "CCCCCCCCCC",
        "B2": "CTACTGGT",
    }
    # Ne doit pas lever d'erreur même sans cible ni si F1/B1 sont invalides
    domains = reconstruct_fip_bip("Panel1", primers, [], linker="TTTT")
    assert primers["FIP"] == full_fip
    assert primers["BIP"] == full_bip
    assert "F1" not in primers and "F2" not in primers
    assert "B1" not in primers and "B2" not in primers
    assert domains == {}

    # Vérification via parse_primer_file complet avec >Panel_FIP / >Panel_BIP
    with tempfile.NamedTemporaryFile("w", delete=False, encoding="utf-8") as f:
        f.write(f">Panel1_FIP\n{full_fip}\n>Panel1_BIP\n{full_bip}\n>Panel1_F1\nAAAAAAAA\n")
        path = f.name
    sets = parse_primer_file(path, [])
    assert len(sets) == 1
    assert sets[0].primers["FIP"].seq == full_fip
    assert sets[0].primers["BIP"].seq == full_bip

def test_parse_primer_file_versions():
    with tempfile.NamedTemporaryFile("w", delete=False, encoding="utf-8") as f:
        f.write(">Panel_F3_1\nAT\n>Panel_B3_1\nCG\n>Panel_F3_2\nTA\n>Panel_B3_2\nGC\n")
        path = f.name
        
    with pytest.warns(UserWarning, match="Multiple versions"):
        sets = parse_primer_file(path, [("Panel", "ATCGTA")])
        
    assert len(sets) == 2
    assert sets[0].target == "Panel"
    assert sets[1].target == "Panel"
    assert sets[0].primers["F3"].seq == "AT"
    assert sets[1].primers["F3"].seq == "TA"
