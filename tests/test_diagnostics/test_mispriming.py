import pytest
from labcraft.diagnostics.mispriming import detect_inter_target_mispriming, _revcomp
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.thermo.salt import UnifiedSaltModel
from labcraft.diagnostics.enzyme import ENZYME_REGISTRY

@pytest.fixture
def backend():
    return ViennaSaltShiftBackend(UnifiedSaltModel())

@pytest.fixture
def enzyme():
    return ENZYME_REGISTRY["bst2.0"]

def test_mispriming_perfect_match(backend, enzyme):
    # L'amorce s'apparie parfaitement sur la cible YFV
    p = PhysicalPrimer.from_simple("F3_X", "ATCGATCGATCGATCG", PrimerRole.F3)
    target_seq = "NNNN" + _revcomp("ATCGATCGATCGATCG") + "NNNN"
    
    primer_to_panel = {"F3_X": "X"}
    targets = {"YFV": target_seq}
    
    risks = detect_inter_target_mispriming(
        [p], primer_to_panel, targets, backend, enzyme, temp_celsius=65.0
    )
    
    assert len(risks) > 0
    assert risks[0].primer_name == "F3_X"
    assert risks[0].target_id == "YFV"

def test_mispriming_internal_mismatch(backend, enzyme):
    # L'extrémité 3' (5 derniers nt) est parfaite.
    # On introduit un mismatch interne à la position 7 (en partant du 3').
    seq = "ATCGTTAGCCACAGTA"
    p = PhysicalPrimer.from_simple("F3_X", seq, PrimerRole.F3)
    rc_seq = _revcomp(seq)
    # L'index du 7ème nucléotide depuis le 3' dans rc_seq est 6 (0-indexed).
    target_rc = rc_seq[:6] + ("A" if rc_seq[6] != "A" else "C") + rc_seq[7:]
    
    target_seq = "NNNN" + target_rc + "NNNN"
    
    primer_to_panel = {"F3_X": "X"}
    targets = {"YFV": target_seq}
    
    risks = detect_inter_target_mispriming(
        [p], primer_to_panel, targets, backend, enzyme, temp_celsius=65.0
    )
    
    # Doit être détecté car l'ancrage 3' est conservé et le mismatch interne est toléré
    assert len(risks) > 0

def test_mispriming_terminal_mismatch(backend, enzyme):
    # L'extrémité 3' (position 1 terminale) a un mismatch.
    # L'ancrage sur l'extrémité 3' ne trouvera pas la cible.
    # Ce cas n'est PAS remonté (veto/ancrage).
    seq = "ATCGTTAGCCACAGTA"
    p = PhysicalPrimer.from_simple("F3_X", seq, PrimerRole.F3)
    
    rc_seq = _revcomp(seq)
    # Remplacer le premier nucléotide de rc_seq (complémentaire du 3' terminal de seq)
    target_rc = ("A" if rc_seq[0] != "A" else "C") + rc_seq[1:]
    
    target_seq = "NNNN" + target_rc + "NNNN"
    
    primer_to_panel = {"F3_X": "X"}
    targets = {"YFV": target_seq}
    
    risks = detect_inter_target_mispriming(
        [p], primer_to_panel, targets, backend, enzyme, temp_celsius=65.0
    )
    
    # Ne doit pas être détecté (veto / pas ancré)
    assert len(risks) == 0

def test_mispriming_no_homology(backend, enzyme):
    p = PhysicalPrimer.from_simple("F3_X", "ATCGATCGATCGATCG", PrimerRole.F3)
    target_seq = "AAAAAAAAAAAAAAAAAAAAAAA"
    
    primer_to_panel = {"F3_X": "X"}
    targets = {"YFV": target_seq}
    
    risks = detect_inter_target_mispriming(
        [p], primer_to_panel, targets, backend, enzyme, temp_celsius=65.0
    )
    
    assert len(risks) == 0

def test_mispriming_case_insensitivity(backend, enzyme):
    p = PhysicalPrimer.from_simple("F3_X", "ATCGATCGATCGATCG", PrimerRole.F3)
    target_rc = _revcomp("ATCGATCGATCGATCG")
    
    # Génome majuscules
    target_upper = "NNNN" + target_rc + "NNNN"
    risks_upper = detect_inter_target_mispriming(
        [p], {"F3_X": "X"}, {"YFV": target_upper}, backend, enzyme, temp_celsius=65.0
    )
    
    # Génome minuscules
    target_lower = target_upper.lower()
    risks_lower = detect_inter_target_mispriming(
        [p], {"F3_X": "X"}, {"YFV": target_lower}, backend, enzyme, temp_celsius=65.0
    )
    
    assert len(risks_upper) > 0
    assert len(risks_upper) == len(risks_lower)


def test_part_e_self_target_mispriming_excludes_only_legitimate_locus_via_cli(tmp_path):
    """Partie E : Vérifie via CliRunner sur `labcraft analyze` avec une seule cible où le site de F3
    est recopié une seconde fois ailleurs que :
    1. Le scan de mésamorçage est bien lancé même avec une seule cible.
    2. Le second site copié est signalé comme mésamorçage, mais le site légitime ne l'est pas
       (une seule alerte pour F3 et non deux).
    """
    from typer.testing import CliRunner
    from labcraft.cli.main import app

    f3_seq = "CTCGTGTAGGAATGGGAG"
    b3_seq = "CTGGCTTGAGACATCTTCT"
    # Cible unique contenant le site F3 légitime à l'indice 30, et une copie ectopique à l'indice 150
    target_seq = (
        "A" * 30
        + f3_seq
        + "C" * 50
        + _revcomp(b3_seq)
        + "T" * 33
        + f3_seq
        + "G" * 30
    )

    primers_file = tmp_path / "primers.txt"
    primers_file.write_text(f">P_F3\n{f3_seq}\n>P_B3\n{b3_seq}\n")

    target_file = tmp_path / "single_target.fasta"
    target_file.write_text(f">SynthSingle\n{target_seq}\n")

    out_html = tmp_path / "report_single.html"
    runner = CliRunner()
    res = runner.invoke(app, [
        "analyze", "-p", str(primers_file), "-t", str(target_file), "-o", str(out_html),
        "--temperature", "63.0", "--na", "50.0", "--mg", "8.0", "--dntp", "1.4"
    ])
    assert res.exit_code == 0, f"CLI failed: {res.stdout}"
    html = out_html.read_text()

    # Le mésamorçage doit être signalé dans le rapport HTML même avec une seule cible
    assert "MÉSAMORÇAGE" in html
    assert html.count("F3_P sur la cible SynthSingle") == 1, (
        "Seul le second site copié doit être signalé, le site légitime doit être exclu !"
    )

