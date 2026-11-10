import pytest
from labcraft.cli.config import PanelConfig, build_engine_from_config
from labcraft.diagnostics.enzyme import get_enzyme
import yaml

def test_enzyme_tolerance():
    assert get_enzyme("bst 2.0").name == "Bst 2.0"
    assert get_enzyme("BST2.0").name == "Bst 2.0"
    assert get_enzyme("Bst2.0_WarmStart").name == "Bst 2.0 WarmStart"
    assert get_enzyme("bst 2.0 warmstart").name == "Bst 2.0 WarmStart"
    assert get_enzyme("taq").name == "Taq"

    with pytest.raises(ValueError):
        get_enzyme("unknown_enzyme")

def test_config_buffer_parsing():
    yaml_str = """
experiment:
  name: "Test"
  chemistry: LAMP
  temperature_C: 65.0
  enzyme: bst2.0
  buffer:
    na_mM: 50.0
    mg_mM: 8.0
targets: []
primer_sets: []
    """
    data = yaml.safe_load(yaml_str)
    config = PanelConfig.model_validate(data)
    
    assert config.experiment.buffer.na_mM == 50.0
    assert config.experiment.buffer.mg_mM == 8.0
    
    _, _, backend, backend_kwargs, _, _, _, _ = build_engine_from_config(config, {})
    assert backend_kwargs['mg_mm'] == 8.0


def test_part_g_default_buffer_and_enzyme_coherence_via_cli(tmp_path):
    """Partie G : Vérifie via CliRunner que :
    1. Sans options de tampon (--na, --mg, --dntp), `analyze` et `coverage` utilisent et affichent
       50.0 mM Na, 8.0 mM Mg, 1.4 mM dNTP.
    2. Le profil `gspssd` est disponible et reprend les paramètres de Bst 2.0.
    3. L'option `--enzyme taq` est bien transmise jusqu'au calcul sur cible et aux dimères.
    """
    from typer.testing import CliRunner
    from labcraft.cli.main import app

    gspssd = get_enzyme("gspssd")
    bst2 = get_enzyme("bst2.0")
    assert gspssd.name == "GspSSD"
    assert gspssd.dimer_dg_threshold == bst2.dimer_dg_threshold
    assert gspssd.three_prime_window == bst2.three_prime_window
    assert gspssd.three_prime_absolute_window == bst2.three_prime_absolute_window

    f3_seq = "CTCGTGTAGGAATGGGAG"
    b3_seq = "CTGGCTTGAGACATCTTCT"
    fip_seq = "GGGAAGACGGTGTGGTTCTGAAGCCTTTCCCCAGA"
    bip_seq = "TCATTTCCAGCTTTGATGCACCACAAAGTCCCAATCA"
    primers_file = tmp_path / "primers.txt"
    primers_file.write_text(
        f">P_F3\n{f3_seq}\n>P_B3\n{b3_seq}\n>P_FIP\n{fip_seq}\n>P_BIP\n{bip_seq}\n"
    )
    target_file = tmp_path / "target.fasta"
    target_file.write_text(f">T1\n{'A'*20}{f3_seq}{'C'*20}CTGAAGCCTTTCCCCAGA{'G'*20}TGATTGGGACTTTGTGGT{'T'*20}AGAAGATGTCTCAAGCCAG{'A'*20}\n")

    out_an = tmp_path / "an.html"
    out_cov = tmp_path / "cov.html"

    runner = CliRunner()
    res_an = runner.invoke(app, ["analyze", "-p", str(primers_file), "-t", str(target_file), "-o", str(out_an)])
    assert res_an.exit_code == 0, f"Analyze failed: {res_an.stdout}"
    assert "Na=50.0 mM, Mg=8.0 mM, dNTP=1.4 mM" in res_an.stdout

    res_cov = runner.invoke(app, ["coverage", "-p", str(primers_file), "-f", str(target_file), "-o", str(out_cov)])
    assert res_cov.exit_code == 0, f"Coverage failed: {res_cov.stdout}"
    assert "Na=50.0 mM, Mg=8.0 mM, dNTP=1.4 mM" in res_cov.stdout

    # Test avec --enzyme taq
    out_taq = tmp_path / "an_taq.html"
    res_taq = runner.invoke(app, [
        "analyze", "-p", str(primers_file), "-t", str(target_file), "-o", str(out_taq), "--enzyme", "taq"
    ])
    assert res_taq.exit_code == 0, f"Analyze taq failed: {res_taq.stdout}"
    assert "Enzyme=Taq" in res_taq.stdout
    html_taq = out_taq.read_text()
    assert "Taq" in html_taq


def test_part_i_missing_files_clean_error_exit_code_1_via_cli(tmp_path):
    """Partie I.3 : Vérifie qu'un fichier absent passé à `analyze` ou `coverage`
    affiche un message d'erreur propre et quitte avec le code 1 sans traceback.
    Part I.3: Verify that a missing file passed to `analyze` or `coverage`
    prints a clean error message and exits with code 1 without traceback.
    """
    from typer.testing import CliRunner
    from labcraft.cli.main import app

    runner = CliRunner()
    missing = str(tmp_path / "introuvable.fasta")
    existing = tmp_path / "existant.fasta"
    existing.write_text(">T1\nATGC\n")

    for args in [
        ["analyze", "-c", missing],
        ["analyze", "-p", missing],
        ["analyze", "-p", str(existing), "-t", missing],
        ["coverage", "-p", missing, "-f", str(existing)],
        ["coverage", "-p", str(existing), "-f", missing],
        ["coverage", "-p", str(existing), "-f", str(existing), "-s", missing],
    ]:
        res = runner.invoke(app, args)
        assert res.exit_code == 1, f"Attendu exit_code=1 pour {args}, obtenu {res.exit_code}: {res.stdout}"
        assert "Erreur de fichier introuvable" in res.stdout
        assert "Traceback" not in res.stdout



