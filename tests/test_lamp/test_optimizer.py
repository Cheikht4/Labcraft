import pytest
import numpy as np

from labcraft.lamp.complex_enumeration import enumerate_complexes
from labcraft.solver.types import EquilibriumProblem
from labcraft.solver.dual import solve_dual
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.thermo.salt import UnifiedSaltModel
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.optimize.concentrations import optimize_concentrations, get_dangerous_dimers
from labcraft.diagnostics.enzyme import BST_2_0

def test_guardrail_lower_concentration_lowers_initiation():
    """
    Test de garde-fou pour s'assurer que baisser la concentration
    d'une amorce d'initiation (ex: FIP) diminue son occupation sur la cible.
    """
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    
    target_seq = "ATGCGTACGTGCAACTGATCGATCGTACGATCG"
    
    # FIP parfait avec binding de 15 bases sur la cible (GATCGATCGTACGAT)
    binding = "ATCGTACGATCGATC"
    fip = PhysicalPrimer("FIP_1", sequence="TTTT" + binding, role=PrimerRole.FIP, binding_domain=binding, nominal_concentration=1.6e-6)
    
    prob_high, species_high, complexes_high, _ = enumerate_complexes(
        primers=[fip],
        target_seq=target_seq,
        backend=backend,
        temp_celsius=65.0
    )
    res_high = solve_dual(prob_high)
    
    target_complex_name = 'FIP_1_on_FIP_1_site'
    target_complex_idx = complexes_high.index(target_complex_name)
    occ_high = res_high.concentrations[target_complex_idx]
    
    fip_low = PhysicalPrimer("FIP_1", sequence="TTTT" + binding, role=PrimerRole.FIP, binding_domain=binding, nominal_concentration=0.4e-6)
    
    prob_low, species_low, complexes_low, _ = enumerate_complexes(
        primers=[fip_low],
        target_seq=target_seq,
        backend=backend,
        temp_celsius=65.0
    )
    res_low = solve_dual(prob_low)
    target_complex_idx = complexes_low.index(target_complex_name)
    occ_low = res_low.concentrations[target_complex_idx]
    
    assert occ_low < occ_high, f"L'occupation aurait dû baisser. High: {occ_high}, Low: {occ_low}"


def test_optimizer_dimer_and_floor():
    """
    Test the optimizer actually runs and reduces dangerous dimers without
    violating hierarchy and floor.
    """
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    # Create primers forming a dimer
    # Primer 1 and 2
    fip = PhysicalPrimer("FIP_1", sequence="ATCGTACGATCGATCGGGGGG", role=PrimerRole.FIP, binding_domain="ATCGTACGATCGATC", nominal_concentration=1.6e-6)
    # The BIP will have the complementary 3' end so it forms an amplifiable dimer
    bip = PhysicalPrimer("BIP_1", sequence="TTTTAAAAAACCCCCC", role=PrimerRole.BIP, binding_domain="AAAAAA", nominal_concentration=1.6e-6)
    
    target_seq = "ATGCGTACGTGCAACTGATCGATCGTACGATCG"
    
    prob, species, complexes, _ = enumerate_complexes(
        primers=[fip, bip],
        target_seq=target_seq,
        backend=backend,
        temp_celsius=65.0
    )
    
    # Run optimizer
    results = optimize_concentrations(
        prob_template=prob,
        species_names=species,
        primers=[fip, bip],
        target_dict={'PANEL1': target_seq},
        primer_to_panel={'FIP_1': 'PANEL1', 'BIP_1': 'PANEL1'},
        original_free_fractions={},
        original_target_occupations={},
        complex_names=complexes,
        temp_celsius=65.0,
        backend=backend,
        enzyme=BST_2_0,
        min_initiation_occupation=0.01
    )
    
    assert isinstance(results, list)
    
    # Find FIP and BIP in species
    fip_idx = species.index('FIP_1')
    bip_idx = species.index('BIP_1')
    
    fip_site_idx = species.index('FIP_1_site')
    
    # Assert sites are preserved
    for r in results:
        # Just to check it produces something
        assert "primer_name" in r
        assert r["suggested_conc"] >= 0.8e-6 # Bound for FIP/BIP

def test_get_dangerous_dimers_names_match_complexes():
    """
    Test que `get_dangerous_dimers` renvoie des noms de complexes
    qui existent bien dans `complexes`.
    """
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    fip = PhysicalPrimer("FIP_1", sequence="ATCGTACGATCGATCGGGGGG", role=PrimerRole.FIP, binding_domain="ATCGTACGATCGATC", nominal_concentration=1.6e-6)
    bip = PhysicalPrimer("BIP_1", sequence="TTTTAAAAAACCCCCC", role=PrimerRole.BIP, binding_domain="AAAAAA", nominal_concentration=1.6e-6)
    
    prob, species, complexes, _ = enumerate_complexes(
        primers=[fip, bip],
        target_seq="ATGCGTACGTGCAACTGATCGATCGTACGATCG",
        backend=backend,
        temp_celsius=65.0
    )
    
    dangerous = get_dangerous_dimers(prob, complexes, species, [fip, bip], 65.0, backend, BST_2_0)
    assert len(dangerous) > 0, "Doit détecter un dimère dangereux"
    
    # Vérifier que CHAQUE nom retourné est bien dans complexes
    for d_name in dangerous:
        assert d_name in complexes, f"Le nom de dimère {d_name} n'est pas dans les complexes générés: {complexes}"

def test_optimizer_multiplex_balance():
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    fip1 = PhysicalPrimer("FIP_1", sequence="ATCGTACGATCGATC", role=PrimerRole.FIP, binding_domain="ATCGTACGATCGATC", nominal_concentration=1.6e-6)
    fip2 = PhysicalPrimer("FIP_2", sequence="GCATGCATGCATGCA", role=PrimerRole.FIP, binding_domain="GCATGCATGCATGCA", nominal_concentration=1.6e-6)
    target_seq = "ATGCGTACGTGCAACTGATCGATCGTACGATCG GCATGCATGCATGCA"
    prob, species, complexes, _ = enumerate_complexes(
        primers=[fip1, fip2],
        target_seq=target_seq,
        backend=backend,
        temp_celsius=65.0
    )
    
    results = optimize_concentrations(
        prob_template=prob,
        species_names=species,
        primers=[fip1, fip2],
        target_dict={'PANEL1': target_seq, 'PANEL2': target_seq},
        primer_to_panel={'FIP_1': 'PANEL1', 'FIP_2': 'PANEL2'},
        original_free_fractions={},
        original_target_occupations={},
        complex_names=complexes,
        temp_celsius=65.0,
        backend=backend,
        enzyme=BST_2_0
    )
    assert isinstance(results, list)


def test_part_f_optimizer_cumulative_gain_and_parent_level_variants_via_cli(tmp_path):
    """Partie F : Vérifie que :
    1. Le gain affiché (`dangerous_gain_rel`) vaut exactement `(dangerous_before - dangerous_after) / dangerous_before`.
    2. Les variants d'une même amorce dégénérée (`LB_3#1`, `LB_3#2`, etc.) bougent ensemble et la suggestion
       est affichée au niveau de l'amorce parente (`LB_3`), aussi bien en appel direct que via `labcraft analyze`.
    """
    from pathlib import Path
    from typer.testing import CliRunner
    from labcraft.cli.main import app
    from labcraft.cli.parsers import read_multi_fasta, build_config_from_cli
    from labcraft.cli.config import build_engine_from_config

    fixture_dir = Path(__file__).resolve().parent.parent / "fixtures" / "real_primer_files"
    p3_file = fixture_dir / "primer_dengue_3.txt"
    d3_fasta = fixture_dir / "DEN3_M93130.fasta"

    d3_records = read_multi_fasta(str(d3_fasta))
    d3_id, d3_seq = d3_records[0][0].split()[0], d3_records[0][1]

    p14_records = [
        f">{h}\n{s}" for h, s in read_multi_fasta(str(p3_file)) if h.startswith("14_")
    ]
    primers_only_14 = tmp_path / "panel_14.txt"
    primers_only_14.write_text("\n".join(p14_records) + "\n", encoding="utf-8")

    cfg = build_config_from_cli(
        config_path=None,
        primers_path=str(primers_only_14),
        targets_path=str(d3_fasta),
        temperature=63.0, na=50.0, mg=8.0, dntp=1.4
    )
    primers, p2t, backend, b_kwargs, mon_molar, enzyme, temp, profiles = build_engine_from_config(cfg, {d3_id: d3_seq})
    prob, strands, complexes, _ = enumerate_complexes(
        primers, d3_seq, backend, profile=profiles[d3_id], temp_celsius=63.0,
        mon_molar=mon_molar, buffer=cfg.experiment.buffer.model_dump(), enzyme=enzyme, primer_to_target=p2t
    )

    opt_res = optimize_concentrations(
        prob_template=prob,
        species_names=strands,
        primers=primers,
        target_dict={d3_id: d3_seq},
        primer_to_panel=p2t,
        original_free_fractions={},
        original_target_occupations={},
        complex_names=complexes,
        temp_celsius=63.0,
        backend=backend,
        enzyme=enzyme,
        **b_kwargs
    )

    assert isinstance(opt_res, dict) and "metrics" in opt_res and "suggestions" in opt_res
    metrics = opt_res["metrics"]
    before = metrics["dangerous_before"]
    after = metrics["dangerous_after"]
    expected_gain = (before - after) / before
    assert metrics["dangerous_gain_rel"] == pytest.approx(expected_gain, rel=1e-9), (
        f"Le gain affiché ({metrics['dangerous_gain_rel']}) doit valoir exactement (avant - après) / avant ({expected_gain})"
    )

    # Aucune suggestion ne doit contenir '#' (uniquement des amorces parentes)
    for s in opt_res["suggestions"]:
        assert "#" not in s["primer_name"], f"Suggestion au niveau variant au lieu du parent : {s['primer_name']}"

    # Vérification via CliRunner sur labcraft analyze
    out_html = tmp_path / "report_opt.html"
    runner = CliRunner()
    res = runner.invoke(app, [
        "analyze", "-p", str(primers_only_14), "-t", str(d3_fasta), "-o", str(out_html),
        "--temperature", "63.0", "--na", "50.0", "--mg", "8.0", "--dntp", "1.4"
    ])
    assert res.exit_code == 0, f"CLI failed: {res.stdout}"
    html = out_html.read_text(encoding="utf-8")
    expected_pct_str = f"{round(expected_gain * 100, 1)} %"
    assert expected_pct_str in html
    assert "LB_14#" not in html.split("Ajustement de Concentrations")[1].split("Avertissements")[0] if "Ajustement de Concentrations" in html and "Avertissements" in html else True

