"""Tests unitaires et d'intégration de l'Étape 3 (`labcraft redesign` V2).
Unit and integration tests for Step 3 (`labcraft redesign` V2).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from labcraft.cli.main import app
from labcraft.cli.parsers import read_multi_fasta
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.optimize.redesign import (
    Candidate,
    Locus,
    compute_inner_comparison,
    compute_max_3p_complementarity,
    prescreen_sort_key,
    resolve_role_inner_reference,
    rewrite_primer_file,
    select_top_candidates,
)
from labcraft.report.renderer import render_redesign_report
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.thermo.salt import UnifiedSaltModel
from labcraft.thermo.vienna import dna_session

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
GT_DIR = REPO_ROOT / "tests" / "fixtures" / "ground_truth"
DUPLEX_PRIMERS = GT_DIR / "duplex_YFV_DENV4_primers.txt"
DUPLEX_TARGETS = GT_DIR / "duplex_YFV_DENV4_targets.fasta"


def test_bidirectional_and_self_3p_complementarity_and_tiebreak():
    """Vérifie le calcul des complémentarités 3' inverse (partenaire -> candidat) et homodimère (candidat -> candidat).
    Verifies reverse (partner 3' -> candidate 5') and self (candidate 3' -> candidate 5') 3' complementarities.
    """
    # Partenaire finissant par 5'-...GACTG-3' -> revcomp(GACTG) = CAGTC
    # Partner ending in 5'-...GACTG-3' -> revcomp(GACTG) = CAGTC
    partner = PhysicalPrimer.from_simple("P_FIP", "TTTTAAAGACTG", PrimerRole.FIP)

    # Candidat A : comp. directe = 2, mais offre la matrice 5' 'A' devant 'CAGTC' (comp. inverse = 5)
    # Candidate A: forward 3' comp = 2, but offers 5' template 'A' upstream of 'CAGTC' (reverse comp = 5)
    cand_high_rev = Candidate(-1, 12, "ACAGTCGGGGAA", Locus(0, 12, "+"), dg_duplex=-9.0, gc_percent=50.0)
    # Candidat B : comp. directe = 2, pas de matrice CAGTC (comp. inverse <= 2), homo <= 2
    # Candidate B: forward 3' comp = 2, no CAGTC template (reverse comp <= 2), self <= 2
    cand_clean = Candidate(+1, 12, "ACACACGGGGAA", Locus(0, 12, "+"), dg_duplex=-9.0, gc_percent=50.0)
    # Candidat C : auto-complémentarité 3' homodimérique de 4 nt ('CCGG' à la fin et en interne)
    # Candidate C: 4-nt homodimeric self-3' complementarity ('CCGG' at 3' and internally)
    cand_self4 = Candidate(0, 12, "AACCGGTTCCGG", Locus(0, 12, "+"), dg_duplex=-9.0, gc_percent=50.0)

    compute_max_3p_complementarity(cand_high_rev, [partner])
    compute_max_3p_complementarity(cand_clean, [partner])
    compute_max_3p_complementarity(cand_self4, [partner])

    assert cand_high_rev.max_reverse_3p_comp == 5
    assert cand_high_rev.max_reverse_3p_partner == "P_FIP"
    assert cand_clean.max_reverse_3p_comp < 5
    assert cand_self4.max_self_3p_comp >= 4

    # À comp. directe, longueur et GC identiques, `cand_clean` doit devancer `cand_high_rev`
    # With equal forward 3' comp, length, and GC, `cand_clean` must rank ahead of `cand_high_rev`
    cand_high_rev.max_3p_comp = 2
    cand_clean.max_3p_comp = 2
    cand_high_rev.gc_percent = 50.0
    cand_clean.gc_percent = 50.0
    ranked = sorted([cand_high_rev, cand_clean], key=lambda c: prescreen_sort_key(c, orig_len=12))
    assert ranked[0] is cand_clean


def test_role_specific_inner_reference_b2_f1c_b1c():
    """Vérifie que F3 utilise F2, B3 utilise B2, LF utilise F1c et LB utilise B1c comme référence interne.
    Verifies F3 uses F2, B3 uses B2, LF uses F1c, and LB uses B1c as inner stability reference.
    """
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    # FIP avec F2 court (faible) et F1c long/GC-riche ; BIP avec B2 très stable (GC-riche)
    # FIP with short/weak F2 and GC-rich F1c; BIP with much more stable GC-rich B2
    fip = PhysicalPrimer(
        "FIP_P1", "GCGCGCGCGCGCGCGCGC" + "ATATATATATATATAT", PrimerRole.FIP,
        binding_domain="ATATATATATATATAT", tail_domain="GCGCGCGCGCGCGCGCGC", parent_name="FIP_P1",
    )
    bip = PhysicalPrimer(
        "BIP_P1", "CGCGCGCGCGCGCGCGCG" + "GGCCGGCCGGCCGGCC", PrimerRole.BIP,
        binding_domain="GGCCGGCCGGCCGGCC", tail_domain="CGCGCGCGCGCGCGCGCG", parent_name="BIP_P1",
    )
    b_kwargs = {"na_mm": 50.0, "mg_mm": 8.0, "dntp_mm": 1.4}

    with dna_session(63.0):
        dom_f3, dg_ref_f3, dg_f2 = resolve_role_inner_reference("F3", [fip, bip], "P1", backend, 63.0, b_kwargs)
        dom_b3, dg_ref_b3, _ = resolve_role_inner_reference("B3", [fip, bip], "P1", backend, 63.0, b_kwargs)
        dom_lf, dg_ref_lf, _ = resolve_role_inner_reference("LF", [fip, bip], "P1", backend, 63.0, b_kwargs)
        dom_lb, dg_ref_lb, _ = resolve_role_inner_reference("LB", [fip, bip], "P1", backend, 63.0, b_kwargs)

    assert dom_f3 == "F2" and dg_ref_f3 == dg_f2
    assert dom_b3 == "B2" and dg_ref_b3 is not None and dg_f2 is not None and dg_ref_b3 < dg_f2 - 5.0
    assert dom_lf == "F1c" and dg_ref_lf is not None
    assert dom_lb == "B1c" and dg_ref_lb is not None

    # Un candidat B3 de ΔG = -10 kcal/mol est plus stable que F2 (~ -4), mais MOINS stable que B2 (~ -18) :
    # il doit être accepté pour B3 !
    # A B3 candidate with ΔG = -10 kcal/mol is more stable than F2 (~ -4) yet LESS stable than B2 (~ -18):
    # it must be accepted for B3!
    cand_b3_ok = Candidate(0, 18, "ATGC" * 4 + "AT", Locus(0, 18, "-"), dg_duplex=-10.0)
    cand_b3_excess = Candidate(1, 20, "GCGC" * 5, Locus(0, 20, "-"), dg_duplex=dg_ref_b3 - 2.0)
    compute_inner_comparison(cand_b3_ok, dom_b3, dg_ref_b3, dg_f2)
    compute_inner_comparison(cand_b3_excess, dom_b3, dg_ref_b3, dg_f2)

    assert cand_b3_ok.more_stable_than_f2 is True
    assert cand_b3_ok.more_stable_than_ref_inner is False
    assert cand_b3_excess.more_stable_than_ref_inner is True

    kept, n_disc = select_top_candidates([cand_b3_ok, cand_b3_excess], top=2, role="B3")
    assert n_disc == 1
    assert kept == [cand_b3_ok]


def test_incremental_evaluator_exact_parity_with_analyze_and_b3_lf_cli(tmp_path: Path):
    """Vérifie la parité numérique exacte entre l'évaluateur incrémental de `redesign` et `labcraft analyze`,
    ainsi que le redessin CLI sur B3 et LF avec template Jinja2.
    Verifies exact numerical parity between `redesign`'s incremental evaluator and `labcraft analyze`,
    plus CLI redesign on B3 and LF with Jinja2 template rendering.
    """
    runner = CliRunner()

    # 1. Redessin de ParidaD4_B3 sur le duplex complet (shift réduit ±2 pour test rapide)
    # 1. Redesign ParidaD4_B3 on the full duplex (compact ±2 shift for fast verification)
    out_b3 = tmp_path / "redesign_b3.html"
    res_b3 = runner.invoke(app, [
        "redesign",
        "-p", str(DUPLEX_PRIMERS),
        "-t", str(DUPLEX_TARGETS),
        "--temperature", "63.0",
        "--conc-fip-bip", "0.4", "--conc-loop", "0.2", "--conc-f3-b3", "0.1",
        "--primer", "ParidaD4_B3",
        "--shift", "2",
        "--length", "18-19",
        "--top", "2",
        "-o", str(out_b3),
    ])
    assert res_b3.exit_code == 0, res_b3.stdout
    js_b3 = json.loads(out_b3.with_suffix(".json").read_text(encoding="utf-8"))
    html_b3 = out_b3.read_text(encoding="utf-8")

    assert js_b3["ref_inner_domain"] == "B2"
    assert "ΔG vs B2" in html_b3
    assert 'id="redesign-ranking"' in html_b3

    # Vérifier la parité exacte avec `labcraft analyze` sur le candidat #1 retenu
    # Check exact parity with `labcraft analyze` on the top-1 selected candidate
    top1 = js_b3["ranked"][1]
    rewritten_primers = tmp_path / "primers_cand1.txt"
    records = read_multi_fasta(str(DUPLEX_PRIMERS))
    rewrite_primer_file(records, "ParidaD4_B3", top1["sequence"], str(rewritten_primers))

    out_an = tmp_path / "analyze_cand1.html"
    res_an = runner.invoke(app, [
        "analyze",
        "-p", str(rewritten_primers),
        "-t", str(DUPLEX_TARGETS),
        "--temperature", "63.0",
        "--conc-fip-bip", "0.4", "--conc-loop", "0.2", "--conc-f3-b3", "0.1",
        "-o", str(out_an),
    ])
    assert res_an.exit_code == 0, res_an.stdout
    js_an = json.loads(out_an.with_suffix(".json").read_text(encoding="utf-8"))

    eq_inc = top1["equilibrium"]
    tot_an = js_an["totals"]
    assert eq_inc["amplifiable_total_nm"] == pytest.approx(tot_an["amplifiable_dimers_total_nm"], rel=1e-9, abs=1e-12)
    assert eq_inc["amplifiable_inter_nm"] == pytest.approx(tot_an["amplifiable_dimers_inter_nm"], rel=1e-9, abs=1e-12)
    assert eq_inc["mispriming_alerts"] == len(js_an["mispriming_risks"])
    assert eq_inc["verdict"] == js_an["verdict"]

    # 2. Redessin d'une amorce de boucle (ParidaD4_LF -> référence F1c)
    # 2. Redesign of a loop primer (ParidaD4_LF -> F1c reference)
    out_lf = tmp_path / "redesign_lf.html"
    res_lf = runner.invoke(app, [
        "redesign",
        "-p", str(GT_DIR / "parida_DENV4.txt"),
        "-t", str(GT_DIR / "Dengue_4_S1028.fasta"),
        "--temperature", "63.0",
        "--primer", "ParidaD4_LF",
        "--shift", "2",
        "--length", "17-19",
        "--top", "2",
        "-o", str(out_lf),
    ])
    assert res_lf.exit_code == 0, res_lf.stdout
    js_lf = json.loads(out_lf.with_suffix(".json").read_text(encoding="utf-8"))
    html_lf = out_lf.read_text(encoding="utf-8")
    assert js_lf["ref_inner_domain"] == "F1c"
    assert "ΔG vs F1c" in html_lf


def test_redesign_jinja2_autoescape():
    """Vérifie que le rendu Jinja2 échappe correctement les caractères spéciaux HTML.
    Verifies Jinja2 rendering properly escapes special HTML characters.
    """
    html = render_redesign_report({
        "primer": "Panel<script>_F3",
        "role_name": "F3",
        "panel_id": "Panel<script>",
        "target_id": "Target&1",
        "orig_locus": Locus(10, 28, "+"),
        "ref_inner_domain": "F2",
        "dg_ref_inner": -12.0,
        "dg_f2": -12.0,
        "n_discarded_inner": 0,
        "shift": 5,
        "len_min": 18,
        "len_max": 20,
        "n_valid": 1,
        "n_total": 1,
        "temp_celsius": 63.0,
        "enzyme_name": "Bst 2.0",
        "rows": [],
    })
    assert "<script>" not in html
    assert "&lt;script&gt;" in html
    assert "Target&amp;1" in html
