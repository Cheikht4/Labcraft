import pytest
from labcraft.metrics.fractions import PrimerFractions
from labcraft.metrics.verdict import generate_verdict

def test_causal_verdict_attribution():
    """
    Vérifie que la fonction de verdict attribue les bonnes causes aux 
    bonnes amorces selon leurs profils de fractions physiques.
    """
    fractions = {
        "F3_A": PrimerFractions(
            free=0.95, # Très libre
            hairpin=0.0,
            homodimer=0.0,
            heterodimer_intra=0.0,
            heterodimer_inter=0.0,
            target_bound=0.05,
            dominant_complex="F3_A_free",
            dominant_fraction=0.95
        ),
        "FIP_A": PrimerFractions(
            free=0.01, # Quasiment aucune amorce libre
            hairpin=0.0,
            homodimer=0.99, # Totalement séquestrée
            heterodimer_intra=0.0,
            heterodimer_inter=0.0,
            target_bound=0.0,
            dominant_complex="FIP_A_homo",
            dominant_fraction=0.99
        ),
        "FIP_B": PrimerFractions(
            free=0.05,
            hairpin=0.0,
            homodimer=0.0,
            heterodimer_intra=0.0,
            heterodimer_inter=0.90, # Hybridation croisée !
            target_bound=0.05,
            dominant_complex="FIP_B_BIP_A",
            dominant_fraction=0.90
        )
    }
    
    # 0.05 (5%) d'occupation pour tous, donc tous sont en difficulté (< 10%)
    target_occupations = {
        "SynthA": {
            "F3_A_site": 0.05,
            "FIP_A_site": 0.0
        },
        "SynthB": {
            "FIP_B_site": 0.05
        }
    }
    
    verdict = generate_verdict(fractions, target_occupations, risks=[])
    
    assert verdict.status == "FAILURE"
    assert len(verdict.issues) == 3
    
    issues_by_primer = {i.primer_name: i.cause for i in verdict.issues}
    
    # F3_A est libre mais non liée -> Inaccessibilité
    assert "Inaccessibilité du site cible pour F3_A" in issues_by_primer["F3_A"]
    
    # FIP_A est dans un homodimère -> Séquestration intrapanel
    assert "Séquestration de FIP_A à 99.0% dans FIP_A_homo" in issues_by_primer["FIP_A"]
    
    # FIP_B est dans un hétérodimère interpanel -> Hybridation croisée
    assert "Hybridation croisée inter-jeux séquestrant FIP_B à 90.0% dans FIP_B_BIP_A" in issues_by_primer["FIP_B"]

from labcraft.metrics.risk import evaluate_risks, RiskItem

def test_risk_item_alignment_columns():
    """Vérifie que les colonnes d'alignement sont bien transmises au RiskItem"""
    complexes = ["dimer_A_B"]
    concs = [1e-6]
    flags = [True]
    details = [{
        "seq_a": "ATGC",
        "seq_b": "GCAT",
        "structure": "((&))",
        "delta_g": -3.0,
        "delta_g_3p": -3.0,
        "alignment_columns": [
            {"top": "A", "bottom": "T", "paired": True, "role": "three_prime"},
            {"top": "-", "bottom": "G", "paired": False, "role": "template"},
            {"top": "-", "bottom": "C", "paired": False, "role": "template"}
        ],
        "arrow_metrics": {"show": True, "margin_cols": 0, "width_cols": 2}
    }]
    
    risks = evaluate_risks(complexes, concs, flags, is_warm_start=False, dimer_details=details)
    assert len(risks) == 1
    assert len(risks[0].alignment_columns) == 3
    assert risks[0].arrow_metrics["show"] is True
    assert risks[0].arrow_metrics["margin_cols"] == 0
    assert risks[0].arrow_metrics["width_cols"] == 2


def test_part_d_dimer_conc_max_and_recommendations_coherence_via_cli(tmp_path):
    """
    Partie D : Vérifie via CliRunner sur `labcraft analyze` :
    1. L'option `--dimer-conc-max` gouverne la frontière WARNING / FAILURE sur la concentration
       totale de dimères amplifiables (Panel 3 sur M93130.1 a ~1.09 nM de dimères amplifiables :
       WARNING avec --dimer-conc-max 2.5, FAILURE avec --dimer-conc-max 0.5).
    2. Un panel en WARNING ou FAILURE ne reçoit jamais le message
       « Aucune recommandation non-séquentielle critique. Le panel semble sain. ».
    3. `is_warm_start` a pour défaut `False` dans `evaluate_risks` et la CLI transmet `enzyme.warm_start`.
    """
    from pathlib import Path
    from typer.testing import CliRunner
    from labcraft.cli.main import app
    from labcraft.cli.parsers import read_multi_fasta
    from labcraft.optimize.recommendations import generate_recommendations
    from labcraft.metrics.verdict import PanelVerdict, PrimerIssue

    # Vérification du défaut is_warm_start=False dans evaluate_risks
    # Verify default is_warm_start=False in evaluate_risks
    r_default = evaluate_risks(["dimer_1"], [1e-7], [True])
    assert "WarmStart" not in r_default[0].description

    # Vérification qu'un verdict WARNING ou FAILURE sans dimères ne dit jamais "Le panel semble sain."
    # Verify that a WARNING or FAILURE verdict without dimers never says "Le panel semble sain."
    for st in ("WARNING", "FAILURE"):
        v_test = PanelVerdict(
            status=st,
            issues=[PrimerIssue("F3", "T1", 0.0, "Site absent", is_critical=(st == "FAILURE"))],
            global_cause="Problème sur cible."
        )
        recs = generate_recommendations(v_test, risks=[])
        assert all("Le panel semble sain" not in rec for rec in recs)

    # Test CLI réel sur Parida DEN-3 + M93130.1 (~7.26 nM de dimères amplifiables à 63 °C avec BIP corrigée)
    # Real CLI test on Parida DEN-3 + M93130.1 (~7.26 nM amplifiable dimers at 63 °C with fixed BIP orientation)
    fixture_dir = Path(__file__).resolve().parent.parent / "fixtures" / "real_primer_files"
    parida_file = fixture_dir / "amorceparida.txt"
    d3_fasta = fixture_dir / "DEN3_M93130.fasta"

    pd3_records = [
        f">{h}\n{s}" for h, s in read_multi_fasta(str(parida_file)) if h.startswith("ParidaD3_")
    ]
    primers_only_pd3 = tmp_path / "parida_d3.txt"
    primers_only_pd3.write_text("\n".join(pd3_records) + "\n", encoding="utf-8")

    runner = CliRunner()
    # 1) Avec --dimer-conc-max 20.0 (au-dessus de 7.26 nM) -> WARNING
    out_warn = tmp_path / "pd3_warn.html"
    res_warn = runner.invoke(app, [
        "analyze", "-p", str(primers_only_pd3), "-t", str(d3_fasta),
        "--temperature", "63.0", "--dimer-conc-max", "20.0", "-o", str(out_warn)
    ])
    assert res_warn.exit_code == 0, f"CLI failed: {res_warn.stdout}"
    html_warn = out_warn.read_text(encoding="utf-8")
    assert "Verdict Global : WARNING" in html_warn
    assert "Le panel semble sain" not in html_warn

    # 2) Avec --dimer-conc-max 0.1 (en-dessous de 7.26 nM) -> FAILURE
    out_fail = tmp_path / "pd3_fail.html"
    res_fail = runner.invoke(app, [
        "analyze", "-p", str(primers_only_pd3), "-t", str(d3_fasta),
        "--temperature", "63.0", "--dimer-conc-max", "0.1", "-o", str(out_fail)
    ])
    assert res_fail.exit_code == 0, f"CLI failed: {res_fail.stdout}"
    html_fail = out_fail.read_text(encoding="utf-8")
    assert "Verdict Global : FAILURE" in html_fail
    assert "Le panel semble sain" not in html_fail

