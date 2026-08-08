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
