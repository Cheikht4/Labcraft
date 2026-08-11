"""Tests d'intégration du rapport HTML / HTML report integration tests.

Vérifie que les sections conditionnelles du rapport apparaissent effectivement
dans le HTML produit, selon le mode d'analyse (cible, sans cible, multiplexe).
Ces tests protègent contre la classe de bugs où une donnée est correctement
calculée mais perdue à l'affichage par un flag incorrect ou un template mal
conditionné.
"""

import pytest
from labcraft.report.renderer import render_report
from labcraft.metrics.verdict import PanelVerdict


def _make_metadata(**overrides):
    """Construit un dict metadata minimal pour le rapport."""
    base = {
        "version": "0.0.1-test",
        "timestamp": "2026-01-01 00:00:00 UTC",
        "file_hash": "deadbeef",
        "max_residual": 1e-12,
        "temperature": 65.0,
        "buffer": "Reference conditions (1.0 M Na+, 0 mM Mg2+)",
        "enzyme": "Bst 2.0",
        "dimer_dg_threshold": -6.0,
        "concentrations_fip_bip": "1.6 µM",
        "concentrations_target": "1 pM",
        "interaction_matrix": {},
        "primer_names": [],
        "unfolding_penalties": {},
        "target_occupations": {},
        "has_true_target": False,
        "primer_to_panel": {},
        "panel_summaries": {},
        "balance_cv": None,
        "mispriming_risks": [],
        "chemistry": "LAMP",
        "optimization_results": {},
        "probe_tm_results": {},
        "recommendations": [],
        "loop_primer_parents": [],
        "primers_parent": [],
        "warnings": [],
    }
    base.update(overrides)
    return base


def _make_verdict(**overrides):
    """Construit un PanelVerdict minimal."""
    defaults = {
        "status": "OK",
        "issues": [],
        "global_cause": "Panel sain.",
    }
    defaults.update(overrides)
    return PanelVerdict(**defaults)


class TestReportTargetMode:
    """Vérifie que les sections cible apparaissent en mode cible."""

    def test_diagnostic_par_cible_present(self):
        """Le titre 'Diagnostic par Cible' est dans le HTML en mode cible."""
        metadata = _make_metadata(has_true_target=True)
        verdict = _make_verdict()
        html = render_report(verdict, {}, [], metadata)
        assert "Diagnostic par Cible" in html

    def test_accessibilite_initiation_present(self):
        """La section d'accessibilité apparaît avec unfolding_penalties renseignées."""
        metadata = _make_metadata(
            has_true_target=True,
            unfolding_penalties={
                "DEN-3": {
                    "F3_DEN3_site": 4.354,
                    "B3_DEN3_site": 3.178,
                    "FIP_DEN3_site": 3.620,
                    "BIP_DEN3_site": 4.270,
                    "LF_DEN3_site": 0.0,
                    "LB_DEN3_site": 0.0,
                }
            },
            target_occupations={
                "DEN-3": {
                    "F3_DEN3_site": 0.12,
                    "B3_DEN3_site": 0.18,
                    "FIP_DEN3_site": 0.15,
                    "BIP_DEN3_site": 0.10,
                    "LF_DEN3_site": 0.0,
                    "LB_DEN3_site": 0.0,
                }
            },
            loop_primer_parents=["LF_DEN3", "LB_DEN3"],
        )
        verdict = _make_verdict()
        html = render_report(verdict, {}, [], metadata)

        # La section doit être présente
        # The section must be present
        assert "Accessibilité à l'initiation" in html

        # Les valeurs de ΔG_unfold pour les amorces d'initiation doivent apparaître
        # ΔG_unfold values for initiation primers must appear
        assert "4.35" in html  # F3
        assert "3.18" in html  # B3
        assert "3.62" in html  # FIP
        assert "4.27" in html  # BIP

        # Les loop primers doivent afficher "n/a"
        # Loop primers must display "n/a"
        assert "n/a" in html


class TestReportNoTargetMode:
    """Vérifie que les sections cible sont absentes en mode sans cible."""

    def test_diagnostic_par_cible_absent(self):
        """Pas de 'Diagnostic par Cible' en mode sans cible."""
        metadata = _make_metadata(has_true_target=False)
        verdict = _make_verdict()
        html = render_report(verdict, {}, [], metadata)
        assert "Diagnostic par Cible" not in html

    def test_accessibilite_absent(self):
        """Pas d'accessibilité en mode sans cible."""
        metadata = _make_metadata(
            has_true_target=False,
            unfolding_penalties={},
        )
        verdict = _make_verdict()
        html = render_report(verdict, {}, [], metadata)
        assert "Accessibilité à l'initiation" not in html


class TestReportMultiplexMode:
    """Vérifie les sections multiplexe."""

    def test_comparaison_panels_present(self):
        """'Comparaison des panels' apparaît en multiplexe."""
        metadata = _make_metadata(
            has_true_target=True,
            balance_cv=0.15,
            panel_summaries={
                "DEN-3": {
                    "mean_occupation": 0.749,
                    "limiting_primer": "F3_DEN3",
                    "min_occupation": 0.082,
                    "mean_free": 0.85,
                },
                "DEN-1": {
                    "mean_occupation": 0.680,
                    "limiting_primer": "B3_DEN1",
                    "min_occupation": 0.050,
                    "mean_free": 0.90,
                },
            },
        )
        verdict = _make_verdict()
        html = render_report(verdict, {}, [], metadata)
        assert "Comparaison des panels" in html
