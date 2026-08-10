"""Tests pour la cohérence du rapport et le traitement des variants dégénérés.

Tests for report coherence and degenerate variant handling.
"""
import pytest
import tempfile
import yaml
import os

from labcraft.optimize.recommendations import generate_recommendations
from labcraft.metrics.verdict import PanelVerdict, PrimerIssue, generate_verdict
from labcraft.metrics.risk import RiskItem
from labcraft.metrics.balance import calculate_multiplex_balance
from labcraft.diagnostics.mispriming import detect_inter_target_mispriming, _revcomp
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole, expand_degenerate
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.thermo.salt import UnifiedSaltModel
from labcraft.diagnostics.enzyme import ENZYME_REGISTRY


# ============================
# Partie 1 : Recommandations cohérentes
# ============================

class TestRecommendationsCoherence:
    """Vérifie que generate_recommendations utilise les risques réels."""
    
    def test_amplifiable_dimers_trigger_recommendations(self):
        """Un panel avec des dimères amplifiables ne doit PAS être qualifié de sain."""
        # Verdict sans issues (les dimères ne remontent pas dans verdict.issues)
        verdict = PanelVerdict(status="OK", issues=[], global_cause="Tout va bien.")
        
        # Mais il y a un dimère amplifiable dans les risques
        risks = [
            RiskItem(
                complex_name="FIP_x:LF_x",
                concentration=5e-9,
                severity=10.0,  # Amplifiable
                description="Dimère amplifiable (3' extensible)"
            )
        ]
        
        recs = generate_recommendations(verdict, risks=risks)
        
        # Ne doit PAS contenir "panel semble sain"
        assert not any("sain" in r.lower() for r in recs), \
            f"Recommandations contradictoires : {recs}"
        # Doit contenir la recommandation WarmStart
        assert any("WarmStart" in r for r in recs)
        
    def test_no_risks_panel_sain(self):
        """Un panel sans aucun dimère dangereux doit afficher 'panel sain'."""
        verdict = PanelVerdict(status="OK", issues=[], global_cause="OK")
        risks = []
        
        recs = generate_recommendations(verdict, risks=risks)
        assert any("sain" in r.lower() for r in recs)
        
    def test_blocking_dimers_trigger_concentration_advice(self):
        """Des dimères bloquants (severity < 10) déclenchent le conseil de concentration."""
        verdict = PanelVerdict(status="WARNING", issues=[], global_cause="Warning")
        risks = [
            RiskItem(
                complex_name="F3_x:B3_x",
                concentration=5e-9,
                severity=1.0,  # Bloquant, pas amplifiable
                description="Dimère bloquant"
            )
        ]
        
        recs = generate_recommendations(verdict, risks=risks)
        assert any("concentration" in r.lower() for r in recs)
        assert not any("sain" in r.lower() for r in recs)

    def test_backward_compatible_without_risks(self):
        """Sans argument risks, fallback sur le verdict (rétrocompatibilité)."""
        verdict = PanelVerdict(status="OK", issues=[], global_cause="OK")
        recs = generate_recommendations(verdict)
        assert any("sain" in r.lower() for r in recs)


# ============================
# Partie 3a : Mésamorçage et variants
# ============================

class TestMisprimingVariants:
    """Vérifie que les variants ne sont pas exclus de la détection de mésamorçage."""
    
    @pytest.fixture
    def backend(self):
        return ViennaSaltShiftBackend(UnifiedSaltModel())
    
    @pytest.fixture
    def enzyme(self):
        return ENZYME_REGISTRY["bst2.0"]
    
    def test_variant_not_excluded_from_mispriming(self, backend, enzyme):
        """Un variant (#1) ne doit PAS être exclu du scan de mésamorçage."""
        seq = "ATCGATCGATCGATCG"
        # Créer le variant avec nom#1
        p = PhysicalPrimer.from_simple("F3_X#1", seq, PrimerRole.F3, parent_name="F3_X")
        
        target_seq = "NNNN" + _revcomp(seq) + "NNNN"
        
        # primer_to_panel contient le parent ET le variant
        primer_to_panel = {"F3_X": "X", "F3_X#1": "X"}
        targets = {"YFV": target_seq}
        
        risks = detect_inter_target_mispriming(
            [p], primer_to_panel, targets, backend, enzyme, temp_celsius=65.0
        )
        
        # Le variant doit être détecté (même ancrage 3' que le parent)
        assert len(risks) > 0, "Le variant a été exclu du scan de mésamorçage !"
        assert risks[0].primer_name == "F3_X#1"
        
    def test_degenerate_base_does_not_lose_mispriming(self, backend, enzyme):
        """Remplacer un A par un R ne doit PAS faire disparaître un risque de mésamorçage."""
        seq_pure = "ATCGATCGATCGATCG"
        
        # Version pure (sans dégénérescence)
        p_pure = PhysicalPrimer.from_simple("F3_X", seq_pure, PrimerRole.F3)
        target_seq = "NNNN" + _revcomp(seq_pure) + "NNNN"
        primer_to_panel_pure = {"F3_X": "X"}
        targets = {"YFV": target_seq}
        
        risks_pure = detect_inter_target_mispriming(
            [p_pure], primer_to_panel_pure, targets, backend, enzyme, temp_celsius=65.0
        )
        
        # Version avec R (A→R, les variants contiennent A et G)
        # Le variant #1 (A) est identique à la séquence pure
        p_variant = PhysicalPrimer.from_simple("F3_X#1", seq_pure, PrimerRole.F3, parent_name="F3_X")
        primer_to_panel_variant = {"F3_X": "X", "F3_X#1": "X"}
        
        risks_variant = detect_inter_target_mispriming(
            [p_variant], primer_to_panel_variant, targets, backend, enzyme, temp_celsius=65.0
        )
        
        # Le variant A doit détecter le MÊME risque que la version pure
        assert len(risks_pure) == len(risks_variant), \
            f"Mésamorçage perdu avec variant : {len(risks_pure)} vs {len(risks_variant)}"


# ============================
# Partie 3b : Balance multiplexe avec variants
# ============================

class TestBalanceVariants:
    """Vérifie que la balance multiplexe fonctionne correctement avec des variants."""
    
    def test_balance_with_variants(self):
        """min_occupation et limiting_primer doivent être corrects avec des variants."""
        # Deux variants de FIP_T1
        primer_to_panel = {
            "FIP_T1": "T1",
            "FIP_T1#1": "T1",
            "FIP_T1#2": "T1",
            "BIP_T1": "T1",
        }
        
        target_occupations = {
            "T1": {
                "FIP_T1_site": 0.15,
                "BIP_T1_site": 0.20,
            }
        }
        
        free_fractions = {
            "FIP_T1#1": 0.50,
            "FIP_T1#2": 0.45,
            "BIP_T1": 0.60,
        }
        
        summaries, cv = calculate_multiplex_balance(
            primer_to_panel, target_occupations, free_fractions
        )
        
        assert "T1" in summaries
        s = summaries["T1"]
        # min_occupation doit être 0.15 (FIP_T1), pas 0.0
        assert s["min_occupation"] > 0.0, f"min_occupation est {s['min_occupation']}, devrait être > 0"
        assert abs(s["min_occupation"] - 0.15) < 0.01
        # L'amorce limitante doit être FIP_T1 (parent), pas "N/A"
        assert s["limiting_primer"] == "FIP_T1"


# ============================
# Partie 3c : Sous-domaines FIP/BIP développés pour variants
# ============================

class TestSubdomainExpansion:
    """Vérifie que les sous-domaines de FIP/BIP sont correctement développés pour les variants."""
    
    def test_subdomain_expanded_by_position(self):
        """Les sous-domaines doivent être extraits par position dans la séquence du variant."""
        from labcraft.cli.config import PanelConfig, build_engine_from_config, PrimerDomains
        
        # Séquence FIP avec un R (= A ou G) dans F2
        # Structure: 5'-F1c-linker-F2-3'
        # F1c = AAAAACCCCC (10 nt)
        # linker = TT (2 nt)
        # F2 = GGGGGAARRG (10 nt, avec 2 R = A|G → 4 variants)
        fip_seq = "AAAAACCCCCTTGGGGGAARRG"
        f1c = "AAAAACCCCC"
        f2 = "GGGGGAARRG"  # contient 2 R
        
        variants = expand_degenerate(fip_seq)
        # Avec 2 R, ça fait 2^2 = 4 variants
        assert len(variants) == 4
        
        # Pour chaque variant, le F2 doit être cohérent (pas de R littéral)
        # et doit être extrait par position
        for v in variants:
            v_f2 = v[12:]  # positions 12..21 (après F1c + linker)
            assert 'R' not in v_f2, f"R littéral trouvé dans F2 du variant : {v_f2}"
            # Les positions 3 et 4 dans F2 (qui étaient R) doivent être A ou G
            assert v_f2[3] in ('A', 'G'), f"Position 3 inattendue : {v_f2[3]}"
            assert v_f2[4] in ('A', 'G'), f"Position 4 inattendue : {v_f2[4]}"
