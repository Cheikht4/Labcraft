import pytest

try:
    import RNA
    HAS_VIENNA = True
except ImportError:
    HAS_VIENNA = False

from labcraft.target.unfolding import calc_unfolding_penalty


@pytest.mark.skipif(not HAS_VIENNA, reason="ViennaRNA non installé")
class TestTargetUnfolding:
    
    def test_negative_control_unstructured(self):
        """Un site dans une région totalement simple brin doit donner ~0."""
        # Poly-A n'a aucune structure secondaire
        seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        # On cible le milieu
        dg_unfold = calc_unfolding_penalty(seq, 10, 20, temp_celsius=65.0)
        
        assert abs(dg_unfold) < 1e-3, f"Attendu proche de 0, obtenu {dg_unfold}"
        assert dg_unfold >= 0.0, "La pénalité doit toujours être positive"

    def test_strong_hairpin(self):
        """Un site enfoui dans une forte épingle donne une pénalité égale à l'énergie de l'épingle."""
        # Tige très forte (14 bp), boucle (6 bp)
        stem = "GGGGGGGGGGGGCC"
        loop = "AAAAAA"
        # Épingle complète
        seq = stem + loop + "GGCCCCCCCCCCCC"
        
        # On calcule l'énergie propre de l'épingle libre pour référence
        from labcraft.thermo.vienna import dna_params
        with dna_params(temp_celsius=65.0):
            fc = RNA.fold_compound(seq)
            _, g_free = fc.pf()
            
        # Maintenant on contraint tout le bras 5'
        dg_unfold = calc_unfolding_penalty(seq, 0, len(stem), temp_celsius=65.0)
        
        # Le coût d'ouverture doit être très proche de -g_free
        # puisque ouvrir le bras casse entièrement la tige.
        expected_cost = -g_free
        # On s'attend à un coût d'environ 12 kcal/mol, avec une faible différence entropique résiduelle.
        assert abs(dg_unfold - expected_cost) < 0.5, f"Le coût {dg_unfold} devrait valoir l'énergie de la tige {expected_cost}"

    def test_temperature_monotony(self):
        """La pénalité doit décroître quand la température monte (fusion)."""
        seq = "GCGCGCGCATCG" + "AAAA" + "CGATGCGCGCGC"
        
        dg_50 = calc_unfolding_penalty(seq, 0, 10, temp_celsius=50.0)
        dg_65 = calc_unfolding_penalty(seq, 0, 10, temp_celsius=65.0)
        dg_80 = calc_unfolding_penalty(seq, 0, 10, temp_celsius=80.0)
        
        assert dg_80 < dg_65 < dg_50, "Le coût d'ouverture doit décroître avec la T° (fusion)."

    def test_stem_strength_monotony(self):
        """La pénalité doit croître quand la tige se renforce à température fixe."""
        # 8 bp
        seq_weak = "CGCGATAT" + "AAAA" + "ATATCGCG"
        # 10 bp
        seq_strong = "GCGCGCATCG" + "AAAA" + "CGATGCGCGC"
        
        dg_weak = calc_unfolding_penalty(seq_weak, 0, 8, temp_celsius=65.0)
        dg_strong = calc_unfolding_penalty(seq_strong, 0, 10, temp_celsius=65.0)
        
        assert dg_strong > dg_weak, "Renforcer la tige doit augmenter le coût d'ouverture."

    def test_dengue_synthetic_anchor(self):
        """Ancre synthétique type Dengue (site F3 dans tige stable) : ~8.5 kcal/mol à 65°C."""
        # Tige de 10 paires de bases riche en GC
        seq = "GCGCGCATCG" + "AAAA" + "CGATGCGCGC"
        # On lie l'amorce sur le bras 5' (indices 0 à 10)
        dg_unfold = calc_unfolding_penalty(seq, 0, 10, temp_celsius=65.0)
        
        # Tolérance large pour encadrer ~8.1-8.5
        assert 8.0 < dg_unfold < 9.0, f"Le cas de validation Dengue est faussé: {dg_unfold} kcal/mol"

    def test_index_mapping(self):
        """Vérifie que la bonne base est convertie entre base-0 (Python) et base-1 (ViennaRNA)."""
        # Si on a une tige avec un mismatch au milieu.
        # "GCGC-A-CGCG" -> la position 4 (index 4 en 0-based) est le 'A' non apparié.
        seq = "GCGCA" + "AAAA" + "TGCGC"
        # "GCGCA" vs "TGCGC" -> tige de 4 bp (GCGC...GCGC)
        # La position 4 (A) est dans la boucle ou fait un mismatch si on force la fermeture.
        # Peu importe, l'objectif est que si on contraint la boucle (déjà ouverte), ça coûte ~0
        # Mais si on contraint la tige, ça coûte cher.
        
        # 1. Contrainte sur la boucle (pos 5 à 8, 0-based : AAAA)
        dg_loop = calc_unfolding_penalty(seq, 5, 9, temp_celsius=37.0)
        assert abs(dg_loop) < 0.2, "Ouvrir une boucle ne coûte rien."
        
        # 2. Contrainte sur la tige 5' (pos 0 à 4, 0-based : GCGCA)
        dg_stem = calc_unfolding_penalty(seq, 0, 4, temp_celsius=37.0)
        assert dg_stem > 2.0, "Ouvrir une tige doit coûter cher."
