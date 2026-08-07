import pytest

try:
    import RNA
    HAS_VIENNA = True
except ImportError:
    HAS_VIENNA = False

from labcraft.thermo.vienna import dna_params, rna_params, _VIENNA_LOCK


@pytest.mark.skipif(not HAS_VIENNA, reason="ViennaRNA non installé")
def test_vienna_state_leak():
    """Vérifie que le gestionnaire de contexte restaure correctement l'état.
    
    Un repliement avec les paramètres ADN ne doit pas affecter un repliement ARN
    subséquent.
    """
    seq = "GGGGCCAAAAAAGGCCCC"
    
    # 1. Énergie ARN par défaut (sans context manager)
    RNA.params_load_RNA_Turner2004()
    RNA.cvar.temperature = 37.0
    fc_rna_ref = RNA.fold_compound(seq)
    _, g_rna_ref = fc_rna_ref.pf()
    
    # 2. Utilisation du context manager pour passer en ADN à 65°C
    with dna_params(temp_celsius=65.0):
        fc_dna = RNA.fold_compound(seq)
        _, g_dna = fc_dna.pf()
        # On vérifie que la température a bien été appliquée au niveau de cvar
        assert RNA.cvar.temperature == 65.0
        
    # 3. À la sortie, la température et les paramètres doivent être restaurés
    assert RNA.cvar.temperature == 37.0
    
    # On refait un repliement pour vérifier que les paramètres ARN sont bien restaurés
    fc_rna_post = RNA.fold_compound(seq)
    _, g_rna_post = fc_rna_post.pf()
    
    assert abs(g_rna_ref - g_rna_post) < 1e-4, "L'état des paramètres a fui après le context manager !"
    assert abs(g_rna_ref - g_dna) > 0.1, "L'ADN et l'ARN ne devraient pas avoir la même énergie libre pour cette séquence."


@pytest.mark.skipif(not HAS_VIENNA, reason="ViennaRNA non installé")
def test_vienna_state_leak_deliberate_failure():
    """Ce test prouve qu'une fuite d'état modifie bien silencieusement les calculs
    s'il n'y a pas de gestionnaire.
    """
    seq = "GGGGCCAAAAAAGGCCCC"
    
    # Énergie ARN pure
    RNA.params_load_RNA_Turner2004()
    fc_rna = RNA.fold_compound(seq)
    _, g_rna = fc_rna.pf()
    
    # Modification sauvage de l'état (pas de context manager)
    RNA.params_load_DNA_Mathews2004()
    
    # Repliement d'une structure censée être ARN, mais qui est désormais ADN
    fc_corrupted = RNA.fold_compound(seq)
    _, g_corrupted = fc_corrupted.pf()
    
    assert abs(g_rna - g_corrupted) > 0.1, "La bibliothèque C aurait dû fuir."
    
    # Restauration propre pour les autres tests
    RNA.params_load_RNA_Turner2004()


@pytest.mark.skipif(not HAS_VIENNA, reason="ViennaRNA non installé")
def test_vienna_dna_params_reentrant():
    """Test anti-regression pour le deadlock de dna_params imbriqué."""
    import threading
    
    result = []
    
    def worker():
        try:
            with dna_params(63.0):
                # Premier appel (depth 1)
                with dna_params(63.0):
                    # Deuxième appel imbriqué (depth 2)
                    result.append("success")
        except Exception as e:
            result.append(str(e))
            
    t = threading.Thread(target=worker)
    t.start()
    t.join(timeout=2.0)
    
    assert not t.is_alive(), "Deadlock détecté ! Le thread est bloqué."
    assert result == ["success"], f"Erreur pendant l'exécution: {result}"
