# Validation / Validation

This directory contains tests validating LabCraft's thermodynamic models.
Ce répertoire contient les tests de validation des modèles thermodynamiques de LabCraft.

## Layers of Validation / Couches de validation

### Layer A / Couche A
Validates the native Python thermodynamic solver against established open-source tools (ViennaRNA, Primer3).
Valide le solveur thermodynamique Python natif par rapport aux outils open-source établis (ViennaRNA, Primer3).

### Layer B / Couche B
Validates multi-strand interacting systems and advanced equilibrium calculations against NUPACK.
Valide les systèmes multi-brins en interaction et les calculs d'équilibre avancés par rapport à NUPACK.

**Note:** NUPACK is an optional dependency used exclusively for Layer B validation and is never required to run LabCraft.
**Note :** NUPACK est une dépendance optionnelle utilisée exclusivement pour la validation de la Couche B et n'est jamais requise pour exécuter LabCraft.
