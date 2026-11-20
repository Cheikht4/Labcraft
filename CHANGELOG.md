# Changelog

Toutes les modifications notables apportées à **LabCraft** sont documentées dans ce fichier.
Toutes les valeurs numériques de bancs d'essai, seuils et comparatifs ci-dessous sont copiées directement depuis la sortie versionnée `validation/reference_benchmarks.md`, générée par `scripts/reference_benchmarks.py` (et `scripts/generate_coverage_cohorts.py`).

Pour chaque modification, l'impact est distingué en :
- **(a) Équilibre, couverture ou verdict** : modifie les énergies libres `ΔG` injectées dans le solveur multi-espèces, les concentrations à l'équilibre, la couverture multi-souches ou la classification globale du panel (`OK` / `INFO` / `WARNING` / `FAILURE`).
- **(b) Affichage, `Tm`/`ΔH`/`ΔS`, JSON, CLI ou outillage seulement** : enrichit les métadonnées rapportées, la documentation, le typage ou l'empaquetage sans modifier ni l'équilibre thermodynamique ni le verdict.

---

## [Unreleased] — Étapes 1 à 6 & Corrections de relecture (`a044d3f` → `HEAD`)

### Corrections de relecture post-Étape 6 (`fa39fff` → `HEAD`)

- **Point 1 — Provenance heuristique des seuils enzymatiques `dimer_dg_threshold` et intervalle de séparation `--dimer-conc-max`** (`fa39fff`)
  - **(a) Équilibre / Verdict** : aucun changement des seuils codés (`Bst = -2.0`, `Bst 2.0 = -2.0`, `Bst 2.0 WarmStart = -2.0`, `GspSSD = -2.0`, `Bst 3.0 = -1.5`, `Taq = -3.0` kcal/mol dans `src/labcraft/diagnostics/enzyme.py:61-117` ; `--dimer-conc-max = 0.230 nM` dans `src/labcraft/cli/main.py:131, 495`). Sur les cinq cas de référence de paillasse à `63.0 °C` (`validation/reference_benchmarks.md:27-36` et test de régression `@pytest.mark.regression` dans `tests/test_diagnostics/test_point1_dg3p_threshold_coherence.py:98-154`), les concentrations totales de dimères amplifiables sont identiques ou quasi identiques entre `-2.0` et `-3.5` kcal/mol :
    - `singleplex_denv3` : `0.0000` nM (`-2.0`) vs `0.0000` nM (`-3.5`), verdict `INFO` / `INFO` ;
    - `singleplex_denv4` : `0.1536` nM (`-2.0`) vs `0.1536` nM (`-3.5`), verdict `WARNING` / `WARNING` ;
    - `singleplex_yfv` : `0.0011` nM (`-2.0`) vs `0.0000` nM (`-3.5`), verdict `WARNING` / `WARNING` ;
    - `duplex_half` : `0.3527` nM dont `0.3109` nM inter-panels (`-2.0`) vs `0.3363` nM dont `0.2948` nM inter-panels (`-3.5`), verdict `FAILURE` / `FAILURE` ;
    - `duplex_full` : `1.3547` nM dont `1.2010` nM inter-panels (`-2.0`) vs `1.2938` nM dont `1.1412` nM inter-panels (`-3.5`), verdict `FAILURE` / `FAILURE`.
  - **(b) Affichage / Documentation / CLI** : documentation explicite en heuristique interne non calibrée (`src/labcraft/diagnostics/enzyme.py:36-58`, `src/labcraft/cli/main.py:72-90`, `README.md:12-13`), décalage de `local_dg` entre (`37.0 °C, 1 M Na+`) et (`63.0 °C`, tampon LAMP `50.0 mM Na+`, `10.0 mM K+`, `20.0 mM Tris`, `6.0 mM Mg2+`, `1.4 mM dNTP`) sans pénalité d'initiation (`-3.14` → `-1.63` soit `+1.51` kcal/mol sur `3` pb `5'-GGA-3'`, `-4.56` → `-2.19` soit `+2.38` kcal/mol sur `4` pb `5'-GCAT-3'`, `-6.02` → `-2.85` soit `+3.18` kcal/mol sur `5` pb `5'-GCATG-3'`, `-8.26` → `-4.24` soit `+4.02` kcal/mol sur `6` pb `5'-GCATGC-3'` ; `validation/reference_benchmarks.md:16-24`), retrait de la comparaison avec PrimerExplorer, et précision que `duplex_half` et `duplex_full` sont deux dilutions du même panel dont la séparation repose sur l'unique paire inter-panels `FIP_Nunes_YFV × F3_ParidaD4` (`0.2931 nM` sur `duplex_half`, `1.1351 nM` sur `duplex_full` ; `validation/reference_benchmarks.md:34-35`).

- **Point 2 & Point A — Alignement des boucles internes et bulges sur SantaLucia & Hicks et statut `couple_loops_in_salt_dg=False`** (`8cf5bd8`)
  - **(a) Équilibre / Couverture / Verdict** :
    - *Équilibre `analyze` / `redesign`* : **aucun effet**. `couple_loops_in_salt_dg` vaut `False` par défaut (`src/labcraft/thermo/backends/vienna_salt.py:375`) et `build_engine_from_config` instancie `ViennaSaltShiftBackend(UnifiedSaltModel())` sans l'activer (`src/labcraft/cli/config.py:141, 325`) afin d'éviter un double comptage des boucles déjà incluses par ViennaRNA. Sur cible parfaitement appariée, `calculate_hybridization_dg` renvoie directement `res_hyb.dg_kcal` (`src/labcraft/thermo/mismatch.py:399-407`).
    - *Couverture `labcraft coverage` (`58ae5d9` vs `8cf5bd8`)* : `nn_duplex_energy` (`src/labcraft/thermo/mismatch.py:132-270`) utilise désormais la table exacte des boucles internes et l'extrapolation de Jacobson-Stockmayer de SantaLucia & Hicks (`src/labcraft/thermo/mismatch.py:31-66`). Le comparatif reproductible (`validation/reference_benchmarks.md:91-135`, généré via `scripts/generate_coverage_cohorts.py` avec la graine fixe `SEED = 20261119` à `63.0 °C`, `50.0 mM Na+`, `8.0 mM Mg2+`, `1.4 mM dNTP`, `--ddg-max 3.0 kcal/mol`) porte sur les trois génomes GenBank de référence (`tests/fixtures/ground_truth/DEN3_M93130.fasta`, `tests/fixtures/ground_truth/Dengue_4_S1028.fasta`, `tests/fixtures/ground_truth/YFV_17D_17DD.fasta` : `0 / 1` changement chacun, `0 / 6` évaluations modifiées) et sur deux cohortes de non-régression de **mutants synthétiques** (`120` mutants DENV3 + `120` mutants DENV4, et non des souches cliniques de terrain) :
      - Cohorte synthétique DENV3 (`120` mutants) : `0 / 67 / 53` (`55.8 %`) sous `58ae5d9` comme sous `8cf5bd8`, soit `0 / 120` changement de classe et `0` changement de `SiteVerdict` (`26 / 712` évaluations modifiées à `--errors 2`, `34 / 720` à `--errors 3`) ;
      - Cohorte synthétique DENV4 (`120` mutants) : `0 / 69 / 51` (`57.5 %`) sous `58ae5d9` comme sous `8cf5bd8`, soit `0 / 120` changement de classe et `0` changement de `SiteVerdict` (`25 / 710` évaluations modifiées à `--errors 2`, `35 / 720` à `--errors 3`) ;
      - Toutes les évaluations modifiées (boucles internes contiguës `2×2` et `3×3`) présentent le décalage uniforme `Δ(ΔΔG) = +0.5419 kcal/mol` à `63.0 °C` (`validation/reference_benchmarks.md:107-135`).
  - **(b) Affichage / `Tm` / `ΔH` / `ΔS`** : alignement exact de `_BULGE_LOOP_TABLE4_DG37_KCAL` et `_TERMINAL_AT_PENALTY` dans `src/labcraft/thermo/backends/vienna_salt.py:28-80`, validé contre Primer3 dans `tests/test_thermo/test_step5_bulges_and_boltzmann_3p.py:346-410`.

- **Point C — Déplacement de `primer3-py` dans les dépendances optionnelles `dev`** (`fc97587`)
  - **(a) Équilibre / Verdict** : aucun effet (`primer3` n'est importé par aucun module de `src/labcraft/`).
  - **(b) Empaquetage / Tests** : déplacement de `primer3-py` hors de `[project].dependencies` vers `[project.optional-dependencies].dev` (`pyproject.toml:27-53`) pour préserver la licence MIT d'exécution, avec garde `pytest.importorskip("primer3")` vérifiée par `tests/test_thermo/test_primer3_optional.py:1-63`.

- **Point 3, filtre relatif, dédoublonnage 3' `--boltzmann-dedup-3p-window` & défaut illimité `--boltzmann-max-dimers` — Probabilité de Boltzmann 3' conditionnelle intermoléculaire** (`526ddba`, `e52306b`, `51dd563`, `e3b8d4c`, `81919b2`, `824eb09`, `604cd50`)
  - **(a) Équilibre / Verdict** : **aucun effet**. `calc_boltzmann_3p_pairing` et `enrich_dimer_boltzmann_3p` (`src/labcraft/diagnostics/amplifiable_dimer.py:49-229, 459-532`) n'entrent jamais dans `solve_dual` ni dans `generate_verdict`. Comme `2` dimères (`LB_3#1_homo` à `0.1012 nM`, `P_Boltzmann = 42.76 %`, et `F3_3#2_BIP_3#1` à `0.0404 nM`, `P_Boltzmann = 20.53 %`, soit `6` variants de `0.0398` à `0.1012 nM` avant fusion des doublons 5'/milieu ; `validation/reference_benchmarks.md:47-66`) sont signalés `suboptimal_3p_masked = True` sur `singleplex_denv3` alors que ce panel fonctionne en paillasse, la documentation (`README.md:14`, `src/labcraft/diagnostics/amplifiable_dimer.py:15-23`, `src/labcraft/report/templates/report.html.j2:1048`) précise explicitement que cet indicateur ne doit pas entrer dans le verdict tant qu'il n'est pas validé expérimentalement.
  - **(b) Affichage / JSON / Performance** :
    - Séparation entre auto-appariements intramoléculaires (`prob_3p_intra_a`, `prob_3p_intra_b`) et probabilité intermoléculaire conditionnelle sur le sous-ensemble connecté `A:B` avec respect du veto ARMS (`src/labcraft/diagnostics/amplifiable_dimer.py:49-229`).
    - Seuil relatif d'éligibilité `0.1 × 0.230 = 0.023 nM` dérivé de `--dimer-conc-max`, dédoublonnage par clé canonique symétrique `(couple parent, 3'_A, 3'_B)` sur `--boltzmann-dedup-3p-window = 5` nt par défaut (fusionnant sur leur variant le plus concentré les variants ne différant qu'en 5' ou au milieu tout en évaluant séparément ceux qui diffèrent dans les `5` dernières bases 3'), défaut `--boltzmann-max-dimers = None` illimité (`0` ou `None` = sans plafond ; `K > 0` conservé comme garde-fou optionnel pour les très gros panels), et rapport de la `P_Boltzmann` maximale par couple parent (`parent_pair_max_boltzmann_3p_prob` et `boltzmann_max_prob_by_parent_pair` ; `src/labcraft/pipeline/analysis.py:114-131, 272-375`, `src/labcraft/cli/main.py:96-106, 134-135`, `src/labcraft/metrics/risk.py:29`, `src/labcraft/report/json_export.py:113-117, 141-146`, `src/labcraft/report/templates/report.html.j2:1055`, `tests/test_thermo/test_step5_bulges_and_boltzmann_3p.py:524-873`).
    - Réconciliation complète sur les cinq cas de référence (`validation/reference_benchmarks.md:45-87`) : avec les paramètres par défaut, LabCraft retrouve les `7` dimères `suboptimal_3p_masked` distincts en 3' (`2 + 0 + 1 + 0 + 4 = 7`) issus des `12` variants bruts sans dédoublonnage (`5` doublons ne différant qu'en 5'/milieu fusionnés) :
      - `singleplex_denv3` (`2` dimères sur `6` groupes 3' éligibles) : `LB_3#1_homo` (`0.1012 nM`, `42.76 %`, fusionnant `2` variants identiques en `TAGGA × TAGGA`) et `F3_3#2_BIP_3#1` (`0.0404 nM`, `20.53 %`, fusionnant `3` variants identiques en `AATCA × GGGAG`) ;
      - `singleplex_yfv` (`1` dimère sur `9` groupes 3' éligibles) : `F3_Nunes_YFV#6_LB_Nunes_YFV#2` (`0.0496 nM` / `0.04963376 nM`, `CACTG × AGTGA`, `P_Boltzmann = 15.17 %`, `max P_Boltzmann` couple parent = `15.17 %`), évalué séparément de `F3_Nunes_YFV#5_LB_Nunes_YFV#2` (`0.04963523 nM`, `CACTA × AGTGA`, `P_Boltzmann = 0.02 %`) ;
      - `duplex_full` (`4` dimères sur `46` groupes 3' éligibles) : `F3_Nunes_YFV#6_LB_Nunes_YFV#2` (`0.0494 nM`, `15.17 %`), `LF_Nunes_YFV_LB_ParidaD4` (`0.0324 nM`, `10.27 %`), `F3_Nunes_YFV#15_FIP_ParidaD4` (`0.0262 nM` / `0.02615681 nM`, `CATTA × AAACC`, `P_Boltzmann = 10.74 %`, évalué séparément de `#16` `CATTG × AAACC` à `0.02615694 nM` et `0.46 %`, et fusionnant le doublon 5' `F3_Nunes_YFV#7_FIP_ParidaD4` à `0.0261 nM` et `10.31 %`) et `LF_Nunes_YFV_LF_ParidaD4` (`0.0253 nM`, `16.17 %`).
    - Sur `81` espèces (`3 321` paires, `N = 3` répétitions entrelacées ; `validation/reference_benchmarks.md:139-146`), le défaut illimité (`77` appels Boltzmann, `46` groupes 3' éligibles, `0` non évalué, `40` couples parents éligibles) donne un `pf_time` médian de `0.456` s (étendue `0.015` s, `[0.455–0.470]` s, soit `3.25 %` du total médian `14.04` s d'étendue `0.15` s et `+3.37 %` rapporté à `13.52` s d'étendue `2.25` s hors Boltzmann) sur fenêtres locales, et de `0.452` s (étendue `0.005` s, `[0.451–0.456]` s, soit `3.22 %` de `14.03` s d'étendue `0.10` s et `+3.34 %` rapporté à `13.53` s d'étendue `0.10` s) sur `duplex_full` ; avec le garde-fou `--boltzmann-max-dimers 25` (`56` appels, `21` groupes non évalués), `pf_time` vaut `0.392` s (`+2.90 %`) et `0.391` s (`+2.89 %`).

- **Point B & formatage à 3 décimales — Niveau `INFO`, déséquilibre relatif entre rôles jumeaux (Option B2) et plancher relatif de dimères** (`cf947db`, `e52306b`)
  - **(a) Verdict uniquement (concentrations d'équilibre inchangées)** :
    - Ajout du niveau `INFO` (`src/labcraft/metrics/verdict.py:76, 335-346`, sortie CLI de code `0` dans `src/labcraft/cli/main.py:293-309`), du seuil relatif `--role-imbalance-ratio = 20.0` sans normalisation par le nombre de variants (`src/labcraft/metrics/verdict.py:21, 134-188`) distinguant la dilution par dégénérescence (`singleplex_yfv` : `B3_Nunes_YFV : 16 variants, ~12,5 nM par variant ; déséquilibre 55,8× face à F3_Nunes_YFV` ; `validation/reference_benchmarks.md:41`, `src/labcraft/metrics/verdict.py:246-267`) du blocage thermodynamique, et du plancher relatif `--dimer-warn-ratio = 0.500` (`WARNING = 0.115 nM` pour `0.230 nM` ; `validation/reference_benchmarks.md:14`, `src/labcraft/metrics/verdict.py:22, 300-305`).
  - **(b) Affichage / HTML / CLI** : affichage à trois décimales (`0.115 nM`, `0.230 nM`, `0.154 nM`) dans `global_cause` (`src/labcraft/metrics/verdict.py:54-56, 327-341`), `src/labcraft/report/templates/report.html.j2:419-422` et `src/labcraft/cli/main.py:81-90` (`tests/test_metrics/test_verdict_info_and_imbalance.py:90-146`).

- **Point 4 — Avertissements explicites sur `calc_hairpin` et LNA** (`2a177fb`)
  - **(a) Équilibre / Verdict** : inchangé en fonctionnement nominal.
  - **(b) Diagnostics / Avertissements** : `RuntimeWarning` explicite sur toute erreur inattendue hors `NotImplementedError` (`src/labcraft/lamp/complex_enumeration.py:397-409`, `src/labcraft/optimize/redesign.py:884-896, 1063-1075`) et `UserWarning` lorsque `ViennaRNABackend.calc_hairpin` reçoit `lna_positions` non vide (`src/labcraft/thermo/backends/vienna.py:106-114`).

---

### Étapes 1 à 6 initiales (`a044d3f` → `58ae5d9`)

- **Pré-étape / Phase 3 — Session ADN ViennaRNA unique et départage LAMP** (`a044d3f`)
  - **(a) Sélection `redesign`** : critères de départage LAMP dans `src/labcraft/optimize/redesign.py` ; équilibre `analyze` inchangé.
  - **(b) Performance** : gestionnaire de contexte réentrant `dna_session` dans `src/labcraft/thermo/vienna.py:173-181`.

- **Étape 1 — Veto ARMS du second brin, ensemencement bilatéral et normalisation multiplex** (`4a36e17`)
  - **(a) Équilibre / Couverture / Verdict** : correction de l'indice partenaire antiparallèle du second brin dans `is_amplifiable_dimer` (`src/labcraft/diagnostics/amplifiable_dimer.py:390-412`), découpage symétrique inverse-complémentaire et saut des blocs de `N` dans `src/labcraft/target/seeding.py:82-125`, normalisation par cible dans `src/labcraft/optimize/concentrations.py` et propagation de `PolymeraseProfile` dans `src/labcraft/optimize/constraints.py` (`tests/test_step1_reliability_fixes.py`).
  - **(b) Affichage** : alignement du second brin dans `src/labcraft/report/alignment.py`.

- **Étape 2 — Mésappariements terminaux, boucles internes symétriques et option wobble G·T** (`6afcf8e`, affinée en `8cf5bd8`)
  - **(a) Couverture / Hybridation imparfaite** : branchement de `_MISMATCH_TERMINAL` et des boucles internes symétriques dans `nn_duplex_energy` (`src/labcraft/thermo/mismatch.py:21, 81-90, 132-270`) et option `allow_penultimate_gt_wobble` (`False` par défaut dans `src/labcraft/diagnostics/enzyme.py:33`).
  - **(b) Qualification** : fonction `classify_mismatch_pair` dans `src/labcraft/thermo/mismatch.py:109-129`.

- **Étape 3 — Évaluateur incrémental en mémoire et gabarit Jinja2 pour `labcraft redesign`** (`0e243be`)
  - **(a) Sélection `redesign`** : `IncrementalRedesignEvaluator` résolvant l'équilibre couplé en mémoire (`src/labcraft/optimize/redesign.py:724-1160`).
  - **(b) Affichage HTML** : extraction du rapport vers `src/labcraft/report/templates/redesign.html.j2` (`src/labcraft/report/renderer.py:53-64`).

- **Étape 4 — Pipeline d'analyse découplé, optimiseur salin et `mypy --strict` en CI** (`f8ddaea`, `7dfe331`, `fbe9955`)
  - **(a) Optimiseur local** : support de `enzyme`, `backend` et `**backend_kwargs` dans `optimize_primer` (`src/labcraft/optimize/optimizer.py:15-24, 77-85`) ; équilibre `analyze` inchangé.
  - **(b) Architecture / CI** : extraction de `src/labcraft/pipeline/analysis.py` et activation stricte de `mypy` (`pyproject.toml:85-94`, `.github/workflows/ci.yml`, `tests/test_step4_pipeline_and_optimizer.py`).

- **Étape 5 — Thermodynamique des boucles et probabilités de Boltzmann** (`94579a3`, corrigée en `8cf5bd8`, `526ddba`, `e52306b`, `51dd563`, `81919b2`, `824eb09`, `604cd50`)
  - **(a) Équilibre / Verdict** : **aucun effet** (`couple_loops_in_salt_dg=False` dans `src/labcraft/thermo/backends/vienna_salt.py:375` ; probabilités de Boltzmann purement informatives dans `src/labcraft/pipeline/analysis.py:272-375`).
  - **(b) Affichage / `Tm` / `ΔH` / `ΔS` / JSON** : cohérence de Gibbs-Helmholtz sur les grandeurs rapportées (`src/labcraft/thermo/backends/vienna_salt.py:390-455`) et champs Boltzmann dans `src/labcraft/report/json_export.py:110-125, 141-146`.

- **Étape 6 — Épingles intramoléculaires monomères et connexité intermoléculaire des dimères** (`58ae5d9`, affinée en `2a177fb`)
  - **(a) Équilibre / Verdict** : ajout des espèces monomoléculaires `{p.name}_hairpin` dans le bilan de masse (`src/labcraft/lamp/complex_enumeration.py:387-413`, `src/labcraft/optimize/redesign.py:875-899, 1052-1079`) et imposition d'au moins une paire intermoléculaire via `_has_intermolecular_pair` et `_duplexfold_connected` (`src/labcraft/thermo/backends/vienna.py:19-93`, `src/labcraft/metrics/risk.py:53`).
  - **(b) Affichage / JSON** : exposition de la fraction `"hairpin"` dans `src/labcraft/metrics/fractions.py:9, 88` et `src/labcraft/report/json_export.py:48`.
