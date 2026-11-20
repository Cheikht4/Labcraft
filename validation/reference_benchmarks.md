# Bancs d'essai reproductibles de référence (`reference_benchmarks`)

Document généré automatiquement par `python scripts/reference_benchmarks.py`.
Toute valeur numérique citée dans `CHANGELOG.md` doit provenir exclusivement de ce fichier.

## Métadonnées et conditions expérimentales

- **Commit de code évalué (`src/`)** : `f5cfc96` (`f5cfc9635f45a2f5c4ea68612691a3ff25c7d563`)
- **HEAD lors de l'exécution** : `f5cfc96`
- **Température de réaction (`analyze` & `coverage`)** : `63.0 °C`
- **Tampon LAMP par défaut (`analyze`)** : `50.0 mM Na+`, `10.0 mM K+`, `20.0 mM Tris`, `6.0 mM Mg2+`, `1.4 mM dNTP`
- **Tampon explicite du comparatif `coverage` (Section B)** : `50.0 mM Na+`, `0.0 mM K+`, `0.0 mM Tris`, `8.0 mM Mg2+`, `1.4 mM dNTP`, `--ddg-max 3.0 kcal/mol`
- **Profils enzymatiques (`dimer_dg_threshold`)** : `Bst = -2.0`, `Bst 2.0 = -2.0`, `Bst 2.0 WarmStart = -2.0`, `GspSSD = -2.0`, `Bst 3.0 = -1.5`, `Taq = -3.0` kcal/mol
- **Seuils de verdict et Boltzmann (`analyze`)** : `--dimer-conc-max = 0.230 nM`, `--dimer-warn-ratio = 0.500` (plancher `WARNING = 0.115 nM`), `--role-imbalance-ratio = 20.0`, `boltzmann_mode = exhaustive` par défaut (`--boltzmann-dedup-3p-window = 0`, `--boltzmann-max-dimers = None`, plancher d'éligibilité `0.1 × 0.230 = 0.023 nM`)

### Décalage de `local_dg` (`37.0 °C, 1000 mM Na+` → `63.0 °C, tampon LAMP par défaut`)

| Motif 3' contigu (sans initiation) | `37.0 °C, 1 M Na+` (`kcal/mol`) | `63.0 °C, 1 M Na+` (`kcal/mol`) | `63.0 °C, tampon LAMP` (`kcal/mol`) | Décalage total (`kcal/mol`) |
|---|---:|---:|---:|---:|
| 3 pb (`5'-GGA-3'`) | `-3.14` | `-2.05` | `-1.63` | `+1.51` |
| 4 pb (`5'-GCAT-3'`) | `-4.56` | `-2.81` | `-2.19` | `+2.38` |
| 5 pb (`5'-GCATG-3'`) | `-6.02` | `-3.68` | `-2.85` | `+3.18` |
| 6 pb (`5'-GCATGC-3'`) | `-8.26` | `-5.28` | `-4.24` | `+4.02` |

---

## A. Cinq cas de référence de paillasse (`-2.0` vs `-3.5 kcal/mol`, génomes complets)

| Cas | Paillasse | Dimères amp. `-2.0` (`nM`) | dont inter `-2.0` (`nM`) | Dimères amp. `-3.5` (`nM`) | dont inter `-3.5` (`nM`) | Famille dominante (`-2.0`) | Verdict (`-2.0` / `-3.5`) | `boltzmann_mode` | `suboptimal_3p_masked` (défaut) | Non-amp éligibles évalués | Non-amp éligibles non évalués | Couples parents éligibles |
|---|---|---:|---:|---:|---:|---|---|---|---:|---:|---:|---:|
| `singleplex_denv3` | SUCCÈS paillasse | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `—` | `INFO` / `INFO` | `exhaustive` | `6` | `22` | `0` | `5` |
| `singleplex_denv4` | SUCCÈS paillasse | `0.1536` | `0.0000` | `0.1536` | `0.0000` | `BIP_ParidaD4 × LF_ParidaD4 (0.1536 nM)` | `WARNING` / `WARNING` | `exhaustive` | `0` | `13` | `0` | `13` |
| `singleplex_yfv` | SUCCÈS paillasse | `0.0011` | `0.0000` | `0.0000` | `0.0000` | `B3_Nunes_YFV × LF_Nunes_YFV (0.0011 nM)` | `WARNING` / `WARNING` | `exhaustive` | `1` | `223` | `0` | `6` |
| `duplex_half` | ÉCHEC paillasse (faux +) | `0.3527` | `0.3109` | `0.3363` | `0.2948` | `FIP_Nunes_YFV × F3_ParidaD4 (0.2931 nM)` | `FAILURE` / `FAILURE` | `exhaustive` | `0` | `105` | `0` | `19` |
| `duplex_full` | ÉCHEC paillasse (faux +) | `1.3547` | `1.2010` | `1.2938` | `1.1412` | `FIP_Nunes_YFV × F3_ParidaD4 (1.1351 nM)` | `FAILURE` / `FAILURE` | `exhaustive` | `5` | `353` | `0` | `40` |

### Causes globales (`verdict.global_cause` à `-2.0 kcal/mol`)

- `singleplex_denv3` (`INFO`) : Occupations d'initiation basses mais équilibrées entre amorces de même rôle (comportement attendu sur génome complet), sans dimère amplifiable significatif.
- `singleplex_denv4` (`WARNING`) : Dimères amplifiables modérés (0.154 nM >= plancher WARNING 0.115 nM, < seuil FAILURE 0.230 nM).
- `singleplex_yfv` (`WARNING`) : Occupation réduite par dégénérescence (B3_Nunes_YFV : 16 variants, ~12,5 nM par variant ; déséquilibre 55,8× face à F3_Nunes_YFV).
- `duplex_half` (`FAILURE`) : Présence de dimères d'amorces amplifiables majeurs ou séquestration critique interdisant l'amplification.
- `duplex_full` (`FAILURE`) : Présence de dimères d'amorces amplifiables majeurs ou séquestration critique interdisant l'amplification.

### Détail des dimères `suboptimal_3p_masked = True` (défaut `exhaustive` : `--boltzmann-dedup-3p-window 0`, `--boltzmann-max-dimers None`)

- `singleplex_denv3` : `6` dimère(s) `suboptimal_3p_masked` :
  - `LB_3#1_homo` : `conc = 0.1012 nM`, `P_Boltzmann = 42.76 %` (`max P_Boltzmann` couple parent = `43.20 %`), `ΔG = -5.24 kcal/mol`
  - `LB_3#1_LB_3#2` : `conc = 0.0722 nM`, `P_Boltzmann = 43.20 %` (`max P_Boltzmann` couple parent = `43.20 %`), `ΔG = -5.01 kcal/mol`
  - `F3_3#2_BIP_3#1` : `conc = 0.0404 nM`, `P_Boltzmann = 20.53 %` (`max P_Boltzmann` couple parent = `22.85 %`), `ΔG = -5.09 kcal/mol`
  - `F3_3#2_BIP_3#2` : `conc = 0.0403 nM`, `P_Boltzmann = 17.37 %` (`max P_Boltzmann` couple parent = `22.85 %`), `ΔG = -5.09 kcal/mol`
  - `F3_3#2_BIP_3#3` : `conc = 0.0398 nM`, `P_Boltzmann = 22.85 %` (`max P_Boltzmann` couple parent = `22.85 %`), `ΔG = -5.09 kcal/mol`
  - `F3_3#2_BIP_3#4` : `conc = 0.0398 nM`, `P_Boltzmann = 19.59 %` (`max P_Boltzmann` couple parent = `22.85 %`), `ΔG = -5.09 kcal/mol`
- `singleplex_denv4` : `0` dimère `suboptimal_3p_masked`.
- `singleplex_yfv` : `1` dimère(s) `suboptimal_3p_masked` :
  - `F3_Nunes_YFV#6_LB_Nunes_YFV#2` : `conc = 0.0496 nM`, `P_Boltzmann = 15.17 %` (`max P_Boltzmann` couple parent = `15.17 %`), `ΔG = -6.62 kcal/mol`
- `duplex_half` : `0` dimère `suboptimal_3p_masked`.
- `duplex_full` : `5` dimère(s) `suboptimal_3p_masked` :
  - `F3_Nunes_YFV#6_LB_Nunes_YFV#2` : `conc = 0.0494 nM`, `P_Boltzmann = 15.17 %` (`max P_Boltzmann` couple parent = `15.17 %`), `ΔG = -6.62 kcal/mol`
  - `LF_Nunes_YFV_LB_ParidaD4` : `conc = 0.0324 nM`, `P_Boltzmann = 10.27 %` (`max P_Boltzmann` couple parent = `10.27 %`), `ΔG = -3.55 kcal/mol`
  - `F3_Nunes_YFV#15_FIP_ParidaD4` : `conc = 0.0262 nM`, `P_Boltzmann = 10.74 %` (`max P_Boltzmann` couple parent = `10.74 %`), `ΔG = -5.27 kcal/mol`
  - `F3_Nunes_YFV#7_FIP_ParidaD4` : `conc = 0.0261 nM`, `P_Boltzmann = 10.31 %` (`max P_Boltzmann` couple parent = `10.74 %`), `ΔG = -5.27 kcal/mol`
  - `LF_Nunes_YFV_LF_ParidaD4` : `conc = 0.0253 nM`, `P_Boltzmann = 16.17 %` (`max P_Boltzmann` couple parent = `16.17 %`), `ΔG = -3.51 kcal/mol`

### Comparaison `exhaustive` (défaut, fenêtre `0`) vs `dedup_approx` (fenêtre `5` nt) sur les 5 cas

| Cas | Non-amp éligibles (`exhaustive`) | `suboptimal_3p_masked` (`exhaustive`) | Groupes 3' (`dedup_approx`, `5` nt) | `suboptimal_3p_masked` (`dedup_approx`, `5` nt) | Clés `5` nt où `ΔP_Boltzmann > 10.0` points |
|---|---:|---:|---:|---:|---:|
| `singleplex_denv3` | `22` | `6` | `6` | `2` | `1` |
| `singleplex_denv4` | `13` | `0` | `13` | `0` | `0` |
| `singleplex_yfv` | `223` | `1` | `9` | `1` | `0` |
| `duplex_half` | `105` | `0` | `19` | `0` | `0` |
| `duplex_full` | `353` | `5` | `46` | `4` | `0` |

#### Groupes `(couple parent, 5 nt 3'_A, 5 nt 3'_B)` où `P_Boltzmann` varie de plus de `10.0` points entre variants d'une même clé

- `singleplex_denv3` : `1` groupe(s) avec `ΔP_Boltzmann > 10.0` points entre variants d'une même clé `5` nt :
  - Clé `LB_3[3':TAGGA] × LB_3[3':TAGGA]` (`3` variants) : écart `ΔP = 43.20` points (`LB_3#2_homo` à `0.00 %` → `LB_3#1_LB_3#2` à `43.20 %`) ; représentant retenu par concentration = `LB_3#1_homo` (`0.1012 nM`, `P_Boltzmann = 42.76 %`, `masked = True`). Variants : `LB_3#1_homo` (`0.1012 nM`, `P = 42.76 %`, `masked = True`), `LB_3#1_LB_3#2` (`0.0722 nM`, `P = 43.20 %`, `masked = True`), `LB_3#2_homo` (`0.0510 nM`, `P = 0.00 %`, `masked = False`).
- `singleplex_denv4` : `0` groupe `(couple parent, 5 nt 3'_A, 5 nt 3'_B)` avec `ΔP_Boltzmann > 10.0` points.
- `singleplex_yfv` : `0` groupe `(couple parent, 5 nt 3'_A, 5 nt 3'_B)` avec `ΔP_Boltzmann > 10.0` points.
- `duplex_half` : `0` groupe `(couple parent, 5 nt 3'_A, 5 nt 3'_B)` avec `ΔP_Boltzmann > 10.0` points.
- `duplex_full` : `0` groupe `(couple parent, 5 nt 3'_A, 5 nt 3'_B)` avec `ΔP_Boltzmann > 10.0` points.

---

## B. Comparatif `labcraft coverage` : commit `58ae5d9` vs commit `8cf5bd8`

Cohortes synthétiques de non-régression générées par `scripts/generate_coverage_cohorts.py` (graine fixe `SEED = 20261119` — mutants synthétiques et non souches cliniques de terrain) :
- `synthetic_cohort_denv3_120.fasta` : `120` séquences, `sha256 = 7bd32b881e44e19acd3f73971fc46bfcc2a42ccc31160185c20bc6f295c8b74f`
- `synthetic_cohort_denv4_120.fasta` : `120` séquences, `sha256 = f1b65f222ae76bdcab2f38137a5cc8da012c408f607b08cde9d69e5a5bc41125`

| Jeu évalué | `58ae5d9` (`COUVERT / A_RISQUE / ECHEC`) | `8cf5bd8` (`COUVERT / A_RISQUE / ECHEC`) | Changements de classe | Changements `SiteVerdict` | Évaluations `|Δ(ΔΔG)| > 0` |
|---|---|---|---|---:|---:|
| Génome GenBank réf. DENV3 (`DEN3_M93130.fasta`, panel 3, `--errors 2`) | `0 / 1 / 0` (`100.0 %`) | `0 / 1 / 0` (`100.0 %`) | `0 / 1` | `0` | `0 / 6` |
| Génome GenBank réf. DENV4 (`Dengue_4_S1028.fasta`, Parida, `--errors 2`) | `0 / 1 / 0` (`100.0 %`) | `0 / 1 / 0` (`100.0 %`) | `0 / 1` | `0` | `0 / 6` |
| Génome GenBank réf. YFV (`YFV_17D_17DD.fasta`, Nunes, `--errors 2`) | `1 / 0 / 0` (`100.0 %`) | `1 / 0 / 0` (`100.0 %`) | `0 / 1` | `0` | `0 / 6` |
| Cohorte synthétique DENV3 (`120` mutants, `--errors 2`, défaut CLI) | `0 / 67 / 53` (`55.8 %`) | `0 / 67 / 53` (`55.8 %`) | `0 / 120` | `0` | `26 / 712` |
| Cohorte synthétique DENV4 (`120` mutants, `--errors 2`, défaut CLI) | `0 / 69 / 51` (`57.5 %`) | `0 / 69 / 51` (`57.5 %`) | `0 / 120` | `0` | `25 / 710` |
| Cohorte synthétique DENV3 (`120` mutants, `--errors 3`, boucles `3×3` incluses) | `0 / 67 / 53` (`55.8 %`) | `0 / 67 / 53` (`55.8 %`) | `0 / 120` | `0` | `34 / 720` |
| Cohorte synthétique DENV4 (`120` mutants, `--errors 3`, boucles `3×3` incluses) | `0 / 69 / 51` (`57.5 %`) | `0 / 69 / 51` (`57.5 %`) | `0 / 120` | `0` | `35 / 720` |

### Top 5 des plus grands `|Δ(ΔΔG)|` (`8cf5bd8 - 58ae5d9`, à `63.0 °C`) sur les cohortes synthétiques

#### Cohorte synthétique DENV3 (`120` mutants, `--errors 2`, défaut CLI)
- `M93130.1_2mm_adj2x2_02_3_B2` [`BIP`, `2` mm] : `ddG(58ae5d9) = +6.1682` → `ddG(8cf5bd8) = +6.7102 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)
- `M93130.1_2mm_adj2x2_04_3_F3` [`F3`, `2` mm] : `ddG(58ae5d9) = +3.9498` → `ddG(8cf5bd8) = +4.4918 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)
- `M93130.1_2mm_adj2x2_10_3_B2` [`BIP`, `2` mm] : `ddG(58ae5d9) = +4.9792` → `ddG(8cf5bd8) = +5.5211 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)
- `M93130.1_2mm_adj2x2_13_3_F2` [`FIP`, `2` mm] : `ddG(58ae5d9) = +4.0238` → `ddG(8cf5bd8) = +4.5657 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)
- `M93130.1_2mm_adj2x2_19_3_FLOOP` [`LF`, `2` mm] : `ddG(58ae5d9) = +5.6976` → `ddG(8cf5bd8) = +6.2395 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)

#### Cohorte synthétique DENV4 (`120` mutants, `--errors 2`, défaut CLI)
- `MG272272.1_2mm_adj2x2_19_ParidaD4_B3` [`B3`, `2` mm] : `ddG(58ae5d9) = +4.2733` → `ddG(8cf5bd8) = +4.8152 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)
- `MG272272.1_2mm_adj2x2_00_ParidaD4_F3` [`F3`, `2` mm] : `ddG(58ae5d9) = +3.5346` → `ddG(8cf5bd8) = +4.0765 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)
- `MG272272.1_2mm_adj2x2_03_ParidaD4_F3` [`F3`, `2` mm] : `ddG(58ae5d9) = +4.7615` → `ddG(8cf5bd8) = +5.3034 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)
- `MG272272.1_2mm_adj2x2_09_ParidaD4_F3` [`F3`, `2` mm] : `ddG(58ae5d9) = +5.0574` → `ddG(8cf5bd8) = +5.5994 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)
- `MG272272.1_2mm_adj2x2_11_ParidaD4_FLP` [`LF`, `2` mm] : `ddG(58ae5d9) = +12.4879` → `ddG(8cf5bd8) = +13.0298 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)

#### Cohorte synthétique DENV3 (`120` mutants, `--errors 3`, boucles `3×3` incluses)
- `M93130.1_2mm_adj2x2_02_3_B2` [`BIP`, `2` mm] : `ddG(58ae5d9) = +6.1682` → `ddG(8cf5bd8) = +6.7102 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)
- `M93130.1_2mm_adj2x2_04_3_F3` [`F3`, `2` mm] : `ddG(58ae5d9) = +3.9498` → `ddG(8cf5bd8) = +4.4918 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)
- `M93130.1_2mm_adj2x2_10_3_B2` [`BIP`, `2` mm] : `ddG(58ae5d9) = +4.9792` → `ddG(8cf5bd8) = +5.5211 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)
- `M93130.1_2mm_adj2x2_13_3_F2` [`FIP`, `2` mm] : `ddG(58ae5d9) = +4.0238` → `ddG(8cf5bd8) = +4.5657 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)
- `M93130.1_2mm_adj2x2_19_3_FLOOP` [`LF`, `2` mm] : `ddG(58ae5d9) = +5.6976` → `ddG(8cf5bd8) = +6.2395 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)

#### Cohorte synthétique DENV4 (`120` mutants, `--errors 3`, boucles `3×3` incluses)
- `MG272272.1_2mm_adj2x2_19_ParidaD4_B3` [`B3`, `2` mm] : `ddG(58ae5d9) = +4.2733` → `ddG(8cf5bd8) = +4.8152 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)
- `MG272272.1_3mm_adj3x3_01_ParidaD4_B3` [`B3`, `3` mm] : `ddG(58ae5d9) = +6.7385` → `ddG(8cf5bd8) = +7.2804 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)
- `MG272272.1_3mm_adj3x3_03_ParidaD4_F3` [`F3`, `3` mm] : `ddG(58ae5d9) = +5.9586` → `ddG(8cf5bd8) = +6.5005 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)
- `MG272272.1_3mm_adj3x3_09_ParidaD4_F3` [`F3`, `3` mm] : `ddG(58ae5d9) = +13.8831` → `ddG(8cf5bd8) = +14.4250 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)
- `MG272272.1_2mm_adj2x2_00_ParidaD4_F3` [`F3`, `2` mm] : `ddG(58ae5d9) = +3.5346` → `ddG(8cf5bd8) = +4.0765 kcal/mol` (`Δ(ΔΔG) = +0.5419 kcal/mol`, `non_viable` → `non_viable`)

---

## C. Performance `enrich_dimer_boltzmann_3p` sur 81 espèces (`3 321` paires, `N = 3` répétitions : médiane et étendue)

| Scénario | Mode | `boltzmann_mode` | Appels Boltzmann | Groupes/dimères non-amp éligibles | Dimères éligibles non évalués | `pf_time` (`s`, médiane / étendue) | Total hors Boltzmann (`s`, médiane / étendue) | Total avec Boltzmann (`s`, médiane / étendue) | Part médiane `pf / avec` | Surcoût médian `pf / hors` |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 81 espèces — fenêtres locales (`test_step3_redesign_v2`) | défaut `exhaustive` (`--boltzmann-dedup-3p-window 0`, `--boltzmann-max-dimers None`) | `exhaustive` | `384` | `353` | `0` | `4.689 (étendue 0.019, [4.676–4.695])` | `13.63 (étendue 3.28, [13.60–16.88])` | `18.41 (étendue 0.43, [18.35–18.78])` | `25.47 %` | `+34.41 %` |
| 81 espèces — fenêtres locales (`test_step3_redesign_v2`) | approx. `dedup_approx` 5 nt (`--boltzmann-dedup-3p-window 5`, `--boltzmann-max-dimers None`) | `dedup_approx` | `77` | `46` | `307` | `0.474 (étendue 0.027, [0.457–0.484])` | `13.63 (étendue 3.28, [13.60–16.88])` | `14.48 (étendue 0.20, [14.36–14.56])` | `3.27 %` | `+3.48 %` |
| 81 espèces — fenêtres locales (`test_step3_redesign_v2`) | approx. `dedup_approx` 5 nt + K=25 (`--boltzmann-dedup-3p-window 5`, `--boltzmann-max-dimers 25`) | `dedup_approx` | `56` | `46` | `328` | `0.398 (étendue 0.019, [0.393–0.412])` | `13.63 (étendue 3.28, [13.60–16.88])` | `14.12 (étendue 0.67, [14.03–14.70])` | `2.82 %` | `+2.92 %` |
| 81 espèces — `duplex_full` génomes viraux complets | défaut `exhaustive` (`--boltzmann-dedup-3p-window 0`, `--boltzmann-max-dimers None`) | `exhaustive` | `384` | `353` | `0` | `4.653 (étendue 0.066, [4.634–4.700])` | `13.82 (étendue 0.38, [13.57–13.95])` | `18.38 (étendue 0.93, [18.35–19.28])` | `25.31 %` | `+33.66 %` |
| 81 espèces — `duplex_full` génomes viraux complets | approx. `dedup_approx` 5 nt (`--boltzmann-dedup-3p-window 5`, `--boltzmann-max-dimers None`) | `dedup_approx` | `77` | `46` | `307` | `0.452 (étendue 0.059, [0.452–0.511])` | `13.82 (étendue 0.38, [13.57–13.95])` | `14.56 (étendue 1.51, [14.36–15.86])` | `3.10 %` | `+3.27 %` |
| 81 espèces — `duplex_full` génomes viraux complets | approx. `dedup_approx` 5 nt + K=25 (`--boltzmann-dedup-3p-window 5`, `--boltzmann-max-dimers 25`) | `dedup_approx` | `56` | `46` | `328` | `0.390 (étendue 0.009, [0.390–0.398])` | `13.82 (étendue 0.38, [13.57–13.95])` | `14.21 (étendue 0.92, [13.99–14.92])` | `2.75 %` | `+2.82 %` |
