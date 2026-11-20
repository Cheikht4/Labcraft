# Bancs d'essai reproductibles de référence (`reference_benchmarks`)

Document généré automatiquement par `python scripts/reference_benchmarks.py`.
Toute valeur numérique citée dans `CHANGELOG.md` doit provenir exclusivement de ce fichier.

## Métadonnées et conditions expérimentales

- **Commit de code évalué (`src/`)** : `604cd50` (`604cd5092123f927dd97d6782709fbd8ae2d2db4`)
- **HEAD lors de l'exécution** : `604cd50`
- **Température de réaction (`analyze` & `coverage`)** : `63.0 °C`
- **Tampon LAMP par défaut (`analyze`)** : `50.0 mM Na+`, `10.0 mM K+`, `20.0 mM Tris`, `6.0 mM Mg2+`, `1.4 mM dNTP`
- **Tampon explicite du comparatif `coverage` (Section B)** : `50.0 mM Na+`, `0.0 mM K+`, `0.0 mM Tris`, `8.0 mM Mg2+`, `1.4 mM dNTP`, `--ddg-max 3.0 kcal/mol`
- **Profils enzymatiques (`dimer_dg_threshold`)** : `Bst = -2.0`, `Bst 2.0 = -2.0`, `Bst 2.0 WarmStart = -2.0`, `GspSSD = -2.0`, `Bst 3.0 = -1.5`, `Taq = -3.0` kcal/mol
- **Seuils de verdict et Boltzmann (`analyze`)** : `--dimer-conc-max = 0.230 nM`, `--dimer-warn-ratio = 0.500` (plancher `WARNING = 0.115 nM`), `--role-imbalance-ratio = 20.0`, `--boltzmann-dedup-3p-window = 5` nt, `--boltzmann-max-dimers = None` (illimité par défaut ; plancher d'éligibilité `0.1 × 0.230 = 0.023 nM`)

### Décalage de `local_dg` (`37.0 °C, 1000 mM Na+` → `63.0 °C, tampon LAMP par défaut`)

| Motif 3' contigu (sans initiation) | `37.0 °C, 1 M Na+` (`kcal/mol`) | `63.0 °C, 1 M Na+` (`kcal/mol`) | `63.0 °C, tampon LAMP` (`kcal/mol`) | Décalage total (`kcal/mol`) |
|---|---:|---:|---:|---:|
| 3 pb (`5'-GGA-3'`) | `-3.14` | `-2.05` | `-1.63` | `+1.51` |
| 4 pb (`5'-GCAT-3'`) | `-4.56` | `-2.81` | `-2.19` | `+2.38` |
| 5 pb (`5'-GCATG-3'`) | `-6.02` | `-3.68` | `-2.85` | `+3.18` |
| 6 pb (`5'-GCATGC-3'`) | `-8.26` | `-5.28` | `-4.24` | `+4.02` |

---

## A. Cinq cas de référence de paillasse (`-2.0` vs `-3.5 kcal/mol`, génomes complets)

| Cas | Paillasse | Dimères amp. `-2.0` (`nM`) | dont inter `-2.0` (`nM`) | Dimères amp. `-3.5` (`nM`) | dont inter `-3.5` (`nM`) | Famille dominante (`-2.0`) | Verdict (`-2.0` / `-3.5`) | `suboptimal_3p_masked` (défaut) | Groupes 3' éligibles | Groupes 3' non évalués | Couples parents éligibles |
|---|---|---:|---:|---:|---:|---|---|---:|---:|---:|---:|
| `singleplex_denv3` | SUCCÈS paillasse | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `—` | `INFO` / `INFO` | `2` | `6` | `0` | `5` |
| `singleplex_denv4` | SUCCÈS paillasse | `0.1536` | `0.0000` | `0.1536` | `0.0000` | `BIP_ParidaD4 × LF_ParidaD4 (0.1536 nM)` | `WARNING` / `WARNING` | `0` | `13` | `0` | `13` |
| `singleplex_yfv` | SUCCÈS paillasse | `0.0011` | `0.0000` | `0.0000` | `0.0000` | `B3_Nunes_YFV × LF_Nunes_YFV (0.0011 nM)` | `WARNING` / `WARNING` | `1` | `9` | `0` | `6` |
| `duplex_half` | ÉCHEC paillasse (faux +) | `0.3527` | `0.3109` | `0.3363` | `0.2948` | `FIP_Nunes_YFV × F3_ParidaD4 (0.2931 nM)` | `FAILURE` / `FAILURE` | `0` | `19` | `0` | `19` |
| `duplex_full` | ÉCHEC paillasse (faux +) | `1.3547` | `1.2010` | `1.2938` | `1.1412` | `FIP_Nunes_YFV × F3_ParidaD4 (1.1351 nM)` | `FAILURE` / `FAILURE` | `4` | `46` | `0` | `40` |

### Causes globales (`verdict.global_cause` à `-2.0 kcal/mol`)

- `singleplex_denv3` (`INFO`) : Occupations d'initiation basses mais équilibrées entre amorces de même rôle (comportement attendu sur génome complet), sans dimère amplifiable significatif.
- `singleplex_denv4` (`WARNING`) : Dimères amplifiables modérés (0.154 nM >= plancher WARNING 0.115 nM, < seuil FAILURE 0.230 nM).
- `singleplex_yfv` (`WARNING`) : Occupation réduite par dégénérescence (B3_Nunes_YFV : 16 variants, ~12,5 nM par variant ; déséquilibre 55,8× face à F3_Nunes_YFV).
- `duplex_half` (`FAILURE`) : Présence de dimères d'amorces amplifiables majeurs ou séquestration critique interdisant l'amplification.
- `duplex_full` (`FAILURE`) : Présence de dimères d'amorces amplifiables majeurs ou séquestration critique interdisant l'amplification.

### Détail des dimères `suboptimal_3p_masked = True` (défaut : `--boltzmann-dedup-3p-window 5`, `--boltzmann-max-dimers None`)

- `singleplex_denv3` : `2` dimère(s) `suboptimal_3p_masked` :
  - `LB_3#1_homo` : `conc = 0.1012 nM`, `P_Boltzmann = 42.76 %` (`max P_Boltzmann` couple parent = `42.76 %`), `ΔG = -5.24 kcal/mol`
  - `F3_3#2_BIP_3#1` : `conc = 0.0404 nM`, `P_Boltzmann = 20.53 %` (`max P_Boltzmann` couple parent = `20.53 %`), `ΔG = -5.09 kcal/mol`
- `singleplex_denv4` : `0` dimère `suboptimal_3p_masked`.
- `singleplex_yfv` : `1` dimère(s) `suboptimal_3p_masked` :
  - `F3_Nunes_YFV#6_LB_Nunes_YFV#2` : `conc = 0.0496 nM`, `P_Boltzmann = 15.17 %` (`max P_Boltzmann` couple parent = `15.17 %`), `ΔG = -6.62 kcal/mol`
- `duplex_half` : `0` dimère `suboptimal_3p_masked`.
- `duplex_full` : `4` dimère(s) `suboptimal_3p_masked` :
  - `F3_Nunes_YFV#6_LB_Nunes_YFV#2` : `conc = 0.0494 nM`, `P_Boltzmann = 15.17 %` (`max P_Boltzmann` couple parent = `15.17 %`), `ΔG = -6.62 kcal/mol`
  - `LF_Nunes_YFV_LB_ParidaD4` : `conc = 0.0324 nM`, `P_Boltzmann = 10.27 %` (`max P_Boltzmann` couple parent = `10.27 %`), `ΔG = -3.55 kcal/mol`
  - `F3_Nunes_YFV#15_FIP_ParidaD4` : `conc = 0.0262 nM`, `P_Boltzmann = 10.74 %` (`max P_Boltzmann` couple parent = `10.74 %`), `ΔG = -5.27 kcal/mol`
  - `LF_Nunes_YFV_LF_ParidaD4` : `conc = 0.0253 nM`, `P_Boltzmann = 16.17 %` (`max P_Boltzmann` couple parent = `16.17 %`), `ΔG = -3.51 kcal/mol`

### Réconciliation des `7` dimères retenus (`2 + 0 + 1 + 0 + 4 = 7`) avec les `12` variants bruts sans dédoublonnage

- `singleplex_denv3` : sans dédoublonnage = `6` dimère(s) ; dédoublonnage `(couple parent, 5 nt 3'_A, 5 nt 3'_B)` par défaut (illimité) = `2` dimère(s) sur `6` groupe(s) 3' éligible(s) (`5` couple(s) parent(s) éligible(s), `0` groupe non évalué ; avec garde-fou `K = 25` = `2`).
  - Couple parent `LB_3 × LB_3` (`1` groupe(s) 3' distinct(s), `max P_Boltzmann` couple parent = `42.76 %`) :
    - Clé `LB_3[3':TAGGA] × LB_3[3':TAGGA]` (rang `3/6`) : représentant `LB_3#1_homo` (`0.10122321 nM` / `0.1012 nM`, `ΔG = -5.2371 kcal/mol`, `P_Boltzmann = 42.76 %`, `masked = True`) ; fusionne `2` variant(s) identique(s) sur les `5` nt 3' ne différant qu'en 5'/milieu (LB_3#1_LB_3#2 à 0.0722 nM (P = 43.20 %), LB_3#2_homo à 0.0510 nM (P = 0.00 %)).
  - Couple parent `BIP_3 × F3_3` (`1` groupe(s) 3' distinct(s), `max P_Boltzmann` couple parent = `20.53 %`) :
    - Clé `BIP_3[3':AATCA] × F3_3[3':GGGAG]` (rang `4/6`) : représentant `F3_3#2_BIP_3#1` (`0.04035122 nM` / `0.0404 nM`, `ΔG = -5.0864 kcal/mol`, `P_Boltzmann = 20.53 %`, `masked = True`) ; fusionne `3` variant(s) identique(s) sur les `5` nt 3' ne différant qu'en 5'/milieu (F3_3#2_BIP_3#2 à 0.0403 nM (P = 17.37 %), F3_3#2_BIP_3#3 à 0.0398 nM (P = 22.85 %), F3_3#2_BIP_3#4 à 0.0398 nM (P = 19.59 %)).
- `singleplex_denv4` : sans dédoublonnage = `0` dimère(s) ; dédoublonnage `(couple parent, 5 nt 3'_A, 5 nt 3'_B)` par défaut (illimité) = `0` dimère(s) sur `13` groupe(s) 3' éligible(s) (`13` couple(s) parent(s) éligible(s), `0` groupe non évalué ; avec garde-fou `K = 25` = `0`).
- `singleplex_yfv` : sans dédoublonnage = `1` dimère(s) ; dédoublonnage `(couple parent, 5 nt 3'_A, 5 nt 3'_B)` par défaut (illimité) = `1` dimère(s) sur `9` groupe(s) 3' éligible(s) (`6` couple(s) parent(s) éligible(s), `0` groupe non évalué ; avec garde-fou `K = 25` = `1`).
  - Couple parent `F3_Nunes_YFV × LB_Nunes_YFV` (`4` groupe(s) 3' distinct(s), `max P_Boltzmann` couple parent = `15.17 %`) :
    - Clé `F3_Nunes_YFV[3':CACTA] × LB_Nunes_YFV[3':AGTGA]` (rang `5/9`) : représentant `F3_Nunes_YFV#5_LB_Nunes_YFV#2` (`0.04963523 nM` / `0.0496 nM`, `ΔG = -6.6209 kcal/mol`, `P_Boltzmann = 0.02 %`, `masked = False`).
    - Clé `F3_Nunes_YFV[3':CACTG] × LB_Nunes_YFV[3':AGTGA]` (rang `6/9`) : représentant `F3_Nunes_YFV#6_LB_Nunes_YFV#2` (`0.04963376 nM` / `0.0496 nM`, `ΔG = -6.6209 kcal/mol`, `P_Boltzmann = 15.17 %`, `masked = True`).
    - Clé `F3_Nunes_YFV[3':CATTG] × LB_Nunes_YFV[3':AGTGA]` (rang `8/9`) : représentant `F3_Nunes_YFV#8_LB_Nunes_YFV#2` (`0.02388413 nM` / `0.0239 nM`, `ΔG = -6.1309 kcal/mol`, `P_Boltzmann = 1.93 %`, `masked = False`).
    - Clé `F3_Nunes_YFV[3':CATTA] × LB_Nunes_YFV[3':AGTGA]` (rang `9/9`) : représentant `F3_Nunes_YFV#7_LB_Nunes_YFV#2` (`0.02388402 nM` / `0.0239 nM`, `ΔG = -6.1309 kcal/mol`, `P_Boltzmann = 0.00 %`, `masked = False`).
- `duplex_half` : sans dédoublonnage = `0` dimère(s) ; dédoublonnage `(couple parent, 5 nt 3'_A, 5 nt 3'_B)` par défaut (illimité) = `0` dimère(s) sur `19` groupe(s) 3' éligible(s) (`19` couple(s) parent(s) éligible(s), `0` groupe non évalué ; avec garde-fou `K = 25` = `0`).
- `duplex_full` : sans dédoublonnage = `5` dimère(s) ; dédoublonnage `(couple parent, 5 nt 3'_A, 5 nt 3'_B)` par défaut (illimité) = `4` dimère(s) sur `46` groupe(s) 3' éligible(s) (`40` couple(s) parent(s) éligible(s), `0` groupe non évalué ; avec garde-fou `K = 25` = `0`).
  - Couple parent `F3_Nunes_YFV × LB_Nunes_YFV` (`4` groupe(s) 3' distinct(s), `max P_Boltzmann` couple parent = `15.17 %`) :
    - Clé `F3_Nunes_YFV[3':CACTA] × LB_Nunes_YFV[3':AGTGA]` (rang `27/46`) : représentant `F3_Nunes_YFV#5_LB_Nunes_YFV#2` (`0.04936754 nM` / `0.0494 nM`, `ΔG = -6.6209 kcal/mol`, `P_Boltzmann = 0.02 %`, `masked = False`).
    - Clé `F3_Nunes_YFV[3':CACTG] × LB_Nunes_YFV[3':AGTGA]` (rang `28/46`) : représentant `F3_Nunes_YFV#6_LB_Nunes_YFV#2` (`0.04936436 nM` / `0.0494 nM`, `ΔG = -6.6209 kcal/mol`, `P_Boltzmann = 15.17 %`, `masked = True`).
    - Clé `F3_Nunes_YFV[3':CATTG] × LB_Nunes_YFV[3':AGTGA]` (rang `45/46`) : représentant `F3_Nunes_YFV#8_LB_Nunes_YFV#2` (`0.02374708 nM` / `0.0237 nM`, `ΔG = -6.1309 kcal/mol`, `P_Boltzmann = 1.93 %`, `masked = False`).
    - Clé `F3_Nunes_YFV[3':CATTA] × LB_Nunes_YFV[3':AGTGA]` (rang `46/46`) : représentant `F3_Nunes_YFV#7_LB_Nunes_YFV#2` (`0.02374697 nM` / `0.0237 nM`, `ΔG = -6.1309 kcal/mol`, `P_Boltzmann = 0.00 %`, `masked = False`).
  - Couple parent `LB_ParidaD4 × LF_Nunes_YFV` (`1` groupe(s) 3' distinct(s), `max P_Boltzmann` couple parent = `10.27 %`) :
    - Clé `LB_ParidaD4[3':CGCAG] × LF_Nunes_YFV[3':GTATG]` (rang `36/46`) : représentant `LF_Nunes_YFV_LB_ParidaD4` (`0.03241857 nM` / `0.0324 nM`, `ΔG = -3.5502 kcal/mol`, `P_Boltzmann = 10.27 %`, `masked = True`).
  - Couple parent `F3_Nunes_YFV × FIP_ParidaD4` (`2` groupe(s) 3' distinct(s), `max P_Boltzmann` couple parent = `10.74 %`) :
    - Clé `F3_Nunes_YFV[3':CATTG] × FIP_ParidaD4[3':AAACC]` (rang `42/46`) : représentant `F3_Nunes_YFV#16_FIP_ParidaD4` (`0.02615694 nM` / `0.0262 nM`, `ΔG = -5.2681 kcal/mol`, `P_Boltzmann = 0.46 %`, `masked = False`) ; fusionne `1` variant(s) identique(s) sur les `5` nt 3' ne différant qu'en 5'/milieu (F3_Nunes_YFV#8_FIP_ParidaD4 à 0.0261 nM (P = 0.39 %)).
    - Clé `F3_Nunes_YFV[3':CATTA] × FIP_ParidaD4[3':AAACC]` (rang `43/46`) : représentant `F3_Nunes_YFV#15_FIP_ParidaD4` (`0.02615681 nM` / `0.0262 nM`, `ΔG = -5.2681 kcal/mol`, `P_Boltzmann = 10.74 %`, `masked = True`) ; fusionne `1` variant(s) identique(s) sur les `5` nt 3' ne différant qu'en 5'/milieu (F3_Nunes_YFV#7_FIP_ParidaD4 à 0.0261 nM (P = 10.31 %)).
  - Couple parent `LF_Nunes_YFV × LF_ParidaD4` (`1` groupe(s) 3' distinct(s), `max P_Boltzmann` couple parent = `16.17 %`) :
    - Clé `LF_Nunes_YFV[3':GTATG] × LF_ParidaD4[3':GGCAG]` (rang `44/46`) : représentant `LF_Nunes_YFV_LF_ParidaD4` (`0.02526543 nM` / `0.0253 nM`, `ΔG = -3.5061 kcal/mol`, `P_Boltzmann = 16.17 %`, `masked = True`).

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

| Scénario | Mode | Appels Boltzmann | Groupes 3' éligibles | Groupes 3' non évalués | Couples parents éligibles | `pf_time` (`s`, médiane / étendue) | Total hors Boltzmann (`s`, médiane / étendue) | Total avec Boltzmann (`s`, médiane / étendue) | Part médiane `pf / avec` | Surcoût médian `pf / hors` |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 81 espèces — fenêtres locales (`test_step3_redesign_v2`) | défaut illimité, fenêtre 3' = 5 nt (`--boltzmann-max-dimers None`) | `77` | `46` | `0` | `40` | `0.456 (étendue 0.015, [0.455–0.470])` | `13.52 (étendue 2.25, [13.49–15.74])` | `14.04 (étendue 0.15, [13.96–14.11])` | `3.25 %` | `+3.37 %` |
| 81 espèces — fenêtres locales (`test_step3_redesign_v2`) | garde-fou K=25, fenêtre 3' = 5 nt (`--boltzmann-max-dimers 25`) | `56` | `46` | `21` | `40` | `0.392 (étendue 0.003, [0.391–0.394])` | `13.52 (étendue 2.25, [13.49–15.74])` | `13.99 (étendue 0.04, [13.96–14.00])` | `2.80 %` | `+2.90 %` |
| 81 espèces — `duplex_full` génomes viraux complets | défaut illimité, fenêtre 3' = 5 nt (`--boltzmann-max-dimers None`) | `77` | `46` | `0` | `40` | `0.452 (étendue 0.005, [0.451–0.456])` | `13.53 (étendue 0.10, [13.52–13.62])` | `14.03 (étendue 0.10, [13.99–14.09])` | `3.22 %` | `+3.34 %` |
| 81 espèces — `duplex_full` génomes viraux complets | garde-fou K=25, fenêtre 3' = 5 nt (`--boltzmann-max-dimers 25`) | `56` | `46` | `21` | `40` | `0.391 (étendue 0.011, [0.390–0.401])` | `13.53 (étendue 0.10, [13.52–13.62])` | `14.01 (étendue 0.14, [13.95–14.09])` | `2.79 %` | `+2.89 %` |
