# Bancs d'essai reproductibles de référence (`reference_benchmarks`)

Document généré automatiquement par `python scripts/reference_benchmarks.py`.
Toute valeur numérique citée dans `CHANGELOG.md` doit provenir exclusivement de ce fichier.

## Métadonnées et conditions expérimentales

- **Commit de code évalué (`src/`)** : `81919b2` (`81919b2bdda4d2428002628c622a0370cf3d04f8`)
- **HEAD lors de l'exécution** : `81919b2`
- **Température de réaction (`analyze` & `coverage`)** : `63.0 °C`
- **Tampon LAMP par défaut (`analyze`)** : `50.0 mM Na+`, `10.0 mM K+`, `20.0 mM Tris`, `6.0 mM Mg2+`, `1.4 mM dNTP`
- **Tampon explicite du comparatif `coverage` (Section B)** : `50.0 mM Na+`, `0.0 mM K+`, `0.0 mM Tris`, `8.0 mM Mg2+`, `1.4 mM dNTP`, `--ddg-max 3.0 kcal/mol`
- **Profils enzymatiques (`dimer_dg_threshold`)** : `Bst = -2.0`, `Bst 2.0 = -2.0`, `Bst 2.0 WarmStart = -2.0`, `GspSSD = -2.0`, `Bst 3.0 = -1.5`, `Taq = -3.0` kcal/mol
- **Seuils de verdict et Boltzmann (`analyze`)** : `--dimer-conc-max = 0.230 nM`, `--dimer-warn-ratio = 0.500` (plancher `WARNING = 0.115 nM`), `--role-imbalance-ratio = 20.0`, `--boltzmann-max-dimers = 25` (plancher d'éligibilité `0.1 × 0.230 = 0.023 nM`)

### Décalage de `local_dg` (`37.0 °C, 1000 mM Na+` → `63.0 °C, tampon LAMP par défaut`)

| Motif 3' contigu (sans initiation) | `37.0 °C, 1 M Na+` (`kcal/mol`) | `63.0 °C, 1 M Na+` (`kcal/mol`) | `63.0 °C, tampon LAMP` (`kcal/mol`) | Décalage total (`kcal/mol`) |
|---|---:|---:|---:|---:|
| 3 pb (`5'-GGA-3'`) | `-3.14` | `-2.05` | `-1.63` | `+1.51` |
| 4 pb (`5'-GCAT-3'`) | `-4.56` | `-2.81` | `-2.19` | `+2.38` |
| 5 pb (`5'-GCATG-3'`) | `-6.02` | `-3.68` | `-2.85` | `+3.18` |
| 6 pb (`5'-GCATGC-3'`) | `-8.26` | `-5.28` | `-4.24` | `+4.02` |

---

## A. Cinq cas de référence de paillasse (`-2.0` vs `-3.5 kcal/mol`, génomes complets)

| Cas | Paillasse | Dimères amp. `-2.0` (`nM`) | dont inter `-2.0` (`nM`) | Dimères amp. `-3.5` (`nM`) | dont inter `-3.5` (`nM`) | Famille dominante (`-2.0`) | Verdict (`-2.0` / `-3.5`) | `suboptimal_3p_masked` (`cap=25`) | Couples parents éligibles | Couples parents non évalués (`cap=25`) |
|---|---|---:|---:|---:|---:|---|---|---:|---:|---:|
| `singleplex_denv3` | SUCCÈS paillasse | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `—` | `INFO` / `INFO` | `2` | `5` | `0` |
| `singleplex_denv4` | SUCCÈS paillasse | `0.1536` | `0.0000` | `0.1536` | `0.0000` | `BIP_ParidaD4 × LF_ParidaD4 (0.1536 nM)` | `WARNING` / `WARNING` | `0` | `13` | `0` |
| `singleplex_yfv` | SUCCÈS paillasse | `0.0011` | `0.0000` | `0.0000` | `0.0000` | `B3_Nunes_YFV × LF_Nunes_YFV (0.0011 nM)` | `WARNING` / `WARNING` | `0` | `6` | `0` |
| `duplex_half` | ÉCHEC paillasse (faux +) | `0.3527` | `0.3109` | `0.3363` | `0.2948` | `FIP_Nunes_YFV × F3_ParidaD4 (0.2931 nM)` | `FAILURE` / `FAILURE` | `0` | `19` | `0` |
| `duplex_full` | ÉCHEC paillasse (faux +) | `1.3547` | `1.2010` | `1.2938` | `1.1412` | `FIP_Nunes_YFV × F3_ParidaD4 (1.1351 nM)` | `FAILURE` / `FAILURE` | `0` | `40` | `15` |

### Causes globales (`verdict.global_cause` à `-2.0 kcal/mol`)

- `singleplex_denv3` (`INFO`) : Occupations d'initiation basses mais équilibrées entre amorces de même rôle (comportement attendu sur génome complet), sans dimère amplifiable significatif.
- `singleplex_denv4` (`WARNING`) : Dimères amplifiables modérés (0.154 nM >= plancher WARNING 0.115 nM, < seuil FAILURE 0.230 nM).
- `singleplex_yfv` (`WARNING`) : Occupation réduite par dégénérescence (B3_Nunes_YFV : 16 variants, ~12,5 nM par variant ; déséquilibre 55,8× face à F3_Nunes_YFV).
- `duplex_half` (`FAILURE`) : Présence de dimères d'amorces amplifiables majeurs ou séquestration critique interdisant l'amplification.
- `duplex_full` (`FAILURE`) : Présence de dimères d'amorces amplifiables majeurs ou séquestration critique interdisant l'amplification.

### Détail des dimères `suboptimal_3p_masked = True` (`--boltzmann-max-dimers 25`)

- `singleplex_denv3` : `2` dimère(s) `suboptimal_3p_masked` :
  - `LB_3#1_homo` : `conc = 0.1012 nM`, `P_Boltzmann = 42.76 %`, `ΔG = -5.24 kcal/mol`
  - `F3_3#2_BIP_3#1` : `conc = 0.0404 nM`, `P_Boltzmann = 20.53 %`, `ΔG = -5.09 kcal/mol`
- `singleplex_denv4` : `0` dimère `suboptimal_3p_masked`.
- `singleplex_yfv` : `0` dimère `suboptimal_3p_masked`.
- `duplex_half` : `0` dimère `suboptimal_3p_masked`.
- `duplex_full` : `0` dimère `suboptimal_3p_masked`.

### Réconciliation avec les `12` dimères `suboptimal_3p_masked` sans dédoublonnage ni plafond

- `singleplex_denv3` : sans dédoublonnage/plafond = `6` dimère(s) ; dédoublonné par couple parent sans plafond = `2` ; dédoublonné à `K = 25` = `2` (`5` couples parents éligibles).
  - Couple `LB_3 × LB_3` (rang `3/5`) : retenu via `LB_3#1_homo` (`0.1012 nM`, `P_Boltzmann = 42.76 %`) ; `1` autre(s) variant(s) du même couple écarté(s) par dédoublonnage (LB_3#1_LB_3#2 à 0.0722 nM).
  - Couple `BIP_3 × F3_3` (rang `4/5`) : retenu via `F3_3#2_BIP_3#1` (`0.0404 nM`, `P_Boltzmann = 20.53 %`) ; `3` autre(s) variant(s) du même couple écarté(s) par dédoublonnage (F3_3#2_BIP_3#2 à 0.0403 nM, F3_3#2_BIP_3#3 à 0.0398 nM, F3_3#2_BIP_3#4 à 0.0398 nM).
- `singleplex_denv4` : sans dédoublonnage/plafond = `0` dimère(s) ; dédoublonné par couple parent sans plafond = `0` ; dédoublonné à `K = 25` = `0` (`13` couples parents éligibles).
- `singleplex_yfv` : sans dédoublonnage/plafond = `1` dimère(s) ; dédoublonné par couple parent sans plafond = `0` ; dédoublonné à `K = 25` = `0` (`6` couples parents éligibles).
  - Couple `F3_Nunes_YFV × LB_Nunes_YFV` (rang `5/6`) : représentant le plus concentré = `F3_Nunes_YFV#5_LB_Nunes_YFV#2` (`0.04963523 nM`, `ΔG = -6.6209 kcal/mol`, `P_Boltzmann = 0.02 %`, `masked = False`), devant `F3_Nunes_YFV#6_LB_Nunes_YFV#2` (`0.04963376 nM`, `ΔG = -6.6209 kcal/mol`, `P_Boltzmann = 15.17 %`, `masked = True`, total `1` variant(s) masked dans ce couple).
- `duplex_half` : sans dédoublonnage/plafond = `0` dimère(s) ; dédoublonné par couple parent sans plafond = `0` ; dédoublonné à `K = 25` = `0` (`19` couples parents éligibles).
- `duplex_full` : sans dédoublonnage/plafond = `5` dimère(s) ; dédoublonné par couple parent sans plafond = `2` ; dédoublonné à `K = 25` = `0` (`40` couples parents éligibles).
  - Couple `F3_Nunes_YFV × LB_Nunes_YFV` (rang `27/40`) : représentant le plus concentré = `F3_Nunes_YFV#5_LB_Nunes_YFV#2` (`0.04936754 nM`, `ΔG = -6.6209 kcal/mol`, `P_Boltzmann = 0.02 %`, `masked = False`), devant `F3_Nunes_YFV#6_LB_Nunes_YFV#2` (`0.04936436 nM`, `ΔG = -6.6209 kcal/mol`, `P_Boltzmann = 15.17 %`, `masked = True`, total `1` variant(s) masked dans ce couple).
  - Couple `LB_ParidaD4 × LF_Nunes_YFV` (rang `33/40 > 25`) : représentant `LF_Nunes_YFV_LB_ParidaD4` (`0.0324 nM`, `P_Boltzmann = 10.27 %`, `masked = True`) hors plafond `K = 25`.
  - Couple `F3_Nunes_YFV × FIP_ParidaD4` (rang `39/40`) : représentant le plus concentré = `F3_Nunes_YFV#16_FIP_ParidaD4` (`0.02615694 nM`, `ΔG = -5.2681 kcal/mol`, `P_Boltzmann = 0.46 %`, `masked = False`), devant `F3_Nunes_YFV#15_FIP_ParidaD4` (`0.02615681 nM`, `ΔG = -5.2681 kcal/mol`, `P_Boltzmann = 10.74 %`, `masked = True`, total `2` variant(s) masked dans ce couple).
  - Couple `LF_Nunes_YFV × LF_ParidaD4` (rang `40/40 > 25`) : représentant `LF_Nunes_YFV_LF_ParidaD4` (`0.0253 nM`, `P_Boltzmann = 16.17 %`, `masked = True`) hors plafond `K = 25`.

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

| Scénario | Mode | Appels Boltzmann | Couples parents éligibles | Couples parents non évalués | `pf_time` (`s`, médiane / étendue) | Total hors Boltzmann (`s`, médiane / étendue) | Total avec Boltzmann (`s`, médiane / étendue) | Part médiane `pf / avec` | Surcoût médian `pf / hors` |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 81 espèces — fenêtres locales (`test_step3_redesign_v2`) | défaut parent-dédoublonné (`--boltzmann-max-dimers 25`) | `56` | `40` | `15` | `0.416 (étendue 0.014, [0.414–0.428])` | `16.64 (étendue 1.47, [15.69–17.16])` | `15.83 (étendue 1.12, [15.14–16.26])` | `2.63 %` | `+2.50 %` |
| 81 espèces — fenêtres locales (`test_step3_redesign_v2`) | tous couples parents (`--boltzmann-max-dimers 999999`) | `71` | `40` | `0` | `0.469 (étendue 0.135, [0.456–0.591])` | `16.64 (étendue 1.47, [15.69–17.16])` | `16.10 (étendue 2.27, [15.63–17.90])` | `2.91 %` | `+2.82 %` |
| 81 espèces — `duplex_full` génomes viraux complets | défaut parent-dédoublonné (`--boltzmann-max-dimers 25`) | `56` | `40` | `15` | `0.422 (étendue 0.016, [0.411–0.427])` | `15.33 (étendue 1.82, [14.14–15.96])` | `15.51 (étendue 1.24, [14.91–16.15])` | `2.72 %` | `+2.75 %` |
| 81 espèces — `duplex_full` génomes viraux complets | tous couples parents (`--boltzmann-max-dimers 999999`) | `71` | `40` | `0` | `0.463 (étendue 0.003, [0.461–0.465])` | `15.33 (étendue 1.82, [14.14–15.96])` | `15.62 (étendue 1.04, [14.85–15.90])` | `2.96 %` | `+3.02 %` |
