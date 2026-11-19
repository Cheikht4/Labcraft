# Bancs d'essai reproductibles de référence (`reference_benchmarks`)

Document généré automatiquement par `python scripts/reference_benchmarks.py`.
Toute valeur numérique citée dans `CHANGELOG.md` doit provenir exclusivement de ce fichier.

## Métadonnées et conditions expérimentales

- **Commit de code évalué (`src/`)** : `e3b8d4c` (`e3b8d4cb3e142b491f49cb478a27bceaafdb46de`)
- **HEAD lors de l'exécution** : `e3b8d4c`
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

| Cas | Paillasse | Dimères amp. `-2.0` (`nM`) | dont inter `-2.0` (`nM`) | Dimères amp. `-3.5` (`nM`) | dont inter `-3.5` (`nM`) | Famille dominante (`-2.0`) | Verdict (`-2.0` / `-3.5`) | `suboptimal_3p_masked` (`cap=25`) | Éligibles non évalués (`cap=25`) |
|---|---|---:|---:|---:|---:|---|---|---:|---:|
| `singleplex_denv3` | SUCCÈS paillasse | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `—` | `INFO` / `INFO` | `6` | `0` |
| `singleplex_denv4` | SUCCÈS paillasse | `0.1536` | `0.0000` | `0.1536` | `0.0000` | `BIP_ParidaD4 × LF_ParidaD4 (0.1536 nM)` | `WARNING` / `WARNING` | `0` | `0` |
| `singleplex_yfv` | SUCCÈS paillasse | `0.0011` | `0.0000` | `0.0000` | `0.0000` | `B3_Nunes_YFV × LF_Nunes_YFV (0.0011 nM)` | `WARNING` / `WARNING` | `0` | `198` |
| `duplex_half` | ÉCHEC paillasse (faux +) | `0.3527` | `0.3109` | `0.3363` | `0.2948` | `FIP_Nunes_YFV × F3_ParidaD4 (0.2931 nM)` | `FAILURE` / `FAILURE` | `0` | `80` |
| `duplex_full` | ÉCHEC paillasse (faux +) | `1.3547` | `1.2010` | `1.2938` | `1.1412` | `FIP_Nunes_YFV × F3_ParidaD4 (1.1351 nM)` | `FAILURE` / `FAILURE` | `0` | `328` |

### Causes globales (`verdict.global_cause` à `-2.0 kcal/mol`)

- `singleplex_denv3` (`INFO`) : Occupations d'initiation basses mais équilibrées entre amorces de même rôle (comportement attendu sur génome complet), sans dimère amplifiable significatif.
- `singleplex_denv4` (`WARNING`) : Dimères amplifiables modérés (0.154 nM >= plancher WARNING 0.115 nM, < seuil FAILURE 0.230 nM).
- `singleplex_yfv` (`WARNING`) : Occupation réduite par dégénérescence (B3_Nunes_YFV : 16 variants, ~12,5 nM par variant ; déséquilibre 55,8× face à F3_Nunes_YFV).
- `duplex_half` (`FAILURE`) : Présence de dimères d'amorces amplifiables majeurs ou séquestration critique interdisant l'amplification.
- `duplex_full` (`FAILURE`) : Présence de dimères d'amorces amplifiables majeurs ou séquestration critique interdisant l'amplification.

### Détail des dimères `suboptimal_3p_masked = True` (`--boltzmann-max-dimers 25`)

- `singleplex_denv3` : `6` dimère(s) `suboptimal_3p_masked` :
  - `LB_3#1_homo` : `conc = 0.1012 nM`, `P_Boltzmann = 42.76 %`, `ΔG = -5.24 kcal/mol`
  - `LB_3#1_LB_3#2` : `conc = 0.0722 nM`, `P_Boltzmann = 43.20 %`, `ΔG = -5.01 kcal/mol`
  - `F3_3#2_BIP_3#1` : `conc = 0.0404 nM`, `P_Boltzmann = 20.53 %`, `ΔG = -5.09 kcal/mol`
  - `F3_3#2_BIP_3#2` : `conc = 0.0403 nM`, `P_Boltzmann = 17.37 %`, `ΔG = -5.09 kcal/mol`
  - `F3_3#2_BIP_3#3` : `conc = 0.0398 nM`, `P_Boltzmann = 22.85 %`, `ΔG = -5.09 kcal/mol`
  - `F3_3#2_BIP_3#4` : `conc = 0.0398 nM`, `P_Boltzmann = 19.59 %`, `ΔG = -5.09 kcal/mol`
- `singleplex_denv4` : `0` dimère `suboptimal_3p_masked`.
- `singleplex_yfv` : `0` dimère `suboptimal_3p_masked`.
- `duplex_half` : `0` dimère `suboptimal_3p_masked`.
- `duplex_full` : `0` dimère `suboptimal_3p_masked`.

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

## C. Performance `enrich_dimer_boltzmann_3p` sur 81 espèces (`3 321` paires)

| Scénario | Mode | Appels Boltzmann | Éligibles non évalués | `pf_time` (`s`) | Total hors Boltzmann (`s`) | Total avec Boltzmann (`s`) | Part `pf / avec` | Surcoût `pf / hors` |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 81 espèces — fenêtres locales (`test_step3_redesign_v2`) | défaut (`--boltzmann-max-dimers 25`) | `56` | `328` | `0.645` | `19.35` | `17.80` | `3.62 %` | `+3.33 %` |
| 81 espèces — fenêtres locales (`test_step3_redesign_v2`) | sans plafond (`999999`) | `384` | `0` | `5.848` | `16.82` | `22.90` | `25.54 %` | `+34.77 %` |
| 81 espèces — `duplex_full` génomes viraux complets | défaut (`--boltzmann-max-dimers 25`) | `56` | `328` | `0.655` | `17.13` | `18.58` | `3.53 %` | `+3.83 %` |
| 81 espèces — `duplex_full` génomes viraux complets | sans plafond (`999999`) | `384` | `0` | `5.634` | `18.07` | `22.26` | `25.31 %` | `+31.17 %` |
