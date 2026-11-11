# Vérité terrain paillasse — Validation RT-LAMP (LabCraft)

Ce répertoire regroupe les séquences d'amorces et les génomes complets de référence utilisés lors des essais expérimentaux en paillasse (kit Optigene ISO-004-RT, volume réactionnel de 25 µL, 63 °C).

## Fichiers et sources bibliographiques

### 1. Amorces RT-LAMP
- **`primer_dengue_3_panel3.txt`** : Panel 3 ciblant le virus de la dengue sérotype 3 (DENV-3), conçu avec Premier Biosoft LAMP Designer.
  - **Résultat paillasse** : **Succès en singleplex** (FIP/BIP 0,8 µM, LF/LB 0,4 µM, F3/B3 0,2 µM).
- **`parida_DENV4.txt`** : Panel RT-LAMP ciblant la région 3'-UTR du virus de la dengue sérotype 4 (DENV-4), saisi en séquences complètes commandées (devis TIB MOLBIOL 1000850, identiques à `validation/reference_data/parida2005_lamp_primers.csv` avec linker inclus pour `ParidaD4_FIP` et `ParidaD4_BIP`).
  - **Source** : Parida, M., Horioke, K., Ishida, H., Dash, P. K., Saxena, P., Jana, A. M., Islam, M. A., Inoue, S., Hosaka, N., & Morita, K. (2005). *Rapid detection and differentiation of dengue virus serotypes by a real-time reverse transcription-loop-mediated isothermal amplification assay.* **Journal of Clinical Microbiology**, 43(6), 2895–2903. https://doi.org/10.1128/JCM.43.6.2895-2903.2005
  - **Résultat paillasse** : **Succès en singleplex** (FIP/BIP 0,8 µM, LF/LB 0,4 µM, F3/B3 0,2 µM).
- **`nunes_YFV.txt`** : Panel RT-LAMP ciblant le gène NS1 du virus de la fièvre jaune (YFV), saisi en séquences complètes commandées (devis TIB MOLBIOL 1000850, linker `TTT` inclus pour `Nunes_YFV_FIP` et `Nunes_YFV_BIP`).
  - **Source** : Nunes, M. R. T., Viana, S., Teixeira, K. N., et al. (2015). *Reverse transcription loop-mediated isothermal amplification for yellow fever virus detection.* **Journal of Virological Methods**.
  - **Résultat paillasse** : **Succès en singleplex** (FIP/BIP 0,8 µM, LF/LB 0,4 µM, F3/B3 0,2 µM).
- **`duplex_YFV_DENV4_primers.txt`** : Réunion des panels `Nunes_YFV` et `ParidaD4` pour l'essai duplex fièvre jaune + dengue 4.
  - **Résultat paillasse** : **Échec** en demi-concentration par amorce (0,4 / 0,2 / 0,1 µM) comme en pleine concentration par amorce (0,8 / 0,4 / 0,2 µM), par perte du pic de fusion DENV-4 (89,5 °C) au profit de dimères parasites.

### 2. Génomes complets de référence (NCBI GenBank)
- **`DEN3_M93130.fasta`** : *Dengue virus type 3 strain H87, complete genome* — GenBank **`M93130.1`**.
- **`Dengue_4_S1028.fasta`** : *Dengue virus type 4 isolate D4/IND/PUNE/IRSHA-FG-01 (1028), complete genome* — GenBank **`MG272272.1`**.
- **`YFV_17D_17DD.fasta`** : *Yellow fever virus vaccine strain 17DD, complete genome* — GenBank **`U17066.1`**.
- **`duplex_YFV_DENV4_targets.fasta`** : Multi-FASTA réunissant **`U17066.1`** (YFV 17DD) et **`MG272272.1`** (DENV-4 S1028) pour l'analyse duplex conjointe.
