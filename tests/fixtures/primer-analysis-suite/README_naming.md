# README — Convention de Nommage des Amorces / Primer Naming Convention

> Ce document décrit la convention de nommage obligatoire pour utiliser correctement
> le Primer Analysis Suite (lamp_coverage.py, primer_analyse.py).

---

## Format général / General format

```
>{SetID}_{PrimerType}
```

**Exemples :**
```fasta
>Dengue_2023_F3
ATGCATGCATGC...
>Dengue_2023_B3
GCATGCATGCAT...
```

Le programme découpe le nom sur le **dernier underscore** pour extraire :
- `SetID` = tout ce qui précède le dernier `_` (peut contenir des `_`)
- `PrimerType` = le type d'amorce (voir tableaux ci-dessous)

---

## Types d'amorces reconnus / Recognised primer types

### Mode LAMP (sans --pcr)

| Nom canonique | Variantes acceptées       | Description               |
|---|---|---|
| `F3`    | `FWD`, `FORWARD`         | Outer Forward             |
| `B3`    | `REV`, `REVERSE`         | Outer Backward            |
| `F2`    | —                        | Inner Forward             |
| `B2`    | —                        | Inner Backward            |
| `F1`    | `F1C`                    | Inner Forward complement  |
| `B1`    | `B1C`                    | Inner Backward complement |
| `FIP`   | —                        | Forward Inner Primer      |
| `BIP`   | —                        | Backward Inner Primer     |
| `FLOOP` | `LF`, `FLP`, `LOOPF`    | Forward Loop              |
| `BLOOP` | `LB`, `BLP`, `LOOPB`    | Backward Loop             |
| `STEMF` | —                        | Stem Forward (optionnel)  |
| `STEMB` | —                        | Stem Backward (optionnel) |

### Mode PCR (--pcr)

| Nom canonique | Variantes acceptées            | Description |
|---|---|---|
| `F`  | `FWD`, `FORWARD`, `F3`, `FP`  | Forward primer   |
| `R`  | `REV`, `REVERSE`, `B3`, `RP`  | Reverse primer   |
| `P`  | `PROBE`, `SONDE`, `Probe`     | Sonde TaqMan     |

---

## Amorces avec plusieurs versions / Multi-version primers

Lorsque plusieurs séquences candidates existent pour un même type,
utilisez le suffixe `_{N}` où N est un **entier positif** :

```
>{SetID}_{PrimerType}_{VersionNumber}
```

> **Règle stricte** : le suffixe de version doit être un entier pur (1, 2, 3...).
> Pas de lettre ajoutée : pas de `_1a`, `_v1`, `_alt`.

**Exemples LAMP :**
```fasta
>Ebola_2024_F3_1
ATGCATGCATGCATGCATGCATGCATGCA
>Ebola_2024_F3_2
ATGCATGCATGCATGGATGCATGCATGCA
>Ebola_2024_F3_3
ATGCATGCATGCATGCATGCATGGATGCA
>Ebola_2024_B3
GCATGCATGCATGCATGCATGCATGCATG
>Ebola_2024_FLOOP_1
TTTTTTTTTTGCGCGCGCGCGCGCGCG
>Ebola_2024_FLOOP_2
TTTTTTTTTTGCGCGCGCGCGCGAGCG
```

**Exemples PCR (sondes multiples) :**
```fasta
>CCHF_Ibrahim_FP
AAACAGGGGTGGTGATGAGA
>CCHF_Ibrahim_RP
GAACGGCCTGACTTGTTGAT
>CCHF_Ibrahim_Probe_1
TGAACATGTGGAGTGGTGTAGGGAATT
>CCHF_Ibrahim_Probe_2
CATGTGGACTGGTGCAGGGAGTT
```

---

## Ce que le programme calcule / What the program computes

Pour un set avec des amorces multi-versions, le rapport affiche :

```
Matchs individuels par amorce :
  - F31 : 82.00% (657/801)
  - F32 : 75.00% (601/801)
  - F33 : 68.00% (545/801)
  -> Union F3 (F31+F32+F33) : 95.00% (761/801)
  - B3  : 98.00% (785/801)

Match Global du Set (Intersection Brute) : 62.00% (497/801)
Couverture Validation (Union des types essentiels) : 93.10% (746/801)
Match Global Valide (Validation + Ordre Correct) : 92.50% (741/801)

Meilleure Combinaison (3 combinaisons testées, recherche exhaustive) :
  Combinaison optimale : B3=B3 + F3=F31
  Couverture de la combinaison optimale : 82.00% (657/801)
```

- **Couverture Union** = si toutes les versions sont utilisées dans la même réaction
- **Meilleure combinaison** = la combinaison de 1 seule version par type qui maximise la couverture

---

## Résumé des règles / Rules summary

| Règle                           | Correct ✅                   | Incorrect ❌               |
|---------------------------------|------------------------------|----------------------------|
| Format de base                  | `MySet_F3`                   | `MySet_F3primer`           |
| Multi-version avec entier pur   | `MySet_F3_1`                 | `MySet_F3_v1`, `MySet_F31` |
| SetID peut contenir des `_`     | `Dengue_2024_WestNile_F3`   | —                          |
| Pas d'espace dans le nom        | `MySet_F3`                   | `My Set_F3`                |
| Même SetID pour tous les membres| `SetA_F3`, `SetA_B3`         | `SetA_F3`, `SetB_B3`       |
