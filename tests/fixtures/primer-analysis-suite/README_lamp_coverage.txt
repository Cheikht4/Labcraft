---------------------------------------------------------
      README - Outil d'Analyse de Couverture LAMP
---------------------------------------------------------

DESCRIPTION
Le script `lamp_coverage.py` est un outil d'évaluation de la couverture théorique d'un jeu de séquences génomiques par un ou plusieurs sets d'amorces LAMP (Loop-Mediated Isothermal Amplification). 

Il simule l'hybridation des amorces en appliquant une double logique :
1. Une tolérance globale aux mésappariements (mismatches) paramétrable.
2. Une zone 3' strictement conservée sur la molécule physique (0 erreur permise) afin d'assurer l'élongation par la polymérase.

L'outil gère nativement le code IUPAC (ex: Y, R, N, W) pour les amorces dégénérées, et vérifie également la cohérence structurelle de l'amplicon final généré (ordre des amorces, présence, et calcul de la taille de l'amplicon).

---------------------------------------------------------
PRINCIPALES FONCTIONNALITÉS
---------------------------------------------------------
- Tolérance aux Mésappariements : Vous définissez le nombre max d'erreurs (par ex: 2 erreurs).
- Ancrage 3' Strict : Validation biologique vérifiant que l'extrémité 3' physique de chaque amorce matche à 100% sur un nombre défini de nucléotides (par défaut 3 nt).
- Filtrage Intelligent : Le script filtre automatiquement les séquences trop courtes (< 100 nt) présentes dans les alignements (MSA) pour ne pas fausser le pourcentage en testant les amorces sur elles-mêmes.
- Vérification Structurelle : Vérifie l'ordre d'hybridation des amorces selon l'ordre maître universel LAMP (inclut les configurations Loop et Stem) :
  [ F3 -> F2 -> FLOOP -> F1 -> STEMF -> STEMB -> B1 -> BLOOP -> B2 -> B3 ]
- Calcul de l'Amplicon : Estime la taille de l'amplicon réel (de l'extrémité 5' de F3 à l'extrémité 5' de B3 sur le brin cible).
- Format Exportable Excel : Exporte un rapport contenant les statistiques globales et une liste de séquences validées, séparées par des tabulations pour faciliter l'importation vers Excel ou Google Sheets.

---------------------------------------------------------
PRÉREQUIS
---------------------------------------------------------
Python 3.x
Modules nécessaires : biopython, regex
Installation : 
  pip install biopython regex

---------------------------------------------------------
FORMATS DES FICHIERS
---------------------------------------------------------
1. Fichier Cible (FASTA) :
   Génomes ou alignements (FASTA). Les tirets de gaps (-) sont automatiquement nettoyés lors de la lecture.

2. Fichier Amorces (FASTA ou TXT) :
   Il est impératif d'utiliser la nomenclature standardisée `Set_Amorce` pour le nom des séquences afin que le script groupe correctement les sets et vérifie la position de l'amorce.
   Exemple de fichier TXT :
   >1_F3
   GTGCGAAAYCCACTYTCAA
   >1_B3
   GTRGTRTCAGTCATRGCT
   >1_FLOOP
   CCTARGTCCACGTCTCTTT
   
   Noms d'amorces reconnus par le système d'ordre structurel : 
   F3, F2, FLOOP, F1, STEMF, STEMB, B1, BLOOP, B2, B3.

---------------------------------------------------------
UTILISATION (LIGNE DE COMMANDE)
---------------------------------------------------------
Commande de base :
  python3 lamp_coverage.py -t <fichier_cible.fasta> -p <fichier_amorces.txt> -o <rapport_sortie.txt> -e <max_erreurs> -s <strict_3_prime>

Options de base :
  -t, --target        (Obligatoire) Fichier FASTA cible contenant les séquences génomiques.
  -p, --primers       (Obligatoire) Fichier FASTA ou TXT des amorces.
  -o, --output        (Obligatoire) Fichier de rapport en sortie (texte).
  -e, --errors        (Optionnel) Nombre maximum d'erreurs tolérées par amorce (hors zone 3'). Défaut : 0.
  -s, --strict-3prime (Optionnel) Taille de la zone 3' stricte (0 erreur tolérée). Défaut : 3 nucléotides.
  --lng               (Optionnel) Langue du rapport généré (en pour anglais, fr pour français). Défaut : en.

Nouvelles Options de Sortie (Combinables) :
  --summary-only      (Option 2) N'affiche que les statistiques (masque l'affichage des séquences amplifiées dans le rapport principal pour l'alléger).
  --combine           (Option 3) Calcule la couverture combinatoire (multiplexage 2 à 2) des sets d'amorces et l'ajoute à la fin du rapport. Utile pour savoir si Set A + Set B améliorent la couverture globale.
  --export-seqs       (Option 4) Exporte la liste brute des noms de séquences validées dans des fichiers textes séparés pour chaque set (ex: `rapport_Set1_seqs.txt`).

---------------------------------------------------------
EXEMPLE D'EXÉCUTION
---------------------------------------------------------
python3 lamp_coverage.py -t dengue_aligned.fasta -p primers.txt -o rapport_dengue.txt -e 2 -s 3

Ceci testera la couverture sur `dengue_aligned.fasta` avec les amorces contenues dans `primers.txt`, en tolérant jusqu'à 2 erreurs par amorce, tout en gardant une contrainte forte (0 erreur) sur les 3 derniers nucléotides en 3' de chaque molécule.
Le résultat sera sauvegardé dans `rapport_dengue.txt`.

---------------------------------------------------------
COMPRÉHENSION DES RÉSULTATS
---------------------------------------------------------
Dans le rapport final généré, vous trouverez :
1. Matchs individuels : Le pourcentage de couverture théorique de chaque amorce prise indépendamment.
2. Match Global (Intersection Brute) : Le pourcentage des séquences du jeu de données hybrées par TOUTES les amorces du set.
3. Match Global Valide : Même chose, mais vérifiant que les amorces sont accrochées dans le "bon ordre" logique d'une réaction LAMP pour générer l'amplicon.
4. Tableau Détail : La liste complète des souches amplifiées, séparée par des tabulations, incluant la taille estimée de l'amplicon (pb) et son statut d'ordre structurel.
