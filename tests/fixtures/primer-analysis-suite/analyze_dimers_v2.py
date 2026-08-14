#!/usr/bin/env python3
# =============================================================================
# analyze_dimers_v2.py - Thermodynamic analysis of primer dimers and hairpins
#                        Analyse thermodynamique des dimères et hairpins
# =============================================================================
# Author / Auteur  : Talibouya
# Version          : 2.2 (macOS – multilingual / multilingue)
# Description (EN) : Analyzes hairpins and dimerizations of primers via ntthal.
#                    Outputs .txt, .csv, and interactive .html reports.
# Description (FR) : Analyse les hairpins et dimérisations de primers via ntthal.
#                    Génère des rapports .txt, .csv et .html interactif.
# =============================================================================

import subprocess
import itertools
import os
import re
import csv
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

# --- Optional tqdm import / Importation optionnelle de tqdm ---
try:
    from tqdm import tqdm
except ImportError:
    # Fallback if tqdm is not installed / Fallback si tqdm n'est pas installé
    def tqdm(iterable, **kwargs):
        desc = kwargs.get('desc', '')
        if desc:
            print(f"{desc}...")
        for item in iterable:
            yield item

try:
    from Bio import SeqIO
except ImportError:
    # This message is shown before language selection, kept bilingual
    # Ce message s'affiche avant le choix de langue, donc bilingue
    print("❌ Error / Erreur : Biopython is not installed / n'est pas installé.")
    print("   Run / Lancez : pip3 install biopython")
    sys.exit(1)

# =============================================================================
# LANGUAGE DICTIONARY / DICTIONNAIRE DE LANGUES
# =============================================================================
# All user-facing strings are stored here.
# Toutes les chaînes affichées à l'utilisateur sont stockées ici.

LANG = {
    # ---- English ----
    'en': {
        # Language selection / Sélection de la langue
        'lang_prompt'         : "Select language / Choisissez la langue  [en / fr] (default: en) : ",
        'lang_invalid'        : "⚠  Unknown language. Defaulting to English.",

        # Header / En-tête
        'header_title'        : "  🧬 Thermodynamic Analysis of Primer Dimers",

        # ntthal detection / Détection de ntthal
        'detecting_ntthal'    : "\n🔍 Detecting ntthal...",
        'ntthal_found'        : "   ✅ ntthal found",
        'config_found'        : "   ✅ Primer3 config",
        'ntthal_missing'      : "❌ Error: 'ntthal' (Primer3) not found in your PATH.",
        'ntthal_install'      : "   Install it with: brew install primer3",
        'config_missing'      : "⚠️  primer3_config folder not found. Calculations may fail.",
        'config_auto'         : "(auto)",

        # FASTA input / Entrée FASTA
        'fasta_format'        : "\nExpected path format: /path/to/file.fasta",
        'fasta_prompt'        : "Enter the path to the FASTA file: ",
        'fasta_not_found'     : "❌ Error: File not found",
        'fasta_read_error'    : "❌ Error reading FASTA file",
        'fasta_empty'         : "❌ No sequences found in the FASTA file.",
        'primers_loaded'      : "✅ primer(s) loaded",

        # Sequence validation / Validation des séquences
        'invalid_bases_title' : "⚠️  Sequences with non-IUPAC characters:",
        'continue_prompt'     : "Continue anyway? (y/n): ",
        'continue_yes'        : ('y', 'yes'),

        # Analysis / Analyse
        'launching_analysis'  : "⚙️  Launching analysis",
        'hairpins_label'      : "hairpin(s)",
        'dimers_label'        : "dimer(s)",
        'progress_hairpin'    : "Hairpins",
        'progress_dimer'      : "Dimers  ",

        # Alert labels / Labels d'alerte
        'alert_critical'      : "🔴 PROBLEM",
        'alert_warning'       : "🟡 WARNING",
        'alert_ok'            : "🟢 OK",
        'alert_na'            : "❓ N/A",
        # Tm threshold descriptions / Descriptions des seuils Tm
        'tm_thresh_ok'        : "Tm < 0°C",
        'tm_thresh_warn'      : "20°C ≤ Tm ≤ 40°C",
        'tm_thresh_crit'      : "Tm > 50°C",

        # Result types / Types de résultats
        'type_hairpin'        : "Hairpin",
        'type_homo'           : "Homo-dimer",
        'type_hetero'         : "Hetero-dimer",

        # Summary / Résumé
        'summary_label'       : "📊 Summary",
        'summary_problems'    : "problem(s)",
        'summary_warnings'    : "warning(s)",

        # Output files / Fichiers de sortie
        'reports_generated'   : "✅ Reports generated:",
        'txt_label'           : "📄 Text ",
        'csv_label'           : "📊 CSV  ",
        'html_label'          : "🌐 HTML ",
        'html_hint'           : "Open the HTML report in your browser for interactive visualization.",
        'write_error'         : "❌ Error writing output files",

        # Text report / Rapport texte
        'txt_title'           : "Thermodynamic Analysis of Primers",
        'txt_source'          : "Source file",
        'txt_generated'       : "Generated on",
        'txt_summary_header'  : "📊 SUMMARY TABLE (sorted by descending Tm)",
        'txt_col_pair'        : "Pair",
        'txt_col_type'        : "Type",
        'txt_col_dg'          : "dG (cal/mol)",
        'txt_col_tm'          : "Tm (°C)",
        'txt_col_alert'       : "Alert",
        'txt_hairpin_section' : "I. HAIRPIN ANALYSIS (SELF-FOLDING)",
        'txt_dimer_section'   : "II. DIMERIZATION ANALYSIS (HETERO & HOMO-DIMERS)",
        'txt_primer'          : "Primer",
        'txt_sequence'        : "Sequence",
        'txt_pair'            : "Pair",
        'txt_alert'           : "Alert",
        'txt_structure'       : "Structure",

        # CSV headers / En-têtes CSV
        'csv_primer1'         : 'Primer1',
        'csv_primer2'         : 'Primer2',
        'csv_type'            : 'Type',
        'csv_dh'              : 'dH (cal/mol)',
        'csv_ds'              : 'dS (cal/K/mol)',
        'csv_dg'              : 'dG (cal/mol)',
        'csv_tm'              : 'Tm (°C)',
        'csv_alert'           : 'Alert',
        'csv_structure'       : 'Structure',

        # HTML report / Rapport HTML
        'html_lang'           : 'en',
        'html_page_title'     : "Dimer Analysis",
        'html_h1'             : "🧬 Thermodynamic Analysis of Primers",
        'html_source'         : "Source file",
        'html_generated'      : "Generated on",
        'html_legend_ok'      : "🟢 OK — Tm < 0°C",
        'html_legend_warn'    : "🟡 WARNING — 20°C ≤ Tm ≤ 40°C",
        'html_legend_crit'    : "🔴 PROBLEM — Tm > 50°C",
        'html_legend_lt'      : "",
        'html_legend_between' : "",
        'html_stat_total'     : "Total analyses",
        'html_stat_crit'      : "Critical problems",
        'html_stat_warn'      : "Warnings",
        'html_stat_ok'        : "Results OK",
        'html_h2_hairpin'     : "I. Hairpins (Self-folding)",
        'html_h2_dimer'       : "II. Dimerizations (Homo & Hetero-dimers)",
        'html_th_p1'          : "Primer 1",
        'html_th_p2'          : "Primer 2",
        'html_th_type'        : "Type",
        'html_th_dh'          : "dH (cal/mol)",
        'html_th_ds'          : "dS (cal/K/mol)",
        'html_th_dg'          : "dG (cal/mol)",
        'html_th_tm'          : "Tm (°C)",
        'html_th_alert'       : "Alert",
    },

    # ---- Français ----
    'fr': {
        # Language selection / Sélection de la langue
        'lang_prompt'         : "Select language / Choisissez la langue  [en / fr] (défaut: en) : ",
        'lang_invalid'        : "⚠  Langue inconnue. Anglais utilisé par défaut.",

        # Header / En-tête
        'header_title'        : "  🧬 Analyse Thermodynamique des Dimères de Primers",

        # ntthal detection / Détection de ntthal
        'detecting_ntthal'    : "\n🔍 Détection de ntthal...",
        'ntthal_found'        : "   ✅ ntthal trouvé",
        'config_found'        : "   ✅ Config Primer3",
        'ntthal_missing'      : "❌ Erreur : 'ntthal' (Primer3) est introuvable dans votre PATH.",
        'ntthal_install'      : "   Installez-le avec : brew install primer3",
        'config_missing'      : "⚠️  Dossier primer3_config introuvable. Les calculs peuvent échouer.",
        'config_auto'         : "(auto)",

        # FASTA input / Entrée FASTA
        'fasta_format'        : "\nFormat du chemin attendu : /chemin/vers/fichier.fasta",
        'fasta_prompt'        : "Entrez le chemin du fichier FASTA : ",
        'fasta_not_found'     : "❌ Erreur : Fichier introuvable",
        'fasta_read_error'    : "❌ Erreur de lecture du fichier FASTA",
        'fasta_empty'         : "❌ Aucune séquence trouvée dans le fichier FASTA.",
        'primers_loaded'      : "✅ primer(s) chargé(s)",

        # Sequence validation / Validation des séquences
        'invalid_bases_title' : "⚠️  Séquences contenant des caractères non-IUPAC :",
        'continue_prompt'     : "Continuer quand même ? (o/n) : ",
        'continue_yes'        : ('o', 'oui', 'y', 'yes'),

        # Analysis / Analyse
        'launching_analysis'  : "⚙️  Lancement de l'analyse",
        'hairpins_label'      : "hairpin(s)",
        'dimers_label'        : "dimère(s)",
        'progress_hairpin'    : "Hairpins",
        'progress_dimer'      : "Dimères  ",

        # Alert labels / Labels d'alerte
        'alert_critical'      : "🔴 PROBLÈME",
        'alert_warning'       : "🟡 ATTENTION",
        'alert_ok'            : "🟢 OK",
        'alert_na'            : "❓ N/A",
        # Tm threshold descriptions / Descriptions des seuils Tm
        'tm_thresh_ok'        : "Tm < 0°C",
        'tm_thresh_warn'      : "20°C ≤ Tm ≤ 40°C",
        'tm_thresh_crit'      : "Tm > 50°C",

        # Result types / Types de résultats
        'type_hairpin'        : "Hairpin",
        'type_homo'           : "Homo-dimère",
        'type_hetero'         : "Hétéro-dimère",

        # Summary / Résumé
        'summary_label'       : "📊 Résumé",
        'summary_problems'    : "problème(s)",
        'summary_warnings'    : "attention(s)",

        # Output files / Fichiers de sortie
        'reports_generated'   : "✅ Rapports générés :",
        'txt_label'           : "📄 Texte",
        'csv_label'           : "📊 CSV  ",
        'html_label'          : "🌐 HTML ",
        'write_error'         : "❌ Erreur lors de l'écriture des fichiers",
        'html_hint'           : "Ouvrez le rapport HTML dans votre navigateur pour une visualisation interactive.",

        # Text report / Rapport texte
        'txt_title'           : "Analyse thermodynamique des primers",
        'txt_source'          : "Fichier source",
        'txt_generated'       : "Généré le",
        'txt_summary_header'  : "📊 TABLEAU RÉCAPITULATIF (trié par Tm décroissant)",
        'txt_col_pair'        : "Paire",
        'txt_col_type'        : "Type",
        'txt_col_dg'          : "dG (cal/mol)",
        'txt_col_tm'          : "Tm (°C)",
        'txt_col_alert'       : "Alerte",
        'txt_hairpin_section' : "I. ANALYSE DES HAIRPINS (AUTO-REPLIEMENT)",
        'txt_dimer_section'   : "II. ANALYSE DES DIMÉRISATIONS (HÉTÉRO ET HOMO-DIMÈRES)",
        'txt_primer'          : "Primer",
        'txt_sequence'        : "Séquence",
        'txt_pair'            : "Paire",
        'txt_alert'           : "Alerte",
        'txt_structure'       : "Structure",

        # CSV headers / En-têtes CSV
        'csv_primer1'         : 'Primer1',
        'csv_primer2'         : 'Primer2',
        'csv_type'            : 'Type',
        'csv_dh'              : 'dH (cal/mol)',
        'csv_ds'              : 'dS (cal/K/mol)',
        'csv_dg'              : 'dG (cal/mol)',
        'csv_tm'              : 'Tm (°C)',
        'csv_alert'           : 'Alerte',
        'csv_structure'       : 'Structure',

        # HTML report / Rapport HTML
        'html_lang'           : 'fr',
        'html_page_title'     : "Analyse Dimères",
        'html_h1'             : "🧬 Analyse Thermodynamique des Primers",
        'html_source'         : "Fichier source",
        'html_generated'      : "Généré le",
        'html_legend_ok'      : "🟢 OK — Tm < 0°C",
        'html_legend_warn'    : "🟡 ATTENTION — 20°C ≤ Tm ≤ 40°C",
        'html_legend_crit'    : "🔴 PROBLÈME — Tm > 50°C",
        'html_legend_lt'      : "",
        'html_legend_between' : "",
        'html_stat_total'     : "Analyses totales",
        'html_stat_crit'      : "Problèmes critiques",
        'html_stat_warn'      : "Avertissements",
        'html_stat_ok'        : "Résultats OK",
        'html_h2_hairpin'     : "I. Hairpins (Auto-repliement)",
        'html_h2_dimer'       : "II. Dimérisations (Homo & Hétéro-dimères)",
        'html_th_p1'          : "Primer 1",
        'html_th_p2'          : "Primer 2",
        'html_th_type'        : "Type",
        'html_th_dh'          : "dH (cal/mol)",
        'html_th_ds'          : "dS (cal/K/mol)",
        'html_th_dg'          : "dG (cal/mol)",
        'html_th_tm'          : "Tm (°C)",
        'html_th_alert'       : "Alerte",
    },
}

# =============================================================================
# CONSTANTS / CONSTANTES
# =============================================================================
# Valid IUPAC characters for a nucleic acid sequence
# Caractères IUPAC valides pour une séquence d'acide nucléique
VALID_BASES = set("ATGCNRYSWKMBDHVatgcnryswkmbdhv")

# Thermodynamic alert thresholds based on melting temperature Tm (°C)
# Seuils d'alerte thermodynamique basés sur la température de fusion Tm (°C)
# Tm < 0      → OK      (structure trop instable pour être problématique)
# 20 ≤ Tm ≤ 40 → WARNING (structure modérément stable, à surveiller)
# Tm > 50     → PROBLEM  (structure très stable, risque d'inhibition)
TM_OK_MAX      =  0.0   # Below this Tm → OK / En dessous → OK
TM_WARNING_MIN = 20.0   # Start of warning zone / Début zone attention
TM_WARNING_MAX = 40.0   # End of warning zone / Fin zone attention
TM_CRITICAL    = 50.0   # Above this Tm → PROBLEM / Au-dessus → PROBLÈME

# =============================================================================
# UTILITY FUNCTIONS / FONCTIONS UTILITAIRES
# =============================================================================

def select_language() -> dict:
    """
    Prompts user to select language and returns the corresponding translation dict.
    Demande à l'utilisateur de choisir la langue et retourne le dictionnaire correspondant.
    Default / Défaut : English (en)
    """
    # This prompt is always bilingual since language isn't chosen yet
    # Cette invite est toujours bilingue car la langue n'est pas encore choisie
    choice = input(LANG['en']['lang_prompt']).strip().lower()

    if choice == 'fr':
        return LANG['fr']
    elif choice in ('en', ''):
        return LANG['en']
    else:
        T = LANG['en']
        print(T['lang_invalid'])
        return T


def detect_ntthal(T: dict) -> tuple:
    """
    Automatically detects ntthal binary and primer3 config path.
    Détecte automatiquement le binaire ntthal et le chemin de configuration primer3.
    Returns (ntthal_path, config_path) or exits if not found.
    Retourne (ntthal_path, config_path) ou quitte si introuvable.
    """
    ntthal_path = shutil.which("ntthal")
    if not ntthal_path:
        print(T['ntthal_missing'])
        print(T['ntthal_install'])
        sys.exit(1)

    # Find primer3_config relative to binary / Trouver primer3_config relatif au binaire
    bin_dir = os.path.dirname(os.path.realpath(ntthal_path))
    base    = os.path.dirname(bin_dir)

    search_roots = [
        os.path.join(base, "share", "primer3", "primer3_config"),
        os.path.join(base, "Cellar"),
    ]

    config_path = None
    if os.path.isdir(search_roots[0]):
        config_path = search_roots[0]
    else:
        # Walk Cellar to find primer3_config / Parcourir Cellar pour trouver primer3_config
        for root, dirs, _ in os.walk(search_roots[1]):
            if "primer3_config" in dirs:
                config_path = os.path.join(root, "primer3_config")
                break

    if not config_path:
        print(T['config_missing'])
        config_path = ""

    return ntthal_path, config_path


def validate_sequences(primers: dict) -> list:
    """
    Validates that each sequence contains only valid IUPAC bases.
    Valide que chaque séquence ne contient que des bases IUPAC valides.
    Returns a list of error strings / Retourne la liste des erreurs.
    """
    errors = []
    for name, seq in primers.items():
        invalid = set(seq) - VALID_BASES
        if invalid:
            errors.append(f"  - {name} : {''.join(sorted(invalid))}")
    return errors


def get_alert_label(tm_str: str, T: dict) -> str:
    """
    Returns a localized alert label based on the Tm value (melting temperature in °C).
    Retourne un label d'alerte localisé selon la valeur de Tm (température de fusion en °C).

    Classification:
        Tm < 0°C         → OK      (structure trop instable / too unstable)
        20°C ≤ Tm ≤ 40°C → WARNING (structure modérément stable / moderately stable)
        Tm > 50°C        → PROBLEM  (structure très stable / very stable, inhibition risk)
        0°C ≤ Tm < 20°C  → OK      (zone intermédiaire basse / low intermediate zone)
        40°C < Tm ≤ 50°C → WARNING (zone intermédiaire haute / high intermediate zone)
    """
    try:
        tm = float(tm_str)
    except (ValueError, TypeError):
        return T['alert_na']
    if tm > TM_CRITICAL:
        return T['alert_critical']
    elif tm >= TM_WARNING_MIN:
        return T['alert_warning']
    return T['alert_ok']


def read_primers(fasta_file: str) -> dict:
    """
    Reads primer sequences from a FASTA file.
    Lit les séquences d'amorces depuis un fichier FASTA.
    """
    primers = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        primers[record.id] = str(record.seq).upper()
    return primers


# =============================================================================
# THERMODYNAMIC ANALYSIS / ANALYSE THERMODYNAMIQUE
# =============================================================================

def analyze_dimer(seq1: str, seq2: str, name1: str, name2: str,
                  mode: str = 'ANY',
                  ntthal_path: str = 'ntthal',
                  config_path: str = '',
                  T: dict = None) -> dict:
    """
    Runs ntthal to compute thermodynamic parameters for a primer pair.
    Lance ntthal pour calculer les paramètres thermodynamiques d'une paire.

    Parameters / Paramètres:
        seq1, seq2   : Primer sequences / Séquences des amorces
        name1, name2 : Identifiers / Identifiants
        mode         : 'ANY' (dimer) or 'HAIRPIN' (self-folding)
        ntthal_path  : Path to ntthal binary / Chemin vers le binaire ntthal
        config_path  : Path to primer3_config / Chemin vers primer3_config
        T            : Translation dict / Dictionnaire de traduction
    """
    if T is None:
        T = LANG['en']

    path_arg = f' -path "{config_path}/"' if config_path else ''

    if mode == 'HAIRPIN':
        cmd = f'"{ntthal_path}" -s1 {seq1} -a HAIRPIN -t 37{path_arg}'
    else:
        cmd = f'"{ntthal_path}" -s1 {seq1} -s2 {seq2} -a ANY -t 37{path_arg}'

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    # Parse output / Analyser la sortie
    output_lines = result.stdout.split('\n')
    dimer_info   = {}
    structure    = []

    for line in output_lines:
        if "Calculated thermodynamical parameters" in line:
            # Extract values via regex / Extraction des valeurs via regex
            for key, pattern in [('dS', r'dS = ([-\d.]+)'),
                                  ('dH', r'dH = ([-\d.]+)'),
                                  ('dG', r'dG = ([-\d.]+)'),
                                  ('Tm', r't = ([-\d.]+)')]:
                m = re.search(pattern, line)
                if m:
                    dimer_info[key] = m.group(1)
        elif line.strip().startswith(('SEQ', 'STR')):
            # Preserve structure spacing / Conserver les espaces de la structure
            structure.append(line.rstrip())

    dg_val = dimer_info.get('dG', 'N/A')
    tm_val = dimer_info.get('Tm', 'N/A')

    # Determine type label / Déterminer le label de type
    if mode == 'HAIRPIN':
        type_label = T['type_hairpin']
    elif name1 == name2:
        type_label = T['type_homo']
    else:
        type_label = T['type_hetero']

    return {
        'primer1'  : name1,
        'primer2'  : name2,
        'mode'     : mode,
        'type'     : type_label,
        'sequence1': seq1,
        'sequence2': seq2,
        'dH'       : dimer_info.get('dH', 'N/A'),
        'dS'       : dimer_info.get('dS', 'N/A'),
        'dG'       : dg_val,
        'Tm'       : tm_val,
        # Alert is now based on Tm / L'alerte est maintenant basée sur la Tm
        'alert'    : get_alert_label(tm_val, T),
        'structure': '\n'.join(structure) if structure else 'N/A',
    }


def run_analyses_parallel(tasks: list, ntthal_path: str, config_path: str,
                           T: dict, description: str = "Analysis") -> list:
    """
    Runs analyses in parallel to reduce computation time.
    Exécute les analyses en parallèle pour réduire le temps de calcul.
    tasks: list of (seq1, seq2, name1, name2, mode) tuples.
    """
    results = [None] * len(tasks)

    def worker(idx_and_task):
        i, (seq1, seq2, name1, name2, mode) = idx_and_task
        return i, analyze_dimer(seq1, seq2, name1, name2, mode, ntthal_path, config_path, T)

    with ThreadPoolExecutor() as executor:
        futures = {executor.submit(worker, (i, task)): i for i, task in enumerate(tasks)}
        for future in tqdm(as_completed(futures), total=len(futures), desc=description):
            i, result = future.result()
            results[i] = result

    return results


# =============================================================================
# RESULT EXPORT / EXPORTATION DES RÉSULTATS
# =============================================================================

def write_text_report(hairpin_results: list, dimer_results: list,
                      output_file: str, fasta_name: str, T: dict):
    """
    Writes the text report with alert thresholds and summary table.
    Écrit le rapport texte avec seuils d'alerte et tableau récapitulatif.
    """
    all_results = hairpin_results + dimer_results

    # Sort by descending Tm (most critical first) / Trier par Tm décroissant (les plus critiques en premier)
    def sort_key(r):
        try:
            return float(r['Tm'])
        except (ValueError, TypeError):
            return -999.0
    sorted_results = sorted(all_results, key=sort_key, reverse=True)

    with open(output_file, "w", encoding='utf-8') as f:
        f.write(f"{T['txt_title']} ({T['txt_source']}: {fasta_name})\n")
        f.write(f"{T['txt_generated']}: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 70 + "\n\n")

        # ---- Summary table / Tableau récapitulatif ----
        f.write(f"{T['txt_summary_header']}\n")
        f.write("-" * 70 + "\n")
        col_pair  = T['txt_col_pair']
        col_type  = T['txt_col_type']
        col_dg    = T['txt_col_dg']
        col_tm    = T['txt_col_tm']
        col_alert = T['txt_col_alert']
        f.write(f"{col_pair:<35} {col_type:<16} {col_dg:<14} {col_tm:<10} {col_alert}\n")
        f.write("-" * 70 + "\n")
        for r in sorted_results:
            pair = f"{r['primer1']} / {r['primer2']}"
            f.write(f"{pair:<35} {r['type']:<16} {r['dG']:<14} {r['Tm']:<10} {r['alert']}\n")
        f.write("=" * 70 + "\n\n")

        # ---- Hairpin section / Section hairpins ----
        f.write(f"{T['txt_hairpin_section']}\n")
        f.write("-" * 50 + "\n")
        for r in hairpin_results:
            f.write(f"{T['txt_primer']}: {r['primer1']}\n")
            f.write(f"{T['txt_sequence']}: {r['sequence1']}\n")
            f.write(f"dH: {r['dH']} cal/mol | dS: {r['dS']} cal/K/mol | dG: {r['dG']} cal/mol | Tm: {r['Tm']} °C\n")
            f.write(f"{T['txt_alert']}: {r['alert']}\n")
            f.write(f"{T['txt_structure']}:\n{r['structure']}\n")
            f.write("-" * 30 + "\n")
        f.write("\n" + "=" * 70 + "\n\n")

        # ---- Dimer section / Section dimères ----
        f.write(f"{T['txt_dimer_section']}\n")
        f.write("-" * 50 + "\n")
        for r in dimer_results:
            f.write(f"{T['txt_pair']}: {r['primer1']} - {r['primer2']} ({r['type']})\n")
            f.write(f"dH: {r['dH']} cal/mol | dS: {r['dS']} cal/K/mol | dG: {r['dG']} cal/mol | Tm: {r['Tm']} °C\n")
            f.write(f"{T['txt_alert']}: {r['alert']}\n")
            f.write(f"{T['txt_structure']}:\n{r['structure']}\n")
            f.write("-" * 30 + "\n")


def write_csv_report(hairpin_results: list, dimer_results: list,
                     output_file: str, T: dict):
    """
    Exports results to a CSV file compatible with Excel/R/pandas.
    Exporte les résultats dans un fichier CSV compatible Excel/R/pandas.
    """
    fieldnames = [T['csv_primer1'], T['csv_primer2'], T['csv_type'],
                  T['csv_dh'], T['csv_ds'], T['csv_dg'],
                  T['csv_tm'], T['csv_alert'], T['csv_structure']]

    with open(output_file, "w", newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in hairpin_results + dimer_results:
            writer.writerow({
                T['csv_primer1']  : r['primer1'],
                T['csv_primer2']  : r['primer2'],
                T['csv_type']     : r['type'],
                T['csv_dh']       : r['dH'],
                T['csv_ds']       : r['dS'],
                T['csv_dg']       : r['dG'],
                T['csv_tm']       : r['Tm'],
                T['csv_alert']    : r['alert'],
                T['csv_structure']: r['structure'].replace('\n', ' | '),
            })


def write_html_report(hairpin_results: list, dimer_results: list,
                      output_file: str, fasta_name: str, T: dict):
    """
    Generates an interactive HTML report with sorting and color-coding by alert level.
    Génère un rapport HTML interactif avec tri et code couleur selon les alertes.
    """
    all_results = hairpin_results + dimer_results

    def tm_class(tm_str):
        # CSS class based on Tm alert level / Classe CSS selon le niveau d'alerte Tm
        try:
            tm = float(tm_str)
        except (ValueError, TypeError):
            return "na"
        if tm > TM_CRITICAL:
            return "critical"
        elif tm >= TM_WARNING_MIN:
            return "warning"
        return "ok"

    def make_rows(results):
        html = ""
        for r in results:
            cls = tm_class(r['Tm'])
            # Escape < > in structure / Échapper les < > dans la structure
            struct = r['structure'].replace('<', '&lt;').replace('>', '&gt;').replace('\n', '<br>')
            html += f"""
            <tr class="{cls}">
                <td>{r['primer1']}</td>
                <td>{r['primer2']}</td>
                <td>{r['type']}</td>
                <td>{r['dH']}</td>
                <td>{r['dS']}</td>
                <td><strong>{r['dG']}</strong></td>
                <td>{r['Tm']}</td>
                <td class="alert-cell">{r['alert']}</td>
            </tr>
            <tr class="struct-row {cls}">
                <td colspan="8"><pre class="structure">{struct}</pre></td>
            </tr>"""
        return html

    hairpin_rows = make_rows(hairpin_results)
    dimer_rows   = make_rows(dimer_results)
    generated_at = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    n_crit  = sum(1 for r in all_results if r['alert'].startswith('🔴'))
    n_warn  = sum(1 for r in all_results if r['alert'].startswith('🟡'))
    n_ok    = sum(1 for r in all_results if r['alert'].startswith('🟢'))

    # Localized HTML strings / Chaînes HTML localisées
    hl    = T['html_lang']
    h1    = T['html_h1']
    src   = T['html_source']
    gen   = T['html_generated']
    # Legend now shows Tm thresholds / La légende affiche maintenant les seuils Tm
    leg_ok   = T['html_legend_ok']
    leg_warn = T['html_legend_warn']
    leg_crit = T['html_legend_crit']

    html = f"""<!DOCTYPE html>
<html lang="{hl}">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{T['html_page_title']} – {fasta_name}</title>
<style>
  :root {{
    --bg: #0f1117; --surface: #1a1d27; --surface2: #22263a;
    --border: #2e3350; --text: #e2e8f0; --muted: #8892a4;
    --accent: #6366f1; --ok: #22c55e; --warn: #f59e0b; --crit: #ef4444;
    --ok-bg: rgba(34,197,94,0.07); --warn-bg: rgba(245,158,11,0.10);
    --crit-bg: rgba(239,68,68,0.12);
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ background: var(--bg); color: var(--text); font-family: 'Segoe UI', system-ui, sans-serif; padding: 2rem; }}
  h1 {{ font-size: 1.8rem; color: var(--accent); margin-bottom: .25rem; }}
  .subtitle {{ color: var(--muted); font-size: .9rem; margin-bottom: 2rem; }}
  h2 {{ font-size: 1.2rem; margin: 2rem 0 .75rem; color: var(--text); border-left: 3px solid var(--accent); padding-left: .75rem; }}
  .legend {{ display: flex; gap: 1.5rem; margin-bottom: 1.5rem; flex-wrap: wrap; }}
  .legend-item {{ display: flex; align-items: center; gap: .4rem; font-size: .85rem; }}
  .legend-dot {{ width: 12px; height: 12px; border-radius: 50%; }}
  .dot-ok {{ background: var(--ok); }} .dot-warn {{ background: var(--warn); }} .dot-crit {{ background: var(--crit); }}
  .table-wrap {{ overflow-x: auto; border-radius: 10px; border: 1px solid var(--border); margin-bottom: 2rem; }}
  table {{ width: 100%; border-collapse: collapse; font-size: .875rem; }}
  thead tr {{ background: var(--surface2); }}
  thead th {{
    padding: .7rem 1rem; text-align: left; color: var(--muted);
    font-weight: 600; cursor: pointer; user-select: none; white-space: nowrap;
    border-bottom: 1px solid var(--border);
  }}
  thead th:hover {{ color: var(--accent); }}
  thead th::after {{ content: ' ↕'; font-size: .7rem; opacity: .5; }}
  tbody tr {{ border-bottom: 1px solid var(--border); transition: background .15s; }}
  tbody tr:hover {{ background: var(--surface2); }}
  tbody td {{ padding: .6rem 1rem; vertical-align: top; }}
  tr.ok       {{ background: var(--ok-bg); }}
  tr.warning  {{ background: var(--warn-bg); }}
  tr.critical {{ background: var(--crit-bg); }}
  .alert-cell {{ font-weight: 600; white-space: nowrap; }}
  tr.struct-row td {{ padding: .25rem 1rem .6rem; }}
  pre.structure {{
    font-family: 'Courier New', monospace; font-size: .8rem;
    color: var(--muted); background: var(--surface); border-radius: 6px;
    padding: .5rem .75rem; white-space: pre; overflow-x: auto;
  }}
  .stats {{ display: flex; gap: 1rem; flex-wrap: wrap; margin-bottom: 2rem; }}
  .stat-card {{
    flex: 1; min-width: 140px; background: var(--surface); border: 1px solid var(--border);
    border-radius: 10px; padding: 1rem 1.25rem;
  }}
  .stat-num {{ font-size: 2rem; font-weight: 700; }}
  .stat-label {{ color: var(--muted); font-size: .8rem; margin-top: .15rem; }}
  .num-crit {{ color: var(--crit); }} .num-warn {{ color: var(--warn); }} .num-ok {{ color: var(--ok); }}
</style>
</head>
<body>
<h1>{h1}</h1>
<p class="subtitle">{src} : <strong>{fasta_name}</strong> &nbsp;·&nbsp; {gen} {generated_at}</p>

<div class="legend">
  <div class="legend-item"><div class="legend-dot dot-ok"></div> {leg_ok}</div>
  <div class="legend-item"><div class="legend-dot dot-warn"></div> {leg_warn}</div>
  <div class="legend-item"><div class="legend-dot dot-crit"></div> {leg_crit}</div>
</div>

<div class="stats">
  <div class="stat-card"><div class="stat-num">{len(all_results)}</div><div class="stat-label">{T['html_stat_total']}</div></div>
  <div class="stat-card"><div class="stat-num num-crit">{n_crit}</div><div class="stat-label">{T['html_stat_crit']}</div></div>
  <div class="stat-card"><div class="stat-num num-warn">{n_warn}</div><div class="stat-label">{T['html_stat_warn']}</div></div>
  <div class="stat-card"><div class="stat-num num-ok">{n_ok}</div><div class="stat-label">{T['html_stat_ok']}</div></div>
</div>

<h2>{T['html_h2_hairpin']}</h2>
<div class="table-wrap">
<table id="table-hairpin">
  <thead><tr>
    <th>{T['html_th_p1']}</th><th>{T['html_th_p2']}</th><th>{T['html_th_type']}</th>
    <th>{T['html_th_dh']}</th><th>{T['html_th_ds']}</th><th>{T['html_th_dg']}</th>
    <th>{T['html_th_tm']}</th><th>{T['html_th_alert']}</th>
  </tr></thead>
  <tbody>{hairpin_rows}</tbody>
</table>
</div>

<h2>{T['html_h2_dimer']}</h2>
<div class="table-wrap">
<table id="table-dimers">
  <thead><tr>
    <th>{T['html_th_p1']}</th><th>{T['html_th_p2']}</th><th>{T['html_th_type']}</th>
    <th>{T['html_th_dh']}</th><th>{T['html_th_ds']}</th><th>{T['html_th_dg']}</th>
    <th>{T['html_th_tm']}</th><th>{T['html_th_alert']}</th>
  </tr></thead>
  <tbody>{dimer_rows}</tbody>
</table>
</div>

<script>
// Table sorting / Tri de tableau
document.querySelectorAll('table').forEach(table => {{
  table.querySelectorAll('thead th').forEach((th, idx) => {{
    th.addEventListener('click', () => {{
      const tbody = table.querySelector('tbody');
      // Group data rows + structure rows / Regrouper rangées données + structure
      const dataRows = [...tbody.querySelectorAll('tr:not(.struct-row)')];
      const pairs = dataRows.map(dr => [dr, dr.nextElementSibling]);
      const asc = th.dataset.asc !== 'true';
      th.dataset.asc = asc;
      pairs.sort((a, b) => {{
        const va = a[0].cells[idx]?.textContent.trim() || '';
        const vb = b[0].cells[idx]?.textContent.trim() || '';
        const na = parseFloat(va), nb = parseFloat(vb);
        if (!isNaN(na) && !isNaN(nb)) return asc ? na - nb : nb - na;
        return asc ? va.localeCompare(vb) : vb.localeCompare(va);
      }});
      pairs.forEach(([dr, sr]) => {{ tbody.appendChild(dr); if(sr) tbody.appendChild(sr); }});
    }});
  }});
}});
</script>
</body>
</html>"""

    with open(output_file, "w", encoding='utf-8') as f:
        f.write(html)


# =============================================================================
# MAIN PROGRAM / PROGRAMME PRINCIPAL
# =============================================================================

def main():
    # ---- Language selection / Sélection de la langue ----
    # The prompt is always bilingual at this stage / L'invite est toujours bilingue ici
    T = select_language()

    print("\n" + "=" * 60)
    print(T['header_title'])
    print("=" * 60)

    # ---- Detect ntthal / Détecter ntthal ----
    print(T['detecting_ntthal'])
    ntthal_path, config_path = detect_ntthal(T)
    print(f"{T['ntthal_found']}: {ntthal_path}")
    print(f"{T['config_found']}: {config_path or T['config_auto']}")

    # ---- FASTA file / Fichier FASTA ----
    print(T['fasta_format'])
    fasta_file = input(T['fasta_prompt']).strip().strip("'\"")
    fasta_file = os.path.normpath(fasta_file)

    if not os.path.exists(fasta_file):
        print(f"\n{T['fasta_not_found']}: '{fasta_file}'")
        return

    fasta_name = os.path.splitext(os.path.basename(fasta_file))[0]

    # ---- Read primers / Lire les primers ----
    try:
        primers = read_primers(fasta_file)
    except Exception as e:
        print(f"\n{T['fasta_read_error']}: {e}")
        return

    if not primers:
        print(f"\n{T['fasta_empty']}")
        return

    print(f"\n{len(primers)} {T['primers_loaded']}: {', '.join(primers.keys())}")

    # ---- Validate sequences / Valider les séquences ----
    errors = validate_sequences(primers)
    if errors:
        print(f"\n{T['invalid_bases_title']}")
        for e in errors:
            print(e)
        resp = input(T['continue_prompt']).strip().lower()
        if resp not in T['continue_yes']:
            return

    # ---- Prepare tasks / Préparer les tâches ----
    combinations  = list(itertools.combinations_with_replacement(primers.items(), 2))
    hairpin_tasks = [(seq, "", name, "", 'HAIRPIN') for name, seq in primers.items()]
    dimer_tasks   = [(seq1, seq2, n1, n2, 'ANY') for (n1, seq1), (n2, seq2) in combinations]

    # ---- Parallel analysis / Analyse en parallèle ----
    print(f"\n{T['launching_analysis']} "
          f"({len(hairpin_tasks)} {T['hairpins_label']} + {len(dimer_tasks)} {T['dimers_label']})...")

    hairpin_results = run_analyses_parallel(
        hairpin_tasks, ntthal_path, config_path, T, description=T['progress_hairpin']
    )
    dimer_results = run_analyses_parallel(
        dimer_tasks, ntthal_path, config_path, T, description=T['progress_dimer']
    )

    # ---- Console summary / Résumé console ----
    all_results = hairpin_results + dimer_results
    n_crit = sum(1 for r in all_results if r['alert'].startswith('🔴'))
    n_warn = sum(1 for r in all_results if r['alert'].startswith('🟡'))
    n_ok   = sum(1 for r in all_results if r['alert'].startswith('🟢'))
    print(f"\n{T['summary_label']}: 🔴 {n_crit} {T['summary_problems']}  "
          f"🟡 {n_warn} {T['summary_warnings']}  🟢 {n_ok} OK")

    # ---- Write output files / Écrire les fichiers de sortie ----
    txt_file  = f"dimer_analysis_{fasta_name}.txt"
    csv_file  = f"dimer_analysis_{fasta_name}.csv"
    html_file = f"dimer_analysis_{fasta_name}.html"

    try:
        write_text_report(hairpin_results, dimer_results, txt_file, fasta_name, T)
        write_csv_report(hairpin_results, dimer_results, csv_file, T)
        write_html_report(hairpin_results, dimer_results, html_file, fasta_name, T)

        print(f"\n{T['reports_generated']}")
        print(f"   {T['txt_label']}: {txt_file}")
        print(f"   {T['csv_label']}: {csv_file}")
        print(f"   {T['html_label']}: {html_file}")
        print(f"\n{T['html_hint']}\n")

    except Exception as e:
        print(f"\n{T['write_error']}: {e}")


if __name__ == "__main__":
    main()