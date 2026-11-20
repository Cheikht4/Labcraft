import typer
import hashlib
import time
import os
import sys

# Désactiver les codes ANSI forcés en CI/pytest pour préserver la sortie --help sans séquences d'échappement
# Disable forced ANSI codes under CI/pytest so --help output stays free of escape sequences
if not sys.stdout.isatty() and (
    "pytest" in sys.modules or os.getenv("GITHUB_ACTIONS") or os.getenv("CI")
):
    import typer.rich_utils as _rich_utils

    _rich_utils.FORCE_TERMINAL = False

from labcraft import __version__
from labcraft.lamp.domains import PrimerRole
from labcraft.report.renderer import render_report
from labcraft.report.json_export import write_analyze_json
from labcraft.cli.config import PanelConfig, build_engine_from_config
from labcraft.pipeline.analysis import (
    EmptyPanelError,
    TargetFileError,
    evaluate_dimer_complexes as _evaluate_dimer_complexes,
    load_targets_from_config,
    run_single_panel_analysis,
    split_colliding_panel_configs,
)
from labcraft.thermo.vienna import dna_session
from pydantic import ValidationError

from labcraft.cli.parsers import (
    build_config_from_cli,
    parse_primer_conc_args,
    read_multi_fasta,
    validate_target_sequence,
    ParseError,
)
from typing import Any, Dict, List, Optional, Tuple

app = typer.Typer(help="LabCraft CLI - Thermodynamic multiplex primer panel design engine.")


def _version_callback(value: bool) -> None:
    """Affiche la version de LabCraft et quitte / Print LabCraft version and exit."""
    if value:
        typer.echo(f"LabCraft {__version__}")
        raise typer.Exit()


@app.callback()
def main_callback(
    version: Optional[bool] = typer.Option(
        None,
        "--version",
        callback=_version_callback,
        is_eager=True,
        help="Show the LabCraft version and exit.",
    ),
) -> None:
    """LabCraft CLI entry point."""


# Ré-exporté pour compatibilité interne / Re-exported for internal compatibility
__all__ = ["app", "_evaluate_dimer_complexes"]


def hash_file(filepath: str) -> str:
    hasher = hashlib.sha256()
    with open(filepath, 'rb') as f:
        buf = f.read()
        hasher.update(buf)
    return hasher.hexdigest()

_ENZYME_HELP = (
    "Profil enzymatique / Enzyme profile (bst2.0 [défaut], bst, bst2.0ws, gspssd : "
    "ΔG3' = -2.0 kcal/mol ; bst3.0 : -1.5 kcal/mol ; taq : -3.0 kcal/mol — "
    "heuristique non calibrée ; les données de paillasse actuelles ne permettent pas "
    "de discriminer entre -2.0 et -3.5)."
)
_DIMER_CONC_MAX_HELP = (
    "Seuil maximal de concentration totale de dimères amplifiables en nM déclenchant FAILURE "
    "(défaut : 0.230 nM — provisoire, calibré sur 5 cas de paillasse ; intervalle de séparation réel : "
    "0.154 – 0.336 nM ; duplex_half et duplex_full sont deux dilutions du même panel donc non indépendants, "
    "et la séparation repose sur une seule paire inter-panels : FIP_Nunes_YFV × F3_ParidaD4)."
)
_DIMER_WARN_RATIO_HELP = (
    "Fraction de --dimer-conc-max au-delà de laquelle les dimères amplifiables passent d'INFO à WARNING "
    "(heuristique non calibrée, défaut : 0.500 soit un plancher WARNING de 0.115 nM pour un seuil FAILURE de 0.230 nM)."
)
_ROLE_IMBALANCE_RATIO_HELP = (
    "Rapport maximal d'occupation brute toléré entre amorces d'initiation de même rôle (F3 vs B3, FIP vs BIP) "
    "avant passage d'INFO à WARNING (heuristique non calibrée, défaut : 20.0 ≈ ΔΔG de 2,0 kcal/mol à 63 °C, "
    "sans normalisation par le nombre de variants dégénérés)."
)
_BOLTZMANN_MAX_DIMERS_HELP = (
    "Garde-fou optionnel pour les très gros panels : plafond K du nombre de groupes non amplifiables "
    "dédoublonnés par clé (couple parent, 3'_A, 3'_B) (parmi ceux >= 0.1 × --dimer-conc-max, triés par "
    "concentration décroissante) évalués par calc_boltzmann_3p_pairing, en plus de 100 % des dimères "
    "amplifiables (défaut : illimité ; 0 ou non spécifié = pas de plafond)."
)
_BOLTZMANN_DEDUP_3P_WINDOW_HELP = (
    "Longueur en nt de la fenêtre 3' par variant dans la clé de dédoublonnage Boltzmann "
    "(couple parent, 3'_A, 3'_B) : fusionne les variants ne différant qu'en 5'/milieu et évalue "
    "séparément ceux différant dans les N dernières bases 3' (défaut : 5)."
)


@app.command()
def analyze(
    config: Optional[str] = typer.Option(None, "-c", "--config", help="Path to the YAML panel configuration file."),
    primers: Optional[str] = typer.Option(None, "-p", "--primers", help="Path to primers FASTA/TXT file."),
    targets: Optional[str] = typer.Option(None, "-t", "--targets", help="Path to targets FASTA/TXT file."),
    output: str = typer.Option("report.html", "-o", "--output", help="Output HTML report path."),
    temperature: Optional[float] = typer.Option(None, "--temperature", help="Experiment temperature in Celsius."),
    enzyme: Optional[str] = typer.Option(None, "--enzyme", help=_ENZYME_HELP),
    na: Optional[float] = typer.Option(None, "--na", help="Na+ concentration in mM."),
    k: Optional[float] = typer.Option(None, "--k", help="K+ concentration in mM."),
    tris: Optional[float] = typer.Option(None, "--tris", help="Tris concentration in mM."),
    mg: Optional[float] = typer.Option(None, "--mg", help="Mg2+ concentration in mM."),
    dntp: Optional[float] = typer.Option(None, "--dntp", help="dNTP concentration in mM."),
    conc_fip_bip: Optional[float] = typer.Option(None, "--conc-fip-bip", help="FIP/BIP concentration in uM."),
    conc_f3_b3: Optional[float] = typer.Option(None, "--conc-f3-b3", help="F3/B3 concentration in uM."),
    conc_loop: Optional[float] = typer.Option(None, "--conc-loop", help="Loop primer concentration in uM."),
    conc: Optional[List[str]] = typer.Option(None, "--conc", metavar="NOM=µM", help="Per-primer concentration in uM (e.g. --conc Nunes_YFV_FIP=0.1). Repeatable."),  # noqa: B008
    copies: Optional[float] = typer.Option(None, "--copies", help="Target copies per uL."),
    max_variants: int = typer.Option(64, "--max-variants", help="Maximum variants allowed per degenerate primer (default: 64)."),
    errors: int = typer.Option(2, "--errors", help="Max errors for candidate site seeding outside 3' zone."),
    strict_3prime: int = typer.Option(3, "--strict-3prime", help="Length of strict 3' zone for candidate site seeding."),
    strict_3prime_tolerate: int = typer.Option(2, "--strict-3prime-tolerate", help="Number of tolerated substitutions in strict 3' zone for seeding."),
    dimer_conc_max: float = typer.Option(0.23, "--dimer-conc-max", help=_DIMER_CONC_MAX_HELP),
    dimer_warn_ratio: float = typer.Option(0.5, "--dimer-warn-ratio", help=_DIMER_WARN_RATIO_HELP),
    role_imbalance_ratio: float = typer.Option(20.0, "--role-imbalance-ratio", help=_ROLE_IMBALANCE_RATIO_HELP),
    boltzmann_max_dimers: Optional[int] = typer.Option(None, "--boltzmann-max-dimers", min=0, help=_BOLTZMANN_MAX_DIMERS_HELP),
    boltzmann_dedup_3p_window: int = typer.Option(5, "--boltzmann-dedup-3p-window", min=0, help=_BOLTZMANN_DEDUP_3P_WINDOW_HELP),
    separate_targets: bool = typer.Option(False, "--separate-targets", "--per-target", help="Solve targets independently instead of joint multi-target equilibrium."),
    force_joint: bool = typer.Option(False, "--force-joint", help="Force joint analysis of multiple panels for the same target.")
) -> None:
    """
    Analyzes a multiplex primer panel and generates a comprehensive thermodynamic diagnostic report.
    """
    if not config and not primers:
        typer.echo("Erreur: Vous devez spécifier au moins un fichier de configuration YAML (-c) ou un fichier d'amorces (-p).")
        raise typer.Exit(code=1)

    # Vérification de l'existence des fichiers en amont (avant hash_file)
    # Check file existence upfront (before hash_file)
    for fpath in (config, primers, targets):
        if fpath and not os.path.isfile(fpath):
            typer.echo(f"Erreur de fichier introuvable : {fpath}")
            raise typer.Exit(code=1)

    typer.echo("Loading configuration...")
    start_time = time.time()

    # Empreinte des entrées (avant réaffectation de `primers`) / Input fingerprint (before `primers` is reassigned)
    inputs_fp = {
        "config": {"path": config, "sha256": hash_file(config)} if config else None,
        "primers": {"path": primers, "sha256": hash_file(primers)} if primers else None,
        "targets": {"path": targets, "sha256": hash_file(targets)} if targets else None,
    }
    
    # Use config file for hash if available, otherwise primers file
    file_hash = hash_file(config) if config else (hash_file(primers) if primers else "no-file")

    try:
        primer_concs = parse_primer_conc_args(conc)
        config_obj = build_config_from_cli(
            config_path=config,
            primers_path=primers,
            targets_path=targets,
            temperature=temperature,
            enzyme=enzyme,
            na=na, k=k, tris=tris, mg=mg, dntp=dntp,
            conc_fip_bip=conc_fip_bip,
            conc_f3_b3=conc_f3_b3,
            conc_loop=conc_loop,
            primer_concs=primer_concs,
            copies=copies
        )
    except FileNotFoundError as e:
        typer.echo(f"Erreur de fichier introuvable : {e.filename}")
        raise typer.Exit(code=1) from e
    except ParseError as e:
        typer.echo(f"Erreur de lecture de fichier :\n{e}")
        raise typer.Exit(code=1) from e
    except ValidationError as e:
        typer.echo(f"Erreur de validation de la configuration :\n{e}")
        raise typer.Exit(code=1) from e

    # Session ADN unique pour toute l'analyse (chargement unique Mathews 2004, restauration finale Turner 2004)
    # Single DNA session for the whole analysis (single Mathews 2004 load, final Turner 2004 restore)
    with dna_session(temp_celsius=config_obj.experiment.temperature_C):
        _run_analyze_with_config(
            config=config,
            config_obj=config_obj,
            output=output,
            start_time=start_time,
            inputs_fp=inputs_fp,
            file_hash=file_hash,
            max_variants=max_variants,
            errors=errors,
            strict_3prime=strict_3prime,
            strict_3prime_tolerate=strict_3prime_tolerate,
            dimer_conc_max=dimer_conc_max,
            dimer_warn_ratio=dimer_warn_ratio,
            role_imbalance_ratio=role_imbalance_ratio,
            boltzmann_max_dimers=boltzmann_max_dimers,
            boltzmann_dedup_3p_window=boltzmann_dedup_3p_window,
            separate_targets=separate_targets,
            force_joint=force_joint,
        )


def _run_analyze_with_config(
    *,
    config: Optional[str],
    config_obj: PanelConfig,
    output: str,
    start_time: float,
    inputs_fp: Dict[str, Any],
    file_hash: str,
    max_variants: int,
    errors: int,
    strict_3prime: int,
    strict_3prime_tolerate: int,
    dimer_conc_max: float,
    dimer_warn_ratio: float = 0.5,
    role_imbalance_ratio: float = 20.0,
    boltzmann_max_dimers: Optional[int] = None,
    boltzmann_dedup_3p_window: int = 5,
    separate_targets: bool,
    force_joint: bool,
) -> None:
    """Délègue l'orchestration scientifique au module `labcraft.pipeline.analysis` et écrit les rapports.
    Delegates scientific orchestration to `labcraft.pipeline.analysis` and writes HTML/JSON reports.
    """
    try:
        targets_dict, target_warnings = load_targets_from_config(
            config_obj,
            explicit_config=bool(config),
            progress_callback=typer.echo,
        )
    except TargetFileError as e:
        typer.echo(f"Erreur lors de la lecture du fichier cible {e.filepath}:\n{e.cause}")
        raise typer.Exit(code=1) from e

    config_list, had_collisions = split_colliding_panel_configs(
        config_obj, force_joint=force_joint
    )
    if had_collisions:
        typer.echo(
            "Attention: Plusieurs panels alternatifs ciblant la même séquence ont été détectés."
        )
        typer.echo("Analyse séparée des panels...")

    for i_cfg, current_config in enumerate(config_list):
        if len(config_list) > 1:
            panel_n = current_config.primer_sets[0].panel_name or str(i_cfg)
            current_output = output.replace(".html", f"_{panel_n}.html")
            typer.echo(f"--- Analyse du panel {panel_n} ---")
        else:
            panel_n = None
            current_output = output

        try:
            panel_res = run_single_panel_analysis(
                current_config,
                targets_dict,
                panel_label=panel_n,
                target_warnings=target_warnings,
                inputs_fp=inputs_fp,
                file_hash=file_hash,
                max_variants=max_variants,
                errors=errors,
                strict_3prime=strict_3prime,
                strict_3prime_tolerate=strict_3prime_tolerate,
                dimer_conc_max=dimer_conc_max,
                dimer_warn_ratio=dimer_warn_ratio,
                role_imbalance_ratio=role_imbalance_ratio,
                boltzmann_max_dimers=boltzmann_max_dimers,
                boltzmann_dedup_3p_window=boltzmann_dedup_3p_window,
                separate_targets=separate_targets,
                progress_callback=typer.echo,
            )
        except EmptyPanelError as e:
            typer.echo(str(e))
            raise typer.Exit(code=1) from e
        except (ParseError, ValueError) as e:
            typer.echo(f"Erreur de validation des séquences :\n{e}")
            raise typer.Exit(code=1) from e

        typer.echo("Generating HTML report...")
        html = render_report(
            panel_res.verdict,
            panel_res.fractions,
            panel_res.risks,
            panel_res.metadata,
        )
        with open(current_output, "w", encoding="utf-8") as f:
            f.write(html)

        # Sortie machine JSON à côté du HTML / Machine-readable JSON next to the HTML
        json_out = write_analyze_json(current_output, panel_res.json_payload)
        typer.echo(f"JSON saved to {json_out}.")
        typer.echo(
            f"Analysis complete in {time.time() - start_time:.2f}s. Report saved to {current_output}."
        )
        


@app.command()
def coverage(
    primers: str = typer.Option(..., "-p", "--primers", help="Path to primers FASTA/TXT file."),
    fasta: str = typer.Option(..., "-f", "--fasta", help="Path to multi-FASTA file of strains."),
    sites: Optional[str] = typer.Option(None, "-s", "--sites", help="Path to candidate sites CSV (optional, automatic screening if omitted)."),
    output: str = typer.Option("coverage_report.html", "-o", "--output", help="Output HTML report path."),
    panel: Optional[str] = typer.Option(None, "--panel", help="Panel name to analyze if multiple panels are present in primers file."),
    temperature: float = typer.Option(65.0, "--temperature", help="Experiment temperature in Celsius."),
    enzyme: Optional[str] = typer.Option(None, "--enzyme", help=_ENZYME_HELP),
    na: float = typer.Option(50.0, "--na", help="Na+ concentration in mM."),
    mg: float = typer.Option(8.0, "--mg", help="Mg2+ concentration in mM."),
    dntp: float = typer.Option(1.4, "--dntp", help="dNTP concentration in mM."),
    ddg_max: float = typer.Option(3.0, "--ddg-max", help="Maximum allowable delta-delta-G penalty (kcal/mol) relative to perfect duplex."),
    dg_threshold: Optional[float] = typer.Option(None, "--dg-threshold", help="Optional absolute viability dG ceiling in kcal/mol (disabled by default)."),
    max_mismatches: int = typer.Option(2, "--max-mismatches", help="Counting rule threshold."),
    errors: int = typer.Option(2, "--errors", help="Max errors for seeding outside 3' zone."),
    strict_3prime: int = typer.Option(3, "--strict-3prime", help="Length of strict 3' zone for seeding."),
    strict_3prime_tolerate: int = typer.Option(2, "--strict-3prime-tolerate", help="Number of tolerated substitutions in strict 3' zone for seeding."),
    export_sites: Optional[str] = typer.Option(None, "--export-sites", help="Export candidate sites to CSV.")
) -> None:
    """Analyse de la couverture thermodynamique multi-souches."""
    import csv
    import time
    import tempfile
    from labcraft.cli.parsers import build_config_from_cli, read_multi_fasta
    from labcraft.cli.config import build_engine_from_config
    from labcraft.lamp.coverage import CoverageAnalyzer
    from labcraft.target.seeding import find_candidate_sites
    from labcraft.report.coverage_report import render_coverage_report

    # Vérification de l'existence des fichiers en amont en mode coverage
    # Upfront file existence check in coverage mode
    for fpath in (primers, fasta, sites):
        if fpath and not os.path.isfile(fpath):
            typer.echo(f"Erreur de fichier introuvable : {fpath}")
            raise typer.Exit(code=1)

    first_strain_path: Optional[str] = None
    try:
        print("Lecture des entrées...")
        fasta_list: List[Tuple[str, str]] = []
        for name, seq in read_multi_fasta(fasta):
            name = name.split()[0]
            clean_seq, w_msg = validate_target_sequence(name, seq, warn_dashes=True)
            if w_msg:
                typer.echo(f"Avertissement : {w_msg}")
            fasta_list.append((name, clean_seq))
        fasta_dict = {name: seq for name, seq in fasta_list}
        
        if fasta_list:
            with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", delete=False, suffix=".fasta") as tmp:
                tmp.write(f">{fasta_list[0][0]}\n{fasta_list[0][1]}\n")
                first_strain_path = tmp.name

        config_obj = build_config_from_cli(
            config_path=None,
            primers_path=primers,
            targets_path=first_strain_path,
            temperature=temperature,
            enzyme=enzyme,
            na=na,
            mg=mg,
            dntp=dntp,
            allow_unmatched_targets=True
        )
    except FileNotFoundError as e:
        typer.echo(f"Erreur de fichier introuvable : {e.filename}")
        raise typer.Exit(code=1) from e
    except (ParseError, ValueError) as e:
        typer.echo(f"Erreur de lecture de fichier :\n{e}")
        raise typer.Exit(code=1) from e
    except ValidationError as e:
        typer.echo(f"Erreur de validation de la configuration :\n{e}")
        raise typer.Exit(code=1) from e
    finally:
        # Nettoyage du fichier temporaire de première souche / Clean up temporary first-strain FASTA
        if first_strain_path and os.path.exists(first_strain_path):
            try:
                os.remove(first_strain_path)
            except OSError:
                pass
    
    # Sélection du panel si plusieurs sont présents
    if panel:
        filtered_sets = [ps for ps in config_obj.primer_sets if ps.panel_name == panel or ps.target == panel]
        if not filtered_sets:
            available = [ps.panel_name for ps in config_obj.primer_sets]
            raise ValueError(f"Panel '{panel}' non trouvé. Panels disponibles : {available}")
        config_obj.primer_sets = filtered_sets
    elif len(config_obj.primer_sets) > 1:
        available = [ps.panel_name for ps in config_obj.primer_sets]
        print(f"Plusieurs panels détectés ({', '.join(available)}). Analyse du premier panel '{config_obj.primer_sets[0].panel_name}'. Utilisez --panel pour en choisir un autre.")
        config_obj.primer_sets = [config_obj.primer_sets[0]]

    selected_panel_name = config_obj.primer_sets[0].panel_name if config_obj.primer_sets else "DefaultPanel"

    dummy_targets = {fasta_list[0][0]: fasta_list[0][1]} if fasta_list else {}
    physical_primers, _, backend, backend_kwargs, _, enzyme_obj, _, _ = build_engine_from_config(config_obj, dummy_targets)
    buf_cfg = config_obj.experiment.buffer
    if buf_cfg:
        print(
            f"Conditions utilisées : T={temperature:.1f} °C, Enzyme={enzyme_obj.name}, "
            f"Na={buf_cfg.na_mM} mM, Mg={buf_cfg.mg_mM} mM, dNTP={buf_cfg.dntp_mM} mM"
        )
    
    print(f"Chargement de {len(fasta_dict)} souches.")
    
    csv_records = []
    if sites:
        print(f"Lecture des sites candidats depuis : {sites}")
        with open(sites, 'r', encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                csv_records.append(row)
    else:
        print(f"Criblage interne des sites candidats (erreurs <= {errors}, zone 3' stricte = {strict_3prime}, tolérance 3' = {strict_3prime_tolerate})...")
        t_seed_0 = time.time()
        csv_records = find_candidate_sites(
            fasta_dict,
            physical_primers,
            max_errors=errors,
            strict_3prime_len=strict_3prime,
            strict_3prime_tolerate=strict_3prime_tolerate,
            panel_name=selected_panel_name
        )
        t_seed_1 = time.time()
        print(f"Criblage terminé en {t_seed_1 - t_seed_0:.3f} s ({len(csv_records)} sites trouvés).")

    if export_sites:
        with open(export_sites, 'w', newline='', encoding="utf-8") as f:
            fieldnames = ["strain_id", "primer_role", "primer_name", "position", "strand", "site_seq", "n_mismatches", "panel"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in csv_records:
                writer.writerow({k: r.get(k, "") for k in fieldnames})
        print(f"Sites candidats exportés dans : {export_sites}")

    analyzer = CoverageAnalyzer(
        physical_primers, fasta_dict, backend, enzyme_obj,
        temp_celsius=temperature,
        ddg_max=ddg_max,
        dg_threshold=dg_threshold,
        max_mismatches_count=max_mismatches,
        **backend_kwargs,
    )
    
    print("Analyse thermodynamique des souches...")
    t0 = time.time()
    verdicts = analyzer.analyze_strains(csv_records)
    t1 = time.time()
    
    print(f"Analyse terminée en {t1 - t0:.2f} s.")
    
    render_coverage_report(
        verdicts, output, selected_panel_name,
        ddg_max=ddg_max, dg_threshold=dg_threshold,
        max_mismatches_count=max_mismatches, temperature_C=temperature
    )
    print(f"Rapport généré : {output}")


@app.command()
def redesign(
    primers: str = typer.Option(..., "-p", "--primers", help="Path to primers FASTA/TXT file."),
    targets: str = typer.Option(..., "-t", "--targets", help="Path to targets FASTA/TXT file."),
    primer: str = typer.Option(..., "--primer", help="Name of the primer to redesign (F3, B3, LF or LB), e.g. ParidaD4_F3."),
    shift: int = typer.Option(12, "--shift", help="Maximum shift (nt) of the 3' end around the legitimate locus."),
    length: str = typer.Option("17-22", "--length", help="Candidate length range, e.g. 17-22."),
    top: int = typer.Option(5, "--top", help="Number of best candidates re-evaluated by the full coupled equilibrium."),
    output: str = typer.Option("redesign_report.html", "-o", "--output", help="Output HTML report path (JSON written alongside)."),
    temperature: Optional[float] = typer.Option(None, "--temperature", help="Experiment temperature in Celsius."),
    enzyme: Optional[str] = typer.Option(None, "--enzyme", help=_ENZYME_HELP),
    na: Optional[float] = typer.Option(None, "--na", help="Na+ concentration in mM."),
    k: Optional[float] = typer.Option(None, "--k", help="K+ concentration in mM."),
    tris: Optional[float] = typer.Option(None, "--tris", help="Tris concentration in mM."),
    mg: Optional[float] = typer.Option(None, "--mg", help="Mg2+ concentration in mM."),
    dntp: Optional[float] = typer.Option(None, "--dntp", help="dNTP concentration in mM."),
    conc_fip_bip: Optional[float] = typer.Option(None, "--conc-fip-bip", help="FIP/BIP concentration in uM."),
    conc_f3_b3: Optional[float] = typer.Option(None, "--conc-f3-b3", help="F3/B3 concentration in uM."),
    conc_loop: Optional[float] = typer.Option(None, "--conc-loop", help="Loop primer concentration in uM."),
    conc: Optional[List[str]] = typer.Option(None, "--conc", metavar="NOM=µM", help="Per-primer concentration in uM (NAME=VALUE). Repeatable."),  # noqa: B008
    copies: Optional[float] = typer.Option(None, "--copies", help="Target copies per uL."),
    max_variants: int = typer.Option(64, "--max-variants", help="Maximum variants allowed per degenerate primer."),
    dimer_conc_max: float = typer.Option(0.23, "--dimer-conc-max", help=_DIMER_CONC_MAX_HELP),
    dimer_warn_ratio: float = typer.Option(0.5, "--dimer-warn-ratio", help=_DIMER_WARN_RATIO_HELP),
    role_imbalance_ratio: float = typer.Option(20.0, "--role-imbalance-ratio", help=_ROLE_IMBALANCE_RATIO_HELP),
) -> None:
    """Propose des séquences de remplacement pour une amorce externe ou de boucle.
    Proposes replacement sequences for an outer or loop primer.
    """
    import json
    from labcraft.optimize.redesign import (
        IncrementalRedesignEvaluator,
        RedesignError,
        check_role_redesignable,
        compute_duplex_energy,
        compute_inner_comparison,
        compute_max_3p_complementarity,
        compute_unfolding_energy,
        equilibrium_sort_key,
        generate_candidates,
        geometry_satisfied,
        locate_on_target,
        panel_geometry,
        prescreen_dimers,
        prescreen_sort_key,
        resolve_role_inner_reference,
        select_top_candidates,
    )
    from labcraft.lamp.domains import PhysicalPrimer
    from labcraft.cli.parsers import parse_primer_name
    from labcraft.report.json_export import json_path_for_output
    from labcraft.report.renderer import render_redesign_report

    for fpath in (primers, targets):
        if not os.path.isfile(fpath):
            typer.echo(f"Erreur de fichier introuvable : {fpath}")
            raise typer.Exit(code=1)

    # 1. Rôle demandé : refus explicite des FIP/BIP / Requested role: explicit FIP/BIP refusal
    try:
        _pan, role_req, _ver = parse_primer_name(primer)
        check_role_redesignable(role_req, primer)
        lo_s, _, hi_s = length.partition("-")
        len_min = int(lo_s)
        len_max = int(hi_s) if hi_s else len_min
        if len_min <= 0 or len_max < len_min:
            raise RedesignError(f"Plage de longueur invalide : '{length}'.")
    except (RedesignError, ParseError, ValueError) as e:
        typer.echo(f"Erreur : {e}")
        raise typer.Exit(code=1) from e

    # 2. Configuration et moteur, identiques à analyze / Config and engine, identical to analyze
    try:
        primer_concs = parse_primer_conc_args(conc)
        config_obj = build_config_from_cli(
            primers_path=primers, targets_path=targets, temperature=temperature, enzyme=enzyme,
            na=na, k=k, tris=tris, mg=mg, dntp=dntp,
            conc_fip_bip=conc_fip_bip, conc_f3_b3=conc_f3_b3, conc_loop=conc_loop,
            primer_concs=primer_concs, copies=copies,
        )
        targets_dict: Dict[str, str] = {}
        for h, seq in read_multi_fasta(targets):
            t_id = h.split()[0]
            clean_seq, _w = validate_target_sequence(t_id, seq, warn_dashes=False)
            targets_dict[t_id] = clean_seq
        engine_primers, primer_to_target, backend, backend_kwargs, mon_molar, enzyme_obj, temp_celsius, profiles = (
            build_engine_from_config(config_obj, targets_dict, max_variants=max_variants)
        )
    except (ParseError, ValueError, ValidationError) as e:
        typer.echo(f"Erreur de lecture de fichier :\n{e}")
        raise typer.Exit(code=1) from e

    primer_to_panel_map: Dict[str, str] = {}
    for ps in config_obj.primer_sets:
        p_id = ps.panel_name if ps.panel_name else ps.target
        for r_n in ps.primers:
            primer_to_panel_map[f"{r_n}_{p_id}"] = p_id
    for ep in engine_primers:
        ep_par = ep.parent_name or ep.name.split("#")[0]
        if ep_par in primer_to_panel_map:
            primer_to_panel_map[ep.name] = primer_to_panel_map[ep_par]
        elif ep.name not in primer_to_panel_map:
            primer_to_panel_map[ep.name] = primer_to_target.get(ep.name, "")

    # 3. Retrouver l'amorce visée / Locate the requested primer
    sel = None
    for pset in config_obj.primer_sets:
        panel_id = pset.panel_name or pset.target
        for role_name, p_data in pset.primers.items():
            names = {p_data.raw_name, f"{role_name}_{panel_id}", f"{panel_id}_{role_name}"}
            if primer in names:
                sel = (pset, panel_id, role_name, p_data)
                break
        if sel is not None:
            break
    if sel is None:
        typer.echo(f"Erreur : amorce '{primer}' introuvable dans {primers}.")
        raise typer.Exit(code=1)
    pset, panel_id, role_name, p_data = sel
    internal_name = f"{role_name}_{panel_id}"
    target_id = pset.target
    target_seq = targets_dict.get(target_id, "")
    from labcraft.thermo.lna import parse_lna_sequence
    orig_seq, _lna = parse_lna_sequence(p_data.seq)
    orig_len = len(orig_seq)

    try:
        orig_locus = locate_on_target(orig_seq, target_seq)
        if orig_locus is None:
            raise RedesignError(f"Locus légitime de '{primer}' introuvable sur la cible '{target_id}'.")
        panel_parents: Dict[str, PhysicalPrimer] = {}
        for p in engine_primers:
            par_id = p.parent_name or p.name.split("#")[0]
            if par_id == f"{p.role.value}_{panel_id}" and p.role.value not in panel_parents:
                panel_parents[p.role.value] = p
        geo = panel_geometry(panel_parents, target_seq)
        cands = generate_candidates(
            orig_seq, orig_locus, target_seq, shift, len_min, len_max, raw_original_seq=p_data.seq
        )
        for c in cands:
            c.geometry_ok = geometry_satisfied(role_name, c.locus, geo)
    except RedesignError as e:
        typer.echo(f"Erreur : {e}")
        raise typer.Exit(code=1) from e

    ref = cands[0]
    valid = [c for c in cands[1:] if c.geometry_ok]
    typer.echo(
        f"Amorce {primer} ({role_name}) : locus {orig_locus.start}-{orig_locus.end} brin {orig_locus.strand} ; "
        f"{len(valid)} candidats géométriquement valides sur {len(cands) - 1}."
    )

    # 4. & 5. Session ADN unique pour le criblage et l'évaluation incrémentale en mémoire
    # Single DNA session across prescreen and in-memory incremental coupled equilibrium
    total_conc_uM = sum(
        float(p.nominal_concentration or 0.0) * 1e6 for p in engine_primers
        if (p.parent_name or p.name.split("#")[0]) == internal_name
    )
    partners = [p for p in engine_primers if (p.parent_name or p.name.split("#")[0]) != internal_name]
    role_enum = PrimerRole[role_name.upper()]
    pair_cache: Dict[Tuple[Any, ...], Tuple[bool, float, Dict[str, Any]]] = {}
    buf_dump = config_obj.experiment.buffer.model_dump() if config_obj.experiment.buffer else None

    with dna_session(temp_celsius=temp_celsius):
        ref_inner_domain, dg_ref_inner, dg_f2 = resolve_role_inner_reference(
            role_name, engine_primers, panel_id, backend, temp_celsius, backend_kwargs
        )
        for c in [ref] + valid:
            compute_duplex_energy(c, backend, temp_celsius, backend_kwargs)
            compute_inner_comparison(c, ref_inner_domain, dg_ref_inner, dg_f2)
            compute_max_3p_complementarity(c, partners)
            prescreen_dimers(
                c, internal_name, role_enum, total_conc_uM * 1e-6, partners,
                backend, enzyme_obj, temp_celsius, backend_kwargs, pair_cache=pair_cache,
            )

        # Tri préliminaire sans ouverture génomique coûteuse ; `dg_open` n'intervient qu'au dernier
        # niveau de départage : on ne le calcule que pour la référence, le top retenu et tout ex aequo frontière.
        # Preliminary sort without costly genomic pf() unfolding; `dg_open` only enters the final
        # tie-breaker: compute it only for reference, selected top candidates, and any boundary ties.
        valid.sort(key=lambda c: prescreen_sort_key(c, orig_len)[:-1])
        prelim_best, _ = select_top_candidates(valid, top, role_name)
        tied_prefixes = {prescreen_sort_key(c, orig_len)[:-1] for c in prelim_best}
        needs_unfold = [ref] + [
            c for c in valid if c in prelim_best or prescreen_sort_key(c, orig_len)[:-1] in tied_prefixes
        ]
        seen_unfold: set[int] = set()
        for c in needs_unfold:
            if id(c) not in seen_unfold:
                seen_unfold.add(id(c))
                compute_unfolding_energy(c, role_name, target_seq, temp_celsius, mon_molar)

        valid.sort(key=lambda c: prescreen_sort_key(c, orig_len))
        best, n_discarded_inner = select_top_candidates(valid, top, role_name)
        for c in [ref] + best:
            if id(c) not in seen_unfold:
                seen_unfold.add(id(c))
                compute_unfolding_energy(c, role_name, target_seq, temp_celsius, mon_molar)
        if dg_ref_inner is not None:
            typer.echo(
                f"Référence interne {ref_inner_domain} ({panel_id}) : ΔG({ref_inner_domain}) = {dg_ref_inner:.2f} kcal/mol "
                f"({n_discarded_inner} candidat(s) plus stable(s) que {ref_inner_domain} écarté(s))."
            )

        evaluator = IncrementalRedesignEvaluator(
            engine_primers=engine_primers,
            internal_name=internal_name,
            target_id=target_id,
            role_enum=role_enum,
            total_conc_molar=total_conc_uM * 1e-6,
            targets_dict=targets_dict,
            primer_to_target=primer_to_target,
            primer_to_panel_map=primer_to_panel_map,
            profiles=profiles,
            buffer_dump=buf_dump,
            backend=backend,
            backend_kwargs=backend_kwargs,
            mon_molar=mon_molar,
            enzyme=enzyme_obj,
            temp_celsius=temp_celsius,
            max_variants=max_variants,
            dimer_conc_max=dimer_conc_max,
            dimer_warn_ratio=dimer_warn_ratio,
            role_imbalance_ratio=role_imbalance_ratio,
            pair_cache=pair_cache,
        )
        for i_c, c in enumerate([ref] + best):
            typer.echo(f"Équilibre couplé : candidat {i_c} (décalage {c.shift:+d}, longueur {c.length})...")
            c.equilibrium = evaluator.evaluate_candidate(c)

    # 6. Sortie HTML (template Jinja2) et JSON / HTML (Jinja2 template) and JSON output
    ranked_sorted = [ref] + sorted(best, key=lambda c: equilibrium_sort_key(c, orig_len))
    n_discarded_f2 = n_discarded_inner if ref_inner_domain == "F2" else sum(1 for c in valid if c.more_stable_than_f2)
    payload = {
        "primer": primer,
        "role": role_name,
        "panel": panel_id,
        "target": target_id,
        "original_sequence": orig_seq,
        "original_locus": {"start": orig_locus.start, "end": orig_locus.end, "strand": orig_locus.strand},
        "dg_f2": dg_f2,
        "n_discarded_more_stable_than_f2": n_discarded_f2,
        "ref_inner_domain": ref_inner_domain,
        "dg_ref_inner": dg_ref_inner,
        "n_discarded_more_stable_than_ref_inner": n_discarded_inner,
        "shift": shift,
        "length_range": [len_min, len_max],
        "n_candidates_total": len(cands) - 1,
        "n_candidates_geometry_ok": len(valid),
        "ranked": [c.to_dict() for c in ranked_sorted],
        "prescreen": [c.to_dict() for c in valid],
    }
    json_out = json_path_for_output(output)
    with open(json_out, "w", encoding="utf-8") as jf:
        json.dump(payload, jf, indent=2, ensure_ascii=False)

    template_rows = []
    for idx_rank, c in enumerate(ranked_sorted):
        eq = c.equilibrium or {}
        acc = eq.get("initiation_access_by_panel", {})
        acc_txt = "; ".join(
            f"{pan}: {(v.get('min_occupation') or 0) * 100:.1f} %"
            for pan, v in acc.items()
            if isinstance(v, dict)
        )
        fams = ", ".join(eq.get("families", [])) or "—"
        template_rows.append({
            "is_reference": c.is_reference,
            "label": "Référence" if c.is_reference else f"#{idx_rank}",
            "shift": c.shift,
            "length": c.length,
            "gc_percent": c.gc_percent,
            "sequence": c.sequence,
            "locus_start": c.locus.start,
            "locus_end": c.locus.end,
            "locus_strand": c.locus.strand,
            "dg_duplex": c.dg_duplex,
            "ddg_vs_ref": c.ddg_vs_ref_inner if c.ddg_vs_ref_inner is not None else c.ddg_vs_f2,
            "more_stable_than_ref": c.more_stable_than_ref_inner if c.dg_ref_inner is not None else c.more_stable_than_f2,
            "dg_open": c.dg_open,
            "max_3p_comp": c.max_3p_comp,
            "max_3p_partner": c.max_3p_partner,
            "max_rev_3p_comp": c.max_reverse_3p_comp,
            "max_rev_3p_partner": c.max_reverse_3p_partner,
            "max_self_3p_comp": c.max_self_3p_comp,
            "n_amplifiable_pairs": c.n_amplifiable_pairs,
            "amplifiable_total_nm": float(eq.get("amplifiable_total_nm", float("nan"))),
            "amplifiable_inter_nm": float(eq.get("amplifiable_inter_nm", float("nan"))),
            "initiation_access_txt": acc_txt,
            "mispriming_alerts": eq.get("mispriming_alerts", 0),
            "verdict": eq.get("verdict", ""),
            "families_txt": fams,
        })

    html = render_redesign_report({
        "primer": primer,
        "role_name": role_name,
        "panel_id": panel_id,
        "target_id": target_id,
        "orig_locus": orig_locus,
        "ref_inner_domain": ref_inner_domain,
        "dg_ref_inner": dg_ref_inner,
        "dg_f2": dg_f2,
        "n_discarded_inner": n_discarded_inner,
        "shift": shift,
        "len_min": len_min,
        "len_max": len_max,
        "n_valid": len(valid),
        "n_total": len(cands) - 1,
        "temp_celsius": temp_celsius,
        "enzyme_name": enzyme_obj.name,
        "rows": template_rows,
    })
    with open(output, "w", encoding="utf-8") as hf:
        hf.write(html)

    for c in ranked_sorted:
        eq = c.equilibrium or {}
        ddg_disp = c.ddg_vs_ref_inner if c.ddg_vs_ref_inner is not None else c.ddg_vs_f2
        exc_disp = c.more_stable_than_ref_inner if c.dg_ref_inner is not None else c.more_stable_than_f2
        ref_cmp = (
            f"ΔΔG({ref_inner_domain})={ddg_disp:+.2f}{f'(>{ref_inner_domain}!)' if exc_disp else ''}"
            if ddg_disp is not None
            else ""
        )
        typer.echo(
            f"{'REF ' if c.is_reference else '    '}décalage {c.shift:+d} L={c.length} GC={c.gc_percent:.1f}% "
            f"{c.locus.start}-{c.locus.end}({c.locus.strand}) {c.sequence} | "
            f"ΔG={c.dg_duplex:.2f} {ref_cmp} 3'max={c.max_3p_comp}nt | "
            f"dimères amplifiables {eq.get('amplifiable_total_nm', float('nan')):.4f} nM "
            f"(inter {eq.get('amplifiable_inter_nm', float('nan')):.4f}) | verdict {eq.get('verdict', '')}"
        )
    typer.echo(f"Rapport de redessin : {output} ; JSON : {json_out}")


if __name__ == "__main__":
    app()
