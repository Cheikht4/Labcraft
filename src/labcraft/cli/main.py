import typer
import yaml
import hashlib
import time
import importlib.metadata
import os

from labcraft.lamp.domains import PrimerRole
from labcraft.lamp.complex_enumeration import enumerate_complexes
from labcraft.solver.dual import solve_dual
from labcraft.metrics.fractions import compute_fractions
from labcraft.diagnostics.amplifiable_dimer import is_amplifiable_dimer, evaluate_pair_amplifiable
from labcraft.diagnostics.probe_tm import check_probes_tm
from labcraft.metrics.risk import evaluate_risks
from labcraft.metrics.verdict import generate_verdict
from labcraft.report.renderer import render_report
from labcraft.cli.config import PanelConfig, build_engine_from_config
from pydantic import ValidationError

from labcraft.cli.parsers import build_config_from_cli, read_multi_fasta, ParseError
from typing import Optional

app = typer.Typer(help="LabCraft CLI - Thermodynamic multiplex primer panel design engine.")

def hash_file(filepath: str) -> str:
    hasher = hashlib.sha256()
    with open(filepath, 'rb') as f:
        buf = f.read()
        hasher.update(buf)
    return hasher.hexdigest()

@app.command()
def analyze(
    config: Optional[str] = typer.Option(None, "-c", "--config", help="Path to the YAML panel configuration file."),
    primers: Optional[str] = typer.Option(None, "-p", "--primers", help="Path to primers FASTA/TXT file."),
    targets: Optional[str] = typer.Option(None, "-t", "--targets", help="Path to targets FASTA/TXT file."),
    output: str = typer.Option("report.html", "-o", "--output", help="Output HTML report path."),
    temperature: Optional[float] = typer.Option(None, "--temperature", help="Experiment temperature in Celsius."),
    enzyme: Optional[str] = typer.Option(None, "--enzyme", help="Enzyme name (e.g. bst2.0)."),
    na: Optional[float] = typer.Option(None, "--na", help="Na+ concentration in mM."),
    k: Optional[float] = typer.Option(None, "--k", help="K+ concentration in mM."),
    tris: Optional[float] = typer.Option(None, "--tris", help="Tris concentration in mM."),
    mg: Optional[float] = typer.Option(None, "--mg", help="Mg2+ concentration in mM."),
    dntp: Optional[float] = typer.Option(None, "--dntp", help="dNTP concentration in mM."),
    conc_fip_bip: Optional[float] = typer.Option(None, "--conc-fip-bip", help="FIP/BIP concentration in uM."),
    conc_f3_b3: Optional[float] = typer.Option(None, "--conc-f3-b3", help="F3/B3 concentration in uM."),
    conc_loop: Optional[float] = typer.Option(None, "--conc-loop", help="Loop primer concentration in uM."),
    copies: Optional[float] = typer.Option(None, "--copies", help="Target copies per uL."),
    max_variants: int = typer.Option(64, "--max-variants", help="Maximum variants allowed per degenerate primer (default: 64)."),
    separate_targets: bool = typer.Option(False, "--separate-targets", "--per-target", help="Solve targets independently instead of joint multi-target equilibrium."),
    force_joint: bool = typer.Option(False, "--force-joint", help="Force joint analysis of multiple panels for the same target.")
):
    """
    Analyzes a multiplex primer panel and generates a comprehensive thermodynamic diagnostic report.
    """
    if not config and not primers:
        typer.echo("Erreur: Vous devez spécifier au moins un fichier de configuration YAML (-c) ou un fichier d'amorces (-p).")
        raise typer.Exit(code=1)
        
    typer.echo("Loading configuration...")
    start_time = time.time()
    
    # Use config file for hash if available, otherwise primers file
    file_hash = hash_file(config) if config else (hash_file(primers) if primers else "no-file")

    try:
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
            copies=copies
        )
    except FileNotFoundError as e:
        typer.echo(f"Erreur de fichier introuvable : {e.filename}")
        raise typer.Exit(code=1)
    except ParseError as e:
        typer.echo(f"Erreur de lecture de fichier :\n{e}")
        raise typer.Exit(code=1)
    except ValidationError as e:
        typer.echo(f"Erreur de validation de la configuration :\n{e}")
        raise typer.Exit(code=1)

    # 1. Parse Targets
    targets_dict = {}
    for t in config_obj.targets:
        try:
            records = read_multi_fasta(t.sequence_file)
            for h, seq in records:
                t_id = h.split()[0]
                # In build_config_from_cli, we might have multiple targets in one file.
                # If we parsed it via CLI -t, they are all inside config_obj.targets.
                # We can just match the one we need, or collect them all.
                if t_id == t.id or not config: # if not config, we just collect them all
                    targets_dict[t_id] = seq
        except Exception as e:
            typer.echo(f"Erreur lors de la lecture du fichier cible {t.sequence_file}:\n{e}")
            raise typer.Exit(code=1)
        
    if not targets_dict:
        targets_dict = {"Compétition Sans Cible": ""}
        
    has_true_target = any(k != "Compétition Sans Cible" for k in targets_dict.keys())

    target_to_panels = {}
    for pset in config_obj.primer_sets:
        target_to_panels.setdefault(pset.target, []).append(pset)
    
    has_collisions = any(len(panels) > 1 for panels in target_to_panels.values())
    config_list = []
    if has_collisions and not force_joint:
        typer.echo("Attention: Plusieurs panels alternatifs ciblant la même séquence ont été détectés.")
        typer.echo("Analyse séparée des panels...")
        import copy
        for pset in config_obj.primer_sets:
            new_config = copy.deepcopy(config_obj)
            new_config.primer_sets = [pset]
            config_list.append(new_config)
    else:
        config_list.append(config_obj)
        
    for i_cfg, current_config in enumerate(config_list):
        if len(config_list) > 1:
            panel_n = current_config.primer_sets[0].panel_name or str(i_cfg)
            current_output = output.replace('.html', f'_{panel_n}.html')
            typer.echo(f'--- Analyse du panel {panel_n} ---')
        else:
            current_output = output
            
        # 2. & 3. Build Internal Engine Components
        primers, primer_to_panel, backend, backend_kwargs, mon_molar, enzyme, temp_celsius, profiles = build_engine_from_config(current_config, targets_dict, max_variants=max_variants)
        
        if not primers:
            typer.echo("Erreur : Le panel d'amorces est vide. Veuillez définir au moins un jeu d'amorces (primer_sets).")
            raise typer.Exit(code=1)
        
        # 4. Simulation
        target_occupations = {}
        all_fractions = {}
        all_risks = []
        max_residual_global = 0.0
        all_unfolding_penalties = {}
        interaction_matrix = {}
        
        prob = None
        strands = None
        complexes = None
        
        import warnings
        with warnings.catch_warnings(record=True) as captured_warnings:
            warnings.simplefilter("always")
            
            if not separate_targets and len(targets_dict) > 1:
                # Mode CONJOINT (par défaut en multiplexe multi-cibles)
                # JOINT mode (default in multi-target multiplex)
                typer.echo(f"Building joint complex network and solving coupled thermodynamic equilibrium for {len(targets_dict)} targets ({', '.join(targets_dict.keys())})...")
                
                prob, strands, complexes, unfolding_penalties = enumerate_complexes(
                    primers, targets_dict if has_true_target else None, backend, profile=profiles, temp_celsius=temp_celsius, mon_molar=mon_molar, buffer=config_obj.experiment.buffer.model_dump() if config_obj.experiment.buffer else None
                )
                
                res = solve_dual(prob)
                max_residual_global = res.max_residual
                
                # Diagnostics
                amplifiable_flags = []
                concs = []
                dimer_details = []
                
                for i, c_name in enumerate(complexes):
                    stoich = prob.stoichiometry[i]
                    conc = res.concentrations[i]
                    concs.append(conc)
                    
                    is_amp = False
                    details = {}
                    has_target = any(k >= len(primers) and val > 0 for k, val in enumerate(prob.stoichiometry[i]))
                    
                    if not has_target:
                        p_counts = {k: prob.stoichiometry[i][k] for k, val in enumerate(prob.stoichiometry[i]) if val > 0 and k < len(primers)}
                        total_primers = sum(p_counts.values())
                        
                        if total_primers == 2:
                            if len(p_counts) == 1:
                                k = list(p_counts.keys())[0]
                                p_a = primers[k]
                                is_amp, dg_3p, details = evaluate_pair_amplifiable(p_a, p_a, backend, enzyme, temp_celsius, **backend_kwargs)
                            elif len(p_counts) == 2:
                                k1, k2 = list(p_counts.keys())
                                p_a = primers[k1]
                                p_b = primers[k2]
                                is_amp, dg_3p, details = evaluate_pair_amplifiable(p_a, p_b, backend, enzyme, temp_celsius, **backend_kwargs)
                                
                            from labcraft.report.alignment import dotbracket_to_alignment, get_alignment_columns
                            if details and 'alignment' not in details:
                                details['alignment'] = dotbracket_to_alignment(details['seq_a'], details['seq_b'], details['structure'])
                            if details and 'alignment_columns' not in details:
                                cols = get_alignment_columns(details['seq_a'], details['seq_b'], details['structure'], details.get('extensible_strand'))
                                details['alignment_columns'] = cols
                                
                                three_prime_idx = next((i for i, c in enumerate(cols) if c.get('role') == 'three_prime'), -1)
                                template_count = sum(1 for c in cols if c.get('role') == 'template')
                                before_arrow = sum(1 for c in cols[:three_prime_idx] if not c.get('is_truncation')) if three_prime_idx != -1 else 0
                                
                                details['arrow_metrics'] = {
                                    "show": three_prime_idx != -1 and template_count > 0,
                                    "margin_cols": before_arrow,
                                    "width_cols": template_count
                                }
                    amplifiable_flags.append(is_amp)
                    dimer_details.append(details)
                    
                risks = evaluate_risks(complexes, concs, amplifiable_flags, is_warm_start=False, dimer_details=dimer_details)
                fractions = compute_fractions(strands, complexes, prob.stoichiometry, res.free_concentrations, prob.delta_g, temp_celsius, primer_to_panel)
                
                all_risks = risks
                all_fractions = fractions
                
                for t_id in targets_dict.keys():
                    target_occupations[t_id] = {}
                    
                for s in strands:
                    if s.endswith("_site"):
                        idx = strands.index(s)
                        occ = (prob.total_concentrations[idx] - res.free_concentrations[idx]) / prob.total_concentrations[idx] if prob.total_concentrations[idx] > 0 else 0.0
                        if "@" in s:
                            p_part, rest_part = s.split("@", 1)
                            t_part = rest_part.replace("_site", "")
                            if t_part in target_occupations:
                                target_occupations[t_part][s] = occ
                                target_occupations[t_part][f"{p_part}_site"] = occ
                        else:
                            for t_id in targets_dict.keys():
                                target_occupations[t_id][s] = occ
                                
                if isinstance(unfolding_penalties, dict) and any(isinstance(v, dict) for v in unfolding_penalties.values()):
                    all_unfolding_penalties = unfolding_penalties
                else:
                    for t_id in targets_dict.keys():
                        all_unfolding_penalties[t_id] = unfolding_penalties
            else:
                # Mode Séparé (ou cible unique)
                # Separate mode (or single target)
                for t_id, t_seq in targets_dict.items():
                    typer.echo(f"Solving for target {t_id}...")
                    profile = profiles.get(t_id)
                    if profile is None and profiles:
                        profile = next(iter(profiles.values()))
                        
                    prob, strands, complexes, unfolding_penalties = enumerate_complexes(
                        primers, t_seq if t_seq else None, backend, profile=profile, temp_celsius=temp_celsius, mon_molar=mon_molar, buffer=config_obj.experiment.buffer.model_dump() if config_obj.experiment.buffer else None
                    )
                    
                    res = solve_dual(prob)
                    max_residual_global = max(max_residual_global, res.max_residual)
                    
                    # Diagnostics
                    amplifiable_flags = []
                    concs = []
                    dimer_details = []
                    
                    for i, c_name in enumerate(complexes):
                        stoich = prob.stoichiometry[i]
                        conc = res.concentrations[i]
                        concs.append(conc)
                        
                        is_amp = False
                        details = {}
                        has_target = any(k >= len(primers) and val > 0 for k, val in enumerate(prob.stoichiometry[i]))
                        
                        if not has_target:
                            p_counts = {k: prob.stoichiometry[i][k] for k, val in enumerate(prob.stoichiometry[i]) if val > 0 and k < len(primers)}
                            total_primers = sum(p_counts.values())
                            
                            if total_primers == 2:
                                if len(p_counts) == 1:
                                    k = list(p_counts.keys())[0]
                                    p_a = primers[k]
                                    is_amp, dg_3p, details = evaluate_pair_amplifiable(p_a, p_a, backend, enzyme, temp_celsius, **backend_kwargs)
                                elif len(p_counts) == 2:
                                    k1, k2 = list(p_counts.keys())
                                    p_a = primers[k1]
                                    p_b = primers[k2]
                                    is_amp, dg_3p, details = evaluate_pair_amplifiable(p_a, p_b, backend, enzyme, temp_celsius, **backend_kwargs)
                                    
                                from labcraft.report.alignment import dotbracket_to_alignment, get_alignment_columns
                                if details and 'alignment' not in details:
                                    details['alignment'] = dotbracket_to_alignment(details['seq_a'], details['seq_b'], details['structure'])
                                if details and 'alignment_columns' not in details:
                                    cols = get_alignment_columns(details['seq_a'], details['seq_b'], details['structure'], details.get('extensible_strand'))
                                    details['alignment_columns'] = cols
                                    
                                    three_prime_idx = next((i for i, c in enumerate(cols) if c.get('role') == 'three_prime'), -1)
                                    template_count = sum(1 for c in cols if c.get('role') == 'template')
                                    before_arrow = sum(1 for c in cols[:three_prime_idx] if not c.get('is_truncation')) if three_prime_idx != -1 else 0
                                    
                                    details['arrow_metrics'] = {
                                        "show": three_prime_idx != -1 and template_count > 0,
                                        "margin_cols": before_arrow,
                                        "width_cols": template_count
                                    }
                        amplifiable_flags.append(is_amp)
                        dimer_details.append(details)
                        
                    risks = evaluate_risks(complexes, concs, amplifiable_flags, is_warm_start=False, dimer_details=dimer_details)
                    fractions = compute_fractions(strands, complexes, prob.stoichiometry, res.free_concentrations, prob.delta_g, temp_celsius, primer_to_panel)
                    
                    for r in risks:
                        if not any(ar.complex_name == r.complex_name for ar in all_risks):
                            all_risks.append(r)
                            
                    target_occupations[t_id] = {}
                    for s in strands:
                        if s.endswith("_site"):
                            idx = strands.index(s)
                            occ = (prob.total_concentrations[idx] - res.free_concentrations[idx]) / prob.total_concentrations[idx] if prob.total_concentrations[idx] > 0 else 0.0
                            target_occupations[t_id][s] = occ
                            
                    all_fractions.update(fractions)
                    if t_id not in all_unfolding_penalties:
                        all_unfolding_penalties[t_id] = unfolding_penalties
                
            # Matrice d'interaction globale pour la heatmap
            def get_parent(name: str) -> str:
                return name.split('#')[0] if '#' in name else name
                
            unique_parents = list(dict.fromkeys(get_parent(p.name) for p in primers))
            
            for parent1 in unique_parents:
                interaction_matrix[parent1] = {}
                for parent2 in unique_parents:
                    interaction_matrix[parent1][parent2] = 0.0
                    
            if complexes and prob is not None:
                for i_c, c_name in enumerate(complexes):
                    if "_on_" in c_name or "_free" in c_name:
                        continue
                    dg = float(prob.delta_g[i_c])
                    stoich = prob.stoichiometry[i_c]
                    p_indices = [k for k, val in enumerate(stoich) if val > 0 and k < len(primers)]
                    if len(p_indices) == 1:
                        k = p_indices[0]
                        par = get_parent(primers[k].name)
                        if interaction_matrix[par][par] == 0.0 or dg < interaction_matrix[par][par]:
                            interaction_matrix[par][par] = dg
                    elif len(p_indices) == 2:
                        k1, k2 = p_indices
                        par1 = get_parent(primers[k1].name)
                        par2 = get_parent(primers[k2].name)
                        if interaction_matrix[par1][par2] == 0.0 or dg < interaction_matrix[par1][par2]:
                            interaction_matrix[par1][par2] = dg
                            interaction_matrix[par2][par1] = dg
                        
        # 5. Verdict
        def get_parent(name: str) -> str:
            return name.split('#')[0] if '#' in name else name
    
        loop_primer_parents = {get_parent(p.name) for p in primers if p.role in (PrimerRole.LF, PrimerRole.LB)}
        
        verdict = generate_verdict(all_fractions, target_occupations, all_risks, loop_primer_parents=loop_primer_parents, primer_to_panel=primer_to_panel)
        
        # Recommandations
        from labcraft.optimize.recommendations import generate_recommendations
        recommendations = generate_recommendations(verdict, risks=all_risks)
        
        # Contrôle des Sondes TaqMan : déclenché par la présence de sondes, pas par la chimie
        # TaqMan probe check: triggered by presence of probes, not by chemistry field
        probe_tm_results = []
        has_probes = any(p.role == PrimerRole.PROBE for p in primers)
        if has_probes:
            probe_tm_results = check_probes_tm(primers, backend, temp_celsius, **backend_kwargs)
            
        # Optimisation des concentrations
        optimization_results = []
        # On n'optimise que si le génome est présent et qu'il y a un vrai problème
        # Mais le rapport peut toujours afficher les résultats.
        # We need the last `prob`, `complexes`, `strands`
        # Let's run it.
        from labcraft.optimize.concentrations import optimize_concentrations
        free_fractions = {p: f.free for p, f in all_fractions.items()}
        if prob and complexes and strands:
            optimization_results = optimize_concentrations(
                prob_template=prob,
                species_names=strands,
                primers=primers,
                target_dict=targets,
                primer_to_panel=primer_to_panel,
                original_free_fractions=free_fractions,
                original_target_occupations={},
                complex_names=complexes,
                temp_celsius=temp_celsius,
                backend=backend,
                enzyme=enzyme
            )
        
        
        # 5.5 Multiplex Balance & Mispriming
        from labcraft.metrics.balance import calculate_multiplex_balance
        from labcraft.diagnostics.mispriming import detect_inter_target_mispriming
        
        panel_summaries, balance_cv = {}, None
        if target_occupations:
            panel_summaries, balance_cv = calculate_multiplex_balance(primer_to_panel, target_occupations, free_fractions, loop_primer_parents=loop_primer_parents)
        
        mispriming_risks = []
        if config_obj.targets and len(config_obj.targets) > 1:
            target_dict = targets_dict  # targets est le dictionnaire {t.id: t.sequence} déjà chargé
            mispriming_risks = detect_inter_target_mispriming(
                primers, primer_to_panel, target_dict, backend, enzyme, temp_celsius, **backend_kwargs
            )
        
        
        # 6. Generate Report
        typer.echo("Generating HTML report...")
        
        metadata = {
            "version": importlib.metadata.version("labcraft") if importlib.metadata.version else "0.0.1",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
            "file_hash": file_hash,
            "max_residual": max_residual_global,
            "temperature": temp_celsius,
            "buffer": config_obj.experiment.buffer.model_dump() if config_obj.experiment.buffer else "Reference conditions (1.0 M Na+, 0 mM Mg2+)",
            "enzyme": enzyme.name,
            "dimer_dg_threshold": enzyme.dimer_dg_threshold,
            "concentrations_fip_bip": "Per target (check config)",
            "concentrations_target": "Per target (check config)",
            "interaction_matrix": interaction_matrix,
            "primer_names": list(interaction_matrix.keys()),
            "unfolding_penalties": all_unfolding_penalties,
            "target_occupations": target_occupations,
            "has_true_target": has_true_target,
            "primer_to_panel": primer_to_panel,
            "panel_summaries": panel_summaries,
            "balance_cv": balance_cv,
            "mispriming_risks": mispriming_risks,
            "chemistry": config_obj.experiment.chemistry,
            "optimization_results": optimization_results,
            "probe_tm_results": probe_tm_results,
            "recommendations": recommendations,
            "loop_primer_parents": list(loop_primer_parents),
            "primers_parent": primers,
            "multiplex_mode": "Conjoint (joint multi-target equilibrium)" if not separate_targets else "Séparé (cible par cible)",
            "warnings": [str(w.message) for w in captured_warnings]
        }
        
        html = render_report(verdict, all_fractions, all_risks, metadata)
        
        with open(current_output, "w") as f:
            f.write(html)
            
        typer.echo(f"Analysis complete in {time.time() - start_time:.2f}s. Report saved to {current_output}.")
        


@app.command()
def coverage(
    primers: str = typer.Option(..., "-p", "--primers", help="Path to primers FASTA/TXT file."),
    fasta: str = typer.Option(..., "-f", "--fasta", help="Path to multi-FASTA file of strains."),
    sites: Optional[str] = typer.Option(None, "-s", "--sites", help="Path to candidate sites CSV (optional, automatic screening if omitted)."),
    output: str = typer.Option("coverage_report.html", "-o", "--output", help="Output HTML report path."),
    panel: Optional[str] = typer.Option(None, "--panel", help="Panel name to analyze if multiple panels are present in primers file."),
    temperature: float = typer.Option(65.0, "--temperature", help="Experiment temperature in Celsius."),
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
):
    """Analyse de la couverture thermodynamique multi-souches."""
    import csv
    import time
    import tempfile
    from labcraft.cli.parsers import build_config_from_cli, read_multi_fasta
    from labcraft.cli.config import build_engine_from_config
    from labcraft.lamp.coverage import CoverageAnalyzer
    from labcraft.target.seeding import find_candidate_sites
    from labcraft.report.coverage_report import render_coverage_report

    print("Lecture des entrées...")
    fasta_list = read_multi_fasta(fasta)
    fasta_dict = {name: seq for name, seq in fasta_list}
    
    first_strain_path = None
    if fasta_list:
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as tmp:
            tmp.write(f">{fasta_list[0][0]}\n{fasta_list[0][1]}\n")
            first_strain_path = tmp.name

    config_obj = build_config_from_cli(
        config_path=None,
        primers_path=primers,
        targets_path=first_strain_path,
        temperature=temperature,
        na=na,
        mg=mg,
        dntp=dntp,
        allow_unmatched_targets=True
    )
    
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
    physical_primers, _, backend, _, _, enzyme, _, _ = build_engine_from_config(config_obj, dummy_targets)
    
    print(f"Chargement de {len(fasta_dict)} souches.")
    
    csv_records = []
    if sites:
        print(f"Lecture des sites candidats depuis : {sites}")
        with open(sites, 'r') as f:
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
        with open(export_sites, 'w', newline='') as f:
            fieldnames = ["strain_id", "primer_role", "primer_name", "position", "strand", "site_seq", "n_mismatches", "panel"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for r in csv_records:
                writer.writerow({k: r.get(k, "") for k in fieldnames})
        print(f"Sites candidats exportés dans : {export_sites}")

    analyzer = CoverageAnalyzer(
        physical_primers, fasta_dict, backend, enzyme,
        temp_celsius=temperature,
        ddg_max=ddg_max,
        dg_threshold=dg_threshold,
        max_mismatches_count=max_mismatches
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


if __name__ == "__main__":
    app()
