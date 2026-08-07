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
from labcraft.metrics.risk import evaluate_risks
from labcraft.metrics.verdict import generate_verdict
from labcraft.report.renderer import render_report
from labcraft.cli.config import PanelConfig, build_engine_from_config
from pydantic import ValidationError

app = typer.Typer(help="LabCraft CLI - Thermodynamic multiplex primer panel design engine.")

def read_fasta(filepath: str) -> str:
    with open(filepath, 'r') as f:
        lines = f.read().splitlines()
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    return seq

def hash_file(filepath: str) -> str:
    hasher = hashlib.sha256()
    with open(filepath, 'rb') as f:
        buf = f.read()
        hasher.update(buf)
    return hasher.hexdigest()

@app.command()
def analyze(
    config: str = typer.Argument(..., help="Path to the YAML panel configuration file."),
    output: str = typer.Option("report.html", "-o", "--output", help="Output HTML report path.")
):
    """
    Analyzes a multiplex primer panel and generates a comprehensive thermodynamic diagnostic report.
    """
    typer.echo(f"Loading configuration from {config}...")
    
    start_time = time.time()
    file_hash = hash_file(config)
    
    with open(config, "r") as f:
        data = yaml.safe_load(f)

    try:
        config_obj = PanelConfig.model_validate(data)
    except ValidationError as e:
        typer.echo(f"Erreur de validation de la configuration :\n{e}")
        raise typer.Exit(code=1)

    # 1. Parse Targets
    targets = {}
    for t in config_obj.targets:
        seq = read_fasta(t.sequence_file)
        targets[t.id] = seq

    # 2. & 3. Build Internal Engine Components
    primers, primer_to_panel, backend, backend_kwargs, mon_molar, enzyme, temp_celsius, profiles = build_engine_from_config(config_obj, targets)
    
    # 4. Simulation
    typer.echo("Building complex network and solving thermodynamic equilibrium...")
    
    # Target occupation dictionary to collect
    target_occupations = {}
    all_fractions = {}
    all_risks = []
    max_residual_global = 0.0
    all_unfolding_penalties = {}
    
    # Interaction matrix for heatmap
    interaction_matrix = {}
    
    # For now, we simulate everything together against each target to get the target occupation.
    # Actually, a full multiplex solve computes everything in one matrix!
    # To keep it simple, we solve the full matrix with Target A, then Target B, and aggregate.
    
    results_by_target = {}
    
    for t_id, t_seq in targets.items():
        typer.echo(f"Solving for target {t_id}...")
        prob, strands, complexes, unfolding_penalties = enumerate_complexes(
            primers, t_seq, backend, profile=profiles[t_id], temp_celsius=temp_celsius, mon_molar=mon_molar, buffer=config_obj.experiment.buffer.model_dump() if config_obj.experiment.buffer else None
        )
        
        res = solve_dual(prob)
        max_residual_global = max(max_residual_global, res.max_residual)
        
        # Diagnostics
        R = 0.00198720425864083
        RT = R * (273.15 + temp_celsius)
        import numpy as np
        u = np.log(res.free_concentrations)
        
        amplifiable_flags = []
        concs = []
        dimer_details = []
        
        for i, c_name in enumerate(complexes):
            stoich = prob.stoichiometry[i]
            conc = np.exp(-prob.delta_g[i] / RT + np.dot(stoich, u))
            concs.append(conc)
            
            is_amp = False
            details = {}
            if "_on_" not in c_name and "_free" not in c_name:
                parts = c_name.split('_')
                p1_name = "_".join(parts[:2]) if len(parts) >= 2 else parts[0]
                if "homo" in c_name:
                    p_a_seq = next(p.sequence for p in primers if p.name == p1_name)
                    is_amp, dg_3p, details = evaluate_pair_amplifiable(p_a_seq, p_a_seq, backend, enzyme, temp_celsius, **backend_kwargs)
                    from labcraft.report.alignment import dotbracket_to_alignment
                    if 'alignment' not in details:
                        details['alignment'] = dotbracket_to_alignment(details['seq_a'], details['seq_b'], details['structure'])
                elif len(parts) >= 4:
                    p2_name = "_".join(parts[2:4])
                    try:
                        p_a_seq = next(p.sequence for p in primers if p.name == p1_name)
                        p_b_seq = next(p.sequence for p in primers if p.name == p2_name)
                        is_amp, dg_3p, details = evaluate_pair_amplifiable(p_a_seq, p_b_seq, backend, enzyme, temp_celsius, **backend_kwargs)
                        from labcraft.report.alignment import dotbracket_to_alignment
                        if 'alignment' not in details:
                            details['alignment'] = dotbracket_to_alignment(details['seq_a'], details['seq_b'], details['structure'])
                    except StopIteration:
                        pass
            amplifiable_flags.append(is_amp)
            dimer_details.append(details)
            
        risks = evaluate_risks(complexes, concs, amplifiable_flags, is_warm_start=False, dimer_details=dimer_details)
        fractions = compute_fractions(strands, complexes, prob.stoichiometry, res.free_concentrations, prob.delta_g, temp_celsius, primer_to_panel)
        
        # Update globals
        # Note: In a real system, the solver would include ALL targets simultaneously.
        # We merge risks (take the max concentration).
        for r in risks:
            if not any(ar.complex_name == r.complex_name for ar in all_risks):
                all_risks.append(r)
                
        # Target occupation
        target_occupations[t_id] = {}
        for s in strands:
            if s.endswith("_site"):
                idx = strands.index(s)
                occ = (prob.total_concentrations[idx] - res.free_concentrations[idx]) / prob.total_concentrations[idx]
                target_occupations[t_id][s] = occ
                
        all_fractions.update(fractions)
        if t_id not in all_unfolding_penalties:
            all_unfolding_penalties[t_id] = unfolding_penalties
            
        # Extract interaction matrix (only need to do it once, as it's target-independent for dimer interactions)
        if not interaction_matrix:
            for p1 in primers:
                interaction_matrix[p1.name] = {}
                for p2 in primers:
                    if p1.name == p2.name:
                        res = backend.calc_homodimer(p1.sequence, temp_celsius=temp_celsius, **backend_kwargs)
                    else:
                        res = backend.calc_heterodimer(p1.sequence, p2.sequence, temp_celsius=temp_celsius, **backend_kwargs)
                    interaction_matrix[p1.name][p2.name] = res.dg_kcal

    # 5. Verdict
    verdict = generate_verdict(all_fractions, target_occupations, all_risks)
    
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
        "primer_names": [p.name for p in primers],
        "unfolding_penalties": all_unfolding_penalties
    }
    
    html = render_report(verdict, all_fractions, all_risks, metadata)
    
    with open(output, "w") as f:
        f.write(html)
        
    typer.echo(f"Analysis complete in {time.time() - start_time:.2f}s. Report saved to {output}.")
    
if __name__ == "__main__":
    app()
