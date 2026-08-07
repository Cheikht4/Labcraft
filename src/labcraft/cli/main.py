import typer
import yaml
import hashlib
import time
import importlib.metadata
import os

from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.lamp.stoichiometry import ConcentrationProfile, target_copies_to_molar
from labcraft.lamp.complex_enumeration import enumerate_complexes
from labcraft.thermo.backends.vienna import ViennaRNABackend
from labcraft.solver.dual import solve_dual
from labcraft.metrics.fractions import compute_fractions
from labcraft.diagnostics.enzyme import BST_2_0, TAQ, PolymeraseProfile
from labcraft.diagnostics.amplifiable_dimer import is_amplifiable_dimer
from labcraft.metrics.risk import evaluate_risks
from labcraft.metrics.verdict import generate_verdict
from labcraft.report.renderer import render_report

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

    # 1. Parse Targets
    targets = {}
    for t in data.get("targets", []):
        seq = read_fasta(t["sequence_file"])
        targets[t["id"]] = seq

    # 2. Parse Primers & Panel Mapping
    primers = []
    primer_to_panel = {}
    
    # Simple hardcoded parser for the YAML structure (In production, use Pydantic models)
    for pset in data.get("primer_sets", []):
        t_id = pset["target"]
        p_dict = pset["primers"]
        
        for role_name, p_data in p_dict.items():
            name = f"{role_name}_{t_id.replace('Synth', '')}" # e.g. F3_A
            seq = p_data["seq"]
            
            primer_to_panel[name] = t_id
            
            if role_name in ("FIP", "BIP") and "domains" in p_data and p_data["domains"] != "auto":
                d = p_data["domains"]
                primers.append(PhysicalPrimer(name, seq, PrimerRole[role_name], d["F2"], d["F1c"], d.get("linker", "")))
            elif role_name in ("FIP", "BIP"):
                primers.append(PhysicalPrimer.from_alignment(name, seq, PrimerRole[role_name], targets[t_id]))
            else:
                try:
                    role_enum = PrimerRole[role_name]
                except KeyError:
                    role_enum = PrimerRole.F3
                primers.append(PhysicalPrimer(name, seq, role_enum, seq))
                
    # 3. Setup Thermodynamic Engine
    temp_celsius = data.get("experiment", {}).get("temperature_C", 65.0)
    backend = ViennaRNABackend()
    
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
        # Enumerate with the full primer multiplex but only one target
        profile = ConcentrationProfile(
            target=target_copies_to_molar(1000), 
            fip_bip=1.6e-6, 
            f3_b3=0.2e-6, 
            lf_lb=0.8e-6
        )
        
        prob, strands, complexes, unfolding_penalties = enumerate_complexes(
            primers, t_seq, backend, profile=profile, temp_celsius=temp_celsius
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
                    res_homo = backend.calc_homodimer(p_a_seq, temp_celsius=temp_celsius)
                    struct, mfe = res_homo.structure, res_homo.dg_kcal
                    is_amp, dg_3p = is_amplifiable_dimer(p_a_seq, p_a_seq, struct, mfe, BST_2_0, temp_celsius)
                    from labcraft.report.alignment import dotbracket_to_alignment
                    details = {
                        "seq_a": p_a_seq, "seq_b": p_a_seq, "structure": struct, 
                        "delta_g": mfe, "delta_g_3p": dg_3p,
                        "alignment": dotbracket_to_alignment(p_a_seq, p_a_seq, struct)
                    }
                elif len(parts) >= 4:
                    p2_name = "_".join(parts[2:4])
                    try:
                        p_a_seq = next(p.sequence for p in primers if p.name == p1_name)
                        p_b_seq = next(p.sequence for p in primers if p.name == p2_name)
                        res_hetero = backend.calc_heterodimer(p_a_seq, p_b_seq, temp_celsius=temp_celsius)
                        struct, mfe = res_hetero.structure, res_hetero.dg_kcal
                        is_amp, dg_3p = is_amplifiable_dimer(p_a_seq, p_b_seq, struct, mfe, BST_2_0, temp_celsius)
                        from labcraft.report.alignment import dotbracket_to_alignment
                        details = {
                            "seq_a": p_a_seq, "seq_b": p_b_seq, "structure": struct, 
                            "delta_g": mfe, "delta_g_3p": dg_3p,
                            "alignment": dotbracket_to_alignment(p_a_seq, p_b_seq, struct)
                        }
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
                        res = backend.calc_homodimer(p1.sequence, temp_celsius=temp_celsius)
                    else:
                        res = backend.calc_heterodimer(p1.sequence, p2.sequence, temp_celsius=temp_celsius)
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
