import csv
import matplotlib.pyplot as plt
import numpy as np
import hashlib

from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.thermo.backends.vienna import ViennaRNABackend
from labcraft.pipeline.titration import simulate_titration

def load_primers(csv_path: str) -> tuple[list[PhysicalPrimer], dict[str, float]]:
    primers = []
    concentrations_1x = {}
    
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            role_str = row["role"]
            # Map FLOOP/BLOOP to LF/LB
            if role_str == "FLOOP":
                role_str = "LF"
            elif role_str == "BLOOP":
                role_str = "LB"
                
            p = PhysicalPrimer(
                name=row["name"],
                sequence=row["sequence"],
                role=PrimerRole(role_str),
                binding_domain=row["sequence"]
            )
            primers.append(p)
            concentrations_1x[p.name] = float(row["conc_std_uM"]) * 1e-6
            
    return primers, concentrations_1x

def run_validation(name: str, csv_path: str, out_svg: str):
    print(f"--- Validation {name} ---")
    primers, base_concs = load_primers(csv_path)
    print(f"Loaded {len(primers)} primers.")
    
    backend = ViennaRNABackend()
    dilutions = [1.0, 0.5, 0.25, 0.125]
    
    res = simulate_titration(
        primers,
        backend,
        base_concs,
        dilutions=dilutions,
        temp_celsius=65.0
    )
    
    idx_1x = dilutions.index(1.0)
    idx_1to4 = dilutions.index(0.25)
    
    free_1x = res.free_fractions[idx_1x, :]
    free_1to4 = res.free_fractions[idx_1to4, :]
    
    print(f"Moyenne libre à 1x: {np.mean(free_1x):.1f}%")
    print(f"Amorces < 20% à 1x: {np.sum(free_1x < 20.0)}")
    
    print(f"Moyenne libre à 1/4x: {np.mean(free_1to4):.1f}%")
    print(f"Amorces < 20% à 1/4x: {np.sum(free_1to4 < 20.0)}")
    
    print("\nClassement des amorces séquestrées (Top pires à 1x) :")
    sorted_indices = np.argsort(free_1x)
    for i in sorted_indices:
        if free_1x[i] < 20.0:
            p_name = res.primer_names[i]
            dom = res.dominant_complexes.get(p_name, "Inconnu")
            print(f"- {p_name} : {free_1x[i]:.1f}% libre (Séquestrée par {dom})")
            
    # Plotting
    plt.figure(figsize=(10, 6))
    
    colors = plt.cm.tab10.colors
    color_idx = 0
    
    for i in range(len(res.primer_names)):
        if free_1x[i] < 20.0:
            # Séquestrée = couleur + épais
            c = colors[color_idx % len(colors)]
            plt.plot(dilutions, res.free_fractions[:, i], color=c, linewidth=2.5, alpha=0.9, marker='o', label=res.primer_names[i])
            color_idx += 1
        else:
            # Libre = gris fin
            plt.plot(dilutions, res.free_fractions[:, i], color='gray', linewidth=1.0, alpha=0.3)
            
    plt.xscale('log', base=2)
    plt.xticks(dilutions, [f"1/{int(1/d)}" if d < 1 else "1x" for d in dilutions])
    plt.gca().invert_xaxis() # Pour avoir 1x à gauche, puis 1/2, 1/4, 1/8
    
    plt.title(f"Titration de la fraction libre - {name} ({len(primers)} amorces)")
    plt.xlabel("Facteur de dilution")
    plt.ylabel("Fraction libre de l'amorce (%)")
    plt.ylim(0, 105)
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # Placer la légende s'il y a des amorces séquestrées
    if color_idx > 0:
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(out_svg)
    print(f"Graphique sauvegardé: {out_svg}\n")
    plt.close()

if __name__ == "__main__":
    run_validation("DENV4 (18plex)", "validation/reference_data/DENV4_primers.csv", "validation/titration_DENV4.svg")
    run_validation("DENV2 (72plex)", "validation/reference_data/DENV2_primers.csv", "validation/titration_DENV2.svg")

    # Génération du rapport HTML complet pour DENV2
    print("\nGénération du rapport HTML pour DENV2 (avec heatmap 72x72)...")
    
    import yaml
    from labcraft.lamp.complex_enumeration import enumerate_complexes
    from labcraft.thermo.backends.vienna import ViennaRNABackend
    from labcraft.lamp.stoichiometry import ConcentrationProfile
    from labcraft.solver.dual import solve_dual
    from labcraft.metrics.fractions import compute_fractions
    from labcraft.metrics.verdict import generate_verdict
    from labcraft.report.renderer import render_report
    import dataclasses

    primers, base_concs = load_primers("validation/reference_data/DENV2_primers.csv")
    
    backend = ViennaRNABackend()
    profile = ConcentrationProfile(target=0.0, fip_bip=1.6, f3_b3=0.2, lf_lb=0.8)
    prob, strand_names, complex_names, _ = enumerate_complexes(primers, "", backend, profile, 65.0)
    c_tot_arr = np.zeros(prob.n_strands)
    for i, s in enumerate(strand_names):
        if s in base_concs:
            c_tot_arr[i] = base_concs[s]
    prob = dataclasses.replace(prob, total_concentrations=c_tot_arr)
    
    sol = solve_dual(prob)
    fractions = compute_fractions(
        primer_names=[p.name for p in primers],
        complex_names=complex_names,
        stoichiometry=prob.stoichiometry,
        free_concentrations=sol.free_concentrations,
        delta_g=prob.delta_g,
        temp_celsius=65.0
    )
    risks = []
    verdict = generate_verdict(fractions, {}, risks)
    
    interaction_matrix = {}
    for p1 in primers:
        interaction_matrix[p1.name] = {}
        for p2 in primers:
            if p1.name == p2.name:
                c_name = f"{p1.name}_homo"
            else:
                p1_idx = primers.index(p1)
                p2_idx = primers.index(p2)
                c_name = f"{primers[min(p1_idx, p2_idx)].name}_{primers[max(p1_idx, p2_idx)].name}"
            
            try:
                c_idx = complex_names.index(c_name)
                interaction_matrix[p1.name][p2.name] = prob.delta_g[c_idx]
            except ValueError:
                interaction_matrix[p1.name][p2.name] = 0.0

    import time
    
    with open("validation/reference_data/DENV2_primers.csv", "rb") as f:
        file_hash = hashlib.sha256(f.read()).hexdigest()
        
    metadata = {
        "experiment_name": "DENV2 LAVA 72plex",
        "description": "Validation sur les amorces DENV2 de Lopez-Jimena 2018",
        "temperature_C": 65.0,
        "temperature": 65.0,
        "buffer": "Reference conditions (1.0 M Na+, 0 mM Mg2+)",
        "enzyme": "Bst 2.0 WarmStart",
        "dimer_dg_threshold": -2.0,
        "primer_names": [p.name for p in primers],
        "interaction_matrix": interaction_matrix,
        "max_residual": sol.max_residual,
        "unfolding_penalties": {},
        "file_hash": file_hash,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        "version": "0.0.1"
    }
    
    html = render_report(verdict, fractions, risks, metadata)
    
    with open("validation/report_DENV2.html", "w") as f:
        f.write(html)
        
    print("Rapport enregistré dans validation/report_DENV2.html")

