import csv
import matplotlib.pyplot as plt
import numpy as np

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
    primers, base_concs = load_primers("validation/reference_data/DENV2_primers.csv")
    
    yaml_data = {
        "experiment_name": "DENV2 LAVA 72plex",
        "description": "Validation sur les amorces DENV2 de Lopez-Jimena 2018",
        "temperature_celsius": 65.0,
        "primers": [
            {
                "name": p.name,
                "sequence": p.sequence,
                "role": p.role.value
            } for p in primers
        ],
        "targets": []
    }
    
    with open("validation/DENV2.yaml", "w") as f:
        yaml.dump(yaml_data, f, sort_keys=False)
        
    import os
    os.system("python3 src/labcraft/cli/main.py validation/DENV2.yaml -o validation/report_DENV2.html")

