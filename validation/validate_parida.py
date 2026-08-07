import csv
import itertools
from collections import defaultdict
import RNA
import hashlib

from labcraft.diagnostics.amplifiable_dimer import is_amplifiable_dimer
from labcraft.diagnostics.enzyme import BST
from labcraft.thermo.vienna import dna_params
from labcraft.metrics.risk import RiskItem
from labcraft.metrics.verdict import PanelVerdict
from labcraft.report.renderer import render_report

def load_parida_csv(path):
    serotypes = defaultdict(list)
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            serotypes[row['serotype']].append(row)
    return serotypes

def format_alignment(seq_a, seq_b, structure):
    # Un format basique pour l'affichage ascii
    return f"{seq_a}\n{structure}\n{seq_b}"

def run_analysis(serotype, primers, skip_flp=False):
    if skip_flp:
        primers = [p for p in primers if p['role'] != 'FLP']
        
    amplifiable_dimers = []
    
    with dna_params(63.0):
        for i in range(len(primers)):
            for j in range(i, len(primers)):
                p1 = primers[i]
                p2 = primers[j]
                seq1 = p1['sequence']
                seq2 = p2['sequence']
                
                # cofold
                seq_concat = f"{seq1}&{seq2}"
                structure, mfe = RNA.cofold(seq_concat)
                struct_clean = structure.replace('&', '')
                
                is_amp, min_dg_3p = is_amplifiable_dimer(
                    seq1, seq2, struct_clean, mfe, BST, temp_celsius=63.0
                )
                
                if is_amp:
                    amplifiable_dimers.append({
                        'primer1': p1['name'],
                        'primer2': p2['name'],
                        'seq1': seq1,
                        'seq2': seq2,
                        'structure': structure,
                        'mfe': mfe,
                        'dg_3p': min_dg_3p,
                        'concentration': max(float(p1['conc_uM']), float(p2['conc_uM'])) * 1e-6
                    })
                    
    # Sort by dg_3p (most negative first)
    amplifiable_dimers.sort(key=lambda x: x['dg_3p'])
    return amplifiable_dimers

def main():
    path = "validation/reference_data/parida2005_lamp_primers.csv"
    serotypes = load_parida_csv(path)
    
    results = {}
    print("--- Validation Parida 2005 ---")
    for s_name, primers in serotypes.items():
        dimers = run_analysis(s_name, primers)
        results[s_name] = dimers
        
        strongest = dimers[0] if dimers else None
        strongest_str = f"{strongest['primer1']}/{strongest['primer2']} (dG 3' = {strongest['dg_3p']:.2f})" if strongest else "Aucun"
        print(f"{s_name}: {len(dimers)} dimères amplifiables. Pire: {strongest_str}")
        if s_name == "DEN-3":
            for d in dimers:
                print(f"  - {d['primer1']}/{d['primer2']}: dG 3' = {d['dg_3p']:.2f}")
        
    # Table récapitulative
    print("\n--- Récapitulatif ---")
    sorted_s = sorted(results.keys(), key=lambda k: len(results[k]), reverse=True)
    table_html = "<table border='1'><tr><th>Sérotype</th><th>Nb Dimères Amplifiables</th><th>Pire Dimère</th><th>ΔG 3' Pire</th></tr>"
    for s in sorted_s:
        dimers = results[s]
        nb = len(dimers)
        pire = f"{dimers[0]['primer1']}/{dimers[0]['primer2']}" if dimers else "-"
        dg = f"{dimers[0]['dg_3p']:.2f}" if dimers else "-"
        print(f"{s:<6} | {nb} dimères | {pire:<20} | {dg}")
        table_html += f"<tr><td>{s}</td><td>{nb}</td><td>{pire}</td><td>{dg}</td></tr>"
    table_html += "</table>"
    
    # Contrôle sans FLP
    dimers_no_flp = run_analysis("DEN-3", serotypes["DEN-3"], skip_flp=True)
    print(f"\nContrôle DEN-3 sans FLP: {len(dimers_no_flp)} dimères amplifiables.")
    
    # Génération du rapport DEN-3
    risks = []
    for d in results["DEN-3"]:
        risks.append(RiskItem(
            complex_name=f"{d['primer1']}_{d['primer2']}",
            concentration=d['concentration'],
            severity=10.0 + max(0.0, -(d['dg_3p'] + 3.0)),
            description=f"Dimère amplifiable (3' extensible).",
            seq_a=d['seq1'],
            seq_b=d['seq2'],
            structure=d['structure'],
            delta_g=d['mfe'],
            delta_g_3p=d['dg_3p'],
            alignment_ascii=format_alignment(d['seq1'], d['seq2'], d['structure'])
        ))
        
    verdict = PanelVerdict("FAILURE", [], "Présence de dimères amplifiables multiples, notamment avec FLP.")
    
    import time
    
    with open(path, "rb") as f:
        file_hash = hashlib.sha256(f.read()).hexdigest()
        
    metadata = {
        "experiment_name": "Parida 2005 - DEN-3",
        "description": f"Validation sur le jeu DEN-3. Limite: le tampon de Parida contient de la bétaïne qui atténue les structures GC (non modélisé).<br><br><b>Comparaison des 4 sérotypes:</b><br>{table_html}",
        "temperature_C": 63.0,
        "primer_names": [p['name'] for p in serotypes["DEN-3"]],
        "interaction_matrix": {},
        "max_residual": 0.0,
        "unfolding_penalties": {},
        "file_hash": file_hash,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        "version": "0.0.1"
    }
    
    html = render_report(verdict, {}, risks, metadata)
    
    with open("validation/report_Parida_DEN3.html", "w") as f:
        f.write(html)
    print("Rapport enregistré dans validation/report_Parida_DEN3.html")

if __name__ == "__main__":
    main()
