#!/usr/bin/env python3
"""Run the synthetic 2-plex LAMP panel and demonstrate the artifacts.

Exécute le panel synthétique et démontre les trois artefacts.
"""
import yaml
import sys
import os

from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.lamp.stoichiometry import ConcentrationProfile, target_copies_to_molar
from labcraft.lamp.complex_enumeration import enumerate_complexes
from labcraft.thermo.backends.vienna import ViennaRNABackend
from labcraft.solver.dual import solve_dual

def read_fasta(filepath):
    with open(filepath, 'r') as f:
        lines = f.read().splitlines()
    seq = ""
    for line in lines:
        if not line.startswith(">"):
            seq += line.strip()
    return seq

def main():
    yaml_path = "examples/synthetic_2plex.yaml"
    with open(yaml_path, "r") as f:
        data = yaml.safe_load(f)

    # 1. Parser les cibles
    targets = {}
    for t in data["targets"]:
        seq = read_fasta(t["sequence_file"])
        targets[t["id"]] = seq

    # 2. Parser les amorces
    primers = []
    # Panel A
    p_a = data["primer_sets"][0]["primers"]
    primers.append(PhysicalPrimer("F3_A", p_a["F3"]["seq"], PrimerRole.F3, p_a["F3"]["seq"]))
    # FIP_A domains: F1c is first 20, linker TTTT, F2 is last 20
    fip_a_seq = p_a["FIP"]["seq"]
    primers.append(PhysicalPrimer("FIP_A", fip_a_seq, PrimerRole.FIP, p_a["FIP"]["domains"]["F2"], fip_a_seq[:20], "TTTT"))
    bip_a_seq = p_a["BIP"]["seq"]
    # BIP auto detect or explicit:
    # TTGCAATGCAATGCAATGCATTTTGACCTTAGGCATTCGAACT (43 nt)
    # F2 size ~20, F1c ~20. Wait, for BIP_A I will just say binding_domain is the last 20.
    primers.append(PhysicalPrimer.from_alignment("BIP_A", bip_a_seq, PrimerRole.BIP, targets["SynthA"]))

    # Panel B
    p_b = data["primer_sets"][1]["primers"]
    primers.append(PhysicalPrimer("F3_B", p_b["F3"]["seq"], PrimerRole.F3, p_b["F3"]["seq"]))
    fip_b_seq = p_b["FIP"]["seq"]
    primers.append(PhysicalPrimer("FIP_B", fip_b_seq, PrimerRole.FIP, p_b["FIP"]["domains"]["F2"], fip_b_seq[:20], "TTTT"))
    bip_b_seq = p_b["BIP"]["seq"]
    primers.append(PhysicalPrimer.from_alignment("BIP_B", bip_b_seq, PrimerRole.BIP, targets["SynthB"]))

    # 3. Paramétrer l'expérience
    profile = ConcentrationProfile(
        target=target_copies_to_molar(1000), 
        fip_bip=1.6e-6, 
        f3_b3=0.2e-6, 
        lf_lb=0.8e-6
    )
    backend = ViennaRNABackend()

    print("=========================================================")
    print(" DÉMONSTRATION DES ARTEFACTS DU PANEL 2-PLEX SYNTHÉTIQUE")
    print("=========================================================\n")

    # Résolution globale sur SynthA (avec tous les primers du multiplex)
    # On ajoute la cible A
    prob_a, strands_a, complexes_a = enumerate_complexes(
        primers, targets["SynthA"], backend, profile=profile, temp_celsius=65.0
    )
    res_a = solve_dual(prob_a)

    print("--- 1. STEM-LOOP SUR F3_A (Blocage thermodynamique) ---")
    # Recherchons l'occupation de F3_A_site
    # Nom du site: F3_A_site
    f3a_site_idx = strands_a.index("F3_A_site")
    f3a_site_free = res_a.free_concentrations[f3a_site_idx]
    f3a_site_tot = prob_a.total_concentrations[f3a_site_idx]
    occupation_a = (f3a_site_tot - f3a_site_free) / f3a_site_tot
    print(f"Occupation du site cible par F3_A : {occupation_a*100:.4f}%")
    if occupation_a < 0.1:
        print("-> Artefact confirmé : le site est enfoui, F3_A ne peut pas s'hybrider de manière compétitive.")
    print("")

    # Affichons les complexes majeurs (artefacts dimériques)
    print("--- 2. DIMÈRES ET HYBRIDATIONS CROISÉES MAJEURES ---")
    # Trier par concentration
    complex_concs = []
    for i, name in enumerate(complexes_a):
        stoich = prob_a.stoichiometry[i]
        # concentration = exp(-dg/RT) * prod(free_i ** stoich_i)
        # Mais le solveur nous donne-t-il les concentrations des complexes ?
        # x_i = exp(-dg/RT + stoich \cdot u)
        import numpy as np
        R = 0.00198720425864083
        RT = R * (273.15 + 65.0)
        u = np.log(res_a.free_concentrations)
        conc = np.exp(-prob_a.delta_g[i] / RT + np.dot(stoich, u))
        if "_free" not in name and conc > 1e-12: # Ne pas afficher les monomères libres
            complex_concs.append((name, prob_a.delta_g[i], conc))

    complex_concs.sort(key=lambda x: x[2], reverse=True)
    for name, dg, conc in complex_concs[:5]:
        print(f"Complexe : {name:20s} | ΔG° : {dg:6.2f} kcal/mol | Conc : {conc:.2e} M")
    
    # Vérifier la présence de FIP_A_BIP_B (cross-hybridization)
    fip_a_bip_b = any(n == "FIP_A_BIP_B" or n == "BIP_B_FIP_A" for n, dg, conc in complex_concs)
    if fip_a_bip_b:
        print("-> Artefact confirmé : forte hybridation croisée entre FIP_A et BIP_B.")
        
    # Vérifier la présence de FIP_B_homo ou FIP_B_BIP_B
    fip_b_homo = any("FIP_B_homo" in n for n, dg, conc in complex_concs)
    fip_b_bip_b = any((n == "FIP_B_BIP_B" or n == "BIP_B_FIP_B") for n, dg, conc in complex_concs)
    if fip_b_homo or fip_b_bip_b:
        print("-> Artefact confirmé : dimères impliquant FIP_B formés à l'équilibre.")

if __name__ == "__main__":
    main()
