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
from labcraft.metrics.fractions import compute_fractions
from labcraft.diagnostics.enzyme import BST_2_0
from labcraft.diagnostics.amplifiable_dimer import is_amplifiable_dimer
from labcraft.metrics.risk import evaluate_risks
from labcraft.metrics.verdict import generate_verdict
import numpy as np

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
        print("-> Artefact confirmé : le site F3_A est enfoui dans la cible, empêchant la liaison.")
    
    print("\n--- 2. DIAGNOSTIC DES ARTEFACTS (DIMÈRES) ---")
    R = 0.00198720425864083
    RT = R * (273.15 + 65.0)
    u = np.log(res_a.free_concentrations)
    
    amplifiable_flags = []
    concs = []
    
    for i, name in enumerate(complexes_a):
        stoich = prob_a.stoichiometry[i]
        conc = np.exp(-prob_a.delta_g[i] / RT + np.dot(stoich, u))
        concs.append(conc)
        
        # Test amplifiabilité
        is_amp = False
        if "_on_" not in name and "_free" not in name:
            # C'est un dimère.
            # Pour l'instant, on lance un avertissement simple (simplification d'usage pour la démo)
            # Dans le vrai cas, on devrait récupérer la structure exacte. 
            # Comme enumerate_complexes ne renvoie pas la structure, on refait un cofold rapide.
            parts = name.split('_')
            # Extraire les noms des amorces
            p1_name = "_".join(parts[:2]) if len(parts) >= 2 else parts[0]
            if "homo" in name:
                p_a_seq = next(p.sequence for p in primers if p.name == p1_name)
                struct, mfe = backend.calc_homodimer(p_a_seq, temp_celsius=65.0).structure, backend.calc_homodimer(p_a_seq, temp_celsius=65.0).dg_kcal
                is_amp, _ = is_amplifiable_dimer(p_a_seq, p_a_seq, struct, mfe, BST_2_0)
            elif len(parts) >= 4: # Heterodimer ex: FIP_A_BIP_B
                p2_name = "_".join(parts[2:4])
                try:
                    p_a_seq = next(p.sequence for p in primers if p.name == p1_name)
                    p_b_seq = next(p.sequence for p in primers if p.name == p2_name)
                    struct, mfe = backend.calc_heterodimer(p_a_seq, p_b_seq, temp_celsius=65.0).structure, backend.calc_heterodimer(p_a_seq, p_b_seq, temp_celsius=65.0).dg_kcal
                    is_amp, _ = is_amplifiable_dimer(p_a_seq, p_b_seq, struct, mfe, BST_2_0)
                except StopIteration:
                    pass
        
        amplifiable_flags.append(is_amp)
        
    risks = evaluate_risks(complexes_a, concs, amplifiable_flags, is_warm_start=False)
    for risk in risks[:5]:
        print(f"Complexe : {risk.complex_name:20s} | Conc : {risk.concentration:.2e} M | Verdict : {risk.description}")

    print("\n--- 3. MÉTRIQUES ET DÉCOMPOSITION (LOI DE CONSERVATION) ---")
    fractions = compute_fractions(strands_a, complexes_a, prob_a.stoichiometry, res_a.free_concentrations, prob_a.delta_g, 65.0)
    
    for p_name in ["F3_A", "FIP_A", "BIP_A", "FIP_B", "BIP_B"]:
        if p_name in fractions:
            f = fractions[p_name]
            print(f"Amorce {p_name:6s} | Libre: {f.free*100:5.2f}% | Dimères: {(f.homodimer+f.heterodimer)*100:5.2f}% | Somme: {f.sum*100:6.2f}%")
            
    verdict = generate_verdict(fractions, risks)
    print("\n--- 4. VERDICT DU PANEL ---")
    print(f"Statut : {verdict.status}")
    print(f"Cause  : {verdict.dominant_cause}")

if __name__ == "__main__":
    main()
