import sys
import yaml
from labcraft.lamp.complex_enumeration import enumerate_complexes
from labcraft.thermo.salt import UnifiedSaltModel
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.diagnostics.enzyme import BST_2_0
from labcraft.optimize.concentrations import optimize_concentrations
from labcraft.solver.dual import solve_dual
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole

def main():
    with open('den3_test.yaml', 'r') as f:
        config = yaml.safe_load(f)
        
    primers = []
    # parse from config
    for pset in config.get('primer_sets', []):
        if pset['target'] == 'DEN3_M93130':
            for role_str, pdata in pset['primers'].items():
                name = f"{role_str}_DEN-3"
                seq = pdata['seq']
                conc = pdata.get('conc_uM', 1.0) * 1e-6
                role_map = {'F3': PrimerRole.F3, 'B3': PrimerRole.B3, 'FIP': PrimerRole.FIP, 'BIP': PrimerRole.BIP, 'LF': PrimerRole.LF, 'LB': PrimerRole.LB}
                role = role_map.get(role_str, PrimerRole.FIP)
                if role in (PrimerRole.F3, PrimerRole.LF):
                    conc *= 20.0
                

                
                bd = seq
                if 'domains' in pdata:
                    # Find the longest domain as binding domain, or use known keys
                    if role == PrimerRole.FIP:
                        bd = pdata['domains'].get('F2', seq)
                    elif role == PrimerRole.BIP:
                        bd = pdata['domains'].get('B2', seq)
                    else:
                        bd = list(pdata['domains'].values())[0]
                        
                primers.append(PhysicalPrimer(name, seq, role=role, binding_domain=bd, nominal_concentration=conc))
                
    target_seq = "GCCACCTTAAGCCACAGTACGGAAGAAGCTGTGCAGCCTGTGAGCCCCGTCCAAGGACGTTAAAAGAAGAAGTCAGGCCCAAAAGCCACGGTTTGAGCAAACCGTGCTGCCTGTAGCTCCGTCGTGGGGACGTAAAGCCTGGGAGGCTGCAAACCGTGGAAGCTGTACGCACGGTGTAGCAGACTAGTGGTTAGAGGAGACCCCTCCCATGACACAAC"
    
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    buffer_cond = {'na_mm': 50.0, 'mg_mm': 8.0, 'dntp_mm': 1.4}
    
    prob, species, complexes, _ = enumerate_complexes(
        primers=primers,
        target_seq=target_seq,
        backend=backend,
        temp_celsius=63.0,
        buffer=buffer_cond
    )
    
    primer_to_panel = {p.name: "DEN-3" for p in primers}
    target_dict = {"DEN-3": target_seq}
    
    from labcraft.optimize.concentrations import get_dangerous_dimers, optimize_concentrations
    dangerous_complexes = get_dangerous_dimers(prob, complexes, species, primers, 63.0, backend, BST_2_0)
    print("Dangerous complexes identified:", dangerous_complexes)
    
    res_before = solve_dual(prob)
    
    # We can calculate the terms directly using the optimizer's evaluate func logic
    # But it's easier to just run the optimizer on the original prob and it will log if needed, or we can just print the terms.
    # I will modify optimize_concentrations to return the initial metrics as well, or I can just re-evaluate here.
    total_primer_conc = sum(prob.total_concentrations)
    
    dangerous_indices = [i for i, c in enumerate(complexes) if c in dangerous_complexes]
    dang_frac = sum(res_before.concentrations[dangerous_indices])
    normalized_dang = dang_frac / total_primer_conc if total_primer_conc > 0 else 0
    
    print("\nBEFORE OPTIMIZATION SCORE TERMS:")
    print(f"Dangerous term: {normalized_dang * 1e6:.2f}")
    
    print("\nBEFORE OPTIMIZATION:")
    for c_idx, c_name in enumerate(complexes):
        if c_name in dangerous_complexes:
            print(f"Dimer {c_name}: {res_before.concentrations[c_idx]*1e9:.2f} nM")
            
    for c_idx, c_name in enumerate(complexes):
        if "_on_" in c_name and c_name.endswith("_site"):
            site_name = c_name.split("_on_")[1]
            try:
                site_idx = species.index(site_name)
            except ValueError:
                continue
            site_conc = prob.total_concentrations[site_idx]
            occ = res_before.concentrations[c_idx] / site_conc if site_conc > 0 else 0
            print(f"Init {c_name}: {occ*100:.2f}%")
            
    print("\n--- RUNNING OPTIMIZER ---")
    results = optimize_concentrations(
        prob_template=prob,
        species_names=species,
        primers=primers,
        target_dict=target_dict,
        primer_to_panel=primer_to_panel,
        original_free_fractions={},
        original_target_occupations={},
        complex_names=complexes,
        temp_celsius=63.0,
        backend=backend,
        enzyme=BST_2_0
    )
    
    print("\nAFTER OPTIMIZATION:")
    for r in results:
        print(f"{r['primer_name']}: {r['original_conc']*1e6:.2f} uM -> {r['suggested_conc']*1e6:.2f} uM ({r['reason']})")
        
    for r in results:
        p_name = r['primer_name']
        idx = species.index(p_name)
        prob.total_concentrations[idx] = r['suggested_conc']
        
    res_after = solve_dual(prob)
    
    print("\n--- METRICS AFTER OPTIMIZATION ---")
    for c_idx, c_name in enumerate(complexes):
        if c_name in dangerous_complexes:
            print(f"Dimer {c_name}: {res_after.concentrations[c_idx]*1e9:.2f} nM")
            
    for c_idx, c_name in enumerate(complexes):
        if "_on_" in c_name and c_name.endswith("_site"):
            site_name = c_name.split("_on_")[1]
            try:
                site_idx = species.index(site_name)
            except ValueError:
                continue
            site_conc = prob.total_concentrations[site_idx]
            occ = res_after.concentrations[c_idx] / site_conc if site_conc > 0 else 0
            print(f"Init {c_name}: {occ*100:.2f}%")
        
    if not results:
        print("No changes suggested.")

if __name__ == '__main__':
    main()
