from validation.validate_parida import load_parida_csv
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.lamp.complex_enumeration import enumerate_complexes
from labcraft.thermo.salt import UnifiedSaltModel
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.diagnostics.enzyme import BST_2_0
from labcraft.optimize.concentrations import optimize_concentrations
from labcraft.solver.dual import solve_dual
from labcraft.optimize.concentrations import get_dangerous_dimers

def test_abundant_dimer():
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    fip = PhysicalPrimer("FIP_1", sequence="ATCGTACGATCGATCGGGGGG", role=PrimerRole.FIP, binding_domain="ATCGTACGATCGATC", nominal_concentration=1.6e-6)
    bip = PhysicalPrimer("BIP_1", sequence="CCCCCCGATCGATCGTACGAT", role=PrimerRole.BIP, binding_domain="GATCGTACGAT", nominal_concentration=1.6e-6)
    
    target_seq = "ATGCGTACGTGCAACTGATCGATCGTACGATCG"
    
    prob, species, complexes, _ = enumerate_complexes(
        primers=[fip, bip],
        target_seq=target_seq,
        backend=backend,
        temp_celsius=65.0
    )
    
    dangerous = get_dangerous_dimers(prob, complexes, species, [fip, bip], 65.0, backend, BST_2_0)
    print("Dangerous complexes:", dangerous)
    
    res_before = solve_dual(prob)
    total_primer_conc = sum(prob.total_concentrations)
    dangerous_indices = [i for i, c in enumerate(complexes) if c in dangerous]
    dang_frac = sum(res_before.concentrations[dangerous_indices])
    normalized_dang = dang_frac / total_primer_conc if total_primer_conc > 0 else 0
    print(f"Normalized dang before: {normalized_dang * 1e6}")
    
    results = optimize_concentrations(
        prob_template=prob,
        species_names=species,
        primers=[fip, bip],
        target_dict={'PANEL1': target_seq},
        primer_to_panel={'FIP_1': 'PANEL1', 'BIP_1': 'PANEL1'},
        original_free_fractions={},
        original_target_occupations={},
        complex_names=complexes,
        temp_celsius=65.0,
        backend=backend,
        enzyme=BST_2_0
    )
    
    print("Optimization results:", results)
    
    if results:
        for r in results:
            p_name = r['primer_name']
            idx = species.index(p_name)
            prob.total_concentrations[idx] = r['suggested_conc']
            
        res_after = solve_dual(prob)
        dang_frac_after = sum(res_after.concentrations[dangerous_indices])
        normalized_dang_after = dang_frac_after / total_primer_conc if total_primer_conc > 0 else 0
        print(f"Normalized dang after: {normalized_dang_after * 1e6}")

test_abundant_dimer()
