import pytest
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.lamp.complex_enumeration import enumerate_complexes
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.thermo.salt import UnifiedSaltModel

def test_stoichiometry_target_binding_is_correct():
    primers = [
        PhysicalPrimer.from_simple("F3_DEN3", "ATCG", PrimerRole.F3),
        PhysicalPrimer.from_simple("B3_DEN3", "GCTA", PrimerRole.B3),
        PhysicalPrimer.from_simple("FIP_DEN3", "CGAT", PrimerRole.FIP),
    ]
    target_seq = "ATCGNNNNGCTANNNNCGAT"
    backend = ViennaSaltShiftBackend(UnifiedSaltModel())
    
    prob, strands, complex_names, unfold = enumerate_complexes(primers, target_seq, backend, temp_celsius=65.0, mon_molar=0.05)
    
    strand_names = strands
    
    failures = []
    for idx, c_name in enumerate(complex_names):
        if "_on_" in c_name and c_name.endswith("_site"):
            primer_name = c_name.split("_on_")[0]
            # Find which primer index has stoich == 1
            # We expect the primer strand to be in the first len(primers) slots
            primer_indices = [i for i, val in enumerate(prob.stoichiometry[idx][:len(primers)]) if val == 1]
            if not primer_indices:
                failures.append(f"{c_name}: No primer consumed")
                continue
            
            consumed_primer_idx = primer_indices[0]
            consumed_primer_name = strand_names[consumed_primer_idx]
            
            if consumed_primer_name != primer_name:
                failures.append(f"{c_name}: claims to be {primer_name} but consumes {consumed_primer_name}")
                
    assert len(failures) == 0, f"Stoichiometry mismatch: {failures}"
