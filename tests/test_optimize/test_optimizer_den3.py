import pathlib
import pytest
from labcraft.optimize.optimizer import optimize_primer
from validation.validate_parida import load_parida_csv, run_analysis

def test_optimizer_den3_flp():
    path = pathlib.Path(__file__).parents[2] / "validation/reference_data/parida2005_lamp_primers.csv"
    serotypes = load_parida_csv(path)
    den3_primers = serotypes["DEN-3"]
    
    # Run optimizer
    top_variants = optimize_primer(
        primer_name="FLP_DEN-3",
        primers=den3_primers,
        temp_celsius=63.0,
        max_mutations=2,
        window_3p=6
    )
    
    assert len(top_variants) > 0
    best_variant = top_variants[0]
    
    # Should have at most 2 mutations
    assert best_variant['num_mutations'] <= 2
    
    # Should have a dG_3p > -2.0
    assert best_variant['worst_dg_3p'] > -2.0
    
    # Replace primer in panel and verify zero amplifiable dimers
    modified_primers = list(den3_primers)
    for i, p in enumerate(modified_primers):
        if p['name'] == "FLP_DEN-3":
            p_copy = dict(p)
            p_copy['sequence'] = best_variant['sequence']
            modified_primers[i] = p_copy
            break
            
    dimers_after = run_analysis("DEN-3", modified_primers, skip_flp=False)
    assert len(dimers_after) == 0
