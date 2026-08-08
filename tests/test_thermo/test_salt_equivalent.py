
def test_sodium_equivalent_for_folding_anchor():
    from labcraft.thermo.salt import sodium_equivalent_for_folding
    # Ancrage : sodium_equivalent_for_folding(0.05, 0.008, 0.0) vaut ~0,200 M (tolérance ± 0,01 M). C'est LA valeur de référence calibrée.
    na_eq = sodium_equivalent_for_folding(0.05, 0.008, 0.0)
    assert 0.19 < na_eq < 0.21, f"Expected around 0.200 M, got {na_eq}"

def test_sodium_equivalent_for_folding_zero_mg():
    from labcraft.thermo.salt import sodium_equivalent_for_folding
    # Mg nul : sodium_equivalent_for_folding(0.05, 0.0, 0.0) == 0.05 exactement.
    na_eq = sodium_equivalent_for_folding(0.05, 0.0, 0.0)
    assert abs(na_eq - 0.05) < 1e-6, f"Expected 0.05 M, got {na_eq}"

def test_sodium_equivalent_for_folding_dntp():
    from labcraft.thermo.salt import sodium_equivalent_for_folding
    # dNTP : le Mg libre baisse quand les dNTP montent, donc na_eq_dntp < na_eq_no_dntp
    na_eq_no_dntp = sodium_equivalent_for_folding(0.05, 0.008, 0.0)
    na_eq_dntp = sodium_equivalent_for_folding(0.05, 0.008, 0.006)
    assert na_eq_dntp < na_eq_no_dntp, "dNTPs should lower free Mg2+ and thus lower Na_eq"

def test_sodium_equivalent_for_folding_monotony():
    from labcraft.thermo.salt import sodium_equivalent_for_folding
    # Monotonie : plus de Mg -> Na_eq plus élevé.
    na_eq_low = sodium_equivalent_for_folding(0.05, 0.002, 0.0)
    na_eq_high = sodium_equivalent_for_folding(0.05, 0.008, 0.0)
    assert na_eq_high > na_eq_low, "Higher Mg2+ should result in higher Na_eq"

def test_sodium_equivalent_for_folding_cap():
    from labcraft.thermo.salt import sodium_equivalent_for_folding
    # Plafond : un Mg énorme ne dépasse pas cap_molar.
    na_eq = sodium_equivalent_for_folding(0.05, 1.0, 0.0, cap_molar=1.0)
    assert na_eq <= 1.0, "Na_eq should not exceed cap_molar"
