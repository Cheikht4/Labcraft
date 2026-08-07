"""LAMP/PCR stoichiometry and concentration profiles / Stœchiométrie.

Gère les profils de concentration asymétriques (ex: LAMP).
"""
from dataclasses import dataclass
from .domains import PrimerRole

@dataclass(frozen=True)
class ConcentrationProfile:
    """Profil de concentrations en Molaire (M)."""
    target: float
    fip_bip: float
    f3_b3: float
    lf_lb: float
    fwd_rev: float = 0.5e-6
    probe: float = 0.2e-6
    
    def get_concentration(self, role: PrimerRole) -> float:
        if role in (PrimerRole.FIP, PrimerRole.BIP):
            return self.fip_bip
        if role in (PrimerRole.F3, PrimerRole.B3):
            return self.f3_b3
        if role in (PrimerRole.LF, PrimerRole.LB):
            return self.lf_lb
        if role in (PrimerRole.FWD, PrimerRole.REV):
            return self.fwd_rev
        if role == PrimerRole.PROBE:
            return self.probe
        raise ValueError(f"Rôle {role} non géré.")

# 1 copie/µL = 1e6 copies/L.
# 1 mole = 6.022e23 copies.
# Concentration = (1e6 / 6.022e23) M = 1.6605e-18 M
COPIES_PER_UL_TO_MOLAR = 1e6 / 6.02214076e23

def target_copies_to_molar(copies_per_ul: float) -> float:
    return copies_per_ul * COPIES_PER_UL_TO_MOLAR

# Valeurs par défaut typiques pour LAMP
LAMP_DEFAULT_PROFILE = ConcentrationProfile(
    target=target_copies_to_molar(1e5), # 100 000 copies/µL par défaut
    fip_bip=1.6e-6,
    f3_b3=0.2e-6,
    lf_lb=0.4e-6
)

PCR_DEFAULT_PROFILE = ConcentrationProfile(
    target=target_copies_to_molar(1e5),
    fip_bip=0.0,
    f3_b3=0.0,
    lf_lb=0.0,
    fwd_rev=0.5e-6,
    probe=0.2e-6
)
