import pytest
import numpy as np

from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.lamp.stoichiometry import ConcentrationProfile, target_copies_to_molar
from labcraft.lamp.complex_enumeration import enumerate_complexes
from labcraft.thermo.backends.native import NativeBackend
from labcraft.solver.dual import solve_dual

from labcraft.thermo.backends.base import DuplexEnergyBackend, DuplexResult

class MockBackend(DuplexEnergyBackend):
    def calc_heterodimer(self, seq1, seq2, **kwargs):
        return DuplexResult(0, 0, -1.0, 0, "", 65.0)
    def calc_homodimer(self, seq, **kwargs):
        return DuplexResult(0, 0, -1.0, 0, "", 65.0)
    def calc_hairpin(self, seq, **kwargs):
        return DuplexResult(0, 0, -1.0, 0, "", 65.0)
    def calc_duplex(self, seq1, seq2, **kwargs):
        return DuplexResult(0, 0, -10.0, 0, "", 65.0)

def test_complex_enumeration_counts():
    """Vérifie le dénombrement exact pour un panel LAMP 2-plex.
    
    12 oligos -> 12 monomers + 66 heterodimers + 12 homodimers = 90 espèces oligo.
    Plus les espèces "site libre" et les duplex "amorce-cible".
    """
    backend = MockBackend()
    
    # On crée 12 oligos distincts (2 panels de 6)
    primers = []
    for p_id in range(1, 3): # Panel 1 et 2
        for role_name in ["F3", "B3", "FIP", "BIP", "LF", "LB"]:
            role = PrimerRole[role_name]
            # Des séquences purement bidons, assez distinctes
            # Utilisons l'index pour générer un domaine de liaison unique
            idx = (p_id - 1) * 6 + ["F3", "B3", "FIP", "BIP", "LF", "LB"].index(role_name)
            base_motif = ["A", "C", "G", "T"][idx % 4] * 5 + ["T", "G", "C", "A"][idx % 4] * 5 + f"{idx:05d}"
            seq = base_motif + f"_tail_{p_id}_{role_name}"
            # On force le domaine de liaison à un petit bout
            binding = base_motif
            primers.append(PhysicalPrimer(
                name=f"{role_name}_{p_id}",
                sequence=seq,
                role=role,
                binding_domain=binding,
                tail_domain=seq[15:],
                linker=None
            ))
            
    # Séquence cible (très basique, ne matche rien de manière pertinente)
    target = "ATGC" * 50
    
    # Énumération
    profile = ConcentrationProfile(target=1e-12, fip_bip=1.6e-6, f3_b3=0.2e-6, lf_lb=0.4e-6)
    prob, strand_names, complex_names, _ = enumerate_complexes(
        primers, target, backend, profile=profile, temp_celsius=65.0
    )
    
    # Validation du nombre de brins (composants de base)
    # 12 amorces physiques + 0 site (car la cible bidon ne matche pas les domaines)
    # Attends, si elle ne matche pas, on n'a que 12 brins.
    # On va tricher et rajouter la cible qui matche le premier F3
    target_match = primers[0].binding_domain
    target = target_match + "ATGC" * 20
    
    prob, strand_names, complex_names, _ = enumerate_complexes(
        primers, target, backend, profile=profile, temp_celsius=65.0
    )
    
    # 12 amorces + 1 site cible = 13 brins
    assert prob.n_strands == 13
    
    # Nombre de complexes:
    # 13 formes libres (12 monomères amorces + 1 monomère site cible)
    # 66 hétérodimères
    # 12 homodimères
    # 1 complexe amorce-cible
    # Total = 13 + 66 + 12 + 1 = 92 complexes
    
    # Compte manuel des types dans complex_names
    frees = [n for n in complex_names if n.endswith("_free")]
    homos = [n for n in complex_names if n.endswith("_homo")]
    ons = [n for n in complex_names if "_on_" in n]
    heteros = [n for n in complex_names if n not in frees and n not in homos and n not in ons]
    
    assert len(frees) == 13, f"Attendu 13, eu {len(frees)}"
    assert len(homos) == 12, f"Attendu 12, eu {len(homos)}"
    assert len(heteros) == 66, f"Attendu 66, eu {len(heteros)}"
    assert len(ons) == 1, f"Attendu 1, eu {len(ons)}"


def test_target_trace_consistency():
    """Une cible en trace ne doit pas déplacer l'équilibre des amorces."""
    backend = NativeBackend()
    
    # Amorce propre
    p_clean = PhysicalPrimer("P_clean", "ATGCATGCATGCATGC", PrimerRole.F3, "ATGCATGCATGCATGC")
    target = "GCATGCATGCATGCAT" # Son reverse complement
    
    profile_no_target = ConcentrationProfile(target=1e-25, fip_bip=0, f3_b3=1e-6, lf_lb=0)
    profile_trace_target = ConcentrationProfile(target=target_copies_to_molar(1e5), fip_bip=0, f3_b3=1e-6, lf_lb=0) # 1e5 copies/µL
    
    prob_no, _, _, _ = enumerate_complexes([p_clean], target, backend, profile_no_target)
    prob_with, _, _, _ = enumerate_complexes([p_clean], target, backend, profile_trace_target)
    
    res_no = solve_dual(prob_no)
    res_with = solve_dual(prob_with)
    
    free_p_no = res_no.free_concentrations[0]
    free_p_with = res_with.free_concentrations[0]
    
    # L'écart relatif doit être infinitésimal (la cible est autour de 1e-13 M,
    # l'amorce à 1e-6 M).
    assert abs(free_p_with - free_p_no) / free_p_no < 1e-5


class DimerMockBackend(DuplexEnergyBackend):
    def calc_heterodimer(self, seq1, seq2, **kwargs):
        return DuplexResult(0, 0, 1.0, 0, "", 65.0)
    def calc_homodimer(self, seq, **kwargs):
        if "CGCGCGCG" in seq:
            return DuplexResult(0, 0, -15.0, 0, "", 65.0)
        return DuplexResult(0, 0, 1.0, 0, "", 65.0)
    def calc_hairpin(self, seq, **kwargs):
        return DuplexResult(0, 0, 1.0, 0, "", 65.0)
    def calc_duplex(self, seq1, seq2, **kwargs):
        return DuplexResult(0, 0, -10.0, 0, "", 65.0)

def test_strong_dimer_depletion():
    """Une amorce avec un fort homodimère est appauvrie et se lie moins à la cible."""
    backend = DimerMockBackend()
    
    # Amorce propre
    p_clean = PhysicalPrimer("P_clean", "ATGCATGCATGCATGC", PrimerRole.F3, "ATGCATGCATGCATGC")
    
    # Amorce auto-complémentaire (CGCGCGCGCGCGCGCG)
    p_dimer = PhysicalPrimer("P_dimer", "CGCGCGCGCGCGCGCG", PrimerRole.F3, "CGCGCGCGCGCGCGCG")
    
    target_clean = "GCATGCATGCATGCAT"
    target_dimer = "CGCGCGCGCGCGCGCG"
    
    profile = ConcentrationProfile(target=target_copies_to_molar(1e5), fip_bip=0, f3_b3=1e-6, lf_lb=0)
    
    prob_clean, _, _, _ = enumerate_complexes([p_clean], target_clean, backend, profile)
    prob_dimer, _, _, _ = enumerate_complexes([p_dimer], target_dimer, backend, profile)
    
    res_clean = solve_dual(prob_clean)
    res_dimer = solve_dual(prob_dimer)
    
    free_clean = res_clean.free_concentrations[0]
    free_dimer = res_dimer.free_concentrations[0]
    
    # L'amorce sale forme massivement un homodimère, donc sa concentration libre est très faible
    assert free_dimer < free_clean * 0.1, "Le fort dimère devrait appauvrir l'amorce libre"
    
    # Vérifions l'occupation de la cible
    # L'espèce cible est l'indice 1, le complexe est l'indice 3 (0=A_free, 1=T_free, 2=A_homo, 3=A_on_T)
    # Complexe concentration = exp(-dg/RT) * A_free * T_free
    # ... ou simplement on regarde T_free
    t_free_clean = res_clean.free_concentrations[1]
    t_free_dimer = res_dimer.free_concentrations[1]
    
    # T_free sera plus HAUT pour l'amorce sale (car elle ne se lie pas à la cible)
    assert t_free_dimer > t_free_clean, "La cible de l'amorce sale est moins occupée"
