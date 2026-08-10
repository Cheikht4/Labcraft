"""Configuration models for LabCraft CLI.

Modèles Pydantic pour la configuration YAML/JSON.
"""
from typing import Dict, List, Optional, Any, Union
from pydantic import BaseModel, Field, model_validator

class BufferConfig(BaseModel):
    na_mM: float = 50.0
    k_mM: float = 0.0
    tris_mM: float = 0.0
    mg_mM: float = 0.0
    dntp_mM: float = 0.0
    
    model_config = {"extra": "forbid"}

class ExperimentConfig(BaseModel):
    name: Optional[str] = None
    chemistry: str = "LAMP"
    temperature_C: float = 65.0
    enzyme: str = "bst2.0"
    buffer: Optional[BufferConfig] = None
    
    model_config = {"extra": "forbid"}

class TargetConfig(BaseModel):
    id: str
    sequence_file: str
    copies_per_uL: float = 1000.0
    
    model_config = {"extra": "forbid"}

class PrimerDomains(BaseModel):
    F2: Optional[str] = None
    B2: Optional[str] = None
    F1c: Optional[str] = None
    B1c: Optional[str] = None
    linker: Optional[str] = ""

class PrimerConfig(BaseModel):
    seq: str
    conc_uM: Optional[float] = None
    domains: Union[PrimerDomains, str, None] = None
    blocked_3prime: bool = False
    mod_3prime: Optional[str] = None
    label_5prime: Optional[str] = None
    
    model_config = {"extra": "forbid"}

class PrimerSetConfig(BaseModel):
    target: str
    primers: Dict[str, PrimerConfig]
    
    model_config = {"extra": "forbid"}

class PanelConfig(BaseModel):
    experiment: ExperimentConfig = Field(default_factory=ExperimentConfig)
    targets: List[TargetConfig] = Field(default_factory=list)
    primer_sets: List[PrimerSetConfig] = Field(default_factory=list)
    
    model_config = {"extra": "forbid"}
    
    @model_validator(mode="after")
    def validate_targets_referenced(self) -> 'PanelConfig':
        target_ids = {t.id for t in self.targets}
        for pset in self.primer_sets:
            if pset.target not in target_ids:
                raise ValueError(
                    f"Target '{pset.target}' referenced by primer set not found in targets list. "
                    f"Available targets: {list(target_ids)}"
                )
        return self
        
    @model_validator(mode="after")
    def check_mg_dntp(self) -> 'PanelConfig':
        if self.experiment.buffer:
            mg = self.experiment.buffer.mg_mM
            dntp = self.experiment.buffer.dntp_mM
            if mg <= dntp:
                import warnings
                warnings.warn("mg_mM is less than or equal to dntp_mM.")
        return self

from typing import Tuple
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.lamp.stoichiometry import ConcentrationProfile, target_copies_to_molar
from labcraft.thermo.backends.vienna import ViennaRNABackend
from labcraft.diagnostics.enzyme import get_enzyme, PolymeraseProfile
from labcraft.thermo.backends.base import DuplexEnergyBackend
from labcraft.thermo.salt import UnifiedSaltModel, SaltCorrectedBackend, sodium_equivalent_for_folding
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.buffer.monovalent import get_total_monovalent

def build_engine_from_config(
    config: PanelConfig, 
    targets: Dict[str, str]
) -> Tuple[List[PhysicalPrimer], Dict[str, str], DuplexEnergyBackend, Dict[str, Any], float, PolymeraseProfile, float, Dict[str, ConcentrationProfile]]:
    """Construit les structures internes à partir d'une configuration.
    
    Returns:
        primers: Liste d'amorces physiques.
        primer_to_panel: Dictionnaire associant chaque nom d'amorce à l'ID de sa cible.
        backend: Backend thermodynamique (avec ou sans sel).
        backend_kwargs: Dictionnaire des paramètres de sel (na_mm, mg_mm, etc.) si applicable.
        mon_molar: Molarité des monovalents si applicable.
        enzyme: Profil de la polymérase.
        temp_celsius: Température de l'expérience.
        profiles: Dictionnaire des profils de concentration par cible.
    """
    temp_celsius = config.experiment.temperature_C
    
    # 1. Tampon Salin
    buffer_conf = config.experiment.buffer
    backend_kwargs = {}
    if buffer_conf:
        backend = ViennaSaltShiftBackend(UnifiedSaltModel())
        backend_kwargs = {
            'na_mm': buffer_conf.na_mM,
            'k_mm': buffer_conf.k_mM,
            'tris_mm': buffer_conf.tris_mM,
            'mg_mm': buffer_conf.mg_mM,
            'dntp_mm': buffer_conf.dntp_mM
        }
        mon_true = get_total_monovalent(backend_kwargs['na_mm'], backend_kwargs['k_mm'], backend_kwargs['tris_mm']) / 1000.0
        mg_molar = backend_kwargs['mg_mm'] / 1000.0
        dntp_molar = (backend_kwargs['dntp_mm'] or 0.0) / 1000.0
        mon_molar = sodium_equivalent_for_folding(mon_true, mg_molar, dntp_molar)
    else:
        backend = ViennaRNABackend()
        mon_molar = None

    # 2. Enzyme
    enzyme = get_enzyme(config.experiment.enzyme)
    
    # 3. Amorces et profils de concentration
    primers = []
    primer_to_panel = {}
    profiles = {}
    
    # Création du dictionnaire de copies par cible
    target_copies = {t.id: t.copies_per_uL for t in config.targets}
    
    for pset in config.primer_sets:
        t_id = pset.target
        p_dict = pset.primers
        
        copies_uL = target_copies.get(t_id, 1000.0)
        target_molar = target_copies_to_molar(copies_uL)
        
        # On garde une trace des concentrations pour le profile de cette cible
        conc_map = {}
        
        for role_name, p_data in p_dict.items():
            name = f"{role_name}_{t_id.replace('Synth', '')}"
            seq = p_data.seq
            primer_to_panel[name] = t_id
            
            try:
                role_enum = PrimerRole[role_name.upper()]
            except KeyError:
                role_enum = PrimerRole.F3
                
            
            
            # Détermination de la concentration nominale totale de l'amorce
            total_conc = p_data.conc_uM * 1e-6 if p_data.conc_uM is not None else (
                1.6e-6 if role_name.upper() in ("FIP", "BIP") else
                0.2e-6 if role_name.upper() in ("F3", "B3", "PROBE") else
                0.8e-6
            )
            
            # Gestion du blocage 3'
            is_blocked = p_data.blocked_3prime
            if p_data.mod_3prime and p_data.mod_3prime.upper() in ['BHQ', 'TAMRA', 'QUENCHER', 'DSPACER', 'C3', 'INVDT', 'PHOS', 'DDC']:
                is_blocked = True
                
            label = p_data.label_5prime
            
            from labcraft.lamp.domains import expand_degenerate
            from labcraft.thermo.lna import parse_lna_sequence
            variants = expand_degenerate(seq)
            variant_conc = total_conc / len(variants)
            
            # Enregistrer le nom parent dans primer_to_panel (pour les cas sans dégénérescence)
            # Register parent name in primer_to_panel (for non-degenerate cases)
            primer_to_panel[name] = t_id
            
            # Parse la séquence parente pour localiser les sous-domaines par position
            # Parse the parent sequence to locate sub-domains by position
            parent_bare, _ = parse_lna_sequence(seq)
            
            for v_idx, variant_seq in enumerate(variants):
                bare_seq, lna_pos = parse_lna_sequence(variant_seq)
                v_name = f"{name}#{v_idx+1}" if len(variants) > 1 else name
                
                # Enregistrer chaque variant dans primer_to_panel
                # Register each variant in primer_to_panel
                if len(variants) > 1:
                    primer_to_panel[v_name] = t_id
                
                if role_name.upper() in ("FIP", "BIP") and isinstance(p_data.domains, PrimerDomains):
                    d = p_data.domains
                    # Développer les sous-domaines par position dans le variant complet
                    # Expand sub-domains by position in the full variant sequence
                    d_f2_raw = d.F2 or d.B2 or ""
                    d_f1c_raw = d.F1c or d.B1c or ""
                    d_linker_raw = d.linker or ""
                    
                    d_f2_bare, _ = parse_lna_sequence(d_f2_raw)
                    d_f1c_bare, _ = parse_lna_sequence(d_f1c_raw)
                    d_linker_bare, _ = parse_lna_sequence(d_linker_raw)
                    
                    if len(variants) > 1 and d_f2_bare and d_f1c_bare:
                        # Localiser les sous-domaines dans la séquence parente nue
                        # Locate sub-domains in the bare parent sequence
                        # Structure FIP: 5'-F1c-linker-F2-3'
                        # Structure BIP: 5'-B1c-linker-B2-3'
                        f1c_start = 0
                        f1c_end = len(d_f1c_bare)
                        linker_end = f1c_end + len(d_linker_bare)
                        f2_start = linker_end
                        f2_end = f2_start + len(d_f2_bare)
                        
                        # Extraire les sous-domaines depuis le variant par position
                        # Extract sub-domains from variant by position
                        variant_f1c = bare_seq[f1c_start:f1c_end]
                        variant_linker = bare_seq[f1c_end:linker_end]
                        variant_f2 = bare_seq[f2_start:f2_end]
                        
                        primers.append(PhysicalPrimer(
                            v_name, bare_seq, role_enum, variant_f2, variant_f1c, variant_linker,
                            nominal_concentration=variant_conc, parent_name=name,
                            blocked_3prime=is_blocked, label_5prime=label,
                            lna_positions=tuple(lna_pos)
                        ))
                    else:
                        primers.append(PhysicalPrimer(
                            v_name, bare_seq, role_enum, d_f2_bare, d_f1c_bare, d_linker_bare,
                            nominal_concentration=variant_conc, parent_name=name,
                            blocked_3prime=is_blocked, label_5prime=label,
                            lna_positions=tuple(lna_pos)
                        ))
                elif role_name.upper() in ("FIP", "BIP"):
                    target_seq = targets.get(t_id, "")
                    primers.append(PhysicalPrimer.from_alignment(
                        v_name, variant_seq, role_enum, target_seq, 
                        nominal_concentration=variant_conc, parent_name=name,
                        blocked_3prime=is_blocked, label_5prime=label
                    ))
                else:
                    primers.append(PhysicalPrimer(
                        v_name, bare_seq, role_enum, bare_seq, 
                        nominal_concentration=variant_conc, parent_name=name,
                        blocked_3prime=is_blocked, label_5prime=label,
                        lna_positions=tuple(lna_pos)
                    ))
                
        # Création du profil pour la cible, avec fallback par défaut si absent
        profile = ConcentrationProfile(
            target=target_molar,
            fip_bip=conc_map.get(PrimerRole.FIP, conc_map.get(PrimerRole.BIP, 1.6e-6)),
            f3_b3=conc_map.get(PrimerRole.F3, conc_map.get(PrimerRole.B3, 0.2e-6)),
            lf_lb=conc_map.get(PrimerRole.LF, conc_map.get(PrimerRole.LB, 0.8e-6))
        )
        profiles[t_id] = profile

    # Vérification présence LNA pour forcer le backend SaltShift (même sans tampon)
    has_lna = any(len(p.lna_positions) > 0 for p in primers)
    if has_lna and not buffer_conf:
        backend = ViennaSaltShiftBackend(UnifiedSaltModel())
        backend_kwargs = {
            'na_mm': 1000.0, # Neutral salt pour annuler ddG_sel
            'k_mm': 0.0,
            'tris_mm': 0.0,
            'mg_mm': 0.0,
            'dntp_mm': 0.0
        }
        # mon_molar reste None pour ne pas perturber calc_unfolding_penalty (qui utilise 1M par défaut)
        
    return primers, primer_to_panel, backend, backend_kwargs, mon_molar, enzyme, temp_celsius, profiles
