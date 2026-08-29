"""Parsers for LabCraft CLI inputs (FASTA, TXT)."""
import re
from typing import List, Tuple, Dict, Optional
from collections import defaultdict

from labcraft.cli.config import PrimerConfig, PrimerSetConfig, PrimerDomains
from labcraft.lamp.domains import PrimerRole
from labcraft.lamp.domains import _find_iupac_substring

class ParseError(Exception):
    pass

def read_multi_fasta(filepath: str) -> List[Tuple[str, str]]:
    """
    Reads a FASTA or mixed TXT file.
    Returns a list of tuples: (header, sequence).
    Raises ParseError on duplicate headers or bad format.
    
    Supports:
    - Standard FASTA: >Header\nSequence
    - One-line format: >Header   Sequence
    """
    records = []
    seen_headers = set()
    
    with open(filepath, 'r') as f:
        lines = f.read().splitlines()
        
    current_header = None
    current_seq_parts = []
    
    def save_current():
        if current_header is not None:
            if current_header in seen_headers:
                raise ParseError(f"Duplicate header found in {filepath}: '{current_header}'")
            seq = "".join(current_seq_parts).upper().replace("U", "T")
            records.append((current_header, seq))
            seen_headers.add(current_header)
            
    for line in lines:
        line = line.strip()
        if not line:
            continue
            
        if line.startswith(">"):
            save_current()
            # Check for one-line format: >Header Sequence
            parts = line[1:].split(None, 1)
            current_header = parts[0]
            current_seq_parts = []
            if len(parts) > 1:
                # Sequence is on the same line
                current_seq_parts.append(parts[1])
        else:
            if current_header is None:
                # If no header yet, maybe it's just raw sequence or a bad file
                raise ParseError(f"File {filepath} doesn't start with a '>' header.")
            current_seq_parts.append(line)
            
    save_current()
    
    return records


def parse_primer_name(name: str) -> Tuple[str, str, Optional[str]]:
    """
    Parses a primer name into (PanelName, Role, Version).
    Valid roles: F3, B3, FIP, BIP, LF, LB, F1, F2, B1, B2, FWD, REV, PROBE.
    Aliases: FLOOP, FLP -> LF; BLOOP, BLP -> LB.
    
    The role is expected to be one of the underscore-separated components.
    Usually it's the last or second to last (if version is present).
    """
    aliases = {
        "FLOOP": "LF", "FLP": "LF", "LOOPF": "LF",
        "BLOOP": "LB", "BLP": "LB", "LOOPB": "LB"
    }
    valid_roles = {"F3", "B3", "FIP", "BIP", "LF", "LB", "F1", "F2", "B1", "B2", "FWD", "REV", "PROBE"}
    
    parts = name.split("_")
    
    # Try to find the role from right to left
    role = None
    role_idx = -1
    for i in range(len(parts)-1, -1, -1):
        p_upper = parts[i].upper()
        p_resolved = aliases.get(p_upper, p_upper)
        if p_resolved in valid_roles:
            role = p_resolved
            role_idx = i
            break
            
    if role is None:
        raise ParseError(f"Primer '{name}' contains no valid role. Accepted roles: {', '.join(sorted(valid_roles))} (and aliases like FLOOP/BLOOP). Note: STEMF is not a valid alias for F1c.")
        
    version = None
    if role_idx < len(parts) - 1:
        # Anything after the role is considered the version (e.g. DENV3_F3_1 -> 1)
        version = "_".join(parts[role_idx+1:])
        
    panel_name = "_".join(parts[:role_idx])
    if not panel_name:
        # Fallback if the name was just the role, e.g. "F3"
        panel_name = "DefaultPanel"
        
    return panel_name, role, version


def reverse_complement(seq: str) -> str:
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
            'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 
            'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}
    return "".join(comp.get(c, c) for c in reversed(seq))


def reconstruct_fip_bip(
    panel_name: str,
    primers: Dict[str, str],
    targets: List[Tuple[str, str]],
    linker: str = ""
) -> Dict[str, PrimerDomains]:
    """
    Reconstructs FIP from F1/F2 and BIP from B1/B2 in-place within `primers`.
    Verifies orientation using `_find_iupac_substring` against the target sequence.
    Returns a dictionary of PrimerDomains for reconstructed composite primers.
    """
    # Find matching target
    target_seq = None
    if len(targets) == 1:
        target_seq = targets[0][1]
    else:
        for t_name, t_seq in targets:
            if t_name == panel_name:
                target_seq = t_seq
                break

    def check_orientation(p1: str, p2: str, role_prefix: str) -> Tuple[str, str]:
        # p1 is F1 or B1. p2 is F2 or B2.
        # Returns (composite_seq, tail_domain_seq)
        if not target_seq:
            p1c = reverse_complement(p1)
            return p1c + linker + p2, p1c
            
        p1_matches = _find_iupac_substring(p1, target_seq) >= 0
        p1rc_matches = _find_iupac_substring(reverse_complement(p1), target_seq) >= 0
        
        # Si la recherche exacte échoue, tenter un criblage tolérant (max 2 erreurs)
        # If exact search fails, try tolerant screening (up to 2 errors)
        if not p1_matches and not p1rc_matches:
            from labcraft.target.seeding import primer_matches_sequence
            match_res = primer_matches_sequence(target_seq, p1, max_errors=2, strict_3prime_len=0)
            if match_res is not None:
                _, _, strand = match_res
                if strand == '+':
                    p1_matches = True
                else:
                    p1rc_matches = True
        
        if p1_matches and not p1rc_matches:
            # p1 is in the same sense as the target. We need its RC.
            p1c = reverse_complement(p1)
            return p1c + linker + p2, p1c
        elif p1rc_matches and not p1_matches:
            # p1 was provided already as the complement.
            return p1 + linker + p2, p1
        elif p1_matches and p1rc_matches:
            raise ParseError(
                f"Panel '{panel_name}': Unable to determine orientation for {role_prefix}1. "
                f"Both orientations match the target '{panel_name}'."
            )
        else:
            import warnings
            warnings.warn(f"L'orientation de {role_prefix}1 pour le panel '{panel_name}' n'a pas pu être vérifiée sur la cible (aucune correspondance). Convention standard (RC) appliquée.")
            p1c = reverse_complement(p1)
            return p1c + linker + p2, p1c

    domains_created: Dict[str, PrimerDomains] = {}

    if 'F1' in primers and 'F2' in primers:
        f1_seq = primers['F1']
        f2_seq = primers['F2']
        if 'FIP' not in primers:
            fip_seq, f1c_seq = check_orientation(f1_seq, f2_seq, 'F')
            primers['FIP'] = fip_seq
            domains_created['FIP'] = PrimerDomains(F2=f2_seq, F1c=f1c_seq, linker=linker)
        del primers['F1']
        del primers['F2']
        
    if 'B1' in primers and 'B2' in primers:
        b1_seq = primers['B1']
        b2_seq = primers['B2']
        if 'BIP' not in primers:
            bip_seq, b1c_seq = check_orientation(b1_seq, b2_seq, 'B')
            primers['BIP'] = bip_seq
            domains_created['BIP'] = PrimerDomains(B2=b2_seq, B1c=b1c_seq, linker=linker)
        del primers['B1']
        del primers['B2']

    return domains_created


def parse_primer_file(
    filepath: str,
    targets: List[Tuple[str, str]],
    linker: str = "",
    allow_unmatched_targets: bool = False
) -> List[PrimerSetConfig]:
    """
    Parses a primer file and constructs PrimerSetConfig objects.
    Groups primers by PanelName and Version.
    """
    records = read_multi_fasta(filepath)
    
    # Nested dict: PanelName -> Version -> {Role: Sequence}
    panels = defaultdict(lambda: defaultdict(dict))
    
    for header, seq in records:
        # Allow spaces in header, but typically the first word is the name
        name = header.split()[0]
        panel_name, role, version = parse_primer_name(name)
        v_key = version or "default"
        
        if role in panels[panel_name][v_key]:
            raise ParseError(f"Duplicate role '{role}' for panel '{panel_name}' version '{v_key}' in {filepath}.")
            
        panels[panel_name][v_key][role] = seq
        
    config_sets = []
    
    for panel_name, versions in panels.items():
        if len(versions) > 1:
            import warnings
            warnings.warn(f"Multiple versions detected for panel '{panel_name}'. They will be treated as separate alternative panels.")
            
        # Match target
        target_name = None
        if len(targets) == 0:
            target_name = panel_name
        elif len(targets) == 1:
            target_name = targets[0][0]
        else:
            t_names = [t[0] for t in targets]
            if panel_name in t_names:
                target_name = panel_name
            else:
                for t_name in t_names:
                    if panel_name.lower() in t_name.lower() or t_name.lower() in panel_name.lower():
                        target_name = t_name
                        break
                if target_name is None:
                    if allow_unmatched_targets:
                        target_name = targets[0][0] if targets else panel_name
                    else:
                        raise ParseError(f"Panel '{panel_name}' does not match any target name. Available targets: {', '.join(t_names)}. If there are multiple targets, the panel name must match the target name exactly or by inclusion.")
                
        for v_key, roles_dict in versions.items():
            domains_map = reconstruct_fip_bip(panel_name, roles_dict, targets, linker)
            
            primer_configs = {}
            for role, seq in roles_dict.items():
                dom = domains_map.get(role)
                primer_configs[role] = PrimerConfig(seq=seq, domains=dom)
                
            p_name = panel_name if v_key == "default" else f"{panel_name}_{v_key}"
            config_sets.append(PrimerSetConfig(
                target=target_name,
                panel_name=p_name,
                primers=primer_configs
            ))
            
    return config_sets


from labcraft.cli.config import PanelConfig, ExperimentConfig, BufferConfig, TargetConfig
import yaml

def build_config_from_cli(
    config_path: Optional[str] = None,
    primers_path: Optional[str] = None,
    targets_path: Optional[str] = None,
    temperature: Optional[float] = None,
    enzyme: Optional[str] = None,
    na: Optional[float] = None,
    k: Optional[float] = None,
    tris: Optional[float] = None,
    mg: Optional[float] = None,
    dntp: Optional[float] = None,
    conc_fip_bip: Optional[float] = None,
    conc_f3_b3: Optional[float] = None,
    conc_loop: Optional[float] = None,
    copies: Optional[float] = None,
    linker: str = "",
    allow_unmatched_targets: bool = False
) -> PanelConfig:
    # 1. Load base config from YAML if provided
    base_data = {}
    if config_path:
        with open(config_path, "r") as f:
            base_data = yaml.safe_load(f) or {}

    config = PanelConfig.model_validate(base_data) if base_data else PanelConfig()

    # 2. Override Experiment settings
    if temperature is not None: config.experiment.temperature_C = temperature
    if enzyme is not None: config.experiment.enzyme = enzyme
    
    if any(x is not None for x in (na, k, tris, mg, dntp)):
        if not config.experiment.buffer:
            config.experiment.buffer = BufferConfig()
        if na is not None: config.experiment.buffer.na_mM = na
        if k is not None: config.experiment.buffer.k_mM = k
        if tris is not None: config.experiment.buffer.tris_mM = tris
        if mg is not None: config.experiment.buffer.mg_mM = mg
        if dntp is not None: config.experiment.buffer.dntp_mM = dntp

    # 3. Handle Targets
    target_list = []
    if targets_path:
        records = read_multi_fasta(targets_path)
        for h, seq in records:
            t_id = h.split()[0]
            target_list.append((t_id, seq))
            
            c = copies if copies is not None else 1000.0
            config.targets.append(TargetConfig(id=t_id, sequence_file=targets_path, copies_per_uL=c))
            
    # 4. Handle Primers
    if primers_path:
        new_sets = parse_primer_file(
            primers_path,
            target_list,
            linker=linker,
            allow_unmatched_targets=allow_unmatched_targets
        )
        
        # Override concentrations
        for pset in new_sets:
            for role, p_data in pset.primers.items():
                r_upper = role.upper()
                if r_upper in ("FIP", "BIP") and conc_fip_bip is not None:
                    p_data.conc_uM = conc_fip_bip
                elif r_upper in ("F3", "B3") and conc_f3_b3 is not None:
                    p_data.conc_uM = conc_f3_b3
                elif r_upper in ("LF", "LB") and conc_loop is not None:
                    p_data.conc_uM = conc_loop
                    
        config.primer_sets.extend(new_sets)

    # Re-valider le modèle pour déclencher les validateurs (ex: check_mg_dntp)
    # Re-validate model to trigger post-validation checks (e.g. check_mg_dntp)
    return PanelConfig.model_validate(config.model_dump())
