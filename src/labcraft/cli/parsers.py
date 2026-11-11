"""Parsers for LabCraft CLI inputs (FASTA, TXT)."""
import re
from typing import List, Tuple, Dict, Optional
from collections import defaultdict

from labcraft.cli.config import PrimerConfig, PrimerSetConfig, PrimerDomains
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
    
    with open(filepath, 'r', encoding="utf-8") as f:
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
            parts = line[1:].strip().split(None, 1)
            current_header = parts[0]
            current_seq_parts = []
            if len(parts) > 1:
                # Sequence is on the same line only if it is a single valid IUPAC token without spaces
                rest = parts[1].strip()
                if " " not in rest and "\t" not in rest and re.match(r'^[A-Za-z]+$', rest) and len(rest) > 3 and set(rest.upper()).issubset(set("ACGTURYMKSWBDHVN")):
                    current_seq_parts.append(rest)
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
        # p1 est F1 ou B1 ; p2 est F2 ou B2.
        # Retourne (sequence_composite, domaine_queue_5p)
        # p1 is F1 or B1; p2 is F2 or B2.
        # Returns (composite_seq, tail_domain_seq)
        if not target_seq:
            raise ParseError(
                f"Panel '{panel_name}': impossible de vérifier l'orientation de {role_prefix}1 "
                f"sans séquence cible. Veuillez fournir la séquence complète de l'amorce "
                f"(>{panel_name}_{role_prefix}IP)."
            )

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

        if p1_matches and p1rc_matches:
            raise ParseError(
                f"Panel '{panel_name}': Unable to determine orientation for {role_prefix}1. "
                f"Both orientations match the target '{panel_name}'. "
                f"Veuillez fournir la séquence complète de l'amorce (>{panel_name}_{role_prefix}IP)."
            )

        if not p1_matches and not p1rc_matches:
            raise ParseError(
                f"Panel '{panel_name}': l'orientation de {role_prefix}1 n'a pas pu être vérifiée "
                f"sur la cible (aucune correspondance). Veuillez fournir la séquence complète de "
                f"l'amorce (>{panel_name}_{role_prefix}IP) plutôt que les sous-domaines."
            )

        # Pour F (FIP = 5'-F1c-F2-3') : F2 est sur le brin +, donc la queue 5' F1c doit être sur le brin - (RC).
        # For F (FIP = 5'-F1c-F2-3'): F2 is on + strand, so 5' tail F1c must be on - strand (RC).
        # Pour B (BIP = 5'-B1c-B2-3') : B2 est sur le brin -, donc la queue 5' B1c doit être sur le brin + (opposé à B2).
        # For B (BIP = 5'-B1c-B2-3'): B2 is on - strand, so 5' tail B1c must be on + strand (opposite to B2).
        if role_prefix == 'F':
            if p1_matches:
                # p1 sur brin + -> inverser-complémenter pour obtenir F1c (-)
                p1c = reverse_complement(p1)
                return p1c + linker + p2, p1c
            # p1 déjà sur brin - (F1c)
            return p1 + linker + p2, p1

        # role_prefix == 'B' : inversé par rapport à F / opposite rule to F
        if p1_matches:
            # p1 déjà sur brin + (B1c, opposé à B2 sur -)
            return p1 + linker + p2, p1
        # p1 sur brin - -> inverser-complémenter pour obtenir B1c (+)
        p1c = reverse_complement(p1)
        return p1c + linker + p2, p1c

    domains_created: Dict[str, PrimerDomains] = {}

    if 'FIP' in primers:
        # La séquence complète >Panel_FIP fait foi sans aucune reconstruction
        # Full >Panel_FIP sequence takes precedence without reconstruction
        primers.pop('F1', None)
        primers.pop('F2', None)
    elif 'F1' in primers and 'F2' in primers:
        f1_seq = primers['F1']
        f2_seq = primers['F2']
        fip_seq, f1c_seq = check_orientation(f1_seq, f2_seq, 'F')
        primers['FIP'] = fip_seq
        domains_created['FIP'] = PrimerDomains(F2=f2_seq, F1c=f1c_seq, linker=linker)
        del primers['F1']
        del primers['F2']

    if 'BIP' in primers:
        # La séquence complète >Panel_BIP fait foi sans aucune reconstruction
        # Full >Panel_BIP sequence takes precedence without reconstruction
        primers.pop('B1', None)
        primers.pop('B2', None)
    elif 'B1' in primers and 'B2' in primers:
        b1_seq = primers['B1']
        b2_seq = primers['B2']
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
            warnings.warn(
                f"Multiple versions detected for panel '{panel_name}'. They will be treated as separate alternative panels.",
                stacklevel=2,
            )
            
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
                    # Fallback par alignement de séquence des domaines d'initiation (F3, B3, F2, B2)
                    # Sequence-alignment fallback using initiation domains (F3, B3, F2, B2)
                    sample_roles = next(iter(versions.values()), {})
                    probe_seqs = [
                        sample_roles[r] for r in ("F3", "B3", "F2", "B2") if r in sample_roles
                    ]
                    best_t_name = None
                    best_hits = 0
                    for t_name, t_seq in targets:
                        hits = sum(
                            1 for s in probe_seqs
                            if _find_iupac_substring(s, t_seq) >= 0
                            or _find_iupac_substring(reverse_complement(s), t_seq) >= 0
                        )
                        if hits > best_hits:
                            best_hits = hits
                            best_t_name = t_name
                    if best_hits >= 2:
                        target_name = best_t_name
                if target_name is None:
                    if allow_unmatched_targets:
                        target_name = targets[0][0] if targets else panel_name
                    else:
                        raise ParseError(f"Panel '{panel_name}' does not match any target name. Available targets: {', '.join(t_names)}. If there are multiple targets, the panel name must match the target name exactly or by inclusion.")
                
        matched_targets = [
            (t_name, t_seq) for t_name, t_seq in targets if t_name == target_name
        ] or targets
        for v_key, roles_dict in versions.items():
            domains_map = reconstruct_fip_bip(panel_name, roles_dict, matched_targets, linker)
            
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


from labcraft.cli.config import PanelConfig, BufferConfig, TargetConfig
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
        with open(config_path, "r", encoding="utf-8") as f:
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
    elif not config.experiment.buffer:
        # Aligner analyze sur les mêmes conditions par défaut que coverage (50 mM Na, 8 mM Mg, 1,4 mM dNTP)
        # Align analyze with the same default buffer conditions as coverage (50 mM Na, 8 mM Mg, 1.4 mM dNTP)
        config.experiment.buffer = BufferConfig(na_mM=50.0, mg_mM=8.0, dntp_mM=1.4)

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
