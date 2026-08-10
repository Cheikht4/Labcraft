from typing import List, Dict, Tuple
from pydantic import BaseModel
from labcraft.lamp.domains import PhysicalPrimer
from labcraft.thermo.backends.base import DuplexEnergyBackend
from labcraft.diagnostics.enzyme import PolymeraseProfile
from labcraft.diagnostics.amplifiable_dimer import is_amplifiable_dimer

def _revcomp(seq: str) -> str:
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'U': 'A',
        'M': 'K', 'R': 'Y', 'W': 'W', 'S': 'S', 'Y': 'R', 'K': 'M',
        'V': 'B', 'H': 'D', 'D': 'H', 'B': 'V',
        'N': 'N'
    }
    try:
        return "".join(complement[c] for c in reversed(seq.upper()))
    except KeyError as e:
        raise ValueError(f"Caractère non reconnu dans la séquence : {e.args[0]}")

class MisprimingRisk(BaseModel):
    primer_name: str
    target_id: str
    delta_g: float
    structure: str
    site_sequence: str

def detect_inter_target_mispriming(
    primers: List[PhysicalPrimer],
    primer_to_panel: Dict[str, str],
    targets: Dict[str, str],
    backend: DuplexEnergyBackend,
    enzyme: PolymeraseProfile,
    temp_celsius: float,
    **backend_kwargs
) -> List[MisprimingRisk]:
    """
    Détecte le mésamorçage inter-cible.
    Ancre sur les 3 nt 3' de l'amorce, extrait le site candidat sur la cible,
    calcule le duplex complet avec calc_heterodimer, puis vérifie le veto 3' ARMS.
    """
    risks = []
    
    # Pour l'ancrage, on utilise les 3 derniers nucléotides
    K_LEN = 3
    
    for p in primers:
        # Résoudre le nom parent pour le lookup (les variants portent nom#N)
        # Resolve parent name for lookup (variants carry name#N)
        p_target = primer_to_panel.get(p.name)
        if not p_target and '#' in p.name:
            p_target = primer_to_panel.get(p.name.split('#')[0])
            
        if len(p.sequence) < K_LEN:
            continue
            
        p_3p = p.sequence[-K_LEN:]
        k_mer = _revcomp(p_3p)
        
        for t_id, t_seq in targets.items():
            t_seq = t_seq.upper()
            if p_target and t_id == p_target:
                continue # Seulement inter-cible
                
            # Chercher k_mer dans les deux brins
            strands = [
                ("+", t_seq),
                ("-", _revcomp(t_seq))
            ]
            
            for strand_dir, s_seq in strands:
                idx = 0
                while True:
                    idx = s_seq.find(k_mer, idx)
                    if idx == -1:
                        break
                        
                    # Site candidat trouvé. Le k-mer de la cible s'apparie aux 3 derniers nt de p.
                    # P s'étend vers son 5', ce qui correspond au 3' du site candidat.
                    # On extrait une fenêtre large autour de l'ancrage.
                    # P (5'->3') est antiparallèle au site.
                    # Le 3' de P s'apparie à s_seq[idx : idx+K_LEN].
                    # Le 5' de P s'appariera aux nucléotides *suivants* dans s_seq.
                    # On extrait donc s_seq de (idx - 5) à (idx + len(p.sequence) + 5).
                    start = max(0, idx - 5)
                    end = min(len(s_seq), idx + len(p.sequence) + 5)
                    site = s_seq[start:end]
                    
                    try:
                        res = backend.calc_heterodimer(
                            p.sequence, site,
                            temp_celsius=temp_celsius,
                            lna_positions_a=p.lna_positions,
                            lna_positions_b=(),
                            **backend_kwargs
                        )
                        
                        # Vérifier si c'est amplifiable par P (ext_strand == 'a')
                        # et si le duplex complet est suffisamment stable
                        struct = res.structure.replace('&', '')
                        mfe = res.dg_kcal
                        
                        is_amp, dg_3p, ext_strand, blocked = is_amplifiable_dimer(
                            p.sequence, site, struct, mfe,
                            enzyme, temp_celsius,
                            blocked_a=p.blocked_3prime, blocked_b=False,
                            lna_positions=p.lna_positions # On passe lna_a, lna_b est vide
                        )
                        
                        if is_amp and ext_strand == 'a' and mfe <= enzyme.dimer_dg_threshold:
                            risks.append(MisprimingRisk(
                                primer_name=p.name,
                                target_id=t_id,
                                delta_g=mfe,
                                structure=res.structure,
                                site_sequence=site
                            ))
                    except Exception:
                        pass
                        
                    idx += 1
                    
    return risks
