from typing import List, Dict, Any
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole, _revcomp
from labcraft.thermo.backends.base import DuplexEnergyBackend

def check_probes_tm(
    primers: List[PhysicalPrimer],
    backend: DuplexEnergyBackend, # Keep signature but we will override it for Tm
    temp_celsius: float,
    **backend_kwargs
) -> List[Dict[str, Any]]:
    """
    Calcule le Tm des amorces et sondes, et vérifie la règle de conception TaqMan :
    Tm_sonde >= max(Tm_amorces) + 5 °C.
    
    Args:
        primers: Liste des oligos (amorces + sondes) du panel
        backend: Moteur de calcul thermodynamique
        temp_celsius: Température de la réaction (pour information/warning, mais Tm est intrinsèque)
        
    Returns:
        Une liste de dictionnaires avec les infos pour l'affichage :
        [{
            "probe_name": "...",
            "probe_tm": float,
            "max_primer_tm": float,
            "delta_tm": float,
            "is_ok": bool,
            "probe_conc": float
        }]
    """
    
    # 1. Grouper par cible/panel
    # Note : actuellement, on reçoit une liste simple, souvent un seul panel.
    # Pour simplifier, on traite tous les oligos passés comme un seul panel.
    
    probes = []
    other_primers = []
    
    for p in primers:
        if p.role == PrimerRole.PROBE:
            probes.append(p)
        else:
            other_primers.append(p)
            
    if not probes:
        return []
        
    # FORCE usage of NativeBackend to ensure LNA rules apply for melting temperature.
    from labcraft.thermo.backends.native import NativeBackend
    from labcraft.thermo.salt import SaltCorrectedBackend, UnifiedSaltModel
    native_backend = SaltCorrectedBackend(NativeBackend(), UnifiedSaltModel())
        
    # 2. Calculer le Tm des amorces
    primer_tms = []
    for p in other_primers:
        # Concentration individuelle
        conc = p.nominal_concentration if p.nominal_concentration else 0.8e-6
        # kwargs locaux
        kw = dict(backend_kwargs)
        kw['ct_molar'] = conc
        # Duplexe parfait (oligo contre complément inverse)
        comp = _revcomp(p.sequence)
        # LNA positions are taken directly from the PhysicalPrimer
        kw['lna_positions'] = p.lna_positions
        res = native_backend.calc_duplex(p.sequence, comp, temp_celsius=temp_celsius, **kw)
        if res.tm_celsius is not None:
            primer_tms.append(res.tm_celsius)
            
    max_primer_tm = max(primer_tms) if primer_tms else 0.0
    
    # 3. Calculer le Tm des sondes et comparer
    results = []
    for p in probes:
        conc = p.nominal_concentration if p.nominal_concentration else 0.2e-6
        kw = dict(backend_kwargs)
        kw['ct_molar'] = conc
        comp = _revcomp(p.sequence)
        kw['lna_positions'] = p.lna_positions
        res = native_backend.calc_duplex(p.sequence, comp, temp_celsius=temp_celsius, **kw)
        probe_tm = res.tm_celsius or 0.0
        
        delta = probe_tm - max_primer_tm
        
        results.append({
            "probe_name": p.name,
            "probe_tm": probe_tm,
            "max_primer_tm": max_primer_tm,
            "delta_tm": delta,
            "is_ok": delta >= 5.0,
            "probe_conc": conc
        })
        
    return results
