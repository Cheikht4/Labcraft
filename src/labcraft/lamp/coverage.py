import csv
from dataclasses import dataclass
from typing import Dict, List, Optional
from enum import Enum

from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.thermo.mismatch import calculate_hybridization_dg, three_prime_extensible
from labcraft.diagnostics.enzyme import PolymeraseProfile

class SiteVerdict(Enum):
    PARFAIT = "parfait"
    TOLERABLE = "tolerable"
    VETO_3P = "veto_3p"
    ABSENT = "absent"

@dataclass
class SiteEvaluation:
    strain_id: str
    primer_name: str
    primer_role: PrimerRole
    verdict: SiteVerdict
    n_mismatches_count: int
    dg_hyb: Optional[float] = None
    first_bad_pos: Optional[int] = None
    severity: Optional[str] = None
    
@dataclass
class StrainVerdict:
    strain_id: str
    evaluations: Dict[str, SiteEvaluation] # key = primer_name
    is_amplifiable_thermo: bool
    is_amplifiable_count: bool
    limiting_primer_thermo: Optional[str] = None

class CoverageAnalyzer:
    def __init__(
        self,
        primers: List[PhysicalPrimer],
        fasta_dict: Dict[str, str],
        backend,
        enzyme: PolymeraseProfile,
        temp_celsius: float = 65.0,
        dg_threshold: float = -6.0,
        max_mismatches_count: int = 2
    ):
        self.primers = primers
        self.primers_by_name = {p.name: p for p in primers}
        self.fasta_dict = fasta_dict
        self.backend = backend
        self.enzyme = enzyme
        self.temp_celsius = temp_celsius
        self.dg_threshold = dg_threshold
        self.max_mismatches_count = max_mismatches_count
        
        # Cache for perfectly matched duplex energy
        self._perfect_dg_cache: Dict[str, float] = {}

    def _get_perfect_dg(self, primer: PhysicalPrimer) -> float:
        if primer.name not in self._perfect_dg_cache:
            comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            perfect_bottom = "".join(comp.get(c, 'N') for c in primer.binding_domain)
            # Use calculate_hybridization_dg which automatically caches base energy via the backend
            dg_hyb, _ = calculate_hybridization_dg(
                primer.binding_domain, perfect_bottom, self.temp_celsius, self.backend
            )
            self._perfect_dg_cache[primer.name] = dg_hyb
        return self._perfect_dg_cache[primer.name]

    def evaluate_site(self, primer: PhysicalPrimer, site_seq: str, strand: str) -> SiteEvaluation:
        # site_seq is extracted from + strand.
        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        if strand == '+':
            bottom_under_top = "".join(comp.get(c, 'N') for c in site_seq)
        else:
            bottom_under_top = site_seq[::-1]
        n_mismatches = sum(1 for a, b in zip(primer.binding_domain, bottom_under_top) if comp.get(a, '') != b)

        if n_mismatches == 0:
            dg_hyb = self._get_perfect_dg(primer)
        else:
            dg_perfect_backend = self._get_perfect_dg(primer)
            
            from labcraft.thermo.mismatch import nn_duplex_energy
            perfect_bottom = "".join(comp.get(c, c) for c in primer.binding_domain)
            _, _, dg_perfect_nn = nn_duplex_energy(primer.binding_domain, perfect_bottom, self.temp_celsius)
            _, _, dg_mismatched_nn = nn_duplex_energy(primer.binding_domain, bottom_under_top, self.temp_celsius)
            
            ddg_mismatch = dg_mismatched_nn - dg_perfect_nn
            dg_hyb = dg_perfect_backend + ddg_mismatch

        extensible, first_bad_pos, severity = three_prime_extensible(
            primer.binding_domain, bottom_under_top, self.enzyme
        )

        if n_mismatches == 0:
            verdict = SiteVerdict.PARFAIT
        elif not extensible:
            verdict = SiteVerdict.VETO_3P
        elif dg_hyb > self.dg_threshold:
            verdict = SiteVerdict.ABSENT
        else:
            verdict = SiteVerdict.TOLERABLE

        return SiteEvaluation(
            strain_id="", # filled by caller
            primer_name=primer.name,
            primer_role=primer.role,
            verdict=verdict,
            n_mismatches_count=n_mismatches,
            dg_hyb=dg_hyb,
            first_bad_pos=first_bad_pos,
            severity=severity
        )

    def analyze_strains(self, csv_records: List[Dict]) -> List[StrainVerdict]:
        # Validate CSV columns
        if not csv_records:
            return []
        
        first_row = csv_records[0]
        expected_cols = {"strain_id", "position", "strand", "n_mismatches"}
        found_cols = set(first_row.keys())
        missing = expected_cols - found_cols
        if missing:
            raise ValueError(f"Le CSV est incomplet. Colonnes attendues : {expected_cols}, trouvées : {found_cols}. Manquant : {missing}")
        if "primer_role" not in found_cols and "primer_name" not in found_cols:
            raise ValueError("Le CSV doit contenir 'primer_role' ou 'primer_name'.")

        # Group records by strain
        records_by_strain = {}
        for r in csv_records:
            s_id = r["strain_id"]
            if s_id not in records_by_strain:
                records_by_strain[s_id] = []
            records_by_strain[s_id].append(r)

        strain_verdicts = []
        for s_id, records in records_by_strain.items():
            if s_id not in self.fasta_dict:
                continue
            genome = self.fasta_dict[s_id].upper()
            
            evaluations = {}
            # We also store the csv n_mismatches to calculate count rule correctly
            csv_mismatches_dict = {}
            for r in records:
                p_name = r.get("primer_role", r.get("primer_name"))
                primer = None
                if p_name in self.primers_by_name:
                    primer = self.primers_by_name[p_name]
                else:
                    for p in self.primers:
                        if p.role.name == p_name or p.role.value == p_name:
                            primer = p
                            break
                
                if not primer:
                    continue
                    
                if "site_seq" in r and r["site_seq"]:
                    site_seq = r["site_seq"]
                else:
                    pos = int(r["position"])
                    site_seq = genome[pos : pos + len(primer.binding_domain)]
                    
                eval_obj = self.evaluate_site(primer, site_seq, r["strand"])
                eval_obj.strain_id = s_id
                evaluations[primer.name] = eval_obj
                csv_mismatches_dict[primer.name] = int(r["n_mismatches"])

            # Check if all init primers are present
            init_primers = [p for p in self.primers if p.role in (PrimerRole.F3, PrimerRole.B3, PrimerRole.FIP, PrimerRole.BIP)]
            
            is_amp_thermo = True
            is_amp_count = True
            limiting = None
            
            for p in init_primers:
                if p.name not in evaluations:
                    is_amp_thermo = False
                    is_amp_count = False
                    if not limiting: limiting = p.name
                    continue
                    
                ev = evaluations[p.name]
                csv_mms = csv_mismatches_dict[p.name]
                
                # Thermo rule
                if ev.verdict in (SiteVerdict.VETO_3P, SiteVerdict.ABSENT):
                    is_amp_thermo = False
                    if not limiting: limiting = p.name
                    
                # Count rule (USING CSV COLUMN)
                if csv_mms > self.max_mismatches_count:
                    is_amp_count = False
                    
            strain_verdicts.append(StrainVerdict(
                strain_id=s_id,
                evaluations=evaluations,
                is_amplifiable_thermo=is_amp_thermo,
                is_amplifiable_count=is_amp_count,
                limiting_primer_thermo=limiting
            ))
            
        return strain_verdicts