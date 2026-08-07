"""ViennaRNA backend (duplexfold/cofold, no GPL) / Backend ViennaRNA"""
import RNA
from labcraft.thermo.vienna import dna_params, _ensure_vienna
from labcraft.thermo.backends.base import DuplexEnergyBackend, DuplexResult

class ViennaRNABackend(DuplexEnergyBackend):
    def __init__(self):
        _ensure_vienna()
        
    def calc_heterodimer(self, seq1: str, seq2: str, *, temp_celsius: float = 65.0, **kwargs) -> DuplexResult:
        with dna_params(temp_celsius):
            seq = f"{seq1}&{seq2}"
            structure, mfe = RNA.cofold(seq)
            return DuplexResult(0.0, 0.0, float(mfe), 0.0, structure, temp_celsius)
            
    def calc_homodimer(self, seq: str, *, temp_celsius: float = 65.0, **kwargs) -> DuplexResult:
        with dna_params(temp_celsius):
            seq2 = f"{seq}&{seq}"
            structure, mfe = RNA.cofold(seq2)
            return DuplexResult(0.0, 0.0, float(mfe), 0.0, structure, temp_celsius)
            
    def calc_hairpin(self, seq: str, *, temp_celsius: float = 65.0, **kwargs) -> DuplexResult:
        with dna_params(temp_celsius):
            structure, mfe = RNA.fold(seq)
            return DuplexResult(0.0, 0.0, float(mfe), 0.0, structure, temp_celsius)
            
    def calc_duplex(self, seq1: str, seq2: str, *, temp_celsius: float = 65.0, **kwargs) -> DuplexResult:
        return self.calc_heterodimer(seq1, seq2, temp_celsius=temp_celsius, **kwargs)
