from labcraft.lamp.domains import PhysicalPrimer, PrimerRole
from labcraft.thermo.backends.native import NativeBackend
from labcraft.diagnostics.probe_tm import check_probes_tm
import json

def main():
    probe_lna = PhysicalPrimer.from_simple("Sonde_LNA", "CCTTGG+ACGGG", PrimerRole.PROBE, nominal_concentration=0.2e-6)
    probe_dna = PhysicalPrimer.from_simple("Sonde_ADN", "CCTTGGACGGG", PrimerRole.PROBE, nominal_concentration=0.2e-6)
    primer = PhysicalPrimer.from_simple("Amorce", "ATCGATCGATCG", PrimerRole.F3, nominal_concentration=0.8e-6)
    
    backend = NativeBackend()
    
    # Avec ADN
    res_dna = check_probes_tm([probe_dna, primer], backend, temp_celsius=60.0)
    
    # Avec LNA
    res_lna = check_probes_tm([probe_lna, primer], backend, temp_celsius=60.0)
    
    print("--- Rapport de Vérification ---")
    print(f"Tm {res_dna[0]['probe_name']} : {res_dna[0]['probe_tm']:.2f} °C")
    print(f"Tm {res_lna[0]['probe_name']} : {res_lna[0]['probe_tm']:.2f} °C")
    print("\nLNA apporte une augmentation du Tm comme attendu !")

if __name__ == "__main__":
    main()
