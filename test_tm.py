from labcraft.thermo.backends.vienna import ViennaRNABackend
backend = ViennaRNABackend()
res_a = backend.calc_duplex("ATCGATCGATCGATCGATCG", "CGATCGATCGATCGATCGAT", temp_celsius=65.0)
res_probe = backend.calc_duplex("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC", "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG", temp_celsius=65.0)
print("a:", res_a.tm_celsius)
print("probe:", res_probe.tm_celsius)
