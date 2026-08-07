import os
import csv
import sys

def main():
    try:
        from Bio.SeqUtils import MeltingTemp as mt
    except ImportError:
        print("Biopython is not installed. Please install it to generate the tables.")
        sys.exit(1)

    out_dir = os.path.join(os.path.dirname(__file__), "..", "src", "labcraft", "data", "thermo")
    os.makedirs(out_dir, exist_ok=True)

    internal_csv = os.path.join(out_dir, "mismatch_nn_internal.csv")
    terminal_csv = os.path.join(out_dir, "mismatch_nn_terminal.csv")

    # DNA_IMM1: internal mismatches
    with open(internal_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["nn_step_top_bottom", "dH_kcal_per_mol", "dS_cal_per_molK", "source"])
        for key, (dH, dS) in mt.DNA_IMM1.items():
            writer.writerow([key, dH, dS, "Allawi and SantaLucia 1997-1998, Peyret 1999"])

    # DNA_TMM1: terminal mismatches
    with open(terminal_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["nn_step_top_bottom", "dH_kcal_per_mol", "dS_cal_per_molK", "source"])
        for key, (dH, dS) in mt.DNA_TMM1.items():
            writer.writerow([key, dH, dS, "SantaLucia and Hicks 2004"])

    print(f"Generated {internal_csv}")
    print(f"Generated {terminal_csv}")

if __name__ == "__main__":
    main()
