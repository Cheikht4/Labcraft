"""Générateur déterministe (graine fixe) des cohortes synthétiques de non-régression pour `labcraft coverage`.
Deterministic (fixed-seed) generator of synthetic non-regression cohorts for `labcraft coverage`.
"""
from __future__ import annotations

import argparse
import hashlib
import random
import sys
from pathlib import Path
from typing import Dict, List, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from labcraft.cli.parsers import read_multi_fasta, reverse_complement  # noqa: E402
from labcraft.lamp.domains import _find_iupac_substring  # noqa: E402

# Graine fixe garantissant la reproductibilité octet par octet des cohortes synthétiques
# Fixed PRNG seed guaranteeing byte-for-byte reproducibility of synthetic cohorts
DEFAULT_SEED: int = 20261119


def build_synthetic_cohort(
    out_fasta: Path,
    panel_tag: str,
    ref_fasta: Path,
    primer_file: Path,
    *,
    seed: int = DEFAULT_SEED,
) -> Tuple[int, str]:
    """Construit de façon déterministe 120 génomes mutants synthétiques dérivés d'un génome de référence.
    Deterministically builds 120 synthetic mutant genomes derived from a reference viral genome.

    Répartition par cohorte (120 séquences) / Cohort breakdown (120 sequences):
    - 10 souches sauvages ou mutées hors amplicon (`WT` + 9 `offsite`)
    - 30 souches à 1 mésappariement interne isolé (`1×1`, hors zone 3' stricte)
    - 15 souches à 1 mésappariement terminal 5'
    - 25 souches à 2 mésappariements internes non adjacents
    - 25 souches à 2 mésappariements internes contigus (boucle interne `2×2`, `n=4`)
    - 10 souches à 3 mésappariements internes contigus (boucle interne `3×3`, `n=6`)
    -  5 souches à 1 mésappariement en position 3' terminale (veto ARMS)

    Returns:
        (n_records, sha256_hex)
    """
    ref_id, ref_seq = read_multi_fasta(str(ref_fasta))[0]
    ref_seq = ref_seq.upper()
    primers = [(h, s.upper()) for h, s in read_multi_fasta(str(primer_file))]
    if panel_tag:
        primers = [
            (h, s)
            for h, s in primers
            if h.startswith(f"{panel_tag}_") or h.endswith(f"_{panel_tag}")
        ]

    # Localiser les sites de liaison sur le génome de référence / Locate binding sites on reference genome
    loci: List[Tuple[str, int, int, str]] = []
    for name, seq in primers:
        sub = seq[-20:] if ("FIP" in name or "BIP" in name) else seq
        pos_p = _find_iupac_substring(sub, ref_seq)
        if pos_p >= 0:
            loci.append((name, pos_p, pos_p + len(sub), "+"))
        else:
            rc_sub = reverse_complement(sub)
            pos_m = _find_iupac_substring(rc_sub, ref_seq)
            if pos_m >= 0:
                loci.append((name, pos_m, pos_m + len(sub), "-"))

    rng = random.Random(seed)
    alt_bases: Dict[str, Tuple[str, str, str]] = {
        "A": ("C", "G", "T"),
        "C": ("A", "G", "T"),
        "G": ("A", "C", "T"),
        "T": ("A", "C", "G"),
    }

    def mutate_pos(seq_list: List[str], idx: int) -> None:
        b = seq_list[idx]
        if b in alt_bases:
            seq_list[idx] = rng.choice(alt_bases[b])

    def safe_internal_positions(start: int, end: int, strand: str) -> List[int]:
        # Éviter les 3 nt côté 3' et les 2 nt côté 5' de l'amorce / Avoid 3'-end 3 nt and 5'-end 2 nt
        if strand == "+":
            return list(range(start + 2, end - 3))
        return list(range(start + 3, end - 2))

    records: List[Tuple[str, str]] = [(f"{ref_id}_WT", ref_seq)]
    for k in range(1, 10):
        s_l = list(ref_seq)
        for _ in range(15):
            mutate_pos(s_l, rng.randrange(100, 1000))
        records.append((f"{ref_id}_offsite_{k:02d}", "".join(s_l)))

    for k in range(30):
        s_l = list(ref_seq)
        lname, st, en, strand = loci[k % len(loci)]
        pos_choices = safe_internal_positions(st, en, strand)
        mutate_pos(s_l, rng.choice(pos_choices))
        records.append((f"{ref_id}_1mm_int_{k:02d}_{lname}", "".join(s_l)))

    for k in range(15):
        s_l = list(ref_seq)
        lname, st, en, strand = loci[k % len(loci)]
        p = st if strand == "+" else (en - 1)
        mutate_pos(s_l, p)
        records.append((f"{ref_id}_1mm_5p_{k:02d}_{lname}", "".join(s_l)))

    for k in range(25):
        s_l = list(ref_seq)
        lname, st, en, strand = loci[k % len(loci)]
        pos_choices = safe_internal_positions(st, en, strand)
        p1 = pos_choices[k % (len(pos_choices) - 3)]
        p2 = pos_choices[min(len(pos_choices) - 1, (k % (len(pos_choices) - 3)) + 3)]
        mutate_pos(s_l, p1)
        mutate_pos(s_l, p2)
        records.append((f"{ref_id}_2mm_sep_{k:02d}_{lname}", "".join(s_l)))

    for k in range(25):
        s_l = list(ref_seq)
        lname, st, en, strand = loci[k % len(loci)]
        pos_choices = safe_internal_positions(st, en, strand)
        p1 = pos_choices[k % (len(pos_choices) - 1)]
        mutate_pos(s_l, p1)
        mutate_pos(s_l, p1 + 1)
        records.append((f"{ref_id}_2mm_adj2x2_{k:02d}_{lname}", "".join(s_l)))

    for k in range(10):
        s_l = list(ref_seq)
        lname, st, en, strand = loci[k % len(loci)]
        pos_choices = safe_internal_positions(st, en, strand)
        p1 = pos_choices[k % (len(pos_choices) - 2)]
        mutate_pos(s_l, p1)
        mutate_pos(s_l, p1 + 1)
        mutate_pos(s_l, p1 + 2)
        records.append((f"{ref_id}_3mm_adj3x3_{k:02d}_{lname}", "".join(s_l)))

    for k in range(5):
        s_l = list(ref_seq)
        lname, st, en, strand = loci[k % len(loci)]
        p = (en - 1) if strand == "+" else st
        mutate_pos(s_l, p)
        records.append((f"{ref_id}_1mm_3pveto_{k:02d}_{lname}", "".join(s_l)))

    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    payload_str = "".join(f">{sid}\n{sq}\n" for sid, sq in records)
    out_fasta.write_text(payload_str, encoding="utf-8")
    digest = hashlib.sha256(payload_str.encode("utf-8")).hexdigest()
    return len(records), digest


def generate_default_cohorts(out_dir: Path, *, seed: int = DEFAULT_SEED) -> Dict[str, Tuple[Path, int, str]]:
    """Génère les deux cohortes synthétiques DENV3 (120) et DENV4 (120) dans `out_dir`.
    Generates both synthetic cohorts DENV3 (120) and DENV4 (120) inside `out_dir`.
    """
    gt_dir = REPO_ROOT / "tests" / "fixtures" / "ground_truth"
    real_dir = REPO_ROOT / "tests" / "fixtures" / "real_primer_files"
    d3_path = out_dir / "synthetic_cohort_denv3_120.fasta"
    d4_path = out_dir / "synthetic_cohort_denv4_120.fasta"

    n_d3, sha_d3 = build_synthetic_cohort(
        d3_path,
        "3",
        gt_dir / "DEN3_M93130.fasta",
        real_dir / "primer_dengue_3.txt",
        seed=seed,
    )
    n_d4, sha_d4 = build_synthetic_cohort(
        d4_path,
        "",
        gt_dir / "Dengue_4_S1028.fasta",
        gt_dir / "parida_DENV4.txt",
        seed=seed,
    )
    return {
        "denv3": (d3_path, n_d3, sha_d3),
        "denv4": (d4_path, n_d4, sha_d4),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate deterministic synthetic mutant FASTA cohorts for coverage benchmarks."
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=REPO_ROOT / "validation" / "generated_cohorts",
        help="Output directory for generated FASTA files.",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Fixed PRNG seed.")
    args = parser.parse_args()

    res = generate_default_cohorts(args.out_dir, seed=args.seed)
    for key, (path, count, sha) in res.items():
        print(f"{key}: {count} sequences -> {path} (sha256={sha})")


if __name__ == "__main__":
    main()
