"""Tests structurels obligatoires sur l'architecture LAMP des FIP/BIP.

Mandatory structural tests for LAMP FIP/BIP domain orientation and amplicon ordering.
"""
from pathlib import Path
import pytest

from labcraft.cli.parsers import parse_primer_file, read_multi_fasta
from labcraft.target.seeding import primer_matches_sequence

GROUND_TRUTH_DIR = Path(__file__).parents[1] / "fixtures" / "ground_truth"
REAL_PRIMER_DIR = Path(__file__).parents[1] / "fixtures" / "real_primer_files"


def _assert_panel_lamp_geometry(primer_set, target_seq: str) -> None:
    """Vérifie la géométrie LAMP des FIP et BIP d'un panel sur son génome cible.

    Checks LAMP strand opposition and amplicon ordering of FIP/BIP on the target genome.
    """
    primers = primer_set.primers
    # Localiser F2 et B2 pour délimiter l'amplicon interne F2..B2 / Locate F2 and B2 bounds
    f2_domain = (
        primers["FIP"].domains.F2
        if "FIP" in primers and primers["FIP"].domains
        else primers["F2"].seq if "F2" in primers else None
    )
    b2_domain = (
        primers["BIP"].domains.B2
        if "BIP" in primers and primers["BIP"].domains
        else primers["B2"].seq if "B2" in primers else None
    )
    assert f2_domain is not None and b2_domain is not None, f"{primer_set.target}: F2/B2 manquants"

    f2_hit = primer_matches_sequence(target_seq, f2_domain, max_errors=2)
    b2_hit = primer_matches_sequence(target_seq, b2_domain, max_errors=2)
    assert f2_hit is not None, f"{primer_set.target}: F2 non localisé sur la cible"
    assert b2_hit is not None, f"{primer_set.target}: B2 non localisé sur la cible"

    f2_start, f2_end, f2_strand = f2_hit
    b2_start, b2_end, b2_strand = b2_hit
    amplicon_left = min(f2_start, b2_start)
    amplicon_right = max(f2_end, b2_end)

    for role, tail_key, bind_key in [("FIP", "F1c", "F2"), ("BIP", "B1c", "B2")]:
        if role not in primers:
            continue
        p_cfg = primers[role]
        assert p_cfg.domains is not None, f"{primer_set.target} {role}: domaines manquants"
        tail_seq = getattr(p_cfg.domains, tail_key, None)
        bind_seq = getattr(p_cfg.domains, bind_key, None)
        assert tail_seq and bind_seq, (
            f"{primer_set.target} {role}: domaines {tail_key}/{bind_key} manquants"
        )

        tail_hit = primer_matches_sequence(target_seq, tail_seq, max_errors=2)
        bind_hit = primer_matches_sequence(target_seq, bind_seq, max_errors=2)
        assert tail_hit is not None, f"{primer_set.target} {role}: queue 5' {tail_key} non localisée"
        assert bind_hit is not None, f"{primer_set.target} {role}: domaine 3' {bind_key} non localisé"

        t_start, t_end, t_strand = tail_hit
        b_start, b_end, b_strand = bind_hit

        # 1) La queue 5' et le domaine 3' doivent être sur des brins opposés
        # 1) 5' tail and 3' binding domain must locate on opposite strands
        assert t_strand != b_strand, (
            f"{primer_set.target} {role}: queue 5' {tail_key} ({t_strand}) et 3' {bind_key} ({b_strand}) "
            f"sur le même brin"
        )

        # 2) La région 1 doit se trouver entre la région 2 et la région 2 de l'autre côté
        # 2) Region 1 must lie strictly between F2 and B2 across the amplicon
        r1_mid = 0.5 * (t_start + t_end)
        assert amplicon_left < r1_mid < amplicon_right, (
            f"{primer_set.target} {role}: région {tail_key} [{t_start}:{t_end}] hors de "
            f"l'amplicon F2..B2 [{amplicon_left}:{amplicon_right}]"
        )
        assert f2_strand == "+" and b2_strand == "-", (
            f"{primer_set.target} {role}: géométrie attendue F2(+) / B2(-), observé {f2_strand}/{b2_strand}"
        )


@pytest.mark.parametrize(
    "primer_file,target_fasta,expected_panels",
    [
        (GROUND_TRUTH_DIR / "primer_dengue_3_panel3.txt", GROUND_TRUTH_DIR / "DEN3_M93130.fasta", 1),
        (REAL_PRIMER_DIR / "primer_dengue_3.txt", REAL_PRIMER_DIR / "DEN3_M93130.fasta", 2),
        (GROUND_TRUTH_DIR / "nunes_YFV.txt", GROUND_TRUTH_DIR / "YFV_17D_17DD.fasta", 1),
        (GROUND_TRUTH_DIR / "parida_DENV4.txt", GROUND_TRUTH_DIR / "Dengue_4_S1028.fasta", 1),
        (GROUND_TRUTH_DIR / "duplex_YFV_DENV4_primers.txt", GROUND_TRUTH_DIR / "duplex_YFV_DENV4_targets.fasta", 2),
    ],
)
def test_ground_truth_and_panel3_lamp_structure(primer_file: Path, target_fasta: Path, expected_panels: int):
    """Vérifie l'architecture LAMP de chaque panel de vérité terrain et du panel 3.

    Verifies LAMP strand opposition and domain order on every ground-truth panel and Panel 3.
    """
    targets = read_multi_fasta(str(target_fasta))
    sets = parse_primer_file(str(primer_file), targets)
    assert len(sets) == expected_panels
    target_map = dict(targets)
    for p_set in sets:
        _assert_panel_lamp_geometry(p_set, target_map[p_set.target])
