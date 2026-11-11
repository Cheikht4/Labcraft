"""Tests structurels obligatoires sur l'architecture LAMP des FIP/BIP.

Mandatory structural tests for LAMP FIP/BIP domain orientation and amplicon ordering.
"""
from pathlib import Path
from typing import Optional, Tuple
import pytest
import regex

from labcraft.cli.parsers import parse_primer_file, read_multi_fasta
from labcraft.lamp.domains import PhysicalPrimer, PrimerRole, _revcomp
from labcraft.target.seeding import primer_matches_sequence, seq_to_regex

GROUND_TRUTH_DIR = Path(__file__).parents[1] / "fixtures" / "ground_truth"
REAL_PRIMER_DIR = Path(__file__).parents[1] / "fixtures" / "real_primer_files"


def _locate_domain(target_seq: str, domain_seq: str) -> Optional[Tuple[int, int, str]]:
    """Localise un domaine sur la cible (substitutions <= 2, puis édition <= 3 pour indels publiés).

    Locates a domain on the target genome (substitutions <= 2, then edit distance <= 3 for published indels).
    """
    hit = primer_matches_sequence(target_seq, domain_seq, max_errors=2, strict_3prime_len=0)
    if hit is not None:
        return hit

    # Fallback édition floue (ex. ParidaD4_F1 avec 1 base surnuméraire publiée sur MG272272.1)
    # Fuzzy edit fallback (e.g. published extra-base ParidaD4_F1 against MG272272.1)
    best_hit: Optional[Tuple[int, int, str]] = None
    best_score = None
    pat_sense = regex.compile(f"(?be)({seq_to_regex(domain_seq)}){{e<=3}}")
    m_s = pat_sense.search(target_seq)
    if m_s is not None:
        score = sum(m_s.fuzzy_counts)
        best_score = score
        best_hit = (m_s.start(), m_s.end(), "+")

    target_rc = _revcomp(target_seq)
    m_rc = pat_sense.search(target_rc)
    if m_rc is not None:
        score = sum(m_rc.fuzzy_counts)
        if best_score is None or score < best_score:
            n_t = len(target_seq)
            best_hit = (n_t - m_rc.end(), n_t - m_rc.start(), "-")

    return best_hit


def _extract_fip_bip_domains(primer_cfg, role: str, target_seq: str) -> Tuple[str, str]:
    """Extrait (queue 5', domaine 3') depuis domains explicites ou par PhysicalPrimer.from_alignment.

    Extracts (5' tail, 3' binding domain) from explicit domains or via PhysicalPrimer.from_alignment.
    """
    tail_key = "F1c" if role == "FIP" else "B1c"
    bind_key = "F2" if role == "FIP" else "B2"
    if primer_cfg.domains is not None:
        tail_seq = getattr(primer_cfg.domains, tail_key, None)
        bind_seq = getattr(primer_cfg.domains, bind_key, None)
        if tail_seq and bind_seq:
            return tail_seq, bind_seq

    role_enum = PrimerRole.FIP if role == "FIP" else PrimerRole.BIP
    pp = PhysicalPrimer.from_alignment(f"{role}_struct", primer_cfg.seq, role_enum, target_seq)
    assert pp.tail_domain and pp.binding_domain, (
        f"{role}: impossible d'extraire queue 5' / domaine 3' depuis '{primer_cfg.seq}'"
    )
    return pp.tail_domain, pp.binding_domain


def _assert_panel_lamp_geometry(primer_set, target_seq: str) -> None:
    """Vérifie la géométrie LAMP des FIP et BIP d'un panel sur son génome cible.

    Checks LAMP strand opposition and amplicon ordering of FIP/BIP on the target genome.
    """
    primers = primer_set.primers
    extracted = {}
    for role in ("FIP", "BIP"):
        if role in primers:
            extracted[role] = _extract_fip_bip_domains(primers[role], role, target_seq)

    # Localiser F2 et B2 pour délimiter l'amplicon interne F2..B2 / Locate F2 and B2 bounds
    f2_domain = extracted["FIP"][1] if "FIP" in extracted else (primers["F2"].seq if "F2" in primers else None)
    b2_domain = extracted["BIP"][1] if "BIP" in extracted else (primers["B2"].seq if "B2" in primers else None)
    assert f2_domain is not None and b2_domain is not None, f"{primer_set.target}: F2/B2 manquants"

    f2_hit = _locate_domain(target_seq, f2_domain)
    b2_hit = _locate_domain(target_seq, b2_domain)
    assert f2_hit is not None, f"{primer_set.target}: F2 non localisé sur la cible"
    assert b2_hit is not None, f"{primer_set.target}: B2 non localisé sur la cible"

    f2_start, f2_end, f2_strand = f2_hit
    b2_start, b2_end, b2_strand = b2_hit
    amplicon_left = min(f2_start, b2_start)
    amplicon_right = max(f2_end, b2_end)

    for role, tail_key, bind_key in [("FIP", "F1c", "F2"), ("BIP", "B1c", "B2")]:
        if role not in extracted:
            continue
        tail_seq, bind_seq = extracted[role]

        tail_hit = _locate_domain(target_seq, tail_seq)
        bind_hit = _locate_domain(target_seq, bind_seq)
        assert tail_hit is not None, f"{primer_set.target} {role}: queue 5' {tail_key} ({tail_seq}) non localisée"
        assert bind_hit is not None, f"{primer_set.target} {role}: domaine 3' {bind_key} ({bind_seq}) non localisé"

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
