"""Redessin d'une amorce externe ou de boucle LAMP (commande `labcraft redesign`).
Redesign of an outer or loop LAMP primer (`labcraft redesign` command).

Périmètre : F3, B3, LF, LB. Les amorces internes FIP/BIP sont refusées.
Scope: F3, B3, LF, LB. Inner primers FIP/BIP are refused.
"""
from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from labcraft.diagnostics.amplifiable_dimer import evaluate_pair_amplifiable
from labcraft.diagnostics.mispriming import detect_inter_target_mispriming
from labcraft.lamp.complex_enumeration import ComplexInfo, enumerate_complexes
from labcraft.lamp.domains import (
    PhysicalPrimer,
    PrimerRole,
    _find_iupac_substring,
    expand_degenerate,
)
from labcraft.lamp.stoichiometry import ConcentrationProfile, LAMP_DEFAULT_PROFILE
from labcraft.metrics.balance import calculate_multiplex_balance
from labcraft.metrics.fractions import compute_fractions
from labcraft.metrics.risk import evaluate_risks, summarize_amplifiable_dimers
from labcraft.metrics.verdict import (
    DEFAULT_DIMER_WARN_RATIO,
    DEFAULT_ROLE_IMBALANCE_RATIO,
    generate_verdict,
)
from labcraft.solver.dual import solve_dual
from labcraft.solver.types import EquilibriumProblem
from labcraft.target.seeding import primer_matches_sequence
from labcraft.thermo.lna import parse_lna_sequence

REDESIGNABLE_ROLES = ("F3", "B3", "LF", "LB")
REFUSED_INNER_ROLES = ("FIP", "BIP", "F1", "F2", "B1", "B2")
UNFOLDING_WINDOW = 150  # identique à enumerate_complexes / same as enumerate_complexes
R_KCAL = 0.0019872

_IUPAC_GC_WEIGHT: Dict[str, float] = {
    "G": 1.0, "C": 1.0, "S": 1.0,
    "A": 0.0, "T": 0.0, "U": 0.0, "W": 0.0,
    "R": 0.5, "Y": 0.5, "K": 0.5, "M": 0.5, "N": 0.5,
    "B": 2.0 / 3.0, "V": 2.0 / 3.0,
    "D": 1.0 / 3.0, "H": 1.0 / 3.0,
}


class RedesignError(ValueError):
    """Erreur de redessin (rôle refusé, locus introuvable, géométrie non vérifiable).
    Redesign error (refused role, locus not found, unverifiable geometry).
    """


def revcomp(seq: str) -> str:
    """Complément inverse IUPAC / IUPAC reverse complement."""
    comp = {
        "A": "T", "T": "A", "C": "G", "G": "C",
        "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
        "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
    }
    return "".join(comp.get(c, c) for c in reversed(seq.upper()))


def gc_percentage(seq: str) -> float:
    """Calcule le pourcentage GC (0–100) d'une séquence nucléotidique (supporte IUPAC).
    Computes GC percentage (0–100) of a nucleotide sequence (supports IUPAC codes).
    """
    s = seq.upper()
    if not s:
        return 0.0
    return 100.0 * sum(_IUPAC_GC_WEIGHT.get(c, 0.0) for c in s) / len(s)


@dataclass
class Locus:
    """Locus sur le brin + de la cible : [start, end) et brin de l'amorce.
    Locus on the + strand of the target: [start, end) and primer strand.
    """
    start: int
    end: int
    strand: str


@dataclass
class Candidate:
    """Candidat de redessin et ses métriques.
    Redesign candidate and its metrics.
    """
    shift: int
    length: int
    sequence: str
    locus: Locus
    raw_sequence: Optional[str] = None
    is_reference: bool = False
    geometry_ok: bool = True
    gc_percent: float = 0.0
    dg_duplex: float = 0.0
    dg_f2: Optional[float] = None
    ddg_vs_f2: Optional[float] = None
    more_stable_than_f2: bool = False
    ref_inner_domain: str = "F2"
    dg_ref_inner: Optional[float] = None
    ddg_vs_ref_inner: Optional[float] = None
    more_stable_than_ref_inner: bool = False
    dg_open: float = 0.0
    max_3p_comp: int = 0
    max_3p_partner: str = ""
    max_reverse_3p_comp: int = 0
    max_reverse_3p_partner: str = ""
    max_self_3p_comp: int = 0
    n_amplifiable_pairs: int = 0
    prescreen_score: float = 0.0
    amplifiable_partners: List[str] = field(default_factory=list)
    equilibrium: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        """Sérialisation JSON / JSON serialization."""
        return {
            "shift": self.shift,
            "length": self.length,
            "gc_percent": round(self.gc_percent, 2),
            "sequence": self.sequence,
            "start": self.locus.start,
            "end": self.locus.end,
            "strand": self.locus.strand,
            "is_reference": self.is_reference,
            "geometry_ok": self.geometry_ok,
            "dg_duplex": self.dg_duplex,
            "dg_f2": self.dg_f2,
            "ddg_vs_f2": self.ddg_vs_f2,
            "more_stable_than_f2": self.more_stable_than_f2,
            "ref_inner_domain": self.ref_inner_domain,
            "dg_ref_inner": self.dg_ref_inner,
            "ddg_vs_ref_inner": self.ddg_vs_ref_inner,
            "more_stable_than_ref_inner": self.more_stable_than_ref_inner,
            "dg_open": self.dg_open,
            "max_3p_complementarity": self.max_3p_comp,
            "max_3p_partner": self.max_3p_partner,
            "max_reverse_3p_complementarity": self.max_reverse_3p_comp,
            "max_reverse_3p_partner": self.max_reverse_3p_partner,
            "max_self_3p_complementarity": self.max_self_3p_comp,
            "n_amplifiable_pairs": self.n_amplifiable_pairs,
            "prescreen_score": self.prescreen_score,
            "amplifiable_partners": self.amplifiable_partners,
            "equilibrium": self.equilibrium,
        }


def check_role_redesignable(role: str, primer_name: str) -> None:
    """Refuse clairement les amorces hors périmètre (FIP/BIP en particulier).
    Clearly refuses out-of-scope primers (FIP/BIP in particular).
    """
    r = role.upper()
    if r in REFUSED_INNER_ROLES:
        raise RedesignError(
            f"Redessin refusé pour '{primer_name}' (rôle {r}) : les amorces internes FIP/BIP "
            f"ne sont pas prises en charge dans cette version. Seules les amorces externes "
            f"(F3, B3) et de boucle (LF, LB) peuvent être redessinées."
        )
    if r not in REDESIGNABLE_ROLES:
        raise RedesignError(
            f"Redessin refusé pour '{primer_name}' (rôle {r}) : seules les amorces externes "
            f"(F3, B3) et de boucle (LF, LB) peuvent être redessinées."
        )


def locate_on_target(
    seq: str,
    target: str,
    anchor_locus: Optional[Locus] = None,
) -> Optional[Locus]:
    """Localise une séquence (amorce ou domaine) sur la cible : exacte IUPAC, substitutions <= 2, puis édition locale <= 3.
    Locates a sequence (primer or domain) on the target: exact IUPAC, <= 2 substitutions, then local <= 3 edits.

    Brin '+' : la séquence est identique au brin + ; '-' : son complément inverse est sur le brin +.
    Strand '+': the sequence equals the + strand; '-': its reverse complement lies on the + strand.
    """
    if not seq:
        return None
    pos = _find_iupac_substring(seq, target)
    if pos >= 0:
        return Locus(pos, pos + len(seq), "+")
    pos = _find_iupac_substring(revcomp(seq), target)
    if pos >= 0:
        return Locus(pos, pos + len(seq), "-")
    res = primer_matches_sequence(target, seq, max_errors=2, strict_3prime_len=0)
    if res is not None:
        start, end, strand = res
        return Locus(start, end, strand)

    # Repli local à distance d'édition <= 3 autour de l'ancre d'amplicon (identique à PhysicalPrimer.from_alignment,
    # ex. polymorphisme + 1 indel de souche clinique dans F1/B1 comme Dengue_4_S1028).
    # Local edit-distance <= 3 fallback around amplicon anchor (matching PhysicalPrimer.from_alignment,
    # e.g. clinical strain polymorphism + 1 indel in F1/B1 such as Dengue_4_S1028).
    if anchor_locus is not None and target:
        import regex as _regex
        from labcraft.target.seeding import seq_to_regex

        w_start = max(0, anchor_locus.start - 500)
        w_end = min(len(target), anchor_locus.end + 500)
        local_sense = target[w_start:w_end].upper()
        pat_plus = seq_to_regex(seq.upper())
        pat_minus = seq_to_regex(revcomp(seq.upper()))
        # En LAMP, la queue F1c/B1c est sur le brin opposé au domaine 3' F2/B2 : tester ce brin en priorité
        # In LAMP, the F1c/B1c tail lies on the opposite strand from F2/B2: probe opposite strand first
        probe_order = (
            [("-", pat_minus), ("+", pat_plus)]
            if anchor_locus.strand == "+"
            else [("+", pat_plus), ("-", pat_minus)]
        )
        for cand_strand, pat in probe_order:
            m = _regex.search(f"(?be)({pat}){{e<=3}}", local_sense)
            if m is not None:
                return Locus(w_start + m.start(), w_start + m.end(), cand_strand)

    return None


def panel_geometry(panel_primers: Dict[str, PhysicalPrimer], target: str) -> Dict[str, Optional[Locus]]:
    """Localise F2, F1 (queue FIP), B2, B1 (queue BIP) du panel sur la cible.
    Locates the panel's F2, F1 (FIP tail), B2, B1 (BIP tail) on the target.
    """
    geo: Dict[str, Optional[Locus]] = {"F2": None, "F1": None, "B2": None, "B1": None}
    fip = panel_primers.get("FIP")
    bip = panel_primers.get("BIP")
    if fip is not None:
        geo["F2"] = locate_on_target(fip.binding_domain, target)
        geo["F1"] = locate_on_target(fip.tail_domain or "", target, anchor_locus=geo["F2"])
    if bip is not None:
        geo["B2"] = locate_on_target(bip.binding_domain, target)
        geo["B1"] = locate_on_target(bip.tail_domain or "", target, anchor_locus=geo["B2"])
    return geo


def geometry_satisfied(role: str, cand: Locus, geo: Dict[str, Optional[Locus]]) -> bool:
    """Vérifie les contraintes géométriques LAMP en coordonnées du brin +.
    Checks LAMP geometric constraints in + strand coordinates.

    Orientation directe : F3 < F2 < F1 < B1 < B2 < B3 ; l'orientation inverse est le miroir.
    Forward orientation: F3 < F2 < F1 < B1 < B2 < B3; reverse orientation is the mirror image.
    """
    f2, b2, f1, b1 = geo["F2"], geo["B2"], geo["F1"], geo["B1"]
    if f2 is None or b2 is None:
        raise RedesignError(
            "Géométrie LAMP non vérifiable : F2 ou B2 introuvable sur la cible. "
            "Fournissez les séquences complètes FIP/BIP du panel."
        )
    forward = f2.start < b2.start
    r = role.upper()
    if r == "F3":
        return cand.end <= f2.start if forward else cand.start >= f2.end
    if r == "B3":
        return cand.start >= b2.end if forward else cand.end <= b2.start
    if r == "LF":
        if f1 is None:
            raise RedesignError("Géométrie LAMP non vérifiable : région F1 (queue FIP) introuvable sur la cible.")
        return (cand.start >= f2.end and cand.end <= f1.start) if forward else (cand.end <= f2.start and cand.start >= f1.end)
    if r == "LB":
        if b1 is None:
            raise RedesignError("Géométrie LAMP non vérifiable : région B1 (queue BIP) introuvable sur la cible.")
        return (cand.start >= b1.end and cand.end <= b2.start) if forward else (cand.end <= b1.start and cand.start >= b2.end)
    return False


def generate_candidates(
    original_seq: str,
    original_locus: Locus,
    target: str,
    shift: int,
    length_min: int,
    length_max: int,
    raw_original_seq: Optional[str] = None,
) -> List[Candidate]:
    """Génère les fenêtres glissées sur le même brin (décalage de l'extrémité 3').
    Generates sliding windows on the same strand (3' end shift).

    Décalage positif : l'extrémité 3' avance dans le sens de l'amorce.
    Positive shift: the 3' end moves forward in the primer's direction.
    La ligne de référence (décalage 0, longueur d'origine, séquence d'origine) est toujours incluse.
    The reference row (shift 0, original length, original sequence) is always included.
    """
    orig_len = len(original_seq)
    lengths = sorted(set(range(length_min, length_max + 1)) | {orig_len})
    orig_up = original_seq.upper()
    out: List[Candidate] = [
        Candidate(
            0, orig_len, orig_up, original_locus,
            raw_sequence=raw_original_seq or orig_up,
            is_reference=True, gc_percent=gc_percentage(orig_up),
        )
    ]
    n = len(target)
    for d in range(-shift, shift + 1):
        for L in lengths:
            if d == 0 and L == orig_len:
                continue
            if original_locus.strand == "+":
                end = original_locus.end + d
                start = end - L
            else:
                start = original_locus.start - d
                end = start + L
            if start < 0 or end > n:
                continue
            window = target[start:end].upper()
            if any(c not in "ACGT" for c in window):
                # Positions ambiguës ou inconnues dans la cible : candidat non commandable
                # Ambiguous or unknown target positions: non-orderable candidate
                continue
            seq = window if original_locus.strand == "+" else revcomp(window)
            out.append(
                Candidate(
                    d, L, seq, Locus(start, end, original_locus.strand),
                    raw_sequence=seq, gc_percent=gc_percentage(seq),
                )
            )
    return out


def compute_panel_domain_dg(
    engine_primers: List[PhysicalPrimer],
    panel_id: str,
    domain_name: str,
    backend: Any,
    temp_celsius: float,
    backend_kwargs: Dict[str, Any],
) -> Optional[float]:
    """Calcule le ΔG du duplex parfait d'un sous-domaine interne (F2, B2, F1c, B1c) du même jeu d'amorces.
    Computes the perfect-duplex ΔG of an inner sub-domain (F2, B2, F1c, B1c) of the same primer set.
    """
    dom_upper = domain_name.upper()
    if dom_upper in ("F2", "F1C"):
        target_parent = f"FIP_{panel_id}"
        use_tail = dom_upper == "F1C"
    elif dom_upper in ("B2", "B1C"):
        target_parent = f"BIP_{panel_id}"
        use_tail = dom_upper == "B1C"
    else:
        return None

    seen: set[Tuple[str, Tuple[int, ...]]] = set()
    dgs: List[float] = []
    for p in engine_primers:
        par = p.parent_name or p.name.split("#")[0]
        if par != target_parent:
            continue
        sub_seq = (p.tail_domain or "") if use_tail else p.binding_domain
        if not sub_seq:
            continue
        offset = p.sequence.find(sub_seq)
        sub_lna = (
            tuple(pos - offset for pos in p.lna_positions if offset <= pos < offset + len(sub_seq))
            if offset >= 0
            else ()
        )
        key = (sub_seq, sub_lna)
        if key in seen:
            continue
        seen.add(key)
        res = backend.calc_duplex(
            sub_seq, revcomp(sub_seq), temp_celsius=temp_celsius,
            lna_positions_a=sub_lna, lna_positions_b=(), **backend_kwargs,
        )
        dgs.append(float(res.dg_kcal))
    return min(dgs) if dgs else None


def compute_panel_f2_dg(
    engine_primers: List[PhysicalPrimer],
    panel_id: str,
    backend: Any,
    temp_celsius: float,
    backend_kwargs: Dict[str, Any],
) -> Optional[float]:
    """Calcule le ΔG du duplex parfait du domaine F2 (liaison de FIP) du même jeu d'amorces.
    Computes the perfect-duplex ΔG of the F2 domain (FIP binding domain) of the same primer set.
    """
    return compute_panel_domain_dg(engine_primers, panel_id, "F2", backend, temp_celsius, backend_kwargs)


def resolve_role_inner_reference(
    role: str,
    engine_primers: List[PhysicalPrimer],
    panel_id: str,
    backend: Any,
    temp_celsius: float,
    backend_kwargs: Dict[str, Any],
) -> Tuple[str, Optional[float], Optional[float]]:
    """Détermine le sous-domaine interne de référence adapté au rôle redessiné.
    Determines the homologous inner reference domain appropriate for the redesigned role.

    - F3 -> F2 (domaine 3' de FIP : F3 ne doit pas précéder FIP)
    - B3 -> B2 (domaine 3' de BIP : B3 ne doit pas précéder BIP ; repli F2 si B2 absent)
    - LF -> F1c (queue 5' de FIP : la boucle LF reste sous l'ancrage de tige F1c ; repli F2)
    - LB -> B1c (queue 5' de BIP : la boucle LB reste sous l'ancrage de tige B1c ; repli B2/F2)

    Retourne (nom_domaine_reference, dg_reference, dg_f2).
    Returns (ref_domain_name, dg_ref, dg_f2).
    """
    dg_f2 = compute_panel_domain_dg(engine_primers, panel_id, "F2", backend, temp_celsius, backend_kwargs)
    r = role.upper()
    if r == "F3":
        return "F2", dg_f2, dg_f2
    if r == "B3":
        dg_b2 = compute_panel_domain_dg(engine_primers, panel_id, "B2", backend, temp_celsius, backend_kwargs)
        return ("B2", dg_b2, dg_f2) if dg_b2 is not None else ("F2", dg_f2, dg_f2)
    if r == "LF":
        dg_f1c = compute_panel_domain_dg(engine_primers, panel_id, "F1C", backend, temp_celsius, backend_kwargs)
        return ("F1c", dg_f1c, dg_f2) if dg_f1c is not None else ("F2", dg_f2, dg_f2)
    if r == "LB":
        dg_b1c = compute_panel_domain_dg(engine_primers, panel_id, "B1C", backend, temp_celsius, backend_kwargs)
        if dg_b1c is not None:
            return "B1c", dg_b1c, dg_f2
        dg_b2 = compute_panel_domain_dg(engine_primers, panel_id, "B2", backend, temp_celsius, backend_kwargs)
        return ("B2", dg_b2, dg_f2) if dg_b2 is not None else ("F2", dg_f2, dg_f2)
    return "F2", dg_f2, dg_f2


def compute_f2_comparison(cand: Candidate, dg_f2: Optional[float]) -> None:
    """Compare le ΔG du candidat à celui du domaine F2 du même jeu (compatibilité unitaire).
    Compares candidate duplex ΔG against the same panel's F2 domain ΔG (unit-test compatibility).
    """
    compute_inner_comparison(cand, "F2", dg_f2, dg_f2)


def compute_inner_comparison(
    cand: Candidate,
    ref_domain: str,
    dg_ref: Optional[float],
    dg_f2: Optional[float],
) -> None:
    """Compare le ΔG du candidat à son domaine interne homologue (`ref_domain`) ainsi qu'à `F2`.
    Compares candidate duplex ΔG against its homologous inner domain (`ref_domain`) as well as `F2`.
    """
    cand.dg_f2 = dg_f2
    if dg_f2 is None:
        cand.ddg_vs_f2 = None
        cand.more_stable_than_f2 = False
    else:
        cand.ddg_vs_f2 = float(cand.dg_duplex - dg_f2)
        cand.more_stable_than_f2 = bool(cand.dg_duplex < dg_f2 - 1e-9)

    cand.ref_inner_domain = ref_domain
    cand.dg_ref_inner = dg_ref
    if dg_ref is None:
        cand.ddg_vs_ref_inner = None
        cand.more_stable_than_ref_inner = False
    else:
        cand.ddg_vs_ref_inner = float(cand.dg_duplex - dg_ref)
        cand.more_stable_than_ref_inner = bool(cand.dg_duplex < dg_ref - 1e-9)


def compute_max_3p_complementarity(cand: Candidate, partners: List[PhysicalPrimer]) -> None:
    """Calcule les complémentarités 3' parfaites avec matrice 5' (directe, inverse et homodimère).
    Computes 5'-templated perfect 3' complementarities (forward, reverse partner->cand, and self).

    - `max_3p_comp` : plus long segment 3' du candidat s'hybridant à l'intérieur (`P[1:]`) d'un partenaire.
    - `max_reverse_3p_comp` : plus long segment 3' d'un partenaire s'hybridant à l'intérieur (`cand[1:]`) du candidat.
    - `max_self_3p_comp` : plus long segment 3' du candidat s'hybridant à l'intérieur (`cand[1:]`) de lui-même.
    """
    cand.gc_percent = gc_percentage(cand.sequence)
    variants = expand_degenerate(cand.sequence)

    # 1. Sens direct : 3' du candidat -> matrice 5' du partenaire
    # 1. Forward direction: candidate 3' end -> partner 5' template
    best_k = 0
    best_partner = ""
    for v in variants:
        v_len = len(v)
        for p in partners:
            p_seq = p.sequence
            max_k = min(v_len, len(p_seq) - 1)
            if max_k <= best_k:
                continue
            p_tail_space = p_seq[1:]
            if revcomp(v[-(best_k + 1):]) not in p_tail_space:
                continue
            for k in range(best_k + 1, max_k + 1):
                if revcomp(v[-k:]) in p_tail_space:
                    best_k = k
                    best_partner = p.parent_name or p.name.split("#")[0]
                else:
                    break

    # 2. Sens inverse : 3' d'un partenaire -> matrice 5' offerte par le candidat
    # 2. Reverse direction: partner 3' end -> 5' template offered by the candidate
    best_rev_k = 0
    best_rev_partner = ""
    for v in variants:
        if len(v) <= 1:
            continue
        v_tail_space = v[1:]
        for p in partners:
            p_seq = p.sequence
            max_k = min(len(p_seq), len(v) - 1)
            if max_k <= best_rev_k:
                continue
            if revcomp(p_seq[-(best_rev_k + 1):]) not in v_tail_space:
                continue
            for k in range(best_rev_k + 1, max_k + 1):
                if revcomp(p_seq[-k:]) in v_tail_space:
                    best_rev_k = k
                    best_rev_partner = p.parent_name or p.name.split("#")[0]
                else:
                    break

    # 3. Auto-amorçage homodimérique : 3' du candidat -> matrice 5' du candidat
    # 3. Homodimeric self-priming: candidate 3' end -> candidate 5' template
    best_self_k = 0
    for v1 in variants:
        for v2 in variants:
            if len(v2) <= 1:
                continue
            v2_tail = v2[1:]
            max_k = min(len(v1), len(v2) - 1)
            if max_k <= best_self_k:
                continue
            if revcomp(v1[-(best_self_k + 1):]) not in v2_tail:
                continue
            for k in range(best_self_k + 1, max_k + 1):
                if revcomp(v1[-k:]) in v2_tail:
                    best_self_k = k
                else:
                    break

    cand.max_3p_comp = best_k
    cand.max_3p_partner = best_partner
    cand.max_reverse_3p_comp = best_rev_k
    cand.max_reverse_3p_partner = best_rev_partner
    cand.max_self_3p_comp = best_self_k


def compute_duplex_energy(
    cand: Candidate,
    backend: Any,
    temp_celsius: float,
    backend_kwargs: Dict[str, Any],
) -> None:
    """Calcule uniquement le ΔG du duplex parfait du candidat sur cible.
    Computes only the candidate's perfect-duplex ΔG on target.
    """
    raw_seq = cand.raw_sequence or cand.sequence
    dgs = []
    for v_raw in expand_degenerate(raw_seq):
        v_bare, lna_pos = parse_lna_sequence(v_raw)
        res = backend.calc_duplex(
            v_bare, revcomp(v_bare), temp_celsius=temp_celsius,
            lna_positions_a=tuple(lna_pos), lna_positions_b=(), **backend_kwargs,
        )
        dgs.append(res.dg_kcal)
    cand.dg_duplex = float(min(dgs))


def compute_unfolding_energy(
    cand: Candidate,
    role: str,
    target: str,
    temp_celsius: float,
    mon_molar: Optional[float],
) -> None:
    """Calcule le ΔG d'ouverture locale de la cible (même fenêtre que `analyze`).
    Computes local target unfolding ΔG (same window as `analyze`).
    """
    from labcraft.target.unfolding import calc_unfolding_penalty

    if role.upper() in ("LF", "LB"):
        # Comme analyze : pas de pénalité d'ouverture pour les boucles
        # Same as analyze: no opening penalty for loop primers
        cand.dg_open = 0.0
        return
    s0, e0 = cand.locus.start, cand.locus.end
    w_start = max(0, s0 - UNFOLDING_WINDOW)
    w_end = min(len(target), e0 + UNFOLDING_WINDOW)
    cand.dg_open = float(calc_unfolding_penalty(
        target[w_start:w_end], s0 - w_start, e0 - w_start,
        temp_celsius=temp_celsius, mon_molar=mon_molar,
    ))


def compute_site_energies(
    cand: Candidate,
    role: str,
    target: str,
    backend: Any,
    temp_celsius: float,
    mon_molar: Optional[float],
    backend_kwargs: Dict[str, Any],
) -> None:
    """ΔG du duplex parfait et ΔG d'ouverture de la cible (même fenêtre que `analyze`).
    Perfect-duplex ΔG and target opening ΔG (same window as `analyze`).
    """
    compute_duplex_energy(cand, backend, temp_celsius, backend_kwargs)
    compute_unfolding_energy(cand, role, target, temp_celsius, mon_molar)


def _eval_pair_cached(
    p_a: PhysicalPrimer,
    p_b: PhysicalPrimer,
    backend: Any,
    enzyme: Any,
    temp_celsius: float,
    backend_kwargs: Dict[str, Any],
    pair_cache: Optional[Dict[Tuple[Any, ...], Tuple[bool, float, Dict[str, Any]]]],
) -> Tuple[bool, float, Dict[str, Any]]:
    """Évalue `evaluate_pair_amplifiable` avec cache indexé par séquences et modifications physiques.
    Evaluates `evaluate_pair_amplifiable` with a cache keyed by sequences and physical modifications.
    """
    pa_parent = p_a.parent_name or p_a.name.split("#")[0]
    pb_parent = p_b.parent_name or p_b.name.split("#")[0]
    key = (
        p_a.sequence, p_b.sequence,
        p_a.lna_positions, p_b.lna_positions,
        p_a.blocked_3prime, p_b.blocked_3prime,
    )
    if pair_cache is not None and key in pair_cache:
        is_amp, dg_3p, base_det = pair_cache[key]
    else:
        is_amp, dg_3p, raw_det = evaluate_pair_amplifiable(
            p_a, p_b, backend, enzyme, temp_celsius, **backend_kwargs
        )
        base_det = {
            "seq_a": raw_det["seq_a"],
            "seq_b": raw_det["seq_b"],
            "structure": raw_det["structure"],
            "delta_g": raw_det["delta_g"],
            "delta_g_3p": raw_det["delta_g_3p"],
            "order": raw_det["order"],
            "extensible_strand": raw_det["extensible_strand"],
            "is_blocked_veto": raw_det["is_blocked_veto"],
        }
        if pair_cache is not None:
            pair_cache[key] = (is_amp, dg_3p, base_det)

    order_str = base_det["order"]
    ext_strand = base_det["extensible_strand"]
    if order_str == "a&a":
        ext_par, tpl_par = (pa_parent, pa_parent) if ext_strand else ("", "")
    else:
        s_a = pa_parent if order_str == "a&b" else pb_parent
        s_b = pb_parent if order_str == "a&b" else pa_parent
        if ext_strand == "a":
            ext_par, tpl_par = s_a, s_b
        elif ext_strand == "b":
            ext_par, tpl_par = s_b, s_a
        else:
            ext_par, tpl_par = "", ""

    det = dict(base_det)
    det["parent_a"] = pa_parent
    det["parent_b"] = pb_parent
    det["extensible_primer_parent"] = ext_par
    det["template_primer_parent"] = tpl_par
    return is_amp, dg_3p, det


def prescreen_dimers(
    cand: Candidate,
    cand_name: str,
    role: PrimerRole,
    total_conc_molar: float,
    partners: List[PhysicalPrimer],
    backend: Any,
    enzyme: Any,
    temp_celsius: float,
    backend_kwargs: Dict[str, Any],
    pair_cache: Optional[Dict[Tuple[Any, ...], Tuple[bool, float, Dict[str, Any]]]] = None,
) -> None:
    """Évalue les dimères amplifiables du candidat face aux autres amorces (critère 3' de `analyze`).
    Evaluates the candidate's amplifiable dimers against other primers (same 3' criterion as `analyze`).

    Score de criblage : somme de c_partenaire × exp(−ΔG/RT) sur les paires amplifiables,
    proportionnelle à la concentration de dimère en régime dilué.
    Prescreen score: sum of c_partner × exp(−ΔG/RT) over amplifiable pairs,
    proportional to dimer concentration in the dilute regime.
    """
    rt = R_KCAL * (temp_celsius + 273.15)
    raw_seq = cand.raw_sequence or cand.sequence
    variants = expand_degenerate(raw_seq)
    v_conc = total_conc_molar / len(variants)
    cand_primers = []
    for i, v_raw in enumerate(variants):
        v_bare, lna_pos = parse_lna_sequence(v_raw)
        v_name = f"{cand_name}#{i+1}" if len(variants) > 1 else cand_name
        cand_primers.append(
            PhysicalPrimer(
                v_name, v_bare, role, v_bare,
                nominal_concentration=v_conc, parent_name=cand_name,
                parent_binding_domain=cand.sequence,
                lna_positions=tuple(lna_pos),
            )
        )
    score = 0.0
    n_amp = 0
    partner_set: Dict[str, None] = {}
    pairs: List[Tuple[PhysicalPrimer, PhysicalPrimer]] = []
    for cp in cand_primers:
        for other in partners:
            pairs.append((cp, other))
    for i, cp in enumerate(cand_primers):
        for cp2 in cand_primers[i:]:
            pairs.append((cp, cp2))
    for p_a, p_b in pairs:
        is_amp, _dg3, details = _eval_pair_cached(
            p_a, p_b, backend, enzyme, temp_celsius, backend_kwargs, pair_cache
        )
        if not is_amp:
            continue
        n_amp += 1
        dg = float(details.get("delta_g", 0.0) or 0.0)
        c_b = float(p_b.nominal_concentration or 0.0)
        score += c_b * math.exp(min(-dg / rt, 700.0))
        partner_set[(p_b.parent_name or p_b.name.split("#")[0])] = None
    cand.n_amplifiable_pairs = n_amp
    cand.prescreen_score = score * v_conc
    cand.amplifiable_partners = list(partner_set.keys())


class IncrementalRedesignEvaluator:
    """Évaluateur incrémental en mémoire de l'équilibre couplé pour `labcraft redesign`.
    In-memory incremental coupled-equilibrium evaluator for `labcraft redesign`.

    Pré-calcule une seule fois le sous-réseau invariant des partenaires (dimères partenaire×partenaire,
    sites cibles des partenaires, ouverture génomique et alertes de mésamorçage des partenaires),
    puis ne met à jour que les espèces impliquant le candidat lors de chaque évaluation.
    Precomputes the invariant partner subnetwork once (partner×partner dimers, partner target sites,
    genomic unfolding, and partner mispriming alerts), then updates only candidate-related species.
    """

    def __init__(
        self,
        *,
        engine_primers: List[PhysicalPrimer],
        internal_name: str,
        target_id: str,
        role_enum: PrimerRole,
        total_conc_molar: float,
        targets_dict: Dict[str, str],
        primer_to_target: Dict[str, str],
        primer_to_panel_map: Dict[str, str],
        profiles: Dict[str, ConcentrationProfile],
        buffer_dump: Optional[Dict[str, Any]],
        backend: Any,
        backend_kwargs: Dict[str, Any],
        mon_molar: Optional[float],
        enzyme: Any,
        temp_celsius: float,
        max_variants: int = 64,
        errors: int = 2,
        strict_3prime: int = 3,
        strict_3prime_tolerate: int = 2,
        dimer_conc_max: Optional[float] = None,
        dimer_warn_ratio: float = DEFAULT_DIMER_WARN_RATIO,
        role_imbalance_ratio: float = DEFAULT_ROLE_IMBALANCE_RATIO,
        pair_cache: Optional[Dict[Tuple[Any, ...], Tuple[bool, float, Dict[str, Any]]]] = None,
    ) -> None:
        self.internal_name = internal_name
        self.target_id = target_id
        self.role_enum = role_enum
        self.total_conc_molar = total_conc_molar
        self.targets_dict = targets_dict
        self.is_multi_target = len(targets_dict) > 1
        self.base_primer_to_target = dict(primer_to_target)
        self.base_primer_to_panel_map = dict(primer_to_panel_map)
        self.profiles = profiles
        self.buffer_dump = buffer_dump
        self.backend = backend
        self.backend_kwargs = backend_kwargs
        self.mon_molar = mon_molar
        self.enzyme = enzyme
        self.temp_celsius = temp_celsius
        self.max_variants = max_variants
        self.errors = errors
        self.strict_3prime = strict_3prime
        self.strict_3prime_tolerate = strict_3prime_tolerate
        self.dimer_conc_max = dimer_conc_max
        self.dimer_warn_ratio = dimer_warn_ratio
        self.role_imbalance_ratio = role_imbalance_ratio
        self.pair_cache = pair_cache if pair_cache is not None else {}

        # Conserver la position exacte de l'amorce remplacée dans la liste des brins
        # Preserve the exact position of the replaced primer within the strand list
        self.primers_before: List[PhysicalPrimer] = []
        self.primers_after: List[PhysicalPrimer] = []
        self.orig_blocked = False
        self.orig_label: Optional[str] = None
        seen_target_primer = False
        for p in engine_primers:
            par = p.parent_name or p.name.split("#")[0]
            if par == internal_name:
                seen_target_primer = True
                self.orig_blocked = p.blocked_3prime
                self.orig_label = p.label_5prime
            elif not seen_target_primer:
                self.primers_before.append(p)
            else:
                self.primers_after.append(p)
        self.partners = self.primers_before + self.primers_after

        # 1. Pré-calculer une seule fois les sites cibles et complexes `_on_` des partenaires
        # 1. Precompute partner target sites and `_on_` complexes once
        partner_target_arg: Any = (
            targets_dict if self.is_multi_target else next(iter(targets_dict.values()))
        )
        partner_prof_arg: Any = (
            profiles if self.is_multi_target else profiles.get(target_id, next(iter(profiles.values()), LAMP_DEFAULT_PROFILE))
        )
        p_prob, p_strands, p_complexes, p_unfold = enumerate_complexes(
            self.partners,
            partner_target_arg,
            backend,
            profile=partner_prof_arg,
            temp_celsius=temp_celsius,
            mon_molar=mon_molar,
            buffer=buffer_dump,
            enzyme=enzyme,
            max_errors=errors,
            strict_3prime_len=strict_3prime,
            strict_3prime_tolerate=strict_3prime_tolerate,
            primer_to_target=primer_to_target,
        )
        self.partner_unfolding = p_unfold
        n_part = len(self.partners)
        self.partner_sites: List[Tuple[str, float, List[Tuple[str, float]]]] = []
        for s_rel_idx, s_name in enumerate(p_strands[n_part:]):
            s_col = n_part + s_rel_idx
            s_conc = float(p_prob.total_concentrations[s_col])
            on_bindings: List[Tuple[str, float]] = []
            for c_idx, c_name in enumerate(p_complexes):
                if "_on_" in c_name and p_prob.stoichiometry[c_idx, s_col] > 0:
                    for p_i in range(n_part):
                        if p_prob.stoichiometry[c_idx, p_i] > 0:
                            on_bindings.append((self.partners[p_i].name, float(p_prob.delta_g[c_idx])))
                            break
            self.partner_sites.append((s_name, s_conc, on_bindings))

        # Déterminer l'indice d'insertion du site candidat parmi les sites cibles (ordre de enumerate_complexes)
        # Determine candidate site insertion index among target sites (matching enumerate_complexes order)
        self.cand_site_insert_idx = len(self.partner_sites)
        if self.is_multi_target:
            # Dans enumerate_complexes, les sites sont groupés par cible dans l'ordre de targets_dict,
            # et au sein de target_id selon l'ordre d'apparition des parents dans primers.
            # In enumerate_complexes, sites are grouped by target in targets_dict order,
            # and within target_id by parent appearance order in primers.
            after_parents_on_target = {
                (p.parent_name or p.name.split("#")[0])
                for p in self.primers_after
                if primer_to_target.get(p.parent_name or p.name.split("#")[0]) == target_id
            }
            target_order = list(targets_dict.keys())
            my_t_pos = target_order.index(target_id)
            later_targets = set(target_order[my_t_pos + 1:])
            for idx_s, (s_name, _sc, _ob) in enumerate(self.partner_sites):
                if "@" in s_name:
                    par_s, rest_s = s_name.split("@", 1)
                    t_s = rest_s.removesuffix("_site")
                    if t_s in later_targets or (t_s == target_id and par_s in after_parents_on_target):
                        self.cand_site_insert_idx = idx_s
                        break
        else:
            after_parents = {(p.parent_name or p.name.split("#")[0]) for p in self.primers_after}
            for idx_s, (s_name, _sc, _ob) in enumerate(self.partner_sites):
                par_s = s_name.removesuffix("_site")
                if par_s in after_parents:
                    self.cand_site_insert_idx = idx_s
                    break

        # 2. Pré-calculer une seule fois les énergies d'épingles et des dimères partenaire×partenaire
        # 2. Precompute partner hairpin energies and partner×partner dimer energies once
        self.partner_hairpin_dg: Dict[str, float] = {}
        for p1 in self.partners:
            try:
                res_hp = backend.calc_hairpin(
                    p1.sequence,
                    temp_celsius=temp_celsius,
                    lna_positions=p1.lna_positions,
                    mon_molar=mon_molar,
                    **backend_kwargs,
                )
                dg_hp = float(res_hp.dg_kcal)
            except NotImplementedError:
                dg_hp = 1.0
            except Exception:
                dg_hp = 1.0
            if dg_hp < 0.0:
                self.partner_hairpin_dg[p1.name] = dg_hp

        self.partner_dimer_dg: Dict[Tuple[str, str], float] = {}
        for i, p1 in enumerate(self.partners):
            for j in range(i, len(self.partners)):
                p2 = self.partners[j]
                try:
                    if i == j:
                        res_d = backend.calc_homodimer(
                            p1.sequence, temp_celsius=temp_celsius,
                            lna_positions=p1.lna_positions, **backend_kwargs,
                        )
                    else:
                        res_d = backend.calc_heterodimer(
                            p1.sequence, p2.sequence, temp_celsius=temp_celsius,
                            lna_positions_a=p1.lna_positions, lna_positions_b=p2.lna_positions,
                            **backend_kwargs,
                        )
                    dg_val = float(res_d.dg_kcal)
                except Exception:
                    dg_val = 1.0
                if dg_val < 0.0:
                    self.partner_dimer_dg[(p1.name, p2.name)] = dg_val
                    _eval_pair_cached(
                        p1, p2, backend, enzyme, temp_celsius, backend_kwargs, self.pair_cache
                    )

        # 3. Pré-calculer une seule fois les alertes de mésamorçage des partenaires fixes
        # 3. Precompute mispriming alerts of fixed partners once
        partner_loci_for_mispriming: Dict[str, Any] = {}
        if self.is_multi_target and isinstance(p_unfold, dict):
            partner_loci_for_mispriming = p_unfold
        else:
            for t_k in targets_dict:
                partner_loci_for_mispriming[t_k] = p_unfold
        self.partner_mispriming_count = len(
            detect_inter_target_mispriming(
                self.partners,
                primer_to_target,
                targets_dict,
                backend,
                enzyme,
                temp_celsius,
                primer_to_target=primer_to_target,
                legitimate_loci=partner_loci_for_mispriming,
                **backend_kwargs,
            )
        )

    def evaluate_candidate(self, cand: Candidate) -> Dict[str, Any]:
        """Résout l'équilibre couplé complet d'un candidat en mémoire et retourne ses métriques.
        Solves full coupled equilibrium of a candidate in memory and returns its metrics.
        """
        raw_seq = cand.raw_sequence or cand.sequence
        bare_parent_seq, _ = parse_lna_sequence(raw_seq)
        variants = expand_degenerate(raw_seq, max_variants=self.max_variants)
        v_conc = self.total_conc_molar / len(variants)

        p2t = dict(self.base_primer_to_target)
        p2pan = dict(self.base_primer_to_panel_map)
        panel_id = p2pan.get(self.internal_name, self.target_id)
        p2t[self.internal_name] = self.target_id
        p2pan[self.internal_name] = panel_id

        cand_primers: List[PhysicalPrimer] = []
        for idx_v, v_raw in enumerate(variants):
            v_bare, lna_pos = parse_lna_sequence(v_raw)
            v_name = f"{self.internal_name}#{idx_v + 1}" if len(variants) > 1 else self.internal_name
            p2t[v_name] = self.target_id
            p2pan[v_name] = panel_id
            cand_primers.append(
                PhysicalPrimer(
                    v_name, v_bare, self.role_enum, v_bare,
                    nominal_concentration=v_conc, parent_name=self.internal_name,
                    parent_binding_domain=bare_parent_seq,
                    blocked_3prime=self.orig_blocked, label_5prime=self.orig_label,
                    lna_positions=tuple(lna_pos),
                )
            )

        all_primers = self.primers_before + cand_primers + self.primers_after
        cand_name_set = {cp.name for cp in cand_primers}

        # Évaluer uniquement le locus cible du candidat sur sa cible homologue
        # Evaluate only the candidate's target locus on its homologous target
        t_prof = self.profiles.get(self.target_id, next(iter(self.profiles.values()), LAMP_DEFAULT_PROFILE))
        c_prob, c_strands, c_complexes, c_unfold = enumerate_complexes(
            cand_primers,
            self.targets_dict[self.target_id],
            self.backend,
            profile=t_prof,
            temp_celsius=self.temp_celsius,
            mon_molar=self.mon_molar,
            buffer=self.buffer_dump,
            enzyme=self.enzyme,
            max_errors=self.errors,
            strict_3prime_len=self.strict_3prime,
            strict_3prime_tolerate=self.strict_3prime_tolerate,
            primer_to_target=p2t,
        )

        all_sites = list(self.partner_sites)
        all_unfolding: Dict[str, Any] = {}
        if self.is_multi_target:
            for t_k in self.targets_dict:
                all_unfolding[t_k] = dict(self.partner_unfolding.get(t_k, {})) if isinstance(self.partner_unfolding, dict) else {}
        else:
            for t_k in self.targets_dict:
                all_unfolding[t_k] = copy.deepcopy(self.partner_unfolding)

        n_cp = len(cand_primers)
        if len(c_strands) > n_cp:
            raw_s_name = c_strands[n_cp]
            eff_s_name = f"{self.internal_name}@{self.target_id}_site" if self.is_multi_target else raw_s_name
            s_conc = float(c_prob.total_concentrations[n_cp])
            c_on_bindings: List[Tuple[str, float]] = []
            for c_idx, c_name in enumerate(c_complexes):
                if "_on_" in c_name and c_prob.stoichiometry[c_idx, n_cp] > 0:
                    for p_i in range(n_cp):
                        if c_prob.stoichiometry[c_idx, p_i] > 0:
                            c_on_bindings.append((cand_primers[p_i].name, float(c_prob.delta_g[c_idx])))
                            break
            all_sites.insert(self.cand_site_insert_idx, (eff_s_name, s_conc, c_on_bindings))
            if isinstance(c_unfold, dict) and raw_s_name in c_unfold:
                all_unfolding[self.target_id][eff_s_name] = c_unfold[raw_s_name]

        # Construire le problème d'équilibre complet (brins, monomères, dimères, complexes sur cible)
        # Assemble the full EquilibriumProblem (strands, monomers, dimers, primer-target complexes)
        n_primers = len(all_primers)
        n_sites = len(all_sites)
        n_strands = n_primers + n_sites
        strand_names = [p.name for p in all_primers] + [s[0] for s in all_sites]
        primer_idx_map = {p.name: idx for idx, p in enumerate(all_primers)}

        concentrations = np.zeros(n_strands, dtype=np.float64)
        for i, p in enumerate(all_primers):
            concentrations[i] = float(p.nominal_concentration or 0.0)
        for i, (_sn, sc, _ob) in enumerate(all_sites):
            concentrations[n_primers + i] = sc

        complexes: List[ComplexInfo] = []
        dimer_pairs_meta: List[Optional[Tuple[PhysicalPrimer, PhysicalPrimer]]] = []

        for i, p in enumerate(all_primers):
            stoich = [0] * n_strands
            stoich[i] = 1
            complexes.append(ComplexInfo(f"{p.name}_free", stoich, 0.0))
            dimer_pairs_meta.append(None)

        for i, (s_name, _sc, _ob) in enumerate(all_sites):
            stoich = [0] * n_strands
            stoich[n_primers + i] = 1
            complexes.append(ComplexInfo(f"{s_name}_free", stoich, 0.0))
            dimer_pairs_meta.append(None)

        for i, p1 in enumerate(all_primers):
            if p1.name not in cand_name_set:
                dg_hp = self.partner_hairpin_dg.get(p1.name, 1.0)
            else:
                try:
                    res_hp = self.backend.calc_hairpin(
                        p1.sequence,
                        temp_celsius=self.temp_celsius,
                        lna_positions=p1.lna_positions,
                        mon_molar=self.mon_molar,
                        **self.backend_kwargs,
                    )
                    dg_hp = float(res_hp.dg_kcal)
                except NotImplementedError:
                    dg_hp = 1.0
                except Exception:
                    dg_hp = 1.0
            if dg_hp < 0.0:
                stoich = [0] * n_strands
                stoich[i] = 1
                complexes.append(ComplexInfo(f"{p1.name}_hairpin", stoich, dg_hp))
                dimer_pairs_meta.append(None)

        for i, p1 in enumerate(all_primers):
            p1_is_cand = p1.name in cand_name_set
            for j in range(i, n_primers):
                p2 = all_primers[j]
                if not p1_is_cand and p2.name not in cand_name_set:
                    dg_d = self.partner_dimer_dg.get((p1.name, p2.name), 1.0)
                else:
                    try:
                        if i == j:
                            res_d = self.backend.calc_homodimer(
                                p1.sequence, temp_celsius=self.temp_celsius,
                                lna_positions=p1.lna_positions, **self.backend_kwargs,
                            )
                        else:
                            res_d = self.backend.calc_heterodimer(
                                p1.sequence, p2.sequence, temp_celsius=self.temp_celsius,
                                lna_positions_a=p1.lna_positions, lna_positions_b=p2.lna_positions,
                                **self.backend_kwargs,
                            )
                        dg_d = float(res_d.dg_kcal)
                    except Exception:
                        dg_d = 1.0
                if dg_d < 0.0:
                    stoich = [0] * n_strands
                    stoich[i] += 1
                    stoich[j] += 1
                    cname = f"{p1.name}_{p2.name}" if i != j else f"{p1.name}_homo"
                    complexes.append(ComplexInfo(cname, stoich, dg_d))
                    dimer_pairs_meta.append((p1, p2))

        for s_rel_i, (s_name, _sc, on_bindings) in enumerate(all_sites):
            site_col = n_primers + s_rel_i
            for p_name, dg_eff in on_bindings:
                p_col = primer_idx_map[p_name]
                stoich = [0] * n_strands
                stoich[p_col] = 1
                stoich[site_col] = 1
                complexes.append(ComplexInfo(f"{p_name}_on_{s_name}", stoich, dg_eff))
                dimer_pairs_meta.append(None)

        stoich_matrix = np.array([c.stoichiometry for c in complexes], dtype=np.float64)
        dg_vector = np.array([c.delta_g for c in complexes], dtype=np.float64)
        complex_names = [c.name for c in complexes]
        prob = EquilibriumProblem(
            n_strands=n_strands,
            n_complexes=len(complexes),
            stoichiometry=stoich_matrix,
            delta_g=dg_vector,
            total_concentrations=concentrations,
            temperature_kelvin=273.15 + self.temp_celsius,
        )

        res = solve_dual(prob)

        amplifiable_flags: List[bool] = []
        concs = [float(x) for x in res.concentrations]
        dimer_details: List[Dict[str, Any]] = []
        for pair_meta in dimer_pairs_meta:
            if pair_meta is None:
                amplifiable_flags.append(False)
                dimer_details.append({})
            else:
                pa, pb = pair_meta
                is_amp, _dg3, det = _eval_pair_cached(
                    pa, pb, self.backend, self.enzyme, self.temp_celsius, self.backend_kwargs, self.pair_cache
                )
                pa_par = det["parent_a"]
                pb_par = det["parent_b"]
                det["panel_a"] = p2pan.get(pa_par, p2pan.get(pa.name, ""))
                det["panel_b"] = p2pan.get(pb_par, p2pan.get(pb.name, ""))
                amplifiable_flags.append(is_amp)
                dimer_details.append(det)

        risks = evaluate_risks(
            complex_names, concs, amplifiable_flags,
            is_warm_start=self.enzyme.warm_start, dimer_details=dimer_details,
        )
        fractions = compute_fractions(
            strand_names, complex_names, stoich_matrix,
            res.free_concentrations, dg_vector, self.temp_celsius, p2t,
        )

        target_occupations: Dict[str, Dict[str, float]] = {t_k: {} for t_k in self.targets_dict}
        site_occupations_json: Dict[str, Dict[str, Optional[float]]] = {t_k: {} for t_k in self.targets_dict}
        for s_rel_i, (s_name, _sc, _ob) in enumerate(all_sites):
            idx_col = n_primers + s_rel_i
            tot_s = float(prob.total_concentrations[idx_col])
            occ = (tot_s - float(res.free_concentrations[idx_col])) / tot_s if tot_s > 0.0 else 0.0
            if "@" in s_name:
                p_part, rest_part = s_name.split("@", 1)
                t_part = rest_part.removesuffix("_site")
                if t_part in target_occupations:
                    target_occupations[t_part][s_name] = occ
                    target_occupations[t_part][f"{p_part}_site"] = occ
            else:
                for t_k in self.targets_dict:
                    target_occupations[t_k][s_name] = occ

        for t_k, s_map in all_unfolding.items():
            if not isinstance(s_map, dict):
                continue
            occ_map = target_occupations.get(t_k, {})
            for s_key, pen in s_map.items():
                if not isinstance(pen, dict):
                    continue
                occ_v = occ_map.get(s_key)
                if occ_v is None and "@" in s_key:
                    occ_v = occ_map.get(f"{s_key.split('@', 1)[0]}_site")
                site_occupations_json[t_k][s_key] = occ_v

        loop_parents = {
            (p.parent_name or p.name.split("#")[0])
            for p in all_primers
            if p.role in (PrimerRole.LF, PrimerRole.LB)
        }
        parent_nominal_concs: Dict[str, float] = {}
        parent_variant_counts: Dict[str, int] = {}
        for p in all_primers:
            par_k = p.parent_name or p.name.split("#")[0]
            parent_nominal_concs[par_k] = parent_nominal_concs.get(par_k, 0.0) + float(
                p.nominal_concentration or 0.0
            )
            parent_variant_counts[par_k] = parent_variant_counts.get(par_k, 0) + 1
        verdict = generate_verdict(
            fractions, target_occupations, risks,
            loop_primer_parents=loop_parents,
            primer_to_panel=p2pan, primer_to_target=p2t,
            unfolding_penalties=all_unfolding,
            dimer_conc_max_nm=self.dimer_conc_max,
            role_imbalance_ratio=self.role_imbalance_ratio,
            dimer_warn_ratio=self.dimer_warn_ratio,
            primer_nominal_concs=parent_nominal_concs,
            primer_variant_counts=parent_variant_counts,
        )
        decomp = summarize_amplifiable_dimers(risks, primer_to_panel=p2pan)
        free_frac_map = {p: f.free for p, f in fractions.items()}
        panel_summaries, _ = calculate_multiplex_balance(
            p2t, target_occupations, free_frac_map,
            loop_primer_parents=loop_parents, primer_to_target=p2t,
        )

        cand_mispriming = detect_inter_target_mispriming(
            cand_primers, p2t, self.targets_dict,
            self.backend, self.enzyme, self.temp_celsius,
            primer_to_target=p2t, legitimate_loci=all_unfolding,
            **self.backend_kwargs,
        )
        fam_labels = [f.get("label", "") for f in decomp.get("families", [])]

        return {
            "amplifiable_total_nm": float(decomp["total_nm"]),
            "amplifiable_inter_nm": float(decomp["inter_nm"]),
            "amplifiable_intra_by_panel_nm": dict(decomp["intra_by_panel_nm"]),
            "dominant_family": decomp.get("dominant_family"),
            "families": fam_labels,
            "initiation_access_by_panel": panel_summaries,
            "site_occupations": site_occupations_json,
            "mispriming_alerts": self.partner_mispriming_count + len(cand_mispriming),
            "verdict": verdict.status,
            "verdict_cause": verdict.global_cause,
        }


def _tiebreak_tail(cand: Candidate, orig_len: int) -> Tuple[int, int, int, int, float, int, int, float]:
    """Critères de départage appliqués après les paires extensibles.
    Tie-breaking criteria applied after extensible pairs:
    1) candidat non plus stable que sa référence interne (F2/B2/F1c/B1c),
    2) plus longue complémentarité 3' parfaite avec matrice 5' (minimisée),
    3) proximité de la longueur d'origine,
    4) GC entre 40 et 60 %,
    5) complémentarité 3' inverse (partenaire -> candidat) et homodimérique,
    6) décalage minimal |shift| puis énergie effective sur cible.
    """
    gc_outside = 0 if 40.0 <= cand.gc_percent <= 60.0 else 1
    gc_dist = round(max(0.0, 40.0 - cand.gc_percent, cand.gc_percent - 60.0), 6)
    exceeds_inner = (
        cand.more_stable_than_ref_inner if cand.dg_ref_inner is not None else cand.more_stable_than_f2
    )
    return (
        1 if exceeds_inner else 0,
        cand.max_3p_comp,
        abs(cand.length - orig_len),
        gc_outside,
        gc_dist,
        max(cand.max_reverse_3p_comp, cand.max_self_3p_comp),
        abs(cand.shift),
        round(cand.dg_duplex + cand.dg_open, 6),
    )


def prescreen_sort_key(cand: Candidate, orig_len: int) -> Tuple[Any, ...]:
    """Clé de tri de présélection / Prescreen ranking key."""
    return (
        1 if cand.n_amplifiable_pairs > 0 else 0,
        round(cand.prescreen_score, 15),
        *_tiebreak_tail(cand, orig_len),
    )


def equilibrium_sort_key(cand: Candidate, orig_len: int) -> Tuple[Any, ...]:
    """Clé de tri après équilibre couplé complet / Full coupled-equilibrium ranking key.

    Les concentrations (affichées au dix-millième de nM, ex. 0,0586 nM) sont comparées à 1e-4 nM
    afin que les ex aequo à l'équilibre soient départagés par les règles de conception.
    Concentrations (displayed to 1e-4 nM, e.g. 0.0586 nM) are compared at 1e-4 nM precision
    so equilibrium ties are broken by the design tie-breaking rules.
    """
    eq = cand.equilibrium
    if not eq:
        return (1e9, 1e9, 1, *(1, 999, 999, 1, 100.0, 999, 999, 1e9))
    return (
        round(float(eq["amplifiable_total_nm"]), 4),
        round(float(eq["amplifiable_inter_nm"]), 4),
        1 if cand.n_amplifiable_pairs > 0 else 0,
        *_tiebreak_tail(cand, orig_len),
    )


def select_top_candidates(valid_sorted: List[Candidate], top: int, role: str) -> Tuple[List[Candidate], int]:
    """Sélectionne les `top` meilleurs candidats en écartant ceux plus stables que leur référence interne.
    Selects `top` best candidates, discarding candidates more stable than their inner reference when possible.

    Retourne (liste_retenue, nombre_ecartes_plus_stables_que_reference_interne).
    Returns (selected_list, n_discarded_more_stable_than_inner_ref).
    """
    k = max(0, top)
    if k == 0 or not valid_sorted:
        return [], 0
    if role.upper() not in REDESIGNABLE_ROLES:
        return valid_sorted[:k], 0

    def _is_too_stable(c: Candidate) -> bool:
        if c.dg_ref_inner is not None:
            return c.more_stable_than_ref_inner
        return c.more_stable_than_f2

    n_more_stable = sum(1 for c in valid_sorted if _is_too_stable(c))
    has_zero_amp_compliant = any(c.n_amplifiable_pairs == 0 and not _is_too_stable(c) for c in valid_sorted)
    has_any_zero_amp = any(c.n_amplifiable_pairs == 0 for c in valid_sorted)

    if has_zero_amp_compliant or (not has_any_zero_amp and any(not _is_too_stable(c) for c in valid_sorted)):
        filtered = [c for c in valid_sorted if not _is_too_stable(c)]
        return filtered[:k], n_more_stable
    # Repli : tous les meilleurs candidats sans dimères dépassent la référence -> on les conserve en les signalant
    # Fallback: all zero-dimer candidates exceed reference stability -> keep and flag them
    return valid_sorted[:k], 0


def rewrite_primer_file(records: List[Tuple[str, str]], raw_name: str, new_seq: str, path: str) -> None:
    """Réécrit le fichier d'amorces en remplaçant la séquence de l'amorce visée.
    Rewrites the primer file, replacing the target primer sequence.
    """
    with open(path, "w", encoding="utf-8") as f:
        for header, seq in records:
            name = header.split()[0]
            f.write(f">{header}\n{new_seq if name == raw_name else seq}\n")
