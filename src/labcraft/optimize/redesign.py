"""Redessin d'une amorce externe ou de boucle LAMP (commande `labcraft redesign`).
Redesign of an outer or loop LAMP primer (`labcraft redesign` command).

Périmètre : F3, B3, LF, LB. Les amorces internes FIP/BIP sont refusées.
Scope: F3, B3, LF, LB. Inner primers FIP/BIP are refused.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

from labcraft.diagnostics.amplifiable_dimer import evaluate_pair_amplifiable
from labcraft.lamp.domains import (
    PhysicalPrimer,
    PrimerRole,
    _find_iupac_substring,
    expand_degenerate,
)
from labcraft.target.seeding import primer_matches_sequence

REDESIGNABLE_ROLES = ("F3", "B3", "LF", "LB")
REFUSED_INNER_ROLES = ("FIP", "BIP", "F1", "F2", "B1", "B2")
UNFOLDING_WINDOW = 150  # identique à enumerate_complexes / same as enumerate_complexes
R_KCAL = 0.0019872


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
    is_reference: bool = False
    geometry_ok: bool = True
    dg_duplex: float = 0.0
    dg_open: float = 0.0
    n_amplifiable_pairs: int = 0
    prescreen_score: float = 0.0
    amplifiable_partners: List[str] = field(default_factory=list)
    equilibrium: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        """Sérialisation JSON / JSON serialization."""
        return {
            "shift": self.shift,
            "length": self.length,
            "sequence": self.sequence,
            "start": self.locus.start,
            "end": self.locus.end,
            "strand": self.locus.strand,
            "is_reference": self.is_reference,
            "geometry_ok": self.geometry_ok,
            "dg_duplex": self.dg_duplex,
            "dg_open": self.dg_open,
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


def locate_on_target(seq: str, target: str) -> Optional[Locus]:
    """Localise une séquence (amorce ou domaine) sur la cible, recherche exacte IUPAC puis tolérante.
    Locates a sequence (primer or domain) on the target, exact IUPAC search then tolerant.

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
        geo["F1"] = locate_on_target(fip.tail_domain or "", target)
    if bip is not None:
        geo["B2"] = locate_on_target(bip.binding_domain, target)
        geo["B1"] = locate_on_target(bip.tail_domain or "", target)
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
    out: List[Candidate] = [
        Candidate(0, orig_len, original_seq.upper(), original_locus, is_reference=True)
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
            out.append(Candidate(d, L, seq, Locus(start, end, original_locus.strand)))
    return out


def compute_site_energies(
    cand: Candidate,
    role: str,
    target: str,
    backend,
    temp_celsius: float,
    mon_molar: Optional[float],
    backend_kwargs: Dict[str, Any],
) -> None:
    """ΔG du duplex parfait et ΔG d'ouverture de la cible (même fenêtre que `analyze`).
    Perfect-duplex ΔG and target opening ΔG (same window as `analyze`).
    """
    from labcraft.target.unfolding import calc_unfolding_penalty

    dgs = []
    for v in expand_degenerate(cand.sequence):
        res = backend.calc_duplex(v, revcomp(v), temp_celsius=temp_celsius, **backend_kwargs)
        dgs.append(res.dg_kcal)
    cand.dg_duplex = float(min(dgs))
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


def prescreen_dimers(
    cand: Candidate,
    cand_name: str,
    role: PrimerRole,
    total_conc_molar: float,
    partners: List[PhysicalPrimer],
    backend,
    enzyme,
    temp_celsius: float,
    backend_kwargs: Dict[str, Any],
) -> None:
    """Évalue les dimères amplifiables du candidat face aux autres amorces (critère 3' de `analyze`).
    Evaluates the candidate's amplifiable dimers against other primers (same 3' criterion as `analyze`).

    Score de criblage : somme de c_partenaire × exp(−ΔG/RT) sur les paires amplifiables,
    proportionnelle à la concentration de dimère en régime dilué.
    Prescreen score: sum of c_partner × exp(−ΔG/RT) over amplifiable pairs,
    proportional to dimer concentration in the dilute regime.
    """
    rt = R_KCAL * (temp_celsius + 273.15)
    variants = expand_degenerate(cand.sequence)
    v_conc = total_conc_molar / len(variants)
    cand_primers = [
        PhysicalPrimer.from_simple(
            f"{cand_name}#{i+1}" if len(variants) > 1 else cand_name, v, role,
            nominal_concentration=v_conc, parent_name=cand_name,
        )
        for i, v in enumerate(variants)
    ]
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
        is_amp, _dg3, details = evaluate_pair_amplifiable(
            p_a, p_b, backend, enzyme, temp_celsius, **backend_kwargs
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


def rewrite_primer_file(records: List[Tuple[str, str]], raw_name: str, new_seq: str, path: str) -> None:
    """Réécrit le fichier d'amorces en remplaçant la séquence de l'amorce visée.
    Rewrites the primer file, replacing the target primer sequence.
    """
    with open(path, "w", encoding="utf-8") as f:
        for header, seq in records:
            name = header.split()[0]
            f.write(f">{header}\n{new_seq if name == raw_name else seq}\n")
