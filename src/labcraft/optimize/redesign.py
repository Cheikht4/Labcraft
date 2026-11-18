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
    is_reference: bool = False
    geometry_ok: bool = True
    gc_percent: float = 0.0
    dg_duplex: float = 0.0
    dg_f2: Optional[float] = None
    ddg_vs_f2: Optional[float] = None
    more_stable_than_f2: bool = False
    dg_open: float = 0.0
    max_3p_comp: int = 0
    max_3p_partner: str = ""
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
            "dg_open": self.dg_open,
            "max_3p_complementarity": self.max_3p_comp,
            "max_3p_partner": self.max_3p_partner,
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
    orig_up = original_seq.upper()
    out: List[Candidate] = [
        Candidate(0, orig_len, orig_up, original_locus, is_reference=True, gc_percent=gc_percentage(orig_up))
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
                Candidate(d, L, seq, Locus(start, end, original_locus.strand), gc_percent=gc_percentage(seq))
            )
    return out


def compute_panel_f2_dg(
    engine_primers: List[PhysicalPrimer],
    panel_id: str,
    backend,
    temp_celsius: float,
    backend_kwargs: Dict[str, Any],
) -> Optional[float]:
    """Calcule le ΔG du duplex parfait du domaine F2 (liaison de FIP) du même jeu d'amorces.
    Computes the perfect-duplex ΔG of the F2 domain (FIP binding domain) of the same primer set.
    """
    fip_parent = f"FIP_{panel_id}"
    seen: set[Tuple[str, Tuple[int, ...]]] = set()
    dgs: List[float] = []
    for p in engine_primers:
        par = p.parent_name or p.name.split("#")[0]
        if par != fip_parent or not p.binding_domain:
            continue
        bd = p.binding_domain
        offset = p.sequence.find(bd)
        bd_lna = (
            tuple(pos - offset for pos in p.lna_positions if offset <= pos < offset + len(bd))
            if offset >= 0
            else ()
        )
        key = (bd, bd_lna)
        if key in seen:
            continue
        seen.add(key)
        res = backend.calc_duplex(
            bd, revcomp(bd), temp_celsius=temp_celsius,
            lna_positions_a=bd_lna, lna_positions_b=(), **backend_kwargs,
        )
        dgs.append(float(res.dg_kcal))
    return min(dgs) if dgs else None


def compute_f2_comparison(cand: Candidate, dg_f2: Optional[float]) -> None:
    """Compare le ΔG du candidat à celui du domaine F2 du même jeu.
    Compares candidate duplex ΔG against the same panel's F2 domain ΔG.
    """
    cand.dg_f2 = dg_f2
    if dg_f2 is None:
        cand.ddg_vs_f2 = None
        cand.more_stable_than_f2 = False
        return
    cand.ddg_vs_f2 = float(cand.dg_duplex - dg_f2)
    # Plus négatif que F2 => plus stable que F2 (contraire à la règle de conception des externes)
    # More negative than F2 => more stable than F2 (violates LAMP outer-primer design rule)
    cand.more_stable_than_f2 = bool(cand.dg_duplex < dg_f2 - 1e-9)


def compute_max_3p_complementarity(cand: Candidate, partners: List[PhysicalPrimer]) -> None:
    """Calcule la plus longue complémentarité 3' parfaite avec une autre séquence offrant une matrice en 5'.
    Computes longest perfect 3' complementarity with another mixture sequence offering a 5' template.

    Une séquence partenaire P offre une matrice en 5' si le complément inverse des k dernières bases
    en 3' du candidat apparaît dans P[1:] (au moins 1 nucléotide en amont côté 5' de P).
    A partner sequence P offers a 5' template iff revcomp(cand[-k:]) occurs in P[1:] (at least 1 nt
    upstream on the 5' side of P).
    """
    cand.gc_percent = gc_percentage(cand.sequence)
    variants = expand_degenerate(cand.sequence)
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
            # Tester k croissant à partir de best_k + 1 / Probe increasing k starting at best_k + 1
            if revcomp(v[-(best_k + 1):]) not in p_tail_space:
                continue
            for k in range(best_k + 1, max_k + 1):
                if revcomp(v[-k:]) in p_tail_space:
                    best_k = k
                    best_partner = p.parent_name or p.name.split("#")[0]
                else:
                    break
    cand.max_3p_comp = best_k
    cand.max_3p_partner = best_partner


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


def _tiebreak_tail(cand: Candidate, orig_len: int) -> Tuple[int, int, int, int, float, int, float]:
    """Critères de départage appliqués après les paires extensibles.
    Tie-breaking criteria applied after extensible pairs:
    1) candidat non plus stable que F2 / not more stable than F2,
    2) plus longue complémentarité 3' parfaite avec matrice 5' (minimisée) / max 3' match with 5' template,
    3) proximité de la longueur d'origine / proximity to original length,
    4) GC entre 40 et 60 % / GC within [40%, 60%].
    """
    gc_outside = 0 if 40.0 <= cand.gc_percent <= 60.0 else 1
    gc_dist = round(max(0.0, 40.0 - cand.gc_percent, cand.gc_percent - 60.0), 6)
    return (
        1 if cand.more_stable_than_f2 else 0,
        cand.max_3p_comp,
        abs(cand.length - orig_len),
        gc_outside,
        gc_dist,
        abs(cand.shift),
        round(cand.dg_duplex + cand.dg_open, 6),
    )


def prescreen_sort_key(cand: Candidate, orig_len: int) -> Tuple:
    """Clé de tri de présélection / Prescreen ranking key."""
    return (
        1 if cand.n_amplifiable_pairs > 0 else 0,
        round(cand.prescreen_score, 15),
        *_tiebreak_tail(cand, orig_len),
    )


def equilibrium_sort_key(cand: Candidate, orig_len: int) -> Tuple:
    """Clé de tri après équilibre couplé complet / Full coupled-equilibrium ranking key.

    Les concentrations (affichées au dix-millième de nM, ex. 0,0586 nM) sont comparées à 1e-4 nM
    afin que les ex aequo à l'équilibre soient départagés par les règles de conception.
    Concentrations (displayed to 1e-4 nM, e.g. 0.0586 nM) are compared at 1e-4 nM precision
    so equilibrium ties are broken by the design tie-breaking rules.
    """
    eq = cand.equilibrium
    if not eq:
        return (1e9, 1e9, 1, *(1, 999, 999, 1, 100.0, 999, 1e9))
    return (
        round(float(eq["amplifiable_total_nm"]), 4),
        round(float(eq["amplifiable_inter_nm"]), 4),
        1 if cand.n_amplifiable_pairs > 0 else 0,
        *_tiebreak_tail(cand, orig_len),
    )


def select_top_candidates(valid_sorted: List[Candidate], top: int, role: str) -> Tuple[List[Candidate], int]:
    """Sélectionne les `top` meilleurs candidats en écartant ceux plus stables que F2 quand possible.
    Selects `top` best candidates, discarding candidates more stable than F2 when compliant ones exist.

    Retourne (liste_retenue, nombre_ecartes_plus_stables_que_f2).
    Returns (selected_list, n_discarded_more_stable_than_f2).
    """
    k = max(0, top)
    if k == 0 or not valid_sorted:
        return [], 0
    if role.upper() not in ("F3", "B3"):
        return valid_sorted[:k], 0

    n_more_stable = sum(1 for c in valid_sorted if c.more_stable_than_f2)
    has_zero_amp_compliant = any(c.n_amplifiable_pairs == 0 and not c.more_stable_than_f2 for c in valid_sorted)
    has_any_zero_amp = any(c.n_amplifiable_pairs == 0 for c in valid_sorted)

    if has_zero_amp_compliant or (not has_any_zero_amp and any(not c.more_stable_than_f2 for c in valid_sorted)):
        filtered = [c for c in valid_sorted if not c.more_stable_than_f2]
        return filtered[:k], n_more_stable
    # Repli : tous les meilleurs candidats sans dimères dépassent F2 -> on les conserve en les signalant
    # Fallback: all zero-dimer candidates exceed F2 stability -> keep and flag them
    return valid_sorted[:k], 0


def rewrite_primer_file(records: List[Tuple[str, str]], raw_name: str, new_seq: str, path: str) -> None:
    """Réécrit le fichier d'amorces en remplaçant la séquence de l'amorce visée.
    Rewrites the primer file, replacing the target primer sequence.
    """
    with open(path, "w", encoding="utf-8") as f:
        for header, seq in records:
            name = header.split()[0]
            f.write(f">{header}\n{new_seq if name == raw_name else seq}\n")
