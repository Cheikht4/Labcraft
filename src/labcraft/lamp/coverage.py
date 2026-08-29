"""Thermodynamic coverage analysis for multi-strain LAMP panels.
Analyse thermodynamique de couverture multi-souches pour panels LAMP.
"""

from __future__ import annotations

import csv
import warnings
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
from enum import Enum

from labcraft.lamp.domains import PhysicalPrimer, PrimerRole, IUPAC_MATCHABLE
from labcraft.thermo.backends.vienna_salt import ViennaSaltShiftBackend
from labcraft.thermo.mismatch import calculate_hybridization_dg, three_prime_extensible, nn_duplex_energy
from labcraft.diagnostics.enzyme import PolymeraseProfile


class SiteVerdict(Enum):
    PARFAIT = "parfait"
    TOLERABLE = "tolerable"
    VETO_3P = "veto_3p"
    NON_VIABLE = "non_viable"
    ABSENT = "absent"


@dataclass
class SiteEvaluation:
    strain_id: str
    primer_name: str
    primer_role: PrimerRole
    verdict: SiteVerdict
    n_mismatches_count: int
    dg_hyb: Optional[float] = None
    ddg: Optional[float] = None
    first_bad_pos: Optional[int] = None
    severity: Optional[str] = None
    position: Optional[int] = None
    strand: Optional[str] = None


@dataclass
class StrainVerdict:
    strain_id: str
    evaluations: Dict[str, SiteEvaluation]  # key = role name (e.g. "F3", "B3", "FIP", "BIP", "LF", "LB")
    is_amplifiable_thermo: bool
    is_amplifiable_count: bool
    limiting_primer_thermo: Optional[str] = None


def _evaluation_rank(eval_obj: SiteEvaluation) -> Tuple[int, int, float]:
    """Calcule le rang de qualité d'une évaluation (plus élevé = meilleur).
    Computes quality rank of a site evaluation (higher = better).

    Hiérarchie / Hierarchy:
    1. Verdict : PARFAIT (4) > TOLERABLE (3) > NON_VIABLE (2) > VETO_3P (1) > ABSENT (0)
    2. Moins de mésappariements / Fewer mismatches
    3. Pénalité relative ddG minimale (ou dG le plus favorable)
    """
    verdict_scores = {
        SiteVerdict.PARFAIT: 4,
        SiteVerdict.TOLERABLE: 3,
        SiteVerdict.NON_VIABLE: 2,
        SiteVerdict.VETO_3P: 1,
        SiteVerdict.ABSENT: 0,
    }
    v_score = verdict_scores.get(eval_obj.verdict, 0)
    mm_score = -eval_obj.n_mismatches_count
    # ddg plus faible = meilleur score
    score_ddg = -eval_obj.ddg if eval_obj.ddg is not None else (-eval_obj.dg_hyb if eval_obj.dg_hyb is not None else -999.0)
    return (v_score, mm_score, score_ddg)


class CoverageAnalyzer:
    """Analyseur de couverture multi-souches basé sur la thermodynamique et les règles d'élongation.
    Multi-strain coverage analyzer based on thermodynamics and elongation rules.
    """

    def __init__(
        self,
        primers: List[PhysicalPrimer],
        fasta_dict: Dict[str, str],
        backend: ViennaSaltShiftBackend,
        enzyme: PolymeraseProfile,
        temp_celsius: float = 65.0,
        ddg_max: float = 3.0,
        dg_threshold: Optional[float] = None,
        max_mismatches_count: int = 2
    ):
        if not isinstance(backend, ViennaSaltShiftBackend):
            warnings.warn(
                "Aucun modèle de sel/tampon (ViennaSaltShiftBackend) fourni à CoverageAnalyzer. "
                "Repli sur les conditions de référence standard du backend (1 M NaCl / sans correction de sel LAMP). "
                "Les énergies de duplexe peuvent différer de la force ionique réelle.",
                UserWarning
            )

        self.primers = primers
        self.primers_by_name = {p.name: p for p in primers}
        self.fasta_dict = fasta_dict
        self.backend = backend
        self.enzyme = enzyme
        self.temp_celsius = temp_celsius
        self.ddg_max = ddg_max
        self.dg_threshold = dg_threshold
        self.max_mismatches_count = max_mismatches_count

        # Cache pour l'énergie des duplexes parfaits
        # Cache for perfectly matched duplex energy
        self._perfect_dg_cache: Dict[str, float] = {}

    def _get_perfect_dg_seq(self, seq: str) -> float:
        """Calcule et met en cache le dG d'un duplexe parfait pour une séquence donnée.
        Computes and caches the perfect duplex dG for a given sequence.
        """
        if seq not in self._perfect_dg_cache:
            comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            perfect_bottom = "".join(comp.get(c, 'N') for c in seq)
            dg_hyb, _ = calculate_hybridization_dg(
                seq, perfect_bottom, self.temp_celsius, self.backend
            )
            self._perfect_dg_cache[seq] = dg_hyb
        return self._perfect_dg_cache[seq]

    def evaluate_site(
        self,
        primer: PhysicalPrimer,
        site_seq: str,
        strand: str,
        strain_id: str = "",
        position: Optional[int] = None
    ) -> SiteEvaluation:
        """Évalue thermodynamiquement un site d'hybridation avec résolution IUPAC.
        Thermodynamically evaluates a binding site with IUPAC base resolution.

        Convention de brin / Strand convention:
        - '+' : la séquence de l'amorce correspond au brin sens (+) de la cible.
          '+' : primer sequence matches sense (+) strand of target.
        - '-' : le complément inverse de l'amorce apparaît sur le brin sens (+).
          '-' : reverse complement of primer appears on sense (+) strand.
        """
        if len(site_seq) != len(primer.binding_domain):
            raise ValueError(
                f"Longueur de site ({len(site_seq)} nt) incompatible avec le domaine de liaison "
                f"de l'amorce {primer.name} ({len(primer.binding_domain)} nt) pour la souche '{strain_id}'."
            )

        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        if strand == '+':
            bottom_under_top = "".join(comp.get(c, 'N') for c in site_seq)
        else:
            bottom_under_top = site_seq[::-1]

        # Résolution IUPAC et comptage des mésappariements
        # IUPAC resolution and mismatch counting
        resolved_primer_chars = []
        mismatch_count = 0

        for a, b in zip(primer.binding_domain, bottom_under_top):
            a_set = IUPAC_MATCHABLE.get(a.upper())
            b_comp = comp.get(b.upper(), '')
            b_comp_set = IUPAC_MATCHABLE.get(b_comp)

            if a_set is not None and b_comp_set is not None and a_set.intersection(b_comp_set):
                # Appariement valide / Valid base pair
                matched_bases = a_set.intersection(b_comp_set)
                if b_comp in matched_bases:
                    resolved_primer_chars.append(b_comp)
                else:
                    resolved_primer_chars.append(sorted(matched_bases)[0])
            else:
                # Mésappariement / Mismatch
                mismatch_count += 1
                if a.upper() in ('A', 'C', 'G', 'T'):
                    resolved_primer_chars.append(a.upper())
                elif a_set:
                    resolved_primer_chars.append(sorted(a_set)[0])
                else:
                    resolved_primer_chars.append('N')

        resolved_primer_seq = "".join(resolved_primer_chars)
        dg_perfect_backend = self._get_perfect_dg_seq(resolved_primer_seq)

        # Calcul thermodynamique / Thermodynamic calculation
        if mismatch_count == 0:
            dg_hyb = dg_perfect_backend
            ddg = 0.0
        else:
            perfect_bottom = "".join(comp.get(c, c) for c in resolved_primer_seq)
            _, _, dg_perfect_nn = nn_duplex_energy(resolved_primer_seq, perfect_bottom, self.temp_celsius)
            _, _, dg_mismatched_nn = nn_duplex_energy(resolved_primer_seq, bottom_under_top, self.temp_celsius)
            ddg_mismatch = dg_mismatched_nn - dg_perfect_nn
            dg_hyb = dg_perfect_backend + ddg_mismatch
            ddg = max(0.0, dg_hyb - dg_perfect_backend)

        # Règle d'extension enzymatique en 3' / Enzymatic 3' extension rule
        extensible, first_bad_pos, severity = three_prime_extensible(
            resolved_primer_seq, bottom_under_top, self.enzyme
        )

        # Règle de viabilité : critère relatif ddG <= ddG_max et filtre absolu optionnel
        is_viable_energy = (ddg <= self.ddg_max)
        if self.dg_threshold is not None and dg_hyb > self.dg_threshold:
            is_viable_energy = False

        if mismatch_count == 0:
            verdict = SiteVerdict.PARFAIT
        elif not extensible:
            verdict = SiteVerdict.VETO_3P
        elif not is_viable_energy:
            verdict = SiteVerdict.NON_VIABLE
        else:
            verdict = SiteVerdict.TOLERABLE

        return SiteEvaluation(
            strain_id=strain_id,
            primer_name=primer.name,
            primer_role=primer.role,
            verdict=verdict,
            n_mismatches_count=mismatch_count,
            dg_hyb=dg_hyb,
            ddg=ddg,
            first_bad_pos=first_bad_pos,
            severity=severity,
            position=position,
            strand=strand
        )

    def analyze_strains(self, csv_records: List[Dict]) -> List[StrainVerdict]:
        """Analyse l'ensemble des souches à partir des enregistrements de sites candidats.
        Analyzes all strains from candidate site records.
        """
        if not csv_records:
            return []

        # Validation stricte du contrat CSV / Strict CSV contract validation
        first_row = csv_records[0]
        found_cols = set(first_row.keys())
        required_cols = {"strain_id", "strand", "n_mismatches"}
        has_primer_col = ("primer_role" in found_cols) or ("primer_name" in found_cols)

        missing = required_cols - found_cols
        if not has_primer_col:
            missing.add("primer_role (ou primer_name)")

        if missing:
            expected_display = {"strain_id", "primer_role", "strand", "position", "site_seq", "n_mismatches"}
            raise ValueError(
                f"Le fichier CSV de sites est incomplet. "
                f"Colonnes attendues : {sorted(expected_display)}, "
                f"colonnes trouvées : {sorted(found_cols)}. "
                f"Colonnes obligatoires manquantes : {sorted(missing)}."
            )

        # Regrouper les enregistrements par souche
        # Group records by strain
        records_by_strain: Dict[str, List[Dict]] = {}
        for r in csv_records:
            s_id = r["strain_id"]
            if s_id not in records_by_strain:
                records_by_strain[s_id] = []
            records_by_strain[s_id].append(r)

        # Amorces d'initiation requises pour l'amplification
        # Required initiation primer roles
        init_roles = {PrimerRole.F3, PrimerRole.B3, PrimerRole.FIP, PrimerRole.BIP}

        strain_verdicts: List[StrainVerdict] = []

        for s_id in self.fasta_dict.keys():
            records = records_by_strain.get(s_id, [])
            genome = self.fasta_dict[s_id].upper()

            # Évaluer tous les sites / variantes par rôle
            # Evaluate all sites / variants grouped by role
            role_evaluations: Dict[str, List[Tuple[SiteEvaluation, int]]] = {}

            for r in records:
                p_ident = r.get("primer_role") or r.get("primer_name")
                if not p_ident:
                    continue

                matching_primers = []
                if p_ident in self.primers_by_name:
                    matching_primers.append(self.primers_by_name[p_ident])
                else:
                    for p in self.primers:
                        if p.role.name == p_ident or p.role.value == p_ident or p.name == p_ident or (p.parent_name and p.parent_name == p_ident):
                            matching_primers.append(p)

                if not matching_primers:
                    continue

                csv_mm = int(r["n_mismatches"]) if "n_mismatches" in r and r["n_mismatches"] != "" else 0
                pos_val = int(r["position"]) if "position" in r and r["position"] != "" else None

                for primer in matching_primers:
                    req_len = len(primer.binding_domain)
                    site_seq = r.get("site_seq", "")
                    if len(site_seq) != req_len:
                        # Ceinture et bretelles : tentative de ré-extraction par coordonnées
                        if pos_val is not None and 0 <= pos_val <= len(genome) - req_len:
                            site_seq = genome[pos_val : pos_val + req_len]
                        if len(site_seq) != req_len:
                            raise ValueError(
                                f"Longueur de site ({len(site_seq)} nt) incompatible avec l'amorce "
                                f"{primer.name} ({req_len} nt) pour la souche '{s_id}'."
                            )

                    eval_obj = self.evaluate_site(primer, site_seq, r["strand"], strain_id=s_id, position=pos_val)

                    role_key = primer.role.value
                    if role_key not in role_evaluations:
                        role_evaluations[role_key] = []
                    role_evaluations[role_key].append((eval_obj, csv_mm))

            # Sélectionner la meilleure variante pour chaque rôle
            # Select the best variant for each role
            best_evaluations: Dict[str, SiteEvaluation] = {}
            min_csv_mms: Dict[str, int] = {}

            for role_key, ev_list in role_evaluations.items():
                best_ev, _ = max(ev_list, key=lambda pair: _evaluation_rank(pair[0]))
                best_evaluations[role_key] = best_ev
                # Règle par comptage indépendante : minimum des mésappariements sur toutes les lignes
                min_csv_mms[role_key] = min(mm for _, mm in ev_list)

            # Déterminer la viabilité d'amplification
            # Determine amplification viability
            is_amp_thermo = True
            is_amp_count = True
            limiting: Optional[str] = None

            # Vérifier toutes les amorces d'initiation
            # Check all initiation primers
            for role in init_roles:
                role_key = role.value
                if role_key not in best_evaluations:
                    # Site totalement introuvable / Completely absent site
                    is_amp_thermo = False
                    is_amp_count = False
                    if not limiting:
                        limiting = role_key
                    best_evaluations[role_key] = SiteEvaluation(
                        strain_id=s_id,
                        primer_name=role_key,
                        primer_role=role,
                        verdict=SiteVerdict.ABSENT,
                        n_mismatches_count=99,
                        dg_hyb=None,
                        ddg=None,
                        position=None,
                        strand=None
                    )
                    continue

                ev = best_evaluations[role_key]
                role_min_mm = min_csv_mms.get(role_key, ev.n_mismatches_count)

                # Règle thermodynamique (indépendante)
                if ev.verdict in (SiteVerdict.VETO_3P, SiteVerdict.ABSENT, SiteVerdict.NON_VIABLE):
                    is_amp_thermo = False
                    if not limiting:
                        limiting = role_key

                # Règle par comptage (indépendante)
                if role_min_mm > self.max_mismatches_count:
                    is_amp_count = False

            strain_verdicts.append(StrainVerdict(
                strain_id=s_id,
                evaluations=best_evaluations,
                is_amplifiable_thermo=is_amp_thermo,
                is_amplifiable_count=is_amp_count,
                limiting_primer_thermo=limiting
            ))

        return strain_verdicts
