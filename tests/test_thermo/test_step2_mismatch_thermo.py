"""Tests unitaires de l'Étape 2 : thermodynamique des mésappariements terminaux, consécutifs et wobble 3'.
Step 2 unit tests: terminal mismatches, consecutive symmetric internal loops, and 3' G·T wobble rules.
"""
from __future__ import annotations

import math
import pytest
import RNA
from Bio.SeqUtils import MeltingTemp as mt

from labcraft.diagnostics.enzyme import BST_2_0, get_enzyme
from labcraft.thermo.constants import R_GAS_KCAL_MOL_K
from labcraft.thermo.mismatch import (
    _MISMATCH_TERMINAL,
    classify_mismatch_pair,
    nn_duplex_energy,
    three_prime_extensible,
)
from labcraft.thermo.vienna import dna_session


def _tm_from_dh_ds(dh_kcal: float, ds_cal: float, ct_molar: float = 5e-8) -> float:
    """Calcule le Tm (°C) pour un hétéroduplexe non auto-complémentaire (convention Biopython dnac1=dnac2=25 nM -> k=12.5 nM).
    Computes Tm (°C) for a non-self-complementary heteroduplex (Biopython dnac1=dnac2=25 nM -> k=12.5 nM).
    """
    r_cal = R_GAS_KCAL_MOL_K * 1000.0
    k = ct_molar / 4.0
    return (dh_kcal * 1000.0) / (ds_cal + r_cal * math.log(k)) - 273.15


def test_terminal_mismatch_table_loaded_and_matches_biopython():
    """Vérifie que `_MISMATCH_TERMINAL` contient les 48 pas SantaLucia & Hicks 2004 et reproduit Biopython `DNA_TMM1`.
    Verifies `_MISMATCH_TERMINAL` loads all 48 SantaLucia & Hicks 2004 steps and matches Biopython `DNA_TMM1`.
    """
    assert len(_MISMATCH_TERMINAL) == 48

    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    seq = "CTCGTGTAGGAATGGGAG"
    bot_perf = "".join(comp[b] for b in seq)

    # 1. Mésappariements terminaux en 5' (pos 0 mutée) : exactitude analytique + accord Biopython (< 0,35 °C)
    # 1. 5'-terminal mismatches (pos 0 mutated): analytical exactness + Biopython agreement (< 0.35 °C)
    dh_core_5p, ds_core_5p, _ = nn_duplex_energy(seq[1:], bot_perf[1:], temp_celsius=63.0)
    for mut_base in ("A", "C", "T"):  # 'G' est le complément parfait de 'C'
        bot_5p = mut_base + bot_perf[1:]
        dh, ds, _ = nn_duplex_energy(seq, bot_5p, temp_celsius=63.0)
        expected_key = f"{bot_5p[1]}{bot_5p[0]}/{seq[1]}{seq[0]}"
        exp_dh, exp_ds = _MISMATCH_TERMINAL[expected_key]
        assert dh == pytest.approx(dh_core_5p + exp_dh, abs=1e-9)
        assert ds == pytest.approx(ds_core_5p + exp_ds, abs=1e-9)

        tm_lc = _tm_from_dh_ds(dh, ds, ct_molar=5e-8)
        tm_bp = mt.Tm_NN(
            seq, c_seq=bot_5p, nn_table=mt.DNA_NN3, tmm_table=mt.DNA_TMM1, saltcorr=0
        )
        assert tm_lc == pytest.approx(tm_bp, abs=0.35), (
            f"Écart Tm sur mismatch 5' ({mut_base}) : LabCraft={tm_lc:.2f} vs Biopython={tm_bp:.2f}"
        )

    # 2. Mésappariements terminaux en 3' (pos -1 mutée) : exactitude analytique + accord Biopython (< 0,35 °C)
    # 2. 3'-terminal mismatches (pos -1 mutated): analytical exactness + Biopython agreement (< 0.35 °C)
    dh_core_3p, ds_core_3p, _ = nn_duplex_energy(seq[:-1], bot_perf[:-1], temp_celsius=63.0)
    for mut_base in ("A", "G", "T"):  # 'C' est le complément parfait de 'G'
        bot_3p = bot_perf[:-1] + mut_base
        dh, ds, _ = nn_duplex_energy(seq, bot_3p, temp_celsius=63.0)
        expected_key = f"{seq[-2:]}/{bot_3p[-2:]}"
        exp_dh, exp_ds = _MISMATCH_TERMINAL[expected_key]
        assert dh == pytest.approx(dh_core_3p + exp_dh, abs=1e-9)
        assert ds == pytest.approx(ds_core_3p + exp_ds, abs=1e-9)

        tm_lc = _tm_from_dh_ds(dh, ds, ct_molar=5e-8)
        tm_bp = mt.Tm_NN(
            seq, c_seq=bot_3p, nn_table=mt.DNA_NN3, tmm_table=mt.DNA_TMM1, saltcorr=0
        )
        assert tm_lc == pytest.approx(tm_bp, abs=0.35), (
            f"Écart Tm sur mismatch 3' ({mut_base}) : LabCraft={tm_lc:.2f} vs Biopython={tm_bp:.2f}"
        )


def test_consecutive_tandem_mismatches_monotonic_and_coherent_with_vienna():
    """Vérifie que des mésappariements internes consécutifs (2×2, 3×3, 4×4) déstabilisent monotonement le duplexe.
    Verifies consecutive internal mismatches (2×2, 3×3, 4×4) monotonically destabilize the duplex.
    """
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    seq = "GCCACCTGGCAGCCACAGTA"
    bot_perf = "".join(comp[b] for b in seq)

    bot_1mm = bot_perf[:7] + "A" + bot_perf[8:]
    bot_2mm = bot_perf[:7] + "AA" + bot_perf[9:]
    bot_3mm = bot_perf[:7] + "AAA" + bot_perf[10:]
    bot_4mm = bot_perf[:7] + "AAAA" + bot_perf[11:]

    _, _, dg_0 = nn_duplex_energy(seq, bot_perf, temp_celsius=63.0)
    _, _, dg_1 = nn_duplex_energy(seq, bot_1mm, temp_celsius=63.0)
    _, _, dg_2 = nn_duplex_energy(seq, bot_2mm, temp_celsius=63.0)
    _, _, dg_3 = nn_duplex_energy(seq, bot_3mm, temp_celsius=63.0)
    _, _, dg_4 = nn_duplex_energy(seq, bot_4mm, temp_celsius=63.0)

    # Monotonie stricte : 0 < 1 < 2 < 3 < 4 mésappariements consécutifs
    # Strict monotonicity: 0 < 1 < 2 < 3 < 4 consecutive mismatches
    assert dg_0 < dg_1 < dg_2 < dg_3 < dg_4
    assert (dg_2 - dg_0) > (dg_1 - dg_0) + 2.0

    # Accord exact (< 0,02 kcal/mol) avec l'implémentation de référence SantaLucia & Hicks 2004 (Primer3 thal)
    # et écart borné vis-à-vis de Mathews 2004 dans ViennaRNA (qui utilise 3.1/3.9/4.2 au lieu de 3.6/4.4/4.8)
    # Exact agreement (< 0.02 kcal/mol) with the SantaLucia & Hicks 2004 reference C implementation (Primer3 thal)
    # and bounded offset vs Mathews 2004 in ViennaRNA (which uses 3.1/3.9/4.2 instead of Table 4's 3.6/4.4/4.8)
    import primer3

    p3_0 = primer3.calc_heterodimer(seq, bot_perf[::-1], mv_conc=1000.0, dv_conc=0.0, dntp_conc=0.0, temp_c=63.0)
    with dna_session(63.0):
        fc_p = RNA.fold_compound(seq + "&" + bot_perf[::-1])
        dg_v0 = fc_p.eval_structure("(" * len(seq) + ")" * len(seq))
        for k_label, bot_k, dg_k in [("2x2", bot_2mm, dg_2), ("3x3", bot_3mm, dg_3), ("4x4", bot_4mm, dg_4)]:
            p3_k = primer3.calc_heterodimer(seq, bot_k[::-1], mv_conc=1000.0, dv_conc=0.0, dntp_conc=0.0, temp_c=63.0)
            ddg_nn = dg_k - dg_0
            ddg_p3 = (p3_k.dg - p3_0.dg) / 1000.0
            assert abs(ddg_nn - ddg_p3) < 0.02, (
                f"Écart ΔΔG SantaLucia & Hicks 2004 sur boucle {k_label} : NN={ddg_nn:.3f} vs Primer3={ddg_p3:.3f}"
            )

            s_top = "".join("(" if comp[a] == c else "." for a, c in zip(seq, bot_k))
            s_bot = "".join(")" if comp[a] == c else "." for a, c in reversed(list(zip(seq, bot_k))))
            fc_k = RNA.fold_compound(seq + "&" + bot_k[::-1])
            dg_vk = fc_k.eval_structure(s_top + s_bot)
            ddg_vrna = dg_vk - dg_v0
            assert abs(ddg_nn - ddg_vrna) < 1.35, (
                f"Écart ΔΔG sur boucle {k_label} : NN={ddg_nn:.2f} vs ViennaRNA(Mathews2004)={ddg_vrna:.2f}"
            )


def test_all_mismatch_alignment_is_unstable():
    """Vérifie qu'un alignement sans aucune paire Watson-Crick retourne une énergie défavorable (> 0).
    Verifies an alignment with zero Watson-Crick pairs returns an unfavorable energy (> 0).
    """
    _, _, dg_none = nn_duplex_energy("AAAAAA", "AAAAAA", temp_celsius=65.0)
    _, _, dg_perf = nn_duplex_energy("AAAAAA", "TTTTTT", temp_celsius=65.0)
    assert dg_none > 0.0
    assert dg_none > dg_perf


def test_classify_mismatch_and_penultimate_gt_wobble():
    """Vérifie la classification chimique des paires et l'option de tolérance wobble G·T en position -2.
    Verifies chemical base-pair classification and penultimate (-2) G·T wobble tolerance option.
    """
    assert classify_mismatch_pair("G", "C") == "match"
    assert classify_mismatch_pair("G", "T") == "wobble"
    assert classify_mismatch_pair("T", "G") == "wobble"
    assert classify_mismatch_pair("A", "G") == "purine_purine"
    assert classify_mismatch_pair("C", "C") == "pyrimidine_pyrimidine"

    primer = "GCCACCTTAAGCCACAGTA"  # ... G(-2) T(-2? wait: -3=G, -2=T, -1=A)
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    bot_perf = "".join(comp[b] for b in primer)

    # En position -2, primer[-2] == 'T' -> template 'G' forme un wobble T:G ; template 'C' forme py:py T:C
    # At position -2, primer[-2] == 'T' -> template 'G' forms a T:G wobble; template 'C' forms T:C py:py
    bot_wobble_pos2 = bot_perf[:-2] + "G" + bot_perf[-1]
    bot_clash_pos2 = bot_perf[:-2] + "C" + bot_perf[-1]
    # En position -1, primer[-1] == 'A' ; testons aussi un primer finissant par G face à T en -1
    # At position -1, test a primer ending in G facing T at -1
    primer_end_g = primer[:-1] + "G"
    bot_wobble_pos1 = bot_perf[:-1] + "T"

    # 1. Comportement ARMS strict par défaut sur BST_2_0
    # 1. Default strict ARMS behavior on BST_2_0
    ext_def, pos_def, sev_def = three_prime_extensible(primer, bot_wobble_pos2, BST_2_0)
    assert (ext_def, pos_def, sev_def) == (False, 2, "block")

    # 2. Profil avec tolérance au wobble G·T en position pénultième (-2)
    # 2. Profile with penultimate (-2) G·T wobble tolerance enabled
    bst_wobble = get_enzyme({"name": "bst2.0", "allow_penultimate_gt_wobble": True})
    ext_w2, pos_w2, sev_w2 = three_prime_extensible(primer, bot_wobble_pos2, bst_wobble)
    assert (ext_w2, pos_w2, sev_w2) == (True, 2, "wobble_weak")

    # Une collision pyrimidine-pyrimidine T:C en -2 reste bloquée
    # A pyrimidine-pyrimidine T:C clash at -2 remains blocked
    ext_c2, pos_c2, sev_c2 = three_prime_extensible(primer, bot_clash_pos2, bst_wobble)
    assert (ext_c2, pos_c2, sev_c2) == (False, 2, "block")

    # Un wobble G:T à l'extrémité terminale stricte (-1) reste bloqué
    # A G:T wobble at the strict 3'-terminal base (-1) remains blocked
    ext_w1, pos_w1, sev_w1 = three_prime_extensible(primer_end_g, bot_wobble_pos1, bst_wobble)
    assert (ext_w1, pos_w1, sev_w1) == (False, 1, "block")
