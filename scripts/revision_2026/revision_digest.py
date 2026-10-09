#!/usr/bin/env python3
"""Collect every revision number into one JSON digest.

data/revision_2026/DIGEST.json is the single source of truth for the macros in
scripts/build_docx_mlst.py. Nothing numeric should be typed into the manuscript
by hand; it should come from here, which in turn comes from the CSVs written by
the revision_*.py scripts.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

OUT = Path("data/revision_2026")


def main() -> int:
    d = {}
    inter = pd.read_csv(OUT / "intersection_delong.csv")
    e2e = pd.read_csv(OUT / "endtoend_delong.csv")
    dec = pd.read_csv(OUT / "decomposition.csv")
    prep = pd.read_csv(OUT / "prep_matched.csv")
    cont = pd.read_csv(OUT / "decoy_contamination.csv")
    sens = pd.read_csv(OUT / "decoy_sensitivity.csv")
    prov = pd.read_csv(OUT / "assay_provenance.csv")
    comp = pd.read_csv(OUT / "filter_composition.csv")
    hit = pd.read_csv(OUT / "hitlist_quality.csv")
    pot = pd.read_csv(OUT / "filtered_active_potency.csv")
    pops = pd.read_csv(OUT / "populations.csv")
    meta_i = json.loads((OUT / "meta_analysis_intersection.json").read_text())
    meta_e = json.loads((OUT / "meta_analysis_endtoend.json").read_text())
    meta_c = json.loads((OUT / "meta_analysis_components.json").read_text())

    d["meta_intersection"] = meta_i
    d["meta_endtoend"] = meta_e
    d["meta_components"] = meta_c

    d["intersection"] = {
        "n_targets": len(inter),
        "n_raw_sig": int((inter.delong_p < 0.05).sum()),
        "n_bh_sig": int((inter.delong_p_adj_BH < 0.05).sum()),
        "actives_min": int(inter.n_actives.min()),
        "actives_max": int(inter.n_actives.max()),
        "actives_median": float(inter.n_actives.median()),
        "n_min": int(inter.intersection_n.min()),
        "n_max": int(inter.intersection_n.max()),
    }
    d["endtoend"] = {
        "n_bh_sig": int((e2e.delong_p_adj_BH < 0.05).sum()),
        "n_bh_sig_negative": int(((e2e.delong_p_adj_BH < 0.05) & (e2e.delong_delta < 0)).sum()),
        "n_bh_sig_positive": int(((e2e.delong_p_adj_BH < 0.05) & (e2e.delong_delta > 0)).sum()),
        "delta_min": float(e2e.delong_delta.min()),
        "delta_max": float(e2e.delong_delta.max()),
        "pooled_zero_convention": float(e2e.delta_zero_convention.mean()),
    }
    d["components"] = {
        "prep_params_n_bh_sig": int((dec.prep_params_p_adj_BH < 0.05).sum()),
        "prep_params_min_p_adj": float(dec.prep_params_p_adj_BH.min()),
        "prep_params_ache": float(dec.loc[dec.target == "ache", "prep_params_delta_C_minus_A"].iloc[0]),
        "prep_params_ache_p": float(dec.loc[dec.target == "ache", "prep_params_p"].iloc[0]),
        "prep_params_cox2": float(dec.loc[dec.target == "cox2", "prep_params_delta_C_minus_A"].iloc[0]),
        "prep_params_cox2_p": float(dec.loc[dec.target == "cox2", "prep_params_p"].iloc[0]),
        "filter_n_sig": int((dec.filter_p_adj_BH < 0.05).sum()),
        "filter_n_sig_pos": int(((dec.filter_p_adj_BH < 0.05) & (dec.filter_delta_C_minus_B > 0)).sum()),
        "filter_n_sig_neg": int(((dec.filter_p_adj_BH < 0.05) & (dec.filter_delta_C_minus_B < 0)).sum()),
        "filter_delta_min": float(dec.filter_delta_C_minus_B.min()),
        "filter_delta_max": float(dec.filter_delta_C_minus_B.max()),
    }
    d["prep_matched"] = {
        "naive_prep_lo": float(prep.naive_prep_success_pct.min()),
        "naive_prep_hi": float(prep.naive_prep_success_pct.max()),
        "skill_prep_lo": float(prep.skill_prep_success_pct.min()),
        "skill_prep_hi": float(prep.skill_prep_success_pct.max()),
        "prep_gap_mean": float((prep.skill_prep_success_pct - prep.naive_prep_success_pct).mean()),
        "naive_dock_lo": float(prep.naive_dock_yield_pct.min()),
        "naive_dock_hi": float(prep.naive_dock_yield_pct.max()),
        "skill_dock_lo": float(prep.skill_dock_yield_pct.min()),
        "skill_dock_hi": float(prep.skill_dock_yield_pct.max()),
        "dock_gap_mean": float((prep.skill_dock_yield_pct - prep.naive_dock_yield_pct).mean()),
        "dock_gap_lo": float((prep.skill_dock_yield_pct - prep.naive_dock_yield_pct).min()),
        "dock_gap_hi": float((prep.skill_dock_yield_pct - prep.naive_dock_yield_pct).max()),
        "naive_lost_min": int(prep.naive_prepared_but_not_scored.min()),
        "naive_lost_max": int(prep.naive_prepared_but_not_scored.max()),
        "skill_lost_min": int(prep.skill_prepared_but_not_scored.min()),
        "skill_lost_max": int(prep.skill_prepared_but_not_scored.max()),
        "fallback_min": int(prep.skill_obabel_fallback_used.min()),
        "fallback_max": int(prep.skill_obabel_fallback_used.max()),
        "fallback_total": int(prep.skill_obabel_fallback_used.sum()),
    }
    d["decoys"] = {
        "same_target_max_pct": float(cont.pct_same_target.max()),
        "other_panel_lo_pct": float(cont.pct_any_other.min()),
        "other_panel_hi_pct": float(cont.pct_any_other.max()),
        "same_family_max_pct": float(cont.pct_same_family.max()),
        "same_family_max_target": str(cont.loc[cont.pct_same_family.idxmax(), "target"]),
        "sensitivity_max_abs_shift": float(sens.delta_shift.abs().max()),
    }
    d["provenance"] = {
        "binding_pct_lo": float(prov.pct_binding_records.min()),
        "binding_pct_hi": float(prov.pct_binding_records.max()),
        "binding_pct_median": float(prov.pct_binding_records.median()),
        "raw_records_lo": int(prov.n_raw_records.min()),
        "raw_records_hi": int(prov.n_raw_records.max()),
        "uniq_lo": int(prov.n_unique_qualifying_molecules.min()),
        "uniq_hi": int(prov.n_unique_qualifying_molecules.max()),
    }
    ca = comp[comp["set"] == "active"]
    cd = comp[comp["set"] == "decoy"]
    d["composition"] = {
        "active_removed_median_pct": float(ca.pct_removed.median()),
        "decoy_removed_median_pct": float(cd.pct_removed.median()),
        "active_removed_lo": float(ca[ca.target != "fpr2"].pct_removed.min()),
        "active_removed_hi": float(ca.pct_removed.max()),
        "active_QED_shift_median": float(ca.QED_shift.median()),
        "decoy_QED_shift_median": float(cd.QED_shift.median()),
    }
    hn = hit[hit.protocol == "naive"].set_index("target")
    hs = hit[hit.protocol == "skill"].set_index("target")
    d["hitlist"] = {
        "naive_alert_pct_median": float(hn.pct_with_alert.median()),
        "naive_alert_pct_lo": float(hn.pct_with_alert.min()),
        "naive_alert_pct_hi": float(hn.pct_with_alert.max()),
        "naive_qed_median": float(hn.median_QED.median()),
        "skill_qed_median": float(hs.median_QED.median()),
        "naive_actives_total": int(hn.n_true_actives_recovered.sum()),
        "skill_actives_total": int(hs.n_true_actives_recovered.sum()),
    }
    d["potency"] = {
        "diff_min": float(pot.difference.min()),
        "diff_max": float(pot.difference.max()),
        "min_p_adj": float(pot.mannwhitney_p_adj_BH.min()),
        "n_sig_after_bh": int((pot.mannwhitney_p_adj_BH < 0.05).sum()),
    }
    p1 = pops[pops.population == "P1_end_to_end"]
    p3 = pops[pops.population == "P3_intersection"]
    d["tie_bias"] = {
        "endtoend_mean_inflation": float((p1.roc_auc_legacy_ties - p1.roc_auc).mean()),
        "endtoend_max_inflation": float((p1.roc_auc_legacy_ties - p1.roc_auc).max()),
        "intersection_max_inflation": float((p3.roc_auc_legacy_ties - p3.roc_auc).max()),
    }

    (OUT / "DIGEST.json").write_text(json.dumps(d, indent=2))
    print(json.dumps(d, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
