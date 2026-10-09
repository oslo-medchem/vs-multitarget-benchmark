#!/usr/bin/env python3
"""Two-component decomposition and matched-input preparation robustness.

Referee 1 point 2 asks for a stepwise ablation that separates preparation-tool
effects from parameter effects. Referee 2 points 3 and 4 ask that mechanistic
claims be limited to what the design can isolate, and that preparation
robustness be evaluated on matched input molecules.

The benchmark contains three executed conditions:

  A  naive          naive preparation  + naive parameters  + no filter
  B  skill          skill preparation  + skill parameters  + PAINS/Brenk filter
  C  skill_nofilter skill preparation  + skill parameters  + no filter
                    (the true filter ablation: the filtered subset re-prepared
                     with apply_filter=False and re-docked at 25 A / exh. 32)

Two components can therefore be isolated exactly:

  FILTER      B vs C. Identical preparation and docking; the only difference is
              whether the substructure filter removed compounds before docking.
              Evaluated end to end on the full library, because the filter acts
              on the population rather than on any score.
  PREP+PARAMS A vs C. Neither condition applies a filter, so the comparison is
              free of the filter. Evaluated as a paired DeLong test on the
              compounds both conditions scored.

What CANNOT be isolated from these three conditions is preparation alone versus
docking parameters alone. That needs the two cross conditions, skill
preparation with naive parameters and naive preparation with skill parameters,
which were not run. scripts/run_stepwise_ablation.py defines them.

Outputs decomposition.csv and prep_matched.csv in data/revision_2026/.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import revision_common as rc

ABL_TARGETS = [t for t in rc.TARGETS if t != "fpr2"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--ablation-root", required=True)
    ap.add_argument("--out-dir", default="data/revision_2026")
    args = ap.parse_args()

    root, abl = Path(args.data_root), Path(args.ablation_root)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    dec_rows, prep_rows = [], []

    for t in ABL_TARGETS:
        td = rc.load_target(t, root, abl)
        # Condition C: skill scores plus the re-docked previously filtered subset
        skill_nf = dict(td.skill)
        n_added = 0
        for cid, s in td.nofilter.items():
            if cid not in skill_nf:
                skill_nf[cid] = s
                n_added += 1

        # ---- FILTER component: B vs C, end to end on the full library -----
        _, y_e2e, r_B = rc.end_to_end_ranks(td.labels, td.skill, "strict")
        _, _, r_C = rc.end_to_end_ranks(td.labels, skill_nf, "strict")
        d_filt, se_filt, p_filt, auc_B, auc_C = rc.delong_paired(y_e2e, r_B, r_C)

        # ---- PREP+PARAMS component: A vs C, paired on shared scored set ----
        common = sorted(set(td.naive) & set(skill_nf) & set(td.labels))
        y_c = np.array([td.labels[c] for c in common])
        rA = np.array([-td.naive[c] for c in common])
        rC2 = np.array([-skill_nf[c] for c in common])
        d_pp, se_pp, p_pp, auc_A, auc_C2 = rc.delong_paired(y_c, rA, rC2)

        dec_rows.append({
            "target": t,
            "n_library": len(td.labels),
            "n_skill_scored_with_filter": len(td.skill),
            "n_skill_scored_no_filter": len(skill_nf),
            "n_redocked_from_ablation": n_added,
            "auc_B_skill_filter_end2end": round(auc_B, 6),
            "auc_C_skill_nofilter_end2end": round(auc_C, 6),
            "filter_delta_C_minus_B": round(d_filt, 6),
            "filter_se": round(se_filt, 6),
            "filter_p": round(p_filt, 6),
            "intersection_n_A_vs_C": len(common),
            "n_actives_A_vs_C": int(y_c.sum()),
            "auc_A_naive_intersection": round(auc_A, 6),
            "auc_C_skill_nofilter_intersection": round(auc_C2, 6),
            "prep_params_delta_C_minus_A": round(d_pp, 6),
            "prep_params_se": round(se_pp, 6),
            "prep_params_ci_lo": round(d_pp - 1.96 * se_pp, 6),
            "prep_params_ci_hi": round(d_pp + 1.96 * se_pp, 6),
            "prep_params_p": round(p_pp, 6),
        })

        # ---- matched-input preparation robustness -------------------------
        # The molecules the skill protocol actually attempted: library minus
        # the compounds its substructure filter removed. Both protocols are
        # scored on this identical molecule set.
        attempted = set(td.labels) - set(td.filtered)
        n_att = len(attempted)
        naive_prep_ok = sum(1 for c in attempted
                            if td.prep_naive.get(c, ("", ""))[0] == "True")
        naive_prep_fail = sum(1 for c in attempted
                              if td.prep_naive.get(c, ("", ""))[0] == "False")
        skill_prep_ok = sum(1 for c in attempted
                            if td.prep_skill.get(c, ("", ""))[0] in ("ok", "ok_obabel_fallback"))
        skill_fallback = sum(1 for c in attempted
                             if td.prep_skill.get(c, ("", ""))[0] == "ok_obabel_fallback")
        skill_prep_fail = sum(1 for c in attempted
                              if td.prep_skill.get(c, ("", ""))[0] == "failed")
        naive_docked = len(attempted & set(td.naive))
        skill_docked = len(attempted & set(td.skill))
        prep_rows.append({
            "target": t,
            "n_library": len(td.labels),
            "n_filtered_out_by_skill": len(td.filtered),
            "n_matched_attempted_set": n_att,
            "naive_prep_ok": naive_prep_ok,
            "naive_prep_failed": naive_prep_fail,
            "naive_prep_success_pct": round(100 * naive_prep_ok / n_att, 2) if n_att else "",
            "skill_prep_ok": skill_prep_ok,
            "skill_prep_failed": skill_prep_fail,
            "skill_prep_success_pct": round(100 * skill_prep_ok / n_att, 2) if n_att else "",
            "skill_obabel_fallback_used": skill_fallback,
            "naive_docked_on_matched": naive_docked,
            "skill_docked_on_matched": skill_docked,
            "naive_dock_yield_pct": round(100 * naive_docked / n_att, 2) if n_att else "",
            "skill_dock_yield_pct": round(100 * skill_docked / n_att, 2) if n_att else "",
            "naive_prepared_but_not_scored": naive_prep_ok - naive_docked,
            "skill_prepared_but_not_scored": skill_prep_ok - skill_docked,
        })

    padj_pp = rc.bh_adjust([r["prep_params_p"] for r in dec_rows])
    padj_f = rc.bh_adjust([r["filter_p"] for r in dec_rows])
    for r, a, b in zip(dec_rows, padj_pp, padj_f):
        r["prep_params_p_adj_BH"] = round(float(a), 6)
        r["filter_p_adj_BH"] = round(float(b), 6)

    meta_pp = rc.dersimonian_laird([r["prep_params_delta_C_minus_A"] for r in dec_rows],
                                   [r["prep_params_se"] for r in dec_rows])
    meta_f = rc.dersimonian_laird([r["filter_delta_C_minus_B"] for r in dec_rows],
                                  [r["filter_se"] for r in dec_rows])

    for name, rows in (("decomposition.csv", dec_rows), ("prep_matched.csv", prep_rows)):
        with (out / name).open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    (out / "meta_analysis_components.json").write_text(json.dumps(
        {"prep_plus_params_C_minus_A": meta_pp, "filter_C_minus_B": meta_f}, indent=2))

    print("PREP+PARAMS (skill prep+params vs naive, no filter in either arm):")
    for r in dec_rows:
        print(f"  {r['target']:9s} d={r['prep_params_delta_C_minus_A']:+.4f} "
              f"[{r['prep_params_ci_lo']:+.3f},{r['prep_params_ci_hi']:+.3f}] "
              f"p={r['prep_params_p']:.4f} p_adj={r['prep_params_p_adj_BH']:.4f} "
              f"(n={r['intersection_n_A_vs_C']}, act={r['n_actives_A_vs_C']})")
    print(json.dumps(meta_pp, indent=2))
    print("\nFILTER (removing the filter, end to end):")
    for r in dec_rows:
        print(f"  {r['target']:9s} d={r['filter_delta_C_minus_B']:+.4f} p={r['filter_p']:.4g}")
    print(json.dumps(meta_f, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
