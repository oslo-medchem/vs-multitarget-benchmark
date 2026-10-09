#!/usr/bin/env python3
"""Per-arm scoring rates for actives and decoys, and their relation to the end-to-end AUC.

The end-to-end population ranks every compound a protocol never scored below every
compound it did. Whether that helps or hurts an arm depends on which class of
compound goes unscored. This script writes, for every target and protocol, the
fraction of actives and the fraction of decoys that received a score, their
difference, and the tie-corrected end-to-end ROC AUC, and reports the correlation
between the difference and the AUC across all 22 arms.

Writes attrition_by_arm.csv.

Usage:
    python attrition.py --data-root <deposit root> --out-dir data/revision_2026
"""
import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import roc_auc_score

from revision_common import PRETTY, TARGETS, end_to_end_ranks, load_target

LABEL = dict(PRETTY)
LABEL.update({"hsp90": "HSP90α", "p38": "p38 MAPK"})


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--out-dir", default="data/revision_2026")
    args = ap.parse_args()
    root, out = Path(args.data_root), Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    rows = []
    for t in TARGETS:
        d = load_target(t, root)
        n_act = sum(d.labels.values())
        n_dec = len(d.labels) - n_act
        for cond, scores in (("naive", d.naive), ("skill", d.skill)):
            act = sum(1 for c, v in d.labels.items() if v == 1 and c in scores)
            dec = sum(1 for c, v in d.labels.items() if v == 0 and c in scores)
            _, y, rank = end_to_end_ranks(d.labels, scores, "strict")
            rows.append({
                "target": t, "target_label": LABEL[t], "protocol": cond,
                "n_actives": n_act, "n_decoys": n_dec,
                "actives_scored": act, "decoys_scored": dec,
                "pct_actives_scored": round(100 * act / n_act, 1),
                "pct_decoys_scored": round(100 * dec / n_dec, 1),
                "attrition_gap_pp": round(100 * act / n_act - 100 * dec / n_dec, 1),
                "auc_end_to_end": round(float(roc_auc_score(y, rank)), 4),
            })

    with (out / "attrition_by_arm.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    gap = np.array([r["attrition_gap_pp"] for r in rows])
    auc = np.array([r["auc_end_to_end"] for r in rows])
    pr, pp = pearsonr(gap, auc)
    sr, sp = spearmanr(gap, auc)
    below = [r for r in rows if r["auc_end_to_end"] < 0.5]
    print(f"wrote {out / 'attrition_by_arm.csv'} ({len(rows)} arms)")
    print(f"Pearson r = {pr:.3f} (p = {pp:.2e}); Spearman rho = {sr:.3f} (p = {sp:.2e})")
    print(f"arms with AUC below 0.5: {len(below)}; "
          f"of these, {sum(1 for r in below if r['attrition_gap_pp'] < 0)} have a negative gap")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
