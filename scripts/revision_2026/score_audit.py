#!/usr/bin/env python3
"""Audit of released docking-score ranges and the non-physical-score sensitivity analysis.

Writes:
  score_range_audit.csv             per target and score file: number scored, min, max,
                                    and counts below -20, below -50, above 0 and above 5 kcal/mol
  nonphysical_score_sensitivity.csv end-to-end ROC AUC per target and protocol as published,
                                    after treating every score below -20 kcal/mol as a
                                    preparation failure (unscored, ranked last), and after
                                    treating both tails outside -20 to +5 kcal/mol that way

The threshold of -20 kcal/mol is the value the naive route never reaches on any
target. Routes audited: naive, skill, the filter-ablation re-dock (skill_nofilter)
where present.

Usage:
    python score_audit.py --data-root <deposit root> [--ablation-root <dir>] --out-dir data/revision_2026
"""

import argparse
import csv
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

from revision_common import TARGETS, end_to_end_ranks, load_target

THRESH = -20.0
THRESH_HI = 5.0


def audit(route: str, target: str, scores: dict) -> dict:
    v = np.array(list(scores.values()), dtype=float)
    return {
        "target": target,
        "route": route,
        "n_scored": int(v.size),
        "min": round(float(v.min()), 1) if v.size else "",
        "max": round(float(v.max()), 1) if v.size else "",
        "lt_m20": int((v < -20).sum()),
        "lt_m50": int((v < -50).sum()),
        "gt_0": int((v > 0).sum()),
        "gt_5": int((v > 5).sum()),
    }


def e2e_auc(labels: dict, scores: dict) -> float:
    _, y, rank = end_to_end_ranks(labels, scores, "strict")
    return float(roc_auc_score(y, rank))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--ablation-root", default=None)
    ap.add_argument("--out-dir", default="data/revision_2026")
    args = ap.parse_args()
    root = Path(args.data_root)
    abl = Path(args.ablation_root) if args.ablation_root else None
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    audit_rows, sens_rows = [], []
    for t in TARGETS:
        d = load_target(t, root, abl)
        audit_rows.append(audit("naive", t, d.naive))
        audit_rows.append(audit("skill", t, d.skill))
        if d.nofilter:
            audit_rows.append(audit("skill_nofilter", t, d.nofilter))

        routes = (("naive", d.naive), ("skill", d.skill))
        lo = {r: {c: s for c, s in sc.items() if s >= THRESH} for r, sc in routes}
        both = {
            r: {c: s for c, s in sc.items() if THRESH <= s <= THRESH_HI}
            for r, sc in routes
        }
        pub = {r: e2e_auc(d.labels, sc) for r, sc in routes}
        clp = {r: e2e_auc(d.labels, lo[r]) for r in lo}
        bot = {r: e2e_auc(d.labels, both[r]) for r in both}
        sens_rows.append(
            {
                "target": t,
                "naive_as_published": round(pub["naive"], 5),
                "naive_clipped": round(clp["naive"], 5),
                "skill_as_published": round(pub["skill"], 5),
                "skill_clipped": round(clp["skill"], 5),
                "delta_published": round(pub["skill"] - pub["naive"], 5),
                "delta_clipped": round(clp["skill"] - clp["naive"], 5),
                "naive_clipped_both_tails": round(bot["naive"], 5),
                "skill_clipped_both_tails": round(bot["skill"], 5),
                "delta_clipped_both_tails": round(bot["skill"] - bot["naive"], 5),
            }
        )

    for name, rows in (
        ("score_range_audit.csv", audit_rows),
        ("nonphysical_score_sensitivity.csv", sens_rows),
    ):
        with (out / name).open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {out / name} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
