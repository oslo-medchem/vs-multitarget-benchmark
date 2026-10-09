#!/usr/bin/env python3
"""Early-recognition metric panel on both evaluation populations.

Writes the two tables behind Additional file 1, Tables S12 and S13:

  TableS_early_recognition_metrics.csv               end-to-end population
  TableS_early_recognition_metrics_intersection.csv  intersection population

For every target and protocol: tie-corrected ROC AUC, BEDROC at alpha 20 and
80, enrichment factor at 1, 5 and 10 percent, LogAUC and AUPR, together with
the chance expectation of each metric for that library. Chance values are
analytic: AUPR under a random ordering equals the active fraction, BEDROC
under a random ordering (RIE = 1) follows Truchon and Bayly (2007), and the
enrichment factor is 1 by construction.

End-to-end rankings use the strict convention of revision_common
(unscored compounds strictly below every scored compound) and the pessimistic
tie policy (actives last within a tied block). Intersection rankings contain
only compounds both protocols scored.

Usage:
    python early_recognition.py --data-root <deposit root> --out-dir data/revision_2026
"""

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from rdkit.ML.Scoring.Scoring import CalcBEDROC, CalcEnrichment
from sklearn.metrics import average_precision_score, roc_auc_score

from revision_common import (
    PRETTY,
    TARGETS,
    calc_logauc,
    end_to_end_ranks,
    load_target,
    ranked_list,
)

LABEL = dict(PRETTY)
LABEL.update({"hsp90": "HSP90α", "p38": "p38 MAPK"})


def bedroc_chance(ra: float, alpha: float) -> float:
    """BEDROC of a random ordering: RIE = 1 in equation 36 of Truchon & Bayly."""
    num = ra * math.sinh(alpha / 2.0)
    den = math.cosh(alpha / 2.0) - math.cosh(alpha / 2.0 - alpha * ra)
    return num / den + 1.0 / (1.0 - math.exp(alpha * (1.0 - ra)))


def panel(y: np.ndarray, rank: np.ndarray) -> dict:
    lst = ranked_list(y, rank, "pessimistic")
    ra = float(y.mean())
    ef = CalcEnrichment(lst, 1, [0.01, 0.05, 0.10])
    return {
        "roc_auc": round(float(roc_auc_score(y, rank)), 4),
        "bedroc_a20": round(float(CalcBEDROC(lst, 1, 20.0)), 4),
        "bedroc_a20_chance": round(bedroc_chance(ra, 20.0), 4),
        "bedroc_a80": round(float(CalcBEDROC(lst, 1, 80.0)), 4),
        "bedroc_a80_chance": round(bedroc_chance(ra, 80.0), 4),
        "ef_1pct": round(float(ef[0]), 4),
        "ef_5pct": round(float(ef[1]), 4),
        "ef_10pct": round(float(ef[2]), 4),
        "logauc": round(calc_logauc(y, rank), 4),
        "aupr": round(float(average_precision_score(y, rank)), 4),
        "aupr_chance": round(ra, 4),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--out-dir", default="data/revision_2026")
    args = ap.parse_args()
    root, out = Path(args.data_root), Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    e2e_rows, int_rows = [], []
    for t in TARGETS:
        d = load_target(t, root)
        for cond, scores in (("naive", d.naive), ("skill", d.skill)):
            cids, y, rank = end_to_end_ranks(d.labels, scores, "strict")
            row = {
                "target": t,
                "target_label": LABEL[t],
                "condition": cond,
                "n_library": len(cids),
                "n_actives": int(y.sum()),
                "n_scored": sum(1 for c in cids if c in scores),
            }
            row.update(panel(y, rank))
            e2e_rows.append(row)

        both = sorted(set(d.naive) & set(d.skill))
        y = np.array([d.labels[c] for c in both])
        for cond, scores in (("naive", d.naive), ("skill", d.skill)):
            rank = np.array([-scores[c] for c in both])
            row = {
                "target": t,
                "target_label": LABEL[t],
                "condition": cond,
                "n_intersection": len(both),
                "n_actives": int(y.sum()),
            }
            row.update(panel(y, rank))
            int_rows.append(row)

    for name, rows in (
        ("TableS_early_recognition_metrics.csv", e2e_rows),
        ("TableS_early_recognition_metrics_intersection.csv", int_rows),
    ):
        with (out / name).open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {out / name} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
