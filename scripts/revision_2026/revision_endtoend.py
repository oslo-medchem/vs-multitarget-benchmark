#!/usr/bin/env python3
"""Paired end-to-end (intention-to-screen) DeLong test and meta-analysis.

Referee 2 point 2 notes that ranking deliberately excluded or failed compounds
last is a legitimate end-to-end retrieval endpoint rather than an artefact, and
asks for the primary evaluation population to be declared. This script runs the
same paired DeLong machinery used for the intersection population on the full
input library, with every compound a protocol did not score placed at worst
rank. Both protocols are therefore evaluated on an identical population of
identical size, which makes the test paired and the comparison like-for-like.

Outputs data/revision_2026/endtoend_delong.csv and
data/revision_2026/meta_analysis_endtoend.json.
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--out-dir", default="data/revision_2026")
    args = ap.parse_args()

    root = Path(args.data_root)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    rows = []
    for t in rc.TARGETS:
        td = rc.load_target(t, root)
        cids, y, rn = rc.end_to_end_ranks(td.labels, td.naive, "strict")
        _, _, rs = rc.end_to_end_ranks(td.labels, td.skill, "strict")
        delta, se, p, auc_n, auc_s = rc.delong_paired(y, rn, rs)
        # sensitivity: the simpler rank-key-zero convention
        _, _, rn0 = rc.end_to_end_ranks(td.labels, td.naive, "zero")
        _, _, rs0 = rc.end_to_end_ranks(td.labels, td.skill, "zero")
        d0, se0, p0, an0, as0 = rc.delong_paired(y, rn0, rs0)
        rows.append({
            "target": t, "n_library": len(cids),
            "n_actives": int(y.sum()), "n_decoys": int(len(y) - y.sum()),
            "n_naive_unscored": int(np.sum(rn == 0.0)),
            "n_skill_unscored": int(np.sum(rs == 0.0)),
            "auc_naive_end2end": round(auc_n, 6),
            "auc_skill_end2end": round(auc_s, 6),
            "delong_delta": round(delta, 6), "delong_se": round(se, 6),
            "ci_lo": round(delta - 1.96 * se, 6), "ci_hi": round(delta + 1.96 * se, 6),
            "delong_p": round(p, 6),
            "delta_zero_convention": round(d0, 6),
            "p_zero_convention": round(p0, 6),
        })

    padj = rc.bh_adjust([r["delong_p"] for r in rows])
    for r, a in zip(rows, padj):
        r["delong_p_adj_BH"] = round(float(a), 6)

    meta = rc.dersimonian_laird([r["delong_delta"] for r in rows],
                                [r["delong_se"] for r in rows])

    with (out / "endtoend_delong.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    (out / "meta_analysis_endtoend.json").write_text(json.dumps(meta, indent=2))

    for r in rows:
        print(f"{r['target']:9s} dAUC={r['delong_delta']:+.4f} "
              f"[{r['ci_lo']:+.3f},{r['ci_hi']:+.3f}] p={r['delong_p']:.4f} "
              f"p_adj={r['delong_p_adj_BH']:.4f}")
    print(json.dumps(meta, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
