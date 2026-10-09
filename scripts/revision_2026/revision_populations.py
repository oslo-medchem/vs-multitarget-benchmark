#!/usr/bin/env python3
"""Population-explicit re-analysis for the September 2026 MLST rewrite.

Addresses referee 1 point 1 (recompute the four carried-over targets) and
referee 2 point 2 (separate full-workflow performance from conditional ranking
performance, with explicit populations, counts, exclusions and uncertainty).

Three evaluation populations are defined and reported side by side:

  P1  end-to-end (intention-to-screen): every compound in the input library.
      Compounds a protocol never scored, whether excluded deliberately by the
      substructure filter or lost to a preparation or docking failure, are
      ranked last. This is the full-workflow endpoint.
  P2  docked-only: compounds each protocol scored, evaluated independently.
      The two protocols are on different compound sets here, so this is
      descriptive only and is never used for a paired test.
  P3  intersection: compounds scored by BOTH protocols. This is the
      conditional ranking endpoint and the only population on which the
      paired DeLong test is run.

Outputs data/revision_2026/populations.csv and intersection_delong.csv.
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
    ap.add_argument("--ablation-root", default=None)
    ap.add_argument("--out-dir", default="data/revision_2026")
    args = ap.parse_args()

    root = Path(args.data_root)
    abl = Path(args.ablation_root) if args.ablation_root else None
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    pop_rows, del_rows = [], []

    for t in rc.TARGETS:
        td = rc.load_target(t, root, abl)
        n_lib = len(td.labels)
        n_act_lib = sum(td.labels.values())

        # ---------------- P1 end-to-end, both protocols -------------------
        for prot, sc in (("naive", td.naive), ("skill", td.skill)):
            _, y, rank = rc.end_to_end_ranks(td.labels, sc, "strict")
            m = rc.all_metrics(y, rank)
            n_scored = len(set(sc) & set(td.labels))
            n_filtered = len(set(td.filtered) & set(td.labels)) if prot == "skill" else 0
            n_unscored = n_lib - n_scored
            pop_rows.append({
                "target": t, "population": "P1_end_to_end", "protocol": prot,
                "n_library": n_lib, "n_evaluated": n_lib,
                "n_actives_eval": n_act_lib, "n_decoys_eval": n_lib - n_act_lib,
                "n_scored": n_scored, "n_unscored": n_unscored,
                "n_excluded_by_filter": n_filtered,
                "n_technical_failure": n_unscored - n_filtered,
                **{k: v for k, v in m.items() if k not in ("n_total",)},
            })

        # ---------------- P2 docked-only ----------------------------------
        for prot, sc in (("naive", td.naive), ("skill", td.skill)):
            cids = sorted(c for c in td.labels if c in sc)
            y = np.array([td.labels[c] for c in cids])
            rank = np.array([-sc[c] for c in cids])
            m = rc.all_metrics(y, rank)
            pop_rows.append({
                "target": t, "population": "P2_docked_only", "protocol": prot,
                "n_library": n_lib, "n_evaluated": len(cids),
                "n_actives_eval": int(y.sum()), "n_decoys_eval": int(len(y) - y.sum()),
                "n_scored": len(cids), "n_unscored": n_lib - len(cids),
                "n_excluded_by_filter": "", "n_technical_failure": "",
                **{k: v for k, v in m.items() if k not in ("n_total",)},
            })

        # ---------------- P3 intersection + paired DeLong -----------------
        common = sorted(set(td.naive) & set(td.skill) & set(td.labels))
        y = np.array([td.labels[c] for c in common])
        rn = np.array([-td.naive[c] for c in common])
        rs = np.array([-td.skill[c] for c in common])
        for prot, rank in (("naive", rn), ("skill", rs)):
            m = rc.all_metrics(y, rank)
            pop_rows.append({
                "target": t, "population": "P3_intersection", "protocol": prot,
                "n_library": n_lib, "n_evaluated": len(common),
                "n_actives_eval": int(y.sum()), "n_decoys_eval": int(len(y) - y.sum()),
                "n_scored": len(common), "n_unscored": n_lib - len(common),
                "n_excluded_by_filter": "", "n_technical_failure": "",
                **{k: v for k, v in m.items() if k not in ("n_total",)},
            })

        delta, se, p, auc_n, auc_s = rc.delong_paired(y, rn, rs)
        del_rows.append({
            "target": t, "intersection_n": len(common),
            "n_actives": int(y.sum()), "n_decoys": int(len(y) - y.sum()),
            "auc_naive_intersection": round(auc_n, 6),
            "auc_skill_intersection": round(auc_s, 6),
            "delong_delta": round(delta, 6), "delong_se": round(se, 6),
            "ci_lo": round(delta - 1.96 * se, 6), "ci_hi": round(delta + 1.96 * se, 6),
            "delong_p": round(p, 6),
        })

    padj = rc.bh_adjust([r["delong_p"] for r in del_rows])
    for r, a in zip(del_rows, padj):
        r["delong_p_adj_BH"] = round(float(a), 6)

    meta = rc.dersimonian_laird([r["delong_delta"] for r in del_rows],
                                [r["delong_se"] for r in del_rows])

    with (out / "populations.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(pop_rows[0].keys()))
        w.writeheader()
        w.writerows(pop_rows)
    with (out / "intersection_delong.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(del_rows[0].keys()))
        w.writeheader()
        w.writerows(del_rows)
    (out / "meta_analysis_intersection.json").write_text(json.dumps(meta, indent=2))

    print(f"wrote {out/'populations.csv'} ({len(pop_rows)} rows)")
    print(f"wrote {out/'intersection_delong.csv'} ({len(del_rows)} rows)")
    print(json.dumps(meta, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
