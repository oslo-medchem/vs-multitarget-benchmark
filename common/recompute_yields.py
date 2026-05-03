#!/usr/bin/env python3
"""
Recompute pipeline yield numbers consistently from the source archive.

Reads the authoritative aggregate_metrics.csv and per-target filter_report.csv
files from the source tree, produces a single yields.csv with three
denominators per target:

    naive_full      = naive_n_docked / n_compounds         (full library)
    skill_post      = skill_n_docked / (n_compounds - filt) (post-filter library)
    skill_full      = skill_n_docked / n_compounds         (full library; apples-to-apples)
    filter_pct_lib  = filtered / n_compounds                (full library filter rate)
    filter_pct_act  = filtered_actives / 100                (active-set filter rate)

These are the numbers used in the manuscript abstract, Table 2, and Discussion.
Run from the project root; expects the source tree at:
    /data/CLAUDE_works/TEST_SKILLS/vina_docking_test/vs_benchmark/vs_benchmark_multitarget/

If the source tree path changes, set SRC_ROOT below.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

SRC_ROOT = Path("/data/CLAUDE_works/TEST_SKILLS/vina_docking_test/vs_benchmark/vs_benchmark_multitarget")


def main() -> int:
    agg_path = SRC_ROOT / "06_meta_analysis" / "aggregate_metrics.csv"
    if not agg_path.exists():
        print(f"ERROR: aggregate_metrics.csv not found at {agg_path}", file=sys.stderr)
        return 1

    by_target: dict[str, dict[str, dict]] = {}
    with agg_path.open() as f:
        for row in csv.DictReader(f):
            by_target.setdefault(row["target"], {})[row["protocol"]] = row

    out = []
    for target in sorted(by_target):
        n_lib = int(by_target[target]["naive"]["n_compounds"])
        n_naive = int(by_target[target]["naive"]["n_docked"])
        n_skill = int(by_target[target]["skill"]["n_docked"])
        filter_csv = SRC_ROOT / "targets" / target / "03_library_preparation" / "skill" / "pdbqt" / "filter_report.csv"
        if filter_csv.exists():
            with filter_csv.open() as f:
                rows = list(csv.DictReader(f))
            n_filt = len(rows)
            n_act_filt = sum(1 for r in rows if r.get("compound_id", "").startswith("ACT_"))
        else:
            n_filt = 0
            n_act_filt = 0
        skill_attempt = n_lib - n_filt
        out.append(
            {
                "target": target,
                "n_compounds": n_lib,
                "n_naive_docked": n_naive,
                "n_skill_docked": n_skill,
                "n_filtered": n_filt,
                "n_actives_filtered": n_act_filt,
                "naive_full_pct": round(100 * n_naive / n_lib, 2) if n_lib else None,
                "skill_post_pct": round(100 * n_skill / skill_attempt, 2) if skill_attempt else None,
                "skill_full_pct": round(100 * n_skill / n_lib, 2) if n_lib else None,
                "filter_lib_pct": round(100 * n_filt / n_lib, 2) if n_lib else None,
                "filter_act_pct": round(100 * n_act_filt / 100, 2) if n_act_filt else 0.0,
            }
        )

    out_path = Path(__file__).resolve().parent.parent / "yields.csv"
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(out[0]))
        writer.writeheader()
        writer.writerows(out)

    print(f"Wrote {out_path}")
    print()
    core = [r for r in out if r["target"] != "fpr2"]
    print("Manuscript-claim ranges (excluding fpr2, which used the older single-target setup):")
    print(f"  naive_full_pct  : {min(r['naive_full_pct'] for r in core):.1f}–{max(r['naive_full_pct'] for r in core):.1f}  (manuscript: 46–72%)")
    print(f"  skill_post_pct  : {min(r['skill_post_pct'] for r in core):.1f}–{max(r['skill_post_pct'] for r in core):.1f}  (manuscript: 86–94%)")
    print(f"  skill_full_pct  : {min(r['skill_full_pct'] for r in core):.1f}–{max(r['skill_full_pct'] for r in core):.1f}  (manuscript: 38–51%)")
    print(f"  filter_lib_pct  : {min(r['filter_lib_pct'] for r in core):.1f}–{max(r['filter_lib_pct'] for r in core):.1f}  (manuscript: 45–57%)")
    print(f"  filter_act_pct  : {min(r['filter_act_pct'] for r in core):.0f}–{max(r['filter_act_pct'] for r in core):.0f}  (manuscript: 38–67%)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
