#!/usr/bin/env python3
"""
Build a partial true_ablation_results.json from per-target scores_skill_nofilter.csv
files that have already completed, without waiting for the full multi-target stream
to finish.

Use case: COX-2 ablation finishes in ~3 hours, but the rest of stream B (bace1, dpp4,
egfr, esr1) needs another 8-12 hours. Run this to capture the COX-2 result and update
the manuscript immediately, then re-run later when more targets finish.

Reads: filter_ablation_data/<target>/skill_nofilter/scores_skill_nofilter.csv
       (one per completed target)

Writes: filter_ablation_data/true_ablation_results.json
        — replaced with the union of all currently-complete targets.

Then chain into integrate_ablation_results.py:
    python scripts/build_partial_results.py && python scripts/integrate_ablation_results.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from run_true_filter_ablation import compute_full_auc, parse_scores  # type: ignore


ALL_TARGETS = ["ache", "bace1", "cdk2", "cox2", "dpp4", "egfr", "esr1", "fpr2", "hsp90", "p38", "thrombin"]


def main() -> int:
    abl_root = ROOT / "filter_ablation_data"
    results = []
    for target in ALL_TARGETS:
        scores_csv = abl_root / target / "skill_nofilter" / "scores_skill_nofilter.csv"
        if scores_csv.exists():
            scores = {}
            for r in csv.DictReader(scores_csv.open()):
                try:
                    scores[r["compound_id"]] = float(r["score"])
                except (ValueError, KeyError):
                    pass
            print(f"{target}: {len(scores)} ablation scores loaded")
            try:
                r = compute_full_auc(target, scores)
                results.append(r)
            except Exception as e:
                print(f"  {target}: compute_full_auc failed: {e}")
                results.append({"target": target, "error": str(e)})
        else:
            # Check for partial results from running docking (PDBQTs without scores csv)
            results_dir = abl_root / target / "skill_nofilter" / "results"
            if results_dir.exists():
                pdbqts = [p for p in results_dir.glob("*_out.pdbqt")]
                if pdbqts:
                    print(f"{target}: docking still in progress ({len(pdbqts)} PDBQTs); skipping")
                else:
                    print(f"{target}: not started")

    out = abl_root / "true_ablation_results.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nWrote {out} with {len(results)} target(s)")
    print("Now run: conda run -n vina_dock python scripts/integrate_ablation_results.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
