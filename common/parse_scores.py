#!/usr/bin/env python3
"""
Parse Uni-Dock output PDBQT files and extract best docking score per compound.

Usage:
    python parse_scores.py --target-dir targets/egfr/

Outputs:
    target_dir/04_docking/naive/scores_naive.csv
    target_dir/04_docking/skill/scores_skill.csv

Compounds with no docking result get score=0.0 (worst rank).
"""
import argparse
import csv
import re
import sys
from pathlib import Path

# Patterns for Uni-Dock output PDBQT
SCORE_PATTERN = re.compile(r"REMARK VINA RESULT:\s+([-\d.]+)", re.MULTILINE)
MODEL_PATTERN = re.compile(r"MODEL\s+1.*?ENDMDL", re.DOTALL)


def extract_best_score(pdbqt_path):
    """Extract best (mode 1) Vina score from docked PDBQT."""
    try:
        content = open(pdbqt_path).read()
        # Try REMARK VINA RESULT pattern
        m = SCORE_PATTERN.search(content)
        if m:
            return float(m.group(1))
        # Try "REMARK minimizedAffinity"
        m2 = re.search(r"REMARK minimizedAffinity\s+([-\d.]+)", content)
        if m2:
            return float(m2.group(1))
        # Try plain affinity line
        m3 = re.search(r"affinity.*?([-\d]+\.\d+)", content)
        if m3:
            return float(m3.group(1))
    except Exception:
        pass
    return None


def parse_protocol(protocol, results_dir, output_csv, all_ids):
    """Parse all docked files in results_dir, write CSV."""
    scores = {}

    if not results_dir.exists():
        print(f"  WARNING: {results_dir} does not exist")
    else:
        pdbqt_files = list(results_dir.glob("*.pdbqt"))
        print(f"  Found {len(pdbqt_files)} docked files in {results_dir}")

        for pf in pdbqt_files:
            # Extract compound ID from filename: ACT_001_out.pdbqt -> ACT_001
            stem = pf.stem
            cid = stem.replace("_out", "").replace("_docked", "")
            score = extract_best_score(pf)
            if score is not None:
                scores[cid] = score

    print(f"  Parsed scores: {len(scores)}")

    # Write CSV: all compounds, missing get 0.0
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["compound_id", "score"])
        writer.writeheader()
        for cid in all_ids:
            writer.writerow({
                "compound_id": cid,
                "score": scores.get(cid, 0.0),
            })

    n_docked = sum(1 for cid in all_ids if cid in scores)
    n_missing = len(all_ids) - n_docked
    print(f"  Docked: {n_docked}, Missing (score=0): {n_missing}")
    print(f"  Written: {output_csv}")


def main():
    parser = argparse.ArgumentParser(
        description="Parse Uni-Dock output and extract best docking scores."
    )
    parser.add_argument("--target-dir", required=True,
                        help="Path to target directory (e.g. targets/egfr/)")
    args = parser.parse_args()

    target_dir = Path(args.target_dir).resolve()
    docking_dir = target_dir / "04_docking"
    labels_csv = target_dir / "03_library_preparation" / "library_labels.csv"

    if not labels_csv.exists():
        print(f"ERROR: {labels_csv} not found")
        sys.exit(1)

    # Load ground truth compound IDs
    rows = list(csv.DictReader(open(labels_csv)))
    all_ids = [r["compound_id"] for r in rows]
    print(f"Total compounds: {len(all_ids)}")

    protocols = {
        "naive": {
            "results_dir": docking_dir / "naive" / "results",
            "output_csv": docking_dir / "naive" / "scores_naive.csv",
        },
        "skill": {
            "results_dir": docking_dir / "skill" / "results",
            "output_csv": docking_dir / "skill" / "scores_skill.csv",
        },
    }

    for protocol, cfg in protocols.items():
        print(f"\n=== {protocol.upper()} ===")
        parse_protocol(
            protocol=protocol,
            results_dir=cfg["results_dir"],
            output_csv=cfg["output_csv"],
            all_ids=all_ids,
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
