#!/usr/bin/env python3
"""
True filter ablation: re-prepare and re-dock PAINS/Brenk-filtered compounds
under skill protocol parameters (RDKit ETKDGv3 + Meeko, 25 Å box, exhaustiveness 32),
then merge with existing skill scores to obtain a full-library skill ROC AUC.

This is the definitive answer to the proxy filter ablation: does the COX-2 / AChE
gain survive when filtered compounds are scored under the skill pipeline (not
naive)?

For each target:
  1. Read library_labels.csv and filter_report.csv → list of filtered compounds + SMILES.
  2. Re-prepare them with apply_filter=False (PDBQTs in skill_nofilter/pdbqt/).
  3. Symlink the existing receptor_skill.pdbqt; build ligand list.
  4. Dock with skill parameters (size=25, exhaustiveness=32, vina scoring).
  5. Parse scores → scores_skill_nofilter.csv (just the previously-filtered subset).
  6. Merge with the existing scores_skill.csv to get a full-library "skill+ablation" score table.
  7. Compute ROC AUC on this full library and ΔAUC vs naive.
  8. Run DeLong's test on the intersection (compounds docked by both naive and skill+ablation).

Usage:
    python scripts/run_true_filter_ablation.py --targets ache cox2 dpp4
    python scripts/run_true_filter_ablation.py --targets all
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Iterable

# ---- Paths ----
SRC_ROOT = Path("/data/CLAUDE_works/TEST_SKILLS/vina_docking_test/vs_benchmark/vs_benchmark_multitarget")
PROJECT_ROOT = Path(__file__).resolve().parent.parent
ABL_ROOT = PROJECT_ROOT / "filter_ablation_data"  # outputs go here so source tree is untouched
ABL_ROOT.mkdir(parents=True, exist_ok=True)

# Make common scripts importable
COMMON = SRC_ROOT / "common"
sys.path.insert(0, str(COMMON))

import multiprocessing as mp
from functools import partial


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--targets", nargs="+", required=True,
                   help='Target names or "all" for all 11 targets')
    p.add_argument("--gpu", type=int, default=0, help="CUDA device index")
    p.add_argument("--max-workers", type=int, default=None, help="Prep parallel workers")
    p.add_argument("--skip-prep", action="store_true", help="Skip prep step (assume done)")
    p.add_argument("--skip-dock", action="store_true", help="Skip dock step (assume done)")
    return p.parse_args()


ALL_TARGETS = ["ache", "bace1", "cdk2", "cox2", "dpp4", "egfr", "esr1", "fpr2", "hsp90", "p38", "thrombin"]


# --- Box centers and skill params per target ---
def load_box_center(target: str) -> tuple[float, float, float]:
    """Read centre from existing skill run_unidock.sh."""
    script = SRC_ROOT / "targets" / target / "04_docking" / "skill" / "run_unidock.sh"
    cx = cy = cz = None
    for line in script.read_text().splitlines():
        if line.startswith("CENTER_X="):
            cx = float(line.split("=", 1)[1].strip())
        elif line.startswith("CENTER_Y="):
            cy = float(line.split("=", 1)[1].strip())
        elif line.startswith("CENTER_Z="):
            cz = float(line.split("=", 1)[1].strip())
    assert all(v is not None for v in (cx, cy, cz)), f"failed to parse box centre for {target}"
    return cx, cy, cz


# --- Prep step (no filter) ---
def prep_target(target: str, max_workers: int | None) -> Path:
    """Re-prepare the previously-filtered compounds with apply_filter=False."""
    out_dir = ABL_ROOT / target / "skill_nofilter"
    pdbqt_dir = out_dir / "pdbqt"
    pdbqt_dir.mkdir(parents=True, exist_ok=True)

    # Load filtered compound list + their SMILES
    filt_csv = SRC_ROOT / "targets" / target / "03_library_preparation" / "skill" / "pdbqt" / "filter_report.csv"
    if not filt_csv.exists():
        print(f"[{target}] no filter_report — skipping prep")
        return pdbqt_dir

    filtered_ids = {r["compound_id"] for r in csv.DictReader(filt_csv.open())}
    lib_csv = SRC_ROOT / "targets" / target / "03_library_preparation" / "library_labels.csv"
    rows = [r for r in csv.DictReader(lib_csv.open()) if r["compound_id"] in filtered_ids]
    print(f"[{target}] re-preparing {len(rows)} previously-filtered compounds (apply_filter=False)")

    # Use the existing prepare_one
    from prepare_skill import prepare_one  # type: ignore

    prep = partial(prepare_one, output_dir=pdbqt_dir, apply_filter=False)
    if max_workers is None:
        max_workers = max(1, mp.cpu_count() - 2)

    log_path = out_dir / "preparation_log.csv"
    t0 = time.time()
    with mp.Pool(max_workers) as pool, log_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["compound_id", "status", "pdbqt_size", "note"])
        for i, result in enumerate(pool.imap_unordered(prep, rows, chunksize=20), 1):
            w.writerow(result)
            if i % 200 == 0:
                print(f"  [{target}] prep {i}/{len(rows)} ({time.time()-t0:.0f}s)")
    print(f"[{target}] prep done: {time.time()-t0:.0f}s, log → {log_path}")
    return pdbqt_dir


# --- Dock step ---
def dock_target(target: str, gpu: int) -> Path:
    """Dock the prepared compounds with skill parameters."""
    out_dir = ABL_ROOT / target / "skill_nofilter"
    pdbqt_dir = out_dir / "pdbqt"
    results_dir = out_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    receptor = SRC_ROOT / "targets" / target / "04_docking" / "skill" / "receptor_skill.pdbqt"
    if not receptor.exists():
        raise FileNotFoundError(f"missing skill receptor for {target}: {receptor}")

    # Build ligand list (pdbqt files >50 bytes)
    lig_list = out_dir / "ligand_list.txt"
    valid_ligs = sorted(p for p in pdbqt_dir.glob("*.pdbqt") if p.stat().st_size > 50)
    lig_list.write_text("\n".join(str(p) for p in valid_ligs) + "\n")
    print(f"[{target}] {len(valid_ligs)} ligands to dock (skill+ablation)")

    if not valid_ligs:
        return results_dir

    cx, cy, cz = load_box_center(target)
    # Use the existing batched runner
    cmd = [
        sys.executable, str(COMMON / "run_unidock_batched.py"),
        "--receptor", str(receptor),
        "--ligand-list", str(lig_list),
        "--center-x", str(cx), "--center-y", str(cy), "--center-z", str(cz),
        "--size", "25", "--exhaustiveness", "32",
        "--results-dir", str(results_dir),
        "--gpu", str(gpu),
        "--batch-size", "100",
    ]
    log_path = out_dir / "unidock.log"
    print(f"[{target}] docking → log {log_path}")
    t0 = time.time()
    with log_path.open("w") as f:
        proc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    print(f"[{target}] dock done: {time.time()-t0:.0f}s, exit={proc.returncode}")
    return results_dir


# --- Parse scores from PDBQT outputs ---
def parse_scores(target: str) -> dict[str, float]:
    """Walk results dir and extract first-mode Vina score per compound."""
    out_dir = ABL_ROOT / target / "skill_nofilter"
    results_dir = out_dir / "results"
    scores: dict[str, float] = {}
    for pdbqt in results_dir.glob("*_out.pdbqt"):
        compound_id = pdbqt.stem.removesuffix("_out")
        try:
            with pdbqt.open() as f:
                for line in f:
                    if line.startswith("REMARK VINA RESULT"):
                        score = float(line.split()[3])
                        scores[compound_id] = score
                        break
        except Exception as e:
            print(f"  parse error {pdbqt}: {e}")
    csv_out = out_dir / "scores_skill_nofilter.csv"
    with csv_out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["compound_id", "score"])
        for cid, s in sorted(scores.items()):
            w.writerow([cid, s])
    print(f"[{target}] parsed {len(scores)} scores → {csv_out}")
    return scores


# --- Compute new full-library AUC ---
def compute_full_auc(target: str, ablation_scores: dict[str, float]) -> dict:
    """Merge skill scores with ablation scores, compute full-library AUC and DeLong vs naive."""
    from sklearn.metrics import roc_auc_score
    import numpy as np

    # Load existing skill scores (post-filter); for filtered compounds, score=0
    skill_scores_csv = SRC_ROOT / "targets" / target / "04_docking" / "skill" / "scores_skill.csv"
    skill_full = {}
    for r in csv.DictReader(skill_scores_csv.open()):
        try:
            s = float(r["score"])
            if s != 0.0:
                skill_full[r["compound_id"]] = s
        except (ValueError, KeyError):
            pass

    # Add ablation scores
    n_added = 0
    for cid, s in ablation_scores.items():
        if cid not in skill_full and s < 0:
            skill_full[cid] = s
            n_added += 1

    # Naive scores
    naive_csv = SRC_ROOT / "targets" / target / "04_docking" / "naive" / "scores_naive.csv"
    naive_full = {}
    for r in csv.DictReader(naive_csv.open()):
        try:
            s = float(r["score"])
            if s != 0.0:
                naive_full[r["compound_id"]] = s
        except (ValueError, KeyError):
            pass

    # Labels
    labels_csv = SRC_ROOT / "targets" / target / "03_library_preparation" / "library_labels.csv"
    labels = {r["compound_id"]: r["label"] for r in csv.DictReader(labels_csv.open())}

    def auc(scores: dict[str, float]) -> float:
        y_true, y_score = [], []
        for cid, s in scores.items():
            lab = labels.get(cid)
            if lab not in ("active", "decoy") or s >= 0:
                continue
            y_true.append(1 if lab == "active" else 0)
            y_score.append(-s)
        if sum(y_true) < 2 or len(y_true) - sum(y_true) < 2:
            return float("nan")
        return roc_auc_score(y_true, y_score)

    auc_n = auc(naive_full)
    auc_skill_post = None  # will use aggregate_metrics value
    auc_skill_ablation = auc(skill_full)

    # Read the original docked-only AUC
    agg = SRC_ROOT / "06_meta_analysis" / "aggregate_metrics.csv"
    for r in csv.DictReader(agg.open()):
        if r["target"] == target and r["protocol"] == "skill":
            auc_skill_post = float(r["roc_auc"])

    # DeLong on intersection of compounds in both naive_full and skill_full+ablation
    common_ids = set(naive_full) & set(skill_full)
    common = sorted(common_ids)
    y = np.array([1 if labels.get(c) == "active" else 0 for c in common])
    sn = np.array([-naive_full[c] for c in common])
    sk = np.array([-skill_full[c] for c in common])

    delong_p = delong_delta = delong_se = None
    try:
        from scipy.stats import norm
        delong_delta, delong_se, delong_p = delong_paired(y, sn, sk)
    except Exception as e:
        print(f"[{target}] DeLong failed: {e}")

    return {
        "target": target,
        "auc_naive_full": round(auc_n, 4) if not math.isnan(auc_n) else None,
        "auc_skill_post_filter": round(auc_skill_post, 4) if auc_skill_post else None,
        "auc_skill_ablation_full": round(auc_skill_ablation, 4) if not math.isnan(auc_skill_ablation) else None,
        "delta_ablation_minus_naive": round(auc_skill_ablation - auc_n, 4) if not (math.isnan(auc_skill_ablation) or math.isnan(auc_n)) else None,
        "delta_post_minus_naive": round(auc_skill_post - auc_n, 4) if (auc_skill_post and not math.isnan(auc_n)) else None,
        "n_compounds_skill_ablation": len(skill_full),
        "n_added_from_ablation": n_added,
        "intersection_n": len(common),
        "delong_delta": round(delong_delta, 4) if delong_delta is not None else None,
        "delong_se": round(delong_se, 6) if delong_se is not None else None,
        "delong_p": round(delong_p, 4) if delong_p is not None else None,
    }


# --- DeLong's test (paired) ---
def delong_paired(y: "np.ndarray", scores_a: "np.ndarray", scores_b: "np.ndarray") -> tuple[float, float, float]:
    """Paired DeLong test for two correlated AUCs.
    Returns (auc_b - auc_a, SE of difference, two-sided p-value).
    """
    import numpy as np
    from scipy.stats import norm
    pos = y == 1
    neg = ~pos
    m, n = pos.sum(), neg.sum()
    if m < 2 or n < 2:
        raise ValueError("need >=2 positives and >=2 negatives")

    def structural(scores: "np.ndarray") -> tuple["np.ndarray", "np.ndarray", float]:
        sp = scores[pos]
        sn = scores[neg]
        # V10 = P(Y_pos > X_neg | Y_pos), V01 = P(Y_pos > X_neg | X_neg)
        v10 = np.array([np.mean((sp[i] > sn) + 0.5 * (sp[i] == sn)) for i in range(m)])
        v01 = np.array([np.mean((sp > sn[j]) + 0.5 * (sp == sn[j])) for j in range(n)])
        auc = v10.mean()
        return v10, v01, auc

    v10a, v01a, auc_a = structural(scores_a)
    v10b, v01b, auc_b = structural(scores_b)
    # Covariance matrices
    s10 = np.cov(np.vstack([v10a, v10b]))
    s01 = np.cov(np.vstack([v01a, v01b]))
    s = s10 / m + s01 / n
    diff = auc_b - auc_a
    var_diff = s[0, 0] + s[1, 1] - 2 * s[0, 1]
    se = math.sqrt(var_diff) if var_diff > 0 else 0.0
    z = diff / se if se > 0 else 0.0
    p = 2 * (1 - norm.cdf(abs(z)))
    return diff, se, p


# --- Main ---
def main() -> int:
    args = parse_args()
    targets = ALL_TARGETS if args.targets == ["all"] else args.targets
    print(f"=== True filter ablation ===")
    print(f"Targets: {targets}")
    print(f"GPU: {args.gpu}")
    print(f"Max workers (prep): {args.max_workers or 'auto'}")
    print(f"Output root: {ABL_ROOT}")
    print()

    results = []
    for target in targets:
        print(f"\n{'='*60}\n{target.upper()}\n{'='*60}")
        try:
            if not args.skip_prep:
                prep_target(target, args.max_workers)
            if not args.skip_dock:
                dock_target(target, args.gpu)
            scores = parse_scores(target)
            r = compute_full_auc(target, scores)
            results.append(r)
            print(f"[{target}] result: {json.dumps(r, indent=2)}")
        except Exception as e:
            print(f"[{target}] FAILED: {e}")
            import traceback
            traceback.print_exc()
            results.append({"target": target, "error": str(e)})

    out = ABL_ROOT / "true_ablation_results.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\n=== DONE ===\nResults → {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
