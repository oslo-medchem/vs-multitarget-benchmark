#!/usr/bin/env python3
"""
Single-target orchestrator for the multi-target virtual screening benchmark.

Reads target_config.json and runs the full pipeline (stages 1-5) for one target.

Usage:
    python run_target_pipeline.py --target egfr
    python run_target_pipeline.py --target egfr --stages 1,2,3 --skip-existing
    python run_target_pipeline.py --target cdk2 --stages 4 --gpu 0
"""

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def banner(msg: str) -> None:
    """Print a stage banner."""
    width = max(len(msg) + 4, 60)
    print("\n" + "=" * width)
    print(f"  {msg}")
    print("=" * width, flush=True)


def run(cmd: str, env: dict | None = None) -> None:
    """Run a shell command with check=True, printing the command first."""
    print(f"  >> {cmd}", flush=True)
    subprocess.run(cmd, shell=True, check=True, env=env)


def elapsed(start: float) -> str:
    """Return a human-readable elapsed time string."""
    dt = time.time() - start
    if dt < 60:
        return f"{dt:.1f}s"
    elif dt < 3600:
        return f"{dt / 60:.1f}min"
    else:
        return f"{dt / 3600:.1f}h"


def file_exists_nonempty(path: Path) -> bool:
    """Check if a file exists and is non-empty."""
    return path.exists() and path.stat().st_size > 0


# ---------------------------------------------------------------------------
# Skip-existing checks: key output files per stage
# ---------------------------------------------------------------------------

def stage_complete(stage: int, target_dir: Path) -> bool:
    """Return True if the key output files for a stage already exist."""
    checks = {
        1: [
            target_dir / "01_active_curation" / "actives_100.smi",
            target_dir / "01_active_curation" / "actives_curated.csv",
        ],
        2: [
            target_dir / "02_decoy_generation" / "decoys_local_5000.smi",
        ],
        3: [
            target_dir / "03_library_preparation" / "library_combined.smi",
            target_dir / "03_library_preparation" / "library_labels.csv",
            target_dir / "03_library_preparation" / "naive" / "pdbqt",
            target_dir / "03_library_preparation" / "skill" / "pdbqt",
        ],
        4: [
            target_dir / "04_docking" / "scores_naive.csv",
            target_dir / "04_docking" / "scores_skill.csv",
        ],
        5: [
            target_dir / "05_evaluation" / "results",
        ],
    }
    paths = checks.get(stage, [])
    return all(p.exists() for p in paths)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the VS benchmark pipeline for a single target."
    )
    parser.add_argument(
        "--target", required=True,
        help="Target name, e.g. 'egfr' (must match a name in target_config.json)"
    )
    parser.add_argument(
        "--stages", default="1,2,3,4,5",
        help="Comma-separated list of stages to run (default: 1,2,3,4,5)"
    )
    parser.add_argument(
        "--config", default="target_config.json",
        help="Path to target_config.json (default: target_config.json)"
    )
    parser.add_argument(
        "--gpu", type=int, default=None,
        help="GPU device ID for docking stage (optional)"
    )
    parser.add_argument(
        "--skip-existing", action="store_true",
        help="Skip stages whose key output files already exist"
    )
    args = parser.parse_args()

    # Resolve paths relative to this script's directory
    BASE = Path(__file__).resolve().parent
    COMMON = BASE / "common"

    # If config path is relative, resolve relative to BASE
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = BASE / config_path

    # Load configuration
    with open(config_path) as f:
        config = json.load(f)

    global_cfg = config["global"]

    # Find the requested target
    target_cfg = None
    for t in config["targets"]:
        if t["name"] == args.target:
            target_cfg = t
            break
    if target_cfg is None:
        available = [t["name"] for t in config["targets"]]
        print(f"ERROR: Target '{args.target}' not found in config. "
              f"Available: {', '.join(available)}", file=sys.stderr)
        sys.exit(1)

    TARGET_DIR = BASE / "targets" / args.target
    stages = set(int(s.strip()) for s in args.stages.split(",") if s.strip())

    print(f"Target:     {args.target} ({target_cfg['description']})")
    print(f"ChEMBL ID:  {target_cfg['chembl_id']}")
    print(f"PDB:        {target_cfg['pdb_id']} chain {target_cfg['chain']}")
    print(f"Target dir: {TARGET_DIR}")
    print(f"Stages:     {sorted(stages)}")
    print(f"Skip existing: {args.skip_existing}")

    pipeline_start = time.time()

    # ------------------------------------------------------------------
    # Stage 1: Active curation
    # ------------------------------------------------------------------
    if 1 in stages:
        if args.skip_existing and stage_complete(1, TARGET_DIR):
            print("\n[Stage 1] SKIPPED (output files already exist)")
        else:
            banner(f"Stage 1: Active Curation — {args.target}")
            t0 = time.time()
            out_dir = TARGET_DIR / "01_active_curation"
            out_dir.mkdir(parents=True, exist_ok=True)

            run(
                f"python \"{COMMON}/fetch_actives.py\" "
                f"--chembl-id {target_cfg['chembl_id']} "
                f"--output-dir \"{out_dir}\" "
                f"--pchembl-min {global_cfg['pchembl_min']} "
                f"--target-name {target_cfg['name']}"
            )
            run(
                f"python \"{COMMON}/curate_actives.py\" "
                f"--input-csv \"{out_dir}/actives_raw.csv\" "
                f"--output-dir \"{out_dir}\" "
                f"--n-actives {global_cfg['n_actives_target']} "
                f"--target-name {target_cfg['name']} "
                f"--butina-cutoff {global_cfg['butina_cutoff']} "
                f"--mw-min {global_cfg['mw_min']} "
                f"--mw-max {global_cfg['mw_max']} "
                f"--max-heavy-atoms {global_cfg['max_heavy_atoms']}"
            )
            print(f"  Stage 1 done in {elapsed(t0)}")

    # ------------------------------------------------------------------
    # Stage 2: Decoy generation
    # ------------------------------------------------------------------
    if 2 in stages:
        if args.skip_existing and stage_complete(2, TARGET_DIR):
            print("\n[Stage 2] SKIPPED (output files already exist)")
        else:
            banner(f"Stage 2: Decoy Generation — {args.target}")
            t0 = time.time()
            out_dir = TARGET_DIR / "02_decoy_generation"
            out_dir.mkdir(parents=True, exist_ok=True)

            pool = BASE / global_cfg["pool_smi"]
            run(
                f"python \"{COMMON}/generate_decoys.py\" "
                f"--actives-csv \"{TARGET_DIR}/01_active_curation/actives_curated.csv\" "
                f"--output-dir \"{out_dir}\" "
                f"--pool-file \"{pool}\" "
                f"--n-decoys-per-active {global_cfg['decoy_ratio']} "
                f"--max-tc {global_cfg['max_tc_vs_actives']} "
                f"--seed {global_cfg['random_seed']}"
            )
            print(f"  Stage 2 done in {elapsed(t0)}")

    # ------------------------------------------------------------------
    # Stage 3: Library preparation
    # ------------------------------------------------------------------
    if 3 in stages:
        if args.skip_existing and stage_complete(3, TARGET_DIR):
            print("\n[Stage 3] SKIPPED (output files already exist)")
        else:
            banner(f"Stage 3: Library Preparation — {args.target}")
            t0 = time.time()
            lib_dir = TARGET_DIR / "03_library_preparation"
            lib_dir.mkdir(parents=True, exist_ok=True)

            # 3a: Prepare receptor (downloads PDB, extracts chain, computes box center)
            run(
                f"python \"{COMMON}/prepare_receptor.py\" "
                f"--pdb-id {target_cfg['pdb_id']} "
                f"--chain {target_cfg['chain']} "
                f"--ligand-resname {target_cfg['ligand_resname']} "
                f"--target-dir \"{TARGET_DIR}\""
            )

            # 3b: Combine library
            run(
                f"python \"{COMMON}/prepare_combined_library.py\" "
                f"--actives-smi \"{TARGET_DIR}/01_active_curation/actives_100.smi\" "
                f"--decoys-smi \"{TARGET_DIR}/02_decoy_generation/decoys_local_5000.smi\" "
                f"--output-dir \"{lib_dir}\" "
                f"--seed {global_cfg['random_seed']}"
            )

            # 3c: Prepare naive PDBQT
            run(
                f"python \"{COMMON}/prepare_naive.py\" "
                f"--library-smi \"{lib_dir}/library_combined.smi\" "
                f"--output-dir \"{lib_dir}/naive/pdbqt\""
            )

            # 3d: Prepare skill PDBQT
            run(
                f"python \"{COMMON}/prepare_skill.py\" "
                f"--library-csv \"{lib_dir}/library_labels.csv\" "
                f"--output-dir \"{lib_dir}/skill/pdbqt\""
            )

            # 3e: Generate docking scripts
            # Determine box center: box_center.json > config > error
            box_center_file = TARGET_DIR / "box_center.json"
            if box_center_file.exists():
                with open(box_center_file) as f:
                    bc = json.load(f)
                cx, cy, cz = bc["center_x"], bc["center_y"], bc["center_z"]
            elif target_cfg.get("box_center"):
                cx, cy, cz = target_cfg["box_center"]
            else:
                raise RuntimeError(
                    f"No box center available for {args.target}. "
                    f"Expected box_center.json in {TARGET_DIR} or "
                    f"'box_center' key in target_config.json."
                )

            run(
                f"python \"{COMMON}/generate_docking_scripts.py\" "
                f"--target-dir \"{TARGET_DIR}\" "
                f"--center-x {cx} "
                f"--center-y {cy} "
                f"--center-z {cz} "
                f"--naive-box-size {global_cfg['naive_box_size']} "
                f"--skill-box-size {global_cfg['skill_box_size']} "
                f"--naive-exhaustiveness {global_cfg['naive_exhaustiveness']} "
                f"--skill-exhaustiveness {global_cfg['skill_exhaustiveness']} "
                f"--num-modes {global_cfg['num_modes']}"
            )
            print(f"  Stage 3 done in {elapsed(t0)}")

    # ------------------------------------------------------------------
    # Stage 4: Docking
    # ------------------------------------------------------------------
    if 4 in stages:
        if args.skip_existing and stage_complete(4, TARGET_DIR):
            print("\n[Stage 4] SKIPPED (output files already exist)")
        else:
            banner(f"Stage 4: Docking — {args.target}")
            t0 = time.time()
            dock_dir = TARGET_DIR / "04_docking"

            env = os.environ.copy()
            if args.gpu is not None:
                env["CUDA_VISIBLE_DEVICES"] = str(args.gpu)

            # Run naive docking
            run(f"bash \"{dock_dir}/naive/run_unidock.sh\"", env=env)

            # Run skill docking
            run(f"bash \"{dock_dir}/skill/run_unidock.sh\"", env=env)

            # Parse scores
            run(f"python \"{COMMON}/parse_scores.py\" --target-dir \"{TARGET_DIR}\"")
            print(f"  Stage 4 done in {elapsed(t0)}")

    # ------------------------------------------------------------------
    # Stage 5: Evaluation
    # ------------------------------------------------------------------
    if 5 in stages:
        if args.skip_existing and stage_complete(5, TARGET_DIR):
            print("\n[Stage 5] SKIPPED (output files already exist)")
        else:
            banner(f"Stage 5: Evaluation — {args.target}")
            t0 = time.time()

            run(
                f"python \"{COMMON}/compute_metrics.py\" "
                f"--target-dir \"{TARGET_DIR}\""
            )
            run(
                f"python \"{COMMON}/statistical_analysis.py\" "
                f"--target-dir \"{TARGET_DIR}\""
            )
            run(
                f"python \"{COMMON}/generate_figures.py\" "
                f"--target-dir \"{TARGET_DIR}\" "
                f"--target-name {target_cfg['name']}"
            )
            print(f"  Stage 5 done in {elapsed(t0)}")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    banner(f"Pipeline complete for {args.target}")
    print(f"  Total wall time: {elapsed(pipeline_start)}")
    print(f"  Output directory: {TARGET_DIR}")


if __name__ == "__main__":
    main()
