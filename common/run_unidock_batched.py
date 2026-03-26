#!/usr/bin/env python3
"""
Run Uni-Dock in batches to handle PDBQT parsing errors gracefully.
Uni-Dock can crash on malformed PDBQTs; this script splits the ligand list
into batches and retries failed batches with smaller sizes.
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path


def run_unidock(receptor, ligand_list, center, size, exhaustiveness,
                num_modes, results_dir, gpu=None):
    """Run Uni-Dock and return (success, n_results)."""
    env = os.environ.copy()
    if gpu is not None:
        env["CUDA_VISIBLE_DEVICES"] = str(gpu)

    cmd = [
        "unidock",
        "--receptor", str(receptor),
        "--ligand_index", str(ligand_list),
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size),
        "--size_y", str(size),
        "--size_z", str(size),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--dir", str(results_dir),
        "--scoring", "vina",
    ]

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=3600, env=env
        )
        # Count results
        n_results = len(list(Path(results_dir).glob("*_out.pdbqt")))
        has_batch = "Batch" in result.stdout or "Batch" in result.stderr
        return has_batch, n_results
    except subprocess.TimeoutExpired:
        return False, 0
    except Exception as e:
        print(f"  Error: {e}")
        return False, 0


def main():
    parser = argparse.ArgumentParser(description="Batched Uni-Dock runner")
    parser.add_argument("--receptor", required=True)
    parser.add_argument("--ligand-list", required=True)
    parser.add_argument("--center-x", required=True, type=float)
    parser.add_argument("--center-y", required=True, type=float)
    parser.add_argument("--center-z", required=True, type=float)
    parser.add_argument("--size", required=True, type=int)
    parser.add_argument("--exhaustiveness", required=True, type=int)
    parser.add_argument("--num-modes", type=int, default=9)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--gpu", type=int, default=None)
    parser.add_argument("--batch-size", type=int, default=100)
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    # Load ligand list
    ligands = [l.strip() for l in open(args.ligand_list) if l.strip()]
    print(f"Total ligands: {len(ligands)}")

    center = (args.center_x, args.center_y, args.center_z)
    total_docked = 0
    total_failed_batches = 0

    # Count existing results
    existing = len(list(results_dir.glob("*_out.pdbqt")))
    if existing > 0:
        print(f"Existing results: {existing}")

    # Process in batches
    batch_size = args.batch_size
    for i in range(0, len(ligands), batch_size):
        batch = ligands[i:i + batch_size]
        batch_num = i // batch_size + 1

        # Write temp ligand list
        tmp_list = results_dir / f"_batch_{batch_num}.txt"
        with open(tmp_list, "w") as f:
            f.write("\n".join(batch) + "\n")

        n_before = len(list(results_dir.glob("*_out.pdbqt")))

        print(f"  Batch {batch_num} ({len(batch)} ligands, {i+1}-{i+len(batch)})...",
              end=" ", flush=True)

        ok, n_new = run_unidock(
            args.receptor, tmp_list, center, args.size,
            args.exhaustiveness, args.num_modes, results_dir, args.gpu
        )

        n_after = len(list(results_dir.glob("*_out.pdbqt")))
        docked_this_batch = n_after - n_before

        if ok:
            print(f"OK ({docked_this_batch} docked)")
            total_docked += docked_this_batch
        else:
            # Try smaller sub-batches
            sub_size = max(10, batch_size // 5)
            print(f"FAILED, retrying in sub-batches of {sub_size}")
            for j in range(0, len(batch), sub_size):
                sub_batch = batch[j:j + sub_size]
                sub_list = results_dir / f"_sub_{batch_num}_{j}.txt"
                with open(sub_list, "w") as f:
                    f.write("\n".join(sub_batch) + "\n")

                n_before2 = len(list(results_dir.glob("*_out.pdbqt")))
                ok2, _ = run_unidock(
                    args.receptor, sub_list, center, args.size,
                    args.exhaustiveness, args.num_modes, results_dir, args.gpu
                )
                n_after2 = len(list(results_dir.glob("*_out.pdbqt")))
                sub_docked = n_after2 - n_before2
                total_docked += sub_docked

                if not ok2:
                    total_failed_batches += 1

                sub_list.unlink(missing_ok=True)

        tmp_list.unlink(missing_ok=True)

    # Clean up temp files
    for f in results_dir.glob("_batch_*.txt"):
        f.unlink(missing_ok=True)
    for f in results_dir.glob("_sub_*.txt"):
        f.unlink(missing_ok=True)

    total_results = len(list(results_dir.glob("*_out.pdbqt")))
    print(f"\nDone. Total docked: {total_results}")
    print(f"Failed batches: {total_failed_batches}")
    print(f"Success rate: {total_results}/{len(ligands)} ({100*total_results/max(1,len(ligands)):.1f}%)")


if __name__ == "__main__":
    main()
