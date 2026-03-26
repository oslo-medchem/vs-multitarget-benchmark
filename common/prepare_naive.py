#!/usr/bin/env python3
"""
Naive protocol: Prepare library compounds using OpenBabel defaults.

OpenBabel gen3d -> SDF, then SDF -> PDBQT with Gasteiger charges.
No filtering, no salt removal, no standardization beyond OpenBabel defaults.

Usage:
    python prepare_naive.py \
        --library-smi 03_library_preparation/library_combined.smi \
        --output-dir 03_library_preparation/naive/pdbqt/
"""

import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

TIMEOUT = 60  # seconds per compound


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Prepare ligand PDBQT files using OpenBabel naive protocol."
    )
    p.add_argument(
        "--library-smi",
        required=True,
        type=Path,
        help="Path to library_combined.smi (tab-separated: SMILES, compound_id)",
    )
    p.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Output directory for PDBQT files and preparation_log.csv",
    )
    return p.parse_args(argv)


def prepare_naive(smiles: str, compound_id: str, output_dir: Path) -> tuple[bool, str]:
    """Prepare a single compound using OpenBabel naive protocol.

    Returns (success, message).
    """
    output_pdbqt = output_dir / f"{compound_id}.pdbqt"

    smi_path = None
    sdf_path = None
    try:
        with tempfile.NamedTemporaryFile(
            suffix=".smi", mode="w", delete=False
        ) as smi_f:
            smi_f.write(f"{smiles}\n")
            smi_path = smi_f.name

        with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as sdf_f:
            sdf_path = sdf_f.name

        # Step 1: SMILES -> 3D SDF via gen3d
        result = subprocess.run(
            ["obabel", smi_path, "-O", sdf_path, "--gen3d", "-h"],
            capture_output=True,
            text=True,
            timeout=TIMEOUT,
        )
        if result.returncode != 0 or not Path(sdf_path).exists():
            return False, f"gen3d failed: {result.stderr[:200]}"

        if Path(sdf_path).stat().st_size < 50:
            return False, "Empty SDF"

        # Step 2: SDF -> PDBQT with Gasteiger charges
        result = subprocess.run(
            [
                "obabel", sdf_path, "-O", str(output_pdbqt),
                "-h", "--partialcharge", "gasteiger",
            ],
            capture_output=True,
            text=True,
            timeout=TIMEOUT,
        )
        if result.returncode != 0 or not output_pdbqt.exists():
            return False, f"PDBQT conversion failed: {result.stderr[:200]}"

        if output_pdbqt.stat().st_size < 50:
            return False, "Empty PDBQT"

        return True, "OK"

    except subprocess.TimeoutExpired:
        return False, "Timeout"
    except Exception as e:
        return False, str(e)[:200]
    finally:
        for p in [smi_path, sdf_path]:
            if p is not None:
                try:
                    Path(p).unlink(missing_ok=True)
                except Exception:
                    pass


def main(argv=None):
    args = parse_args(argv)
    output_dir = args.output_dir.resolve()
    pdbqt_dir = output_dir
    pdbqt_dir.mkdir(parents=True, exist_ok=True)

    # Load library
    compounds = []
    with open(args.library_smi) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                compounds.append((parts[0], parts[1]))

    print(f"Preparing {len(compounds)} compounds (naive protocol)")
    print(f"Output: {pdbqt_dir}")

    # Process sequentially
    results = []
    success = 0
    fail = 0
    for i, (smiles, cid) in enumerate(compounds):
        ok, msg = prepare_naive(smiles, cid, pdbqt_dir)
        results.append({"compound_id": cid, "success": ok, "message": msg})
        if ok:
            success += 1
        else:
            fail += 1

        if (i + 1) % 100 == 0:
            print(f"  [{i + 1}/{len(compounds)}] Success: {success}, Failed: {fail}")

    # Save log
    log_file = output_dir / "preparation_log.csv"
    with open(log_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["compound_id", "success", "message"])
        writer.writeheader()
        writer.writerows(results)

    rate = success / len(compounds) * 100 if compounds else 0
    print(f"\nDone. Success: {success}/{len(compounds)} ({rate:.1f}%)")
    print(f"Failed: {fail}")
    print(f"PDBQT files in: {pdbqt_dir}")
    print(f"Log: {log_file}")


if __name__ == "__main__":
    main()
