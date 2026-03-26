#!/usr/bin/env python3
"""
Prepare receptor PDBQT files for both naive and skill protocols.

Downloads a PDB structure from RCSB, extracts the specified chain,
computes the docking box center from a co-crystallised ligand, and
prepares receptor PDBQT files using OpenBabel (naive) and Meeko (skill).

Usage:
    python prepare_receptor.py \
        --pdb-id 1M17 --chain A --ligand-resname AQ4 \
        --target-dir targets/egfr/
"""

import argparse
import json
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Download PDB, extract chain, compute box center, prepare receptor PDBQT."
    )
    p.add_argument("--pdb-id", required=True, help="PDB accession code (e.g. 1M17)")
    p.add_argument("--chain", required=True, help="Chain identifier (e.g. A)")
    p.add_argument(
        "--ligand-resname",
        required=True,
        help="Residue name of co-crystallised ligand for box center (e.g. AQ4)",
    )
    p.add_argument(
        "--target-dir",
        required=True,
        type=Path,
        help="Target output directory (e.g. targets/egfr/)",
    )
    p.add_argument(
        "--work-dir",
        type=Path,
        default=None,
        help="Working directory for temp files (default: <target-dir>/00_setup)",
    )
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def download_pdb(pdb_id: str, work_dir: Path) -> Path:
    """Download PDB file from RCSB, caching if already present."""
    pdb_path = work_dir / f"{pdb_id.upper()}.pdb"
    if pdb_path.exists() and pdb_path.stat().st_size > 1000:
        print(f"Using cached {pdb_path}")
        return pdb_path
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    print(f"Downloading {url} ...")
    urllib.request.urlretrieve(url, pdb_path)
    print(f"Downloaded: {pdb_path.stat().st_size:,} bytes")
    return pdb_path


# ---------------------------------------------------------------------------
# Chain extraction and box center
# ---------------------------------------------------------------------------

def compute_box_center(pdb_path: Path, chain: str, ligand_resname: str) -> dict:
    """Compute centroid of all HETATM atoms matching ligand_resname on the given chain."""
    all_lines = open(pdb_path).readlines()
    xs, ys, zs = [], [], []
    for line in all_lines:
        if line[:6].strip() == "HETATM":
            res = line[17:20].strip()
            ch = line[21]
            if res == ligand_resname and ch == chain:
                try:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
                except ValueError:
                    continue
    if not xs:
        # Try all chains as fallback
        for line in all_lines:
            if line.startswith("HETATM"):
                resname = line[17:20].strip()
                if resname == ligand_resname:
                    try:
                        xs.append(float(line[30:38]))
                        ys.append(float(line[38:46]))
                        zs.append(float(line[46:54]))
                    except (ValueError, IndexError):
                        continue
        if xs:
            print(f"NOTE: Found {len(xs)} {ligand_resname} atoms on other chains (not {chain})")
    if not xs:
        available = sorted(set(
            l[17:20].strip() for l in all_lines
            if l.startswith("HETATM") and l[17:20].strip() != "HOH"
        ))
        raise RuntimeError(
            f"ERROR: No HETATM atoms found for residue '{ligand_resname}' in PDB. "
            f"Available HETATM residues: {', '.join(available)}"
        )

    center = {
        "center_x": round(sum(xs) / len(xs), 3),
        "center_y": round(sum(ys) / len(ys), 3),
        "center_z": round(sum(zs) / len(zs), 3),
    }
    print(
        f"Box center from {len(xs)} {ligand_resname} atoms: "
        f"({center['center_x']}, {center['center_y']}, {center['center_z']})"
    )
    return center


def extract_chain(pdb_path: Path, chain: str, work_dir: Path) -> Path:
    """Extract ATOM records for the specified chain (clean protein, no HETATM)."""
    pdb_id = pdb_path.stem
    chain_path = work_dir / f"{pdb_id}_chain{chain}.pdb"
    kept = []
    with open(pdb_path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec == "ATOM" and line[21] == chain:
                kept.append(line)
            elif rec == "TER":
                # Keep TER lines associated with our chain
                if kept and kept[-1][:6].strip() == "ATOM":
                    kept.append(line)
    kept.append("END\n")
    with open(chain_path, "w") as fh:
        fh.writelines(kept)
    n_atoms = sum(1 for l in kept if l[:6].strip() == "ATOM")
    print(f"Extracted chain {chain}: {n_atoms} ATOM records -> {chain_path}")
    return chain_path


# ---------------------------------------------------------------------------
# Receptor preparation
# ---------------------------------------------------------------------------

def prepare_naive_receptor(chain_path: Path, output_pdbqt: Path) -> bool:
    """OpenBabel: PDB -> PDBQT (add hydrogens, Gasteiger charges)."""
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)

    # Try with -xr (rigid receptor) first
    result = subprocess.run(
        [
            "obabel", str(chain_path), "-O", str(output_pdbqt),
            "-xr", "--partialcharge", "gasteiger",
        ],
        capture_output=True, text=True,
    )
    if output_pdbqt.exists() and output_pdbqt.stat().st_size > 1000:
        print(f"Naive receptor: {output_pdbqt} ({output_pdbqt.stat().st_size:,} bytes)")
        return True

    # Fallback without -xr
    if result.returncode != 0:
        print(f"obabel -xr failed ({result.stderr[:300]}), retrying without -xr ...")
    result2 = subprocess.run(
        [
            "obabel", str(chain_path), "-O", str(output_pdbqt),
            "--partialcharge", "gasteiger",
        ],
        capture_output=True, text=True,
    )
    if output_pdbqt.exists() and output_pdbqt.stat().st_size > 1000:
        print(f"Naive receptor (fallback): {output_pdbqt} ({output_pdbqt.stat().st_size:,} bytes)")
        return True

    print(f"ERROR: Naive receptor preparation failed. stderr: {result2.stderr[:300]}", file=sys.stderr)
    return False


def prepare_skill_receptor(chain_path: Path, output_pdbqt: Path) -> bool:
    """Meeko mk_prepare_receptor.py -> PDBQT, with obabel fallback."""
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)

    result = subprocess.run(
        [
            "mk_prepare_receptor.py", "--read_pdb", str(chain_path),
            "-p", str(output_pdbqt),
            "--default_altloc", "A", "-a",
        ],
        capture_output=True, text=True,
    )
    if result.returncode == 0 and output_pdbqt.exists() and output_pdbqt.stat().st_size > 1000:
        print(f"Skill receptor: {output_pdbqt} ({output_pdbqt.stat().st_size:,} bytes)")
        return True

    # Fallback to obabel (must use -xr for rigid receptor)
    print(f"Meeko failed ({result.stderr[:200]}), falling back to obabel ...")
    result2 = subprocess.run(
        [
            "obabel", str(chain_path), "-O", str(output_pdbqt),
            "-xr", "--partialcharge", "gasteiger",
        ],
        capture_output=True, text=True,
    )
    if output_pdbqt.exists() and output_pdbqt.stat().st_size > 1000:
        print(f"Skill receptor (obabel fallback): {output_pdbqt} ({output_pdbqt.stat().st_size:,} bytes)")
        return True

    print(f"ERROR: Skill receptor preparation failed. stderr: {result2.stderr[:300]}", file=sys.stderr)
    return False


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    args = parse_args(argv)

    pdb_id = args.pdb_id.upper()
    chain = args.chain
    ligand_resname = args.ligand_resname.upper()
    target_dir = args.target_dir.resolve()
    work_dir = (args.work_dir or target_dir / "00_setup").resolve()

    work_dir.mkdir(parents=True, exist_ok=True)
    target_dir.mkdir(parents=True, exist_ok=True)

    print(f"=== Receptor Preparation: {pdb_id} chain {chain} ===")
    print(f"  Ligand resname : {ligand_resname}")
    print(f"  Target dir     : {target_dir}")
    print(f"  Work dir       : {work_dir}")

    # 1. Download PDB
    pdb_path = download_pdb(pdb_id, work_dir)

    # 2. Compute box center from ligand HETATM coordinates
    center = compute_box_center(pdb_path, chain, ligand_resname)
    center_path = target_dir / "box_center.json"
    with open(center_path, "w") as fh:
        json.dump(center, fh, indent=2)
    print(f"Box center saved to {center_path}")

    # 3. Extract clean chain (ATOM only)
    chain_path = extract_chain(pdb_path, chain, work_dir)

    # 4. Prepare naive receptor
    naive_dir = target_dir / "04_docking" / "naive"
    naive_pdbqt = naive_dir / "receptor_naive.pdbqt"
    print(f"\n--- Naive protocol (OpenBabel) ---")
    ok_naive = prepare_naive_receptor(chain_path, naive_pdbqt)

    # 5. Prepare skill receptor
    skill_dir = target_dir / "04_docking" / "skill"
    skill_pdbqt = skill_dir / "receptor_skill.pdbqt"
    print(f"\n--- Skill protocol (Meeko) ---")
    ok_skill = prepare_skill_receptor(chain_path, skill_pdbqt)

    # Summary
    print(f"\n=== Summary ===")
    print(f"  Chain PDB      : {chain_path}")
    print(f"  Box center     : {center_path}")
    print(f"  Naive receptor : {naive_pdbqt} ({'OK' if ok_naive else 'FAILED'})")
    print(f"  Skill receptor : {skill_pdbqt} ({'OK' if ok_skill else 'FAILED'})")

    if not (ok_naive and ok_skill):
        sys.exit(1)


if __name__ == "__main__":
    main()
