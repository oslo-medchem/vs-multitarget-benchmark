#!/usr/bin/env python3
"""Generate property-matched decoys for virtual screening benchmarks.

Loads a pre-existing ChEMBL compound pool, property-matches against curated
actives, and applies Tanimoto similarity filtering to ensure structural
dissimilarity.  Falls back to widened tolerances when too few decoys are found.

All parameters are supplied via command-line arguments so this script can be
called by a pipeline orchestrator for multiple targets.
"""

import argparse
import csv
import random
import sys
import numpy as np
from pathlib import Path
from collections import defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import TanimotoSimilarity

RDLogger.logger().setLevel(RDLogger.ERROR)

# Default property-match tolerances (tight pass)
TOL_MW = 25
TOL_LOGP = 1.0
TOL_HBD = 1
TOL_HBA = 2
TOL_ROTBONDS = 2

# Widened tolerances (fallback pass)
TOL_MW_WIDE = 50
TOL_LOGP_WIDE = 2.0
TOL_HBD_WIDE = 2
TOL_HBA_WIDE = 3
TOL_ROTBONDS_WIDE = 3


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate property-matched decoys from a pre-existing compound pool."
    )
    parser.add_argument(
        "--actives-csv",
        required=True,
        type=Path,
        help="Path to actives_curated.csv from curate_actives.py",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for decoys_local_5000.smi, decoys_curated.csv",
    )
    parser.add_argument(
        "--pool-file",
        required=True,
        type=Path,
        help="Path to shared ChEMBL pool SMILES file (tab-separated: SMILES<tab>CHEMBL_ID)",
    )
    parser.add_argument(
        "--n-decoys-per-active",
        type=int,
        default=50,
        help="Target number of decoys per active (default: 50)",
    )
    parser.add_argument(
        "--max-tc",
        type=float,
        default=0.35,
        help="Maximum Tanimoto coefficient vs any active (default: 0.35)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)",
    )
    return parser.parse_args(argv)


def compute_properties(mol):
    """Compute drug-like molecular properties."""
    return {
        "MW": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "RotBonds": Descriptors.NumRotatableBonds(mol),
        "TPSA": Descriptors.TPSA(mol),
        "NumHeavyAtoms": mol.GetNumHeavyAtoms(),
        "NetCharge": Chem.GetFormalCharge(mol),
    }


def load_actives(actives_csv):
    """Load curated actives from CSV."""
    actives = []
    with open(actives_csv) as f:
        for row in csv.DictReader(f):
            mol = Chem.MolFromSmiles(row["std_smiles"])
            if mol is None:
                continue
            actives.append(
                {
                    "compound_id": row["compound_id"],
                    "smiles": row["std_smiles"],
                    "mol": mol,
                    "MW": float(row["MW"]),
                    "LogP": float(row["LogP"]),
                    "HBD": int(float(row["HBD"])),
                    "HBA": int(float(row["HBA"])),
                    "RotBonds": int(float(row["RotBonds"])),
                    "NetCharge": int(float(row["NetCharge"])),
                }
            )
    return actives


def load_pool(pool_file):
    """Load compound pool from tab-separated SMILES file."""
    pool = []
    with open(pool_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                pool.append({"smiles": parts[0], "chembl_id": parts[1]})
    return pool


def property_match(pool_props, active, tol_mw, tol_logp, tol_hbd, tol_hba, tol_rotbonds, check_charge=True):
    """Check whether pool compound properties match an active within tolerances."""
    if abs(pool_props["MW"] - active["MW"]) > tol_mw:
        return False
    if abs(pool_props["LogP"] - active["LogP"]) > tol_logp:
        return False
    if abs(pool_props["HBD"] - active["HBD"]) > tol_hbd:
        return False
    if abs(pool_props["HBA"] - active["HBA"]) > tol_hba:
        return False
    if abs(pool_props["RotBonds"] - active["RotBonds"]) > tol_rotbonds:
        return False
    if check_charge and pool_props["NetCharge"] != active["NetCharge"]:
        return False
    return True


def main(argv=None):
    args = parse_args(argv)

    # Validate inputs
    if not args.actives_csv.exists():
        print(f"ERROR: Actives file not found: {args.actives_csv}", file=sys.stderr)
        sys.exit(1)
    if not args.pool_file.exists():
        print(f"ERROR: Pool file not found: {args.pool_file}", file=sys.stderr)
        sys.exit(1)

    # Ensure output directory exists
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_smi = args.output_dir / "decoys_local_5000.smi"
    output_csv = args.output_dir / "decoys_curated.csv"

    random.seed(args.seed)
    np.random.seed(args.seed)

    # Load actives
    actives = load_actives(args.actives_csv)
    print(f"Loaded {len(actives)} actives from {args.actives_csv}")
    if not actives:
        print("ERROR: No valid actives loaded", file=sys.stderr)
        sys.exit(1)

    # Load pool
    raw_pool = load_pool(args.pool_file)
    print(f"Pool: {len(raw_pool)} molecules from {args.pool_file}")
    if not raw_pool:
        print("ERROR: Empty compound pool", file=sys.stderr)
        sys.exit(1)

    # Compute fingerprints for actives
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    active_fps = []
    for a in actives:
        a["fp"] = fpgen.GetFingerprint(a["mol"])
        active_fps.append(a["fp"])

    # Process pool: compute properties and fingerprints
    print("Computing properties and fingerprints for pool...")
    pool = []
    for i, entry in enumerate(raw_pool):
        mol = Chem.MolFromSmiles(entry["smiles"])
        if mol is None:
            continue
        props = compute_properties(mol)
        fp = fpgen.GetFingerprint(mol)
        pool.append(
            {
                "chembl_id": entry["chembl_id"],
                "smiles": entry["smiles"],
                "mol": mol,
                "fp": fp,
                "props": props,
            }
        )
        if (i + 1) % 10000 == 0:
            print(f"  Processed {i + 1}/{len(raw_pool)}...")

    print(f"Valid pool: {len(pool)} molecules")

    # Build InChIKey exclusion set from actives
    used_inchikeys = set()
    for a in actives:
        ik = Chem.MolToInchiKey(a["mol"])
        if ik:
            used_inchikeys.add(ik[:14])

    # Phase 1: Tight-tolerance matching
    print(f"\nPhase 1: Matching decoys (tight tolerances, max Tc < {args.max_tc})...")
    all_decoys = []
    decoys_per_active = defaultdict(int)

    for ai, active in enumerate(actives):
        needed = args.n_decoys_per_active

        candidates = []
        for p in pool:
            if not property_match(
                p["props"], active, TOL_MW, TOL_LOGP, TOL_HBD, TOL_HBA, TOL_ROTBONDS
            ):
                continue

            ik = Chem.MolToInchiKey(p["mol"])
            if ik and ik[:14] in used_inchikeys:
                continue

            max_tc = max(TanimotoSimilarity(p["fp"], afp) for afp in active_fps)
            if max_tc >= args.max_tc:
                continue

            candidates.append((p, ik, max_tc))

        # Sort by lowest Tc first (most dissimilar preferred)
        candidates.sort(key=lambda x: x[2])

        picked = 0
        for p, ik, max_tc in candidates:
            if picked >= needed:
                break
            decoy = {
                "matched_to": active["compound_id"],
                "chembl_id": p["chembl_id"],
                "smiles": p["smiles"],
                "max_tc_to_actives": max_tc,
            }
            decoy.update(p["props"])
            all_decoys.append(decoy)
            if ik:
                used_inchikeys.add(ik[:14])
            picked += 1

        decoys_per_active[active["compound_id"]] = picked
        if (ai + 1) % 20 == 0:
            print(f"  Active {ai + 1}/{len(actives)}: total decoys = {len(all_decoys)}")

    print(f"\nPhase 1 total decoys: {len(all_decoys)}")

    # Report under-matched actives
    insufficient = [(aid, n) for aid, n in decoys_per_active.items() if n < args.n_decoys_per_active]
    if insufficient:
        print(f"Actives with < {args.n_decoys_per_active} decoys: {len(insufficient)}")

    # Phase 2: Widened tolerances if too few decoys
    if len(all_decoys) < 3000:
        print("\nPhase 2: Widening tolerances for second pass...")
        for ai, active in enumerate(actives):
            needed = args.n_decoys_per_active - decoys_per_active[active["compound_id"]]
            if needed <= 0:
                continue

            for p in pool:
                if needed <= 0:
                    break
                if not property_match(
                    p["props"],
                    active,
                    TOL_MW_WIDE,
                    TOL_LOGP_WIDE,
                    TOL_HBD_WIDE,
                    TOL_HBA_WIDE,
                    TOL_ROTBONDS_WIDE,
                    check_charge=False,
                ):
                    continue

                ik = Chem.MolToInchiKey(p["mol"])
                if ik and ik[:14] in used_inchikeys:
                    continue

                max_tc = max(TanimotoSimilarity(p["fp"], afp) for afp in active_fps)
                if max_tc >= args.max_tc:
                    continue

                decoy = {
                    "matched_to": active["compound_id"],
                    "chembl_id": p["chembl_id"],
                    "smiles": p["smiles"],
                    "max_tc_to_actives": max_tc,
                }
                decoy.update(p["props"])
                all_decoys.append(decoy)
                if ik:
                    used_inchikeys.add(ik[:14])
                needed -= 1
                decoys_per_active[active["compound_id"]] += 1

        print(f"After Phase 2: {len(all_decoys)} decoys")

    # Trim to 5000 max
    max_decoys = args.n_decoys_per_active * len(actives)
    max_decoys = min(max_decoys, 5000)
    if len(all_decoys) > max_decoys:
        print(f"Trimming from {len(all_decoys)} to {max_decoys} decoys...")
        random.shuffle(all_decoys)
        all_decoys = all_decoys[:max_decoys]

    # Save SMILES file
    with open(output_smi, "w") as f:
        for i, d in enumerate(all_decoys, 1):
            did = f"DEC_{i:04d}"
            d["compound_id"] = did
            f.write(f"{d['smiles']}\t{did}\t{d['chembl_id']}\n")
    print(f"\nSaved {len(all_decoys)} decoys to {output_smi}")

    # Save CSV
    csv_fields = [
        "compound_id",
        "matched_to",
        "chembl_id",
        "smiles",
        "MW",
        "LogP",
        "HBD",
        "HBA",
        "RotBonds",
        "TPSA",
        "NumHeavyAtoms",
        "NetCharge",
        "max_tc_to_actives",
    ]
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_fields, extrasaction="ignore")
        writer.writeheader()
        for d in all_decoys:
            writer.writerow(d)
    print(f"Saved CSV to {output_csv}")

    # Summary statistics
    tc_vals = [d["max_tc_to_actives"] for d in all_decoys]
    print(
        f"\nMax Tc to actives: min={min(tc_vals):.3f}, "
        f"max={max(tc_vals):.3f}, mean={np.mean(tc_vals):.3f}"
    )

    matched_counts = list(decoys_per_active.values())
    print(
        f"Decoys per active: min={min(matched_counts)}, "
        f"max={max(matched_counts)}, mean={np.mean(matched_counts):.1f}"
    )

    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
