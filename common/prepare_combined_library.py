#!/usr/bin/env python3
"""
Merge actives + decoys into a single shuffled library with ground truth labels.

Usage:
    python prepare_combined_library.py \
        --actives-smi 01_active_curation/actives_100.smi \
        --decoys-smi 02_decoy_generation/decoys_local_5000.smi \
        --output-dir 03_library_preparation/ \
        --seed 42
"""

import argparse
import csv
import random
from pathlib import Path


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Merge actives and decoys into a shuffled library with ground truth labels."
    )
    p.add_argument(
        "--actives-smi",
        required=True,
        type=Path,
        help="Path to actives SMILES file (tab-separated: SMILES, compound_id [, source_id])",
    )
    p.add_argument(
        "--decoys-smi",
        required=True,
        type=Path,
        help="Path to decoys SMILES file (tab-separated: SMILES, compound_id [, source_id])",
    )
    p.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Output directory for library_combined.smi and library_labels.csv",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for shuffling (default: 42)",
    )
    return p.parse_args(argv)


def load_smi(path: Path, label: str) -> list[dict]:
    """Load a tab-separated SMILES file into a list of dicts."""
    entries = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            entries.append(
                {
                    "smiles": parts[0],
                    "compound_id": parts[1],
                    "label": label,
                    "source_id": parts[2] if len(parts) > 2 else "",
                }
            )
    return entries


def main(argv=None):
    args = parse_args(argv)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    random.seed(args.seed)

    # Load actives and decoys
    actives = load_smi(args.actives_smi, "active")
    print(f"Loaded {len(actives)} actives from {args.actives_smi}")

    decoys = load_smi(args.decoys_smi, "decoy")
    print(f"Loaded {len(decoys)} decoys from {args.decoys_smi}")

    # Combine and shuffle
    library = actives + decoys
    random.shuffle(library)
    print(f"Combined library: {len(library)} compounds (shuffled, seed={args.seed})")

    # Save labels CSV
    labels_csv = output_dir / "library_labels.csv"
    with open(labels_csv, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["compound_id", "label", "smiles", "source_id"]
        )
        writer.writeheader()
        for entry in library:
            writer.writerow(entry)
    print(f"Labels saved to {labels_csv}")

    # Save combined SMILES
    library_smi = output_dir / "library_combined.smi"
    with open(library_smi, "w") as f:
        for entry in library:
            f.write(f"{entry['smiles']}\t{entry['compound_id']}\n")
    print(f"Library SMILES saved to {library_smi}")

    # Summary
    n_act = sum(1 for e in library if e["label"] == "active")
    n_dec = sum(1 for e in library if e["label"] == "decoy")
    print(f"\nSummary: {n_act} actives + {n_dec} decoys = {len(library)} total")
    print(f"Active ratio: {n_act / len(library) * 100:.1f}%")


if __name__ == "__main__":
    main()
