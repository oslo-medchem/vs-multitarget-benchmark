#!/usr/bin/env python3
"""Fetch active compounds from ChEMBL for any target.

Queries bioactivity data with configurable pChEMBL threshold, saves raw CSV.
All parameters are supplied via command-line arguments so this script can be
called by a pipeline orchestrator for multiple targets.
"""

import argparse
import csv
import sys
from pathlib import Path


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Fetch active compounds from ChEMBL for a given target."
    )
    parser.add_argument(
        "--chembl-id",
        required=True,
        help="ChEMBL target ID, e.g. CHEMBL203",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for actives_raw.csv output",
    )
    parser.add_argument(
        "--pchembl-min",
        type=float,
        default=5.0,
        help="Minimum pChEMBL value filter (default: 5.0)",
    )
    parser.add_argument(
        "--target-name",
        default=None,
        help="Human-readable target name for logging (default: use chembl-id)",
    )
    return parser.parse_args(argv)


def fetch_activities(chembl_id, pchembl_min):
    """Query ChEMBL activities API and return list of record dicts."""
    from chembl_webresource_client.new_client import new_client

    activity = new_client.activity

    results = activity.filter(
        target_chembl_id=chembl_id,
        pchembl_value__gte=pchembl_min,
        standard_relation="=",
        standard_type__in=["EC50", "IC50", "Ki", "Kd"],
    )

    rows = []
    for r in results:
        if r.get("canonical_smiles") and r.get("pchembl_value"):
            rows.append(
                {
                    "molecule_chembl_id": r["molecule_chembl_id"],
                    "canonical_smiles": r["canonical_smiles"],
                    "standard_type": r["standard_type"],
                    "standard_value": r["standard_value"],
                    "standard_units": r["standard_units"],
                    "pchembl_value": float(r["pchembl_value"]),
                    "assay_chembl_id": r["assay_chembl_id"],
                    "assay_type": r.get("assay_type", ""),
                    "assay_description": r.get("assay_description", ""),
                    "document_chembl_id": r.get("document_chembl_id", ""),
                    "document_year": r.get("document_year", ""),
                }
            )

    return rows


def print_summary(rows):
    """Print summary statistics for fetched records."""
    unique_mols = len(set(r["molecule_chembl_id"] for r in rows))
    pchembl_vals = [r["pchembl_value"] for r in rows]
    types = {}
    for r in rows:
        types[r["standard_type"]] = types.get(r["standard_type"], 0) + 1

    print(f"\nSummary:")
    print(f"  Total records: {len(rows)}")
    print(f"  Unique molecules: {unique_mols}")
    print(f"  pChEMBL range: {min(pchembl_vals):.2f} - {max(pchembl_vals):.2f}")
    print(f"  Activity types: {types}")


def main(argv=None):
    args = parse_args(argv)

    target_name = args.target_name or args.chembl_id

    # Ensure output directory exists
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_path = args.output_dir / "actives_raw.csv"

    print(
        f"Querying ChEMBL for {target_name} ({args.chembl_id}) "
        f"activities (pChEMBL >= {args.pchembl_min})..."
    )

    rows = fetch_activities(args.chembl_id, args.pchembl_min)

    print(f"Retrieved {len(rows)} activity records")

    if not rows:
        print(
            f"ERROR: No activities found for {args.chembl_id} "
            f"with pChEMBL >= {args.pchembl_min}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Save to CSV
    fieldnames = list(rows[0].keys())
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved to {output_path}")
    print_summary(rows)

    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
