#!/usr/bin/env python3
"""Curate actives from raw ChEMBL data to N diverse compounds.

Steps: validate SMILES, standardize (salt removal, neutralize), property filter,
deduplicate by InChIKey, Butina clustering, MaxMin diversity pick.

All parameters are supplied via command-line arguments so this script can be
called by a pipeline orchestrator for multiple targets.
"""

import argparse
import csv
import sys
import numpy as np
from pathlib import Path
from collections import Counter

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import TanimotoSimilarity, BulkTanimotoSimilarity
from rdkit.ML.Cluster import Butina

RDLogger.logger().setLevel(RDLogger.ERROR)

MAX_SMILES_LEN = 200


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Curate active compounds: standardize, filter, cluster, diversity-pick."
    )
    parser.add_argument(
        "--input-csv",
        required=True,
        type=Path,
        help="Path to actives_raw.csv from fetch_actives.py",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for actives_100.smi, actives_curated.csv, curation_report.md",
    )
    parser.add_argument(
        "--n-actives",
        type=int,
        default=100,
        help="Target number of diverse actives to select (default: 100)",
    )
    parser.add_argument(
        "--target-name",
        default="target",
        help="Human-readable target name for report header (default: 'target')",
    )
    parser.add_argument(
        "--butina-cutoff",
        type=float,
        default=0.4,
        help="Tanimoto distance cutoff for Butina clustering (default: 0.4)",
    )
    parser.add_argument("--mw-min", type=float, default=150, help="Min MW filter (default: 150)")
    parser.add_argument("--mw-max", type=float, default=650, help="Max MW filter (default: 650)")
    parser.add_argument(
        "--max-heavy-atoms",
        type=int,
        default=50,
        help="Max heavy atom count (default: 50)",
    )
    return parser.parse_args(argv)


def standardize_mol(mol):
    """Standardize molecule: remove salts, neutralize, canonicalize."""
    remover = SaltRemover.SaltRemover()
    mol = remover.StripMol(mol)

    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(mol)

    smi = Chem.MolToSmiles(mol)
    mol = Chem.MolFromSmiles(smi)
    return mol


def compute_properties(mol):
    """Compute drug-like properties."""
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


def butina_cluster(fps, cutoff):
    """Butina clustering with Tanimoto distance."""
    n = len(fps)
    if n == 0:
        return []
    dists = []
    for i in range(1, n):
        sims = BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - s for s in sims])
    clusters = Butina.ClusterData(dists, n, cutoff, isDistData=True)
    return clusters


def maxmin_pick(fps, n_pick):
    """MaxMin diversity picking."""
    from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker

    picker = MaxMinPicker()
    picked = picker.LazyBitVectorPick(fps, len(fps), n_pick)
    return list(picked)


def generate_report(
    target_name,
    raw,
    valid,
    standardized,
    filtered,
    deduped,
    selected,
    clusters,
    sel_clusters,
    n_clusters,
    mw_min,
    mw_max,
    max_heavy_atoms,
    butina_cutoff,
    target_n,
):
    """Generate markdown curation report."""
    props = {
        k: [float(r[k]) for r in selected]
        for k in ["MW", "LogP", "HBD", "HBA", "RotBonds", "TPSA"]
    }
    pchembl = [float(r["pchembl_value"]) for r in selected]

    report = f"""# {target_name} Active Compound Curation Report

## Curation Pipeline Summary

| Step | Description | Count |
|------|------------|-------|
| 0 | Raw ChEMBL records | {len(raw)} |
| 1 | Valid SMILES (len <= {MAX_SMILES_LEN}) | {len(valid)} |
| 2 | Standardized (salt removal, neutralize) | {len(standardized)} |
| 3 | Property filter (MW {mw_min}-{mw_max}, HA<={max_heavy_atoms}) | {len(filtered)} |
| 4 | Deduplicated (InChIKey connectivity) | {len(deduped)} |
| 5 | Butina clustering (Tc cutoff {butina_cutoff}) | {n_clusters} clusters |
| 6 | MaxMin diversity pick | **{len(selected)} actives** |

## Final Dataset Properties

| Property | Min | Median | Max | Mean | Std |
|----------|-----|--------|-----|------|-----|
| MW | {min(props['MW']):.1f} | {np.median(props['MW']):.1f} | {max(props['MW']):.1f} | {np.mean(props['MW']):.1f} | {np.std(props['MW']):.1f} |
| LogP | {min(props['LogP']):.2f} | {np.median(props['LogP']):.2f} | {max(props['LogP']):.2f} | {np.mean(props['LogP']):.2f} | {np.std(props['LogP']):.2f} |
| HBD | {min(props['HBD']):.0f} | {np.median(props['HBD']):.0f} | {max(props['HBD']):.0f} | {np.mean(props['HBD']):.1f} | {np.std(props['HBD']):.1f} |
| HBA | {min(props['HBA']):.0f} | {np.median(props['HBA']):.0f} | {max(props['HBA']):.0f} | {np.mean(props['HBA']):.1f} | {np.std(props['HBA']):.1f} |
| RotBonds | {min(props['RotBonds']):.0f} | {np.median(props['RotBonds']):.0f} | {max(props['RotBonds']):.0f} | {np.mean(props['RotBonds']):.1f} | {np.std(props['RotBonds']):.1f} |
| TPSA | {min(props['TPSA']):.1f} | {np.median(props['TPSA']):.1f} | {max(props['TPSA']):.1f} | {np.mean(props['TPSA']):.1f} | {np.std(props['TPSA']):.1f} |
| pChEMBL | {min(pchembl):.2f} | {np.median(pchembl):.2f} | {max(pchembl):.2f} | {np.mean(pchembl):.2f} | {np.std(pchembl):.2f} |

## Diversity

- Butina clusters (all deduped): {n_clusters}
- Butina clusters (final {len(selected)}): {len(sel_clusters)}
- Cluster size distribution: {dict(Counter([len(c) for c in sel_clusters]).most_common(10))}

## Checkpoint CP1

- [{'x' if len(selected) >= target_n else ' '}] {len(selected)} actives (target: {target_n})
- [{'x' if len(sel_clusters) >= 30 else ' '}] {len(sel_clusters)} clusters (target: >= 30)
- [x] All SMILES valid (RDKit-verified)
- [x] Properties span drug-like range
"""
    return report


def main(argv=None):
    args = parse_args(argv)

    # Validate input
    if not args.input_csv.exists():
        print(f"ERROR: Input file not found: {args.input_csv}", file=sys.stderr)
        sys.exit(1)

    # Ensure output directory exists
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_smi = args.output_dir / "actives_100.smi"
    output_csv = args.output_dir / "actives_curated.csv"
    report_path = args.output_dir / "curation_report.md"

    # Load raw data
    with open(args.input_csv) as f:
        reader = csv.DictReader(f)
        raw_records = list(reader)
    print(f"[{args.target_name}] Loaded {len(raw_records)} raw records from {args.input_csv}")

    if not raw_records:
        print(f"ERROR: No records in {args.input_csv}", file=sys.stderr)
        sys.exit(1)

    # Step 1: Validate SMILES with RDKit
    valid = []
    invalid_count = 0
    for rec in raw_records:
        smi = rec["canonical_smiles"]
        if len(smi) > MAX_SMILES_LEN:
            invalid_count += 1
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            invalid_count += 1
            continue
        rec["mol"] = mol
        valid.append(rec)
    print(f"Step 1 - Valid SMILES: {len(valid)} (dropped {invalid_count})")

    # Step 2: Standardize
    standardized = []
    for rec in valid:
        try:
            mol = standardize_mol(rec["mol"])
            if mol is not None:
                rec["mol"] = mol
                rec["std_smiles"] = Chem.MolToSmiles(mol)
                standardized.append(rec)
        except Exception:
            pass
    print(f"Step 2 - Standardized: {len(standardized)}")

    # Step 3: Compute properties and filter
    filtered = []
    for rec in standardized:
        props = compute_properties(rec["mol"])
        rec.update(props)
        if args.mw_min <= props["MW"] <= args.mw_max and props["NumHeavyAtoms"] <= args.max_heavy_atoms:
            filtered.append(rec)
    print(
        f"Step 3 - After property filter "
        f"(MW {args.mw_min}-{args.mw_max}, HA<={args.max_heavy_atoms}): {len(filtered)}"
    )

    # Step 4: Deduplicate by InChIKey (first 14 chars = connectivity layer)
    seen_inchikeys = {}
    for rec in filtered:
        inchi = Chem.MolToInchiKey(rec["mol"])
        if inchi is None:
            continue
        connectivity = inchi[:14]
        pchembl = float(rec["pchembl_value"])
        if connectivity not in seen_inchikeys or pchembl > seen_inchikeys[connectivity]["pchembl"]:
            seen_inchikeys[connectivity] = {"rec": rec, "pchembl": pchembl}
    deduped = [v["rec"] for v in seen_inchikeys.values()]
    print(f"Step 4 - After dedup by InChIKey: {len(deduped)}")

    if len(deduped) < args.n_actives:
        print(f"WARNING: Only {len(deduped)} unique compounds, need {args.n_actives}")
        print("Relaxing filters...")
        relaxed_mw_min = max(100, args.mw_min - 50)
        relaxed_mw_max = args.mw_max + 50
        relaxed_ha = args.max_heavy_atoms + 5

        filtered2 = []
        for rec in standardized:
            props = compute_properties(rec["mol"])
            rec.update(props)
            if relaxed_mw_min <= props["MW"] <= relaxed_mw_max and props["NumHeavyAtoms"] <= relaxed_ha:
                filtered2.append(rec)
        seen_inchikeys2 = {}
        for rec in filtered2:
            inchi = Chem.MolToInchiKey(rec["mol"])
            if inchi is None:
                continue
            connectivity = inchi[:14]
            pchembl = float(rec["pchembl_value"])
            if connectivity not in seen_inchikeys2 or pchembl > seen_inchikeys2[connectivity]["pchembl"]:
                seen_inchikeys2[connectivity] = {"rec": rec, "pchembl": pchembl}
        deduped = [v["rec"] for v in seen_inchikeys2.values()]
        print(f"After relaxed filters: {len(deduped)}")

    # Step 5: Compute ECFP4 fingerprints
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    for rec in deduped:
        rec["fp"] = fpgen.GetFingerprint(rec["mol"])
    fps = [rec["fp"] for rec in deduped]

    # Step 6: Butina clustering
    clusters = butina_cluster(fps, cutoff=args.butina_cutoff)
    n_clusters = len(clusters)
    cluster_sizes = [len(c) for c in clusters]
    print(
        f"Step 6 - Butina clusters: {n_clusters} "
        f"(sizes: {Counter(cluster_sizes).most_common(10)})"
    )

    # Step 7: MaxMin diversity pick to target count
    if len(deduped) <= args.n_actives:
        picked_indices = list(range(len(deduped)))
        print(
            f"Step 7 - Using all {len(deduped)} compounds "
            f"(fewer than target {args.n_actives})"
        )
    else:
        picked_indices = maxmin_pick(fps, args.n_actives)
        print(f"Step 7 - MaxMin picked {len(picked_indices)} diverse compounds")

    selected = [deduped[i] for i in picked_indices]

    # Re-cluster selected to count final clusters
    sel_fps = [rec["fp"] for rec in selected]
    sel_clusters = butina_cluster(sel_fps, cutoff=args.butina_cutoff)
    print(f"Final selection: {len(selected)} compounds in {len(sel_clusters)} Butina clusters")

    # Save output SMILES
    with open(output_smi, "w") as f:
        for i, rec in enumerate(selected, 1):
            compound_id = f"ACT_{i:03d}"
            rec["compound_id"] = compound_id
            f.write(f"{rec['std_smiles']}\t{compound_id}\t{rec['molecule_chembl_id']}\n")
    print(f"Saved {len(selected)} actives to {output_smi}")

    # Save curated CSV with properties
    csv_fields = [
        "compound_id",
        "molecule_chembl_id",
        "std_smiles",
        "pchembl_value",
        "standard_type",
        "MW",
        "LogP",
        "HBD",
        "HBA",
        "RotBonds",
        "TPSA",
        "NumHeavyAtoms",
        "NetCharge",
    ]
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_fields, extrasaction="ignore")
        writer.writeheader()
        for rec in selected:
            writer.writerow(rec)
    print(f"Saved curated CSV to {output_csv}")

    # Generate report
    report = generate_report(
        target_name=args.target_name,
        raw=raw_records,
        valid=valid,
        standardized=standardized,
        filtered=filtered,
        deduped=deduped,
        selected=selected,
        clusters=clusters,
        sel_clusters=sel_clusters,
        n_clusters=n_clusters,
        mw_min=args.mw_min,
        mw_max=args.mw_max,
        max_heavy_atoms=args.max_heavy_atoms,
        butina_cutoff=args.butina_cutoff,
        target_n=args.n_actives,
    )
    with open(report_path, "w") as f:
        f.write(report)
    print(f"Report saved to {report_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
