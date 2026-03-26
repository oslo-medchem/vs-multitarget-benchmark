#!/usr/bin/env python3
"""
Skill-guided protocol: RDKit ETKDGv3 + Meeko ligand preparation.

Applies PAINS + Brenk filtering, generates 3D conformers with ETKDGv3,
optimises with MMFF, and converts to PDBQT via Meeko (obabel fallback).
Uses multiprocessing with (cpu_count - 2) workers.

Usage:
    python prepare_skill.py \
        --library-csv 03_library_preparation/library_labels.csv \
        --output-dir 03_library_preparation/skill/
"""

import argparse
import csv
import multiprocessing as mp
import os
import subprocess
import sys
import tempfile
from functools import partial
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

TIMEOUT = 90  # seconds per compound
MAX_WORKERS = max(1, mp.cpu_count() - 2)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Prepare ligand PDBQT files using skill-guided protocol "
        "(PAINS/Brenk filter, RDKit ETKDGv3, Meeko)."
    )
    p.add_argument(
        "--library-csv",
        required=True,
        type=Path,
        help="Path to library_labels.csv (must contain 'smiles' and 'compound_id' columns)",
    )
    p.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Output directory for pdbqt/, preparation_log.csv, filter_report.csv",
    )
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# PAINS + Brenk filter catalog (built per-worker to avoid pickling issues)
# ---------------------------------------------------------------------------

def build_filter() -> FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog(params)


# ---------------------------------------------------------------------------
# Per-compound preparation
# ---------------------------------------------------------------------------

def prepare_one(
    row: dict,
    output_dir: Path,
    apply_filter: bool = True,
) -> tuple[str, str, int, str]:
    """Prepare a single compound.

    Returns (compound_id, status, pdbqt_size, note).
    Status is one of: ok, ok_obabel_fallback, filtered, failed, skipped.
    """
    compound_id = row["compound_id"]
    smiles = row["smiles"]
    output_pdbqt = output_dir / f"{compound_id}.pdbqt"

    # Skip if already prepared and non-trivial
    if output_pdbqt.exists() and output_pdbqt.stat().st_size > 100:
        try:
            first4 = open(output_pdbqt, "rb").read(4)
            if first4 != b"\x00\x00\x00\x00":
                return compound_id, "skipped", 0, ""
        except Exception:
            pass

    # Build filter catalog (per-worker, avoids pickle issues)
    catalog = build_filter()

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return compound_id, "failed", 0, "invalid_smiles"

    # PAINS / Brenk filter
    if apply_filter:
        entry = catalog.GetFirstMatch(mol)
        if entry is not None:
            filter_name = entry.GetDescription()
            return compound_id, "filtered", 0, filter_name

    # Add explicit hydrogens
    mol_h = Chem.AddHs(mol, addCoords=False)

    # Generate 3D conformer (ETKDGv3)
    params_etkdg = AllChem.ETKDGv3()
    params_etkdg.randomSeed = 42
    params_etkdg.numThreads = 1
    result = AllChem.EmbedMolecule(mol_h, params_etkdg)
    if result == -1:
        # Fallback: retry with different seed
        params_etkdg2 = AllChem.ETKDGv3()
        params_etkdg2.randomSeed = 0
        result = AllChem.EmbedMolecule(mol_h, params_etkdg2)
        if result == -1:
            return compound_id, "failed", 0, "embed_failed"

    # Optimise geometry with MMFF
    try:
        AllChem.MMFFOptimizeMolecule(mol_h, maxIters=2000)
    except Exception:
        pass  # Non-fatal; unoptimised conformer is still usable

    # Write 3D SDF to temp file
    sdf_path = None
    try:
        with tempfile.NamedTemporaryFile(
            suffix=".sdf", delete=False, mode="w"
        ) as sdf_f:
            sdf_path = sdf_f.name

        writer = Chem.SDWriter(sdf_path)
        writer.write(mol_h)
        writer.close()

        # Meeko: SDF -> PDBQT
        result_meeko = subprocess.run(
            ["mk_prepare_ligand.py", "-i", sdf_path, "-o", str(output_pdbqt)],
            capture_output=True,
            text=True,
            timeout=TIMEOUT,
        )

        if (
            result_meeko.returncode == 0
            and output_pdbqt.exists()
            and output_pdbqt.stat().st_size > 50
        ):
            size = output_pdbqt.stat().st_size
            return compound_id, "ok", size, ""

        # Fallback: obabel
        result_ob = subprocess.run(
            [
                "obabel", sdf_path, "-O", str(output_pdbqt),
                "--partialcharge", "gasteiger",
            ],
            capture_output=True,
            text=True,
            timeout=TIMEOUT,
        )
        if output_pdbqt.exists() and output_pdbqt.stat().st_size > 50:
            return (
                compound_id,
                "ok_obabel_fallback",
                output_pdbqt.stat().st_size,
                "",
            )

        return (
            compound_id,
            "failed",
            0,
            f"meeko: {result_meeko.stderr[:100]}",
        )

    except subprocess.TimeoutExpired:
        return compound_id, "failed", 0, "timeout"
    except Exception as e:
        return compound_id, "failed", 0, str(e)[:100]
    finally:
        if sdf_path is not None:
            try:
                os.unlink(sdf_path)
            except Exception:
                pass


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    args = parse_args(argv)
    output_dir = args.output_dir.resolve()
    pdbqt_dir = output_dir
    pdbqt_dir.mkdir(parents=True, exist_ok=True)

    print(f"Workers: {MAX_WORKERS}")
    print(f"Output : {pdbqt_dir}")

    # Load library
    with open(args.library_csv) as f:
        rows = list(csv.DictReader(f))
    print(f"Total compounds: {len(rows)}")

    # Validate required columns
    if rows and not {"smiles", "compound_id"}.issubset(rows[0].keys()):
        print(
            f"ERROR: library CSV must contain 'smiles' and 'compound_id' columns. "
            f"Found: {list(rows[0].keys())}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Process with multiprocessing
    prep_func = partial(prepare_one, output_dir=pdbqt_dir, apply_filter=True)

    log_rows: list[dict] = []
    filter_rows: list[dict] = []
    ok = failed = filtered = skipped = 0

    with mp.Pool(MAX_WORKERS) as pool:
        for i, result in enumerate(pool.imap(prep_func, rows, chunksize=20)):
            compound_id, status, size, msg = result
            log_rows.append(
                {
                    "compound_id": compound_id,
                    "status": status,
                    "pdbqt_size": size,
                    "note": msg,
                }
            )
            if status in ("ok", "ok_obabel_fallback", "skipped"):
                ok += 1
                if status == "skipped":
                    skipped += 1
            elif status == "filtered":
                filtered += 1
                filter_rows.append({"compound_id": compound_id, "filter_hit": msg})
            else:
                failed += 1

            if (i + 1) % 200 == 0:
                print(
                    f"  {i + 1}/{len(rows)}  ok={ok} filtered={filtered} "
                    f"failed={failed} skipped={skipped}"
                )

    print(
        f"\nDone: {ok} ok ({skipped} skipped), "
        f"{filtered} filtered, {failed} failed"
    )

    # Write preparation log
    log_file = output_dir / "preparation_log.csv"
    with open(log_file, "w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["compound_id", "status", "pdbqt_size", "note"]
        )
        w.writeheader()
        w.writerows(log_rows)
    print(f"Log: {log_file}")

    # Write filter report
    filter_file = output_dir / "filter_report.csv"
    with open(filter_file, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["compound_id", "filter_hit"])
        w.writeheader()
        w.writerows(filter_rows)
    print(f"Filter report: {filter_file}")

    # Final count
    pdbqt_count = len(list(pdbqt_dir.glob("*.pdbqt")))
    print(f"PDBQT files in output dir: {pdbqt_count}")


if __name__ == "__main__":
    main()
