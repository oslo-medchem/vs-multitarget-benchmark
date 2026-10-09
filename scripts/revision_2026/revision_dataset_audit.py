#!/usr/bin/env python3
"""Dataset-construction audit: assay provenance and decoy cross-reactivity.

Referee 1 point 3 asks about cross-reactivity risk in the shared decoy pool for
the kinase and protease groups. Referee 2 point 5 asks whether binding and
functional assays were distinguished, whether ChEMBL target-assignment
confidence was restricted, how conflicting activity records were resolved,
whether the decoy similarity restriction was applied against all qualifying
actives or only the 100 selected, and whether other known actives were
explicitly excluded.

The answers are read out of the benchmark's own code and data rather than
asserted:

  * common/fetch_actives.py queries ChEMBL with target_chembl_id,
    pchembl_value >= 5, standard_relation "=" and standard_type in
    {EC50, IC50, Ki, Kd}. It does NOT restrict assay_type and does NOT restrict
    the target-assignment confidence score. assay_type is preserved in
    actives_raw.csv, so the binding/functional mix is quantified here.
  * common/curate_actives.py deduplicates on the InChIKey connectivity layer
    and keeps the record with the HIGHEST pChEMBL value, so conflicting records
    are resolved by maximum reported potency.
  * common/generate_decoys.py builds its InChIKey exclusion set and its
    Tanimoto reference set from actives_curated.csv, i.e. from the 100 SELECTED
    actives only. Qualifying actives that were not selected were therefore not
    excluded from the decoy pool, and neither were actives of other targets.

This script quantifies the consequences of the last point.

Outputs assay_provenance.csv, decoy_contamination.csv and
decoy_contamination_pairs.csv in data/revision_2026/.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import revision_common as rc

KINASES = {"egfr", "cdk2", "p38"}
PROTEASES = {"dpp4", "thrombin", "bace1"}


def read_raw_actives(root: Path, t: str):
    """All ChEMBL records that met the query, before diversity selection."""
    p = root / "targets" / t / "01_active_curation" / "actives_raw.csv"
    rows = []
    if p.exists():
        with p.open() as f:
            rows = list(csv.DictReader(f))
    return rows


def read_selected(root: Path, t: str):
    p = root / "targets" / t / "01_active_curation" / "actives_curated.csv"
    if not p.exists():
        return []
    with p.open() as f:
        return list(csv.DictReader(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--out-dir", default="data/revision_2026")
    args = ap.parse_args()
    root, out = Path(args.data_root), Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # qualifying ChEMBL molecule ids per target (pChEMBL >= 5 by construction)
    qualifying: dict[str, set] = {}
    best_p: dict[str, dict] = {}
    prov_rows = []
    for t in rc.TARGETS:
        raw = read_raw_actives(root, t)
        sel = read_selected(root, t)
        ids = {r["molecule_chembl_id"] for r in raw}
        qualifying[t] = ids
        bp = {}
        for r in raw:
            try:
                v = float(r["pchembl_value"])
            except (ValueError, KeyError, TypeError):
                continue
            mid = r["molecule_chembl_id"]
            if mid not in bp or v > bp[mid]:
                bp[mid] = v
        best_p[t] = bp
        at = Counter(r.get("assay_type", "") or "unspecified" for r in raw)
        sel_ids = {r["molecule_chembl_id"] for r in sel}
        at_sel = Counter()
        for r in raw:
            if r["molecule_chembl_id"] in sel_ids:
                at_sel[r.get("assay_type", "") or "unspecified"] += 1
        st = Counter(r.get("standard_type", "") for r in raw)
        prov_rows.append({
            "target": t,
            "n_raw_records": len(raw),
            "n_unique_qualifying_molecules": len(ids),
            "n_selected_actives": len(sel),
            "assay_type_B_binding": at.get("B", 0),
            "assay_type_F_functional": at.get("F", 0),
            "assay_type_other": sum(v for k, v in at.items() if k not in ("B", "F")),
            "pct_binding_records": round(100 * at.get("B", 0) / len(raw), 1) if raw else "",
            "selected_from_B": at_sel.get("B", 0),
            "selected_from_F": at_sel.get("F", 0),
            "std_type_IC50": st.get("IC50", 0), "std_type_Ki": st.get("Ki", 0),
            "std_type_Kd": st.get("Kd", 0), "std_type_EC50": st.get("EC50", 0),
            "n_molecules_with_conflicting_records":
                sum(1 for mid in ids if sum(1 for r in raw if r["molecule_chembl_id"] == mid) > 1),
        })

    # ---- decoy contamination ---------------------------------------------
    cont_rows, pair_rows = [], []
    for t in rc.TARGETS:
        td = rc.load_target(t, root)
        decoy_src = {cid: td.source_id.get(cid, "") for cid, lab in td.labels.items() if lab == 0}
        n_dec = len(decoy_src)

        same_target = {cid for cid, sid in decoy_src.items() if sid in qualifying.get(t, set())}

        other_hits = defaultdict(set)
        for cid, sid in decoy_src.items():
            if not sid:
                continue
            for other in rc.TARGETS:
                if other == t:
                    continue
                if sid in qualifying.get(other, set()):
                    other_hits[other].add(cid)

        any_other = set().union(*other_hits.values()) if other_hits else set()
        fam = KINASES if t in KINASES else (PROTEASES if t in PROTEASES else set())
        within_family = set().union(*[other_hits[o] for o in other_hits if o in fam]) \
            if any(o in fam for o in other_hits) else set()

        cont_rows.append({
            "target": t, "n_decoys": n_dec,
            "decoys_that_are_qualifying_actives_of_same_target": len(same_target),
            "pct_same_target": round(100 * len(same_target) / n_dec, 2) if n_dec else "",
            "decoys_that_are_qualifying_actives_of_any_other_panel_target": len(any_other),
            "pct_any_other": round(100 * len(any_other) / n_dec, 2) if n_dec else "",
            "decoys_active_at_a_same_family_target": len(within_family),
            "pct_same_family": round(100 * len(within_family) / n_dec, 2) if n_dec else "",
            "family": "kinase" if t in KINASES else ("protease" if t in PROTEASES else "other"),
        })
        for other, cids in sorted(other_hits.items()):
            pair_rows.append({"decoy_library_of": t, "also_active_at": other,
                              "n_decoys": len(cids),
                              "same_family": other in fam})

    for name, rows in (("assay_provenance.csv", prov_rows),
                       ("decoy_contamination.csv", cont_rows),
                       ("decoy_contamination_pairs.csv", pair_rows)):
        with (out / name).open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {out/name} ({len(rows)} rows)")

    print("\n=== assay provenance ===")
    for r in prov_rows:
        print(f"  {r['target']:9s} raw={r['n_raw_records']:6d} uniq={r['n_unique_qualifying_molecules']:5d} "
              f"B={r['assay_type_B_binding']:5d} F={r['assay_type_F_functional']:5d} "
              f"({r['pct_binding_records']}% binding)  selected B/F={r['selected_from_B']}/{r['selected_from_F']}")
    print("\n=== decoy contamination ===")
    for r in cont_rows:
        print(f"  {r['target']:9s} decoys={r['n_decoys']:5d} "
              f"same-target={r['decoys_that_are_qualifying_actives_of_same_target']:4d} ({r['pct_same_target']}%) "
              f"other-panel={r['decoys_that_are_qualifying_actives_of_any_other_panel_target']:4d} ({r['pct_any_other']}%) "
              f"same-family={r['decoys_active_at_a_same_family_target']:3d} ({r['pct_same_family']}%)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
