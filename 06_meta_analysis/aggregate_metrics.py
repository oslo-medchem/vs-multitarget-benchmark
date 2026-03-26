#!/usr/bin/env python3
"""Aggregate per-target evaluation metrics into cross-target tables.

Reads metrics_naive.json, metrics_skill.json, and comparison_stats.json
for every target listed in target_config.json.  Produces:
  - aggregate_metrics.csv / .json  (long-form: one row per target x protocol)
  - delong_deltas.csv              (one row per target with DeLong delta AUC)
"""

import argparse
import json
import math
import sys
from pathlib import Path

import pandas as pd


def load_json(path: Path) -> dict | None:
    """Return parsed JSON or None if file missing / broken."""
    try:
        with open(path) as fh:
            return json.load(fh)
    except (FileNotFoundError, json.JSONDecodeError) as exc:
        print(f"  WARNING: cannot load {path}: {exc}", file=sys.stderr)
        return None


def extract_protocol_row(target_name: str, target_class: str,
                         metrics: dict) -> dict:
    """Build one row of the aggregate table from a metrics JSON."""
    ci = metrics.get("bootstrap_ci", {})

    def ci_val(metric, bound):
        entry = ci.get(metric, {})
        v = entry.get(bound, float("nan"))
        # Handle JSON NaN (stored as null or literal NaN)
        if v is None or (isinstance(v, float) and math.isnan(v)):
            return float("nan")
        return v

    return {
        "target": target_name,
        "target_class": target_class,
        "protocol": metrics.get("protocol", "unknown"),
        "n_compounds": metrics.get("n_compounds", 0),
        "n_docked": metrics.get("n_docked", 0),
        "n_actives": metrics.get("n_actives", 0),
        "n_decoys": metrics.get("n_decoys", 0),
        "roc_auc": metrics.get("roc_auc", float("nan")),
        "roc_auc_ci_lo": ci_val("roc_auc", "ci_low"),
        "roc_auc_ci_hi": ci_val("roc_auc", "ci_high"),
        "bedroc_20": metrics.get("bedroc_20", float("nan")),
        "ef_1pct": metrics.get("ef_1pct", float("nan")),
        "ef_5pct": metrics.get("ef_5pct", float("nan")),
        "ef_10pct": metrics.get("ef_10pct", float("nan")),
        "log_auc": metrics.get("log_auc", float("nan")),
        "aupr": metrics.get("aupr", float("nan")),
    }


def extract_delong_row(target_name: str, target_class: str,
                       comparison: dict) -> dict | None:
    """Extract DeLong delta AUC from comparison_stats.json."""
    delong = comparison.get("delong_auc")
    if delong is None:
        return None

    delta = delong.get("delta_auc", float("nan"))
    ci_lo = delong.get("ci_95_low", float("nan"))
    ci_hi = delong.get("ci_95_high", float("nan"))
    p_val = delong.get("p_value", float("nan"))
    z_stat = delong.get("z_stat", float("nan"))

    # Compute SE from 95 % CI: SE = (ci_hi - ci_lo) / (2 * 1.96)
    if not (math.isnan(ci_lo) or math.isnan(ci_hi)):
        se = (ci_hi - ci_lo) / (2 * 1.96)
    else:
        se = float("nan")

    return {
        "target": target_name,
        "target_class": target_class,
        "delta_auc": delta,
        "ci_low": ci_lo,
        "ci_high": ci_hi,
        "se": se,
        "z_stat": z_stat,
        "p_value": p_val,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate per-target VS metrics into cross-target tables.")
    parser.add_argument("--config", default="../target_config.json",
                        help="Path to target_config.json")
    parser.add_argument("--targets-dir", default="../targets",
                        help="Base directory containing per-target results")
    parser.add_argument("--output-dir", default=".",
                        help="Directory for output files")
    args = parser.parse_args()

    config_path = Path(args.config)
    targets_dir = Path(args.targets_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Load target config
    # ------------------------------------------------------------------
    with open(config_path) as fh:
        config = json.load(fh)

    targets = config.get("targets", [])
    if not targets:
        sys.exit("ERROR: no targets found in config")

    # ------------------------------------------------------------------
    # Iterate targets and collect rows
    # ------------------------------------------------------------------
    metric_rows = []
    delong_rows = []
    skipped = []

    for tgt in targets:
        name = tgt["name"]
        tclass = tgt.get("target_class", "unknown")
        results_dir = targets_dir / name / "05_evaluation" / "results"

        print(f"Processing {name} ({tclass}) ...")

        m_naive = load_json(results_dir / "metrics_naive.json")
        m_skill = load_json(results_dir / "metrics_skill.json")
        comp = load_json(results_dir / "comparison_stats.json")

        if m_naive is None or m_skill is None:
            print(f"  SKIPPED (missing metrics files)")
            skipped.append(name)
            continue

        metric_rows.append(extract_protocol_row(name, tclass, m_naive))
        metric_rows.append(extract_protocol_row(name, tclass, m_skill))

        if comp is not None:
            drow = extract_delong_row(name, tclass, comp)
            if drow is not None:
                delong_rows.append(drow)

    if not metric_rows:
        sys.exit("ERROR: no targets had valid metrics files")

    # ------------------------------------------------------------------
    # Build DataFrames and save
    # ------------------------------------------------------------------
    df_metrics = pd.DataFrame(metric_rows)
    df_delong = pd.DataFrame(delong_rows) if delong_rows else pd.DataFrame()

    # CSV
    csv_metrics = output_dir / "aggregate_metrics.csv"
    df_metrics.to_csv(csv_metrics, index=False, float_format="%.6f")
    print(f"\nWrote {csv_metrics}  ({len(df_metrics)} rows, "
          f"{len(df_metrics) // 2} targets)")

    # JSON (records orient, human-readable)
    json_metrics = output_dir / "aggregate_metrics.json"
    with open(json_metrics, "w") as fh:
        json.dump(df_metrics.to_dict(orient="records"), fh, indent=2,
                  default=str)
    print(f"Wrote {json_metrics}")

    # DeLong deltas
    if not df_delong.empty:
        csv_delong = output_dir / "delong_deltas.csv"
        df_delong.to_csv(csv_delong, index=False, float_format="%.6f")
        print(f"Wrote {csv_delong}  ({len(df_delong)} targets)")
    else:
        print("WARNING: no DeLong delta data found; delong_deltas.csv not written")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    if skipped:
        print(f"\nSkipped targets (missing data): {', '.join(skipped)}")

    print("\n--- Aggregate summary ---")
    for proto in ("naive", "skill"):
        sub = df_metrics[df_metrics["protocol"] == proto]
        if sub.empty:
            continue
        print(f"\n{proto.upper()} (n={len(sub)} targets):")
        for col in ("roc_auc", "bedroc_20", "ef_1pct", "ef_5pct",
                     "log_auc", "aupr"):
            vals = sub[col].dropna()
            if len(vals) > 0:
                print(f"  {col:>12s}: mean={vals.mean():.4f}  "
                      f"median={vals.median():.4f}  "
                      f"range=[{vals.min():.4f}, {vals.max():.4f}]")

    if not df_delong.empty:
        print(f"\nDeLong delta AUC (skill - naive):")
        d = df_delong["delta_auc"]
        print(f"  mean={d.mean():.4f}  median={d.median():.4f}  "
              f"range=[{d.min():.4f}, {d.max():.4f}]")
        sig = df_delong[df_delong["p_value"] < 0.05]
        print(f"  Significant (p<0.05): {len(sig)}/{len(df_delong)}")


if __name__ == "__main__":
    main()
