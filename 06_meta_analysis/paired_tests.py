#!/usr/bin/env python3
"""Paired statistical tests across targets: skill vs naive for each metric.

Reads aggregate_metrics.csv (produced by aggregate_metrics.py) and runs:
  - Wilcoxon signed-rank test (two-sided) per metric
  - Sign test per metric
  - Cliff's delta effect size per metric
  - Benjamini-Hochberg correction across all tests

Produces:
  - paired_test_results.json
  - paired_test_summary.csv
  - metrics_heatmap.pdf / .svg  (targets x metrics, colored by delta)
"""

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


# ---------------------------------------------------------------------------
# Statistical helpers
# ---------------------------------------------------------------------------

def sign_test(diffs: np.ndarray):
    """Two-sided sign test on paired differences.

    Returns n_pos, n_neg, n_zero, p_value.
    """
    n_pos = int(np.sum(diffs > 0))
    n_neg = int(np.sum(diffs < 0))
    n_zero = int(np.sum(diffs == 0))
    n_valid = n_pos + n_neg
    if n_valid == 0:
        return n_pos, n_neg, n_zero, 1.0
    # Two-sided binomial test under H0: p = 0.5
    try:
        # scipy >= 1.7 uses binomtest
        from scipy.stats import binomtest
        result = binomtest(n_pos, n_valid, 0.5, alternative="two-sided")
        p = result.pvalue
    except ImportError:
        p = binom_test(n_pos, n_valid, 0.5)
    return n_pos, n_neg, n_zero, float(p)


def cliffs_delta(x: np.ndarray, y: np.ndarray):
    """Cliff's delta for paired observations.

    For paired data, this compares each x_i with each y_j (all pairs),
    computing (concordant - discordant) / total.
    """
    n = len(x)
    if n == 0:
        return 0.0
    count_more = 0
    count_less = 0
    total = 0
    for xi in x:
        for yj in y:
            total += 1
            if xi > yj:
                count_more += 1
            elif xi < yj:
                count_less += 1
    if total == 0:
        return 0.0
    return (count_more - count_less) / total


def cliffs_delta_interpret(d: float) -> str:
    """Interpret Cliff's delta magnitude."""
    ad = abs(d)
    if ad < 0.147:
        return "negligible"
    elif ad < 0.33:
        return "small"
    elif ad < 0.474:
        return "medium"
    else:
        return "large"


def benjamini_hochberg(p_values: list[float], alpha: float = 0.05):
    """Benjamini-Hochberg correction.

    Returns list of (p_adjusted, significant) tuples.
    """
    n = len(p_values)
    if n == 0:
        return []

    # Sort by p-value, keeping track of original indices
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    adjusted = [0.0] * n
    prev_adj = 0.0

    # Walk from largest to smallest p
    for rank_rev, (orig_idx, p) in enumerate(reversed(indexed)):
        rank = n - rank_rev  # 1-based rank
        adj = min(1.0, p * n / rank)
        if rank_rev > 0:
            adj = min(adj, prev_adj)
        adjusted[orig_idx] = adj
        prev_adj = adj

    # Re-walk forward to ensure monotonicity
    sorted_indices = [idx for idx, _ in indexed]
    running_min = 1.0
    for i in range(n - 1, -1, -1):
        idx = sorted_indices[i]
        adjusted[idx] = min(adjusted[idx], running_min)
        running_min = adjusted[idx]

    return [(adj, adj < alpha) for adj in adjusted]


# ---------------------------------------------------------------------------
# Heatmap
# ---------------------------------------------------------------------------

METRICS = ["roc_auc", "bedroc_20", "ef_1pct", "ef_5pct", "log_auc", "aupr"]
METRIC_LABELS = {
    "roc_auc": "ROC AUC",
    "bedroc_20": "BEDROC",
    "ef_1pct": "EF 1%",
    "ef_5pct": "EF 5%",
    "log_auc": "LogAUC",
    "aupr": "AUPR",
}


def _setup_style():
    """Nature-style formatting."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 11,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


# Target class colors (same palette as forest_plot.py)
CLASS_COLORS = {
    "GPCR":              "#E64B35",
    "kinase":            "#4DBBD5",
    "enzyme":            "#00A087",
    "nuclear_receptor":  "#3C5488",
    "serine_protease":   "#F39B7F",
    "hydrolase":         "#8491B4",
    "aspartyl_protease": "#91D1C2",
    "chaperone":         "#DC9A6C",
}
DEFAULT_COLOR = "#7E6148"


def make_heatmap(df_wide: pd.DataFrame, sig_matrix: pd.DataFrame,
                 target_classes: pd.Series, output_dir: Path):
    """Heatmap: rows=targets, cols=metrics, colored by delta (skill - naive).

    sig_matrix has same shape; True where BH-significant.
    """
    _setup_style()

    n_targets = len(df_wide)
    n_metrics = len(df_wide.columns)

    fig_height = max(3.5, 0.45 * n_targets + 1.5)
    fig_width = max(4.5, 0.9 * n_metrics + 2.5)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    data = df_wide.values
    vmax = np.nanmax(np.abs(data)) if not np.all(np.isnan(data)) else 0.1
    if vmax < 0.01:
        vmax = 0.1  # avoid degenerate colorbar

    cmap = plt.get_cmap("RdBu_r")
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    im = ax.imshow(data, cmap=cmap, norm=norm, aspect="auto")

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label("Delta (skill - naive)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    # Cell annotations: value + significance star
    for i in range(n_targets):
        for j in range(n_metrics):
            val = data[i, j]
            if math.isnan(val):
                txt = "n/a"
            else:
                txt = f"{val:+.3f}"
                if sig_matrix.iloc[i, j]:
                    txt += "*"
            # Choose text color for readability
            text_color = "white" if abs(val) > 0.6 * vmax else "black"
            ax.text(j, i, txt, ha="center", va="center", fontsize=6.5,
                    color=text_color)

    # Axis labels
    col_labels = [METRIC_LABELS.get(c, c) for c in df_wide.columns]
    ax.set_xticks(range(n_metrics))
    ax.set_xticklabels(col_labels, fontsize=8, rotation=30, ha="right")

    row_labels = [t.upper() for t in df_wide.index]
    ax.set_yticks(range(n_targets))
    ax.set_yticklabels(row_labels, fontsize=8)

    # Target class color patches on left margin
    for i, tgt in enumerate(df_wide.index):
        tc = target_classes.get(tgt, "unknown")
        c = CLASS_COLORS.get(tc, DEFAULT_COLOR)
        ax.plot(-0.7, i, "s", color=c, markersize=7, clip_on=False,
                transform=ax.transData, zorder=10)

    ax.set_title("Metric deltas: skill-guided vs naive\n"
                 "(* = BH-adjusted p < 0.05)", fontsize=10, pad=10)

    plt.tight_layout()
    for ext in ("pdf", "svg"):
        outpath = output_dir / f"metrics_heatmap.{ext}"
        fig.savefig(outpath, dpi=300)
        print(f"Wrote {outpath}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Paired statistical tests (skill vs naive) across targets.")
    parser.add_argument("--input", default="aggregate_metrics.csv",
                        help="Path to aggregate_metrics.csv")
    parser.add_argument("--output-dir", default=".",
                        help="Directory for output files")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Load and pivot data
    # ------------------------------------------------------------------
    df = pd.read_csv(input_path)
    required = {"target", "protocol", "target_class"} | set(METRICS)
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"ERROR: missing columns in input: {missing}")

    # Build lookup: target -> target_class
    tclass_map = df.groupby("target")["target_class"].first()

    # Separate protocols
    naive = df[df["protocol"] == "naive"].set_index("target")
    skill = df[df["protocol"] == "skill"].set_index("target")

    # Intersect targets present in both
    common_targets = sorted(set(naive.index) & set(skill.index))
    if len(common_targets) < 2:
        sys.exit(f"ERROR: need >= 2 targets with both protocols; "
                 f"found {len(common_targets)}")

    naive = naive.loc[common_targets]
    skill = skill.loc[common_targets]

    print(f"Loaded {len(common_targets)} targets with both protocols.")

    # ------------------------------------------------------------------
    # Paired tests for each metric
    # ------------------------------------------------------------------
    results = {}
    all_p_values = []
    p_labels = []

    # Also build delta matrix for heatmap
    delta_data = {}

    for metric in METRICS:
        n_vals = naive[metric].values.astype(float)
        s_vals = skill[metric].values.astype(float)

        # Drop pairs where either is NaN
        valid = ~(np.isnan(n_vals) | np.isnan(s_vals))
        n_v = n_vals[valid]
        s_v = s_vals[valid]
        diffs = s_v - n_v
        n_valid = len(diffs)

        delta_data[metric] = {t: float(s_vals[i] - n_vals[i])
                              for i, t in enumerate(common_targets)
                              if valid[i]}

        if n_valid < 2:
            results[metric] = {
                "n_targets": n_valid,
                "note": "insufficient data for test",
            }
            all_p_values.append(1.0)
            p_labels.append(metric)
            continue

        # Wilcoxon signed-rank (two-sided)
        # Handle case where all diffs are zero
        if np.all(diffs == 0):
            w_stat, w_p = 0.0, 1.0
        else:
            try:
                w_stat, w_p = wilcoxon(diffs, alternative="two-sided")
            except ValueError:
                # All differences are zero
                w_stat, w_p = 0.0, 1.0

        # Sign test
        n_pos, n_neg, n_zero, sign_p = sign_test(diffs)

        # Cliff's delta (all pairs from skill vs naive)
        cd = cliffs_delta(s_v, n_v)
        cd_interp = cliffs_delta_interpret(cd)

        # Summary stats
        mean_diff = float(np.mean(diffs))
        median_diff = float(np.median(diffs))
        mean_naive = float(np.mean(n_v))
        mean_skill = float(np.mean(s_v))

        results[metric] = {
            "n_targets": n_valid,
            "mean_naive": mean_naive,
            "mean_skill": mean_skill,
            "mean_delta": mean_diff,
            "median_delta": median_diff,
            "wilcoxon_stat": float(w_stat),
            "wilcoxon_p": float(w_p),
            "sign_test_pos": n_pos,
            "sign_test_neg": n_neg,
            "sign_test_zero": n_zero,
            "sign_test_p": float(sign_p),
            "cliffs_delta": float(cd),
            "cliffs_delta_interp": cd_interp,
        }

        all_p_values.append(float(w_p))
        p_labels.append(metric)

    # ------------------------------------------------------------------
    # Benjamini-Hochberg correction
    # ------------------------------------------------------------------
    bh_results = benjamini_hochberg(all_p_values)
    for i, label in enumerate(p_labels):
        p_adj, sig = bh_results[i]
        if label in results and "wilcoxon_p" in results[label]:
            results[label]["wilcoxon_p_bh"] = float(p_adj)
            results[label]["significant_bh"] = bool(sig)

    # ------------------------------------------------------------------
    # Print summary table
    # ------------------------------------------------------------------
    print("\n" + "=" * 100)
    print(f"{'Metric':>12s}  {'n':>3s}  {'Mean naive':>11s}  "
          f"{'Mean skill':>11s}  {'Mean delta':>11s}  "
          f"{'Wilcoxon p':>11s}  {'p (BH)':>8s}  "
          f"{'Cliff d':>8s}  {'Interp':>10s}  {'Sign +/-':>8s}")
    print("-" * 100)

    summary_rows = []
    for metric in METRICS:
        r = results.get(metric, {})
        if "wilcoxon_p" not in r:
            print(f"{metric:>12s}  insufficient data")
            continue

        sig_marker = " *" if r.get("significant_bh", False) else "  "
        print(f"{metric:>12s}  {r['n_targets']:3d}  "
              f"{r['mean_naive']:11.4f}  {r['mean_skill']:11.4f}  "
              f"{r['mean_delta']:+11.4f}  "
              f"{r['wilcoxon_p']:11.4g}  {r['wilcoxon_p_bh']:8.4g}{sig_marker}  "
              f"{r['cliffs_delta']:+8.3f}  {r['cliffs_delta_interp']:>10s}  "
              f"{r['sign_test_pos']:d}/{r['sign_test_neg']:d}")

        summary_rows.append({
            "metric": metric,
            "n_targets": r["n_targets"],
            "mean_naive": r["mean_naive"],
            "mean_skill": r["mean_skill"],
            "mean_delta": r["mean_delta"],
            "median_delta": r["median_delta"],
            "wilcoxon_stat": r["wilcoxon_stat"],
            "wilcoxon_p": r["wilcoxon_p"],
            "wilcoxon_p_bh": r.get("wilcoxon_p_bh", float("nan")),
            "significant_bh": r.get("significant_bh", False),
            "sign_pos": r["sign_test_pos"],
            "sign_neg": r["sign_test_neg"],
            "sign_test_p": r["sign_test_p"],
            "cliffs_delta": r["cliffs_delta"],
            "cliffs_delta_interp": r["cliffs_delta_interp"],
        })

    print("=" * 100)
    print("* = significant after BH correction (alpha=0.05)")

    # ------------------------------------------------------------------
    # Save results
    # ------------------------------------------------------------------
    json_path = output_dir / "paired_test_results.json"
    with open(json_path, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nWrote {json_path}")

    if summary_rows:
        df_summary = pd.DataFrame(summary_rows)
        csv_path = output_dir / "paired_test_summary.csv"
        df_summary.to_csv(csv_path, index=False, float_format="%.6f")
        print(f"Wrote {csv_path}")

    # ------------------------------------------------------------------
    # Generate heatmap
    # ------------------------------------------------------------------
    # Build delta matrix: targets x metrics
    # Sort targets by class then name
    sorted_targets = sorted(common_targets,
                            key=lambda t: (tclass_map.get(t, "zzz"), t))

    delta_matrix = pd.DataFrame(index=sorted_targets, columns=METRICS,
                                dtype=float)
    for metric in METRICS:
        for tgt in sorted_targets:
            delta_matrix.loc[tgt, metric] = delta_data.get(metric, {}).get(
                tgt, float("nan"))

    # Build significance matrix (BH-significant per metric, applied to
    # all targets; we flag cells where that metric is BH-significant)
    sig_matrix = pd.DataFrame(False, index=sorted_targets, columns=METRICS)
    for metric in METRICS:
        is_sig = results.get(metric, {}).get("significant_bh", False)
        sig_matrix[metric] = is_sig

    make_heatmap(delta_matrix, sig_matrix,
                 tclass_map, output_dir)


if __name__ == "__main__":
    main()
