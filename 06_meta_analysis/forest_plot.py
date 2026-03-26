#!/usr/bin/env python3
"""DerSimonian-Laird random-effects meta-analysis with forest and funnel plots.

Reads delong_deltas.csv (produced by aggregate_metrics.py) and generates:
  - forest_plot.pdf / .svg / .png   (per-target + pooled effect)
  - funnel_plot.pdf                 (delta AUC vs SE)
  - meta_analysis_results.json      (numerical results)
"""

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Meta-analysis helpers
# ---------------------------------------------------------------------------

def dersimonian_laird(deltas: np.ndarray, ses: np.ndarray):
    """DerSimonian-Laird random-effects meta-analysis.

    Parameters
    ----------
    deltas : array of effect sizes (delta AUC)
    ses    : array of standard errors

    Returns
    -------
    dict with pooled estimate, CI, tau2, Q, I2, weights, etc.
    """
    k = len(deltas)
    if k == 0:
        return None

    w_fixed = 1.0 / (ses ** 2)
    sum_w = w_fixed.sum()
    pooled_fixed = (w_fixed * deltas).sum() / sum_w

    # Cochran's Q
    Q = (w_fixed * (deltas - pooled_fixed) ** 2).sum()
    df = k - 1

    # Q p-value (chi-squared with k-1 df)
    from scipy.stats import chi2
    Q_pvalue = 1.0 - chi2.cdf(Q, df) if df > 0 else 1.0

    # I-squared
    I2 = max(0.0, (Q - df) / Q * 100.0) if Q > 0 else 0.0

    # tau-squared (DerSimonian-Laird estimator)
    sum_w2 = (w_fixed ** 2).sum()
    C = sum_w - sum_w2 / sum_w
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0

    # Random-effects weights
    w_re = 1.0 / (ses ** 2 + tau2)
    sum_w_re = w_re.sum()
    pooled_re = (w_re * deltas).sum() / sum_w_re
    se_pooled = 1.0 / math.sqrt(sum_w_re)
    ci_low = pooled_re - 1.96 * se_pooled
    ci_high = pooled_re + 1.96 * se_pooled

    # z-test for pooled
    z = pooled_re / se_pooled if se_pooled > 0 else 0.0
    from scipy.stats import norm
    p_pooled = 2.0 * (1.0 - norm.cdf(abs(z)))

    return {
        "k": int(k),
        "pooled_delta": float(pooled_re),
        "se_pooled": float(se_pooled),
        "ci_low": float(ci_low),
        "ci_high": float(ci_high),
        "z": float(z),
        "p_value": float(p_pooled),
        "tau2": float(tau2),
        "tau": float(math.sqrt(tau2)),
        "Q": float(Q),
        "Q_df": int(df),
        "Q_pvalue": float(Q_pvalue),
        "I2": float(I2),
        "weights_re": w_re.tolist(),
        "pooled_fixed": float(pooled_fixed),
    }


# ---------------------------------------------------------------------------
# Plotting style helpers
# ---------------------------------------------------------------------------

# Target-class color palette (colorblind-safe)
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


def _setup_style():
    """Apply Nature-style formatting defaults."""
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


# ---------------------------------------------------------------------------
# Forest plot
# ---------------------------------------------------------------------------

def make_forest_plot(df: pd.DataFrame, meta: dict, output_dir: Path):
    """Generate forest plot: rows = targets, diamond = pooled."""
    _setup_style()

    # Sort by target class then name
    df_sorted = df.sort_values(["target_class", "target"]).reset_index(drop=True)
    n = len(df_sorted)

    fig_height = max(4.0, 1.0 + 0.45 * (n + 2))
    fig, ax = plt.subplots(figsize=(7.5, fig_height))

    y_positions = list(range(n + 1, 0, -1))  # targets at top, pooled at bottom
    y_target = y_positions[:n]
    y_pooled = y_positions[n] - 1.0  # extra gap before pooled

    # Per-target points and error bars
    for i, (_, row) in enumerate(df_sorted.iterrows()):
        y = y_target[i]
        color = CLASS_COLORS.get(row["target_class"], DEFAULT_COLOR)
        delta = row["delta_auc"]
        lo = row["ci_low"]
        hi = row["ci_high"]

        # Horizontal CI bar
        ax.plot([lo, hi], [y, y], color=color, linewidth=1.5, solid_capstyle="round")
        # Point estimate
        ax.plot(delta, y, "o", color=color, markersize=6, markeredgecolor="white",
                markeredgewidth=0.5, zorder=5)

        # Right-side annotation
        ann = f"{delta:+.3f} [{lo:+.3f}, {hi:+.3f}]"
        ax.annotate(ann, xy=(1.0, y), xycoords=("axes fraction", "data"),
                    fontsize=6.5, va="center", ha="left",
                    xytext=(6, 0), textcoords="offset points",
                    family="monospace")

    # Pooled diamond
    pooled = meta["pooled_delta"]
    plo = meta["ci_low"]
    phi = meta["ci_high"]
    diamond_hw = 0.35  # half-height of diamond
    diamond_x = [plo, pooled, phi, pooled]
    diamond_y = [y_pooled, y_pooled + diamond_hw, y_pooled,
                 y_pooled - diamond_hw]
    ax.fill(diamond_x, diamond_y, color="#333333", alpha=0.8, zorder=5)
    ax.plot(diamond_x + [diamond_x[0]], diamond_y + [diamond_y[0]],
            color="black", linewidth=0.8, zorder=6)

    # Pooled annotation
    ann_pooled = f"{pooled:+.3f} [{plo:+.3f}, {phi:+.3f}]"
    ax.annotate(ann_pooled, xy=(1.0, y_pooled),
                xycoords=("axes fraction", "data"),
                fontsize=6.5, va="center", ha="left",
                xytext=(6, 0), textcoords="offset points",
                family="monospace", fontweight="bold")

    # Vertical dashed line at 0
    ax.axvline(0, color="grey", linestyle="--", linewidth=0.8, zorder=1)

    # Horizontal separator before pooled
    sep_y = y_pooled + 0.6
    ax.axhline(sep_y, color="grey", linestyle="-", linewidth=0.5, alpha=0.5)

    # Y-axis labels
    y_labels = list(df_sorted["target"].str.upper()) + ["Pooled (RE)"]
    y_all = y_target + [y_pooled]
    ax.set_yticks(y_all)
    ax.set_yticklabels(y_labels, fontsize=8)

    # Bold the pooled label
    labels = ax.get_yticklabels()
    labels[-1].set_fontweight("bold")

    # Target class color patches on the left margin
    for i, (_, row) in enumerate(df_sorted.iterrows()):
        color = CLASS_COLORS.get(row["target_class"], DEFAULT_COLOR)
        ax.plot(-0.02, y_target[i], "s", color=color, markersize=5,
                transform=ax.get_yaxis_transform(), clip_on=False, zorder=10)

    # Axis labels and limits
    ax.set_xlabel("Delta AUC (skill - naive)", fontsize=9)
    ax.set_ylim(y_pooled - 1.0, y_target[0] + 1.0)

    # Auto x-limits with some padding
    all_lo = list(df_sorted["ci_low"]) + [plo]
    all_hi = list(df_sorted["ci_high"]) + [phi]
    x_pad = 0.02
    ax.set_xlim(min(all_lo) - x_pad, max(all_hi) + x_pad)

    # Title with heterogeneity statistics
    I2 = meta["I2"]
    Q = meta["Q"]
    Q_df = meta["Q_df"]
    Q_p = meta["Q_pvalue"]
    ax.set_title(
        f"Forest plot: delta ROC AUC (skill - naive)\n"
        f"Random-effects (DerSimonian-Laird)   "
        f"$I^2$ = {I2:.1f}%,  Q = {Q:.2f} (df={Q_df}, p={Q_p:.3f})",
        fontsize=10, pad=12)

    # Legend for target classes
    seen_classes = df_sorted["target_class"].unique()
    handles = []
    for tc in sorted(seen_classes):
        c = CLASS_COLORS.get(tc, DEFAULT_COLOR)
        handles.append(mpatches.Patch(color=c, label=tc.replace("_", " ")))
    ax.legend(handles=handles, loc="lower right", fontsize=6.5,
              framealpha=0.9, edgecolor="grey")

    # Spine cleanup
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(left=False)

    plt.tight_layout()

    for ext in ("pdf", "svg", "png"):
        outpath = output_dir / f"forest_plot.{ext}"
        fig.savefig(outpath, dpi=300)
        print(f"Wrote {outpath}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Funnel plot
# ---------------------------------------------------------------------------

def make_funnel_plot(df: pd.DataFrame, meta: dict, output_dir: Path):
    """Funnel plot: delta AUC (x) vs SE (y, inverted)."""
    _setup_style()

    fig, ax = plt.subplots(figsize=(5.5, 4.5))

    deltas = df["delta_auc"].values
    ses = df["se"].values
    classes = df["target_class"].values

    for d, s, tc in zip(deltas, ses, classes):
        c = CLASS_COLORS.get(tc, DEFAULT_COLOR)
        ax.plot(d, s, "o", color=c, markersize=6, markeredgecolor="white",
                markeredgewidth=0.5, zorder=5)

    # Pooled effect line
    pooled = meta["pooled_delta"]
    ax.axvline(pooled, color="#333333", linewidth=1.0, linestyle="-",
               label=f"Pooled = {pooled:.4f}")

    # Pseudo-95% CI funnel
    se_max = max(ses) * 1.15 if len(ses) > 0 else 0.1
    se_range = np.linspace(0, se_max, 200)
    ax.fill_betweenx(se_range,
                     pooled - 1.96 * se_range,
                     pooled + 1.96 * se_range,
                     color="grey", alpha=0.1, label="95% pseudo-CI")

    # No-effect line
    ax.axvline(0, color="grey", linewidth=0.7, linestyle="--", alpha=0.6)

    ax.set_xlabel("Delta AUC (skill - naive)", fontsize=9)
    ax.set_ylabel("Standard error", fontsize=9)
    ax.set_title("Funnel plot", fontsize=11)
    ax.invert_yaxis()
    ax.legend(fontsize=7, loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    outpath = output_dir / "funnel_plot.pdf"
    fig.savefig(outpath, dpi=300)
    print(f"Wrote {outpath}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="DerSimonian-Laird random-effects meta-analysis with "
                    "forest and funnel plots.")
    parser.add_argument("--input", default="delong_deltas.csv",
                        help="Path to delong_deltas.csv")
    parser.add_argument("--output-dir", default="meta_figures",
                        help="Directory for output figures and results")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    df = pd.read_csv(input_path)
    required_cols = {"target", "delta_auc", "ci_low", "ci_high", "se"}
    missing = required_cols - set(df.columns)
    if missing:
        sys.exit(f"ERROR: input CSV missing columns: {missing}")

    # Drop rows with NaN SE (cannot contribute to meta-analysis)
    valid = df.dropna(subset=["delta_auc", "se"])
    invalid = df[~df.index.isin(valid.index)]
    if len(invalid) > 0:
        print(f"WARNING: dropping {len(invalid)} targets with missing "
              f"delta/SE: {list(invalid['target'])}", file=sys.stderr)
    df = valid.reset_index(drop=True)

    if len(df) < 2:
        sys.exit("ERROR: need at least 2 targets with valid data for "
                 "meta-analysis")

    print(f"Loaded {len(df)} targets from {input_path}")

    # Ensure target_class column exists
    if "target_class" not in df.columns:
        df["target_class"] = "unknown"

    # ------------------------------------------------------------------
    # Run meta-analysis
    # ------------------------------------------------------------------
    deltas = df["delta_auc"].values
    ses = df["se"].values
    meta = dersimonian_laird(deltas, ses)

    print(f"\n--- Random-effects meta-analysis ---")
    print(f"  k = {meta['k']} studies")
    print(f"  Pooled delta AUC = {meta['pooled_delta']:.4f}  "
          f"95% CI [{meta['ci_low']:.4f}, {meta['ci_high']:.4f}]")
    print(f"  z = {meta['z']:.3f},  p = {meta['p_value']:.4g}")
    print(f"  tau^2 = {meta['tau2']:.6f},  tau = {meta['tau']:.4f}")
    print(f"  Q = {meta['Q']:.2f} (df={meta['Q_df']}, p={meta['Q_pvalue']:.4f})")
    print(f"  I^2 = {meta['I2']:.1f}%")

    # ------------------------------------------------------------------
    # Generate plots
    # ------------------------------------------------------------------
    make_forest_plot(df, meta, output_dir)
    make_funnel_plot(df, meta, output_dir)

    # ------------------------------------------------------------------
    # Save numerical results
    # ------------------------------------------------------------------
    results = {
        "method": "DerSimonian-Laird random-effects",
        "metric": "delta_ROC_AUC (skill - naive)",
        **meta,
        "per_target": df[["target", "target_class", "delta_auc",
                          "ci_low", "ci_high", "se"]].to_dict(orient="records"),
    }
    json_path = output_dir / "meta_analysis_results.json"
    with open(json_path, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nWrote {json_path}")


if __name__ == "__main__":
    main()
