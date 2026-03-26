#!/usr/bin/env python3
"""
Generate all figures for a single-target VS benchmark evaluation.

Usage:
    python generate_figures.py --target-dir targets/egfr/
    python generate_figures.py --target-dir targets/egfr/ --target-name "EGFR"

Figures:
    fig_main.pdf/svg  -- 3-panel composite (ROC, violin, box)
    figS1.pdf/svg     -- Semi-log ROC (LogAUC)
    figS2.pdf/svg     -- Precision-Recall curves
    figS3.pdf/svg     -- Property QC distributions
    figS4.pdf/svg     -- Score rank plots
    figS5.pdf/svg     -- Full metrics heatmap

All saved to target_dir/05_evaluation/figures/.
"""
import argparse
import json
import csv
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrowPatch
from pathlib import Path
from sklearn.metrics import roc_curve, precision_recall_curve, auc as sklearn_auc

# ---------------------------------------------------------------------------
# Nature color palette
# ---------------------------------------------------------------------------
PALETTE = {
    "naive": "#4477AA",   # blue
    "skill": "#EE6677",   # red/orange
    "naive_ci": "#4477AA44",
    "skill_ci": "#EE667744",
    "grey":   "#AAAAAA",
    "random": "#228833",
}

PROTOCOL_LABELS = {
    "naive": "Naive (OpenBabel)",
    "skill": "Skill-guided (Meeko/RDKit)",
}


# ---------------------------------------------------------------------------
# matplotlib rcParams -- Nature style
# ---------------------------------------------------------------------------
def set_nature_style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "lines.linewidth": 0.8,
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.minor.width": 0.3,
        "ytick.minor.width": 0.3,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "axes.grid": False,
    })


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_data(target_dir):
    """Load labels and scores for both protocols."""
    labels_csv = target_dir / "03_library_preparation" / "library_labels.csv"
    scores_files = {
        "naive": target_dir / "04_docking" / "naive" / "scores_naive.csv",
        "skill": target_dir / "04_docking" / "skill" / "scores_skill.csv",
    }

    labels = {}
    with open(labels_csv) as f:
        for row in csv.DictReader(f):
            labels[row["compound_id"]] = 1 if row["label"] == "active" else 0

    data = {}
    for protocol, score_file in scores_files.items():
        scores = {}
        if score_file.exists():
            with open(score_file) as f:
                for row in csv.DictReader(f):
                    try:
                        scores[row["compound_id"]] = float(row["score"])
                    except (ValueError, KeyError):
                        scores[row["compound_id"]] = 0.0
        cids = sorted(labels.keys())
        y_true = np.array([labels[c] for c in cids])
        y_score_raw = np.array([scores.get(c, 0.0) for c in cids])
        # Negate for ranking (lower docking score = better)
        y_score = -y_score_raw
        data[protocol] = {
            "cids": cids, "y_true": y_true,
            "y_score_raw": y_score_raw, "y_score": y_score,
        }
    return data


def load_metrics(results_dir, protocol):
    path = results_dir / f"metrics_{protocol}.json"
    if path.exists():
        return json.load(open(path))
    return {}


def load_roc(results_dir, protocol):
    path = results_dir / f"roc_{protocol}.json"
    if path.exists():
        return json.load(open(path))
    return {}


# ---------------------------------------------------------------------------
# Bootstrap ROC bands
# ---------------------------------------------------------------------------
def bootstrap_roc_band(y_true, y_score, n_boot=500, fpr_grid=None, seed=42):
    rng = np.random.default_rng(seed)
    n = len(y_true)
    if fpr_grid is None:
        fpr_grid = np.linspace(0, 1, 200)
    tpr_boot = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yt, ys = y_true[idx], y_score[idx]
        if yt.sum() == 0 or yt.sum() == n:
            continue
        fpr, tpr, _ = roc_curve(yt, ys)
        tpr_interp = np.interp(fpr_grid, fpr, tpr)
        tpr_boot.append(tpr_interp)
    tpr_boot = np.array(tpr_boot)
    return fpr_grid, np.percentile(tpr_boot, 2.5, axis=0), np.percentile(tpr_boot, 97.5, axis=0)


# ---------------------------------------------------------------------------
# Panel helper
# ---------------------------------------------------------------------------
def label_panel(ax, letter, x=-0.12, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=10, fontweight="bold", va="top", ha="left")


def save_figure(fig, stem, fig_dir):
    """Save as PDF, SVG, EPS, and high-resolution PNG (600 DPI)."""
    fig.savefig(fig_dir / f"{stem}.pdf")
    fig.savefig(fig_dir / f"{stem}.svg")
    fig.savefig(fig_dir / f"{stem}.eps")
    fig.savefig(fig_dir / f"{stem}.png", dpi=600)
    print(f"  Saved: {stem}.pdf + .svg + .eps + .png (600 DPI)")


# ---------------------------------------------------------------------------
# Main figure: 3-panel composite
# ---------------------------------------------------------------------------
def make_fig_main(data, metrics, results_dir, fig_dir, target_name):
    from scipy.stats import mannwhitneyu

    fig = plt.figure(figsize=(7.09, 3.0))  # 180 mm x 76 mm (Nature full width)
    gs = gridspec.GridSpec(1, 3, figure=fig, wspace=0.45,
                           left=0.065, right=0.98, bottom=0.15, top=0.90,
                           width_ratios=[1.15, 1.0, 1.0])

    fig.suptitle(f"{target_name} -- Virtual Screening Benchmark",
                 fontsize=9, fontweight="bold", y=0.98)

    # -- Panel A: ROC curves --
    ax_a = fig.add_subplot(gs[0])
    ax_a.plot([0, 1], [0, 1], "--", color=PALETTE["grey"], lw=0.6, zorder=1)

    for protocol, d in data.items():
        fpr, tpr, _ = roc_curve(d["y_true"], d["y_score"])
        roc_auc = sklearn_auc(fpr, tpr)
        m = metrics.get(protocol, {})
        ci = m.get("bootstrap_ci", {}).get("roc_auc", {})
        ci_lo = ci.get("ci_low", roc_auc)
        ci_hi = ci.get("ci_high", roc_auc)
        short = "Naive" if protocol == "naive" else "Skill"
        label = f"{short} (AUC = {roc_auc:.3f} [{ci_lo:.2f}\u2013{ci_hi:.2f}])"
        ax_a.plot(fpr, tpr, color=PALETTE[protocol], lw=1.2, label=label, zorder=3)

        # Bootstrap CI band
        fpr_g, tpr_lo, tpr_hi = bootstrap_roc_band(d["y_true"], d["y_score"], n_boot=300)
        ax_a.fill_between(fpr_g, tpr_lo, tpr_hi, color=PALETTE[protocol],
                          alpha=0.12, linewidth=0, zorder=2)

    ax_a.set_xlabel("False Positive Rate")
    ax_a.set_ylabel("True Positive Rate")
    ax_a.legend(loc="lower right", fontsize=6, handlelength=1.2,
                borderpad=0.4, labelspacing=0.3, handletextpad=0.5)
    ax_a.set_xlim(-0.02, 1.02)
    ax_a.set_ylim(-0.02, 1.02)
    label_panel(ax_a, "A")

    # -- Panel B: Bootstrap ROC AUC distributions --
    ax_b = fig.add_subplot(gs[1])

    rng = np.random.default_rng(42)
    n_boot = 2000
    auc_boots = {}
    for protocol, d in data.items():
        n = len(d["y_true"])
        aucs = []
        for _ in range(n_boot):
            idx = rng.integers(0, n, n)
            yt, ys = d["y_true"][idx], d["y_score"][idx]
            if yt.sum() == 0 or yt.sum() == n:
                continue
            fpr_b, tpr_b, _ = roc_curve(yt, ys)
            aucs.append(sklearn_auc(fpr_b, tpr_b))
        auc_boots[protocol] = np.array(aucs)

    # Paired violin plots
    positions = [0, 1]
    parts = ax_b.violinplot([auc_boots["naive"], auc_boots["skill"]],
                             positions=positions, showmedians=False,
                             showextrema=False, widths=0.65)
    for pc, protocol in zip(parts["bodies"], ["naive", "skill"]):
        pc.set_facecolor(PALETTE[protocol])
        pc.set_alpha(0.5)
        pc.set_edgecolor(PALETTE[protocol])
        pc.set_linewidth(0.6)

    # Overlay box plots for quartiles
    bp = ax_b.boxplot([auc_boots["naive"], auc_boots["skill"]],
                       positions=positions, widths=0.18,
                       patch_artist=True, showfliers=False,
                       medianprops=dict(color="white", linewidth=1.2),
                       whiskerprops=dict(color="#333333", linewidth=0.6),
                       capprops=dict(color="#333333", linewidth=0.6),
                       boxprops=dict(linewidth=0.5))
    for patch, protocol in zip(bp["boxes"], ["naive", "skill"]):
        patch.set_facecolor(PALETTE[protocol])
        patch.set_alpha(0.85)

    # Point estimates
    for i, protocol in enumerate(["naive", "skill"]):
        m = metrics.get(protocol, {})
        observed = m.get("roc_auc", np.median(auc_boots[protocol]))
        ax_b.plot(i, observed, "D", color="white", markersize=4,
                  markeredgecolor="#333333", markeredgewidth=0.6, zorder=5)

    # DeLong p-value bracket
    stats_file = results_dir / "comparison_stats.json"
    if stats_file.exists():
        comp_stats = json.load(open(stats_file))
        delong_p = comp_stats.get("delong_auc", {}).get("p_value", 0)
        delta = comp_stats.get("delong_auc", {}).get("delta_auc", 0)
    else:
        delong_p = 0
        delta = 0

    y_bracket = max(auc_boots["naive"].max(), auc_boots["skill"].max()) + 0.008
    ax_b.plot([0, 0, 1, 1], [y_bracket - 0.004, y_bracket, y_bracket, y_bracket - 0.004],
              color="#333333", lw=0.7)
    p_text = f"$\\Delta$AUC = +{delta:.3f}\nDeLong $p$ = {delong_p:.3f}"
    ax_b.text(0.5, y_bracket + 0.003, p_text, ha="center", va="bottom",
              fontsize=6, color="#333333")

    ax_b.set_xticks(positions)
    ax_b.set_xticklabels(["Naive", "Skill"], fontsize=7)
    ax_b.set_ylabel("ROC AUC (bootstrap)")
    ax_b.axhline(0.5, color=PALETTE["grey"], linestyle=":", lw=0.5)
    ax_b.text(1.05, 0.502, "random", color=PALETTE["grey"], fontsize=5,
              va="bottom", ha="left")
    label_panel(ax_b, "B")

    # -- Panel C: Box plots of docking scores --
    ax_c = fig.add_subplot(gs[2])

    box_data = []
    box_colors = []
    box_labels = []
    mw_results = []

    for protocol in ["naive", "skill"]:
        d = data[protocol]
        mask_a = d["y_true"] == 1
        mask_d = d["y_true"] == 0
        scores_a = d["y_score_raw"][mask_a]
        scores_d = d["y_score_raw"][mask_d]
        # Remove unscored (score==0)
        scores_a = scores_a[scores_a != 0]
        scores_d = scores_d[scores_d != 0]
        short = "Naive" if protocol == "naive" else "Skill"

        box_data.append(scores_a)
        box_colors.append(PALETTE[protocol])
        box_labels.append(f"{short}\nActives")

        box_data.append(scores_d)
        box_colors.append(PALETTE[protocol] + "66")
        box_labels.append(f"{short}\nDecoys")

        if len(scores_a) > 0 and len(scores_d) > 0:
            _, pval = mannwhitneyu(scores_a, scores_d, alternative="less")
            mw_results.append((len(box_data) - 2, len(box_data) - 1, pval))

    positions = np.array([0, 0.7, 1.8, 2.5])
    bp = ax_c.boxplot(box_data, positions=positions, widths=0.5,
                      patch_artist=True, showfliers=False,
                      medianprops=dict(color="black", linewidth=1.0),
                      whiskerprops=dict(linewidth=0.6),
                      capprops=dict(linewidth=0.6),
                      boxprops=dict(linewidth=0.5))

    for patch, color in zip(bp["boxes"], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)

    # Significance brackets
    for (i, j, pval) in mw_results:
        p_text = "***" if pval < 0.001 else ("**" if pval < 0.01 else ("*" if pval < 0.05 else "ns"))
        x1, x2 = positions[i], positions[j]
        all_scores = np.concatenate([box_data[i], box_data[j]])
        q1, q3 = np.percentile(all_scores, [25, 75])
        iqr = q3 - q1
        whisker_top = min(all_scores[all_scores <= q3 + 1.5 * iqr].max(),
                          all_scores.max()) if len(all_scores) > 0 else 0
        y = whisker_top - 0.3
        ax_c.plot([x1, x1, x2, x2], [y - 0.15, y, y, y - 0.15],
                  color="#333333", lw=0.6)
        ax_c.text((x1 + x2) / 2, y + 0.05, p_text, ha="center", va="bottom",
                  fontsize=6, color="#333333")

    ax_c.set_xticks(positions)
    ax_c.set_xticklabels(box_labels, fontsize=6)
    ax_c.set_ylabel("Docking Score (kcal/mol)")

    # Sensible y-axis from IQR (exclude extreme outliers)
    all_nonzero = np.concatenate([s for s in box_data if len(s) > 0])
    q1_all, q3_all = np.percentile(all_nonzero, [5, 95])
    margin = (q3_all - q1_all) * 0.3
    ax_c.set_ylim(q1_all - margin, q3_all + margin)
    label_panel(ax_c, "C")

    save_figure(fig, "fig_main", fig_dir)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure S1: Semi-log ROC
# ---------------------------------------------------------------------------
def make_fig_s1(data, metrics, fig_dir, target_name):
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 3.0))

    fpr_log = np.logspace(-3, 0, 1000)
    ax.plot(fpr_log, fpr_log, "--", color=PALETTE["grey"], lw=0.7, label="Random")

    for protocol, d in data.items():
        fpr, tpr, _ = roc_curve(d["y_true"], d["y_score"])
        tpr_interp = np.interp(fpr_log, fpr, tpr)
        m = metrics.get(protocol, {})
        logauc = m.get("log_auc", 0)
        label = f"{PROTOCOL_LABELS[protocol]}\nLogAUC = {logauc:.3f}"
        ax.semilogx(fpr_log, tpr_interp, color=PALETTE[protocol], lw=1.2, label=label)

    ax.set_xlabel("False Positive Rate (log scale)")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(f"Semi-log ROC -- {target_name}")
    ax.legend(fontsize=7)
    ax.set_xlim(1e-3, 1)
    save_figure(fig, "figS1_semilog_roc", fig_dir)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure S2: Precision-Recall
# ---------------------------------------------------------------------------
def make_fig_s2(data, metrics, fig_dir, target_name):
    fig, ax = plt.subplots(1, 1, figsize=(3.5, 3.0))

    for protocol, d in data.items():
        prec, rec, _ = precision_recall_curve(d["y_true"], d["y_score"])
        aupr = sklearn_auc(rec, prec)
        m = metrics.get(protocol, {})
        label = f"{PROTOCOL_LABELS[protocol]}\nAUPR = {aupr:.3f}"
        ax.step(rec, prec, color=PALETTE[protocol], lw=1.2, where="post", label=label)

    baseline = d["y_true"].mean()
    ax.axhline(baseline, color=PALETTE["grey"], linestyle="--", lw=0.7,
               label=f"Random (baseline = {baseline:.3f})")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title(f"Precision-Recall -- {target_name}")
    ax.legend(fontsize=7)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0, 1.05)
    save_figure(fig, "figS2_pr_curves", fig_dir)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure S3: Property QC distributions
# ---------------------------------------------------------------------------
def make_fig_s3(target_dir, fig_dir, target_name):
    """Re-plot property distributions from curated data files."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except ImportError:
        print("  Skipping S3 (rdkit not available)")
        return

    # Try to load curated CSVs with SMILES for property calculation
    actives_csv = target_dir / "01_active_curation" / "actives_curated.csv"
    decoys_csv = target_dir / "02_decoy_generation" / "decoys_curated.csv"
    labels_csv = target_dir / "03_library_preparation" / "library_labels.csv"

    actives_rows = []
    decoys_rows = []

    # Prefer curated CSVs if available; fallback to library_labels.csv
    if actives_csv.exists() and decoys_csv.exists():
        actives_rows = list(csv.DictReader(open(actives_csv)))
        decoys_rows = list(csv.DictReader(open(decoys_csv)))
        # Find the SMILES column name (could be std_smiles or smiles)
        act_smi_key = "std_smiles" if "std_smiles" in actives_rows[0] else "smiles"
        dec_smi_key = "std_smiles" if "std_smiles" in decoys_rows[0] else "smiles"
    elif labels_csv.exists():
        rows = list(csv.DictReader(open(labels_csv)))
        actives_rows = [r for r in rows if r["label"] == "active"]
        decoys_rows = [r for r in rows if r["label"] == "decoy"]
        act_smi_key = "smiles"
        dec_smi_key = "smiles"
    else:
        print("  Skipping S3 (no data files found)")
        return

    def calc_props(rows, smi_key, max_n=500):
        mws, logps, hbds, hbas = [], [], [], []
        for r in rows[:max_n]:
            smi = r.get(smi_key, "")
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            mws.append(Descriptors.MolWt(mol))
            logps.append(Descriptors.MolLogP(mol))
            hbds.append(Descriptors.NumHDonors(mol))
            hbas.append(Descriptors.NumHAcceptors(mol))
        return {"MW": mws, "LogP": logps, "HBD": hbds, "HBA": hbas}

    act_props = calc_props(actives_rows, act_smi_key)
    dec_props = calc_props(decoys_rows, dec_smi_key)

    if not act_props["MW"] or not dec_props["MW"]:
        print("  Skipping S3 (could not parse SMILES)")
        return

    fig, axes = plt.subplots(2, 2, figsize=(5.5, 4.5))
    props = list(act_props.keys())
    units = {"MW": "Da", "LogP": "", "HBD": "", "HBA": ""}
    for ax, prop in zip(axes.flatten(), props):
        ax.hist(act_props[prop], bins=20, color=PALETTE["skill"], alpha=0.6,
                density=True, label="Actives")
        ax.hist(dec_props[prop], bins=20, color=PALETTE["naive"], alpha=0.6,
                density=True, label="Decoys")
        ax.set_xlabel(f"{prop} {units[prop]}" if units[prop] else prop)
        ax.set_ylabel("Density")
        ax.legend(fontsize=6)
    fig.suptitle(f"Property Distributions: Actives vs. Decoys -- {target_name}", fontsize=8)
    plt.tight_layout()
    save_figure(fig, "figS3_property_qc", fig_dir)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure S4: Score rank plots
# ---------------------------------------------------------------------------
def make_fig_s4(data, fig_dir, target_name):
    fig, axes = plt.subplots(1, 2, figsize=(6.0, 3.0))
    for ax, (protocol, d) in zip(axes, data.items()):
        # Sort by descending score (best first)
        order = np.argsort(-d["y_score"])
        y_sorted = d["y_true"][order]
        cumact = np.cumsum(y_sorted)
        total = y_sorted.sum()
        rank_frac = np.arange(1, len(y_sorted) + 1) / len(y_sorted)
        cum_frac = cumact / total if total > 0 else cumact

        ax.plot(rank_frac, cum_frac, color=PALETTE[protocol], lw=1.2,
                label=PROTOCOL_LABELS[protocol])
        ax.plot([0, 1], [0, 1], "--", color=PALETTE["grey"], lw=0.7, label="Random")
        ax.set_xlabel("Fraction of Database Screened")
        ax.set_ylabel("Fraction of Actives Found")
        ax.set_title(f"{PROTOCOL_LABELS[protocol]}")
        ax.legend(fontsize=6)

    fig.suptitle(f"Score Rank -- {target_name}", fontsize=8, y=1.02)
    plt.tight_layout()
    save_figure(fig, "figS4_score_rank", fig_dir)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure S5: Metrics heatmap
# ---------------------------------------------------------------------------
def make_fig_s5(metrics, fig_dir, target_name):
    try:
        import seaborn as sns
    except ImportError:
        print("  Skipping S5 (seaborn not available)")
        return

    metric_keys = ["roc_auc", "bedroc_20", "bedroc_80", "ef_1pct", "ef_5pct",
                   "ef_10pct", "log_auc", "aupr"]
    metric_labels = ["ROC AUC", "BEDROC a=20", "BEDROC a=80", "EF 1%", "EF 5%",
                     "EF 10%", "LogAUC", "AUPR"]
    protocols = list(metrics.keys())

    matrix = []
    for p in protocols:
        row = [metrics[p].get(k, np.nan) for k in metric_keys]
        matrix.append(row)
    matrix = np.array(matrix, dtype=float)

    # Normalize each column to [0,1] for heatmap coloring
    col_min = np.nanmin(matrix, axis=0)
    col_max = np.nanmax(matrix, axis=0)
    matrix_norm = (matrix - col_min) / np.maximum(col_max - col_min, 1e-10)

    fig, ax = plt.subplots(figsize=(6.0, 2.0))
    im = ax.imshow(matrix_norm, aspect="auto", cmap="RdYlGn",
                   vmin=0, vmax=1, interpolation="nearest")
    ax.set_xticks(range(len(metric_labels)))
    ax.set_xticklabels(metric_labels, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(protocols)))
    ax.set_yticklabels([PROTOCOL_LABELS.get(p, p) for p in protocols])
    for i in range(len(protocols)):
        for j in range(len(metric_keys)):
            val = matrix[i, j]
            if not np.isnan(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=6.5, color="k")
    plt.colorbar(im, ax=ax, label="Relative performance\n(min-max normalized)")
    ax.set_title(f"All Metrics -- {target_name}")
    plt.tight_layout()
    save_figure(fig, "figS5_metrics_heatmap", fig_dir)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Generate all figures for a single-target VS benchmark."
    )
    parser.add_argument("--target-dir", required=True,
                        help="Path to target directory (e.g. targets/egfr/)")
    parser.add_argument("--target-name", default=None,
                        help="Human-readable target name for titles (default: extracted from dir)")
    args = parser.parse_args()

    target_dir = Path(args.target_dir).resolve()
    results_dir = target_dir / "05_evaluation" / "results"
    fig_dir = target_dir / "05_evaluation" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Determine target name
    if args.target_name:
        target_name = args.target_name
    else:
        target_name = target_dir.name.upper().replace("_", " ")

    set_nature_style()

    print(f"Target: {target_name}")
    print(f"Target dir: {target_dir}")
    print("Loading data...")
    data = load_data(target_dir)

    # Check if scores are available
    has_data = any(d["y_score"].any() for d in data.values())
    if not has_data:
        print("WARNING: No score data found. Generating placeholder figures.")

    metrics = {p: load_metrics(results_dir, p) for p in data}

    print("Generating figures...")
    print("  Figure 1: Main 3-panel composite")
    make_fig_main(data, metrics, results_dir, fig_dir, target_name)

    print("  Figure S1: Semi-log ROC")
    make_fig_s1(data, metrics, fig_dir, target_name)

    print("  Figure S2: Precision-Recall")
    make_fig_s2(data, metrics, fig_dir, target_name)

    print("  Figure S3: Property QC")
    make_fig_s3(target_dir, fig_dir, target_name)

    print("  Figure S4: Score rank plots")
    make_fig_s4(data, fig_dir, target_name)

    print("  Figure S5: Metrics heatmap")
    make_fig_s5(metrics, fig_dir, target_name)

    print(f"\nAll figures saved to: {fig_dir}")
    for f in sorted(fig_dir.iterdir()):
        print(f"  {f.name:40s} ({f.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
