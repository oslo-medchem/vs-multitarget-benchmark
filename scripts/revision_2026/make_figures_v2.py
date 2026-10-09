#!/usr/bin/env python3
"""Figures for the September 2026 MLST rewrite.

Figure 1  Both endpoints side by side: conditional ranking (intersection) and
          complete workflow (end to end), each as a forest plot with the
          DerSimonian-Laird pooled estimate.
Figure 2  Component decomposition: preparation and parameters with the filter
          removed from both arms, and the effect of the filter itself.
Figure 3  Matched-input throughput: preparation success and docking yield on
          the identical molecule set.

All inputs come from data/revision_2026/. Run from the project root.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "revision_2026"
OUTDIR = ROOT / "figures" / "main"
OUTDIR.mkdir(parents=True, exist_ok=True)

PRETTY = {"fpr2": "FPR2", "egfr": "EGFR", "cdk2": "CDK2", "cox2": "COX-2",
          "esr1": "ESR1", "dpp4": "DPP-4", "ache": "AChE", "bace1": "BACE-1",
          "hsp90": "HSP90α", "p38": "p38 MAPK", "thrombin": "Thrombin"}
ORDER = ["fpr2", "egfr", "cdk2", "cox2", "esr1", "dpp4", "ache",
         "bace1", "hsp90", "p38", "thrombin"]

NAVY = "#1f4e79"
RUST = "#a63603"
GREY = "#6e6e6e"
LIGHT = "#bdbdbd"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 8,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "savefig.dpi": 400,
})


def forest(ax, targets, deltas, los, his, padj, pooled, title, xlabel):
    y = np.arange(len(targets))[::-1]
    for yi, d, lo, hi, pa in zip(y, deltas, los, his, padj):
        sig = pa < 0.05
        c = RUST if sig else GREY
        ax.plot([lo, hi], [yi, yi], color=c, lw=1.2, solid_capstyle="round")
        ax.plot([d], [yi], marker="s", ms=4.5, color=c,
                markeredgecolor="white", markeredgewidth=0.4)
    ax.axvline(0, color=LIGHT, lw=0.8, zorder=0)

    # pooled diamond
    yp = -1.4
    mu, plo, phi = pooled["pooled"], pooled["ci_lo"], pooled["ci_hi"]
    ax.add_patch(plt.Polygon([[plo, yp], [mu, yp + 0.32], [phi, yp],
                              [mu, yp - 0.32]], closed=True,
                             facecolor=NAVY, edgecolor=NAVY))
    ax.set_yticks(list(y) + [yp])
    ax.set_yticklabels([PRETTY.get(t, t) for t in targets] + ["Pooled"], fontsize=7.5)
    for lbl, pa in zip(ax.get_yticklabels()[:-1], padj):
        if pa < 0.05:
            lbl.set_fontweight("bold")
    ax.get_yticklabels()[-1].set_fontweight("bold")
    ax.set_ylim(yp - 2.7, len(targets) - 0.4)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_title(title, fontsize=8.5, fontweight="bold", loc="left", pad=6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7.5)
    txt = (f"pooled {mu:+.4f}\n95% CI [{plo:+.3f}, {phi:+.3f}]\n"
           f"p = {pooled['p']:.2f},  I² = {pooled['I2']:.0f}%")
    ax.text(0.98, 0.02, txt, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=6.8, color=NAVY,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=LIGHT, lw=0.5))


def figure1():
    i = pd.read_csv(DATA / "intersection_delong.csv").set_index("target")
    e = pd.read_csv(DATA / "endtoend_delong.csv").set_index("target")
    mi = json.loads((DATA / "meta_analysis_intersection.json").read_text())
    me = json.loads((DATA / "meta_analysis_endtoend.json").read_text())

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.4))
    forest(axes[0], ORDER,
           [i.loc[t, "delong_delta"] for t in ORDER],
           [i.loc[t, "ci_lo"] for t in ORDER],
           [i.loc[t, "ci_hi"] for t in ORDER],
           [i.loc[t, "delong_p_adj_BH"] for t in ORDER], mi,
           "A  Conditional ranking (intersection)", "ΔAUC (skill − naive)")
    forest(axes[1], ORDER,
           [e.loc[t, "delong_delta"] for t in ORDER],
           [e.loc[t, "ci_lo"] for t in ORDER],
           [e.loc[t, "ci_hi"] for t in ORDER],
           [e.loc[t, "delong_p_adj_BH"] for t in ORDER], me,
           "B  Complete workflow (end to end)", "ΔAUC (skill − naive)")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUTDIR / f"fig1_two_endpoints.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("wrote fig1_two_endpoints")


def figure2():
    d = pd.read_csv(DATA / "decomposition.csv").set_index("target")
    mc = json.loads((DATA / "meta_analysis_components.json").read_text())
    tg = [t for t in ORDER if t in d.index]

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.0))
    forest(axes[0], tg,
           [d.loc[t, "prep_params_delta_C_minus_A"] for t in tg],
           [d.loc[t, "prep_params_ci_lo"] for t in tg],
           [d.loc[t, "prep_params_ci_hi"] for t in tg],
           [d.loc[t, "prep_params_p_adj_BH"] for t in tg],
           mc["prep_plus_params_C_minus_A"],
           "A  Preparation and parameters\n     (filter removed from both arms)",
           "ΔAUC (skill − naive)")

    ax = axes[1]
    y = np.arange(len(tg))[::-1]
    vals = [d.loc[t, "filter_delta_C_minus_B"] for t in tg]
    padj = [d.loc[t, "filter_p_adj_BH"] for t in tg]
    cols = [(NAVY if v > 0 else RUST) if p < 0.05 else LIGHT
            for v, p in zip(vals, padj)]
    ax.barh(y, vals, color=cols, height=0.62, edgecolor="white", linewidth=0.4)
    ax.axvline(0, color=GREY, lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([PRETTY[t] for t in tg], fontsize=7.5)
    ax.set_ylim(-2.7, len(tg) - 0.4)
    for lbl, p in zip(ax.get_yticklabels(), padj):
        if p < 0.05:
            lbl.set_fontweight("bold")
    ax.set_xlabel("ΔAUC from removing the filter (end to end)", fontsize=8)
    ax.set_title("B  Effect of the substructure filter", fontsize=8.5,
                 fontweight="bold", loc="left", pad=6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7.5)
    ax.text(0.98, 0.02,
            "blue: removing the filter helps\nred: removing the filter hurts\n"
            f"I² = {mc['filter_C_minus_B']['I2']:.0f}%",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.8,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=LIGHT, lw=0.5))
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUTDIR / f"fig2_decomposition.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("wrote fig2_decomposition")


def figure3():
    p = pd.read_csv(DATA / "prep_matched.csv").set_index("target")
    tg = [t for t in ORDER if t in p.index]
    x = np.arange(len(tg))
    w = 0.38

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 4.6), sharex=True)
    ax = axes[0]
    ax.bar(x - w / 2, [p.loc[t, "naive_prep_success_pct"] for t in tg], w,
           label="Naive", color=LIGHT, edgecolor="white", linewidth=0.5)
    ax.bar(x + w / 2, [p.loc[t, "skill_prep_success_pct"] for t in tg], w,
           label="Skill-guided", color=NAVY, edgecolor="white", linewidth=0.5)
    ax.set_ylim(95, 100.4)
    ax.set_ylabel("Preparation success (%)", fontsize=8)
    ax.set_title("A  Producing a valid ligand file from the same molecules",
                 fontsize=8.5, fontweight="bold", loc="left", pad=6)
    ax.legend(frameon=False, fontsize=7.5, ncol=2, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    ax.bar(x - w / 2, [p.loc[t, "naive_dock_yield_pct"] for t in tg], w,
           label="Naive", color=LIGHT, edgecolor="white", linewidth=0.5)
    ax.bar(x + w / 2, [p.loc[t, "skill_dock_yield_pct"] for t in tg], w,
           label="Skill-guided", color=NAVY, edgecolor="white", linewidth=0.5)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Yielded a docking score (%)", fontsize=8)
    ax.set_title("B  Ending with a usable docking score from the same molecules",
                 fontsize=8.5, fontweight="bold", loc="left", pad=6)
    ax.set_xticks(x)
    ax.set_xticklabels([PRETTY[t] for t in tg], rotation=30, ha="right", fontsize=7.5)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUTDIR / f"fig3_matched_throughput.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("wrote fig3_matched_throughput")


if __name__ == "__main__":
    figure1()
    figure2()
    figure3()
