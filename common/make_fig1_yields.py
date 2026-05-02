#!/usr/bin/env python3
"""Build the corrected Fig. 1 (yield panels with three denominators).

Panel A: stacked bar per target showing
  - naive successfully docked (% of full library)
  - skill PAINS/Brenk filtered out (% of full library)
  - skill successfully docked (% of full library)
  - skill prep failures (% of full library, ~0.5%, usually invisible)

Panel B: proxy filter-ablation results — three AUCs side-by-side per target:
  - naive (full library)
  - skill (post-filter)
  - skill + naive scores for filtered (proxy full library)

Inputs: yields.csv (from recompute_yields.py) and filter_ablation.csv (from
proxy_filter_ablation.py). Run those first.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def load_csv(p: Path) -> list[dict]:
    with p.open() as f:
        return list(csv.DictReader(f))


def main() -> int:
    yields = load_csv(ROOT / "yields.csv")
    abl = {r["target"]: r for r in load_csv(ROOT / "filter_ablation.csv")}

    # Order: kinases together, proteases together, etc.
    order = ["fpr2", "egfr", "cdk2", "p38", "cox2", "esr1", "dpp4", "ache", "bace1", "hsp90", "thrombin"]
    yields = sorted(yields, key=lambda r: order.index(r["target"]))

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 4.5), gridspec_kw={"width_ratios": [1.1, 1]})

    # ------------------ Panel A: yield stack ------------------
    targets = [r["target"].upper() if r["target"] != "fpr2" else "FPR2" for r in yields]
    naive_doc = np.array([float(r["naive_full_pct"]) for r in yields])
    skill_doc = np.array([float(r["skill_full_pct"]) for r in yields])
    skill_filt = np.array([float(r["filter_lib_pct"]) for r in yields])
    # Skill bar = doc + filtered + (failures = 100 - doc - filtered)
    skill_fail = 100.0 - skill_doc - skill_filt
    skill_fail = np.clip(skill_fail, 0, None)

    x = np.arange(len(targets))
    w = 0.36

    # Naive bar — single colour, two-component (docked vs failed)
    naive_fail = 100.0 - naive_doc
    axA.bar(x - w / 2, naive_doc, w, label="Naive: docked", color="#2c7fb8")
    axA.bar(x - w / 2, naive_fail, w, bottom=naive_doc, label="Naive: prep/dock failure", color="#bdd7e7", edgecolor="white", linewidth=0.5)

    # Skill bar — three components (docked, filtered, true failures)
    axA.bar(x + w / 2, skill_doc, w, label="Skill: docked", color="#31a354")
    axA.bar(x + w / 2, skill_filt, w, bottom=skill_doc, label="Skill: PAINS/Brenk filtered", color="#fdae6b", edgecolor="white", linewidth=0.5)
    axA.bar(x + w / 2, skill_fail, w, bottom=skill_doc + skill_filt, label="Skill: prep failure", color="#bdbdbd", edgecolor="white", linewidth=0.5)

    axA.set_xticks(x)
    axA.set_xticklabels(targets, rotation=30, ha="right", fontsize=9)
    axA.set_ylabel("Fraction of full library (%)")
    axA.set_ylim(0, 100)
    axA.set_title("(A) Pipeline yield: where the library goes", fontsize=11, loc="left")
    axA.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3, fontsize=8, frameon=False)
    axA.grid(axis="y", alpha=0.25, linewidth=0.5)
    axA.set_axisbelow(True)
    for spine in ["top", "right"]:
        axA.spines[spine].set_visible(False)

    # ------------------ Panel B: AUC ablation ------------------
    auc_naive = np.array([float(abl[r["target"]]["naive_auc_full"]) for r in yields])
    auc_skill_post = np.array([float(abl[r["target"]]["skill_auc_post"]) for r in yields])
    auc_proxy = np.array([float(abl[r["target"]]["skill_auc_proxy_full"]) for r in yields])

    w2 = 0.27
    axB.bar(x - w2, auc_naive, w2, label="Naive (full library)", color="#888888")
    axB.bar(x, auc_skill_post, w2, label="Skill (post-filter)", color="#31a354")
    axB.bar(x + w2, auc_proxy, w2, label="Skill + naive-fallback for filtered (proxy full)", color="#e6550d")
    axB.axhline(0.5, color="k", linestyle=":", linewidth=0.7, alpha=0.5)

    axB.set_xticks(x)
    axB.set_xticklabels(targets, rotation=30, ha="right", fontsize=9)
    axB.set_ylabel("ROC AUC")
    axB.set_ylim(0.4, 0.78)
    axB.set_title("(B) Proxy filter ablation: most of the COX-2 / AChE gain is filter-driven", fontsize=11, loc="left")
    axB.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize=8, frameon=False)
    axB.grid(axis="y", alpha=0.25, linewidth=0.5)
    axB.set_axisbelow(True)
    for spine in ["top", "right"]:
        axB.spines[spine].set_visible(False)

    # Annotate AChE and COX-2 with delta values
    for tgt in ("cox2", "ache"):
        i = order.index(tgt)
        d_post = float(abl[tgt]["delta_post_minus_naive"])
        d_proxy = float(abl[tgt]["delta_proxy_minus_naive"])
        axB.annotate(
            f"Δ post: {d_post:+.3f}\nΔ proxy: {d_proxy:+.3f}",
            xy=(i, max(auc_skill_post[i], auc_proxy[i]) + 0.015),
            ha="center",
            fontsize=7,
            color="#333",
        )

    fig.tight_layout()
    png_path = FIG_DIR / "fig1_yields_and_ablation.png"
    pdf_path = FIG_DIR / "fig1_yields_and_ablation.pdf"
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {png_path}")
    print(f"Wrote {pdf_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
