#!/usr/bin/env python3
"""
Proxy filter ablation: estimate whether the AChE/COX-2 AUC gain attributed to
the skill protocol survives if PAINS/Brenk-filtered compounds are NOT excluded.

Method (proxy, not a true ablation):
  - For each target, take the skill-protocol scores on the post-filter library.
  - For PAINS/Brenk-filtered compounds, substitute the naive-protocol score
    (which exists, since the naive protocol did not filter).
  - Recompute ROC AUC on the resulting *full library*.
  - Compare with naive full-library AUC and skill post-filter AUC.

Limitation: the substituted scores were generated with naive docking parameters
(20 Å box, exhaustiveness 8), not skill parameters. A true ablation requires
re-docking the filtered compounds with skill parameters. This script gives a
directional estimate only and is intended to inform whether the planned full
ablation is worth running.

Outputs:
  - filter_ablation.csv  : per-target AUCs (naive_full, skill_post, skill_proxy_full)
  - filter_ablation.png  : bar chart comparing the three AUCs per target
  - filter_ablation_report.md : human-readable summary
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import numpy as np

try:
    from sklearn.metrics import roc_auc_score
except ImportError:
    print("ERROR: scikit-learn required. Install with: mamba install scikit-learn", file=sys.stderr)
    sys.exit(1)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC_ROOT = Path("/data/CLAUDE_works/TEST_SKILLS/vina_docking_test/vs_benchmark/vs_benchmark_multitarget")
OUT_DIR = Path(__file__).resolve().parent.parent

TARGETS = [
    "ache",
    "bace1",
    "cdk2",
    "cox2",
    "dpp4",
    "egfr",
    "esr1",
    "fpr2",
    "hsp90",
    "p38",
    "thrombin",
]


def load_scores(target: str, protocol: str) -> dict[str, float]:
    p = SRC_ROOT / "targets" / target / "04_docking" / protocol / f"scores_{protocol}.csv"
    out = {}
    with p.open() as f:
        for row in csv.DictReader(f):
            try:
                s = float(row["score"])
            except (ValueError, KeyError):
                continue
            out[row["compound_id"]] = s
    return out


def load_labels(target: str) -> dict[str, str]:
    p = SRC_ROOT / "targets" / target / "03_library_preparation" / "library_labels.csv"
    out = {}
    with p.open() as f:
        for row in csv.DictReader(f):
            out[row["compound_id"]] = row["label"]
    return out


def load_filtered(target: str) -> set[str]:
    p = SRC_ROOT / "targets" / target / "03_library_preparation" / "skill" / "pdbqt" / "filter_report.csv"
    if not p.exists():
        return set()
    with p.open() as f:
        return {row["compound_id"] for row in csv.DictReader(f)}


def auc_from_scores(scores: dict[str, float], labels: dict[str, str]) -> tuple[float, int, int]:
    """Compute ROC AUC. Lower (more negative) score = better rank for actives.

    Filters out compounds with score = 0.0 (Uni-Dock failure / not docked) since
    treating these as worst-rank biases the AUC. Reports n_actives and n_decoys.
    """
    y_true, y_score = [], []
    for cid, s in scores.items():
        lab = labels.get(cid)
        if lab not in ("active", "decoy"):
            continue
        if s == 0.0 or math.isnan(s) or s > 0:
            continue
        y_true.append(1 if lab == "active" else 0)
        y_score.append(-s)  # invert: more negative score → higher rank
    n_act = sum(y_true)
    n_dec = len(y_true) - n_act
    if n_act < 2 or n_dec < 2:
        return float("nan"), n_act, n_dec
    return roc_auc_score(y_true, y_score), n_act, n_dec


def main() -> int:
    rows = []
    print(f"{'target':10s} {'naive_AUC':>10s} {'skill_AUC':>10s} {'proxy_AUC':>10s} {'ΔAUC_proxy':>11s} {'ΔAUC_post':>10s}")
    print("-" * 65)

    for t in TARGETS:
        naive = load_scores(t, "naive")
        skill = load_scores(t, "skill")
        labels = load_labels(t)
        filtered = load_filtered(t)

        # Proxy: use skill score where available; for PAINS/Brenk-filtered compounds,
        # OVERWRITE the skill-protocol's 0.0 (filtered = not docked) with the naive score
        # where the naive protocol successfully docked the compound. Compounds missing
        # from naive (or with naive score = 0.0) remain at 0.0 and are excluded by the
        # AUC step.
        proxy = dict(skill)
        n_substituted = 0
        for cid in filtered:
            naive_s = naive.get(cid)
            if naive_s is not None and naive_s < 0.0:  # valid naive docking score
                proxy[cid] = naive_s
                n_substituted += 1

        auc_n, _, _ = auc_from_scores(naive, labels)
        auc_s, _, _ = auc_from_scores(skill, labels)
        auc_p, np_act, np_dec = auc_from_scores(proxy, labels)
        d_proxy = auc_p - auc_n if not (math.isnan(auc_p) or math.isnan(auc_n)) else float("nan")
        d_post = auc_s - auc_n if not (math.isnan(auc_s) or math.isnan(auc_n)) else float("nan")

        rows.append(
            {
                "target": t,
                "naive_auc_full": round(auc_n, 4) if not math.isnan(auc_n) else None,
                "skill_auc_post": round(auc_s, 4) if not math.isnan(auc_s) else None,
                "skill_auc_proxy_full": round(auc_p, 4) if not math.isnan(auc_p) else None,
                "delta_proxy_minus_naive": round(d_proxy, 4) if not math.isnan(d_proxy) else None,
                "delta_post_minus_naive": round(d_post, 4) if not math.isnan(d_post) else None,
                "n_actives_proxy": np_act,
                "n_decoys_proxy": np_dec,
            }
        )
        fmt = lambda x: f"{x:.4f}" if not math.isnan(x) else "  NA  "
        print(f"{t:10s} {fmt(auc_n):>10s} {fmt(auc_s):>10s} {fmt(auc_p):>10s} {fmt(d_proxy):>11s} {fmt(d_post):>10s}")

    # Save CSV
    csv_path = OUT_DIR / "filter_ablation.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {csv_path}")

    # Plot: bar chart per target with three AUCs side-by-side
    fig, ax = plt.subplots(figsize=(9, 4.5))
    targets = [r["target"] for r in rows]
    naive_aucs = [r["naive_auc_full"] or 0 for r in rows]
    skill_aucs = [r["skill_auc_post"] or 0 for r in rows]
    proxy_aucs = [r["skill_auc_proxy_full"] or 0 for r in rows]
    x = np.arange(len(targets))
    w = 0.27
    ax.bar(x - w, naive_aucs, w, label="Naive (full library)", color="#888")
    ax.bar(x, skill_aucs, w, label="Skill (post-filter)", color="#3a7")
    ax.bar(x + w, proxy_aucs, w, label="Skill + naive scores for filtered (proxy full)", color="#c63")
    ax.axhline(0.5, color="k", linestyle=":", linewidth=0.7, alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(targets, rotation=30, ha="right")
    ax.set_ylabel("ROC AUC")
    ax.set_ylim(0.4, 0.8)
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
    ax.set_title("Proxy filter ablation — does the AChE/COX-2 gain survive when filtered compounds are scored?")
    fig.tight_layout()
    png_path = OUT_DIR / "figures" / "filter_ablation_proxy.png"
    pdf_path = OUT_DIR / "figures" / "filter_ablation_proxy.pdf"
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=200)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"Wrote {png_path}")
    print(f"Wrote {pdf_path}")

    # Markdown report
    md = OUT_DIR / "filter_ablation_report.md"
    with md.open("w") as f:
        f.write("# Proxy Filter-Ablation Report\n\n")
        f.write(
            "**Method:** For each target, replaced PAINS/Brenk-filtered compounds' missing skill scores with their naive-protocol scores, then recomputed ROC AUC on the full library. This is a proxy for a true ablation — the substituted scores were generated under naive docking parameters, not skill — but is informative about whether the AChE/COX-2 AUC gain is filter-driven (active reshaping) or docking-driven.\n\n"
        )
        f.write("**Interpretation:**\n\n")
        f.write("- If `delta_proxy_minus_naive` is close to zero, the skill AUC gain *was* primarily driven by filtering out hard-to-dock actives.\n")
        f.write("- If `delta_proxy_minus_naive` is similar to `delta_post_minus_naive`, the gain *survives* filter ablation and is at least partly attributable to skill docking parameters.\n\n")
        f.write("| Target | Naive AUC (full) | Skill AUC (post-filter) | Proxy full AUC | Δ proxy − naive | Δ post − naive |\n")
        f.write("|---|---:|---:|---:|---:|---:|\n")
        for r in rows:
            def fm(v):
                return f"{v:.3f}" if isinstance(v, (int, float)) else "NA"
            f.write(
                f"| {r['target']} | {fm(r['naive_auc_full'])} | {fm(r['skill_auc_post'])} | {fm(r['skill_auc_proxy_full'])} | {fm(r['delta_proxy_minus_naive'])} | {fm(r['delta_post_minus_naive'])} |\n"
            )
        f.write("\n")
        # Spotlight on AChE and COX-2
        for spot in ("ache", "cox2"):
            r = next(x for x in rows if x["target"] == spot)
            f.write(f"## {spot.upper()} spotlight\n\n")
            f.write(f"- Skill post-filter Δ vs naive: **{r['delta_post_minus_naive']:+.3f}** (this was the original BIB-submitted gain)\n")
            f.write(f"- Skill proxy full-library Δ vs naive: **{r['delta_proxy_minus_naive']:+.3f}**\n")
            d_post = r["delta_post_minus_naive"] or 0.0
            d_proxy = r["delta_proxy_minus_naive"] or 0.0
            shrink = (1 - abs(d_proxy) / abs(d_post)) * 100 if d_post != 0 else 0
            direction = "shrinks toward zero" if abs(d_proxy) < abs(d_post) else "is preserved or grows"
            f.write(f"- Effect {direction} by {abs(shrink):.0f}% under proxy ablation.\n\n")
        f.write("**Caveat:** This is a *proxy* analysis. The substituted scores were generated under naive docking parameters (20 Å box, exhaustiveness 8), so any preserved Δ in the proxy may underestimate the true ablation effect. A full ablation — re-docking the PAINS/Brenk-filtered compounds with skill docking parameters — remains required for a definitive answer. This proxy is intended only to inform whether the full ablation is worth running (yes if the proxy shows substantial preservation; less informative if the proxy already collapses the effect to zero).\n")
    print(f"Wrote {md}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
