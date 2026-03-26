#!/usr/bin/env python3
"""
Compute all VS evaluation metrics for naive and skill-guided protocols.

Usage:
    python compute_metrics.py --target-dir targets/egfr/

Metrics: ROC AUC, BEDROC(alpha=20,80), EF1/5/10%, LogAUC, pROC AUC (1%, 5%),
         AUPR, RIE(alpha=20).
Bootstrap 95% CI (BCa, n=10,000) for all metrics.

Outputs to target_dir/05_evaluation/results/:
    metrics_naive.json, metrics_skill.json
    roc_naive.json, roc_skill.json
    y_true_{naive,skill}.npy, y_score_{naive,skill}.npy
"""
import argparse
import csv
import json
import sys
import numpy as np
from pathlib import Path
from scipy.stats import bootstrap as scipy_bootstrap
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve

# Try RDKit ML scoring
try:
    from rdkit.ML.Scoring.Scoring import CalcAUC, CalcBEDROC, CalcEnrichment, CalcRIE
    RDKIT_SCORING = True
except ImportError:
    RDKIT_SCORING = False
    print("WARNING: RDKit ML Scoring not available, using sklearn fallback")

N_BOOTSTRAP = 10_000
CI_LEVEL = 0.95


def load_data(labels_csv, scores_csv):
    """Load labels and scores, return aligned arrays."""
    labels = {}
    for r in csv.DictReader(open(labels_csv)):
        labels[r["compound_id"]] = 1 if r["label"] == "active" else 0

    scores = {}
    for r in csv.DictReader(open(scores_csv)):
        try:
            scores[r["compound_id"]] = float(r["score"])
        except (ValueError, KeyError):
            scores[r["compound_id"]] = 0.0

    # Align
    cids = sorted(labels.keys())
    y_true = np.array([labels[c] for c in cids])
    y_score = np.array([scores.get(c, 0.0) for c in cids])

    # For docking: more negative = better. Negate for rank (higher = better).
    # Treat 0.0 (failed docking) as worst score.
    y_score_rank = -y_score  # negate: less negative (0) -> worst rank

    return cids, y_true, y_score, y_score_rank


def load_data_docked_only(labels_csv, scores_csv):
    """Load only compounds that were successfully docked (score != 0).

    This removes the confound where failed PDBQT preparation or PAINS
    filtering assigns worst-rank scores to compounds.
    """
    labels = {}
    for r in csv.DictReader(open(labels_csv)):
        labels[r["compound_id"]] = 1 if r["label"] == "active" else 0

    scores = {}
    for r in csv.DictReader(open(scores_csv)):
        try:
            s = float(r["score"])
            if s != 0.0:
                scores[r["compound_id"]] = s
        except (ValueError, KeyError):
            pass

    # Only keep docked compounds
    cids = sorted(c for c in labels if c in scores)
    y_true = np.array([labels[c] for c in cids])
    y_score = np.array([scores[c] for c in cids])
    y_score_rank = -y_score

    return cids, y_true, y_score, y_score_rank


def rdkit_sorted_list(y_true, y_score_rank):
    """Build RDKit-style sorted list: [(score, label)] sorted descending."""
    pairs = sorted(zip(y_score_rank, y_true), reverse=True)
    return [[sc, lab] for sc, lab in pairs]


def calc_logauc(y_true, y_score_rank, fpr_min=0.001, fpr_max=1.0):
    """Calculate LogAUC (log-scale ROC AUC)."""
    fpr, tpr, _ = roc_curve(y_true, y_score_rank)
    # Interpolate on log scale
    fpr_log = np.logspace(np.log10(fpr_min), np.log10(fpr_max), 1000)
    tpr_interp = np.interp(fpr_log, fpr, tpr)
    # AUC in log space
    trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
    log_auc = trapz(tpr_interp, np.log10(fpr_log)) / (np.log10(fpr_max) - np.log10(fpr_min))
    return float(log_auc)


def calc_partial_roc_auc(y_true, y_score_rank, max_fpr):
    """Partial ROC AUC normalized to [0, 1] range."""
    try:
        pauc = roc_auc_score(y_true, y_score_rank, max_fpr=max_fpr)
        # sklearn returns the standardized form -- normalize to [0,1]
        pauc_norm = pauc / max_fpr
        return float(pauc_norm)
    except Exception:
        return 0.5


def bedroc_sklearn(y_true, y_score_rank, alpha=20.0):
    """BEDROC calculation from sorted arrays (RDKit-independent)."""
    n = len(y_true)
    n_actives = y_true.sum()
    if n_actives == 0 or n_actives == n:
        return 0.5
    ra = n_actives / n
    order = np.argsort(-y_score_rank)
    y_sorted = y_true[order]

    ranks = np.where(y_sorted == 1)[0] + 1  # 1-indexed
    s = np.sum(np.exp(-alpha * ranks / n))
    s_bar = (1 - np.exp(-alpha * ra)) / (np.exp(alpha / n) - 1)
    ri = np.exp(-alpha * ra) - np.exp(-alpha)
    d_val = np.exp(alpha / n) - 1

    bedroc = (s * d_val * np.sinh(alpha / 2) /
              (np.cosh(alpha / 2) - np.cosh(alpha / 2 - alpha * ra))) * ra + ri / (1 - np.exp(alpha * (1 - ra)))
    # Normalize
    min_bedroc = (1 - np.exp(alpha * ra)) / (1 - np.exp(alpha)) * ra
    max_bedroc = (1 - np.exp(-alpha * ra)) / (1 - np.exp(-alpha)) * ra
    bedroc_norm = (bedroc - min_bedroc) / (max_bedroc - min_bedroc) if (max_bedroc - min_bedroc) > 0 else 0.5
    return float(bedroc_norm)


def calc_ef(y_true, y_score_rank, fraction):
    """Enrichment Factor at given fraction of the library."""
    n = len(y_true)
    n_actives = y_true.sum()
    n_top = max(1, int(n * fraction))
    order = np.argsort(-y_score_rank)
    y_sorted = y_true[order]
    hits_in_top = y_sorted[:n_top].sum()
    expected = n_actives * fraction
    return float(hits_in_top / expected) if expected > 0 else 0.0


def calc_rie(y_true, y_score_rank, alpha=20.0):
    """Robust Initial Enhancement."""
    n = len(y_true)
    n_actives = y_true.sum()
    if n_actives == 0 or n_actives == n:
        return 1.0
    ra = n_actives / n
    order = np.argsort(-y_score_rank)
    y_sorted = y_true[order]
    ranks = np.where(y_sorted == 1)[0] + 1
    s = np.sum(np.exp(-alpha * ranks / n))
    rie_raw = s / (n * ra) * (1 - np.exp(-alpha)) / (np.exp(alpha / n) - 1)
    return float(rie_raw)


def compute_all_metrics(y_true, y_score_rank):
    """Compute all metrics for one protocol."""
    n_actives = int(y_true.sum())
    n_total = len(y_true)

    # RDKit-based metrics
    rdkit_list = rdkit_sorted_list(y_true, y_score_rank)

    if RDKIT_SCORING:
        roc_auc = float(CalcAUC(rdkit_list, 1))
        bedroc_20 = float(CalcBEDROC(rdkit_list, 1, 20.0))
        bedroc_80 = float(CalcBEDROC(rdkit_list, 1, 80.0))
        ef1 = float(CalcEnrichment(rdkit_list, 1, [0.01])[0])
        ef5 = float(CalcEnrichment(rdkit_list, 1, [0.05])[0])
        ef10 = float(CalcEnrichment(rdkit_list, 1, [0.10])[0])
        rie_20 = float(CalcRIE(rdkit_list, 1, 20.0))
    else:
        # sklearn fallback
        roc_auc = float(roc_auc_score(y_true, y_score_rank))
        bedroc_20 = bedroc_sklearn(y_true, y_score_rank, alpha=20.0)
        bedroc_80 = bedroc_sklearn(y_true, y_score_rank, alpha=80.0)
        ef1 = calc_ef(y_true, y_score_rank, 0.01)
        ef5 = calc_ef(y_true, y_score_rank, 0.05)
        ef10 = calc_ef(y_true, y_score_rank, 0.10)
        rie_20 = calc_rie(y_true, y_score_rank, alpha=20.0)

    log_auc = calc_logauc(y_true, y_score_rank)
    proc_1pct = calc_partial_roc_auc(y_true, y_score_rank, 0.01)
    proc_5pct = calc_partial_roc_auc(y_true, y_score_rank, 0.05)
    aupr = float(average_precision_score(y_true, y_score_rank))

    return {
        "n_total": n_total,
        "n_actives": n_actives,
        "n_decoys": n_total - n_actives,
        "roc_auc": roc_auc,
        "bedroc_20": bedroc_20,
        "bedroc_80": bedroc_80,
        "ef_1pct": ef1,
        "ef_5pct": ef5,
        "ef_10pct": ef10,
        "rie_20": rie_20,
        "log_auc": log_auc,
        "proc_auc_1pct": proc_1pct,
        "proc_auc_5pct": proc_5pct,
        "aupr": aupr,
    }


def bootstrap_metric(y_true, y_score_rank, metric_fn, n_resamples=N_BOOTSTRAP):
    """BCa bootstrap CI for a metric."""
    def _fn(y_t, y_s):
        return metric_fn(y_t, y_s)

    try:
        result = scipy_bootstrap(
            (y_true, y_score_rank),
            lambda yt, ys: _fn(yt, ys),
            n_resamples=n_resamples,
            confidence_level=CI_LEVEL,
            method="BCa",
            vectorized=False,
            paired=True,
        )
        return {
            "ci_low": float(result.confidence_interval.low),
            "ci_high": float(result.confidence_interval.high),
        }
    except Exception as e:
        # Fallback: percentile bootstrap
        vals = []
        rng = np.random.default_rng(42)
        n = len(y_true)
        for _ in range(min(1000, n_resamples)):
            idx = rng.integers(0, n, n)
            yt_b = y_true[idx]
            ys_b = y_score_rank[idx]
            if yt_b.sum() > 0 and yt_b.sum() < len(yt_b):
                try:
                    vals.append(_fn(yt_b, ys_b))
                except Exception:
                    pass
        if vals:
            return {
                "ci_low": float(np.percentile(vals, 2.5)),
                "ci_high": float(np.percentile(vals, 97.5)),
            }
        return {"ci_low": None, "ci_high": None}


def get_roc_data(y_true, y_score_rank):
    """Return ROC curve points for plotting."""
    fpr, tpr, thresh = roc_curve(y_true, y_score_rank)
    return {
        "fpr": fpr.tolist(),
        "tpr": tpr.tolist(),
        "thresholds": thresh.tolist(),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Compute VS evaluation metrics for naive and skill protocols."
    )
    parser.add_argument("--target-dir", required=True,
                        help="Path to target directory (e.g. targets/egfr/)")
    args = parser.parse_args()

    target_dir = Path(args.target_dir).resolve()
    labels_csv = target_dir / "03_library_preparation" / "library_labels.csv"
    results_dir = target_dir / "05_evaluation" / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    if not labels_csv.exists():
        print(f"ERROR: {labels_csv} not found")
        sys.exit(1)

    protocols = {
        "naive": target_dir / "04_docking" / "naive" / "scores_naive.csv",
        "skill": target_dir / "04_docking" / "skill" / "scores_skill.csv",
    }

    print("Loading data...")
    labels_data = list(csv.DictReader(open(labels_csv)))
    labels = {r["compound_id"]: r["label"] for r in labels_data}

    metric_fns = {
        "roc_auc": lambda yt, ys: float(roc_auc_score(yt, ys)) if yt.sum() > 0 else 0.5,
        "bedroc_20": lambda yt, ys: bedroc_sklearn(yt, ys, 20.0),
        "bedroc_80": lambda yt, ys: bedroc_sklearn(yt, ys, 80.0),
        "ef_1pct": lambda yt, ys: calc_ef(yt, ys, 0.01),
        "ef_5pct": lambda yt, ys: calc_ef(yt, ys, 0.05),
        "log_auc": lambda yt, ys: calc_logauc(yt, ys),
        "aupr": lambda yt, ys: float(average_precision_score(yt, ys)) if yt.sum() > 0 else 0.0,
    }

    all_results = {}

    for protocol, scores_csv in protocols.items():
        print(f"\n=== {protocol.upper()} ===")

        if not scores_csv.exists():
            print(f"  WARNING: {scores_csv} not found, skipping")
            continue

        # --- Mode 1: Full library (undocked → worst rank) ---
        cids_full, y_true_full, y_score_full, y_rank_full = load_data(labels_csv, scores_csv)
        n_docked = int((y_score_full != 0.0).sum())
        print(f"  Full library: {len(cids_full)} compounds, {n_docked} docked")
        print(f"  Actives: {int(y_true_full.sum())}, Decoys: {int((1-y_true_full).sum())}")

        metrics_full = compute_all_metrics(y_true_full, y_rank_full)
        print(f"  [full]  ROC AUC: {metrics_full['roc_auc']:.4f}  BEDROC20: {metrics_full['bedroc_20']:.4f}  EF1%: {metrics_full['ef_1pct']:.2f}")

        # --- Mode 2: Docked-only (exclude undocked) ---
        cids_do, y_true_do, y_score_do, y_rank_do = load_data_docked_only(labels_csv, scores_csv)
        n_actives_do = int(y_true_do.sum())
        print(f"  Docked-only: {len(cids_do)} compounds ({n_actives_do} actives, {len(cids_do)-n_actives_do} decoys)")

        metrics_do = compute_all_metrics(y_true_do, y_rank_do)
        print(f"  [dock]  ROC AUC: {metrics_do['roc_auc']:.4f}  BEDROC20: {metrics_do['bedroc_20']:.4f}  EF1%: {metrics_do['ef_1pct']:.2f}")

        # Bootstrap CIs (on docked-only, the primary analysis)
        print("  Computing bootstrap CIs (docked-only, 10,000 resamples)...")
        bootstrap_cis = {}
        for mname, mfn in metric_fns.items():
            ci = bootstrap_metric(y_true_do, y_rank_do, mfn, N_BOOTSTRAP)
            bootstrap_cis[mname] = ci
            print(f"    {mname}: {metrics_do[mname]:.4f} [{ci['ci_low']:.4f}, {ci['ci_high']:.4f}]")

        # Save results (include both modes)
        output = {
            "protocol": protocol,
            "n_compounds": len(cids_full),
            "n_docked": n_docked,
            "n_actives_full": int(y_true_full.sum()),
            "n_actives_docked": n_actives_do,
            # Primary: docked-only metrics
            **{f"{k}": v for k, v in metrics_do.items()},
            "bootstrap_ci": bootstrap_cis,
            # Secondary: full-library metrics (with score=0 penalty)
            "full_library": metrics_full,
        }

        metrics_file = results_dir / f"metrics_{protocol}.json"
        with open(metrics_file, "w") as f:
            json.dump(output, f, indent=2)
        print(f"  Saved: {metrics_file}")

        # Save ROC curve data (docked-only)
        roc_data = get_roc_data(y_true_do, y_rank_do)
        roc_file = results_dir / f"roc_{protocol}.json"
        with open(roc_file, "w") as f:
            json.dump(roc_data, f, indent=2)

        # Save score arrays for statistical_analysis.py
        # Full-library arrays (for DeLong on same compound set)
        np.save(results_dir / f"y_true_{protocol}.npy", y_true_full)
        np.save(results_dir / f"y_score_{protocol}.npy", y_rank_full)
        # Docked-only arrays
        np.save(results_dir / f"y_true_{protocol}_docked.npy", y_true_do)
        np.save(results_dir / f"y_score_{protocol}_docked.npy", y_rank_do)
        # Compound IDs for intersection analysis
        np.save(results_dir / f"cids_{protocol}_docked.npy", np.array(cids_do))

        all_results[protocol] = output

    print("\n=== Summary (primary: docked-only) ===")
    for p, res in all_results.items():
        print(f"{p:10s}: AUC={res['roc_auc']:.4f} BEDROC20={res['bedroc_20']:.4f} EF1%={res['ef_1pct']:.2f} AUPR={res['aupr']:.4f} (n={res['n_docked']})")
    print("\n=== Summary (secondary: full-library with score=0 penalty) ===")
    for p, res in all_results.items():
        fl = res["full_library"]
        print(f"{p:10s}: AUC={fl['roc_auc']:.4f} BEDROC20={fl['bedroc_20']:.4f} EF1%={fl['ef_1pct']:.2f} AUPR={fl['aupr']:.4f} (n={res['n_compounds']})")

    print("\nDone.")


if __name__ == "__main__":
    main()
