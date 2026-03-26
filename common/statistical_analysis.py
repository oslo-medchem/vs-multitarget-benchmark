#!/usr/bin/env python3
"""
Statistical analysis comparing naive vs skill-guided VS protocols.

Usage:
    python statistical_analysis.py --target-dir targets/egfr/

Tests: DeLong's AUC comparison, bootstrap permutation, Mann-Whitney U,
       Benjamini-Hochberg FDR, Net Reclassification Improvement.

Inputs:  target_dir/05_evaluation/results/{y_true,y_score}_{naive,skill}.npy,
         metrics_{naive,skill}.json
Outputs: target_dir/05_evaluation/results/comparison_stats.json,
         target_dir/05_evaluation/results/bootstrap_stats.json
"""
import argparse
import json
import sys
import numpy as np
from pathlib import Path
from scipy import stats
from scipy.stats import mannwhitneyu, bootstrap as scipy_bootstrap
from sklearn.metrics import roc_auc_score

N_BOOTSTRAP = 10_000
N_PERM = 10_000
RNG = np.random.default_rng(42)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_data(results_dir, protocol):
    yt = np.load(results_dir / f"y_true_{protocol}.npy")
    ys = np.load(results_dir / f"y_score_{protocol}.npy")
    return yt, ys


def load_data_intersection(results_dir):
    """Load scores for compounds docked in BOTH protocols (intersection).

    This removes the confound where one protocol has more missing compounds
    (e.g., due to PAINS filtering) getting worst-rank scores.
    """
    cids_n = set(np.load(results_dir / "cids_naive_docked.npy", allow_pickle=True))
    cids_s = set(np.load(results_dir / "cids_skill_docked.npy", allow_pickle=True))
    common = sorted(cids_n & cids_s)

    if not common:
        return None, None, None, None, None

    # Reload full-library data and filter to intersection
    import csv
    labels_csv = results_dir.parent.parent / "03_library_preparation" / "library_labels.csv"
    labels = {}
    for r in csv.DictReader(open(labels_csv)):
        labels[r["compound_id"]] = 1 if r["label"] == "active" else 0

    scores_n = {}
    for r in csv.DictReader(open(results_dir.parent.parent / "04_docking" / "naive" / "scores_naive.csv")):
        scores_n[r["compound_id"]] = float(r["score"])
    scores_s = {}
    for r in csv.DictReader(open(results_dir.parent.parent / "04_docking" / "skill" / "scores_skill.csv")):
        scores_s[r["compound_id"]] = float(r["score"])

    y_true = np.array([labels[c] for c in common])
    ys_n = -np.array([scores_n[c] for c in common])
    ys_s = -np.array([scores_s[c] for c in common])

    return common, y_true, ys_n, ys_s, len(common)


def load_metrics(results_dir, protocol):
    with open(results_dir / f"metrics_{protocol}.json") as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# DeLong's test (compare two correlated ROC curves on the same sample)
# Reference: DeLong et al. Biometrika 1988
# ---------------------------------------------------------------------------

def _delong_placement(y_true, y_score):
    """Compute placement values for DeLong's test."""
    n_pos = int(y_true.sum())
    n_neg = int((1 - y_true).sum())
    pos_scores = y_score[y_true == 1]
    neg_scores = y_score[y_true == 0]
    # V10: for each positive, fraction of negatives it outranks
    v10 = np.array([np.mean(ps > neg_scores) + 0.5 * np.mean(ps == neg_scores) for ps in pos_scores])
    # V01: for each negative, fraction of positives it is outranked by
    v01 = np.array([np.mean(ns < pos_scores) + 0.5 * np.mean(ns == pos_scores) for ns in neg_scores])
    return v10, v01, n_pos, n_neg


def delong_test(y_true, y_score_a, y_score_b):
    """
    DeLong's test for comparing two correlated ROC AUC estimates.
    Returns: z_stat, p_value, delta_auc, ci_low, ci_high (95% CI on delta AUC)
    """
    v10_a, v01_a, n_pos, n_neg = _delong_placement(y_true, y_score_a)
    v10_b, v01_b, _, _ = _delong_placement(y_true, y_score_b)

    auc_a = v10_a.mean()
    auc_b = v10_b.mean()
    delta = auc_b - auc_a

    # Covariance matrix of the two AUC estimators
    s_10_aa = np.var(v10_a, ddof=1) / n_pos
    s_10_bb = np.var(v10_b, ddof=1) / n_pos
    s_10_ab = np.cov(v10_a, v10_b, ddof=1)[0, 1] / n_pos

    s_01_aa = np.var(v01_a, ddof=1) / n_neg
    s_01_bb = np.var(v01_b, ddof=1) / n_neg
    s_01_ab = np.cov(v01_a, v01_b, ddof=1)[0, 1] / n_neg

    var_delta = s_10_aa + s_10_bb - 2 * s_10_ab + s_01_aa + s_01_bb - 2 * s_01_ab
    if var_delta <= 0:
        return 0.0, 1.0, float(delta), None, None

    se = np.sqrt(var_delta)
    z = delta / se
    p = 2 * (1 - stats.norm.cdf(abs(z)))

    ci_low = delta - 1.96 * se
    ci_high = delta + 1.96 * se

    return float(z), float(p), float(delta), float(ci_low), float(ci_high)


# ---------------------------------------------------------------------------
# Bootstrap permutation test for a generic metric difference
# ---------------------------------------------------------------------------

def _bedroc_sklearn(y_true, y_score_rank, alpha=20.0):
    n = len(y_true)
    n_actives = y_true.sum()
    if n_actives == 0 or n_actives == n:
        return 0.5
    ra = n_actives / n
    order = np.argsort(-y_score_rank)
    y_sorted = y_true[order]
    ranks = np.where(y_sorted == 1)[0] + 1
    s = np.sum(np.exp(-alpha * ranks / n))
    d_val = np.exp(alpha / n) - 1
    ri = np.exp(-alpha * ra) - np.exp(-alpha)
    bedroc = (s * d_val * np.sinh(alpha / 2) /
              (np.cosh(alpha / 2) - np.cosh(alpha / 2 - alpha * ra))) * ra + ri / (1 - np.exp(alpha * (1 - ra)))
    min_b = (1 - np.exp(alpha * ra)) / (1 - np.exp(alpha)) * ra
    max_b = (1 - np.exp(-alpha * ra)) / (1 - np.exp(-alpha)) * ra
    return float((bedroc - min_b) / (max_b - min_b)) if (max_b - min_b) > 0 else 0.5


def _calc_ef(y_true, y_score_rank, fraction):
    n = len(y_true)
    n_actives = y_true.sum()
    n_top = max(1, int(n * fraction))
    order = np.argsort(-y_score_rank)
    hits = y_true[order[:n_top]].sum()
    expected = n_actives * fraction
    return float(hits / expected) if expected > 0 else 0.0


METRIC_FNS = {
    "roc_auc": lambda yt, ys: float(roc_auc_score(yt, ys)) if 0 < yt.sum() < len(yt) else 0.5,
    "bedroc_20": lambda yt, ys: _bedroc_sklearn(yt, ys, 20.0),
    "ef_1pct": lambda yt, ys: _calc_ef(yt, ys, 0.01),
    "ef_5pct": lambda yt, ys: _calc_ef(yt, ys, 0.05),
}


def permutation_test_delta(y_true_a, y_score_a, y_true_b, y_score_b, metric_fn, n_perm=N_PERM):
    """
    Two-sided permutation test for H0: metric(A) == metric(B).
    Permutes the protocol labels (not compound labels).
    Returns: observed_delta, p_value
    """
    obs_a = metric_fn(y_true_a, y_score_a)
    obs_b = metric_fn(y_true_b, y_score_b)
    observed_delta = obs_b - obs_a

    # Permutation: randomly swap which score vector goes to "a" vs "b"
    all_scores_a = np.concatenate([y_score_a, y_score_b])
    n = len(y_score_a)
    null_deltas = []
    for _ in range(n_perm):
        perm = RNG.permutation(len(all_scores_a))
        s_a = all_scores_a[perm[:n]]
        s_b = all_scores_a[perm[n:]]
        # y_true stays the same for each (same compounds)
        d = metric_fn(y_true_b, s_b) - metric_fn(y_true_a, s_a)
        null_deltas.append(d)

    null_deltas = np.array(null_deltas)
    p_val = float(np.mean(np.abs(null_deltas) >= np.abs(observed_delta)))

    return float(observed_delta), p_val


# ---------------------------------------------------------------------------
# Mann-Whitney U test (actives vs decoys score distributions)
# ---------------------------------------------------------------------------

def mannwhitney_test(y_true, y_score_raw):
    """
    Compare raw docking scores of actives vs decoys.
    Lower (more negative) score = better binding.
    H1: actives have lower (better) docking scores than decoys.
    Returns: U, p_value, r (rank-biserial correlation effect size)
    """
    actives_scores = y_score_raw[y_true == 1]
    decoys_scores = y_score_raw[y_true == 0]

    # Exclude failed dockings (score == 0.0)
    actives_scores = actives_scores[actives_scores != 0.0]
    decoys_scores = decoys_scores[decoys_scores != 0.0]

    if len(actives_scores) == 0 or len(decoys_scores) == 0:
        return None, None, None

    # alternative='less': test that actives scores < decoys scores
    result = mannwhitneyu(actives_scores, decoys_scores, alternative='less')
    U = float(result.statistic)
    p = float(result.pvalue)

    # Rank-biserial correlation as effect size: r = 1 - 2U/(n1*n2)
    n1, n2 = len(actives_scores), len(decoys_scores)
    r = 1 - 2 * U / (n1 * n2)

    return U, p, float(r)


# ---------------------------------------------------------------------------
# Benjamini-Hochberg FDR correction
# ---------------------------------------------------------------------------

def bh_correction(p_values):
    """Benjamini-Hochberg FDR correction. Returns corrected p-values."""
    n = len(p_values)
    sorted_idx = np.argsort(p_values)
    sorted_p = np.array(p_values)[sorted_idx]
    corrected = np.zeros(n)
    for i in range(n - 1, -1, -1):
        if i == n - 1:
            corrected[i] = sorted_p[i]
        else:
            corrected[i] = min(corrected[i + 1], sorted_p[i] * n / (i + 1))
    corrected = np.minimum(corrected, 1.0)
    # Map back to original order
    result = np.zeros(n)
    result[sorted_idx] = corrected
    return result.tolist()


# ---------------------------------------------------------------------------
# Net Reclassification Improvement at EF1% threshold
# ---------------------------------------------------------------------------

def calc_nri_at_ef_threshold(y_true_a, y_score_a, y_true_b, y_score_b, fraction=0.01):
    """
    NRI at EF1% threshold.
    Compounds in top 1% of skill that were NOT in top 1% of naive (up-events for actives).
    Returns: NRI, NRI_actives, NRI_decoys, 95% CI via bootstrap.
    """
    n = len(y_true_a)
    n_top = max(1, int(n * fraction))

    order_a = np.argsort(-y_score_a)
    order_b = np.argsort(-y_score_b)

    in_top_a = np.zeros(n, dtype=bool)
    in_top_b = np.zeros(n, dtype=bool)
    in_top_a[order_a[:n_top]] = True
    in_top_b[order_b[:n_top]] = True

    # Reclassification among actives: moved up - moved down
    act = y_true_a == 1
    up_act = np.sum(in_top_b[act] & ~in_top_a[act])
    down_act = np.sum(~in_top_b[act] & in_top_a[act])
    n_act = act.sum()

    # Reclassification among decoys (up = wrong, down = correct)
    dec = y_true_a == 0
    up_dec = np.sum(in_top_b[dec] & ~in_top_a[dec])
    down_dec = np.sum(~in_top_b[dec] & in_top_a[dec])
    n_dec = dec.sum()

    nri_act = float((up_act - down_act) / n_act) if n_act > 0 else 0.0
    nri_dec = float((down_dec - up_dec) / n_dec) if n_dec > 0 else 0.0
    nri = nri_act + nri_dec

    return {
        "nri": float(nri),
        "nri_actives": float(nri_act),
        "nri_decoys": float(nri_dec),
        "n_actives_up": int(up_act),
        "n_actives_down": int(down_act),
    }


# ---------------------------------------------------------------------------
# Bootstrap CI for delta metric
# ---------------------------------------------------------------------------

def bootstrap_delta_ci(y_true_a, y_score_a, y_true_b, y_score_b, metric_fn,
                        n_resamples=N_BOOTSTRAP):
    """Bootstrap 95% CI for (metric_b - metric_a)."""
    n = len(y_true_a)
    deltas = []
    for _ in range(n_resamples):
        idx = RNG.integers(0, n, n)
        yt_a = y_true_a[idx]
        ys_a = y_score_a[idx]
        yt_b = y_true_b[idx]
        ys_b = y_score_b[idx]
        if yt_a.sum() == 0 or yt_b.sum() == 0:
            continue
        try:
            d = metric_fn(yt_b, ys_b) - metric_fn(yt_a, ys_a)
            deltas.append(d)
        except Exception:
            pass

    if len(deltas) < 100:
        return {"ci_low": None, "ci_high": None, "n_valid": len(deltas)}

    deltas = np.array(deltas)
    return {
        "ci_low": float(np.percentile(deltas, 2.5)),
        "ci_high": float(np.percentile(deltas, 97.5)),
        "n_valid": len(deltas),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Statistical comparison of naive vs skill-guided VS protocols."
    )
    parser.add_argument("--target-dir", required=True,
                        help="Path to target directory (e.g. targets/egfr/)")
    args = parser.parse_args()

    target_dir = Path(args.target_dir).resolve()
    results_dir = target_dir / "05_evaluation" / "results"

    if not results_dir.exists():
        print(f"ERROR: {results_dir} not found. Run compute_metrics.py first.")
        sys.exit(1)

    print("Loading data...")
    yt_n, ys_n = load_data(results_dir, "naive")
    yt_s, ys_s = load_data(results_dir, "skill")
    m_n = load_metrics(results_dir, "naive")
    m_s = load_metrics(results_dir, "skill")

    # Try loading intersection data (compounds docked in both protocols)
    try:
        common_cids, yt_int, ys_n_int, ys_s_int, n_int = load_data_intersection(results_dir)
        has_intersection = common_cids is not None and n_int > 0
        if has_intersection:
            n_act_int = int(yt_int.sum())
            print(f"Intersection: {n_int} compounds ({n_act_int} actives) docked in BOTH protocols")
    except Exception as e:
        has_intersection = False
        print(f"Intersection data not available ({e}), using full-library only")

    # Raw scores (before negation) -- needed for Mann-Whitney
    # y_score saved in compute_metrics is already negated (y_score_rank = -y_score_raw)
    # So raw score = -y_score_rank (more negative raw = better)
    ys_n_raw = -ys_n
    ys_s_raw = -ys_s

    comparison = {}
    bootstrap_stats = {}

    # -----------------------------------------------------------------------
    # 1. DeLong's test -- AUC comparison
    # -----------------------------------------------------------------------
    print("\n=== DeLong's Test (AUC naive vs skill) ===")
    z, p_delong, delta_auc, ci_low, ci_high = delong_test(yt_n, ys_n, ys_s)
    print(f"  AUC naive:  {m_n['roc_auc']:.4f}")
    print(f"  AUC skill:  {m_s['roc_auc']:.4f}")
    print(f"  ΔAUC:       {delta_auc:.4f} [{ci_low:.4f}, {ci_high:.4f}]")
    print(f"  Z={z:.3f}, p={p_delong:.4f}")
    comparison["delong_auc_full"] = {
        "analysis": "full_library_with_score0_penalty",
        "auc_naive": m_n.get("full_library", m_n).get("roc_auc", m_n["roc_auc"]),
        "auc_skill": m_s.get("full_library", m_s).get("roc_auc", m_s["roc_auc"]),
        "delta_auc": delta_auc,
        "ci_95_low": ci_low,
        "ci_95_high": ci_high,
        "z_stat": z,
        "p_value": p_delong,
    }

    # -----------------------------------------------------------------------
    # 1b. DeLong's test on INTERSECTION (compounds docked in both protocols)
    # -----------------------------------------------------------------------
    if has_intersection:
        print("\n=== DeLong's Test (INTERSECTION — docked in both protocols) ===")
        z_int, p_int, delta_int, ci_lo_int, ci_hi_int = delong_test(yt_int, ys_n_int, ys_s_int)
        auc_n_int = float(roc_auc_score(yt_int, ys_n_int)) if yt_int.sum() > 0 else 0.5
        auc_s_int = float(roc_auc_score(yt_int, ys_s_int)) if yt_int.sum() > 0 else 0.5
        print(f"  N compounds: {n_int} ({n_act_int} actives)")
        print(f"  AUC naive:  {auc_n_int:.4f}")
        print(f"  AUC skill:  {auc_s_int:.4f}")
        print(f"  ΔAUC:       {delta_int:.4f} [{ci_lo_int:.4f}, {ci_hi_int:.4f}]")
        print(f"  Z={z_int:.3f}, p={p_int:.4f}")
        comparison["delong_auc"] = {
            "analysis": "intersection_docked_in_both",
            "n_compounds": n_int,
            "n_actives": n_act_int,
            "auc_naive": auc_n_int,
            "auc_skill": auc_s_int,
            "delta_auc": delta_int,
            "ci_95_low": ci_lo_int,
            "ci_95_high": ci_hi_int,
            "z_stat": z_int,
            "p_value": p_int,
        }
    else:
        # Fall back to full-library comparison
        comparison["delong_auc"] = comparison["delong_auc_full"]

    # -----------------------------------------------------------------------
    # 2. Permutation tests for BEDROC and EF differences
    # -----------------------------------------------------------------------
    print("\n=== Permutation Tests (10,000 permutations) ===")
    perm_p_values = {}
    for mname, fn in METRIC_FNS.items():
        print(f"  {mname}...", end=" ", flush=True)
        obs_delta, p_perm = permutation_test_delta(yt_n, ys_n, yt_s, ys_s, fn)
        print(f"Δ={obs_delta:.4f}, p={p_perm:.4f}")
        perm_p_values[mname] = p_perm
        comparison[f"permutation_{mname}"] = {
            "observed_delta": obs_delta,
            "p_value_permutation": p_perm,
            "metric_naive": m_n.get(mname),
            "metric_skill": m_s.get(mname),
        }

    # -----------------------------------------------------------------------
    # 3. Mann-Whitney U test (actives vs decoys within each protocol)
    # -----------------------------------------------------------------------
    print("\n=== Mann-Whitney U (actives vs decoys score distributions) ===")
    for protocol, yt, ys_raw in [("naive", yt_n, ys_n_raw), ("skill", yt_s, ys_s_raw)]:
        U, p_mw, r = mannwhitney_test(yt, ys_raw)
        if U is not None:
            print(f"  {protocol}: U={U:.0f}, p={p_mw:.2e}, r={r:.4f}")
        else:
            print(f"  {protocol}: insufficient data")
        comparison[f"mannwhitney_{protocol}"] = {
            "U_statistic": U,
            "p_value": p_mw,
            "rank_biserial_r": r,
            "interpretation": "actives ranked better than decoys (p < 0.05 = significant)" if p_mw and p_mw < 0.05 else "not significant",
        }

    # -----------------------------------------------------------------------
    # 4. Bootstrap CIs for metric deltas
    # -----------------------------------------------------------------------
    print("\n=== Bootstrap CIs for Metric Deltas ===")
    for mname, fn in METRIC_FNS.items():
        print(f"  {mname}...", end=" ", flush=True)
        ci = bootstrap_delta_ci(yt_n, ys_n, yt_s, ys_s, fn, N_BOOTSTRAP)
        print(f"Δ CI=[{ci['ci_low']:.4f}, {ci['ci_high']:.4f}]" if ci["ci_low"] is not None else "failed")
        bootstrap_stats[mname] = ci

    # -----------------------------------------------------------------------
    # 5. Benjamini-Hochberg FDR correction
    # -----------------------------------------------------------------------
    print("\n=== Benjamini-Hochberg FDR Correction ===")
    all_p_raw = {}
    all_p_raw["delong_auc"] = comparison["delong_auc"]["p_value"]
    for mname in METRIC_FNS:
        all_p_raw[f"perm_{mname}"] = comparison[f"permutation_{mname}"]["p_value_permutation"]
    all_p_raw["mw_naive"] = comparison["mannwhitney_naive"]["p_value"] or 1.0
    all_p_raw["mw_skill"] = comparison["mannwhitney_skill"]["p_value"] or 1.0

    keys = list(all_p_raw.keys())
    raw_vals = [all_p_raw[k] for k in keys]
    corrected_vals = bh_correction(raw_vals)

    bh_results = {}
    for k, p_raw, p_corr in zip(keys, raw_vals, corrected_vals):
        bh_results[k] = {"p_raw": p_raw, "p_bh": p_corr, "significant_bh": p_corr < 0.05}
        sig = "* " if p_corr < 0.05 else "  "
        print(f"  {sig}{k}: p_raw={p_raw:.4f} -> p_BH={p_corr:.4f}")

    comparison["bh_correction"] = bh_results

    # -----------------------------------------------------------------------
    # 6. Net Reclassification Improvement at EF1%
    # -----------------------------------------------------------------------
    print("\n=== NRI at EF1% threshold ===")
    nri_result = calc_nri_at_ef_threshold(yt_n, ys_n, yt_s, ys_s, fraction=0.01)
    print(f"  NRI total: {nri_result['nri']:.4f}")
    print(f"  NRI actives: {nri_result['nri_actives']:.4f}  (actives up: {nri_result['n_actives_up']}, down: {nri_result['n_actives_down']})")
    print(f"  NRI decoys: {nri_result['nri_decoys']:.4f}")
    comparison["nri_ef1pct"] = nri_result

    # -----------------------------------------------------------------------
    # 7. PAINS/Brenk filtering impact analysis
    # -----------------------------------------------------------------------
    print("\n=== PAINS/Brenk Filtering Impact ===")
    filter_report = target_dir / "03_library_preparation" / "skill" / "pdbqt" / "filter_report.csv"
    pains_analysis = {}
    if filter_report.exists():
        import csv as csv_mod
        filtered_cids = set()
        for row in csv_mod.DictReader(open(filter_report)):
            filtered_cids.add(row["compound_id"])

        # Load labels to check which filtered compounds are actives
        labels_csv = target_dir / "03_library_preparation" / "library_labels.csv"
        labels = {r["compound_id"]: r["label"] for r in csv_mod.DictReader(open(labels_csv))}

        n_filtered_total = len(filtered_cids)
        filtered_actives = [c for c in filtered_cids if labels.get(c) == "active"]
        filtered_decoys = [c for c in filtered_cids if labels.get(c) == "decoy"]
        n_actives_total = sum(1 for v in labels.values() if v == "active")
        n_decoys_total = sum(1 for v in labels.values() if v == "decoy")

        print(f"  Filtered by PAINS/Brenk: {n_filtered_total} compounds")
        print(f"    Actives removed: {len(filtered_actives)}/{n_actives_total} ({100*len(filtered_actives)/max(1,n_actives_total):.1f}%)")
        print(f"    Decoys removed:  {len(filtered_decoys)}/{n_decoys_total} ({100*len(filtered_decoys)/max(1,n_decoys_total):.1f}%)")

        # Check if filtered actives would have been hits in naive
        scores_naive_csv = target_dir / "04_docking" / "naive" / "scores_naive.csv"
        if scores_naive_csv.exists():
            naive_scores = {}
            for r in csv_mod.DictReader(open(scores_naive_csv)):
                try:
                    naive_scores[r["compound_id"]] = float(r["score"])
                except (ValueError, KeyError):
                    pass

            # Get naive score distribution for actives
            all_naive_active_scores = [naive_scores.get(c, 0.0) for c in labels if labels[c] == "active" and naive_scores.get(c, 0.0) != 0.0]
            if all_naive_active_scores:
                median_active = np.median(all_naive_active_scores)
                filtered_active_scores = [(c, naive_scores.get(c, 0.0)) for c in filtered_actives if naive_scores.get(c, 0.0) != 0.0]
                n_would_be_hits = sum(1 for _, s in filtered_active_scores if s <= median_active)
                print(f"    Filtered actives docked in naive: {len(filtered_active_scores)}/{len(filtered_actives)}")
                print(f"    Would-be hits (score <= median active): {n_would_be_hits}")

        pains_analysis = {
            "n_filtered_total": n_filtered_total,
            "n_filtered_actives": len(filtered_actives),
            "n_filtered_decoys": len(filtered_decoys),
            "filtered_active_ids": filtered_actives,
            "pct_actives_removed": round(100 * len(filtered_actives) / max(1, n_actives_total), 1),
            "pct_decoys_removed": round(100 * len(filtered_decoys) / max(1, n_decoys_total), 1),
        }
    else:
        print("  No filter_report.csv found")
    comparison["pains_filtering"] = pains_analysis

    # -----------------------------------------------------------------------
    # Save outputs
    # -----------------------------------------------------------------------
    with open(results_dir / "comparison_stats.json", "w") as f:
        json.dump(comparison, f, indent=2, default=lambda x: None if np.isnan(x) else float(x))

    with open(results_dir / "bootstrap_stats.json", "w") as f:
        json.dump(bootstrap_stats, f, indent=2)

    print(f"\nSaved: {results_dir}/comparison_stats.json")
    print(f"Saved: {results_dir}/bootstrap_stats.json")
    print("\nDone.")


if __name__ == "__main__":
    main()
