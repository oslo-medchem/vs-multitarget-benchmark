#!/usr/bin/env python3
"""Four-way decomposition of the skill-file benchmark: preparation and parameters separated.

Reads scores_<target>_<cond>.csv for conditions A, C, D, E (A and C from the published archive,
D and E from run_cross_conditions.py) and produces:

  decomposition_4way.csv    per-target paired DeLong for every contrast, tie-corrected
  decomposition_pooled.csv  DerSimonian-Laird pooled estimate per contrast
  figure_decomposition.png  forest panels, one per contrast

AUC is the midrank (tie-corrected) estimator; the paired test is fast DeLong
(Sun and Xu, IEEE Signal Process Lett 21:1389, 2014), which supplies the covariance
between two correlated AUCs on the same cases.

Usage
  python analyse_cross_conditions.py --indir cross_conditions_out --outdir results
"""
import argparse, csv, itertools, json, os
from pathlib import Path
import numpy as np
from scipy import stats

TARGETS = ["egfr","cdk2","cox2","esr1","dpp4","ache","bace1","hsp90","p38","thrombin"]
LABELS = {"egfr":"EGFR","cdk2":"CDK2","cox2":"COX-2","esr1":"ESR1","dpp4":"DPP-4","ache":"AChE",
          "bace1":"BACE-1","hsp90":"HSP90\u03b1","p38":"p38 MAPK","thrombin":"Thrombin"}

# Contrast: (name, condition_minus, condition_base, what it isolates)
CONTRASTS = [
    ("prep_at_naive_params", "D", "A", "preparation alone, naive parameters"),
    ("prep_at_skill_params", "C", "E", "preparation alone, skill parameters"),
    ("params_at_naive_prep", "E", "A", "parameters alone, naive preparation"),
    ("params_at_skill_prep", "C", "D", "parameters alone, skill preparation"),
    ("joint_published",      "C", "A", "preparation and parameters jointly (published)"),
]


def midrank(x):
    """Midranks, used by the tie-corrected AUC and by fast DeLong."""
    x = np.asarray(x, dtype=float)
    order = np.argsort(x, kind="mergesort")
    xs = x[order]
    n = len(x)
    t = np.zeros(n, dtype=float)
    i = 0
    while i < n:
        j = i
        while j < n - 1 and xs[j + 1] == xs[i]:
            j += 1
        t[i:j + 1] = 0.5 * (i + j) + 1
        i = j + 1
    out = np.empty(n, dtype=float)
    out[order] = t
    return out


def fast_delong(scores, labels):
    """scores: (k, n) array of k predictors on the same n cases. labels: 0/1.

    Returns (aucs, covariance). Higher score must mean more likely positive.
    """
    scores = np.atleast_2d(np.asarray(scores, dtype=float))
    labels = np.asarray(labels)
    pos = labels == 1
    m, n = int(pos.sum()), int((~pos).sum())
    if m == 0 or n == 0:
        raise ValueError("need both classes")
    k = scores.shape[0]
    sp, sn = scores[:, pos], scores[:, ~pos]
    tz = np.empty((k, m + n)); tx = np.empty((k, m)); ty = np.empty((k, n))
    aucs = np.empty(k)
    for r in range(k):
        tx[r] = midrank(sp[r]); ty[r] = midrank(sn[r])
        tz[r] = midrank(np.concatenate([sp[r], sn[r]]))
        aucs[r] = (tz[r, :m].sum() / (m * n) - (m + 1) / (2.0 * n))
    v01 = (tz[:, :m] - tx[:]) / n
    v10 = 1.0 - (tz[:, m:] - ty[:]) / m
    sx = np.cov(v01); sy = np.cov(v10)
    sx = np.atleast_2d(sx); sy = np.atleast_2d(sy)
    cov = sx / m + sy / n
    return aucs, cov


def delong_paired(y, s1, s2):
    """Paired DeLong for two correlated ROC curves. Returns delta, se, ci, p, aucs."""
    aucs, cov = fast_delong(np.vstack([s1, s2]), y)
    delta = aucs[1] - aucs[0]
    var = cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
    se = float(np.sqrt(max(var, 0.0)))
    if se == 0:
        return delta, 0.0, (delta, delta), 1.0, aucs
    z = delta / se
    p = 2 * stats.norm.sf(abs(z))
    return float(delta), se, (delta - 1.96 * se, delta + 1.96 * se), float(p), aucs


def dersimonian_laird(d, se):
    d = np.asarray(d, float); se = np.asarray(se, float)
    w = 1.0 / se ** 2
    fe = (w * d).sum() / w.sum()
    Q = float((w * (d - fe) ** 2).sum())
    k = len(d)
    c = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (Q - (k - 1)) / c) if c > 0 else 0.0
    ws = 1.0 / (se ** 2 + tau2)
    m = float((ws * d).sum() / ws.sum())
    s = float(ws.sum() ** -0.5)
    I2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
    return {"estimate": m, "se": s, "ci_lo": m - 1.96 * s, "ci_hi": m + 1.96 * s,
            "p": float(2 * stats.norm.sf(abs(m / s))), "Q": Q, "Q_p": float(stats.chi2.sf(Q, k - 1)),
            "I2": I2, "tau2": tau2, "k": k}


def benjamini_hochberg(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p)
    q = p[o] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1.0)
    return out


def load(indir, target, cond):
    p = Path(indir) / f"scores_{target}_{cond}.csv"
    if not p.is_file():
        return None
    d = {}
    with open(p, newline="") as f:
        for row in csv.DictReader(f):
            if row.get("status") != "ok" or row.get("score") in ("", None):
                continue
            d[row["mol_id"]] = (float(row["score"]), int(row["label"]))
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="cross_conditions_out")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--population", choices=["all4", "pairwise"], default="all4",
                    help="all4 restricts every contrast to compounds scored in all four conditions")
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

    rows, per_contrast = [], {c[0]: [] for c in CONTRASTS}
    for t in TARGETS:
        data = {c: load(a.indir, t, c) for c in "ACDE"}
        missing = [c for c, v in data.items() if v is None]
        if missing:
            print(f"{t}: missing conditions {missing}, skipped")
            continue
        if a.population == "all4":
            common = set.intersection(*(set(v) for v in data.values()))
        for name, cm, cb, what in CONTRASTS:
            ids = sorted(common) if a.population == "all4" else sorted(set(data[cm]) & set(data[cb]))
            if len(ids) < 20:
                print(f"{t} {name}: only {len(ids)} shared compounds, skipped")
                continue
            y = np.array([data[cb][i][1] for i in ids])
            # more negative docking score ranks higher, so negate for "higher is more active"
            s_b = np.array([-data[cb][i][0] for i in ids])
            s_m = np.array([-data[cm][i][0] for i in ids])
            delta, se, ci, p, aucs = delong_paired(y, s_b, s_m)
            rec = {"target": t, "contrast": name, "isolates": what,
                   "minus": cm, "base": cb, "n": len(ids),
                   "n_actives": int((y == 1).sum()), "n_decoys": int((y == 0).sum()),
                   "auc_base": round(float(aucs[0]), 4), "auc_minus": round(float(aucs[1]), 4),
                   "delta_auc": round(delta, 4), "se": round(se, 4),
                   "ci_lo": round(ci[0], 4), "ci_hi": round(ci[1], 4), "p": p}
            rows.append(rec); per_contrast[name].append(rec)

    for name, recs in per_contrast.items():
        if recs:
            q = benjamini_hochberg([r["p"] for r in recs])
            for r, qq in zip(recs, q):
                r["p_adj_BH"] = float(qq)

    if not rows:
        print("no contrasts computed; check --indir")
        return
    with open(out / "decomposition_4way.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    pooled = []
    for name, cm, cb, what in CONTRASTS:
        recs = per_contrast[name]
        if len(recs) < 2:
            continue
        r = dersimonian_laird([x["delta_auc"] for x in recs], [x["se"] for x in recs])
        r.update({"contrast": name, "isolates": what, "minus": cm, "base": cb,
                  "n_targets": len(recs),
                  "n_sig_BH": int(sum(1 for x in recs if x.get("p_adj_BH", 1) < 0.05)),
                  "n_sig_pos": int(sum(1 for x in recs if x.get("p_adj_BH", 1) < 0.05 and x["delta_auc"] > 0)),
                  "n_sig_neg": int(sum(1 for x in recs if x.get("p_adj_BH", 1) < 0.05 and x["delta_auc"] < 0))})
        pooled.append(r)
    with open(out / "decomposition_pooled.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(pooled[0].keys())); w.writeheader(); w.writerows(pooled)

    # interaction on the pooled scale, reported with a note that it is a derived quantity
    byname = {r["contrast"]: r for r in pooled}
    if {"prep_at_skill_params", "prep_at_naive_params"} <= byname.keys():
        a1, a2 = byname["prep_at_skill_params"], byname["prep_at_naive_params"]
        inter = a1["estimate"] - a2["estimate"]
        se = (a1["se"] ** 2 + a2["se"] ** 2) ** 0.5
        json.dump({"interaction_estimate": inter, "se": se,
                   "ci": [inter - 1.96 * se, inter + 1.96 * se],
                   "p": float(2 * stats.norm.sf(abs(inter / se))),
                   "note": "difference of pooled preparation effects at the two parameter settings; "
                           "treats the two pooled estimates as independent, which is conservative only "
                           "if the per-target correlation is positive"},
                  open(out / "interaction.json", "w"), indent=1)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        show = [c for c in CONTRASTS if per_contrast[c[0]]][:4]
        fig, axes = plt.subplots(1, len(show), figsize=(4.2 * len(show), 5.2), sharey=True)
        axes = np.atleast_1d(axes)
        for ax, (name, cm, cb, what) in zip(axes, show):
            recs = per_contrast[name]
            ys = np.arange(len(recs))[::-1]
            for yy, r in zip(ys, recs):
                sig = r.get("p_adj_BH", 1) < 0.05
                col = "#b4442a" if sig else "#6e6e6e"
                ax.plot([r["ci_lo"], r["ci_hi"]], [yy, yy], color=col, lw=1.4)
                ax.plot([r["delta_auc"]], [yy], "s", color=col, ms=5)
            pr = byname.get(name)
            if pr:
                ax.plot([pr["ci_lo"], pr["ci_hi"]], [-1, -1], color="#1f4e79", lw=2.4)
                ax.plot([pr["estimate"]], [-1], "D", color="#1f4e79", ms=8)
            ax.axvline(0, color="0.6", lw=0.9)
            ax.set_yticks(list(ys) + [-1])
            ax.set_yticklabels([LABELS.get(r["target"], r["target"]) for r in recs] + ["Pooled"])
            ax.set_title(what.replace(", ", "\n"), fontsize=9)
            ax.set_xlabel("$\\Delta$AUC")
        fig.tight_layout()
        fig.savefig(out / "figure_decomposition.png", dpi=300)
    except Exception as e:  # plotting is not load-bearing
        print("figure skipped:", e)

    print(f"wrote {out/'decomposition_4way.csv'} ({len(rows)} rows) and "
          f"{out/'decomposition_pooled.csv'} ({len(pooled)} contrasts)")
    for r in pooled:
        print(f"  {r['contrast']:22s} pooled {r['estimate']:+.4f} "
              f"[{r['ci_lo']:+.3f},{r['ci_hi']:+.3f}] p={r['p']:.3f} I2={r['I2']:.0f}% "
              f"sig {r['n_sig_BH']}/{r['n_targets']}")


if __name__ == "__main__":
    main()
