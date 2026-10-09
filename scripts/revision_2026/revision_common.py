#!/usr/bin/env python3
"""Shared loaders and metrics for the September 2026 MLST rewrite analyses.

Metric definitions are taken verbatim from the benchmark's own
common/compute_metrics.py so that every number recomputed here is on the same
footing as the originally published values. The only additions are the paired
DeLong implementation (also copied from the original ablation script) and the
population bookkeeping the referees asked for.

Data root must point at the extracted vs_benchmark_multitarget tree.
"""
from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.stats import norm
from sklearn.metrics import roc_auc_score, roc_curve, average_precision_score

TARGETS = ["fpr2", "egfr", "cdk2", "cox2", "esr1", "dpp4", "ache",
           "bace1", "hsp90", "p38", "thrombin"]

PRETTY = {"fpr2": "FPR2", "egfr": "EGFR", "cdk2": "CDK2", "cox2": "COX-2",
          "esr1": "ESR1", "dpp4": "DPP-4", "ache": "AChE", "bace1": "BACE-1",
          "hsp90": "HSP90", "p38": "p38", "thrombin": "Thrombin"}

TARGET_CLASS = {"fpr2": "GPCR", "egfr": "kinase", "cdk2": "kinase", "p38": "kinase",
                "cox2": "oxidoreductase", "esr1": "nuclear receptor",
                "dpp4": "serine protease", "thrombin": "serine protease",
                "bace1": "aspartyl protease", "ache": "hydrolase",
                "hsp90": "chaperone"}


# ----------------------------------------------------------------- loaders
@dataclass
class TargetData:
    target: str
    labels: dict = field(default_factory=dict)      # cid -> 1 active / 0 decoy
    smiles: dict = field(default_factory=dict)      # cid -> smiles
    source_id: dict = field(default_factory=dict)   # cid -> ChEMBL id
    naive: dict = field(default_factory=dict)       # cid -> score (<0 only)
    skill: dict = field(default_factory=dict)
    nofilter: dict = field(default_factory=dict)    # ablation re-dock scores
    filtered: dict = field(default_factory=dict)    # cid -> filter_hit
    prep_skill: dict = field(default_factory=dict)  # cid -> (status, note)
    prep_naive: dict = field(default_factory=dict)  # cid -> (success, message)


def _read_scores(path: Path) -> dict:
    """Docked-only convention from the original code: score == 0.0 means no dock."""
    out = {}
    if not path.exists():
        return out
    with path.open() as f:
        for r in csv.DictReader(f):
            try:
                s = float(r["score"])
            except (ValueError, KeyError, TypeError):
                continue
            if s != 0.0:
                out[r["compound_id"]] = s
    return out


def load_target(target: str, data_root: Path, ablation_root: Path | None = None) -> TargetData:
    t = TargetData(target=target)
    tdir = data_root / "targets" / target

    with (tdir / "03_library_preparation" / "library_labels.csv").open() as f:
        for r in csv.DictReader(f):
            cid = r["compound_id"]
            t.labels[cid] = 1 if r["label"] == "active" else 0
            t.smiles[cid] = r.get("smiles", "")
            t.source_id[cid] = r.get("source_id", "")

    t.naive = _read_scores(tdir / "04_docking" / "naive" / "scores_naive.csv")
    t.skill = _read_scores(tdir / "04_docking" / "skill" / "scores_skill.csv")

    filt = tdir / "03_library_preparation" / "skill" / "pdbqt" / "filter_report.csv"
    if filt.exists():
        with filt.open() as f:
            for r in csv.DictReader(f):
                t.filtered[r["compound_id"]] = r.get("filter_hit", "")

    ps = tdir / "03_library_preparation" / "skill" / "pdbqt" / "preparation_log.csv"
    if ps.exists():
        with ps.open() as f:
            for r in csv.DictReader(f):
                t.prep_skill[r["compound_id"]] = (r.get("status", ""), r.get("note", ""))

    pn = tdir / "03_library_preparation" / "naive" / "pdbqt" / "preparation_log.csv"
    if pn.exists():
        with pn.open() as f:
            for r in csv.DictReader(f):
                t.prep_naive[r["compound_id"]] = (r.get("success", ""), r.get("message", ""))

    if ablation_root is not None:
        ab = ablation_root / target / "skill_nofilter" / "scores_skill_nofilter.csv"
        t.nofilter = _read_scores(ab)
    return t


# ----------------------------------------------------------------- metrics
# The published benchmark computed ROC AUC, BEDROC, EF and RIE through RDKit's
# ML.Scoring module (compute_metrics.py takes the RDKIT_SCORING branch when
# RDKit is importable, which it was). The same code path is used here so that
# recomputed values are directly comparable with the published ones. LogAUC and
# AUPR are computed exactly as in the original.
from rdkit.ML.Scoring.Scoring import CalcAUC, CalcBEDROC, CalcEnrichment, CalcRIE


def ranked_list(y_true, y_score_rank, tie_policy="pessimistic"):
    """Order (score, label) pairs by decreasing score for the RDKit scorers.

    RDKit's ML.Scoring functions assume the caller has already sorted the list
    and they do not correct for ties. The benchmark's original helper sorted
    with ``sorted(zip(score, label), reverse=True)``, which breaks ties by
    label in DESCENDING order, so a tied active is always placed ahead of a
    tied decoy. Where ties are rare, as in the docked-only population, the
    effect is below 1e-3. Where a large block of compounds shares one score, as
    in the end-to-end population in which every unscored compound is given the
    same worst rank, that convention inflates the AUC substantially.

    tie_policy:
      "pessimistic" (default) places actives last within a tied block, which is
          the conservative convention for compounds a protocol never scored.
      "legacy" reproduces the original optimistic ordering.
    """
    if tie_policy == "legacy":
        return [[sc, lab] for sc, lab in sorted(zip(y_score_rank, y_true), reverse=True)]
    pairs = sorted(zip(y_score_rank, y_true), key=lambda p: (-p[0], p[1]))
    return [[sc, lab] for sc, lab in pairs]


# Backwards-compatible alias used by the validation script.
def rdkit_sorted_list(y_true, y_score_rank):
    return ranked_list(y_true, y_score_rank, tie_policy="legacy")


def calc_logauc(y_true, y_score_rank, fpr_min=0.001, fpr_max=1.0):
    fpr, tpr, _ = roc_curve(y_true, y_score_rank)
    fpr_log = np.logspace(np.log10(fpr_min), np.log10(fpr_max), 1000)
    tpr_interp = np.interp(fpr_log, fpr, tpr)
    trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))
    return float(trapz(tpr_interp, np.log10(fpr_log)) / (np.log10(fpr_max) - np.log10(fpr_min)))


def all_metrics(y_true: np.ndarray, y_rank: np.ndarray, tie_policy: str = "pessimistic") -> dict:
    """Metric panel for one protocol on one population.

    ROC AUC is always the tie-corrected rank statistic (sklearn's
    roc_auc_score, which gives tied pairs half credit and is identical to the
    AUC implied by the DeLong structural components). The early-recognition
    metrics come from RDKit's ML.Scoring module, as in the original benchmark,
    on a list ordered by ``ranked_list`` under the stated tie policy.
    """
    na, nd = int(y_true.sum()), int(len(y_true) - y_true.sum())
    base = {"n_total": int(len(y_true)), "n_actives": na, "n_decoys": nd,
            "n_tied_at_worst": int(np.sum(y_rank == y_rank.min()))}
    if na < 2 or nd < 2:
        return base
    lst = ranked_list(y_true, y_rank, tie_policy)
    base.update({
        "roc_auc": float(roc_auc_score(y_true, y_rank)),
        "roc_auc_legacy_ties": float(CalcAUC(ranked_list(y_true, y_rank, "legacy"), 1)),
        "bedroc_20": float(CalcBEDROC(lst, 1, 20.0)),
        "ef_1pct": float(CalcEnrichment(lst, 1, [0.01])[0]),
        "ef_5pct": float(CalcEnrichment(lst, 1, [0.05])[0]),
        "rie_20": float(CalcRIE(lst, 1, 20.0)),
        "log_auc": calc_logauc(y_true, y_rank),
        "aupr": float(average_precision_score(y_true, y_rank)),
    })
    return base


def end_to_end_ranks(labels: dict, scores: dict, convention: str = "strict"):
    """Rank keys for the end-to-end (intention-to-screen) population.

    Every library compound is evaluated. Docking scores are negative and more
    negative is better, so the rank key is -score and larger is better.
    Compounds a protocol never scored, whether excluded deliberately by the
    substructure filter or lost to a preparation or docking failure, are placed
    at the bottom.

    A small fraction of compounds receive an unfavourable positive Vina score
    (0.2% of naive and 2.3% of skill scores in this benchmark), whose rank key
    is therefore negative. The convention decides where the unscored block sits
    relative to those:

      "strict" (default) places unscored compounds strictly below every scored
          compound, including the unfavourable ones.
      "zero" reproduces the simpler convention of giving unscored compounds a
          rank key of exactly 0.0, so an unfavourable positive score ranks
          below an unscored compound.

    Both are reported in the revision so that the choice is visible and does
    not silently drive the result.
    """
    cids = sorted(labels)
    y = np.array([labels[c] for c in cids])
    vals = [-scores[c] for c in cids if c in scores]
    floor = (min(vals) - 1.0) if (convention == "strict" and vals) else 0.0
    rank = np.array([-scores[c] if c in scores else floor for c in cids])
    return cids, y, rank


# ----------------------------------------------------------------- DeLong
def delong_paired(y: np.ndarray, scores_a: np.ndarray, scores_b: np.ndarray):
    """Paired DeLong test. Returns (auc_b - auc_a, SE of difference, two-sided p).

    Verbatim from scripts/run_true_filter_ablation.py so that recomputed
    intersection tests match the published ones exactly.
    """
    pos = y == 1
    neg = ~pos
    m, n = int(pos.sum()), int(neg.sum())
    if m < 2 or n < 2:
        raise ValueError("need >=2 positives and >=2 negatives")

    def structural(scores):
        sp, sn = scores[pos], scores[neg]
        v10 = np.array([np.mean((sp[i] > sn) + 0.5 * (sp[i] == sn)) for i in range(m)])
        v01 = np.array([np.mean((sp > sn[j]) + 0.5 * (sp == sn[j])) for j in range(n)])
        return v10, v01, float(v10.mean())

    v10a, v01a, auc_a = structural(scores_a)
    v10b, v01b, auc_b = structural(scores_b)
    s10 = np.cov(np.vstack([v10a, v10b]))
    s01 = np.cov(np.vstack([v01a, v01b]))
    var = (s10[0, 0] + s10[1, 1] - 2 * s10[0, 1]) / m + (s01[0, 0] + s01[1, 1] - 2 * s01[0, 1]) / n
    se = math.sqrt(var) if var > 0 else float("nan")
    delta = auc_b - auc_a
    p = 2 * (1 - norm.cdf(abs(delta / se))) if se and not math.isnan(se) and se > 0 else float("nan")
    return delta, se, p, auc_a, auc_b


# ----------------------------------------------------------------- BH + DL
def bh_adjust(pvals):
    p = np.asarray(pvals, dtype=float)
    ok = ~np.isnan(p)
    out = np.full_like(p, np.nan)
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.minimum(adj, 1.0)
    res = np.empty(n)
    res[order] = adj
    out[ok] = res
    return out


def dersimonian_laird(effects, ses):
    y = np.asarray(effects, dtype=float)
    s = np.asarray(ses, dtype=float)
    ok = ~(np.isnan(y) | np.isnan(s) | (s <= 0))
    y, s = y[ok], s[ok]
    w = 1.0 / s**2
    mu_fixed = np.sum(w * y) / np.sum(w)
    Q = float(np.sum(w * (y - mu_fixed) ** 2))
    k = len(y)
    df = k - 1
    C = np.sum(w) - np.sum(w**2) / np.sum(w)
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    wr = 1.0 / (s**2 + tau2)
    mu = float(np.sum(wr * y) / np.sum(wr))
    se_mu = float(math.sqrt(1.0 / np.sum(wr)))
    z = mu / se_mu
    p = float(2 * (1 - norm.cdf(abs(z))))
    i2 = float(max(0.0, (Q - df) / Q) * 100) if Q > 0 else 0.0
    q_p = float(1 - __import__("scipy.stats", fromlist=["chi2"]).chi2.cdf(Q, df)) if df > 0 else float("nan")
    return {"k": k, "pooled": mu, "se": se_mu, "ci_lo": mu - 1.96 * se_mu,
            "ci_hi": mu + 1.96 * se_mu, "z": float(z), "p": p,
            "Q": Q, "df": df, "Q_p": q_p, "I2": i2, "tau2": float(tau2)}
