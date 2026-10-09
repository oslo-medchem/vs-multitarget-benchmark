#!/usr/bin/env python3
"""Decoy-contamination sensitivity, filtering's effect on composition, and an
independent hit-quality endpoint.

Three referee requests are answered here.

1. Referee 1 point 3 and referee 2 point 5, sensitivity to the constructed
   decoy population. The intersection DeLong test is recomputed after removing
   every decoy that is a qualifying ChEMBL active (pChEMBL >= 5) of any other
   target in the panel. If the conclusions are unchanged, the shared-pool
   concern is answered quantitatively rather than by disclaimer.

2. Referee 2 point 5, whether filtering materially changes the benchmark's
   chemical composition. Property distributions are compared before and after
   the PAINS/Brenk filter, separately for actives and decoys.

3. Referee 2, significance. The claim that the skill file improves library
   curation needs an endpoint beyond demonstrating that the filter runs. Two
   are computed. First, the chemical-liability composition of the nominal hit
   list each protocol would hand to a chemist, taken as the top 1% of its
   end-to-end ranking. Second, following Capuzzi et al., the reported potency
   of the curated actives the filter removes, which tests whether the filter
   discards weak chemistry or well-characterised potent chemistry.

Membership of the PAINS/Brenk alert set is read from the benchmark's own
filter_report.csv rather than re-derived, so it is exactly the filter the skill
protocol applied.

Outputs decoy_sensitivity.csv, filter_composition.csv, hitlist_quality.csv and
filtered_active_potency.csv in data/revision_2026/.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, QED, Crippen

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import revision_common as rc

RDLogger.DisableLog("rdApp.*")


def props(smi: str):
    m = Chem.MolFromSmiles(smi) if smi else None
    if m is None:
        return None
    return {"MW": Descriptors.MolWt(m), "cLogP": Crippen.MolLogP(m),
            "HBD": Descriptors.NumHDonors(m), "HBA": Descriptors.NumHAcceptors(m),
            "RotB": Descriptors.NumRotatableBonds(m), "TPSA": Descriptors.TPSA(m),
            "QED": QED.qed(m), "HeavyAtoms": m.GetNumHeavyAtoms()}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--out-dir", default="data/revision_2026")
    args = ap.parse_args()
    root, out = Path(args.data_root), Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # qualifying actives per target, for the contamination screen
    qualifying = {}
    potency = {}
    for t in rc.TARGETS:
        p = root / "targets" / t / "01_active_curation" / "actives_raw.csv"
        ids, best = set(), {}
        if p.exists():
            with p.open() as f:
                for r in csv.DictReader(f):
                    mid = r["molecule_chembl_id"]
                    ids.add(mid)
                    try:
                        v = float(r["pchembl_value"])
                    except (ValueError, KeyError, TypeError):
                        continue
                    if mid not in best or v > best[mid]:
                        best[mid] = v
        qualifying[t] = ids
        potency[t] = best

    sens_rows, comp_rows, hit_rows, pot_rows = [], [], [], []

    for t in rc.TARGETS:
        td = rc.load_target(t, root)
        alert = set(td.filtered)          # compounds carrying a PAINS/Brenk alert

        # ---------- 1. decoy-contamination sensitivity --------------------
        contaminated = set()
        for cid, lab in td.labels.items():
            if lab != 0:
                continue
            sid = td.source_id.get(cid, "")
            if not sid:
                continue
            if any(sid in qualifying[o] for o in rc.TARGETS if o != t):
                contaminated.add(cid)

        common = sorted(set(td.naive) & set(td.skill) & set(td.labels))
        y = np.array([td.labels[c] for c in common])
        rn = np.array([-td.naive[c] for c in common])
        rs = np.array([-td.skill[c] for c in common])
        d_all, se_all, p_all, _, _ = rc.delong_paired(y, rn, rs)

        keep = [i for i, c in enumerate(common) if c not in contaminated]
        yk, rnk, rsk = y[keep], rn[keep], rs[keep]
        d_cl, se_cl, p_cl, _, _ = rc.delong_paired(yk, rnk, rsk)

        sens_rows.append({
            "target": t, "intersection_n_all": len(common),
            "n_contaminated_decoys_removed": len(common) - len(keep),
            "intersection_n_clean": len(keep),
            "delta_all_decoys": round(d_all, 6), "p_all_decoys": round(p_all, 6),
            "delta_clean_decoys": round(d_cl, 6), "p_clean_decoys": round(p_cl, 6),
            "delta_shift": round(d_cl - d_all, 6),
        })

        # ---------- 2. composition before/after filtering ------------------
        for lab_name, lab_val in (("active", 1), ("decoy", 0)):
            pre, post = [], []
            for cid, lv in td.labels.items():
                if lv != lab_val:
                    continue
                pr = props(td.smiles.get(cid, ""))
                if pr is None:
                    continue
                pre.append(pr)
                if cid not in alert:
                    post.append(pr)
            if not pre:
                continue
            row = {"target": t, "set": lab_name, "n_before": len(pre), "n_after": len(post),
                   "pct_removed": round(100 * (1 - len(post) / len(pre)), 2)}
            for k in ("MW", "cLogP", "HBD", "HBA", "RotB", "TPSA", "QED"):
                a = np.array([x[k] for x in pre], dtype=float)
                b = np.array([x[k] for x in post], dtype=float) if post else np.array([np.nan])
                row[f"{k}_median_before"] = round(float(np.median(a)), 3)
                row[f"{k}_median_after"] = round(float(np.nanmedian(b)), 3)
                row[f"{k}_shift"] = round(float(np.nanmedian(b) - np.median(a)), 3)
            comp_rows.append(row)

        # ---------- 3a. hit-list chemical liability -----------------------
        for prot, sc in (("naive", td.naive), ("skill", td.skill)):
            cids, yv, rank = rc.end_to_end_ranks(td.labels, sc, "strict")
            n_top = max(1, int(0.01 * len(cids)))
            order = np.argsort(-rank, kind="stable")
            top = [cids[i] for i in order[:n_top]]
            n_alert = sum(1 for c in top if c in alert)
            n_act = sum(1 for c in top if td.labels[c] == 1)
            qs = [props(td.smiles.get(c, "")) for c in top]
            qs = [x["QED"] for x in qs if x is not None]
            hit_rows.append({
                "target": t, "protocol": prot, "hitlist_size_top1pct": n_top,
                "n_with_PAINS_or_Brenk_alert": n_alert,
                "pct_with_alert": round(100 * n_alert / n_top, 2),
                "n_true_actives_recovered": n_act,
                "median_QED": round(float(np.median(qs)), 4) if qs else "",
            })

        # ---------- 3b. potency of the actives the filter removes ---------
        act_flagged, act_kept = [], []
        for cid, lv in td.labels.items():
            if lv != 1:
                continue
            sid = td.source_id.get(cid, "")
            v = potency[t].get(sid)
            if v is None:
                continue
            (act_flagged if cid in alert else act_kept).append(v)
        if act_flagged and act_kept:
            from scipy.stats import mannwhitneyu
            u, pu = mannwhitneyu(act_flagged, act_kept, alternative="two-sided")
            pot_rows.append({
                "target": t,
                "n_actives_flagged": len(act_flagged), "n_actives_kept": len(act_kept),
                "median_pChEMBL_flagged": round(float(np.median(act_flagged)), 3),
                "median_pChEMBL_kept": round(float(np.median(act_kept)), 3),
                "difference": round(float(np.median(act_flagged) - np.median(act_kept)), 3),
                "mannwhitney_p": round(float(pu), 5),
            })

    padj = rc.bh_adjust([r["mannwhitney_p"] for r in pot_rows]) if pot_rows else []
    for r, a in zip(pot_rows, padj):
        r["mannwhitney_p_adj_BH"] = round(float(a), 5)

    for name, rows in (("decoy_sensitivity.csv", sens_rows),
                       ("filter_composition.csv", comp_rows),
                       ("hitlist_quality.csv", hit_rows),
                       ("filtered_active_potency.csv", pot_rows)):
        if not rows:
            continue
        with (out / name).open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {out/name} ({len(rows)} rows)")

    print("\n=== decoy-contamination sensitivity (intersection DeLong) ===")
    for r in sens_rows:
        print(f"  {r['target']:9s} removed={r['n_contaminated_decoys_removed']:3d}  "
              f"d_all={r['delta_all_decoys']:+.4f} (p={r['p_all_decoys']:.4f})  "
              f"d_clean={r['delta_clean_decoys']:+.4f} (p={r['p_clean_decoys']:.4f})  "
              f"shift={r['delta_shift']:+.4f}")
    print("\n=== potency of curated actives removed by the filter ===")
    for r in pot_rows:
        print(f"  {r['target']:9s} flagged n={r['n_actives_flagged']:3d} med={r['median_pChEMBL_flagged']:.2f}  "
              f"kept n={r['n_actives_kept']:3d} med={r['median_pChEMBL_kept']:.2f}  "
              f"diff={r['difference']:+.2f} p={r['mannwhitney_p']:.4f} p_adj={r.get('mannwhitney_p_adj_BH')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
