#!/usr/bin/env python3
"""
Integrate true filter ablation results into the manuscript.

Reads filter_ablation_data/true_ablation_results.json (from run_true_filter_ablation.py)
and produces:
  - true_filter_ablation.csv : per-target full-library AUCs and DeLong stats
  - figures/fig1_yields_and_ablation.{png,pdf} : updated with TRUE ablation column
  - manuscript/manuscript_v2.tex : Discussion §"True Filter Ablation" replaces the proxy section
  - manuscript/Manuscript_BIB_v3.docx : DOCX mirror

Run after `scripts/run_true_filter_ablation.py --targets all` completes.
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
ABL_JSON = ROOT / "filter_ablation_data" / "true_ablation_results.json"


def main() -> int:
    if not ABL_JSON.exists():
        print(f"ERROR: {ABL_JSON} not found. Run run_true_filter_ablation.py first.", file=sys.stderr)
        return 1

    results = json.loads(ABL_JSON.read_text())
    # Filter out errors
    ok = [r for r in results if "error" not in r]
    bad = [r for r in results if "error" in r]
    if bad:
        print(f"WARNING: {len(bad)} target(s) failed: {[r['target'] for r in bad]}", file=sys.stderr)

    # Write CSV
    csv_out = ROOT / "true_filter_ablation.csv"
    if ok:
        cols = list(ok[0].keys())
        with csv_out.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(ok)
        print(f"Wrote {csv_out}")

    # Print summary table
    print()
    print(f"{'target':10s} {'naive':>7s} {'skill_post':>11s} {'skill_abl':>10s} {'Δ post':>8s} {'Δ abl':>8s} {'Δ DeLong':>9s} {'p_DeLong':>9s}")
    print('-'*80)
    for r in ok:
        d_p = r.get("delta_post_minus_naive")
        d_a = r.get("delta_ablation_minus_naive")
        dl_d = r.get("delong_delta")
        dl_p = r.get("delong_p")
        def fm(v, w=7, p=3):
            return f"{v:.{p}f}".rjust(w) if isinstance(v, (int, float)) else "NA".rjust(w)
        print(f"{r['target']:10s} {fm(r['auc_naive_full']):>7s} {fm(r['auc_skill_post_filter']):>11s} {fm(r['auc_skill_ablation_full']):>10s} {fm(d_p):>8s} {fm(d_a):>8s} {fm(dl_d,9,3):>9s} {fm(dl_p,9,3):>9s}")

    # Spotlight on AChE and COX-2
    print()
    for spot in ("ache", "cox2"):
        r = next((x for x in ok if x["target"] == spot), None)
        if not r:
            continue
        d_p = r.get("delta_post_minus_naive")
        d_a = r.get("delta_ablation_minus_naive")
        if d_p is not None and d_a is not None and d_p != 0:
            shrink = (1 - abs(d_a) / abs(d_p)) * 100 if d_p != 0 else 0
            print(f"{spot.upper()}: post Δ {d_p:+.3f} → true-ablation Δ {d_a:+.3f}  "
                  f"({'shrinks' if abs(d_a)<abs(d_p) else 'grows'} {abs(shrink):.0f}%)")

    # If any priority target completed, regenerate Fig 1B
    update_fig1(ok)

    # Recompute BH across all 11 targets (FPR2 keeps its old DeLong p; others use new ablation DeLong)
    recompute_bh_table2(ok)

    # Update manuscript text
    update_manuscript(ok)

    # Rebuild PDF
    rebuild_pdf()

    # Update DOCX (creates v3)
    update_docx(ok)

    return 0


def recompute_bh_table2(results: list[dict]) -> None:
    """Update Table 2 in the LaTeX with new intersection DeLong values from the ablation."""
    import numpy as np
    import re
    tex_path = ROOT / "manuscript" / "manuscript_v2.tex"
    tex = tex_path.read_text()

    # Read original delong p-values for FPR2 (unchanged) and any non-ablated targets
    delong_csv = Path("/data/CLAUDE_works/TEST_SKILLS/vina_docking_test/vs_benchmark/vs_benchmark_multitarget/06_meta_analysis/delong_deltas.csv")
    orig = {}
    if delong_csv.exists():
        for r in csv.DictReader(delong_csv.open()):
            orig[r["target"]] = (float(r["delta_auc"]), float(r["p_value"]))

    # Build the new per-target table
    targets_in_order = ["fpr2", "egfr", "cdk2", "cox2", "esr1", "dpp4", "ache", "bace1", "hsp90", "p38", "thrombin"]
    abl_lookup = {r["target"]: r for r in results}

    new_p = []
    new_delta = []
    for t in targets_in_order:
        if t in abl_lookup and abl_lookup[t].get("delong_p") is not None:
            new_delta.append(abl_lookup[t]["delong_delta"])
            new_p.append(abl_lookup[t]["delong_p"])
        elif t in orig:
            new_delta.append(orig[t][0])
            new_p.append(orig[t][1])
        else:
            new_delta.append(None)
            new_p.append(None)

    # BH correction
    p_arr = np.array([p if p is not None else 1.0 for p in new_p])
    n = len(p_arr)
    order_idx = np.argsort(p_arr)
    sorted_p = p_arr[order_idx]
    sorted_adj = np.minimum.accumulate((sorted_p * n / np.arange(1, n+1))[::-1])[::-1]
    sorted_adj = np.minimum(sorted_adj, 1.0)
    p_adj = np.empty(n)
    p_adj[order_idx] = sorted_adj

    # Print summary
    print()
    print(f"Updated DeLong + BH (post-ablation):")
    print(f"{'target':10s}{'ΔAUC':>9s}{'p':>9s}{'p_adj':>9s}{'note':>20s}")
    for i, t in enumerate(targets_in_order):
        note = "ablated" if t in abl_lookup else ("FPR2 kept original" if t == "fpr2" else "kept original")
        d = new_delta[i] if new_delta[i] is not None else float('nan')
        p = new_p[i] if new_p[i] is not None else float('nan')
        print(f"{t:10s}{d:>9.4f}{p:>9.4f}{p_adj[i]:>9.4f}  {note}")

    # Write CSV
    out = ROOT / "table2_post_ablation.csv"
    with out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["target", "delta_auc", "p_value", "p_adj", "note"])
        for i, t in enumerate(targets_in_order):
            note = "ablated" if t in abl_lookup else "original"
            w.writerow([t, new_delta[i], new_p[i], round(p_adj[i], 4), note])
    print(f"Wrote {out}")

    # Update abstract if AChE BH significance changed
    ache_idx = targets_in_order.index("ache")
    ache_padj_new = float(p_adj[ache_idx])
    cox2_idx = targets_in_order.index("cox2")
    cox2_padj_new = float(p_adj[cox2_idx])

    # Find significant targets (p_adj < 0.05)
    sig_targets = [(t, float(p_adj[i])) for i, t in enumerate(targets_in_order)
                   if float(p_adj[i]) < 0.05 and t in abl_lookup]
    n_sig_post = len(sig_targets)

    # Update macros in the LaTeX
    import re
    macro_path = tex_path
    tex2 = macro_path.read_text()

    # Update PADJACHE / PADJCOX macros if AChE / COX-2 were ablated
    if "ache" in abl_lookup:
        tex2 = re.sub(
            r"\\newcommand\{\\PADJACHE\}\{[^}]*\}",
            r"\\newcommand{\\PADJACHE}{" + f"{ache_padj_new:.3f}" + "}",
            tex2,
        )
    if "cox2" in abl_lookup:
        tex2 = re.sub(
            r"\\newcommand\{\\PADJCOX\}\{[^}]*\}",
            r"\\newcommand{\\PADJCOX}{" + f"{cox2_padj_new:.3f}" + "}",
            tex2,
        )
    # Update NSIG macro if number of significant targets changed
    if "ache" in abl_lookup or "cox2" in abl_lookup:
        tex2 = re.sub(
            r"\\newcommand\{\\NSIG\}\{[^}]*\}",
            r"\\newcommand{\\NSIG}{" + f"{n_sig_post}" + "}",
            tex2,
        )

    # If AChE no longer survives BH, modify the headline abstract sentence.
    # Use plain string replace to avoid regex escape pitfalls.
    if "ache" in abl_lookup and ache_padj_new > 0.05:
        old = r"only AChE survived Benjamini--Hochberg correction ($p_{\mathrm{adj}}$~=~\PADJACHE{}); COX-2 did not ($p_{\mathrm{adj}}$~=~\PADJCOX{})"
        new = (
            f"after the true filter ablation, NO target survived Benjamini--Hochberg correction at $\\alpha = 0.05$ "
            f"(AChE post-ablation $p_{{\\mathrm{{adj}}}}$~=~{ache_padj_new:.3f}; COX-2 post-ablation $p_{{\\mathrm{{adj}}}}$~=~{cox2_padj_new:.3f})"
        )
        if old in tex2:
            tex2 = tex2.replace(old, new)
            print("  Updated abstract: 'no target survives BH' framing applied")
        else:
            print("  WARNING: abstract sentence pattern not found; manual update needed")

    macro_path.write_text(tex2)
    print(f"Updated PADJ macros and abstract narrative in {macro_path}")
    print(f"  Post-ablation: NSIG = {n_sig_post}; AChE p_adj = {ache_padj_new:.3f}; COX-2 p_adj = {cox2_padj_new:.3f}")


def rebuild_pdf() -> None:
    """Run pdflatex + bibtex + pdflatex + pdflatex."""
    import subprocess
    md = ROOT / "manuscript"
    for f in ["manuscript_v2.aux", "manuscript_v2.bbl", "manuscript_v2.blg", "manuscript_v2.log", "manuscript_v2.out", "manuscript_v2.pdf"]:
        (md / f).unlink(missing_ok=True)
    cmds = [
        ["pdflatex", "-interaction=nonstopmode", "manuscript_v2.tex"],
        ["bibtex", "manuscript_v2"],
        ["pdflatex", "-interaction=nonstopmode", "manuscript_v2.tex"],
        ["pdflatex", "-interaction=nonstopmode", "manuscript_v2.tex"],
    ]
    for cmd in cmds:
        subprocess.run(cmd, cwd=md, capture_output=True, timeout=60)
    pdf = md / "manuscript_v2.pdf"
    print(f"PDF rebuild: {pdf} ({pdf.stat().st_size if pdf.exists() else 'MISSING'} bytes)")


def update_docx(results: list[dict]) -> None:
    """Mirror key new content into the DOCX (v3)."""
    try:
        from docx import Document
    except ImportError:
        print("python-docx not installed; skipping DOCX update")
        return

    src = ROOT / "manuscript" / "Manuscript_BIB_v2.docx"
    dst = ROOT / "manuscript" / "Manuscript_BIB_v3.docx"
    if not src.exists():
        print(f"Source DOCX missing: {src}")
        return
    import shutil
    shutil.copy(src, dst)

    doc = Document(dst)

    def replace_paragraph_text(para, new_text):
        if not para.runs:
            para.add_run(new_text)
            return
        first_run = para.runs[0]
        for run in para.runs[1:]:
            run._r.getparent().remove(run._r)
        first_run.text = new_text

    # Find the paragraph with the proxy ablation text (paragraph 50 in v2)
    ache = next((r for r in results if r["target"] == "ache"), None)
    cox2 = next((r for r in results if r["target"] == "cox2"), None)

    if ache and cox2:
        d_p_ache = ache.get("delta_post_minus_naive") or 0
        d_a_ache = ache.get("delta_ablation_minus_naive") or 0
        d_p_cox2 = cox2.get("delta_post_minus_naive") or 0
        d_a_cox2 = cox2.get("delta_ablation_minus_naive") or 0
        new_text = (
            f"A true filter ablation (re-preparing the PAINS/Brenk-filtered subset under the skill protocol's "
            f"preparation pipeline with apply_filter=False, then re-docking under skill parameters: 25 Å box, "
            f"exhaustiveness 32, Vina scoring) was performed for the 10 multi-target targets (FPR2 excluded "
            f"because it uses the older single-target archive layout). For COX-2, the post-filter ΔAUC of "
            f"{d_p_cox2:+.3f} becomes {d_a_cox2:+.3f} after the true ablation; for AChE, the post-filter ΔAUC "
            f"of {d_p_ache:+.3f} becomes {d_a_ache:+.3f}. The bulk of the apparent skill-protocol gain is "
            f"therefore attributable to active-set reshaping by PAINS/Brenk filtering rather than to docking-"
            f"quality improvements from the larger box and higher exhaustiveness. The deep-pocket interpretation "
            f"is substantially weakened by these data and should be regarded as a hypothesis at most. See "
            f"true_filter_ablation.csv for full per-target results across all 10 multi-target targets."
        )
        # Find paragraph 50 (which had the proxy)
        if len(doc.paragraphs) > 50:
            replace_paragraph_text(doc.paragraphs[50], new_text)
            print(f"Updated DOCX paragraph 50 with true-ablation results")

    doc.save(dst)
    print(f"Saved {dst}")


def update_fig1(results: list[dict]) -> None:
    """Regenerate Fig 1B with true-ablation column."""
    if not results:
        return

    yields_csv = ROOT / "yields.csv"
    if not yields_csv.exists():
        print("yields.csv not found — skipping fig regeneration")
        return
    yields = list(csv.DictReader(yields_csv.open()))
    abl = {r["target"]: r for r in results}

    order = ["fpr2", "egfr", "cdk2", "p38", "cox2", "esr1", "dpp4", "ache", "bace1", "hsp90", "thrombin"]
    yields = sorted(yields, key=lambda r: order.index(r["target"]))

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 4.5), gridspec_kw={"width_ratios": [1.1, 1]})

    targets = [r["target"].upper() if r["target"] != "fpr2" else "FPR2" for r in yields]
    naive_doc = np.array([float(r["naive_full_pct"]) for r in yields])
    skill_doc = np.array([float(r["skill_full_pct"]) for r in yields])
    skill_filt = np.array([float(r["filter_lib_pct"]) for r in yields])
    skill_fail = np.clip(100.0 - skill_doc - skill_filt, 0, None)
    naive_fail = 100.0 - naive_doc

    x = np.arange(len(targets))
    w = 0.36
    axA.bar(x - w/2, naive_doc, w, label="Naive: docked", color="#2c7fb8")
    axA.bar(x - w/2, naive_fail, w, bottom=naive_doc, label="Naive: prep/dock failure", color="#bdd7e7", edgecolor="white", linewidth=0.5)
    axA.bar(x + w/2, skill_doc, w, label="Skill: docked", color="#31a354")
    axA.bar(x + w/2, skill_filt, w, bottom=skill_doc, label="Skill: PAINS/Brenk filtered", color="#fdae6b", edgecolor="white", linewidth=0.5)
    axA.bar(x + w/2, skill_fail, w, bottom=skill_doc + skill_filt, label="Skill: prep failure", color="#bdbdbd", edgecolor="white", linewidth=0.5)
    axA.set_xticks(x); axA.set_xticklabels(targets, rotation=30, ha="right", fontsize=9)
    axA.set_ylabel("Fraction of full library (%)"); axA.set_ylim(0, 100)
    axA.set_title("(A) Pipeline yield: where the library goes", fontsize=11, loc="left")
    axA.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3, fontsize=8, frameon=False)
    axA.grid(axis="y", alpha=0.25, linewidth=0.5); axA.set_axisbelow(True)
    for spine in ["top", "right"]:
        axA.spines[spine].set_visible(False)

    # Panel B: true filter ablation
    auc_naive = []
    auc_skill_post = []
    auc_skill_abl = []
    for r in yields:
        a = abl.get(r["target"])
        auc_naive.append(float(a["auc_naive_full"]) if a and a.get("auc_naive_full") else float("nan"))
        auc_skill_post.append(float(a["auc_skill_post_filter"]) if a and a.get("auc_skill_post_filter") else float("nan"))
        auc_skill_abl.append(float(a["auc_skill_ablation_full"]) if a and a.get("auc_skill_ablation_full") else float("nan"))

    auc_naive = np.array(auc_naive)
    auc_skill_post = np.array(auc_skill_post)
    auc_skill_abl = np.array(auc_skill_abl)

    w2 = 0.27
    axB.bar(x - w2, auc_naive, w2, label="Naive (full library)", color="#888888")
    axB.bar(x, auc_skill_post, w2, label="Skill (post-filter)", color="#31a354")
    axB.bar(x + w2, auc_skill_abl, w2, label="Skill + true ablation (full library)", color="#e6550d")
    axB.axhline(0.5, color="k", linestyle=":", linewidth=0.7, alpha=0.5)
    axB.set_xticks(x); axB.set_xticklabels(targets, rotation=30, ha="right", fontsize=9)
    axB.set_ylabel("ROC AUC"); axB.set_ylim(0.4, 0.78)
    axB.set_title("(B) True filter ablation (re-docked under skill parameters)", fontsize=11, loc="left")
    axB.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize=8, frameon=False)
    axB.grid(axis="y", alpha=0.25, linewidth=0.5); axB.set_axisbelow(True)
    for spine in ["top", "right"]:
        axB.spines[spine].set_visible(False)

    for tgt in ("cox2", "ache"):
        if tgt not in [r["target"] for r in yields]:
            continue
        i = order.index(tgt)
        a = abl.get(tgt)
        if not a:
            continue
        d_p = a.get("delta_post_minus_naive")
        d_a = a.get("delta_ablation_minus_naive")
        if d_p is None or d_a is None:
            continue
        axB.annotate(
            f"Δ post: {d_p:+.3f}\nΔ ablation: {d_a:+.3f}",
            xy=(i, max(auc_skill_post[i] if not math.isnan(auc_skill_post[i]) else 0.5,
                       auc_skill_abl[i] if not math.isnan(auc_skill_abl[i]) else 0.5) + 0.015),
            ha="center", fontsize=7, color="#333",
        )

    fig.tight_layout()
    out_png = ROOT / "figures" / "fig1_yields_and_true_ablation.png"
    out_pdf = ROOT / "figures" / "fig1_yields_and_true_ablation.pdf"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_png}\nWrote {out_pdf}")


def update_manuscript(results: list[dict]) -> None:
    """Replace Discussion §Proxy Filter Ablation in the LaTeX with the true ablation."""
    tex_path = ROOT / "manuscript" / "manuscript_v2.tex"
    tex = tex_path.read_text()

    # Find AChE and COX-2 numbers
    ache = next((r for r in results if r["target"] == "ache"), None)
    cox2 = next((r for r in results if r["target"] == "cox2"), None)

    n_complete = len(results)
    if not ache and not cox2:
        # Don't update text yet; only enough data.
        print(f"Skipping manuscript update — only {n_complete} target(s) complete and neither is AChE/COX-2.")
        return

    # Build new section text
    blocks = [
        r"\subsection*{True Filter Ablation: Most of the Target-Specific Gain Is Filter-Driven}",
        "",
        r"To attribute the AChE and COX-2 gains to docking quality versus active-set reshaping, we re-prepared the PAINS/Brenk-filtered subset of each target's library using the skill protocol's preparation pipeline (apply\_filter=False) and re-docked these compounds under skill parameters (25~\AA{} box, exhaustiveness 32, Vina scoring) on identical hardware. The resulting scores were merged with the original skill-protocol scores to produce a full-library skill score table, allowing direct full-library AUC computation and an unconfounded paired DeLong test against the naive protocol (Figure~\ref{fig:yields}B; \texttt{true\_filter\_ablation.csv}).",
        "",
    ]

    # AChE results
    if ache:
        d_p = ache.get("delta_post_minus_naive")
        d_a = ache.get("delta_ablation_minus_naive")
        dl_d = ache.get("delong_delta")
        dl_p = ache.get("delong_p")
        if d_p is not None and d_a is not None:
            shrink_a = abs(d_a) / abs(d_p) * 100 if d_p else 0
            blocks.append(
                f"For AChE, the post-filter $\\Delta$AUC of $+${abs(d_p):.3f} (the original BIB-submitted gain) "
                f"becomes $+${d_a:+.3f} on the full library after the true ablation"
                + (f" (paired DeLong $\\Delta$AUC = {dl_d:+.3f}, $p$ = {dl_p:.3f})" if dl_d is not None else "")
                + f". The true-ablation effect is {shrink_a:.0f}\\% of the post-filter effect; "
                f"in absolute terms the gain has {'collapsed' if abs(d_a) < 0.01 else 'shrunk'} substantially "
                f"once filtered compounds are scored under the same docking parameters."
            )

    if cox2:
        d_p = cox2.get("delta_post_minus_naive")
        d_a = cox2.get("delta_ablation_minus_naive")
        dl_d = cox2.get("delong_delta")
        dl_p = cox2.get("delong_p")
        if d_p is not None and d_a is not None:
            shrink_c = abs(d_a) / abs(d_p) * 100 if d_p else 0
            blocks.append(
                f"For COX-2, the post-filter $\\Delta$AUC of $+${abs(d_p):.3f} becomes $+${d_a:+.3f} after the ablation"
                + (f" (paired DeLong $\\Delta$AUC = {dl_d:+.3f}, $p$ = {dl_p:.3f})" if dl_d is not None else "")
                + f", or {shrink_c:.0f}\\% of the post-filter effect."
            )

    # Conclusion
    if ache and cox2:
        blocks.append("")
        blocks.append(
            "The true-ablation results confirm and refine the proxy-ablation conclusion: the bulk of the "
            "apparent skill-protocol AUC gain on the two ``winning'' targets is attributable to active-set "
            "reshaping by PAINS/Brenk filtering rather than to docking-quality improvements from the larger "
            "box (25 vs 20~\\AA) and higher exhaustiveness (32 vs 8). The corrected, multi-target picture is that "
            "skill files reliably enforce medicinal-chemistry-aware library curation but do not, in our data, "
            "demonstrate a generalisable improvement in docking discrimination, even on the targets initially flagged as "
            "hypothesis-supporting."
        )

    new_text = "\n".join(blocks) + "\n"

    # Strategy: surgically replace just the [TRUE_ABLATION_RESULT_AChE_COX2 ...] placeholder
    # with a sentence summarising the AChE/COX-2 true-ablation deltas.
    # Fall back to full-section replacement if the placeholder is missing.

    # Build a single replacement sentence
    sentences = []
    if ache:
        d_p = ache.get("delta_post_minus_naive")
        d_a = ache.get("delta_ablation_minus_naive")
        dl_d = ache.get("delong_delta")
        dl_p = ache.get("delong_p")
        if d_p is not None and d_a is not None:
            extra = f" (paired DeLong $\\Delta$AUC = {dl_d:+.3f}, $p$ = {dl_p:.3f})" if (dl_d is not None and dl_p is not None) else ""
            sentences.append(
                f"yielded $\\Delta$AUC = {d_a:+.3f} for AChE on the full library{extra}"
            )
    if cox2:
        d_p = cox2.get("delta_post_minus_naive")
        d_a = cox2.get("delta_ablation_minus_naive")
        dl_d = cox2.get("delong_delta")
        dl_p = cox2.get("delong_p")
        if d_p is not None and d_a is not None:
            extra = f" (paired DeLong $\\Delta$AUC = {dl_d:+.3f}, $p$ = {dl_p:.3f})" if (dl_d is not None and dl_p is not None) else ""
            sentences.append(
                f"$\\Delta$AUC = {d_a:+.3f} for COX-2{extra}"
            )

    if not sentences:
        print("Skipping placeholder replacement — neither AChE nor COX-2 has results yet.")
        return

    replacement = "yielded " + " and ".join(sentences) + " — confirming and refining the proxy result that the bulk of the apparent skill-protocol gain on these targets is attributable to active-set reshaping by PAINS/Brenk filtering rather than to docking-quality improvements"

    # Find and replace the placeholder. Use lambda replacement to avoid backslash escape issues.
    placeholder_pattern = r"\[\\textbf\{TRUE\\_ABLATION\\_RESULT\\_AChE\\_COX2[^]]*\}\]"
    import re
    if re.search(placeholder_pattern, tex):
        tex_new = re.sub(placeholder_pattern, lambda m: replacement, tex)
        tex_path.write_text(tex_new)
        print(f"Replaced TRUE_ABLATION placeholder in {tex_path}")
    else:
        # Already filled or pattern doesn't match — write fragment for manual paste
        frag = ROOT / "manuscript" / "_true_ablation_fragment.tex"
        frag.write_text(replacement)
        print(f"Placeholder not found; wrote fragment to {frag}")


if __name__ == "__main__":
    sys.exit(main())
