# Proxy Filter-Ablation Report

**Method:** For each target, replaced PAINS/Brenk-filtered compounds' missing skill scores with their naive-protocol scores, then recomputed ROC AUC on the full library. This is a proxy for a true ablation — the substituted scores were generated under naive docking parameters, not skill — but is informative about whether the AChE/COX-2 AUC gain is filter-driven (active reshaping) or docking-driven.

**Interpretation:**

- If `delta_proxy_minus_naive` is close to zero, the skill AUC gain *was* primarily driven by filtering out hard-to-dock actives.
- If `delta_proxy_minus_naive` is similar to `delta_post_minus_naive`, the gain *survives* filter ablation and is at least partly attributable to skill docking parameters.

| Target | Naive AUC (full) | Skill AUC (post-filter) | Proxy full AUC | Δ proxy − naive | Δ post − naive |
|---|---:|---:|---:|---:|---:|
| ache | 0.444 | 0.504 | 0.476 | 0.032 | 0.061 |
| bace1 | 0.626 | 0.528 | 0.580 | -0.047 | -0.098 |
| cdk2 | 0.572 | 0.551 | 0.571 | -0.001 | -0.021 |
| cox2 | 0.607 | 0.653 | 0.610 | 0.003 | 0.046 |
| dpp4 | 0.544 | 0.508 | 0.543 | -0.002 | -0.037 |
| egfr | 0.575 | 0.630 | 0.577 | 0.002 | 0.055 |
| esr1 | 0.564 | 0.602 | 0.547 | -0.016 | 0.038 |
| fpr2 | 0.719 | 0.685 | 0.701 | -0.018 | -0.034 |
| hsp90 | 0.533 | 0.534 | 0.515 | -0.018 | 0.000 |
| p38 | 0.657 | 0.649 | 0.651 | -0.006 | -0.008 |
| thrombin | 0.695 | 0.633 | 0.647 | -0.048 | -0.062 |

## ACHE spotlight

- Skill post-filter Δ vs naive: **+0.061** (this was the original BIB-submitted gain)
- Skill proxy full-library Δ vs naive: **+0.032**
- Effect shrinks toward zero by 47% under proxy ablation.

## COX2 spotlight

- Skill post-filter Δ vs naive: **+0.046** (this was the original BIB-submitted gain)
- Skill proxy full-library Δ vs naive: **+0.003**
- Effect shrinks toward zero by 94% under proxy ablation.

**Caveat:** This is a *proxy* analysis. The substituted scores were generated under naive docking parameters (20 Å box, exhaustiveness 8), so any preserved Δ in the proxy may underestimate the true ablation effect. A full ablation — re-docking the PAINS/Brenk-filtered compounds with skill docking parameters — remains required for a definitive answer. This proxy is intended only to inform whether the full ablation is worth running (yes if the proxy shows substantial preservation; less informative if the proxy already collapses the effect to zero).
