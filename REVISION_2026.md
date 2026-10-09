# Revision analyses, cross conditions and environment (October 2026)

This directory documents what was added to the repository for the Journal of
Cheminformatics version of the manuscript. It supersedes nothing: the original
pipeline in `common/`, the prompts in `prompts/` and the skill file `SKILL.md`
are unchanged and remain the record of what the agent generated.

## What is here

| Path | Contents |
|---|---|
| `scripts/revision_2026/` | The analysis scripts that produce every number in the article. `revision_common.py` holds the loaders, the tie-corrected metric panel, the paired DeLong test, Benjamini-Hochberg adjustment and the DerSimonian-Laird pooling. |
| `cross_conditions/` | Protocol, node-side runner, preparation and docking driver, and analysis script for conditions D and E, the two cross conditions that separate preparation tools from docking parameters. |
| `cross_conditions/results/` | Per-target contrasts and pooled estimates for D minus A, E minus A and C minus A, the interaction term, and the per-condition scored counts. |
| `data/revision_2026/` | Every derived table the article and Additional file 1 are rendered from. |
| `environment/` | Conda environment of the analysis and docking machine, as a package list and as an explicit lockfile. |

## Reproducing the analysis tables

Download the data deposit, unpack it so that `targets/` sits beside this file,
then run, with `scripts/revision_2026` on the Python path:

```bash
cd scripts/revision_2026
python revision_populations.py   --data-root ../.. --ablation-root <ablation dir> --out-dir ../../data/revision_2026
python revision_endtoend.py      --data-root ../.. --out-dir ../../data/revision_2026
python revision_decomposition.py --data-root ../.. --ablation-root <ablation dir> --out-dir ../../data/revision_2026
python revision_dataset_audit.py --data-root ../.. --out-dir ../../data/revision_2026
python revision_hitquality.py    --data-root ../.. --out-dir ../../data/revision_2026
python early_recognition.py      --data-root ../.. --out-dir ../../data/revision_2026
python score_audit.py            --data-root ../.. --ablation-root <ablation dir> --out-dir ../../data/revision_2026
python attrition.py              --data-root ../.. --out-dir ../../data/revision_2026
python make_figures_v2.py
```

`--ablation-root` points at a directory holding, per target, a
`<target>/skill_nofilter/scores_skill_nofilter.csv` file; in the data deposit
those files are at `targets/<target>/05_evaluation/skill_nofilter/`.

Verified on 9 October 2026 against the published data deposit
(10.5281/zenodo.23270845). Reproducing from the deposit alone gives
byte-identical `intersection_delong.csv`, `endtoend_delong.csv`,
`decomposition.csv`, `decoy_contamination.csv`, `decoy_sensitivity.csv`,
`assay_provenance.csv`, `score_range_audit.csv`,
`nonphysical_score_sensitivity.csv`, `attrition_by_arm.csv` and all three meta-analysis JSON
files.

With deposit v1.0 three tables did not reproduce, because it carried no
preparation records:

- `prep_matched.csv` and the filter-exclusion columns of `populations.csv` need
  the per-target preparation logs (`03_library_preparation/{naive,skill}/preparation_log.csv`
  and `skill/filter_report.csv`).
- `filter_composition.csv`, `hitlist_quality.csv` and
  `filtered_active_potency.csv` need those logs plus the curated activity
  tables.

Deposit v2.0 (10.5281/zenodo.23270845) carries those preparation logs and the
prepared receptor and ligand PDBQT files, so every table now reproduces from the
deposit alone, with no file from the author's machine.

## Early-recognition metrics

`early_recognition.py` writes the two tables behind Additional file 1, Tables
S12 and S13, with the analytic chance expectation of every metric alongside its
value. BEDROC chance follows Truchon and Bayly (2007) at RIE = 1 and is near
0.06 at alpha = 20 for these libraries, not 0.5. Enrichment-factor values differ
from the October 2026 working tables in the third decimal because the released
script computes them on the same ranked list as the other metrics; no
comparison, median or count in the article changes.

## Attrition by arm

`attrition.py` writes the per-arm scoring rates behind Additional file 1, Table S18: for every
target and protocol, the fraction of actives and the fraction of decoys that received a score,
their difference, and the end-to-end ROC AUC of the same arm. Across the 22 arms the difference
and the AUC correlate at Pearson r = 0.903 (p = 8.5e-09), and all six arms with an AUC below 0.5
have a negative difference. This is what makes the end-to-end endpoint interpretable: it charges a
workflow for the compounds it never scored, so it is driven by which class of compound goes
unscored as much as by the ordering of the rest.

## Score-range audit

`score_audit.py` records the released score ranges per route and repeats the
end-to-end comparison after treating every score below -20 kcal/mol as a
preparation failure. The naive route reaches no value below -20 on any target,
which is why that threshold was chosen.
