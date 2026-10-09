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
python make_figures_v2.py
```

`--ablation-root` points at a directory holding, per target, a
`<target>/skill_nofilter/scores_skill_nofilter.csv` file; in the data deposit
those files are at `targets/<target>/05_evaluation/skill_nofilter/`.

Verified on 9 October 2026 against the published data deposit
(10.5281/zenodo.20036964). Reproducing from the deposit alone gives
byte-identical `intersection_delong.csv`, `endtoend_delong.csv`,
`decomposition.csv`, `decoy_contamination.csv`, `decoy_sensitivity.csv`,
`assay_provenance.csv`, `score_range_audit.csv`,
`nonphysical_score_sensitivity.csv` and all three meta-analysis JSON files.

Three tables need inputs the deposit does not carry and therefore do not
reproduce from it alone:

- `prep_matched.csv` and the filter-exclusion columns of `populations.csv` need
  the per-target preparation logs (`03_library_preparation/{naive,skill}/preparation_log.csv`
  and `skill/filter_report.csv`).
- `filter_composition.csv`, `hitlist_quality.csv` and
  `filtered_active_potency.csv` need those logs plus the curated activity
  tables.

The preparation logs and the prepared receptor and ligand PDBQT files are added
to the data deposit in the version cited by the article, so that every table
reproduces from the deposit without any file from the author's machine.

## Early-recognition metrics

`early_recognition.py` writes the two tables behind Additional file 1, Tables
S12 and S13, with the analytic chance expectation of every metric alongside its
value. BEDROC chance follows Truchon and Bayly (2007) at RIE = 1 and is near
0.06 at alpha = 20 for these libraries, not 0.5. Enrichment-factor values differ
from the October 2026 working tables in the third decimal because the released
script computes them on the same ranked list as the other metrics; no
comparison, median or count in the article changes.

## Score-range audit

`score_audit.py` records the released score ranges per route and repeats the
end-to-end comparison after treating every score below -20 kcal/mol as a
preparation failure. The naive route reaches no value below -20 on any target,
which is why that threshold was chosen.
