# Data deposit

The per-target benchmark data (`targets/` directory: ChEMBL-curated actives,
property-matched decoy sets, Uni-Dock score tables for both protocols and for
the filter-ablation condition, the original run scripts, box centres and
per-target evaluation outputs across 11 protein targets) are archived as a
separate Zenodo Dataset to keep this code repository lean and under Zenodo's
GitHub auto-archive size threshold.

Deposit v2.0 (10.5281/zenodo.23270845) adds the per-target preparation logs,
the substructure-filter reports and a second archive with the prepared receptor
and ligand PDBQT files for both routes, so every analysis table regenerates from
it and the docking itself can be re-run. It supersedes 10.5281/zenodo.20036964,
which held the score tables only.

This split follows the GigaScience reviewer-preferred pattern:

- **Code repository (this repo, MIT):** scripts, prompts, skill file,
  meta-analysis outputs, figures.
  - GitHub: https://github.com/oslo-medchem/vs-multitarget-benchmark
  - Zenodo concept DOI: [10.5281/zenodo.20036546](https://doi.org/10.5281/zenodo.20036546)
  - Per-version DOI for v1.1-gigascience: [10.5281/zenodo.20036620](https://doi.org/10.5281/zenodo.20036620)

- **Data deposit (Zenodo Dataset, CC-BY 4.0):** ChEMBL actives, decoys,
  Uni-Dock score tables, per-target evaluation outputs.
  - Zenodo DOI: [10.5281/zenodo.23270845](https://doi.org/10.5281/zenodo.23270845)
  - Concept DOI: [10.5281/zenodo.23270844](https://doi.org/10.5281/zenodo.23270844)
  - Cross-linked to the software archive via Zenodo's `isSupplementTo` relation.

To re-execute the analysis end to end, clone this repository, download the
data zip from the Zenodo Dataset, and unzip it into the repository root so
`targets/` is restored alongside `common/`, `prompts/`, `06_meta_analysis/`,
and `07_manuscript/`.
