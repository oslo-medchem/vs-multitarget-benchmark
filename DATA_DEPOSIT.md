# Data deposit

The per-target benchmark data (`targets/` directory: ChEMBL-curated actives,
property-matched decoy sets, Uni-Dock score tables for both protocols and for
the filter-ablation condition, the original run scripts, box centres and
per-target evaluation outputs across 11 protein targets) are archived as a
separate Zenodo Dataset to keep this code repository lean and under Zenodo's
GitHub auto-archive size threshold.

The deposit released as v1.0 does **not** contain the prepared receptor and
ligand PDBQT files, the per-target preparation logs or the decoy-matching
intermediates, so docking and decoy selection cannot be re-run from it as
released. `REVISION_2026.md` lists exactly which analysis tables do and do not
reproduce from it. Those preparation inputs are being added in a new version of
the deposit.

This split follows the GigaScience reviewer-preferred pattern:

- **Code repository (this repo, MIT):** scripts, prompts, skill file,
  meta-analysis outputs, figures.
  - GitHub: https://github.com/oslo-medchem/vs-multitarget-benchmark
  - Zenodo concept DOI: [10.5281/zenodo.20036546](https://doi.org/10.5281/zenodo.20036546)
  - Per-version DOI for v1.1-gigascience: [10.5281/zenodo.20036620](https://doi.org/10.5281/zenodo.20036620)

- **Data deposit (Zenodo Dataset, CC-BY 4.0):** ChEMBL actives, decoys,
  Uni-Dock score tables, per-target evaluation outputs.
  - Zenodo DOI: [10.5281/zenodo.20036964](https://doi.org/10.5281/zenodo.20036964)
  - Concept DOI: [10.5281/zenodo.20036963](https://doi.org/10.5281/zenodo.20036963)
  - Cross-linked to the software archive via Zenodo's `isSupplementTo` relation.

To re-execute the analysis end to end, clone this repository, download the
data zip from the Zenodo Dataset, and unzip it into the repository root so
`targets/` is restored alongside `common/`, `prompts/`, `06_meta_analysis/`,
and `07_manuscript/`.
