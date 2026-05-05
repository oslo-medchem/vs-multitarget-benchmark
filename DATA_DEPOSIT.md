# Data deposit

The full per-target benchmark data (`targets/` directory: ChEMBL-curated
actives, property-matched decoys, receptor and ligand PDBQT files, Uni-Dock
score tables, and per-target evaluation outputs across 11 protein targets) is
archived as a separate Zenodo Dataset to keep this code repository lean and
under Zenodo's GitHub auto-archive size threshold.

This split follows the GigaScience reviewer-preferred pattern:

- **Code repository (this repo, MIT):** scripts, prompts, skill file,
  meta-analysis outputs, figures.
  - GitHub: https://github.com/oslo-medchem/vs-multitarget-benchmark
  - Zenodo concept DOI: [10.5281/zenodo.20036546](https://doi.org/10.5281/zenodo.20036546)
  - Per-version DOI for v1.1-gigascience: [10.5281/zenodo.20036620](https://doi.org/10.5281/zenodo.20036620)

- **Data deposit (Zenodo Dataset, CC-BY 4.0):** ChEMBL actives, decoys,
  receptors, Uni-Dock score tables, per-target evaluation outputs.
  - Zenodo DOI: [10.5281/zenodo.20036964](https://doi.org/10.5281/zenodo.20036964)
  - Concept DOI: [10.5281/zenodo.20036963](https://doi.org/10.5281/zenodo.20036963)
  - Cross-linked to the software archive via Zenodo's `isSupplementTo` relation.

To re-execute the analysis end to end, clone this repository, download the
data zip from the Zenodo Dataset, and unzip it into the repository root so
`targets/` is restored alongside `common/`, `prompts/`, `06_meta_analysis/`,
and `07_manuscript/`.
