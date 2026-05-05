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
  - Per-version DOI for v1.1-gigascience: [10.5281/zenodo.20036547](https://doi.org/10.5281/zenodo.20036547)

- **Data deposit (Zenodo Dataset, CC-BY 4.0):** ChEMBL actives, decoys,
  receptors, Uni-Dock score tables, per-target evaluation outputs.
  - DOI: pending (manual upload at https://zenodo.org/uploads/new); will be
    added here once minted.

To re-execute the analysis end to end, clone this repository, download the
data deposit, and unzip it into the repository root so `targets/` is restored
alongside `common/`, `prompts/`, `06_meta_analysis/`, and `07_manuscript/`.
