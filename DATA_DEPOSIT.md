# Data deposit

The full per-target benchmark data (`targets/` directory: ChEMBL-curated
actives, property-matched decoys, receptor and ligand PDBQT files, Uni-Dock
score tables, and per-target evaluation outputs across 11 protein targets) is
archived as a separate Zenodo Dataset to keep this code repository lean and
under Zenodo's GitHub auto-archive size threshold.

**Data deposit (Zenodo):** DOI to be added in v1.1.1 once minted.

To re-execute the analysis end to end, clone this repository, download the
data deposit, and unzip it into the repository root so `targets/` is restored
alongside `common/`, `prompts/`, `06_meta_analysis/`, and `07_manuscript/`.

This split follows the GigaScience reviewer-preferred pattern: lean code +
small example outputs on GitHub (with a software concept DOI), bulky
benchmark data on a separate Zenodo Dataset deposit (with a data DOI). The
manuscript Data Availability and Code Availability statements cite both DOIs.
