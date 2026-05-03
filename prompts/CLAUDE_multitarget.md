# CLAUDE.md

This file provides guidance to Claude Code when working with this repository.

## Purpose

Multi-target virtual screening benchmark comparing **naive defaults** vs **Claude Code skill-guided protocols** across 11 diverse drug targets. Expansion of the single-target FPR2 study in response to journal reviewer feedback requesting generalization evidence.

## Environment

```bash
eval "$(conda shell.bash hook)" && conda activate vina_dock
# Required: vina, meeko, rdkit, openbabel, unidock, scipy, numpy,
#           seaborn, scikit-learn, chembl_webresource_client, matplotlib, pandas
```

Hardware: 2x NVIDIA RTX 4500 Ada (24 GB VRAM), CUDA 12.9, Uni-Dock 1.1.3

## Targets (11 total)

| Name | ChEMBL | PDB | Class |
|------|--------|-----|-------|
| fpr2 | CHEMBL4227 | 7T6S:R | GPCR |
| egfr | CHEMBL203 | 1M17:A | Kinase |
| cdk2 | CHEMBL301 | 1H1Q:A | Kinase |
| cox2 | CHEMBL230 | 3LN1:A | Enzyme |
| esr1 | CHEMBL206 | 1SJ0:A | Nuclear receptor |
| dpp4 | CHEMBL284 | 2RGU:A | Serine protease |
| ache | CHEMBL220 | 4EY7:A | Hydrolase |
| bace1 | CHEMBL4822 | 4IVS:A | Aspartyl protease |
| hsp90 | CHEMBL3880 | 2WI7:A | Chaperone |
| p38 | CHEMBL260 | 1KV2:A | Kinase |
| thrombin | CHEMBL204 | 1HAH:B | Serine protease |

## Pipeline

All scripts are parameterized in `common/` and driven by `target_config.json`.

### Single target
```bash
# Run full pipeline for one target
python run_target_pipeline.py --target egfr --stages 1,2,3,4,5 --gpu 0

# Run specific stages
python run_target_pipeline.py --target egfr --stages 1,2 --skip-existing
```

### All targets
```bash
# Full pipeline (auto-parallelizes CPU stages, sequential GPU)
bash run_all_targets.sh --stages 1,2,3,4,5 --max-parallel 4
```

### Stage breakdown per target
```
Stage 1: Active curation     → targets/{name}/01_active_curation/actives_100.smi
Stage 2: Decoy generation    → targets/{name}/02_decoy_generation/decoys_local_5000.smi
Stage 3: Library preparation → targets/{name}/03_library_preparation/library_labels.csv
         Receptor preparation → targets/{name}/04_docking/{naive,skill}/receptor_*.pdbqt
         PDBQT generation    → targets/{name}/03_library_preparation/{naive,skill}/pdbqt/
Stage 4: GPU docking         → targets/{name}/04_docking/{naive,skill}/scores_*.csv
Stage 5: Evaluation          → targets/{name}/05_evaluation/results/*.json, figures/*.pdf
```

### Cross-target meta-analysis
```bash
cd 06_meta_analysis
python aggregate_metrics.py    # → aggregate_metrics.csv, delong_deltas.csv
python forest_plot.py          # → meta_figures/forest_plot.pdf
python paired_tests.py         # → paired_test_results.json, metrics_heatmap.pdf
```

## Directory Structure

```
target_config.json           # Master config (all targets + global params)
run_target_pipeline.py       # Single-target orchestrator
run_all_targets.sh           # Master scheduler
common/                      # 12 parameterized pipeline scripts
shared/
  chembl_pool_150k.smi       # Shared decoy pool (download once)
targets/
  {name}/                    # Per-target data (01-05 subdirs)
06_meta_analysis/            # Cross-target statistics + figures
07_manuscript/               # Revised manuscript
```

## Key Notes

- `shared/chembl_pool_150k.smi` is downloaded once and reused for all targets
- `target_config.json` is the single source of truth for target parameters
- Box centers are auto-computed from co-crystallized ligand HETATM coordinates and saved to `targets/{name}/box_center.json`
- FPR2 data can be copied from `../vs_benchmark_SUMISSION/` (status="complete" in config)
- Each target's `library_labels.csv` is the ground truth for evaluation
- Skill protocol removes PAINS/Brenk compounds; evaluation handles size mismatch

## Protocol Differences (same as single-target study)

| Parameter | Naive | Skill-Guided |
|-----------|-------|-------------|
| Receptor prep | OpenBabel | Meeko |
| Ligand 3D | OpenBabel gen3d | RDKit ETKDGv3 |
| Filters | None | PAINS + Brenk |
| Box size | 20 A | 25 A |
| Exhaustiveness | 8 | 32 |

## Meeko/Vina Gotchas

- `mk_prepare_receptor.py` uses `--read_pdb` (not `--pdb`), writes PDBQT via `-p` flag
- `mk_prepare_ligand.py` requires explicit hydrogens — add with RDKit `Chem.AddHs(mol, addCoords=True)`
- Conda `vina` build does not support `--log`; use `2>&1 | tee` instead
- Meeko 0.7.1 requires `gemmi` (installed via pip)
