# dpp4 Active Compound Curation Report

## Curation Pipeline Summary

| Step | Description | Count |
|------|------------|-------|
| 0 | Raw ChEMBL records | 5008 |
| 1 | Valid SMILES (len <= 200) | 5005 |
| 2 | Standardized (salt removal, neutralize) | 5005 |
| 3 | Property filter (MW 150.0-650.0, HA<=50) | 4988 |
| 4 | Deduplicated (InChIKey connectivity) | 3789 |
| 5 | Butina clustering (Tc cutoff 0.4) | 724 clusters |
| 6 | MaxMin diversity pick | **100 actives** |

## Final Dataset Properties

| Property | Min | Median | Max | Mean | Std |
|----------|-----|--------|-----|------|-----|
| MW | 168.2 | 355.1 | 637.7 | 370.3 | 110.8 |
| LogP | -2.35 | 2.43 | 9.56 | 2.32 | 2.01 |
| HBD | 0 | 2 | 5 | 1.8 | 1.2 |
| HBA | 1 | 5 | 12 | 5.0 | 2.0 |
| RotBonds | 0 | 4 | 12 | 4.5 | 2.4 |
| TPSA | 17.1 | 85.2 | 178.1 | 86.6 | 32.8 |
| pChEMBL | 5.10 | 6.98 | 9.89 | 6.98 | 1.18 |

## Diversity

- Butina clusters (all deduped): 724
- Butina clusters (final 100): 100
- Cluster size distribution: {1: 100}

## Checkpoint CP1

- [x] 100 actives (target: 100)
- [x] 100 clusters (target: >= 30)
- [x] All SMILES valid (RDKit-verified)
- [x] Properties span drug-like range
