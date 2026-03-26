# bace1 Active Compound Curation Report

## Curation Pipeline Summary

| Step | Description | Count |
|------|------------|-------|
| 0 | Raw ChEMBL records | 14518 |
| 1 | Valid SMILES (len <= 200) | 14503 |
| 2 | Standardized (salt removal, neutralize) | 14503 |
| 3 | Property filter (MW 150.0-650.0, HA<=50) | 13714 |
| 4 | Deduplicated (InChIKey connectivity) | 8029 |
| 5 | Butina clustering (Tc cutoff 0.4) | 1220 clusters |
| 6 | MaxMin diversity pick | **100 actives** |

## Final Dataset Properties

| Property | Min | Median | Max | Mean | Std |
|----------|-----|--------|-----|------|-----|
| MW | 151.1 | 442.4 | 639.7 | 434.7 | 108.8 |
| LogP | -0.47 | 4.06 | 10.96 | 4.03 | 1.83 |
| HBD | 0 | 2 | 8 | 2.1 | 1.6 |
| HBA | 1 | 5 | 11 | 5.2 | 2.2 |
| RotBonds | 0 | 5 | 24 | 5.4 | 4.0 |
| TPSA | 23.8 | 80.7 | 229.4 | 85.7 | 37.4 |
| pChEMBL | 5.01 | 6.16 | 9.21 | 6.41 | 1.06 |

## Diversity

- Butina clusters (all deduped): 1220
- Butina clusters (final 100): 100
- Cluster size distribution: {1: 100}

## Checkpoint CP1

- [x] 100 actives (target: 100)
- [x] 100 clusters (target: >= 30)
- [x] All SMILES valid (RDKit-verified)
- [x] Properties span drug-like range
