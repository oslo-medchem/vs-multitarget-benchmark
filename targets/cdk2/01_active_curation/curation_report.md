# cdk2 Active Compound Curation Report

## Curation Pipeline Summary

| Step | Description | Count |
|------|------------|-------|
| 0 | Raw ChEMBL records | 3123 |
| 1 | Valid SMILES (len <= 200) | 3123 |
| 2 | Standardized (salt removal, neutralize) | 3123 |
| 3 | Property filter (MW 150.0-650.0, HA<=50) | 3106 |
| 4 | Deduplicated (InChIKey connectivity) | 2292 |
| 5 | Butina clustering (Tc cutoff 0.4) | 716 clusters |
| 6 | MaxMin diversity pick | **100 actives** |

## Final Dataset Properties

| Property | Min | Median | Max | Mean | Std |
|----------|-----|--------|-----|------|-----|
| MW | 172.2 | 365.0 | 564.5 | 370.5 | 81.5 |
| LogP | -0.83 | 3.17 | 5.49 | 3.17 | 1.31 |
| HBD | 0 | 2 | 9 | 2.1 | 1.3 |
| HBA | 2 | 5 | 12 | 5.2 | 1.9 |
| RotBonds | 0 | 3 | 14 | 3.7 | 2.4 |
| TPSA | 32.6 | 84.9 | 217.6 | 89.7 | 31.9 |
| pChEMBL | 5.00 | 6.00 | 8.47 | 6.20 | 0.86 |

## Diversity

- Butina clusters (all deduped): 716
- Butina clusters (final 100): 100
- Cluster size distribution: {1: 100}

## Checkpoint CP1

- [x] 100 actives (target: 100)
- [x] 100 clusters (target: >= 30)
- [x] All SMILES valid (RDKit-verified)
- [x] Properties span drug-like range
