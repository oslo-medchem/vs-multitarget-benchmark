# cox2 Active Compound Curation Report

## Curation Pipeline Summary

| Step | Description | Count |
|------|------------|-------|
| 0 | Raw ChEMBL records | 5319 |
| 1 | Valid SMILES (len <= 200) | 5307 |
| 2 | Standardized (salt removal, neutralize) | 5307 |
| 3 | Property filter (MW 150.0-650.0, HA<=50) | 5267 |
| 4 | Deduplicated (InChIKey connectivity) | 3549 |
| 5 | Butina clustering (Tc cutoff 0.4) | 896 clusters |
| 6 | MaxMin diversity pick | **100 actives** |

## Final Dataset Properties

| Property | Min | Median | Max | Mean | Std |
|----------|-----|--------|-----|------|-----|
| MW | 162.1 | 363.4 | 637.7 | 386.3 | 106.8 |
| LogP | -2.12 | 3.67 | 7.66 | 3.66 | 1.83 |
| HBD | 0 | 1 | 8 | 1.5 | 1.5 |
| HBA | 0 | 5 | 16 | 5.2 | 2.7 |
| RotBonds | 0 | 4 | 17 | 4.7 | 3.0 |
| TPSA | 3.2 | 75.8 | 225.5 | 78.9 | 41.5 |
| pChEMBL | 5.00 | 6.19 | 10.10 | 6.34 | 0.95 |

## Diversity

- Butina clusters (all deduped): 896
- Butina clusters (final 100): 100
- Cluster size distribution: {1: 100}

## Checkpoint CP1

- [x] 100 actives (target: 100)
- [x] 100 clusters (target: >= 30)
- [x] All SMILES valid (RDKit-verified)
- [x] Properties span drug-like range
