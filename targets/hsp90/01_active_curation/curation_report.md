# hsp90 Active Compound Curation Report

## Curation Pipeline Summary

| Step | Description | Count |
|------|------------|-------|
| 0 | Raw ChEMBL records | 2051 |
| 1 | Valid SMILES (len <= 200) | 2049 |
| 2 | Standardized (salt removal, neutralize) | 2049 |
| 3 | Property filter (MW 150.0-650.0, HA<=50) | 1980 |
| 4 | Deduplicated (InChIKey connectivity) | 1478 |
| 5 | Butina clustering (Tc cutoff 0.4) | 330 clusters |
| 6 | MaxMin diversity pick | **100 actives** |

## Final Dataset Properties

| Property | Min | Median | Max | Mean | Std |
|----------|-----|--------|-----|------|-----|
| MW | 170.2 | 378.5 | 620.7 | 380.5 | 96.4 |
| LogP | -1.62 | 3.57 | 7.45 | 3.47 | 1.51 |
| HBD | 0 | 1 | 6 | 1.6 | 1.2 |
| HBA | 1 | 5 | 12 | 5.2 | 2.0 |
| RotBonds | 0 | 4 | 9 | 4.2 | 2.3 |
| TPSA | 12.5 | 78.5 | 187.9 | 81.4 | 28.9 |
| pChEMBL | 5.01 | 6.06 | 10.10 | 6.40 | 1.14 |

## Diversity

- Butina clusters (all deduped): 330
- Butina clusters (final 100): 100
- Cluster size distribution: {1: 100}

## Checkpoint CP1

- [x] 100 actives (target: 100)
- [x] 100 clusters (target: >= 30)
- [x] All SMILES valid (RDKit-verified)
- [x] Properties span drug-like range
