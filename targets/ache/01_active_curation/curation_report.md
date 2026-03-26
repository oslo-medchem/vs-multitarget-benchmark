# ache Active Compound Curation Report

## Curation Pipeline Summary

| Step | Description | Count |
|------|------------|-------|
| 0 | Raw ChEMBL records | 6448 |
| 1 | Valid SMILES (len <= 200) | 6447 |
| 2 | Standardized (salt removal, neutralize) | 6447 |
| 3 | Property filter (MW 150.0-650.0, HA<=50) | 6153 |
| 4 | Deduplicated (InChIKey connectivity) | 4344 |
| 5 | Butina clustering (Tc cutoff 0.4) | 1056 clusters |
| 6 | MaxMin diversity pick | **100 actives** |

## Final Dataset Properties

| Property | Min | Median | Max | Mean | Std |
|----------|-----|--------|-----|------|-----|
| MW | 152.1 | 370.4 | 644.9 | 368.2 | 114.5 |
| LogP | -1.06 | 3.79 | 8.26 | 3.81 | 1.83 |
| HBD | 0 | 1 | 7 | 1.2 | 1.4 |
| HBA | 0 | 4 | 12 | 4.1 | 2.4 |
| RotBonds | 0 | 4 | 17 | 4.4 | 3.7 |
| TPSA | 0.0 | 55.8 | 211.3 | 59.7 | 40.3 |
| pChEMBL | 5.00 | 5.98 | 10.85 | 6.29 | 1.16 |

## Diversity

- Butina clusters (all deduped): 1056
- Butina clusters (final 100): 100
- Cluster size distribution: {1: 100}

## Checkpoint CP1

- [x] 100 actives (target: 100)
- [x] 100 clusters (target: >= 30)
- [x] All SMILES valid (RDKit-verified)
- [x] Properties span drug-like range
