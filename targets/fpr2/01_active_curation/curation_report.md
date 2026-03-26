# FPR2 Active Compound Curation Report

## Curation Pipeline Summary

| Step | Description | Count |
|------|------------|-------|
| 0 | Raw ChEMBL records | 2847 |
| 1 | Valid SMILES (len <= 200) | 2828 |
| 2 | Standardized (salt removal, neutralize) | 2828 |
| 3 | Property filter (MW 150-650, HA<=50) | 2736 |
| 4 | Deduplicated (InChIKey connectivity) | 1901 |
| 5 | Butina clustering (Tc cutoff 0.4) | 288 clusters |
| 6 | MaxMin diversity pick | **100 actives** |

## Final Dataset Properties

| Property | Min | Median | Max | Mean | Std |
|----------|-----|--------|-----|------|-----|
| MW | 288.3 | 428.4 | 633.7 | 436.3 | 73.6 |
| LogP | -3.98 | 3.55 | 6.00 | 3.36 | 1.53 |
| HBD | 0 | 2 | 11 | 2.4 | 1.5 |
| HBA | 2 | 5 | 14 | 4.8 | 2.0 |
| RotBonds | 2 | 6 | 22 | 6.5 | 2.9 |
| TPSA | 43.8 | 88.7 | 295.5 | 96.5 | 40.2 |
| pChEMBL | 5.00 | 6.80 | 9.52 | 7.04 | 1.25 |

## Diversity

- Butina clusters (all deduped): 288
- Butina clusters (final 100): 100
- Cluster size distribution: {1: 100}

## Checkpoint CP1

- [x] 100 actives (target: 100)
- [x] 100 clusters (target: >= 30)
- [x] All SMILES valid (RDKit-verified)
- [x] Properties span drug-like range
