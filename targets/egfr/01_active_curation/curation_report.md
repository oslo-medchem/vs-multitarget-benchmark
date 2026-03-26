# egfr Active Compound Curation Report

## Curation Pipeline Summary

| Step | Description | Count |
|------|------------|-------|
| 0 | Raw ChEMBL records | 19309 |
| 1 | Valid SMILES (len <= 200) | 19306 |
| 2 | Standardized (salt removal, neutralize) | 19306 |
| 3 | Property filter (MW 150.0-650.0, HA<=50) | 18198 |
| 4 | Deduplicated (InChIKey connectivity) | 9301 |
| 5 | Butina clustering (Tc cutoff 0.4) | 1739 clusters |
| 6 | MaxMin diversity pick | **100 actives** |

## Final Dataset Properties

| Property | Min | Median | Max | Mean | Std |
|----------|-----|--------|-----|------|-----|
| MW | 159.2 | 411.0 | 636.7 | 418.0 | 102.8 |
| LogP | -1.04 | 3.92 | 8.95 | 4.01 | 1.77 |
| HBD | 0 | 2 | 8 | 1.9 | 1.4 |
| HBA | 0 | 5 | 12 | 5.5 | 2.3 |
| RotBonds | 0 | 4 | 18 | 5.1 | 3.5 |
| TPSA | 3.9 | 84.4 | 163.2 | 86.8 | 35.1 |
| pChEMBL | 5.00 | 6.17 | 10.00 | 6.44 | 1.16 |

## Diversity

- Butina clusters (all deduped): 1739
- Butina clusters (final 100): 100
- Cluster size distribution: {1: 100}

## Checkpoint CP1

- [x] 100 actives (target: 100)
- [x] 100 clusters (target: >= 30)
- [x] All SMILES valid (RDKit-verified)
- [x] Properties span drug-like range
