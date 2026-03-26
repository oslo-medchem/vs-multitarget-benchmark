# thrombin Active Compound Curation Report

## Curation Pipeline Summary

| Step | Description | Count |
|------|------------|-------|
| 0 | Raw ChEMBL records | 5505 |
| 1 | Valid SMILES (len <= 200) | 5451 |
| 2 | Standardized (salt removal, neutralize) | 5451 |
| 3 | Property filter (MW 150.0-650.0, HA<=50) | 5071 |
| 4 | Deduplicated (InChIKey connectivity) | 4154 |
| 5 | Butina clustering (Tc cutoff 0.4) | 848 clusters |
| 6 | MaxMin diversity pick | **100 actives** |

## Final Dataset Properties

| Property | Min | Median | Max | Mean | Std |
|----------|-----|--------|-----|------|-----|
| MW | 158.2 | 449.0 | 641.8 | 437.7 | 122.4 |
| LogP | -3.16 | 3.21 | 8.02 | 2.89 | 2.03 |
| HBD | 0 | 2 | 8 | 2.4 | 1.7 |
| HBA | 1 | 6 | 12 | 5.7 | 2.1 |
| RotBonds | 0 | 6 | 16 | 6.3 | 3.5 |
| TPSA | 20.2 | 103.7 | 272.7 | 109.1 | 45.3 |
| pChEMBL | 5.02 | 6.08 | 9.42 | 6.44 | 1.09 |

## Diversity

- Butina clusters (all deduped): 848
- Butina clusters (final 100): 100
- Cluster size distribution: {1: 100}

## Checkpoint CP1

- [x] 100 actives (target: 100)
- [x] 100 clusters (target: >= 30)
- [x] All SMILES valid (RDKit-verified)
- [x] Properties span drug-like range
