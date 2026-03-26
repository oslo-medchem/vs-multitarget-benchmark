# esr1 Active Compound Curation Report

## Curation Pipeline Summary

| Step | Description | Count |
|------|------------|-------|
| 0 | Raw ChEMBL records | 7255 |
| 1 | Valid SMILES (len <= 200) | 7215 |
| 2 | Standardized (salt removal, neutralize) | 7215 |
| 3 | Property filter (MW 150.0-650.0, HA<=50) | 6974 |
| 4 | Deduplicated (InChIKey connectivity) | 3544 |
| 5 | Butina clustering (Tc cutoff 0.4) | 798 clusters |
| 6 | MaxMin diversity pick | **100 actives** |

## Final Dataset Properties

| Property | Min | Median | Max | Mean | Std |
|----------|-----|--------|-----|------|-----|
| MW | 240.3 | 393.6 | 647.6 | 392.6 | 90.0 |
| LogP | 0.63 | 4.18 | 8.70 | 4.43 | 1.54 |
| HBD | 0 | 1 | 7 | 1.5 | 1.2 |
| HBA | 0 | 4 | 10 | 4.6 | 2.0 |
| RotBonds | 0 | 4 | 15 | 4.5 | 2.9 |
| TPSA | 9.2 | 68.2 | 239.3 | 71.5 | 33.0 |
| pChEMBL | 5.00 | 5.61 | 10.00 | 6.21 | 1.27 |

## Diversity

- Butina clusters (all deduped): 798
- Butina clusters (final 100): 100
- Cluster size distribution: {1: 100}

## Checkpoint CP1

- [x] 100 actives (target: 100)
- [x] 100 clusters (target: >= 30)
- [x] All SMILES valid (RDKit-verified)
- [x] Properties span drug-like range
