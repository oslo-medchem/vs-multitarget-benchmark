# p38 Active Compound Curation Report

## Curation Pipeline Summary

| Step | Description | Count |
|------|------------|-------|
| 0 | Raw ChEMBL records | 6767 |
| 1 | Valid SMILES (len <= 200) | 6761 |
| 2 | Standardized (salt removal, neutralize) | 6761 |
| 3 | Property filter (MW 150.0-650.0, HA<=50) | 6171 |
| 4 | Deduplicated (InChIKey connectivity) | 4405 |
| 5 | Butina clustering (Tc cutoff 0.4) | 851 clusters |
| 6 | MaxMin diversity pick | **100 actives** |

## Final Dataset Properties

| Property | Min | Median | Max | Mean | Std |
|----------|-----|--------|-----|------|-----|
| MW | 206.6 | 414.3 | 637.9 | 426.5 | 102.3 |
| LogP | -1.08 | 4.44 | 8.95 | 4.34 | 1.53 |
| HBD | 0 | 2 | 8 | 1.8 | 1.4 |
| HBA | 1 | 5 | 11 | 5.0 | 2.0 |
| RotBonds | 0 | 4 | 20 | 4.9 | 3.1 |
| TPSA | 31.2 | 81.1 | 197.4 | 83.3 | 31.9 |
| pChEMBL | 5.06 | 6.20 | 9.35 | 6.47 | 1.02 |

## Diversity

- Butina clusters (all deduped): 851
- Butina clusters (final 100): 100
- Cluster size distribution: {1: 100}

## Checkpoint CP1

- [x] 100 actives (target: 100)
- [x] 100 clusters (target: >= 30)
- [x] All SMILES valid (RDKit-verified)
- [x] Properties span drug-like range
