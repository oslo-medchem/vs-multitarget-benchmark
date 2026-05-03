---
name: virtual-screening-pipeline
description: "Use when designing or executing structure-based virtual screening campaigns. Triggers on virtual screening, docking campaign, hit identification, compound library preparation, or screening pipeline design for any protein target."
allowed-tools: Read Write Edit Bash Agent Glob Grep WebFetch WebSearch
---

# Virtual Screening Pipeline

Guide for designing and executing structure-based virtual screening (SBVS) campaigns, from target preparation through hit triage and experimental validation planning.

## When to Use This Skill

- Designing a virtual screening campaign against a protein target
- Preparing compound libraries for docking
- Setting up docking protocols (Glide, AutoDock Vina, DiffDock)
- Analyzing docking results and enrichment metrics
- Hit triage with ADMET filters
- Any mention of: virtual screening, docking campaign, hit identification, compound library, screening pipeline

## When NOT to Use This Skill

- For lead optimization/FEP -- use fep-spell-workflow
- For diffusion-based docking -- use diffdock

## Performance Evaluation Criteria

### Process Adherence
- Correct pipeline stages followed in proper order (target prep through validation planning)
- Hierarchical docking protocol applied (HTVS/Vina, SP, XP, MM-GBSA)
- Benchmarking performed with known actives before running the full screen

### Output Quality
- Enrichment metrics calculated and meeting quality thresholds (EF1% > 10, ROC AUC > 0.7)
- Hit triage cascade applied with appropriate filters at each stage
- Final hit list contains diverse, purchasable compounds with novelty assessment

### Completeness
- All seven pipeline stages addressed (target prep through validation planning)
- ChEMBL integration used for novelty checking of top hits
- Troubleshooting guidance for poor enrichment, tight score clustering, and ADMET filter attrition

## Pipeline Overview

```
Target Prep → Binding Site → Library Prep → Docking → Post-processing → Hit Triage → Validation
```

## Step 1: Target Preparation

### Structure Selection
- Prefer high-resolution crystal structures (< 2.5 Å) with co-crystallized ligand
- Check for missing loops, alternate conformations, crystallographic artifacts
- For GPCRs: check active vs inactive state, consider ensemble docking

### Protein Preparation
```
Schrodinger Protein Preparation Wizard:
  - Add hydrogens at pH 7.0 ± 2.0
  - Assign protonation states (PROPKA)
  - Optimize H-bond network
  - Restrained minimization (RMSD < 0.3 Å from crystal)

Open-source alternative:
  - PDB2PQR + PROPKA for protonation
  - OpenBabel for hydrogen addition
  - Reduce (Richardson lab) for H-bond optimization
```

### For Homology Models
- Template selection: > 30% sequence identity, same conformational state
- Model with MODELLER, AlphaFold2, or RoseTTAFold
- Validate: Ramachandran plot, DOPE score, ProSA
- Use with caution for VS: binding site RMSD to template < 1.5 Å

## Step 2: Binding Site Definition

### Known Binding Site
- Define grid from co-crystallized ligand (center of mass + 10-15 Å box)
- Schrodinger Glide: Receptor Grid Generation with default van der Waals scaling

### Unknown Binding Site
- SiteMap (Schrodinger), fpocket, DoGSiteScorer
- Conservation analysis (ConSurf) to prioritize sites
- For GPCRs: orthosteric site in TM bundle, allosteric sites at intracellular face or lipid interface

## Step 3: Compound Library Preparation

### Library Sources
| Source | Size | Use Case |
|--------|------|----------|
| ZINC20/22 | 230M+ | General screening, purchasable |
| Enamine REAL | 6B+ | Make-on-demand, diverse |
| ChEMBL actives | Variable | Known actives for benchmarking |
| In-house | Variable | Custom focused libraries |

### Preparation Protocol
```
1. Remove salts, neutralize
2. Generate 3D conformers (OMEGA, RDKit ETKDG)
3. Assign protonation states at pH 7.4 (Epik, Dimorphite-DL)
4. Generate tautomers (retain top 3-5)
5. Filter: MW 150-600, ALogP -2 to 5, rotatable bonds ≤ 10
6. Remove PAINS, Brenk, aggregators
7. Deduplicate by InChIKey
```

## Step 4: Docking

### Docking Programs

| Program | Speed | Accuracy | Best For |
|---------|-------|----------|----------|
| **Glide SP** | Medium | Good | Primary screen (10K-1M) |
| **Glide XP** | Slow | Better | Rescoring top hits (1K-10K) |
| **AutoDock Vina** | Fast | Good | Large libraries, open-source |
| **DiffDock** | Slow | Varies | Blind docking, flexible targets |

### Hierarchical Protocol (Recommended)
1. **Glide HTVS** or **Vina**: Full library → top 10%
2. **Glide SP**: Top 10% → top 1%
3. **Glide XP**: Top 1% → final candidates
4. Optional: MM-GBSA rescoring on top 100-500

## Step 5: Post-Processing & Enrichment Analysis

### Enrichment Metrics
| Metric | Good Value | Description |
|--------|-----------|-------------|
| EF1% | > 10 | Enrichment factor at 1% of ranked database |
| ROC AUC | > 0.7 | Area under ROC curve |
| BEDROC (α=20) | > 0.5 | Early enrichment metric |
| Hit rate | > 5% | Active fraction in tested compounds |

### Benchmarking (Before Live Screen)
- Use DUD-E or LIT-PCBA decoys for target
- Calculate enrichment metrics
- Optimize protocol until satisfactory enrichment

## Step 6: Hit Triage

### Multi-Filter Cascade
1. **Docking score cutoff**: Top N or score < threshold
2. **Visual inspection**: Correct binding mode, key interactions preserved
3. **Interaction fingerprint clustering**: Diverse binding modes
4. **ADMET filters**: Lipinski Ro5, Veber rules, QED > 0.4
5. **PAINS/nuisance filter**: Remove frequent hitters
6. **Chemical diversity**: Cluster by Tanimoto (ECFP4), pick representatives
7. **Synthetic accessibility**: SA score < 4, purchasability check
8. **ChEMBL/PubChem check**: Known activity? Novelty assessment

### Final Selection
- Aim for 20-100 diverse compounds for experimental testing
- Balance: best scores, diverse scaffolds, good ADMET, purchasable

## Step 7: Validation Planning

- **Binding assays**: SPR, ITC, or fluorescence polarization
- **Functional assays**: Target-specific (e.g., cAMP for GPCRs, enzyme activity)
- **Counter-screen**: Related targets for selectivity
- **Dose-response**: IC50/EC50 determination
- **Hit confirmation**: Repurchase from different vendor, test purity

## Examples

### Example 1: Design a full screening campaign
User says: "Design a virtual screening campaign against GPR37 for novel agonists"
Actions:
1. Retrieve GPR37 structure (homology model or AlphaFold), prepare with Protein Prep Wizard
2. Define binding site from known ligand poses or SiteMap
3. Prepare ZINC20 drug-like subset (~10M compounds) with LigPrep
4. Run hierarchical docking: HTVS → SP → XP → MM-GBSA
5. Apply hit triage: ADMET filters, PAINS removal, diversity clustering
6. Check ChEMBL for novelty of top hits
Result: Ranked hit list of 50-100 diverse, drug-like compounds with novelty assessment.

### Example 2: Benchmark docking protocol
User says: "Benchmark our Glide protocol against known NPY5R actives before running the full screen"
Actions:
1. Retrieve known actives from ChEMBL for NPY5R
2. Generate DUD-E decoys or use LIT-PCBA set
3. Run docking on actives + decoys, calculate EF1%, ROC AUC, BEDROC
4. Iterate protocol (grid size, constraints) until EF1% > 10
Result: Validated docking protocol with enrichment metrics.

### Example 3: Triage docking results
User says: "I have 5000 Glide SP hits, help me triage to 50 compounds for testing"
Actions:
1. Filter by docking score threshold (top 500)
2. Visual inspection of binding modes for key interactions
3. Apply ADMET filters (Lipinski, QED > 0.4), remove PAINS
4. Cluster by Tanimoto similarity (ECFP4), pick diverse representatives
5. Check purchasability and synthetic accessibility
Result: Final list of 50 diverse, purchasable compounds ranked by score and diversity.

## Troubleshooting

| Issue | Cause | Fix |
|-------|-------|-----|
| Poor enrichment in benchmarking | Grid too large or wrong site | Refine grid around co-crystallized ligand, add pharmacophore constraints |
| Too few hits pass ADMET filters | Library not pre-filtered | Apply drug-likeness filters before docking, not after |
| Docking scores cluster tightly | Scoring function limitation | Use MM-GBSA rescoring or interaction fingerprints to differentiate |
| Known actives score poorly | Receptor conformation mismatch | Try ensemble docking with multiple receptor conformations |
