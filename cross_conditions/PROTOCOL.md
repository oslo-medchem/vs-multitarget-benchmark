# Cross-condition ablation: protocol

Purpose: answer referee 1's request to separate preparation-tool effects from docking-parameter effects. The published design has three executed conditions and cannot do this. Two further conditions complete the factorial.

Tags: [V] taken from the manuscript or a verified file; [J] design judgement made here; [?] must be checked on the node before the full run.

## 1. The design

[V] Executed conditions in the published work:

| Condition | Preparation | Docking parameters | Substructure filter |
|---|---|---|---|
| A naive | OpenBabel gen3d | 20 A box, exhaustiveness 8 | none |
| B skill | RDKit ETKDGv3 + Meeko | 25 A box, exhaustiveness 32 | PAINS + Brenk |
| C skill, filter ablated | RDKit ETKDGv3 + Meeko | 25 A box, exhaustiveness 32 | none |

A versus C isolates preparation and parameters jointly. B versus C isolates the filter. Preparation and parameters are not separated from each other.

[J] New conditions, both with the filter off so that all four no-filter conditions act on the same molecule population:

| Condition | Preparation | Docking parameters | Substructure filter |
|---|---|---|---|
| D | RDKit ETKDGv3 + Meeko (skill) | 20 A box, exhaustiveness 8 (naive) | none |
| E | OpenBabel gen3d (naive) | 25 A box, exhaustiveness 32 (skill) | none |

This gives a complete 2 x 2 factorial over preparation and parameters on a fixed population:

```
                     naive parameters      skill parameters
naive preparation    A                     E
skill preparation    D                     C
```

Contrasts available after the run:

- **Preparation alone, at naive parameters**: D - A
- **Preparation alone, at skill parameters**: C - E
- **Parameters alone, at naive preparation**: E - A
- **Parameters alone, at skill preparation**: C - D
- **Interaction**: (C - D) - (E - A), equivalently (C - E) - (D - A)
- **Joint effect, already published**: C - A

[J] Reporting both marginal contrasts rather than one is the point. If preparation alone is null at both parameter settings and parameters alone is null at both preparation settings, the published joint null is explained. If either marginal is non-null while the joint is null, the two components cancel, which is a finding in itself and is exactly what the referee suspected might be hiding in the joint comparison.

## 2. What is reused and what is new

[V] The manuscript states that `scripts/run_stepwise_ablation.py` reuses the existing PDBQT files so that only the docking step is new.

- **Reused, not regenerated**: receptor PDBQT per target per preparation route; ligand PDBQT per compound per preparation route; the active and decoy identifier lists; the box centre.
- **New**: the docking runs only.

[J] Reuse is what makes this cheap and also what makes it a clean experiment: condition D uses byte-identical ligand files to condition C, and condition E uses byte-identical ligand files to condition A. Any difference between D and C is then attributable to parameters alone, with no re-preparation noise.

[?] Confirm on the node that the skill-route and naive-route ligand PDBQT sets are both retained for all ten targets, including for compounds the filter removed in condition B. Condition C required re-preparing the filtered compounds under the skill route, so those files should exist. Condition E needs the naive-route files for the same molecules; if the naive run never attempted the filtered compounds, the E population is the naive-attempted set and the comparison must be made on the intersection of all four conditions. The runner records which case applies.

## 3. Population and evaluation

[J] Primary comparison population: compounds scored in all four conditions (A, C, D, E). This is the only population on which the 2 x 2 is balanced. Report its size, active count and decoy count per target.

[J] Secondary: the end-to-end population with unscored compounds ranked last, for continuity with the manuscript's primary endpoint. The filter is off in all four conditions, so no deliberate exclusions enter here; only technical failures do.

[V] All AUC values use the tie-corrected rank statistic. The analysis script implements fast DeLong with midranks, which is the tie-corrected estimator and supplies the covariance for the paired test.

Statistics, matching the manuscript: paired DeLong per target per contrast; DerSimonian-Laird random-effects pooling across targets with Cochran Q and I-squared; Benjamini-Hochberg across the ten targets within each contrast; no correction across contrasts, stated as a choice.

## 4. Targets

[V] Ten multi-target targets: EGFR, CDK2, COX-2, ESR1, DPP-4, AChE, BACE-1, HSP90-alpha, p38 MAPK, thrombin. FPR2 is excluded, as in the published decomposition, because its archive layout has no comparable skill receptor.

## 5. Engine and parameters

[V] Uni-Dock v1.1.3, batches of 100 ligands with a 20-ligand sub-batch retry on failure. Box centre computed per target from the centroid of the co-crystallised ligand HETATM coordinates.

[?] The exact Uni-Dock command line must be checked against `unidock --help` on the node before the full run; the runner prints the command it will issue and refuses to proceed past the smoke test unless `--dry-run` output has been inspected. Do not assume the flag names from the published scripts match the installed build.

Scoring function: Vina, as published. No rescoring.

## 6. Cost

[V] Published runs used two RTX 4500 Ada cards under CUDA 12.9 and completed a 10^4 to 10^5 compound library in tens of minutes to a few hours per GPU. The skill protocol's larger box and higher search effort raise per-ligand cost roughly fourfold.

[J] Estimated work here: about 2,000 to 4,100 ligands per target per condition (the no-filter population), ten targets, two conditions. Condition D is at naive cost (cheap), condition E at skill cost (about four times D). Total is roughly five condition-equivalents of naive docking, which is comparable to the published condition A plus C effort. On 8 x RTX 2080 Ti with no scheduler, running one target per idle GPU, the expectation is well under a day of wall-clock for both conditions. The smoke test measures this rather than assuming it.

[?] The 2080 Ti is compute capability 7.5 with 11 GB. Uni-Dock builds for this architecture, but confirm the installed binary runs on it with the smoke test before committing.

## 7. Controls

[J] Three checks, all run by the scripts, because a new run on different hardware must be shown to be comparable with the published one:

1. **Positive control on reuse.** Re-dock a 500-ligand sample of condition C on this hardware and compare the score distribution and AUC with the archived condition C table for that target. Agreement to within docking stochasticity validates the hardware and build; disagreement beyond that must be explained before D and E are interpreted. Uni-Dock is stochastic, so expect small differences and report them rather than expecting equality.
2. **Input integrity.** SHA-256 every reused PDBQT and record it in the manifest, so the files used here are provably the archived ones.
3. **Submission completeness.** Every ligand in the manifest must appear in the output with either a score or a recorded failure reason. The runner fails loudly on silent drops, which is the defect the audit found in the published failure accounting.

## 8. Outputs

Per target per condition: `scores_<target>_<condition>.csv` with columns `mol_id, chembl_id, label, score, status`. Plus `manifest_<target>_<condition>.json` with input hashes, the exact command line, engine and driver versions, GPU model, timing and failure counts.

These transfer back and `analyse_cross_conditions.py` runs locally to produce the decomposition table and figure.

## 9. Order of operations

1. Resolve node access and the archive location.
2. `python run_cross_conditions.py --check` to verify inputs, hashes and the ligand-set question in section 2.
3. `python run_cross_conditions.py --smoke <target>` for one target, 300 ligands, both conditions. Inspect throughput and output format.
4. Positive control (section 7.1) on one target.
5. Full run, ten targets, conditions D and E.
6. Transfer score tables and manifests.
7. `python analyse_cross_conditions.py` locally.

Nothing in the manuscript changes until step 7 produces numbers.
