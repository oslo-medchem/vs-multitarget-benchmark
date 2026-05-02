# Prompts and Skill Files for the Multi-Target VS Benchmark

This directory holds the verbatim prompts and skill-file content used to drive Claude Code (Claude Opus 4.6) for the manuscript "*A controlled multi-target benchmark of skill-file-guided agentic virtual screening across eleven drug targets*."

## Contents

| File | Purpose | When loaded |
|---|---|---|
| `SKILL.md` | The 202-line virtual-screening-pipeline skill file. The **only** difference between the naive and skill-guided protocols is whether this file was loaded into Claude Code's system context. | Skill-guided protocol only. |
| `CLAUDE_multitarget.md` | Repository-level operational instructions (target list, environment, pipeline stages, conventions). Used by both protocols (it describes the experiment setup, not domain knowledge). | Both protocols. |

The naive protocol was run with **only** `CLAUDE_multitarget.md` (and the per-target user prompts, see below). The skill-guided protocol was run with `CLAUDE_multitarget.md` **plus** `SKILL.md`. No other prompt-engineering, no few-shot examples, no role-play instructions, and no in-context retrieval were used.

## Per-stage user prompts (both protocols)

For each of the 11 targets, a single natural-language prompt of the following form was issued:

> "Build a complete virtual screening pipeline for **{target_name}** (ChEMBL ID: **{chembl_id}**, structure: **{pdb_code}** chain **{chain}**). Implement the five stages: (1) active curation from ChEMBL with pChEMBL ≥ 5; (2) property-matched decoy generation against the shared 60,000-compound ChEMBL pool; (3) library preparation including receptor and ligand preparation; (4) GPU docking with Uni-Dock v1.1.3; (5) statistical evaluation. Generate all required Python and shell scripts. Place results under `targets/{target_name}/`. Do not edit existing scripts in other targets."

The prompt was identical between protocols up to the target-specific tokens (`{target_name}`, `{chembl_id}`, `{pdb_code}`, `{chain}`); the only experimental variable was the presence of `SKILL.md`.

## Skill-file semantics

The skill file is loaded as **advisory** system context, not as a deterministic constraint. The agent retains discretion to deviate (e.g., fall back to AutoDockTools `prepare_ligand4.py` if Meeko fails on a non-standard ligand). In practice, however, the agent followed the skill's tool-selection and parameter directives consistently across all 18 generated scripts; the skill file therefore behaves operationally as a reliable bias rather than a hard rule. See manuscript Methods §Skill File Semantics for further discussion.

## Skill-file content (high-level outline)

The full file is in `SKILL.md`. The directives most relevant to the experimental outcome are:

- **Receptor preparation:** prefer Meeko `mk_prepare_receptor.py` over OpenBabel `obabel`. Default to ` --default_altloc A`.
- **Ligand 3D generation:** prefer RDKit ETKDGv3 over OpenBabel `gen3d`. Set `randomSeed=42`, `numConfs=1`, embed and minimise with MMFF94s.
- **Substructure filtering:** apply PAINS (Baell & Holloway 2010) and Brenk (Brenk et al. 2008) filters via RDKit's `FilterCatalog` before docking. Compounds with any flagged substructure are excluded.
- **Docking box:** centre on the centroid of the co-crystallised ligand HETATM coordinates; size 25 Å (vs. 20 Å naive default).
- **Exhaustiveness:** 32 (vs. 8 naive default).
- **Batched execution:** Uni-Dock with 100 ligands per batch, 20-ligand sub-batch retry on failure.
- **Evaluation:** ROC AUC, BEDROC(α=20), EF₁%, EF₅%, LogAUC, AUPR; bootstrap BCa CIs with n=10,000; DeLong's test on the intersection of compounds docked by both protocols; DerSimonian–Laird random-effects meta-analysis across the 11 targets.

## License

Both files are released under the MIT license, matching the rest of the benchmark archive.
