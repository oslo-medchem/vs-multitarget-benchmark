#!/usr/bin/env bash
# Run all remaining docking (Stage 4) and evaluation (Stage 5) for targets
# that have completed Stage 3 but not yet been docked.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMMON="$SCRIPT_DIR/common"

eval "$(conda shell.bash hook)" && conda activate vina_dock

# Targets that need naive+skill docking
# (skip egfr, cdk2, cox2, dpp4, esr1 which are being handled separately)
TARGETS="ache bace1 hsp90 p38 thrombin"

for target in $TARGETS; do
    TARGET_DIR="$SCRIPT_DIR/targets/$target"

    # Read box center
    CX=$(python3 -c "import json; d=json.load(open('$TARGET_DIR/box_center.json')); print(d['center_x'])")
    CY=$(python3 -c "import json; d=json.load(open('$TARGET_DIR/box_center.json')); print(d['center_y'])")
    CZ=$(python3 -c "import json; d=json.load(open('$TARGET_DIR/box_center.json')); print(d['center_z'])")

    echo ""
    echo "============================================================"
    echo "  $target — Docking (center: $CX, $CY, $CZ)"
    echo "============================================================"

    # Build ligand lists
    find "$TARGET_DIR/03_library_preparation/naive/pdbqt" -name "*.pdbqt" -size +50c | sort \
        > "$TARGET_DIR/04_docking/naive/ligand_list_naive.txt"
    find "$TARGET_DIR/03_library_preparation/skill/pdbqt" -name "*.pdbqt" -size +50c | sort \
        > "$TARGET_DIR/04_docking/skill/ligand_list_skill.txt"

    N_NAIVE=$(wc -l < "$TARGET_DIR/04_docking/naive/ligand_list_naive.txt")
    N_SKILL=$(wc -l < "$TARGET_DIR/04_docking/skill/ligand_list_skill.txt")
    echo "  Naive ligands: $N_NAIVE, Skill ligands: $N_SKILL"

    # Run naive on GPU 0
    echo "  --- Naive docking (GPU 0) ---"
    python3 "$COMMON/run_unidock_batched.py" \
        --receptor "$TARGET_DIR/04_docking/naive/receptor_naive.pdbqt" \
        --ligand-list "$TARGET_DIR/04_docking/naive/ligand_list_naive.txt" \
        --center-x $CX --center-y $CY --center-z $CZ \
        --size 20 --exhaustiveness 8 --num-modes 9 \
        --results-dir "$TARGET_DIR/04_docking/naive/results" \
        --gpu 0 --batch-size 100 2>&1 | grep -E "^(Total|Done|Failed|Success)"

    # Run skill on GPU 1
    echo "  --- Skill docking (GPU 1) ---"
    python3 "$COMMON/run_unidock_batched.py" \
        --receptor "$TARGET_DIR/04_docking/skill/receptor_skill.pdbqt" \
        --ligand-list "$TARGET_DIR/04_docking/skill/ligand_list_skill.txt" \
        --center-x $CX --center-y $CY --center-z $CZ \
        --size 25 --exhaustiveness 32 --num-modes 9 \
        --results-dir "$TARGET_DIR/04_docking/skill/results" \
        --gpu 1 --batch-size 100 2>&1 | grep -E "^(Total|Done|Failed|Success)"

    # Parse scores
    echo "  --- Parse scores ---"
    python3 "$COMMON/parse_scores.py" --target-dir "$TARGET_DIR" 2>&1 | tail -5

    # Compute metrics
    echo "  --- Compute metrics ---"
    python3 "$COMMON/compute_metrics.py" --target-dir "$TARGET_DIR" 2>&1 | grep -E "ROC AUC|Summary|==="

    # Statistical analysis
    echo "  --- Statistical analysis ---"
    python3 "$COMMON/statistical_analysis.py" --target-dir "$TARGET_DIR" 2>&1 | grep -E "ΔAUC|Z=|Done"

    # Generate figures
    echo "  --- Figures ---"
    python3 "$COMMON/generate_figures.py" --target-dir "$TARGET_DIR" --target-name "$target" 2>&1 | tail -3

    echo "  $target COMPLETE"
done

echo ""
echo "============================================================"
echo "  All remaining targets docked and evaluated!"
echo "============================================================"
