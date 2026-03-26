#!/usr/bin/env bash
# Re-dock ALL targets' skill protocol with corrected Meeko receptors.
# Then re-run evaluation (parse_scores, compute_metrics, statistical_analysis, figures).
#
# Background: Initial skill docking used obabel-fallback receptors (identical to naive),
# negating the receptor prep difference. Now using Meeko receptors with --default_altloc A -a.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMMON="$SCRIPT_DIR/common"

eval "$(conda shell.bash hook)" && conda activate vina_dock

# Only targets whose skill docking used wrong (obabel-fallback) receptor
# BACE1/HSP90/p38/Thrombin will get correct Meeko receptor from run_remaining_docking.sh
TARGETS="egfr cdk2 cox2 esr1 dpp4 ache"

for target in $TARGETS; do
    TARGET_DIR="$SCRIPT_DIR/targets/$target"

    # Skip if no box center
    if [ ! -f "$TARGET_DIR/box_center.json" ]; then
        echo "SKIP $target: no box_center.json"
        continue
    fi

    CX=$(python3 -c "import json; d=json.load(open('$TARGET_DIR/box_center.json')); print(d['center_x'])")
    CY=$(python3 -c "import json; d=json.load(open('$TARGET_DIR/box_center.json')); print(d['center_y'])")
    CZ=$(python3 -c "import json; d=json.load(open('$TARGET_DIR/box_center.json')); print(d['center_z'])")

    echo ""
    echo "============================================================"
    echo "  $target — Skill re-dock with Meeko receptor"
    echo "  Center: ($CX, $CY, $CZ)"
    echo "============================================================"

    SKILL_DIR="$TARGET_DIR/04_docking/skill"
    RECEPTOR="$SKILL_DIR/receptor_skill.pdbqt"

    # Verify receptor is Meeko-prepared (not obabel fallback)
    if grep -q "active torsion" "$RECEPTOR" 2>/dev/null; then
        echo "  WARNING: Receptor appears to have torsions (flex). Re-preparing..."
        CHAIN_PDB=$(find "$TARGET_DIR/00_setup" -name "*_chain*.pdb" | head -1)
        mk_prepare_receptor.py --read_pdb "$CHAIN_PDB" -p "$RECEPTOR" --default_altloc A -a 2>&1 | tail -1
    fi

    # Build ligand list
    find "$TARGET_DIR/03_library_preparation/skill/pdbqt" -name "*.pdbqt" -size +50c | sort \
        > "$SKILL_DIR/ligand_list_skill.txt"
    N=$(wc -l < "$SKILL_DIR/ligand_list_skill.txt")
    echo "  Ligands: $N"

    # Clear old results
    rm -f "$SKILL_DIR/results/"*_out.pdbqt 2>/dev/null

    # Dock on GPU 1
    python3 "$COMMON/run_unidock_batched.py" \
        --receptor "$RECEPTOR" \
        --ligand-list "$SKILL_DIR/ligand_list_skill.txt" \
        --center-x $CX --center-y $CY --center-z $CZ \
        --size 25 --exhaustiveness 32 --num-modes 9 \
        --results-dir "$SKILL_DIR/results" \
        --gpu 0 --batch-size 100 2>&1 | grep -E "^(Total|Done|Failed|Success)"

    echo "  --- Evaluation ---"
    # Parse scores
    python3 "$COMMON/parse_scores.py" --target-dir "$TARGET_DIR" 2>&1 | grep -E "Parsed|Docked"

    # Compute metrics (new dual-mode)
    python3 "$COMMON/compute_metrics.py" --target-dir "$TARGET_DIR" 2>&1 | grep -E "\[full\]|\[dock\]|Summary"

    # Statistical analysis (with intersection)
    python3 "$COMMON/statistical_analysis.py" --target-dir "$TARGET_DIR" 2>&1 | grep -E "INTERSECTION|ΔAUC|Z=|N comp"

    # Figures
    python3 "$COMMON/generate_figures.py" --target-dir "$TARGET_DIR" --target-name "$target" 2>&1 | grep "All figures"

    echo "  $target COMPLETE"
done

echo ""
echo "============================================================"
echo "  All skill re-docking and evaluation complete!"
echo "============================================================"
