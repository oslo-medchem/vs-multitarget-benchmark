#!/usr/bin/env bash
# Naive protocol: Uni-Dock docking
# Box 20A, exhaustiveness=8, Vina scoring
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TARGET_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"

RECEPTOR="$SCRIPT_DIR/receptor_naive.pdbqt"
PDBQT_DIR="$TARGET_DIR/03_library_preparation/naive/pdbqt"
RESULTS_DIR="$SCRIPT_DIR/results"
LIGAND_LIST="$SCRIPT_DIR/ligand_list_naive.txt"

mkdir -p "$RESULTS_DIR"

CENTER_X=-14.108
CENTER_Y=-43.833
CENTER_Z=27.67
SIZE=20
EXHAUSTIVENESS=8
NUM_MODES=9

echo "Building ligand list..."
find "$PDBQT_DIR" -name "*.pdbqt" -size +50c | sort > "$LIGAND_LIST"
N=$(wc -l < "$LIGAND_LIST")
echo "Ligands to dock: $N"

echo "Starting Uni-Dock (naive protocol)..."
START_TIME=$(date +%s)

unidock \
  --receptor "$RECEPTOR" \
  --ligand_index "$LIGAND_LIST" \
  --center_x $CENTER_X \
  --center_y $CENTER_Y \
  --center_z $CENTER_Z \
  --size_x $SIZE \
  --size_y $SIZE \
  --size_z $SIZE \
  --exhaustiveness $EXHAUSTIVENESS \
  --num_modes $NUM_MODES \
  --dir "$RESULTS_DIR" \
  --scoring vina \
  2>&1 | tee "$SCRIPT_DIR/unidock_naive.log"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo "Elapsed: ${ELAPSED}s" | tee -a "$SCRIPT_DIR/unidock_naive.log"
echo "Naive docking complete."
