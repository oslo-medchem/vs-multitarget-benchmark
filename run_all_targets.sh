#!/usr/bin/env bash
# =============================================================================
# Master orchestrator: run all targets through the VS benchmark pipeline
#
# Usage:
#   bash run_all_targets.sh                        # all stages, max 4 parallel
#   bash run_all_targets.sh --stages 1,2,3         # prep stages only
#   bash run_all_targets.sh --max-parallel 2       # limit parallelism
#   bash run_all_targets.sh --stages 4 --gpu 0     # docking on GPU 0
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$SCRIPT_DIR/target_config.json"

# ---- Defaults ----
STAGES="1,2,3,4,5"
MAX_PARALLEL=4
GPU=""

# ---- Parse arguments ----
while [[ $# -gt 0 ]]; do
    case "$1" in
        --stages)
            STAGES="$2"; shift 2 ;;
        --max-parallel)
            MAX_PARALLEL="$2"; shift 2 ;;
        --gpu)
            GPU="$2"; shift 2 ;;
        --help|-h)
            echo "Usage: $0 [--stages 1,2,3,4,5] [--max-parallel 4] [--gpu DEVICE_ID]"
            exit 0 ;;
        *)
            echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# ---- Activate conda environment ----
eval "$(conda shell.bash hook)" && conda activate vina_dock

# ---- Discover pending targets ----
TARGETS=$(python3 -c "
import json, sys
c = json.load(open('$CONFIG'))
for t in c['targets']:
    if t.get('status') != 'complete':
        print(t['name'])
")

if [ -z "$TARGETS" ]; then
    echo "No pending targets found in $CONFIG"
    exit 0
fi

TARGET_COUNT=$(echo "$TARGETS" | wc -l)

echo "============================================================"
echo "  Multi-Target Virtual Screening Benchmark"
echo "============================================================"
echo "Config:       $CONFIG"
echo "Targets (${TARGET_COUNT}): $(echo $TARGETS | tr '\n' ' ')"
echo "Stages:       $STAGES"
echo "Max parallel: $MAX_PARALLEL"
[ -n "$GPU" ] && echo "GPU device:   $GPU"
echo ""

PIPELINE_START=$(date +%s)

# Helper: extract stages matching a regex from the STAGES string
filter_stages() {
    local pattern="$1"
    echo "$STAGES" | tr ',' '\n' | grep -E "$pattern" | tr '\n' ',' | sed 's/,$//'
}

# ---- Phase 1: CPU-bound stages (1, 2, 3) ----
CPU_STAGES=$(filter_stages '^[123]$')
if [ -n "$CPU_STAGES" ]; then
    echo "============================================================"
    echo "  Phase 1: Preparation (CPU) — stages $CPU_STAGES"
    echo "============================================================"

    for target in $TARGETS; do
        # Ensure target log directory exists
        mkdir -p "$SCRIPT_DIR/targets/$target"

        echo "[$(date '+%H:%M:%S')] Starting preparation for $target ..."
        python3 "$SCRIPT_DIR/run_target_pipeline.py" \
            --target "$target" \
            --stages "$CPU_STAGES" \
            --config "$CONFIG" \
            --skip-existing \
            2>&1 | tee "$SCRIPT_DIR/targets/$target/pipeline.log" &

        # Throttle parallel jobs
        while [ "$(jobs -r | wc -l)" -ge "$MAX_PARALLEL" ]; do
            wait -n 2>/dev/null || true
        done
    done
    wait
    echo ""
    echo "Phase 1 complete."
    echo ""
fi

# ---- Phase 2: GPU docking (stage 4) ----
if echo "$STAGES" | tr ',' '\n' | grep -qx '4'; then
    echo "============================================================"
    echo "  Phase 2: Docking (GPU) — stage 4"
    echo "============================================================"

    GPU_ARG=""
    [ -n "$GPU" ] && GPU_ARG="--gpu $GPU"

    for target in $TARGETS; do
        echo "[$(date '+%H:%M:%S')] Docking $target ..."
        python3 "$SCRIPT_DIR/run_target_pipeline.py" \
            --target "$target" \
            --stages 4 \
            --config "$CONFIG" \
            $GPU_ARG \
            --skip-existing \
            2>&1 | tee -a "$SCRIPT_DIR/targets/$target/pipeline.log"
    done
    echo ""
    echo "Phase 2 complete."
    echo ""
fi

# ---- Phase 3: Evaluation (stage 5) ----
if echo "$STAGES" | tr ',' '\n' | grep -qx '5'; then
    echo "============================================================"
    echo "  Phase 3: Evaluation — stage 5"
    echo "============================================================"

    for target in $TARGETS; do
        echo "[$(date '+%H:%M:%S')] Evaluating $target ..."
        python3 "$SCRIPT_DIR/run_target_pipeline.py" \
            --target "$target" \
            --stages 5 \
            --config "$CONFIG" \
            2>&1 | tee -a "$SCRIPT_DIR/targets/$target/pipeline.log"
    done
    echo ""
    echo "Phase 3 complete."
    echo ""
fi

# ---- Phase 4: Cross-target meta-analysis ----
echo "============================================================"
echo "  Phase 4: Meta-Analysis"
echo "============================================================"

META_DIR="$SCRIPT_DIR/06_meta_analysis"

if [ -d "$META_DIR" ]; then
    [ -f "$META_DIR/aggregate_metrics.py" ] && \
        python3 "$META_DIR/aggregate_metrics.py" --config "$CONFIG"
    [ -f "$META_DIR/forest_plot.py" ] && \
        python3 "$META_DIR/forest_plot.py"
    [ -f "$META_DIR/paired_tests.py" ] && \
        python3 "$META_DIR/paired_tests.py"
    echo "Meta-analysis complete."
else
    echo "WARNING: $META_DIR not found — skipping meta-analysis."
fi

# ---- Summary ----
PIPELINE_END=$(date +%s)
ELAPSED=$(( PIPELINE_END - PIPELINE_START ))
MINS=$(( ELAPSED / 60 ))
SECS=$(( ELAPSED % 60 ))

echo ""
echo "============================================================"
echo "  All targets complete  (${MINS}m ${SECS}s total)"
echo "============================================================"
echo "Results are in: $SCRIPT_DIR/targets/*/05_evaluation/"
echo "Meta-analysis:  $META_DIR/"
