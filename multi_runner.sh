#!/usr/bin/env bash
# Run the multi-galaxy extension simulations sequentially.
#
# Outputs land in data/multi_*.npz; per-run logs in data/multi_*.log.
# All runs use the N=40k baseline particle count.
#
# CMPH Project 3 -- extension/multi-galaxy.

set -euo pipefail
cd "$(dirname "$0")"
source ../.venv/bin/activate

mkdir -p data
echo "=== multi_runner started at $(date) ==="

run_unequal() {
    local label="$1"; local ratio="$2"
    local out="data/multi_unequal_${label}.npz"
    local log="data/multi_unequal_${label}.log"
    if [[ -f "$out" ]]; then
        echo "[$(date '+%H:%M:%S')] skipping unequal_$label (exists)"
        return
    fi
    echo "[$(date '+%H:%M:%S')] running unequal_$label (ratio = $ratio) -> $out"
    python evolve_unequal_merger.py --mass-ratio "$ratio" --out "$out" \
        > "$log" 2>&1
    echo "[$(date '+%H:%M:%S')] finished unequal_$label"
}

run_three_galaxy() {
    local label="$1"; shift
    local out="data/multi_three_${label}.npz"
    local log="data/multi_three_${label}.log"
    if [[ -f "$out" ]]; then
        echo "[$(date '+%H:%M:%S')] skipping three_$label (exists)"
        return
    fi
    echo "[$(date '+%H:%M:%S')] running three_$label -> $out"
    python evolve_three_galaxy.py "$@" --out "$out" > "$log" 2>&1
    echo "[$(date '+%H:%M:%S')] finished three_$label"
}

# --- Unequal-mass scan ---
run_unequal 1to2  0.5
run_unequal 1to4  0.25
run_unequal 1to8  0.125

# --- Three-galaxy run ---
run_three_galaxy default

echo "=== multi_runner finished at $(date) ==="
