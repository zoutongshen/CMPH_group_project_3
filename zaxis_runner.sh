#!/usr/bin/env bash
# Three z-axis-collision merger runs (--collision-axis z).
# Disks remain in the x-y plane; the orbit lies in the x-z plane, so the
# two disks approach each other face-on (rather than edge-on as in the
# default x-axis collision).
#
# Three inclination variants:
#   i =   0 : both disks fully face-on to each other (clean smash)
#   i =  90 : one disk perpendicular to the other (one face, one edge)
#   i = 180 : disks anti-aligned (one prograde, one retrograde in orbital frame)
#
# CMPH Project 3 -- extension/multi-galaxy.

set -euo pipefail
cd "$(dirname "$0")"
VENV_PATH="${CMPH_VENV:-../.venv}"
source "$VENV_PATH/bin/activate"

mkdir -p data
echo "=== zaxis_runner started at $(date) ==="

run_one() {
    local label="$1"; shift
    local out="data/multi_zaxis_${label}.npz"
    local log="data/multi_zaxis_${label}.log"
    if [[ -f "$out" ]]; then
        echo "[$(date '+%H:%M:%S')] skipping zaxis_$label (exists)"
        return
    fi
    echo "[$(date '+%H:%M:%S')] running zaxis_$label -> $out"
    python evolve_merger.py --collision-axis z "$@" --out "$out" \
        > "$log" 2>&1
    echo "[$(date '+%H:%M:%S')] finished zaxis_$label"
}

run_one i0   --inclination 0.0
run_one i90  --inclination 90.0
run_one i180 --inclination 180.0

echo "=== zaxis_runner finished at $(date) ==="
