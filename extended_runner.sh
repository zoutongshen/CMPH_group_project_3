#!/usr/bin/env bash
# Extended-tstop reruns for three runs whose default tstop = 300-500 cut off
# mid-swing, so it isn't clear whether the orbit decays back to a merger or
# escapes. Re-runs them to tstop = 900 (~11.7 Gyr) and saves to *_long.npz.
#
# CMPH Project 3 -- extension follow-up.

set -euo pipefail
cd "$(dirname "$0")"
VENV_PATH="${CMPH_VENV:-../.venv}"
source "$VENV_PATH/bin/activate"

mkdir -p data
echo "=== extended_runner started at $(date) ==="

run_one() {
    local name="$1"
    local out="data/${name}_long.npz"
    local log="data/${name}_long.log"
    shift
    if [[ -f "$out" ]]; then
        echo "[$(date '+%H:%M:%S')] skipping $name (exists)"
        return
    fi
    echo "[$(date '+%H:%M:%S')] launching $name -> $out"
    python "$@" --tstop 900 --num-dumps 300 --out "$out" > "$log" 2>&1
    echo "[$(date '+%H:%M:%S')] finished $name"
}

# Run the three extended sims in parallel.
run_one scan_peri10        evolve_merger.py         --pericentre 10 &
PID_PERI10=$!
run_one scan_peri20        evolve_merger.py         --pericentre 20 &
PID_PERI20=$!
run_one multi_unequal_1to8 evolve_unequal_merger.py --mass-ratio 0.125 &
PID_1TO8=$!

wait "$PID_PERI10" "$PID_PERI20" "$PID_1TO8"

echo "=== extended_runner finished at $(date) ==="
