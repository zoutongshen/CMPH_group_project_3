#!/usr/bin/env bash
# Generate face-on + edge-on gif animations for every available scan run.
#
# Idempotent: skips outputs that already exist, skips inputs that don't
# exist yet. Re-run as scan_runner.sh finishes more trajectories.
#
# CMPH Project 3 -- extension/orbit-scans.

set -euo pipefail
cd "$(dirname "$0")"
VENV_PATH="${CMPH_VENV:-../.venv}"
source "$VENV_PATH/bin/activate"

mkdir -p figures
echo "=== scan_gif_runner started at $(date) ==="

animate_one() {
    local npz="$1"
    if [[ ! -f "$npz" ]]; then
        echo "[skip] $npz not yet available"
        return
    fi
    # If the file was modified in the last 30 s, the simulator may still
    # be writing to it -- skip rather than risk a partial-read BadZipFile.
    # Re-running scan_gif_runner.sh later will pick it up.
    if find "$npz" -mmin -0.5 2>/dev/null | grep -q .; then
        echo "[skip] $npz modified within the last 30 s, retry later"
        return
    fi
    local stem
    stem="$(basename "$npz" .npz)"
    local gif="figures/${stem}.gif"
    if [[ -f "$gif" ]]; then
        echo "[skip] $gif already exists"
        return
    fi
    echo "[$(date '+%H:%M:%S')] animating $npz -> $gif"
    python merger_animation.py "$npz" --stride 2 --max-points 15000 \
        --view-kpc 90 --fps 15 --out "$gif"
}

# Pericentre scan
animate_one data/scan_peri1.npz
animate_one data/merger_N40k_eps01.npz    # default = peri 5
animate_one data/scan_peri10.npz
animate_one data/scan_peri20.npz

# Inclination scan
animate_one data/scan_incl0.npz
# default = incl 30 (same file as peri 5, already animated)
animate_one data/scan_incl60.npz
animate_one data/scan_incl90.npz
animate_one data/scan_incl180.npz

echo "=== scan_gif_runner finished at $(date) ==="
