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

subfolder_for() {
    # Map npz filename to its visualisation subfolder.
    local stem="$1"
    case "$stem" in
        scan_peri*)        echo "pericentre_scan" ;;
        scan_incl*)        echo "inclination_scan" ;;
        multi_unequal_*)   echo "mass_ratio_scan" ;;
        multi_three_*)     echo "three_galaxy" ;;
        merger_N*|merger_test*) echo "baseline_merger" ;;
        z_*|multi_zaxis_*) echo "z_axis" ;;
        *)                 echo "" ;;
    esac
}

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
    local subfolder
    subfolder="$(subfolder_for "$stem")"
    if [[ -z "$subfolder" ]]; then
        echo "[skip] $stem has no configured subfolder"
        return
    fi
    mkdir -p "figures/$subfolder"
    local gif="figures/${subfolder}/${stem}.gif"
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
animate_one data/merger_N40k_eps01.npz    # default = peri 5 (lands in baseline_merger/)
animate_one data/scan_peri10.npz
animate_one data/scan_peri20.npz

# Inclination scan
animate_one data/scan_incl0.npz
# default = incl 30 (same file as peri 5, already animated above)
animate_one data/scan_incl60.npz
animate_one data/scan_incl90.npz
animate_one data/scan_incl180.npz

# Mass-ratio scan
animate_one data/multi_unequal_1to2.npz
animate_one data/multi_unequal_1to4.npz
animate_one data/multi_unequal_1to8.npz

# Three-galaxy run
animate_one data/multi_three_default.npz

# z-axis collision experiments
animate_one data/multi_zaxis_i0.npz
animate_one data/multi_zaxis_i90.npz
animate_one data/multi_zaxis_i180.npz
animate_one data/multi_zaxis_headon.npz

echo "=== scan_gif_runner finished at $(date) ==="
