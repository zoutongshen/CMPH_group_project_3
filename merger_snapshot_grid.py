"""
4-panel grid showing the merger morphology at key moments.

Renders a single PNG with face-on views of the merger at four times:
``before`` (well-separated), ``bridge`` (close passage forming the
inter-galaxy bridge), ``tail`` (post-pericentre with extended tidal
tails), ``after`` (late-time damped orbit / coalescence). Stellar
particles (disk + bulge) only; halo is excluded so the visible
morphology is what an observer would see.

Run as a script:
    python merger_snapshot_grid.py data/merger_N160k_eps01.npz
    python merger_snapshot_grid.py data/merger_N160k_eps01.npz \\
        --times 0 25.5 45 150 --view-kpc 90

CMPH Project 3 -- Zoutong Shen / Zhaoyang Chu, 2026.
"""

import argparse
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from merger_analysis import stellar_mask
from physical_units import milky_way_unit_system


def robust_stellar_mask(
    galaxy_id: np.ndarray,
    num_disk: int,
    num_bulge: int,
    num_halo: int,
) -> np.ndarray:
    """
    Stellar (disk + bulge) mask that handles any number of galaxies and
    unequal-mass mergers, in which only galaxy 1's component counts are
    stored in the trajectory metadata. Assumes each galaxy's particles
    are contiguous in the array (disk, bulge, halo in that order, as
    produced by make_galaxy_initial_conditions) and that the disk:bulge:
    halo ratio is preserved when total mass is scaled.
    """
    reference_per_galaxy = num_disk + num_bulge + num_halo
    stellar_fraction = (num_disk + num_bulge) / reference_per_galaxy
    mask = np.zeros(galaxy_id.size, dtype=bool)
    for gid in np.unique(galaxy_id):
        indices = np.flatnonzero(galaxy_id == gid)
        n_stellar = int(round(stellar_fraction * indices.size))
        mask[indices[:n_stellar]] = True
    return mask


def nearest_snapshot_index(
    snapshot_times: np.ndarray, target_time: float
) -> int:
    """Return the index of the snapshot closest to ``target_time``."""
    return int(np.argmin(np.abs(snapshot_times - target_time)))


def render_grid(
    trajectory: dict,
    output_path: str,
    *,
    times_code: Sequence[float],
    labels: Sequence[str],
    view_kpc: float,
    max_points_per_galaxy: int,
) -> None:
    """Build the 4-panel figure and write it to ``output_path``."""
    units = milky_way_unit_system()
    positions_history = trajectory["positions_history"]
    galaxy_id = trajectory["galaxy_id"]
    snapshot_times = trajectory["snapshot_times"]
    num_disk = int(trajectory["num_disk"])
    num_bulge = int(trajectory["num_bulge"])
    num_halo = int(trajectory["num_halo"])

    stars = robust_stellar_mask(galaxy_id, num_disk, num_bulge, num_halo)
    star_galaxy = galaxy_id[stars]

    rng = np.random.default_rng(0)
    galaxy_colours = ["#5fa8ff", "#ff7a5f", "#6fdc8c", "#d18cff", "#ffd76f"]
    unique_galaxies = [int(g) for g in np.unique(star_galaxy)]

    def decimate(indices: np.ndarray) -> np.ndarray:
        if indices.size <= max_points_per_galaxy:
            return indices
        return indices[
            rng.choice(indices.size, size=max_points_per_galaxy, replace=False)
        ]

    show_per_galaxy = {
        gid: decimate(np.flatnonzero(star_galaxy == gid))
        for gid in unique_galaxies
    }

    figure, axes = plt.subplots(
        2, len(times_code),
        figsize=(4.5 * len(times_code), 9.5),
    )
    figure.patch.set_facecolor("#0a0a0a")

    for col, (target_t, label) in enumerate(zip(times_code, labels)):
        index = nearest_snapshot_index(snapshot_times, target_t)
        actual_t = float(snapshot_times[index])
        gigayears = actual_t * units.time_year / 1.0e9
        positions = positions_history[index][stars].astype(np.float64) * units.length_kpc

        # Row 0: face-on (x-y); row 1: edge-on (x-z).
        for row, vertical_axis in enumerate((1, 2)):
            axis = axes[row, col]
            for gid in unique_galaxies:
                indices = show_per_galaxy[gid]
                axis.scatter(
                    positions[indices, 0],
                    positions[indices, vertical_axis],
                    s=0.4,
                    c=galaxy_colours[gid % len(galaxy_colours)],
                    alpha=0.45, linewidths=0.0,
                )
            axis.set_xlim(-view_kpc, view_kpc)
            axis.set_ylim(-view_kpc, view_kpc)
            axis.set_aspect("equal")
            axis.set_facecolor("#0a0a0a")
            axis.tick_params(colors="#cccccc", labelsize=8)
            for spine in axis.spines.values():
                spine.set_color("#444444")

        axes[0, col].set_title(
            f"{label}\nt = {actual_t:.1f}  ({gigayears:.2f} Gyr)",
            color="#ffffff", fontsize=11,
        )
        axes[1, col].set_xlabel("x  [kpc]", color="#cccccc", fontsize=9)

    axes[0, 0].set_ylabel("Face-on\ny  [kpc]", color="#cccccc", fontsize=10)
    axes[1, 0].set_ylabel("Edge-on\nz  [kpc]", color="#cccccc", fontsize=10)

    plt.tight_layout()
    figure.savefig(output_path, dpi=130, facecolor=figure.get_facecolor())
    plt.close(figure)
    print(f"Saved {output_path}")


def main() -> None:
    """Parse CLI args and render the grid."""
    parser = argparse.ArgumentParser(
        description="4-panel merger morphology grid (face-on)."
    )
    parser.add_argument("trajectory", type=str)
    parser.add_argument(
        "--times",
        type=float,
        nargs=4,
        default=[0.0, 25.5, 45.0, 150.0],
        help="four code-time targets: before bridge tail after",
    )
    parser.add_argument(
        "--labels",
        type=str,
        nargs=4,
        default=["Before", "Bridge (pericentre)", "Tail", "After"],
    )
    parser.add_argument("--view-kpc", type=float, default=90.0)
    parser.add_argument("--max-points", type=int, default=30000)
    parser.add_argument(
        "--out",
        type=str,
        default="",
        help="output PNG (default: figures/<trajectory-stem>_4panel.png)",
    )
    args = parser.parse_args()

    trajectory = np.load(args.trajectory)
    if args.out:
        output_path = args.out
    else:
        stem = args.trajectory.rsplit("/", 1)[-1].rsplit(".", 1)[0]
        output_path = f"figures/{stem}_4panel.png"

    render_grid(
        trajectory,
        output_path,
        times_code=args.times,
        labels=args.labels,
        view_kpc=args.view_kpc,
        max_points_per_galaxy=args.max_points,
    )


if __name__ == "__main__":
    main()
