"""
Multi-trajectory comparison grid for the orbit-scan extension.

Plots a grid of face-on merger snapshots: one row per trajectory, one
column per time slice. Used to compare merger morphology across the
pericentre and inclination scans run by ``scan_runner.sh``.

Run as a script:
    python scan_grid_plot.py --kind pericentre
    python scan_grid_plot.py --kind inclination --times 0 30 60 120

CMPH Project 3 -- extension/orbit-scans.
"""

import os
import argparse
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from merger_analysis import robust_stellar_mask
from physical_units import milky_way_unit_system


PERICENTRE_TRAJECTORIES = [
    ("data/scan_peri1.npz",            r"$r_p = 1$"),
    ("data/merger_N40k_eps01.npz",     r"$r_p = 5$  (default)"),
    ("data/scan_peri10.npz",           r"$r_p = 10$"),
    ("data/scan_peri20.npz",           r"$r_p = 20$"),
]

INCLINATION_TRAJECTORIES = [
    ("data/scan_incl0.npz",            r"$i = 0\degree$  (coplanar)"),
    ("data/merger_N40k_eps01.npz",     r"$i = 30\degree$  (default)"),
    ("data/scan_incl60.npz",           r"$i = 60\degree$"),
    ("data/scan_incl90.npz",           r"$i = 90\degree$  (polar)"),
    ("data/scan_incl180.npz",          r"$i = 180\degree$  (retrograde)"),
]

MULTI_TRAJECTORIES = [
    ("data/merger_N40k_eps01.npz",     "1:1  (equal mass)"),
    ("data/multi_unequal_1to2.npz",    "1:2"),
    ("data/multi_unequal_1to4.npz",    "1:4"),
    ("data/multi_unequal_1to8.npz",    "1:8"),
]

ZAXIS_TRAJECTORIES = [
    ("data/multi_zaxis_i0.npz",        r"z-axis, $i = 0\degree$ (face-on)"),
    ("data/multi_zaxis_i90.npz",       r"z-axis, $i = 90\degree$"),
    ("data/multi_zaxis_i180.npz",      r"z-axis, $i = 180\degree$ (anti-aligned)"),
]


def nearest_index(snapshot_times: np.ndarray, target_time: float) -> int:
    """Index of the snapshot closest to ``target_time``."""
    return int(np.argmin(np.abs(snapshot_times - target_time)))


def draw_panel(
    axis,
    trajectory: dict,
    target_t: float,
    units,
    view_kpc: float,
    rng,
    max_points_per_galaxy: int,
) -> float:
    """Render one face-on panel; return the actual snapshot time used."""
    positions_history = trajectory["positions_history"]
    galaxy_id = trajectory["galaxy_id"]
    snapshot_times = trajectory["snapshot_times"]
    num_disk = int(trajectory["num_disk"])
    num_bulge = int(trajectory["num_bulge"])
    num_halo = int(trajectory["num_halo"])

    stars = robust_stellar_mask(galaxy_id, num_disk, num_bulge, num_halo)
    star_galaxy = galaxy_id[stars]
    galaxy_colours = ["#5fa8ff", "#ff7a5f", "#6fdc8c", "#d18cff", "#ffd76f"]
    unique_galaxies = [int(gid) for gid in np.unique(star_galaxy)]

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

    index = nearest_index(snapshot_times, target_t)
    actual_t = float(snapshot_times[index])

    positions = positions_history[index][stars].astype(np.float64) * units.length_kpc
    for gid in unique_galaxies:
        indices = show_per_galaxy[gid]
        axis.scatter(
            positions[indices, 0],
            positions[indices, 1],
            s=0.4,
            c=galaxy_colours[gid % len(galaxy_colours)],
            alpha=0.45, linewidths=0.0,
        )
    axis.set_xlim(-view_kpc, view_kpc)
    axis.set_ylim(-view_kpc, view_kpc)
    axis.set_aspect("equal")
    axis.set_facecolor("#0a0a0a")
    axis.tick_params(colors="#888888", labelsize=11)
    for spine in axis.spines.values():
        spine.set_color("#444444")
    return actual_t


def render_grid(
    trajectories: Sequence,
    output_path: str,
    *,
    times_code: Sequence[float],
    view_kpc: float,
    max_points_per_galaxy: int,
) -> None:
    """Build the trajectory x time grid and save it."""
    units = milky_way_unit_system()
    rng = np.random.default_rng(0)

    num_rows = len(trajectories)
    num_cols = len(times_code)
    figure, axes = plt.subplots(
        num_rows, num_cols,
        figsize=(3.2 * num_cols, 3.2 * num_rows),
        squeeze=False,
    )
    figure.patch.set_facecolor("#0a0a0a")

    for row, (path, label) in enumerate(trajectories):
        try:
            trajectory = np.load(path)
        except FileNotFoundError:
            for col in range(num_cols):
                axes[row, col].set_facecolor("#0a0a0a")
                axes[row, col].text(
                    0.5, 0.5, f"missing\n{path}",
                    transform=axes[row, col].transAxes,
                    color="#aa6666", fontsize=14, ha="center", va="center",
                )
                axes[row, col].set_xticks([])
                axes[row, col].set_yticks([])
            axes[row, 0].set_ylabel(label, color="#ffffff", fontsize=18)
            continue
        for col, t in enumerate(times_code):
            actual_t = draw_panel(
                axes[row, col], trajectory, t, units,
                view_kpc, rng, max_points_per_galaxy,
            )
            if row == 0:
                gigayears = actual_t * units.time_year / 1.0e9
                axes[row, col].set_title(
                    f"t = {actual_t:.0f}  ({gigayears:.2f} Gyr)",
                    color="#ffffff", fontsize=16,
                )
        axes[row, 0].set_ylabel(label, color="#ffffff", fontsize=18)

    for col in range(num_cols):
        axes[-1, col].set_xlabel("x  [kpc]", color="#cccccc", fontsize=14)
    for row in range(num_rows):
        axes[row, 0].set_ylabel(
            axes[row, 0].get_ylabel() + "\ny  [kpc]",
            color="#cccccc", fontsize=16,
        )

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    figure.savefig(output_path, dpi=120, facecolor=figure.get_facecolor())
    plt.close(figure)
    print(f"Saved {output_path}")


def main() -> None:
    """Parse CLI args and render."""
    parser = argparse.ArgumentParser(
        description="Scan-comparison grid for orbit extension."
    )
    parser.add_argument(
        "--kind", choices=("pericentre", "inclination", "multi", "zaxis"),
        default="pericentre",
        help="which parameter scan to plot",
    )
    parser.add_argument(
        "--times", type=float, nargs="+",
        default=[0.0, 25.5, 60.0, 150.0],
    )
    parser.add_argument("--view-kpc", type=float, default=120.0)
    parser.add_argument("--max-points", type=int, default=15000)
    parser.add_argument("--out", type=str, default="")
    args = parser.parse_args()

    if args.kind == "pericentre":
        trajectories = PERICENTRE_TRAJECTORIES
        default_out = "figures/pericentre_scan/scan_pericentre_grid.png"
    elif args.kind == "inclination":
        trajectories = INCLINATION_TRAJECTORIES
        default_out = "figures/inclination_scan/scan_inclination_grid.png"
    elif args.kind == "multi":
        trajectories = MULTI_TRAJECTORIES
        default_out = "figures/mass_ratio_scan/multi_unequal_grid.png"
    else:  # "zaxis"
        trajectories = ZAXIS_TRAJECTORIES
        default_out = "figures/z_axis/zaxis_grid.png"

    render_grid(
        trajectories,
        args.out or default_out,
        times_code=args.times,
        view_kpc=args.view_kpc,
        max_points_per_galaxy=args.max_points,
    )


if __name__ == "__main__":
    main()
