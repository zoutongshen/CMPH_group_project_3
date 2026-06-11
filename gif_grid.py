"""
Multi-trajectory animated comparison grid (master gif).

Renders an animation with one face-on panel per trajectory, all
synchronised to the same master simulation time. The master time runs
import os
from 0 to the shortest trajectory's tstop, so every panel shows
meaningful (in-range) data at every frame; each panel's actual frame
index is chosen as the nearest snapshot to the master time.

Run as a script:
    python gif_grid.py --kind pericentre
    python gif_grid.py --kind inclination --view-kpc 90
    python gif_grid.py --kind multi --view-kpc 90

CMPH Project 3 -- extensions.
"""

import argparse
import shutil
from typing import List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

from merger_analysis import robust_stellar_mask
from physical_units import milky_way_unit_system


PERICENTRE_TRAJECTORIES = [
    ("data/scan_peri1.npz",        r"$r_p = 1$"),
    ("data/merger_N40k_eps01.npz", r"$r_p = 5$ (default)"),
    ("data/scan_peri10.npz",       r"$r_p = 10$"),
    ("data/scan_peri20.npz",       r"$r_p = 20$"),
]

INCLINATION_TRAJECTORIES = [
    ("data/scan_incl0.npz",        r"$i = 0\degree$ (coplanar)"),
    ("data/merger_N40k_eps01.npz", r"$i = 30\degree$ (default)"),
    ("data/scan_incl60.npz",       r"$i = 60\degree$"),
    ("data/scan_incl90.npz",       r"$i = 90\degree$ (polar)"),
    ("data/scan_incl180.npz",      r"$i = 180\degree$ (retrograde)"),
]

MULTI_TRAJECTORIES = [
    ("data/merger_N40k_eps01.npz",       "1:1 (equal mass)"),
    ("data/multi_unequal_1to2.npz",      "1:2"),
    ("data/multi_unequal_1to4.npz",      "1:4"),
    ("data/multi_unequal_1to8.npz",      "1:8"),
]

ZAXIS_TRAJECTORIES = [
    ("data/multi_zaxis_i0.npz",          r"z-axis, $i = 0\degree$"),
    ("data/multi_zaxis_i90.npz",         r"z-axis, $i = 90\degree$"),
    ("data/multi_zaxis_i180.npz",        r"z-axis, $i = 180\degree$"),
]


GALAXY_COLOURS = ["#5fa8ff", "#ff7a5f", "#6fdc8c", "#d18cff", "#ffd76f"]


def load_trajectory_for_animation(
    path: str, label: str, max_points_per_galaxy: int, rng: np.random.Generator,
) -> dict:
    """Open one trajectory file and pre-compute the stellar mask + decimation."""
    trajectory = np.load(path)
    galaxy_id = trajectory["galaxy_id"]
    num_disk = int(trajectory["num_disk"])
    num_bulge = int(trajectory["num_bulge"])
    num_halo = int(trajectory["num_halo"])
    stars = robust_stellar_mask(galaxy_id, num_disk, num_bulge, num_halo)
    star_galaxy = galaxy_id[stars]

    unique_galaxies = [int(gid) for gid in np.unique(star_galaxy)]
    show_per_galaxy = {}
    for gid in unique_galaxies:
        indices_this = np.flatnonzero(star_galaxy == gid)
        if indices_this.size > max_points_per_galaxy:
            indices_this = indices_this[
                rng.choice(indices_this.size, size=max_points_per_galaxy, replace=False)
            ]
        show_per_galaxy[gid] = indices_this

    return {
        "label": label,
        "positions_history": trajectory["positions_history"],
        "snapshot_times": trajectory["snapshot_times"],
        "stars": stars,
        "unique_galaxies": unique_galaxies,
        "show_per_galaxy": show_per_galaxy,
    }


def build_master_animation(
    trajectories: List[Tuple[str, str]],
    output_path: str,
    *,
    view_kpc: float,
    max_points_per_galaxy: int,
    frames: int,
    frames_per_second: int,
) -> None:
    """Render the side-by-side master animation as a gif (or mp4)."""
    units = milky_way_unit_system()
    rng = np.random.default_rng(0)

    loaded = []
    for path, label in trajectories:
        try:
            loaded.append(
                load_trajectory_for_animation(path, label, max_points_per_galaxy, rng)
            )
        except FileNotFoundError:
            print(f"warning: skipping missing trajectory {path}")

    if not loaded:
        raise RuntimeError("no input trajectories available")

    # Master clock = uniform sweep over the shortest tstop.
    master_tstop = min(traj["snapshot_times"][-1] for traj in loaded)
    master_times = np.linspace(0.0, master_tstop, frames)

    figure, axes = plt.subplots(
        2, len(loaded),
        figsize=(4.5 * len(loaded), 9.5),
        squeeze=False,
    )
    figure.patch.set_facecolor("#0a0a0a")

    # axes[0, col] is face-on (vertical_axis = 1, y); axes[1, col] is edge-on (vertical_axis = 2, z).
    scatters_by_column = []
    for col, traj in enumerate(loaded):
        per_column = {}
        for row, vertical_axis in enumerate((1, 2)):
            axis = axes[row, col]
            axis.set_xlim(-view_kpc, view_kpc)
            axis.set_ylim(-view_kpc, view_kpc)
            axis.set_aspect("equal")
            axis.set_facecolor("#0a0a0a")
            axis.tick_params(colors="#888888", labelsize=11)
            for spine in axis.spines.values():
                spine.set_color("#444444")
            per_column[vertical_axis] = {}
            for gid in traj["unique_galaxies"]:
                per_column[vertical_axis][gid] = axis.scatter(
                    [], [], s=0.5,
                    c=GALAXY_COLOURS[gid % len(GALAXY_COLOURS)],
                    alpha=0.5, linewidths=0.0,
                )
        axes[0, col].set_title(traj["label"], color="#ffffff", fontsize=16)
        axes[1, col].set_xlabel("x  [kpc]", color="#cccccc", fontsize=13)
        scatters_by_column.append(per_column)

    axes[0, 0].set_ylabel("Face-on\ny  [kpc]", color="#cccccc", fontsize=14)
    axes[1, 0].set_ylabel("Edge-on\nz  [kpc]", color="#cccccc", fontsize=14)
    title = figure.suptitle("", color="#ffffff", fontsize=18)

    def update(frame_index: int):
        master_t = master_times[frame_index]
        gigayears = master_t * units.time_year / 1.0e9
        for traj, per_column in zip(loaded, scatters_by_column):
            snap = int(np.argmin(np.abs(traj["snapshot_times"] - master_t)))
            positions = (
                traj["positions_history"][snap][traj["stars"]].astype(np.float64)
                * units.length_kpc
            )
            for vertical_axis, scatters in per_column.items():
                for gid in traj["unique_galaxies"]:
                    indices = traj["show_per_galaxy"][gid]
                    scatters[gid].set_offsets(
                        np.column_stack([
                            positions[indices, 0],
                            positions[indices, vertical_axis],
                        ])
                    )
        title.set_text(f"t = {master_t:.1f}  ({gigayears:.2f} Gyr)")
        return ()

    movie = animation.FuncAnimation(figure, update, frames=frames, blit=False)
    if shutil.which("ffmpeg") is not None and output_path.endswith(".mp4"):
        writer = animation.FFMpegWriter(fps=frames_per_second, bitrate=4000)
    else:
        if not output_path.endswith(".gif"):
            output_path = output_path.rsplit(".", 1)[0] + ".gif"
        writer = animation.PillowWriter(fps=frames_per_second)

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    movie.save(output_path, writer=writer, dpi=90)
    plt.close(figure)
    print(f"Saved {output_path}")


def main() -> None:
    """Parse CLI args and render the master animation."""
    parser = argparse.ArgumentParser(
        description="Synchronised side-by-side merger animation grid."
    )
    parser.add_argument(
        "--kind", choices=("pericentre", "inclination", "multi", "zaxis"),
        default="pericentre",
    )
    parser.add_argument("--view-kpc", type=float, default=120.0)
    parser.add_argument("--max-points", type=int, default=10000)
    parser.add_argument("--frames", type=int, default=80)
    parser.add_argument("--fps", type=int, default=15)
    parser.add_argument("--out", type=str, default="")
    args = parser.parse_args()

    if args.kind == "pericentre":
        trajectories = PERICENTRE_TRAJECTORIES
        default_out = "figures/pericentre_scan/scan_pericentre_grid.gif"
    elif args.kind == "inclination":
        trajectories = INCLINATION_TRAJECTORIES
        default_out = "figures/inclination_scan/scan_inclination_grid.gif"
    elif args.kind == "multi":
        trajectories = MULTI_TRAJECTORIES
        default_out = "figures/mass_ratio_scan/multi_unequal_grid.gif"
    else:  # "zaxis"
        trajectories = ZAXIS_TRAJECTORIES
        default_out = "figures/z_axis/zaxis_grid.gif"

    build_master_animation(
        trajectories,
        args.out or default_out,
        view_kpc=args.view_kpc,
        max_points_per_galaxy=args.max_points,
        frames=args.frames,
        frames_per_second=args.fps,
    )


if __name__ == "__main__":
    main()
