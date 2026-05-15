"""
Face-on and edge-on images for the section 3.0.4 exercise (Q1, Q2).

Reproduces the manual's ``readdata.pro`` / ``plotvine`` views: the disk,
bulge and halo particles plotted in distinct colours, seen both face-on
(x-y) and edge-on (x-z), at the first and last snapshot so the
morphological difference between the thin rotationally-supported disk and
the round pressure-supported bulge (Q1), and the extended halo (Q2), are
directly visible -- and any structural change over t = 0 -> t = 100 too.

Run as a script:
    python disk_visualization.py data/evolve_N20k_eps01.npz

CMPH Project 3 -- Zoutong Shen / Zhaoyang Chu, 2026.
"""

import argparse
from typing import Dict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from disk_profiles import split_components
from physical_units import milky_way_unit_system

# Distinct colours per component (disk light, bulge red, halo faint blue),
# echoing the manual's "disk white / bulge red" convention.
COMPONENT_STYLE = {
    "halo": dict(color="#6fa8dc", size=0.4, alpha=0.25, zorder=1),
    "disk": dict(color="#222222", size=0.4, alpha=0.5, zorder=3),
    "bulge": dict(color="#d62728", size=0.6, alpha=0.6, zorder=2),
}


def _scatter_projection(
    axis,
    components: Dict[str, np.ndarray],
    horizontal_index: int,
    vertical_index: int,
    length_kpc: float,
    half_width_kpc: float,
) -> None:
    """Scatter the three components onto one 2-D projection axis."""
    for name in ("halo", "disk", "bulge"):
        coordinates = components[name]
        style = COMPONENT_STYLE[name]
        axis.scatter(
            coordinates[:, horizontal_index] * length_kpc,
            coordinates[:, vertical_index] * length_kpc,
            s=style["size"],
            c=style["color"],
            alpha=style["alpha"],
            linewidths=0.0,
            zorder=style["zorder"],
            label=name,
        )
    axis.set_xlim(-half_width_kpc, half_width_kpc)
    axis.set_ylim(-half_width_kpc, half_width_kpc)
    axis.set_aspect("equal")


def make_snapshot_figure(
    trajectory: Dict[str, np.ndarray],
    output_png: str,
    halo_view_kpc: float = 80.0,
    disk_view_kpc: float = 25.0,
) -> None:
    """
    Produce a 2x4 panel image: rows = {t=0, t=final},
    columns = {face-on incl. halo, edge-on incl. halo,
               face-on disk+bulge zoom, edge-on disk+bulge zoom}.
    """
    units = milky_way_unit_system()
    positions_history = trajectory["positions_history"]
    masses = trajectory["masses"]
    num_disk = int(trajectory["num_disk"])
    num_bulge = int(trajectory["num_bulge"])
    times = trajectory["snapshot_times"]

    snapshot_indices = [0, positions_history.shape[0] - 1]
    figure, axes = plt.subplots(
        2, 4, figsize=(20, 11), constrained_layout=True
    )

    for row, snapshot_index in enumerate(snapshot_indices):
        positions = positions_history[snapshot_index].astype(np.float64)
        centre = (positions * masses[:, None]).sum(axis=0) / masses.sum()
        positions = positions - centre
        components = split_components(positions, num_disk, num_bulge)

        _scatter_projection(axes[row, 0], components, 0, 1,
                            units.length_kpc, halo_view_kpc)
        _scatter_projection(axes[row, 1], components, 0, 2,
                            units.length_kpc, halo_view_kpc)
        _scatter_projection(axes[row, 2], components, 0, 1,
                            units.length_kpc, disk_view_kpc)
        _scatter_projection(axes[row, 3], components, 0, 2,
                            units.length_kpc, disk_view_kpc)

        time_label = f"t = {times[snapshot_index]:.0f}"
        axes[row, 0].set_ylabel(
            f"{time_label}\n\ny  [kpc]", fontsize=11
        )
        for column in range(4):
            axes[row, column].set_facecolor("#f4f4f4")

    axes[0, 0].set_title("Face-on (all, incl. halo)")
    axes[0, 1].set_title("Edge-on (all, incl. halo)")
    axes[0, 2].set_title("Face-on (disk + bulge zoom)")
    axes[0, 3].set_title("Edge-on (disk + bulge zoom)")
    for column in range(4):
        axes[1, column].set_xlabel("x  [kpc]")
    axes[1, 1].set_ylabel("z  [kpc]")
    axes[1, 3].set_ylabel("z  [kpc]")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles, labels, loc="outside upper center", ncol=3,
        markerscale=8, frameon=False, fontsize=12,
    )
    figure.suptitle(
        f"Isolated disk -- N = {len(masses)}, "
        f"eps = {float(trajectory['softening']):g}  (MW physical scaling)",
        fontsize=14,
    )
    figure.savefig(output_png, dpi=120)
    plt.close(figure)


def main() -> None:
    """Load a trajectory dump and write the Q1/Q2 image."""
    parser = argparse.ArgumentParser(
        description="Face-on / edge-on images for the isolated disk (Q1/Q2)."
    )
    parser.add_argument("trajectory", type=str)
    parser.add_argument("--out", type=str, default="")
    args = parser.parse_args()

    trajectory = np.load(args.trajectory)
    output_png = args.out or args.trajectory.replace(
        ".npz", "_images.png"
    ).replace("data/", "figures/")

    make_snapshot_figure(trajectory, output_png)
    print(f"Saved {output_png}")


if __name__ == "__main__":
    main()
