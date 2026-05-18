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
from stability_test import fit_disk_scale_height

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


def disk_thickening_figure(
    trajectory: Dict[str, np.ndarray],
    output_png: str,
    edge_zoom_kpc: float = 8.0,
    height_zoom_kpc: float = 3.0,
) -> Dict[str, float]:
    """
    Make the disk vertical-thickening figure and return the fitted
    scale heights.

    The wide Q1/Q2 panels make t = 0 and t = final look identical because
    the finite-N vertical heating is only ~1 kpc -- invisible at an
    80 kpc field of view. This figure isolates it: a tight edge-on zoom
    of the disk particles only (t = 0 blue, t = final red) next to the
    normalised vertical density profile n(|z|) at both times, annotated
    with the sech^2 scale height z_0 = sqrt(12 <z^2> / pi^2).
    """
    units = milky_way_unit_system()
    positions_history = trajectory["positions_history"]
    masses = trajectory["masses"]
    num_disk = int(trajectory["num_disk"])
    num_bulge = int(trajectory["num_bulge"])
    times = trajectory["snapshot_times"]
    snapshot_indices = [0, positions_history.shape[0] - 1]

    figure, (axis_edge, axis_profile) = plt.subplots(
        1, 2, figsize=(15, 6), constrained_layout=True
    )
    colours = ("#1f77b4", "#d62728")
    scale_heights = {}

    height_edges = np.linspace(0.0, height_zoom_kpc, 40)
    height_centres = 0.5 * (height_edges[:-1] + height_edges[1:])

    # Measure the disk vertical structure on the raw particle positions,
    # exactly as stability_test does: the IC is built at the origin with
    # zero net momentum, so the disk midplane stays at z = 0. Recentring
    # on the (halo-dominated) global centroid would tilt/shift the disk
    # frame and artificially dilute the measured heating.
    for colour, snapshot_index in zip(colours, snapshot_indices):
        positions = positions_history[snapshot_index].astype(np.float64)
        disk = split_components(positions, num_disk, num_bulge)["disk"]
        time_label = f"t = {times[snapshot_index]:.0f}"

        scale_height_code = fit_disk_scale_height(disk)
        scale_height_kpc = scale_height_code * units.length_kpc
        scale_heights[time_label] = scale_height_kpc

        axis_edge.scatter(
            disk[:, 0] * units.length_kpc,
            disk[:, 2] * units.length_kpc,
            s=0.5, c=colour, alpha=0.4, linewidths=0.0,
            label=f"{time_label}  (z0 = {scale_height_kpc:.2f} kpc)",
        )

        absolute_height = np.abs(disk[:, 2]) * units.length_kpc
        counts, _ = np.histogram(absolute_height, bins=height_edges)
        normalised = counts / counts.max()
        axis_profile.semilogy(
            height_centres, normalised, color=colour,
            label=f"{time_label}  (z0 = {scale_height_kpc:.2f} kpc)",
        )

    axis_edge.set_xlim(-edge_zoom_kpc, edge_zoom_kpc)
    axis_edge.set_ylim(-height_zoom_kpc, height_zoom_kpc)
    axis_edge.set_aspect("equal")
    axis_edge.set_xlabel("x  [kpc]")
    axis_edge.set_ylabel("z  [kpc]")
    axis_edge.set_title("Disk edge-on (tight zoom)")
    axis_edge.legend(loc="upper right", fontsize=9, markerscale=8)

    axis_profile.set_xlabel("|z|  [kpc]")
    axis_profile.set_ylabel("normalised disk count")
    axis_profile.set_title("Vertical density profile")
    axis_profile.set_ylim(1.0e-3, 1.5)
    axis_profile.legend(fontsize=9)

    figure.suptitle(
        f"Disk vertical thickening -- N = {len(masses)}, "
        f"eps = {float(trajectory['softening']):g}  (MW physical scaling)",
        fontsize=13,
    )
    figure.savefig(output_png, dpi=120)
    plt.close(figure)
    return scale_heights


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

    thickening_png = output_png.replace(".png", "_thickening.png")
    scale_heights = disk_thickening_figure(trajectory, thickening_png)
    print(f"Saved {thickening_png}")
    for time_label, scale_height_kpc in scale_heights.items():
        print(f"  {time_label}: z0 = {scale_height_kpc:.3f} kpc")


if __name__ == "__main__":
    main()
