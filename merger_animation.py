"""
Merger animation for the section 3.0.5 exercise / presentation.

Renders the two-galaxy collision as a movie: stellar particles coloured
by which galaxy they came from, shown face-on (x-y) and edge-on (x-z),
in a fixed frame centred on the system centre of mass (the IC is built
with zero net momentum, so the remnant stays near the origin). This
makes the approach, the tidal distortion at first passage, and the
relaxation into a single spheroid all visible in one clip.

The output is an mp4 if ffmpeg is available, otherwise an animated gif.
For high particle counts the plotted points and the frames are decimated
(``--max-points`` / ``--stride``) so the file stays a reasonable size;
this is a cosmetic plotting choice only and does not touch the data.

Run as a script:
    python merger_animation.py data/merger_N80k_eps01.npz
    python merger_animation.py data/merger_N160k_eps01.npz --stride 1 \\
        --max-points 40000 --out figures/merger_N160k.gif

CMPH Project 3 -- Zoutong Shen / Zhaoyang Chu, 2026.
"""

import argparse
import shutil

import matplotlib

matplotlib.use("Agg")
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

from merger_analysis import robust_stellar_mask, stellar_mask
from physical_units import milky_way_unit_system


def decimated_indices(
    count: int, maximum: int, rng: np.random.Generator
) -> np.ndarray:
    """
    Return up to ``maximum`` random indices into ``count`` items.

    Used purely to keep the rendered point cloud (and the output file
    size) manageable for large N; the simulation data is untouched.
    """
    if count <= maximum:
        return np.arange(count)
    return rng.choice(count, size=maximum, replace=False)


def merger_movie(
    trajectory: dict,
    output_path: str,
    *,
    frame_stride: int = 2,
    max_points_per_galaxy: int = 20000,
    view_kpc: float = 90.0,
    frames_per_second: int = 20,
) -> None:
    """
    Build and write the merger animation.

    Args:
        trajectory:            loaded merger .npz contents.
        output_path:           destination .mp4 or .gif path.
        frame_stride:          render every Nth snapshot, keyword-only.
        max_points_per_galaxy: cap on plotted stellar points per galaxy,
                               keyword-only.
        view_kpc:              half-width of the (square) field of view,
                               keyword-only.
        frames_per_second:     playback rate, keyword-only.
    """
    units = milky_way_unit_system()
    positions_history = trajectory["positions_history"]
    galaxy_id = trajectory["galaxy_id"]
    times = trajectory["snapshot_times"]
    num_disk = int(trajectory["num_disk"])
    num_bulge = int(trajectory["num_bulge"])
    num_halo = int(trajectory["num_halo"])

    stars = robust_stellar_mask(galaxy_id, num_disk, num_bulge, num_halo)
    star_galaxy = galaxy_id[stars]
    frame_indices = range(0, positions_history.shape[0], frame_stride)

    rng = np.random.default_rng(0)
    galaxy_colours = ["#5fa8ff", "#ff7a5f", "#6fdc8c", "#d18cff", "#ffd76f"]
    unique_galaxies = [int(g) for g in np.unique(star_galaxy)]
    show_per_galaxy = {
        gid: np.flatnonzero(star_galaxy == gid)[
            decimated_indices(
                int((star_galaxy == gid).sum()),
                max_points_per_galaxy, rng,
            )
        ]
        for gid in unique_galaxies
    }

    figure, (axis_face, axis_edge) = plt.subplots(1, 2, figsize=(15, 7.5))
    for axis, vertical_label in ((axis_face, "y"), (axis_edge, "z")):
        axis.set_xlim(-view_kpc, view_kpc)
        axis.set_ylim(-view_kpc, view_kpc)
        axis.set_aspect("equal")
        axis.set_xlabel("x  [kpc]")
        axis.set_ylabel(f"{vertical_label}  [kpc]")
        axis.set_facecolor("#0a0a0a")
    axis_face.set_title("Face-on (x-y)")
    axis_edge.set_title("Edge-on (x-z)")

    scatters = {}
    for axis, vertical in ((axis_face, 1), (axis_edge, 2)):
        for gid in unique_galaxies:
            scatters[(id(axis), gid)] = axis.scatter(
                [], [], s=0.6,
                c=galaxy_colours[gid % len(galaxy_colours)],
                alpha=0.5, linewidths=0.0,
            )
        scatters[(id(axis), "vertical")] = vertical
    title = figure.suptitle("")

    def update(frame_index: int):
        """Redraw both projections for one snapshot."""
        positions = positions_history[frame_index][stars].astype(
            np.float64
        ) * units.length_kpc
        for axis in (axis_face, axis_edge):
            vertical = scatters[(id(axis), "vertical")]
            for gid in unique_galaxies:
                indices = show_per_galaxy[gid]
                scatters[(id(axis), gid)].set_offsets(
                    np.column_stack([
                        positions[indices, 0],
                        positions[indices, vertical],
                    ])
                )
        gigayears = times[frame_index] * units.time_year / 1.0e9
        title.set_text(
            f"t = {times[frame_index]:.0f}  ({gigayears:.2f} Gyr)   "
            f"N = {len(galaxy_id)}"
        )
        return ()

    movie = animation.FuncAnimation(
        figure, update, frames=frame_indices, blit=False
    )

    if shutil.which("ffmpeg") is not None and output_path.endswith(".mp4"):
        writer = animation.FFMpegWriter(fps=frames_per_second, bitrate=4000)
    else:
        if not output_path.endswith(".gif"):
            output_path = output_path.rsplit(".", 1)[0] + ".gif"
        writer = animation.PillowWriter(fps=frames_per_second)

    movie.save(output_path, writer=writer, dpi=90)
    plt.close(figure)
    print(f"Saved {output_path}")


def main() -> None:
    """Load a merger trajectory and render the animation."""
    parser = argparse.ArgumentParser(
        description="Animate a two-galaxy merger trajectory."
    )
    parser.add_argument("trajectory", type=str)
    parser.add_argument("--out", type=str, default="")
    parser.add_argument("--stride", type=int, default=2)
    parser.add_argument("--max-points", type=int, default=20000)
    parser.add_argument("--view-kpc", type=float, default=90.0)
    parser.add_argument("--fps", type=int, default=20)
    args = parser.parse_args()

    trajectory = np.load(args.trajectory)
    output_path = args.out or args.trajectory.replace(
        ".npz", ".gif"
    ).replace("data/", "figures/")

    merger_movie(
        trajectory,
        output_path,
        frame_stride=args.stride,
        max_points_per_galaxy=args.max_points,
        view_kpc=args.view_kpc,
        frames_per_second=args.fps,
    )


if __name__ == "__main__":
    main()
