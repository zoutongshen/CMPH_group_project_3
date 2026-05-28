"""
Ring-formation diagnostic for the head-on z-axis collision.

After a perpendicular (face-on) galaxy-galaxy impact the struck disk
develops an outward-propagating density wave -- the "Cartwheel" ring.
This script isolates one galaxy's disk material in the head-on trajectory
and visualises (a) face-on scatter and (b) radial surface-density profile
at a handful of snapshots around the moment of penetration, so we can
tell whether a ring actually formed and how fast it propagated.

Run as a script:
    python ring_diagnostic.py
    python ring_diagnostic.py --times 18 22 25.5 30 36 45 --view-kpc 70

CMPH Project 3 -- extension/multi-galaxy (z-axis variant).
"""

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from physical_units import milky_way_unit_system


def disk_indices_for_galaxy(
    galaxy_id: np.ndarray, gid: int, num_disk: int,
) -> np.ndarray:
    """
    Return indices of the disk particles of galaxy ``gid``.

    Layout: for each galaxy the array is [disk | bulge | halo], so the
    first ``num_disk`` indices of a given galaxy's block are its disk.
    """
    indices = np.flatnonzero(galaxy_id == gid)
    return indices[:num_disk]


def radial_surface_density(
    radii_kpc: np.ndarray, mass_per_particle: float, edges_kpc: np.ndarray,
) -> np.ndarray:
    """Surface density Sigma(R) [M_sun kpc^-2] in concentric annuli."""
    counts, _ = np.histogram(radii_kpc, bins=edges_kpc)
    inner = edges_kpc[:-1]
    outer = edges_kpc[1:]
    annulus_area = np.pi * (outer ** 2 - inner ** 2)
    return counts * mass_per_particle / annulus_area


def nearest_snapshot(snapshot_times: np.ndarray, target_t: float) -> int:
    return int(np.argmin(np.abs(snapshot_times - target_t)))


def build_figure(
    trajectory_path: str,
    output_path: str,
    target_times: list,
    view_kpc: float,
    profile_max_kpc: float,
    num_bins: int,
) -> None:
    """Render the two-row diagnostic figure."""
    units = milky_way_unit_system()
    trajectory = np.load(trajectory_path)

    positions_history = trajectory["positions_history"]
    galaxy_id = trajectory["galaxy_id"]
    masses = trajectory["masses"]
    snapshot_times = trajectory["snapshot_times"]
    num_disk = int(trajectory["num_disk"])

    galaxy_1_disk = disk_indices_for_galaxy(galaxy_id, 0, num_disk)
    galaxy_2_disk = disk_indices_for_galaxy(galaxy_id, 1, num_disk)
    disk_indices = galaxy_1_disk

    # All disk particles share the same mass in code units; convert once.
    mass_per_particle_code = float(masses[disk_indices[0]])
    mass_per_particle_msun = mass_per_particle_code * units.mass_msun

    edges_kpc = np.linspace(0.0, profile_max_kpc, num_bins + 1)
    bin_centres = 0.5 * (edges_kpc[:-1] + edges_kpc[1:])

    # Reference profile -- the initial disk at t = 0 -- so the user can
    # see any bump (ring) clearly above the smooth exponential baseline.
    positions_initial = positions_history[0][galaxy_1_disk].astype(np.float64)
    centre_initial = positions_initial.mean(axis=0)
    positions_initial_kpc = (positions_initial - centre_initial) * units.length_kpc
    radii_initial = np.sqrt(
        positions_initial_kpc[:, 0] ** 2 + positions_initial_kpc[:, 1] ** 2
    )
    sigma_initial = radial_surface_density(
        radii_initial, mass_per_particle_msun, edges_kpc,
    )
    sigma_initial_plot = np.where(sigma_initial > 0, sigma_initial, np.nan)

    num_panels = len(target_times)
    figure, axes = plt.subplots(
        2, num_panels, figsize=(3.6 * num_panels, 7.5), squeeze=False,
    )
    figure.patch.set_facecolor("#0a0a0a")

    actual_times = []
    for col, target_t in enumerate(target_times):
        snap = nearest_snapshot(snapshot_times, target_t)
        actual_t = float(snapshot_times[snap])
        actual_times.append(actual_t)
        positions_code = positions_history[snap]

        # Recentre on galaxy 1's disk centre of mass so the ring is
        # measured in the disk's own frame even after it's drifted.
        positions_g1 = positions_code[galaxy_1_disk].astype(np.float64)
        centre_g1 = positions_g1.mean(axis=0)
        positions_g1_kpc = (positions_g1 - centre_g1) * units.length_kpc

        # Top row -- face-on scatter (looking down z, which is the
        # collision axis; the disk lives in the x-y plane).
        axis_scatter = axes[0, col]
        axis_scatter.scatter(
            positions_g1_kpc[:, 0], positions_g1_kpc[:, 1],
            s=0.6, c="#5fa8ff", alpha=0.55, linewidths=0.0,
        )
        # Show galaxy 2's disk faintly so we can see overlap during impact.
        positions_g2_kpc = (
            positions_code[galaxy_2_disk].astype(np.float64) - centre_g1
        ) * units.length_kpc
        axis_scatter.scatter(
            positions_g2_kpc[:, 0], positions_g2_kpc[:, 1],
            s=0.4, c="#ff7a5f", alpha=0.20, linewidths=0.0,
        )
        axis_scatter.set_xlim(-view_kpc, view_kpc)
        axis_scatter.set_ylim(-view_kpc, view_kpc)
        axis_scatter.set_aspect("equal")
        axis_scatter.set_facecolor("#0a0a0a")
        axis_scatter.tick_params(colors="#888888", labelsize=7)
        for spine in axis_scatter.spines.values():
            spine.set_color("#444444")
        gigayears = actual_t * units.time_year / 1.0e9
        axis_scatter.set_title(
            f"t = {actual_t:.1f}  ({gigayears:.2f} Gyr)",
            color="#ffffff", fontsize=10,
        )

        # Bottom row -- radial Sigma(R) of galaxy 1's disk.
        radii_kpc = np.sqrt(
            positions_g1_kpc[:, 0] ** 2 + positions_g1_kpc[:, 1] ** 2
        )
        sigma = radial_surface_density(
            radii_kpc, mass_per_particle_msun, edges_kpc,
        )
        # Replace zero bins with NaN so log axis doesn't go to -inf.
        sigma_plot = np.where(sigma > 0, sigma, np.nan)
        axis_profile = axes[1, col]
        axis_profile.semilogy(
            bin_centres, sigma_initial_plot,
            color="#888888", lw=1.0, ls="--", label="initial disk",
        )
        axis_profile.semilogy(
            bin_centres, sigma_plot,
            color="#5fa8ff", lw=1.6, label="current",
        )
        axis_profile.set_xlim(0.0, profile_max_kpc)
        axis_profile.set_ylim(1.0e5, 1.0e9)
        axis_profile.set_facecolor("#0a0a0a")
        axis_profile.grid(True, which="both", color="#222222", lw=0.4)
        axis_profile.tick_params(colors="#888888", labelsize=7)
        for spine in axis_profile.spines.values():
            spine.set_color("#444444")
        if col == 0:
            legend = axis_profile.legend(
                fontsize=7, loc="upper right", frameon=False,
            )
            for text in legend.get_texts():
                text.set_color("#cccccc")

        # Locate any density bump beyond the (very steep) inner disk.
        # A ring shows up as a *local* maximum of Sigma well outside
        # R ~ 0, and ideally one that sits above the initial profile.
        if np.isfinite(sigma_plot).any():
            inner_skip = max(3, num_bins // 12)
            outer_part = sigma_plot[inner_skip:]
            if np.isfinite(outer_part).any():
                local_peak = int(np.nanargmax(outer_part)) + inner_skip
                peak_radius = bin_centres[local_peak]
                peak_value = sigma_plot[local_peak]
                # Only flag as a candidate ring if the bump is at least
                # 1.5x the initial profile at the same radius.
                ref_value = sigma_initial_plot[local_peak]
                is_ring_bump = (
                    np.isfinite(ref_value) and peak_value > 1.5 * ref_value
                )
                colour = "#ffd76f" if is_ring_bump else "#666666"
                axis_profile.axvline(
                    peak_radius, color=colour, lw=0.9, ls="--", alpha=0.85,
                )
                label = (
                    f"  ring? R = {peak_radius:.1f} kpc"
                    if is_ring_bump
                    else f"  peak R = {peak_radius:.1f} kpc"
                )
                axis_profile.text(
                    peak_radius, 3.0e8, label,
                    color=colour, fontsize=7, va="center",
                )

    axes[0, 0].set_ylabel(
        "Galaxy 1 disk (xy)\ny  [kpc]", color="#cccccc", fontsize=9,
    )
    axes[1, 0].set_ylabel(
        r"$\Sigma(R)$  [M$_\odot$ / kpc$^2$]", color="#cccccc", fontsize=9,
    )
    for col in range(num_panels):
        axes[1, col].set_xlabel("R  [kpc]", color="#cccccc", fontsize=8)
    figure.suptitle(
        "Ring diagnostic -- head-on z-axis collision (galaxy 1 disk shown,"
        " galaxy 2 disk faded)",
        color="#ffffff", fontsize=11,
    )

    plt.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    figure.savefig(output_path, dpi=120, facecolor=figure.get_facecolor())
    plt.close(figure)
    print(f"Saved {output_path}")
    print("Snapshot times used:", ", ".join(f"{t:.2f}" for t in actual_times))


def main() -> None:
    """Parse CLI args and render."""
    parser = argparse.ArgumentParser(
        description="Cartwheel-ring diagnostic for the head-on z-axis run."
    )
    parser.add_argument(
        "--trajectory", type=str,
        default="data/multi_zaxis_headon.npz",
    )
    parser.add_argument(
        "--times", type=float, nargs="+",
        default=[18.0, 22.5, 25.5, 30.0, 36.0, 45.0],
    )
    parser.add_argument("--view-kpc", type=float, default=70.0)
    parser.add_argument("--profile-max-kpc", type=float, default=60.0)
    parser.add_argument("--num-bins", type=int, default=40)
    parser.add_argument(
        "--out", type=str,
        default="figures/z_axis/multi_zaxis_headon_ring_diagnostic.png",
    )
    args = parser.parse_args()

    build_figure(
        args.trajectory, args.out, args.times,
        args.view_kpc, args.profile_max_kpc, args.num_bins,
    )


if __name__ == "__main__":
    main()
