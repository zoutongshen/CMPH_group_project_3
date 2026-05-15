"""
Quantitative profiles for the section 3.0.4 exercise (Q4, Q5, Q6).

Reproduces the analysis the manual's IDL program ``plotvcirc.pro`` does:
from a trajectory dump it measures, directly from the particles,

  * the cumulative (spherically enclosed) mass M(<r),
  * the circular velocity v_c(r) = sqrt(G M(<r) / r)  with G = 1,
  * the face-on projected surface density Sigma(R) of disk and bulge,

for the total system and per component, at the first (t = 0) and last
(t = 100) snapshot, and plots them in physical units -- km/s, solar
masses and solar masses per square parsec (Q4) -- so the t = 0 -> t = 100
change can be read off (Q5), and the eps = 1e-4 rerun compared (Q6).

Run as a script:
    python disk_profiles.py data/evolve_N20k_eps01.npz

CMPH Project 3 -- Zoutong Shen / Zhaoyang Chu, 2026.
"""

import argparse
from typing import Dict, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from physical_units import milky_way_unit_system


def split_components(
    positions: np.ndarray, num_disk: int, num_bulge: int
) -> Dict[str, np.ndarray]:
    """Slice a stacked particle array into disk / bulge / halo views."""
    end_bulge = num_disk + num_bulge
    return {
        "disk": positions[:num_disk],
        "bulge": positions[num_disk:end_bulge],
        "halo": positions[end_bulge:],
    }


def enclosed_mass_profile(
    positions: np.ndarray,
    masses: np.ndarray,
    radius_edges: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Spherically-enclosed mass M(<r) sampled at ``radius_edges``.

    Args:
        positions:    shape (N, 3), recentred on the mass centroid by caller.
        masses:       shape (N,).
        radius_edges: increasing radii at which to evaluate M(<r).

    Returns:
        (radius_edges, enclosed_mass, enclosed_count) where the i-th entry
        sums the mass / counts the particles with spherical radius
        <= radius_edges[i]. The count lets callers reject bins with too few
        particles to give a meaningful v_c.
    """
    spherical_radius = np.linalg.norm(positions, axis=1)
    order = np.argsort(spherical_radius)
    sorted_radius = spherical_radius[order]
    cumulative_mass = np.cumsum(masses[order])
    indices = np.searchsorted(sorted_radius, radius_edges, side="right") - 1
    enclosed_mass = np.where(indices >= 0, cumulative_mass[indices], 0.0)
    enclosed_count = np.where(indices >= 0, indices + 1, 0).astype(float)
    return radius_edges, enclosed_mass, enclosed_count


def circular_velocity_profile(
    positions: np.ndarray,
    masses: np.ndarray,
    radius_edges: np.ndarray,
    minimum_enclosed_count: int = 50,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Circular speed v_c(r) = sqrt(G M(<r) / r), G = 1, from the particles.

    Bins enclosing fewer than ``minimum_enclosed_count`` particles are set
    to NaN: at very small r only a handful of particles contribute and
    sqrt(M(<r)/r) is dominated by sampling noise, producing a spurious
    central spike that is not physical.
    """
    radius, enclosed, count = enclosed_mass_profile(
        positions, masses, radius_edges
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        circular_speed = np.sqrt(enclosed / radius)
    circular_speed = np.where(
        count >= minimum_enclosed_count, circular_speed, np.nan
    )
    return radius, circular_speed


def surface_density_profile(
    component_positions: np.ndarray,
    mass_per_particle: float,
    radius_edges: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Face-on projected surface density Sigma(R) for one component.

    Bins particles by cylindrical radius R = sqrt(x^2 + y^2) and divides
    the binned mass by the annulus area pi (R_out^2 - R_in^2).

    Returns:
        (bin_centres, surface_density).
    """
    cylindrical_radius = np.sqrt(
        component_positions[:, 0] ** 2 + component_positions[:, 1] ** 2
    )
    counts, _ = np.histogram(cylindrical_radius, bins=radius_edges)
    annulus_area = np.pi * (radius_edges[1:] ** 2 - radius_edges[:-1] ** 2)
    bin_centres = 0.5 * (radius_edges[:-1] + radius_edges[1:])
    surface_density = counts * mass_per_particle / annulus_area
    return bin_centres, surface_density


def _recentre(positions: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """Shift positions so the centre of mass sits at the origin."""
    centre = (positions * masses[:, None]).sum(axis=0) / masses.sum()
    return positions - centre


def central_surface_density(
    component_positions: np.ndarray,
    mass_per_particle: float,
    aperture_radius: float,
) -> float:
    """
    Mean surface density inside a fixed circular aperture.

    Robust central-Sigma estimator for the Q5/Q6 summary: counts all
    particles with cylindrical R < aperture_radius and divides by the
    aperture area, instead of reading a single noisy innermost bin.
    """
    cylindrical_radius = np.sqrt(
        component_positions[:, 0] ** 2 + component_positions[:, 1] ** 2
    )
    inside = int(np.count_nonzero(cylindrical_radius < aperture_radius))
    return inside * mass_per_particle / (np.pi * aperture_radius ** 2)


def analyse_snapshot(
    positions: np.ndarray,
    masses: np.ndarray,
    num_disk: int,
    num_bulge: int,
    spherical_edges: np.ndarray,
    cylindrical_edges: np.ndarray,
    central_aperture: float,
) -> Dict[str, object]:
    """
    Compute all profiles for one recentred snapshot.

    v_c and M(<r) use ``spherical_edges`` (log-spaced, with a min-count
    mask -- they are cumulative and otherwise noisy at small r); the
    surface densities use ``cylindrical_edges`` (linear -- equal-width
    annuli give a clean Sigma); the central-Sigma scalar uses a fixed
    aperture of radius ``central_aperture``.
    """
    positions = _recentre(positions, masses)
    end_bulge = num_disk + num_bulge
    component_slices = {
        "disk": slice(0, num_disk),
        "bulge": slice(num_disk, end_bulge),
        "halo": slice(end_bulge, None),
    }
    disk_mass_per_particle = float(masses[0])
    bulge_mass_per_particle = float(masses[num_disk])

    # Per-component v_c and M(<r): each component's circular-velocity
    # contribution uses *its own* enclosed mass, so that the totals satisfy
    # v_c,total^2 = sum_i v_c,i^2 and M_total = sum_i M_i (conventional
    # rotation-curve decomposition; this is what the manual's plotvcirc
    # produces for "the different components").
    circular_velocity = {"total": circular_velocity_profile(
        positions, masses, spherical_edges
    )}
    enclosed_mass = {"total": enclosed_mass_profile(
        positions, masses, spherical_edges
    )}
    for name, component in component_slices.items():
        circular_velocity[name] = circular_velocity_profile(
            positions[component], masses[component], spherical_edges,
            minimum_enclosed_count=20,
        )
        enclosed_mass[name] = enclosed_mass_profile(
            positions[component], masses[component], spherical_edges
        )

    return {
        "circular_velocity": circular_velocity,
        "enclosed_mass": enclosed_mass,
        "sigma_disk": surface_density_profile(
            positions[component_slices["disk"]],
            disk_mass_per_particle, cylindrical_edges,
        ),
        "sigma_bulge": surface_density_profile(
            positions[component_slices["bulge"]],
            bulge_mass_per_particle, cylindrical_edges,
        ),
        "central_sigma_disk": central_surface_density(
            positions[component_slices["disk"]],
            disk_mass_per_particle, central_aperture,
        ),
    }


def plot_comparison(
    trajectory: Dict[str, np.ndarray], output_png: str
) -> Dict[str, float]:
    """
    Plot t = 0 vs final profiles in physical units and return a numeric
    summary used to answer Q5 / Q6.
    """
    units = milky_way_unit_system()
    positions_history = trajectory["positions_history"]
    masses = trajectory["masses"]
    num_disk = int(trajectory["num_disk"])
    num_bulge = int(trajectory["num_bulge"])
    times = trajectory["snapshot_times"]

    spherical_edges = np.logspace(np.log10(0.05), np.log10(15.0), 60)
    cylindrical_edges = np.linspace(0.05, 12.0, 50)
    central_aperture = 1.0 / units.length_kpc  # 1 kpc, in code units

    initial = analyse_snapshot(
        positions_history[0].astype(np.float64), masses,
        num_disk, num_bulge, spherical_edges, cylindrical_edges,
        central_aperture,
    )
    final = analyse_snapshot(
        positions_history[-1].astype(np.float64), masses,
        num_disk, num_bulge, spherical_edges, cylindrical_edges,
        central_aperture,
    )

    radius_kpc = spherical_edges * units.length_kpc
    centres_kpc = (
        0.5 * (cylindrical_edges[:-1] + cylindrical_edges[1:])
        * units.length_kpc
    )
    component_colour = {
        "total": "k", "disk": "C0", "bulge": "C3", "halo": "C2",
    }

    figure, axes = plt.subplots(1, 3, figsize=(17, 5.0))

    # Panel 1: rotation-curve decomposition, t = 0 (solid) vs final (dashed).
    for name, colour in component_colour.items():
        axes[0].plot(
            radius_kpc, initial["circular_velocity"][name][1]
            * units.velocity_kms,
            colour + "-", label=f"{name} (t=0)",
        )
        axes[0].plot(
            radius_kpc, final["circular_velocity"][name][1]
            * units.velocity_kms,
            colour + "--", label=f"{name} (t={times[-1]:.0f})",
        )
    axes[0].set_xlabel("r  [kpc]")
    axes[0].set_ylabel(r"$v_c$  [km/s]")
    axes[0].set_title("Circular velocity (by component)")
    axes[0].legend(fontsize=7, ncol=2)
    axes[0].set_xlim(0, 30)
    axes[0].set_ylim(0, None)

    # Panel 2: enclosed mass, log y, per component + total.
    for name, colour in component_colour.items():
        axes[1].semilogy(
            radius_kpc, initial["enclosed_mass"][name][1] * units.mass_msun,
            colour + "-", label=f"{name} (t=0)",
        )
        axes[1].semilogy(
            radius_kpc, final["enclosed_mass"][name][1] * units.mass_msun,
            colour + "--", label=f"{name} (t={times[-1]:.0f})",
        )
    axes[1].set_xlabel("r  [kpc]")
    axes[1].set_ylabel(r"$M(<r)$  [$M_\odot$]")
    axes[1].set_title("Cumulative mass (by component)")
    axes[1].legend(fontsize=7, ncol=2)
    axes[1].set_xlim(0, 30)
    axes[1].set_ylim(1e7, 1e12)

    sigma_factor = units.surface_density_msun_pc2
    axes[2].semilogy(centres_kpc, initial["sigma_disk"][1] * sigma_factor,
                     "C0-", label="disk t = 0")
    axes[2].semilogy(centres_kpc, final["sigma_disk"][1] * sigma_factor,
                     "C3--", label=f"disk t = {times[-1]:.0f}")
    axes[2].semilogy(centres_kpc, initial["sigma_bulge"][1] * sigma_factor,
                     "C1-", label="bulge t = 0")
    axes[2].semilogy(centres_kpc, final["sigma_bulge"][1] * sigma_factor,
                     "C4--", label=f"bulge t = {times[-1]:.0f}")
    axes[2].set_xlabel("R  [kpc]")
    axes[2].set_ylabel(r"$\Sigma$  [$M_\odot\,\mathrm{pc}^{-2}$]")
    axes[2].set_title("Projected surface density")
    axes[2].legend(fontsize=8)
    axes[2].set_xlim(0, 20)
    axes[2].set_ylim(1e-1, 1e4)

    figure.suptitle(
        f"Isolated disk: N = {len(masses)}, "
        f"eps = {float(trajectory['softening']):g} "
        f"(physical units, MW scaling)"
    )
    figure.tight_layout()
    figure.savefig(output_png, dpi=130)
    plt.close(figure)

    peak_vc_initial = float(np.nanmax(initial["circular_velocity"]["total"][1]))
    peak_vc_final = float(np.nanmax(final["circular_velocity"]["total"][1]))
    return {
        "peak_vc_initial_kms": peak_vc_initial * units.velocity_kms,
        "peak_vc_final_kms": peak_vc_final * units.velocity_kms,
        "central_sigma_disk_initial":
            float(initial["central_sigma_disk"]) * sigma_factor,
        "central_sigma_disk_final":
            float(final["central_sigma_disk"]) * sigma_factor,
    }


def main() -> None:
    """Load a trajectory dump and produce the Q4/Q5 comparison figure."""
    parser = argparse.ArgumentParser(
        description="Profile analysis for the isolated-disk run (Q4/Q5/Q6)."
    )
    parser.add_argument("trajectory", type=str)
    parser.add_argument(
        "--out", type=str, default="",
        help="output PNG path; default derives from the trajectory name",
    )
    args = parser.parse_args()

    trajectory = np.load(args.trajectory)
    output_png = args.out or args.trajectory.replace(
        ".npz", "_profiles.png"
    ).replace("data/", "figures/")

    summary = plot_comparison(trajectory, output_png)
    print(f"Saved {output_png}")
    print("\nNumeric summary (physical units):")
    for key, value in summary.items():
        print(f"  {key:32s} {value:12.4g}")


if __name__ == "__main__":
    main()
