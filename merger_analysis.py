"""
Merger-remnant analysis for the section 3.0.5 exercise (Q7).

Compares the *stellar* (disk + bulge) structure before and after the
collision: the projected surface-density profile and the cumulative-mass
profile of the relaxed remnant versus the original exponential disks.
The key physical question is whether the violent relaxation of the
merger turns the cold exponential disks into a hot, centrally
concentrated spheroid that follows the de Vaucouleurs R^{1/4} law

    log10 Sigma(R)  is linear in  R^{1/4}

which is the observational signature of an elliptical galaxy.

Run as a script:
    python merger_analysis.py data/merger_N40k_eps01.npz

CMPH Project 3 -- Zoutong Shen / Zhaoyang Chu, 2026.
"""

import argparse
from typing import Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from disk_profiles import enclosed_mass_profile, surface_density_profile
from physical_units import milky_way_unit_system


def stellar_mask(
    num_disk: int, num_bulge: int, num_halo: int
) -> np.ndarray:
    """
    Boolean mask selecting the stellar (disk + bulge) particles of both
    galaxies in the stacked merger array.

    The merger IC stacks galaxy 1 then galaxy 2; within each galaxy the
    order is disk, bulge, halo. The halo is dark and excluded so the
    profiles compare with observable starlight.
    """
    per_galaxy = num_disk + num_bulge + num_halo
    mask = np.zeros(2 * per_galaxy, dtype=bool)
    for galaxy_start in (0, per_galaxy):
        mask[galaxy_start : galaxy_start + num_disk + num_bulge] = True
    return mask


def shrinking_sphere_centre(
    positions: np.ndarray,
    masses: np.ndarray,
    shrink_factor: float = 0.85,
    minimum_particles: int = 200,
) -> np.ndarray:
    """
    Locate the remnant centre with the shrinking-sphere method.

    Starts from the mass centroid and repeatedly recomputes it using only
    the particles inside a sphere that shrinks each iteration, which
    converges on the densest concentration and is robust to tidal debris
    and the surviving companion before full coalescence.
    """
    centre = (positions * masses[:, None]).sum(axis=0) / masses.sum()
    radius = np.linalg.norm(positions - centre, axis=1).max()
    while True:
        offsets = positions - centre
        inside = np.linalg.norm(offsets, axis=1) < radius
        if inside.sum() < minimum_particles:
            break
        centre = (
            (positions[inside] * masses[inside, None]).sum(axis=0)
            / masses[inside].sum()
        )
        radius *= shrink_factor
    return centre


def de_vaucouleurs_fit(
    bin_centres: np.ndarray, surface_density: np.ndarray
) -> Tuple[float, float, np.ndarray]:
    """
    Least-squares fit of log10(Sigma) against R^{1/4}.

    Returns:
        (slope, intercept, model) where ``model`` is the fitted
        log10(Sigma) on the input radii. A good straight-line fit
        (small residuals) means the remnant obeys the R^{1/4} law.
    """
    valid = surface_density > 0.0
    root_radius = bin_centres[valid] ** 0.25
    log_density = np.log10(surface_density[valid])
    slope, intercept = np.polyfit(root_radius, log_density, 1)
    model = slope * bin_centres ** 0.25 + intercept
    return slope, intercept, model


def merger_remnant_figure(
    trajectory: dict, output_png: str
) -> dict:
    """
    Plot the stellar surface density (vs R and vs R^{1/4}) and the
    cumulative mass, initial vs final, and return a numeric summary.
    """
    units = milky_way_unit_system()
    positions_history = trajectory["positions_history"]
    masses = trajectory["masses"]
    num_disk = int(trajectory["num_disk"])
    num_bulge = int(trajectory["num_bulge"])
    num_halo = int(trajectory["num_halo"])
    times = trajectory["snapshot_times"]

    stars = stellar_mask(num_disk, num_bulge, num_halo)
    star_masses = masses[stars]
    initial_positions = positions_history[0][stars].astype(np.float64)
    final_positions = positions_history[-1][stars].astype(np.float64)

    # Centre each epoch on its own stellar concentration: at t = 0 the two
    # galaxies are far apart, so use the global stellar centroid; for the
    # remnant use the shrinking-sphere density centre.
    initial_centre = (
        initial_positions * star_masses[:, None]
    ).sum(axis=0) / star_masses.sum()
    final_centre = shrinking_sphere_centre(final_positions, star_masses)
    initial_positions = initial_positions - initial_centre
    final_positions = final_positions - final_centre

    radius_edges = np.linspace(0.05, 12.0, 50)
    bin_centres = 0.5 * (radius_edges[:-1] + radius_edges[1:])
    bin_centres_kpc = bin_centres * units.length_kpc
    sigma_factor = units.surface_density_msun_pc2
    star_mass_per_particle = float(star_masses[0])

    _, sigma_initial = surface_density_profile(
        initial_positions, star_mass_per_particle, radius_edges
    )
    _, sigma_final = surface_density_profile(
        final_positions, star_mass_per_particle, radius_edges
    )

    spherical_edges = np.logspace(np.log10(0.05), np.log10(15.0), 60)
    _, mass_initial, _ = enclosed_mass_profile(
        initial_positions, star_masses, spherical_edges
    )
    _, mass_final, _ = enclosed_mass_profile(
        final_positions, star_masses, spherical_edges
    )
    radius_kpc = spherical_edges * units.length_kpc

    slope, intercept, model = de_vaucouleurs_fit(
        bin_centres_kpc, sigma_final * sigma_factor
    )

    figure, axes = plt.subplots(1, 3, figsize=(17, 5.0))

    axes[0].semilogy(bin_centres_kpc, sigma_initial * sigma_factor,
                     "C0-", label="progenitor disks (t=0)")
    axes[0].semilogy(bin_centres_kpc, sigma_final * sigma_factor,
                     "C3-", label=f"remnant (t={times[-1]:.0f})")
    axes[0].set_xlabel("R  [kpc]")
    axes[0].set_ylabel(r"$\Sigma$  [$M_\odot\,\mathrm{pc}^{-2}$]")
    axes[0].set_title("Stellar surface density vs R")
    axes[0].legend(fontsize=8)
    axes[0].set_xlim(0, 25)

    axes[1].semilogy(bin_centres_kpc ** 0.25, sigma_final * sigma_factor,
                     "C3o", markersize=3, label="remnant")
    axes[1].semilogy(bin_centres_kpc ** 0.25, 10.0 ** model,
                     "k--", label=f"R^1/4 fit (slope {slope:.2f})")
    axes[1].set_xlabel(r"$R^{1/4}$  [kpc$^{1/4}$]")
    axes[1].set_ylabel(r"$\Sigma$  [$M_\odot\,\mathrm{pc}^{-2}$]")
    axes[1].set_title("de Vaucouleurs test (straight = elliptical)")
    axes[1].legend(fontsize=8)

    axes[2].semilogy(radius_kpc, mass_initial * units.mass_msun,
                     "C0-", label="progenitor disks (t=0)")
    axes[2].semilogy(radius_kpc, mass_final * units.mass_msun,
                     "C3-", label=f"remnant (t={times[-1]:.0f})")
    axes[2].set_xlabel("r  [kpc]")
    axes[2].set_ylabel(r"$M_\star(<r)$  [$M_\odot$]")
    axes[2].set_title("Stellar cumulative mass")
    axes[2].legend(fontsize=8)
    axes[2].set_xlim(0, 30)

    figure.suptitle(
        f"Merger remnant -- N = {len(masses)}, "
        f"eps = {float(trajectory['softening']):g}  (MW physical scaling)"
    )
    figure.tight_layout()
    figure.savefig(output_png, dpi=130)
    plt.close(figure)

    half_mass_radius_final = float(
        np.interp(
            0.5 * mass_final[-1], mass_final, radius_kpc
        )
    )
    return {
        "de_vaucouleurs_slope": float(slope),
        "de_vaucouleurs_intercept": float(intercept),
        "remnant_half_mass_radius_kpc": half_mass_radius_final,
    }


def main() -> None:
    """Load a merger trajectory and produce the Q7 remnant figure."""
    parser = argparse.ArgumentParser(
        description="Merger-remnant profile analysis (Q7)."
    )
    parser.add_argument("trajectory", type=str)
    parser.add_argument("--out", type=str, default="")
    args = parser.parse_args()

    trajectory = np.load(args.trajectory)
    output_png = args.out or args.trajectory.replace(
        ".npz", "_remnant.png"
    ).replace("data/", "figures/")

    summary = merger_remnant_figure(trajectory, output_png)
    print(f"Saved {output_png}")
    print("\nNumeric summary:")
    for key, value in summary.items():
        print(f"  {key:34s} {value:12.4g}")


if __name__ == "__main__":
    main()
