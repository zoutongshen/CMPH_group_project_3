"""
Isolated-disk stability test.

Loads (or builds) the equilibrium IC, evolves it in isolation, and reports
how much the disk's radial scale length h and vertical scale height z_0 have
drifted. A correctly-built IC should hold for at least 5 disk dynamical times
(the manual's tstop = 100, which is ~ 8 rotation periods at the half-mass
radius). Significant puffing-up means the velocity-distribution moments
do not match the density profile and the IC needs revisiting.

Diagnostics:
    h(t):        scale length, exponential fit to the disk surface density.
    z_0(t):      scale height, from <z^2>_disk = pi^2 z_0^2 / 6 (sech^2 dist).
    E(t), L_z(t): total energy and angular momentum (integrator quality).

Run as a script (parameters tunable from the CLI; reduce particle count for
quick sanity-check runs):
    python stability_test.py
    python stability_test.py --n-disk 3000 --n-bulge 1000 --n-halo 6000 --tstop 30

CMPH Project 3 -- Zoutong Shen / Zhaoyang Chu, 2026.
"""

import argparse
import time
from typing import Tuple

import numpy as np

from diagnostics import angular_momentum, total_energy
from initial_conditions import GalaxyParams
from integrator import run_simulation
from velocities import make_galaxy_initial_conditions


def fit_disk_scale_length(
    disk_positions: np.ndarray,
    disk_mass_per_particle: float,
    radial_min: float = 0.5,
    radial_max: float = 5.0,
    num_bins: int = 25,
) -> float:
    """
    Exponential-fit estimate of disk scale length h from sampled positions.

    Bins the surface density Sigma(R), fits ln Sigma vs R in the bulk, and
    returns h = -1/slope. The inner and outer bins are excluded: the inner
    region has the bulge contaminant and the outer is Poisson-noisy.

    Args:
        disk_positions:         shape (N_disk, 3).
        disk_mass_per_particle: mass of one disk particle.
        radial_min, radial_max: fit window in cylindrical radius.
        num_bins:               number of radial bins.

    Returns:
        Best-fit scale length h.
    """
    cylindrical_radius = np.sqrt(
        disk_positions[:, 0] ** 2 + disk_positions[:, 1] ** 2
    )
    bin_edges = np.linspace(0.0, radial_max + radial_min, num_bins + 1)
    bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    counts, _ = np.histogram(cylindrical_radius, bins=bin_edges)
    annulus_areas = np.pi * (bin_edges[1:] ** 2 - bin_edges[:-1] ** 2)
    surface_density = counts * disk_mass_per_particle / annulus_areas

    fit_mask = (
        (bin_centres > radial_min)
        & (bin_centres < radial_max)
        & (surface_density > 0.0)
    )
    radii_fit = bin_centres[fit_mask]
    log_density = np.log(surface_density[fit_mask])
    slope, _ = np.polyfit(radii_fit, log_density, 1)
    return -1.0 / slope


def fit_disk_scale_height(disk_positions: np.ndarray) -> float:
    """
    Sech^2 scale-height estimate from the disk particles' vertical RMS.

    For the normalised distribution rho(z) = (1/(2 z_0)) sech^2(z/z_0),
    direct integration gives <z^2> = pi^2 z_0^2 / 12, hence
    z_0 = sqrt(12 <z^2> / pi^2).
    """
    mean_z_squared = float((disk_positions[:, 2] ** 2).mean())
    return np.sqrt(12.0 * mean_z_squared / np.pi ** 2)


def measure_disk_structural_parameters(
    disk_positions: np.ndarray, disk_mass_per_particle: float
) -> Tuple[float, float]:
    """Return (scale_length_h, scale_height_z0) for the disk particles."""
    return (
        fit_disk_scale_length(disk_positions, disk_mass_per_particle),
        fit_disk_scale_height(disk_positions),
    )


def run_stability_test(
    params: GalaxyParams,
    *,
    timestep: float,
    tstop: float,
    softening: float,
    seed: int,
) -> dict:
    """
    Build a single galaxy IC and integrate in isolation. Return a dict of
    time series of diagnostics plus initial / final structural parameters.
    """
    positions, velocities, masses = make_galaxy_initial_conditions(
        params=params, seed=seed
    )
    num_steps = int(round(tstop / timestep))
    snapshot_interval = max(1, num_steps // 20)

    print(f"Building IC -- {positions.shape[0]} particles")
    print(
        f"Running integrator: dt = {timestep}, tstop = {tstop}, "
        f"steps = {num_steps}, snapshots = {num_steps // snapshot_interval + 1}"
    )

    start = time.time()
    positions_history, velocities_history, snapshot_times = run_simulation(
        positions,
        velocities,
        masses,
        timestep=timestep,
        num_steps=num_steps,
        softening=softening,
        snapshot_interval=snapshot_interval,
    )
    elapsed = time.time() - start
    print(f"Simulation finished in {elapsed:.1f} s")

    num_disk = params.num_disk
    disk_mass_per_particle = params.disk_mass / num_disk

    scale_lengths = np.empty(positions_history.shape[0])
    scale_heights = np.empty(positions_history.shape[0])
    energies = np.empty(positions_history.shape[0])
    angular_momentum_z = np.empty(positions_history.shape[0])

    for snapshot_index in range(positions_history.shape[0]):
        snap_positions = positions_history[snapshot_index]
        snap_velocities = velocities_history[snapshot_index]
        h_fit, z0_fit = measure_disk_structural_parameters(
            snap_positions[:num_disk], disk_mass_per_particle
        )
        scale_lengths[snapshot_index] = h_fit
        scale_heights[snapshot_index] = z0_fit
        energies[snapshot_index] = total_energy(
            snap_positions, snap_velocities, masses, softening
        )
        angular_momentum_z[snapshot_index] = angular_momentum(
            snap_positions, snap_velocities, masses
        )[2]

    return {
        "snapshot_times": snapshot_times,
        "scale_lengths": scale_lengths,
        "scale_heights": scale_heights,
        "energies": energies,
        "angular_momentum_z": angular_momentum_z,
        "h_initial": float(scale_lengths[0]),
        "h_final": float(scale_lengths[-1]),
        "z0_initial": float(scale_heights[0]),
        "z0_final": float(scale_heights[-1]),
    }


def report(results: dict, params: GalaxyParams) -> bool:
    """Pretty-print diagnostics and return True if the disk passed the test."""
    h_target = params.disk_scale_length
    z0_target = params.disk_scale_height
    h_drift = abs(results["h_final"] - results["h_initial"]) / results["h_initial"]
    z0_drift = (
        abs(results["z0_final"] - results["z0_initial"]) / results["z0_initial"]
    )
    energy_drift = (
        abs(results["energies"][-1] - results["energies"][0])
        / abs(results["energies"][0])
    )
    angular_initial = results["angular_momentum_z"][0]
    angular_drift = (
        abs(results["angular_momentum_z"][-1] - angular_initial)
        / abs(angular_initial)
    )

    h_threshold = 0.20
    z0_threshold = 0.40
    energy_threshold = 0.02
    angular_threshold = 0.05

    print()
    print(f"Disk structural drift")
    print(
        f"  h:  {results['h_initial']:.4f} -> {results['h_final']:.4f}  "
        f"(target {h_target:.2f}, drift {h_drift:.1%})"
    )
    print(
        f"  z0: {results['z0_initial']:.4f} -> {results['z0_final']:.4f}  "
        f"(target {z0_target:.2f}, drift {z0_drift:.1%})"
    )
    print(f"  E drift:   {energy_drift:.3e}")
    print(f"  L_z drift: {angular_drift:.3e}")

    h_passed = h_drift < h_threshold
    z0_passed = z0_drift < z0_threshold
    energy_passed = energy_drift < energy_threshold
    angular_passed = angular_drift < angular_threshold
    all_passed = h_passed and z0_passed and energy_passed and angular_passed

    print()
    print(
        f"  h drift          {'PASS' if h_passed else 'FAIL'}  "
        f"(tol {h_threshold:.0%})"
    )
    print(
        f"  z0 drift         {'PASS' if z0_passed else 'FAIL'}  "
        f"(tol {z0_threshold:.0%})"
    )
    print(
        f"  energy drift     {'PASS' if energy_passed else 'FAIL'}  "
        f"(tol {energy_threshold:.0e})"
    )
    print(
        f"  L_z drift        {'PASS' if angular_passed else 'FAIL'}  "
        f"(tol {angular_threshold:.0e})"
    )
    return all_passed


def parse_command_line() -> argparse.Namespace:
    """Parse CLI arguments for the stability test."""
    parser = argparse.ArgumentParser(
        description="Isolated-galaxy stability test for the equilibrium IC."
    )
    parser.add_argument("--n-disk", type=int, default=6000)
    parser.add_argument("--n-bulge", type=int, default=2000)
    parser.add_argument("--n-halo", type=int, default=12000)
    parser.add_argument("--tstop", type=float, default=100.0)
    parser.add_argument("--dt", type=float, default=0.125)
    parser.add_argument("--eps", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--save-snapshots",
        type=str,
        default="",
        help="if set, save snapshot arrays to this .npz file",
    )
    return parser.parse_args()


def main() -> None:
    """Top-level entry point."""
    args = parse_command_line()
    params = GalaxyParams(
        num_disk=args.n_disk,
        num_bulge=args.n_bulge,
        num_halo=args.n_halo,
    )
    results = run_stability_test(
        params,
        timestep=args.dt,
        tstop=args.tstop,
        softening=args.eps,
        seed=args.seed,
    )
    passed = report(results, params)

    if args.save_snapshots:
        np.savez(
            args.save_snapshots,
            snapshot_times=results["snapshot_times"],
            scale_lengths=results["scale_lengths"],
            scale_heights=results["scale_heights"],
            energies=results["energies"],
            angular_momentum_z=results["angular_momentum_z"],
        )

    if not passed:
        raise SystemExit("IC stability test FAILED -- see drift values above.")
    print("\nIC stability test PASSED.")


if __name__ == "__main__":
    main()
