"""
Algorithm comparison: direct N^2 summation vs Barnes-Hut octree, SAME galaxy.

Unlike compare_rotation_curves.py (which compares the two *mass models*),
this script holds the galaxy fixed -- the manual-prescription equilibrium
initial conditions of the reference run -- and changes only the force
algorithm. Both kernels use the identical Plummer softening, so every
difference seen here is the tree's opening-angle approximation, nothing else.

Two modes:

  Static (default) -- one force evaluation on the t = 0 snapshot with both
  algorithms: per-particle relative force error |a_tree - a_N2| / |a_N2|
  at several opening angles theta, plus wall-clock timings. This isolates
  the per-step accuracy/speed trade and runs in under a minute.

  --evolve -- integrate the isolated galaxy to t = 100 with the TREE forces
  (same leapfrog, same dt = 0.125, same eps = 0.1 as the direct-N^2
  reference run), then compare the final disk surface-density profile and
  rotation curve against the stored direct-N^2 trajectory. This is the
  "same setup, different algorithm" cross-validation.

Run from the galaxy_merger_octree/ directory:
    python method_comparison.py            # static force test (fast)
    python method_comparison.py --evolve   # tree-driven evolution to t=100

CMPH Project 3 -- Zoutong Shen / Zhaoyang Chu, 2026.
"""

import argparse
import os
import sys
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# The direct-N^2 kernel and analysis helpers live in the parent project.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from forces import gravitational_accelerations
from nbody_simulation_vectorization import (
    build_flat_tree,
    compute_all_accelerations_parallel,
    get_bounds,
)

REFERENCE_TRAJECTORY = "../data/evolve_N20k_eps01.npz"
SCALING_TRAJECTORY = "../data/evolve_N80k_eps01.npz"
SOFTENING = 0.1
THETA_DEFAULT = 0.5
THETA_SCAN = (0.3, 0.5, 0.8)
LENGTH_KPC = 3.5
VELOCITY_KMS = 262.0


def tree_accelerations(positions, masses, softening, theta):
    """
    Accelerations from the flattened Barnes-Hut octree (one fresh build).

    Allocates the flat tree buffers, builds the tree for the given
    positions, and walks it for every particle -- exactly the per-step
    work of the tree integrator, so timing this call is a fair per-step
    cost comparison against one direct-N^2 evaluation.

    Args:
        positions: (N, 3) particle positions.
        masses:    (N,)   particle masses.
        softening: Plummer softening length.
        theta:     Barnes-Hut opening angle.

    Returns:
        (N, 3) array of accelerations.
    """
    num_particles = positions.shape[0]
    bounds = get_bounds(positions, margin=5.0)
    max_nodes = 2 * num_particles + 1000
    node_mass = np.zeros(max_nodes, dtype=np.float64)
    node_com = np.zeros((max_nodes, 3), dtype=np.float64)
    node_children = np.full((max_nodes, 8), -1, dtype=np.int32)
    node_bounds = np.zeros((max_nodes, 6), dtype=np.float64)
    node_leaf_idx = np.full(max_nodes, -1, dtype=np.int32)
    next_node_ptr = np.array([0], dtype=np.int32)
    indices_all = np.arange(num_particles, dtype=np.int32)
    build_flat_tree(positions, indices_all, bounds, masses,
                    node_mass, node_com, node_children, node_bounds,
                    node_leaf_idx, 0, 40, next_node_ptr)
    return compute_all_accelerations_parallel(
        positions, 0,
        node_mass, node_com, node_children, node_bounds, node_leaf_idx,
        softening ** 2, theta, 1.0)


def load_snapshot(npz_path, snapshot_index):
    """Return (positions, velocities, masses) of one stored snapshot."""
    data = np.load(npz_path)
    positions = data["positions_history"][snapshot_index].astype(np.float64)
    velocities = data["velocities_history"][snapshot_index].astype(np.float64)
    masses = data["masses"].astype(np.float64)
    return positions, velocities, masses


def static_force_test(output_png):
    """
    One-shot force comparison on the t = 0 equilibrium galaxy.

    Computes the direct-N^2 reference accelerations once, then the tree
    accelerations at each theta in THETA_SCAN, recording the per-particle
    relative error and the wall time of each evaluation (after a small
    warm-up call so Numba JIT compilation is not timed).
    """
    positions, _, masses = load_snapshot(REFERENCE_TRAJECTORY, 0)
    num_particles = positions.shape[0]
    print(f"Loaded t=0 snapshot: N = {num_particles}, eps = {SOFTENING}")

    # Warm up both JIT kernels on a small subset so compile time is excluded.
    gravitational_accelerations(positions[:2000], masses[:2000], SOFTENING)
    tree_accelerations(positions[:2000], masses[:2000], SOFTENING,
                       THETA_DEFAULT)

    start = time.perf_counter()
    reference = gravitational_accelerations(positions, masses, SOFTENING)
    direct_seconds = time.perf_counter() - start
    reference_norm = np.linalg.norm(reference, axis=1)
    print(f"direct N^2: {direct_seconds:.2f} s per force evaluation")

    errors_by_theta = {}
    tree_seconds_by_theta = {}
    for theta in THETA_SCAN:
        start = time.perf_counter()
        tree_acc = tree_accelerations(positions, masses, SOFTENING, theta)
        tree_seconds_by_theta[theta] = time.perf_counter() - start
        relative_error = (
            np.linalg.norm(tree_acc - reference, axis=1) / reference_norm
        )
        errors_by_theta[theta] = relative_error
        print(f"tree theta={theta}: "
              f"{tree_seconds_by_theta[theta]:.2f} s,  "
              f"median error {np.median(relative_error):.2e},  "
              f"90th pct {np.percentile(relative_error, 90):.2e},  "
              f"max {relative_error.max():.2e}")

    plt.style.use("dark_background")
    plt.rcParams.update({"font.size": 13})
    figure, (axis_error, axis_time) = plt.subplots(
        1, 2, figsize=(13.0, 5.2)
    )

    bins = np.logspace(-7, 0, 70)
    for theta, colour in zip(THETA_SCAN, ("#4da6ff", "#ffd166", "#ff7b4d")):
        axis_error.hist(
            errors_by_theta[theta], bins=bins, histtype="step", lw=2.2,
            color=colour,
            label=(rf"$\theta={theta}$  "
                   rf"(median {np.median(errors_by_theta[theta]):.1e})"),
        )
    axis_error.set_xscale("log")
    axis_error.set_xlabel(r"relative force error  $|a_\mathrm{tree}-a_{N^2}|\,/\,|a_{N^2}|$")
    axis_error.set_ylabel("particles per bin")
    axis_error.set_title(f"Force accuracy (same galaxy, N = {num_particles})")
    axis_error.legend(fontsize=11)

    # Scaling: time both algorithms at N = 80k as well, so the O(N^2) vs
    # O(N log N) growth is visible (16x more pairs vs ~5x more tree work).
    positions_big, _, masses_big = load_snapshot(SCALING_TRAJECTORY, 0)
    start = time.perf_counter()
    gravitational_accelerations(positions_big, masses_big, SOFTENING)
    direct_seconds_big = time.perf_counter() - start
    start = time.perf_counter()
    tree_accelerations(positions_big, masses_big, SOFTENING, THETA_DEFAULT)
    tree_seconds_big = time.perf_counter() - start
    print(f"N = {positions_big.shape[0]}: direct {direct_seconds_big:.2f} s, "
          f"tree theta={THETA_DEFAULT} {tree_seconds_big:.2f} s")

    group_labels = [f"N = {num_particles}", f"N = {positions_big.shape[0]}"]
    direct_times = [direct_seconds, direct_seconds_big]
    tree_times = [tree_seconds_by_theta[THETA_DEFAULT], tree_seconds_big]
    x_positions = np.arange(2)
    bars_direct = axis_time.bar(x_positions - 0.18, direct_times, 0.34,
                                color="#cccccc", label=r"direct $N^2$")
    bars_tree = axis_time.bar(x_positions + 0.18, tree_times, 0.34,
                              color="#4da6ff",
                              label=rf"tree $\theta={THETA_DEFAULT}$")
    axis_time.bar_label(bars_direct, fmt="%.2f s", fontsize=11)
    axis_time.bar_label(bars_tree, fmt="%.2f s", fontsize=11)
    axis_time.set_xticks(x_positions, group_labels)
    axis_time.set_ylabel("wall time per force evaluation  [s]")
    axis_time.set_title(r"Cost scaling: $O(N^2)$ vs $O(N\log N)$")
    axis_time.legend(fontsize=11)

    figure.suptitle(
        "Direct summation vs Barnes-Hut octree -- identical galaxy & softening"
    )
    figure.tight_layout()
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    figure.savefig(output_png, dpi=150)
    plt.close(figure)
    print(f"Saved {output_png}")


def evolve_with_tree(theta, timestep, t_stop, snapshot_times):
    """
    Integrate the reference initial conditions to t_stop with TREE forces.

    Uses the identical kick-drift-kick leapfrog and timestep as the
    direct-N^2 reference run, so the trajectory differs from the stored
    one only through the force algorithm.

    Returns:
        (snapshot_positions, masses): list of (N, 3) arrays at the
        requested snapshot times, and the particle masses.
    """
    positions, velocities, masses = load_snapshot(REFERENCE_TRAJECTORY, 0)
    num_steps = int(round(t_stop / timestep))
    save_steps = {int(round(t / timestep)) for t in snapshot_times}
    saved = {}

    accelerations = tree_accelerations(positions, masses, SOFTENING, theta)
    start = time.perf_counter()
    for step in range(1, num_steps + 1):
        velocities = velocities + 0.5 * timestep * accelerations
        positions = positions + timestep * velocities
        accelerations = tree_accelerations(positions, masses, SOFTENING,
                                           theta)
        velocities = velocities + 0.5 * timestep * accelerations
        if step in save_steps:
            saved[step * timestep] = positions.copy()
            elapsed = time.perf_counter() - start
            print(f"  t = {step * timestep:6.1f}  ({elapsed:.0f} s elapsed)")
    return saved, masses


def evolution_comparison(output_png, theta=THETA_DEFAULT,
                         timestep=0.125, t_stop=100.0):
    """
    Run the tree-driven evolution and compare t = 100 disk structure
    against the stored direct-N^2 trajectory: surface-density profile
    and rotation curve, same estimators for both.
    """
    from disk_profiles import surface_density_profile

    data = np.load(REFERENCE_TRAJECTORY)
    num_disk = int(data["num_disk"])
    final_direct = data["positions_history"][-1].astype(np.float64)
    masses = data["masses"].astype(np.float64)

    print(f"Evolving tree run: theta={theta}, dt={timestep}, "
          f"t_stop={t_stop} (this takes a while)...")
    saved, _ = evolve_with_tree(theta, timestep, t_stop, [t_stop])
    final_tree = saved[t_stop]

    radius_edges = np.linspace(0.05, 8.0, 45)
    bin_centres_kpc = 0.5 * (radius_edges[:-1] + radius_edges[1:]) * LENGTH_KPC
    disk_mass_per_particle = float(masses[0])

    def disk_sigma(positions):
        """Disk surface density about the disk's own centre of mass."""
        disk = positions[:num_disk]
        centred = disk - disk.mean(axis=0)
        _, sigma = surface_density_profile(
            centred, disk_mass_per_particle, radius_edges
        )
        return sigma

    def rotation_curve(positions):
        """v_c(r) from the cumulative mass profile, in km/s."""
        centred = positions - positions.mean(axis=0)
        radius = np.linalg.norm(centred, axis=1)
        order = np.argsort(radius)
        cumulative = np.cumsum(masses[order])
        r_grid = np.linspace(0.15, 10.0, 200)
        enclosed = np.interp(r_grid, radius[order], cumulative)
        return r_grid * LENGTH_KPC, np.sqrt(enclosed / r_grid) * VELOCITY_KMS

    plt.style.use("dark_background")
    plt.rcParams.update({"font.size": 13})
    figure, (axis_sigma, axis_vc) = plt.subplots(1, 2, figsize=(13.0, 5.2))

    axis_sigma.semilogy(bin_centres_kpc, disk_sigma(final_direct),
                        color="#cccccc", lw=2.5, label=r"direct $N^2$")
    axis_sigma.semilogy(bin_centres_kpc, disk_sigma(final_tree),
                        color="#4da6ff", lw=2.5, ls="--",
                        label=rf"tree $\theta={theta}$")
    axis_sigma.set_xlabel("R  [kpc]")
    axis_sigma.set_ylabel(r"disk $\Sigma(R)$  [code units]")
    axis_sigma.set_title(f"Disk surface density at t = {t_stop:.0f}")
    axis_sigma.legend()

    for positions, colour, ls, label in (
        (final_direct, "#cccccc", "-", r"direct $N^2$"),
        (final_tree, "#4da6ff", "--", rf"tree $\theta={theta}$"),
    ):
        r_kpc, vc_kms = rotation_curve(positions)
        axis_vc.plot(r_kpc, vc_kms, color=colour, ls=ls, lw=2.5, label=label)
    axis_vc.set_xlabel("r  [kpc]")
    axis_vc.set_ylabel(r"$v_c(r)$  [km/s]")
    axis_vc.set_title(f"Rotation curve at t = {t_stop:.0f}")
    axis_vc.legend()

    figure.suptitle(
        "Same galaxy, same leapfrog & softening -- only the force algorithm differs"
    )
    figure.tight_layout()
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    figure.savefig(output_png, dpi=150)
    plt.close(figure)
    print(f"Saved {output_png}")


def main():
    """Run the static force test, and optionally the tree-driven evolution."""
    parser = argparse.ArgumentParser(
        description="Direct N^2 vs Barnes-Hut octree on the same galaxy."
    )
    parser.add_argument("--evolve", action="store_true",
                        help="also run the tree-driven evolution to t=100")
    args = parser.parse_args()

    static_force_test("../figures/comparison/method_force_compare.png")
    if args.evolve:
        evolution_comparison("../figures/comparison/method_evolve_compare.png")


if __name__ == "__main__":
    main()
