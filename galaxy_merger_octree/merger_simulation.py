"""
N-body simulation of a two-galaxy merger, with embedded pre-equilibration.

merger_initialization.py produces a COLD merger IC: bulge and halo have
zero internal velocity dispersion. Running that straight through the
merger would let the bulge and halo collapse during the first few time
units and heat the disk before pericentre, suppressing the tidal
response. This driver therefore equilibrates each galaxy in isolation
first: it splits the loaded state by label, subtracts each galaxy's bulk
position and velocity, evolves it with the same Barnes-Hut leapfrog used
for the merger itself, and adds the bulk motion back. Both galaxies then
end up at their original positions with their original bulk velocities
but with the bulge and halo virialised around the rotating disk, and the
actual merger run proceeds from there.

The full pipeline is just two commands:
    python merger_initialization.py    # cold IC -> data/merger_initial_state.npz
    python merger_simulation.py        # equilibrate + merger -> data/snapshot_merger/

Softening is 0.1 for both the equilibration and merger phases, so the
bulge settles into a configuration consistent with the force kernel that
will be used during the encounter itself.

Author: CMPH 2026 Project 3
Date: May 2026
"""

import os
import time

import numpy as np

import nbody_simulation_vectorization as nb


# Simulation parameters. n_equil_steps is per-galaxy; the equilibration
# phase therefore runs the leapfrog twice for n_equil_steps each (once per
# galaxy, sequentially) before the merger phase begins.
softening = 0.1
theta = nb.theta
G = nb.G
dt = 0.01
n_equil_steps = 1500           # per-galaxy equilibration (t = 15)
n_steps = 10000                # merger run (t = 100); covers first passage,
                               # clear separation, tail formation, and merger
save_interval = 10
output_dir = "data/snapshot_merger"

# Tree root-cube padding. Equilibration uses a tighter cube since each
# galaxy is at the origin; the merger uses a larger cube so the galaxies
# (and any unbound debris / tidal tails) stay inside the root cube for
# the full run -- raise this if particles escape during long encounters.
bounds_margin_eq = 20.0
bounds_margin_merger = 60.0


def equilibrate_one_galaxy(positions, velocities, masses, name):
    """
    Equilibrate one galaxy in isolation, preserving its bulk motion.

    The galaxy's centre-of-mass position and velocity are subtracted
    before integrating and re-added at the end, so the galaxy comes out
    at its original location with its original bulk velocity but with
    the bulge and halo virialised around the rotating disk. The same
    Barnes-Hut leapfrog used by the merger phase is applied here, just
    to the single-galaxy subset.

    Args:
        positions: (N, 3) positions of one galaxy's particles in the
            simulation (lab) frame.
        velocities: (N, 3) velocities in the lab frame, including the
            galaxy's bulk velocity.
        masses: (N,) per-particle masses.
        name: Short label used only for the progress prints.

    Returns:
        Tuple (pos_out, vel_out) of equilibrated phase-space arrays in
        the same lab frame as the inputs.
    """
    print(f"Equilibrating {name}...")
    N = len(masses)

    # Move into the galaxy's rest frame so the bulk velocity does not
    # carry the galaxy out of the (small) equilibration root cube.
    total_mass = np.sum(masses)
    com_pos = np.sum(positions * masses[:, None], axis=0) / total_mass
    com_vel = np.sum(velocities * masses[:, None], axis=0) / total_mass
    pos = positions - com_pos
    vel = velocities - com_vel

    bounds = nb.get_bounds(pos, margin=bounds_margin_eq)
    softening_sq = softening ** 2

    # Tree buffers sized for one galaxy.
    max_nodes = N * 2 + 1000
    node_mass = np.zeros(max_nodes, dtype=np.float64)
    node_com = np.zeros((max_nodes, 3), dtype=np.float64)
    node_children = np.full((max_nodes, 8), -1, dtype=np.int32)
    node_bounds = np.zeros((max_nodes, 6), dtype=np.float64)
    node_leaf_idx = np.full(max_nodes, -1, dtype=np.int32)

    indices = np.arange(N, dtype=np.int32)
    nptr = np.array([0], dtype=np.int32)
    nb.build_flat_tree(pos, indices, bounds, masses,
                       node_mass, node_com, node_children, node_bounds,
                       node_leaf_idx, 0, 40, nptr)

    # Initial half-kick so the per-step update is a uniform full kick.
    acc = nb.compute_all_accelerations_parallel(
        pos, 0, node_mass, node_com, node_children, node_bounds,
        node_leaf_idx, softening_sq, theta, G)
    vel_half = vel + 0.5 * dt * acc

    start = time.time()
    for step in range(1, n_equil_steps + 1):
        pos = pos + vel_half * dt

        nptr[0] = 0
        nb.build_flat_tree(pos, indices, bounds, masses,
                           node_mass, node_com, node_children, node_bounds,
                           node_leaf_idx, 0, 40, nptr)
        acc = nb.compute_all_accelerations_parallel(
            pos, 0, node_mass, node_com, node_children, node_bounds,
            node_leaf_idx, softening_sq, theta, G)
        vel_half = vel_half + acc * dt

        if step % 100 == 0:
            elapsed = time.time() - start
            print(f"  [{name}] step {step:4d}/{n_equil_steps}"
                  f"  |  t = {step * dt:.2f}  |  elapsed = {elapsed:.1f}s")

    # Recover integer-time velocities in the rest frame, then add the
    # bulk motion back to return to the original lab frame.
    vel_rest = vel_half - 0.5 * dt * acc
    pos_out = pos + com_pos
    vel_out = vel_rest + com_vel
    return pos_out, vel_out


def main():
    """Load cold IC, equilibrate each galaxy, then run the merger."""
    os.makedirs(output_dir, exist_ok=True)
    data = np.load("data/merger_initial_state.npz")
    positions = data["pos"].astype(np.float64)
    velocities = data["vel"].astype(np.float64)
    masses = data["mass"].astype(np.float64)
    labels = data["label"]
    N = len(masses)
    print(f"Loaded {N} particles (cold merger IC), total mass = {np.sum(masses):.2f}")

    # ------------------------------------------------------------------
    # Equilibration phase: isolated leapfrog on each galaxy in turn.
    # ------------------------------------------------------------------
    mask1 = labels < 3            # labels 0/1/2 = galaxy 1
    mask2 = labels >= 3           # labels 3/4/5 = galaxy 2

    pos1, vel1 = equilibrate_one_galaxy(
        positions[mask1].copy(), velocities[mask1].copy(),
        masses[mask1], name="galaxy 1")
    pos2, vel2 = equilibrate_one_galaxy(
        positions[mask2].copy(), velocities[mask2].copy(),
        masses[mask2], name="galaxy 2")

    positions[mask1] = pos1
    velocities[mask1] = vel1
    positions[mask2] = pos2
    velocities[mask2] = vel2
    print("Equilibration done. Starting merger run.")

    # ------------------------------------------------------------------
    # Merger phase: full two-galaxy leapfrog on the equilibrated state.
    # ------------------------------------------------------------------
    bounds_global = nb.get_bounds(positions, margin=bounds_margin_merger)
    softening_sq = softening ** 2

    max_nodes = N * 2 + 1000
    node_mass = np.zeros(max_nodes, dtype=np.float64)
    node_com = np.zeros((max_nodes, 3), dtype=np.float64)
    node_children = np.full((max_nodes, 8), -1, dtype=np.int32)
    node_bounds = np.zeros((max_nodes, 6), dtype=np.float64)
    node_leaf_idx = np.full(max_nodes, -1, dtype=np.int32)

    indices_all = np.arange(N, dtype=np.int32)
    next_node_ptr = np.array([0], dtype=np.int32)
    nb.build_flat_tree(positions, indices_all, bounds_global, masses,
                       node_mass, node_com, node_children, node_bounds,
                       node_leaf_idx, 0, 40, next_node_ptr)
    root = 0
    print(f"Merger tree built, used {next_node_ptr[0]} nodes.")

    # Initial half-kick for the merger phase.
    acc0 = nb.compute_all_accelerations_parallel(
        positions, root,
        node_mass, node_com, node_children, node_bounds, node_leaf_idx,
        softening_sq, theta, G)
    velocities_half = velocities + 0.5 * dt * acc0

    # Step-0 snapshot stores the equilibrated, pre-merger state.
    np.savez(os.path.join(output_dir, "galaxy_step_0000.npz"),
             pos=positions, vel=velocities, mass=masses, label=labels)
    print(f"Step    0/{n_steps}  |  time =  0.00  |  elapsed =  0.0s")

    start_time = time.time()
    for step in range(1, n_steps + 1):
        positions = positions + velocities_half * dt

        next_node_ptr[0] = 0
        nb.build_flat_tree(positions, indices_all, bounds_global, masses,
                           node_mass, node_com, node_children, node_bounds,
                           node_leaf_idx, 0, 40, next_node_ptr)
        acc = nb.compute_all_accelerations_parallel(
            positions, root,
            node_mass, node_com, node_children, node_bounds, node_leaf_idx,
            softening_sq, theta, G)
        velocities_half = velocities_half + acc * dt

        if step % save_interval == 0:
            velocities_full = velocities_half - 0.5 * dt * acc
            elapsed = time.time() - start_time
            print(f"Step {step:4d}/{n_steps}  |  time = {step * dt:.2f}"
                  f"  |  elapsed = {elapsed:.1f}s")
            np.savez(os.path.join(output_dir, f"galaxy_step_{step:04d}.npz"),
                     pos=positions, vel=velocities_full,
                     mass=masses, label=labels)

    print(f"Simulation finished in {time.time() - start_time:.1f} seconds.")


if __name__ == "__main__":
    main()
