"""
Reference Barnes-Hut N-body simulation built around a Python OctreeNode class.

This module is the readable, single-threaded counterpart of the flattened,
Numba-parallel implementation in nbody_simulation_vectorization.py. The tree
is built every step as a recursive Python object graph; per-particle force
accumulation traverses that graph and falls back to a direct kernel for
leaves. The integrator is a KDK leapfrog (the half-kick before the first
drift is folded into the per-step kernel, which performs drift then full
kick on the half-step velocities).

The Plummer-softened gravitational force used everywhere is
    a_i = G sum_j  m_j (r_j - r_i) / (|r_j - r_i|^2 + eps^2)^(3/2)
with softening eps = 0.05. The Barnes-Hut opening criterion compares the
x-axis side length of a node against the distance to its centre of mass.

Author: CMPH 2026 Project 3
Date: May 2026
"""

import numpy as np
from numba import njit
import time

# Physical and integration parameters (units: G = 1, length and time scales
# set by the initial conditions in data/galaxy_initial_state.npz).
G = 1.0            # gravitational constant
softening = 0.05   # Plummer softening length to avoid singular forces
theta = 0.5        # Barnes-Hut opening angle
dt = 0.01          # time step
n_steps = 1000     # number of simulation steps
save_interval = 50 # how often to save intermediate snapshots


# ---------------------------------------------------------------------------
# Python octree builder. Cheap to write and easy to inspect, but slow; the
# vectorised module replaces this with flat NumPy arrays.
# ---------------------------------------------------------------------------

class OctreeNode:
    """
    One node of a recursive Barnes-Hut octree.

    A node either holds children (eight slots, some possibly None) or is a
    leaf carrying a single particle index. The centre of mass and total
    mass below the node are precomputed during construction so the force
    walk can short-circuit on distant nodes via the opening-angle test.
    """

    def __init__(self, bounds):
        """
        Args:
            bounds: Triple ((xmin, xmax), (ymin, ymax), (zmin, zmax))
                    giving the cubical region this node covers.
        """
        self.bounds = bounds
        self.children = None          # list of 8 child nodes, or None on leaves
        self.center_of_mass = None    # (cx, cy, cz)
        self.total_mass = 0.0
        self.index = None             # particle index, only meaningful on leaves

    def is_leaf(self):
        """Return True when this node has no children."""
        return self.children is None


def build_tree(particles, indices, bounds, masses_all, depth=0, max_depth=40):
    """
    Build the octree of the particles within the given bounds.

    The centre of mass and total mass below each node are computed during
    the same descent that subdivides the region; this avoids a second pass
    afterwards. Recursion stops when a region contains a single particle or
    when max_depth has been reached, which guards against unbounded recursion
    if two particles share nearly identical coordinates.

    Args:
        particles: (N, 3) array of all particle positions in the simulation
        indices: Indices (into particles) of the particles that fall inside
                 the current node's bounds
        bounds: Triple ((xmin, xmax), (ymin, ymax), (zmin, zmax))
        masses_all: (N,) array of particle masses
        depth: Recursion depth of the current call
        max_depth: Maximum allowed recursion depth

    Returns:
        Root OctreeNode of the subtree covering bounds.
    """
    node = OctreeNode(bounds)
    n = len(indices)

    masses = masses_all[indices]
    total_mass = np.sum(masses)
    if total_mass > 0:
        positions = particles[indices]
        com = np.sum(positions * masses[:, np.newaxis], axis=0) / total_mass
        node.center_of_mass = com
        node.total_mass = total_mass
    else:
        node.center_of_mass = np.zeros(3)
        node.total_mass = 0.0

    if n == 1 or depth >= max_depth:
        node.index = indices[0] if n == 1 else None

        if n > 1 and depth >= max_depth:
            print(f"Warning: reached max depth with {n} particles in node.")
        node.index = indices[0]
        return node

    (xmin, xmax), (ymin, ymax), (zmin, zmax) = bounds
    xmid = 0.5 * (xmin + xmax)
    ymid = 0.5 * (ymin + ymax)
    zmid = 0.5 * (zmin + zmax)

    # Cartesian product of the three axes gives the eight child cubes.
    sub_bounds = []
    for xs in [(xmin, xmid), (xmid, xmax)]:
        for ys in [(ymin, ymid), (ymid, ymax)]:
            for zs in [(zmin, zmid), (zmid, zmax)]:
                sub_bounds.append((xs, ys, zs))

    # Assign each particle to one of the eight octants by midpoint compare.
    child_indices = [[] for _ in range(8)]
    pos = particles[indices]
    for i, idx in enumerate(indices):
        x, y, z = pos[i]
        child_idx = (0 if x < xmid else 1) + \
                    (0 if y < ymid else 2) + \
                    (0 if z < zmid else 4)
        child_indices[child_idx].append(idx)

    # Empty octants are left as None so the force walk can skip them quickly.
    node.children = []
    for i, cb in enumerate(sub_bounds):
        if len(child_indices[i]) > 0:
            child_node = build_tree(particles, child_indices[i], cb,
                                    masses_all, depth + 1, max_depth)
            node.children.append(child_node)
        else:
            node.children.append(None)
    return node


# ---------------------------------------------------------------------------
# Barnes-Hut force calculation.
# ---------------------------------------------------------------------------

@njit
def direct_acceleration_computation(target_pos, other_pos, other_mass,
                                    softening_sq, *, G=1.0):
    """
    Direct Plummer-softened acceleration on target from a list of sources.

    Used only at the leaves of the tree walk, where the source set is
    typically a single particle. The Numba decorator pays off because this
    is called N times per simulation step.

    Args:
        target_pos: (3,) array, position of the particle receiving force
        other_pos: (M, 3) array of source positions
        other_mass: (M,) array of source masses
        softening_sq: Plummer softening squared (epsilon^2)
        G: gravitational constant

    Returns:
        (3,) acceleration vector on the target.
    """
    ax = 0.0
    ay = 0.0
    az = 0.0
    for i in range(other_pos.shape[0]):
        dx = other_pos[i, 0] - target_pos[0]
        dy = other_pos[i, 1] - target_pos[1]
        dz = other_pos[i, 2] - target_pos[2]
        dist_sq = dx * dx + dy * dy + dz * dz + softening_sq
        dist_cube = dist_sq * np.sqrt(dist_sq)
        m = other_mass[i]
        prefactor = G * m / dist_cube
        ax += prefactor * dx
        ay += prefactor * dy
        az += prefactor * dz
    return np.array([ax, ay, az])


def compute_acceleration_octree(particles, node, masses_all, target_idx,
                                softening_sq):
    """
    Acceleration on one particle by recursively walking the octree.

    For each visited node, the opening criterion size / dist < theta is
    used to decide between approximating the node by its centre of mass
    and descending into the children. The x-axis side length is used as
    "size", consistent with the cube subdivision in build_tree.

    Args:
        particles: (N, 3) array of all positions
        node: Current octree node
        masses_all: (N,) array of particle masses
        target_idx: Index of the particle whose acceleration we compute
        softening_sq: Plummer softening squared

    Returns:
        (3,) acceleration vector on particle target_idx.
    """
    if node.total_mass == 0:
        return np.zeros(3)

    dx = node.center_of_mass[0] - particles[target_idx, 0]
    dy = node.center_of_mass[1] - particles[target_idx, 1]
    dz = node.center_of_mass[2] - particles[target_idx, 2]
    dist_sq = dx * dx + dy * dy + dz * dz + softening_sq
    dist = np.sqrt(dist_sq)

    if node.is_leaf():
        # No self-force: target's own leaf contributes zero.
        if node.index == target_idx:
            return np.zeros(3)
        other_pos = particles[node.index].reshape(1, 3)
        other_mass = np.array([masses_all[node.index]])
        return direct_acceleration_computation(particles[target_idx],
                                               other_pos, other_mass,
                                               softening_sq, G)

    else:
        # Opening-angle test. If the node looks small from this target,
        # approximate it by its centre of mass; otherwise descend.
        (xmin, xmax), _, _ = node.bounds
        size = xmax - xmin
        if size / dist < theta:
            prefactor = G * node.total_mass / (dist_sq * dist)
            return np.array([prefactor * dx, prefactor * dy, prefactor * dz])
        else:
            acc = np.zeros(3)
            for child in node.children:
                if child is not None:
                    acc += compute_acceleration_octree(particles, child,
                                                      masses_all, target_idx,
                                                      softening_sq)
            return acc


def compute_all_accelerations(particles, root, masses_all, softening_sq):
    """
    Per-particle acceleration via the Barnes-Hut tree.

    Each particle independently walks the same root tree, so this is the
    expensive loop the vectorised module parallelises with prange.

    Args:
        particles: (N, 3) positions
        root: Root OctreeNode
        masses_all: (N,) masses
        softening_sq: Plummer softening squared

    Returns:
        (N, 3) array of accelerations.
    """
    N = particles.shape[0]
    acc = np.zeros((N, 3))
    for i in range(N):
        acc[i] = compute_acceleration_octree(particles, root, masses_all, i,
                                             softening_sq)
    return acc


# ---------------------------------------------------------------------------
# Leapfrog integrator (drift-then-kick on half-step velocities).
# ---------------------------------------------------------------------------

@njit
def leapfrog_drift_kick(positions, velocities_half, accelerations, dt):
    """
    One leapfrog substep: drift positions then kick half-step velocities.

    With v_half initialised to v(0) + 0.5 dt a(0), repeated applications
    of this kernel advance positions by full steps and velocities by full
    steps offset by dt/2; this is the standard symplectic leapfrog.

    Args:
        positions: (N, 3) positions at time t
        velocities_half: (N, 3) velocities at t - dt/2 (after init: t + dt/2)
        accelerations: (N, 3) accelerations evaluated at the current position
        dt: time step

    Returns:
        Tuple (positions_new, velocities_half_new) at time t + dt.
    """
    positions_new = positions + velocities_half * dt
    velocities_half_new = velocities_half + accelerations * dt
    return positions_new, velocities_half_new


# ---------------------------------------------------------------------------
# Main driver.
# ---------------------------------------------------------------------------

def main():
    """
    Load the initial state, integrate for n_steps, and save snapshots.

    The octree bounds are computed once from the initial particle extent
    plus a 5-unit margin and held fixed for the whole run. Snapshots are
    written every save_interval steps, with velocities re-centred to the
    integer time using a fresh acceleration evaluation.
    """
    data = np.load("data/galaxy_initial_state.npz")
    positions = data["pos"].astype(np.float64)
    velocities = data["vel"].astype(np.float64)
    masses_all = data["mass"].astype(np.float64)
    labels = data["label"]

    N_particles = len(masses_all)
    print(f"Loaded {N_particles} particles from initial state, "
          f"total mass: {np.sum(masses_all)}.")

    # Pad the bounds beyond the initial extent so early-time excursions stay
    # inside the root cube.
    margin = 5.0
    xmin, ymin, zmin = np.min(positions, axis=0)
    xmax, ymax, zmax = np.max(positions, axis=0)
    bounds_global = ((xmin - margin, xmax + margin),
                     (ymin - margin, ymax + margin),
                     (zmin - margin, zmax + margin))

    softening_sq = softening ** 2

    # Initialise half-step velocities with one explicit half-kick so the
    # leapfrog kernel can run uniformly inside the loop.
    root = build_tree(positions, np.arange(N_particles), bounds_global,
                      masses_all)
    acc0 = compute_all_accelerations(positions, root, masses_all, softening_sq)
    velocities_half = velocities + 0.5 * acc0 * dt
    print("Initial accelerations computed and velocities at half step "
          "initialized.")

    start_time = time.time()
    for step in range(n_steps + 1):
        root = build_tree(positions, np.arange(N_particles), bounds_global,
                          masses_all)
        acc = compute_all_accelerations(positions, root, masses_all,
                                        softening_sq)
        positions, velocities_half = leapfrog_drift_kick(
            positions, velocities_half, acc, dt)

        if step % save_interval == 0:
            # Recover full-step velocities at the saved time by a half-kick
            # using the freshly evaluated acceleration at the new positions.
            root_new = build_tree(positions, np.arange(N_particles),
                                  bounds_global, masses_all)
            acc_new = compute_all_accelerations(positions, root_new,
                                                masses_all, softening_sq)
            velocities_full = velocities_half + 0.5 * acc_new * dt
            elapsed = time.time() - start_time
            print(f"Step {step}/{n_steps} - time elapsed: {elapsed:.2f} seconds")
            np.savez(f"data/galaxy_step_{step:04d}.npz",
                     pos=positions, vel=velocities_full,
                     mass=masses_all, label=labels)

        if step == n_steps:
            break

    print("Simulation completed, time elapsed: "
          f"{time.time() - start_time:.2f} seconds.")


if __name__ == "__main__":
    main()
