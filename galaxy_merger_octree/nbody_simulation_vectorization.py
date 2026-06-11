"""
Flattened, Numba-parallel Barnes-Hut N-body simulation.

The tree is stored as a struct-of-arrays in plain NumPy buffers
(node_mass, node_com, node_children, node_bounds, node_leaf_idx). This
representation is what allows the per-particle force walk to live inside
an @njit(parallel=True) kernel: prange threads can iterate over particles
without the GIL because they only read shared array data.

Two leaf flavours are encoded in node_leaf_idx:
    >= 0 : single-particle leaf, value is the particle index
     -1  : internal (non-leaf) node
     -2  : "pseudo-leaf" reached at max_depth with more than one particle;
           force from this node is taken from its centre of mass.

The integrator is a leapfrog in the symplectic (drift-then-full-kick) form,
with an explicit initial half-kick used so the in-loop kernel does not
need a special first step. The Plummer-softened gravitational kernel and
the opening criterion mirror nbody_simulation_class.py exactly so the two
implementations can be cross-checked.

Output snapshots land in data/snapshot_single/galaxy_step_NNNN.npz with
positions, integer-time velocities, masses and species labels.

Author: CMPH 2026 Project 3
Date: May 2026
"""

import numpy as np
from numba import njit, prange
import time
import os

# Physical and integration parameters (G = 1; length/time units set by IC).
G = 1.0
softening = 0.05
theta = 0.5
dt = 0.01
n_steps = 1000
save_interval = 10
output_dir = "data/snapshot_single"


# ---------------------------------------------------------------------------
# Helper: derive a padded root cube from the current particle extent.
# ---------------------------------------------------------------------------

def get_bounds(particles, margin=5.0):
    """
    Axis-aligned bounding box of the particle cloud, padded by margin.

    Used once at startup to seed the root cube of the octree. The margin
    leaves room for early-time excursions before any particle escapes the
    root region.

    Args:
        particles: (N, 3) array of particle positions
        margin: Extra padding added to every face of the bounding box

    Returns:
        (3, 2) array [[xmin, xmax], [ymin, ymax], [zmin, zmax]].
    """
    xmin, ymin, zmin = np.min(particles, axis=0)
    xmax, ymax, zmax = np.max(particles, axis=0)
    return np.array([[xmin - margin, xmax + margin],
                     [ymin - margin, ymax + margin],
                     [zmin - margin, zmax + margin]], dtype=np.float64)


# ---------------------------------------------------------------------------
# Numba-compiled flattened octree builder.
# ---------------------------------------------------------------------------

@njit
def build_flat_tree(particles, indices, bounds, masses,
                    node_mass, node_com, node_children, node_bounds,
                    node_leaf_idx, depth, max_depth, next_node_ptr):
    """
    Recursively populate the flat octree arrays.

    Each call writes one node's fields, then recurses into the (up to
    eight) non-empty octants. Node allocation is a bump-pointer:
    next_node_ptr is a length-1 int array used as a mutable counter so
    every recursive call sees and advances the same value.

    Args:
        particles: (N, 3) all particle positions
        indices: int32 array of particle indices currently belonging to this
                 region
        bounds: (3, 2) bounds of this region (xmin/xmax, ymin/ymax, zmin/zmax)
        masses: (N,) particle masses
        node_mass, node_com, node_children, node_bounds, node_leaf_idx:
            Preallocated tree buffers; this call fills the row at
            next_node_ptr[0] and recursively fills its descendants
        depth: Recursion depth of the current call
        max_depth: Hard cap; nodes deeper than this become pseudo-leaves
        next_node_ptr: Length-1 int array acting as a mutable next-free index

    Returns:
        next_node_ptr (same object), with [0] advanced past every node
        written by this call and its descendants.
    """
    n = len(indices)
    node_idx = next_node_ptr[0]
    next_node_ptr[0] += 1

    # Centre of mass and total mass, accumulated in a single scan.
    total_mass = 0.0
    cx = 0.0
    cy = 0.0
    cz = 0.0
    for i in range(n):
        p_idx = indices[i]
        m = masses[p_idx]
        total_mass += m
        cx += m * particles[p_idx, 0]
        cy += m * particles[p_idx, 1]
        cz += m * particles[p_idx, 2]
    if total_mass > 0:
        cx /= total_mass
        cy /= total_mass
        cz /= total_mass

    # Store the per-node fields. Bounds are kept so the force walk can
    # evaluate the opening criterion (size / dist < theta).
    node_mass[node_idx] = total_mass
    node_com[node_idx, 0] = cx
    node_com[node_idx, 1] = cy
    node_com[node_idx, 2] = cz
    node_bounds[node_idx, 0] = bounds[0, 0]   # xmin
    node_bounds[node_idx, 1] = bounds[0, 1]   # xmax
    node_bounds[node_idx, 2] = bounds[1, 0]   # ymin
    node_bounds[node_idx, 3] = bounds[1, 1]   # ymax
    node_bounds[node_idx, 4] = bounds[2, 0]   # zmin
    node_bounds[node_idx, 5] = bounds[2, 1]   # zmax

    # Leaf termination: single particle, or hit the depth cap with several
    # particles still in the cell (pseudo-leaf, sentinel -2).
    if n == 1 or depth >= max_depth:
        node_leaf_idx[node_idx] = indices[0] if n == 1 else -2
        for c in range(8):
            node_children[node_idx, c] = -1
        return next_node_ptr

    # Subdivide into eight octants. Midpoints are computed once and reused
    # both for the child bounds and for the per-particle assignment below.
    xmin, xmax = bounds[0, 0], bounds[0, 1]
    ymin, ymax = bounds[1, 0], bounds[1, 1]
    zmin, zmax = bounds[2, 0], bounds[2, 1]
    xmid = (xmin + xmax) * 0.5
    ymid = (ymin + ymax) * 0.5
    zmid = (zmin + zmax) * 0.5

    # Bounds of the eight child cubes laid out in (xmin, xmax, ymin, ymax,
    # zmin, zmax) order. The octant index packs three binary axis choices
    # as ix + 2*iy + 4*iz (equivalent to a bit pattern).
    sub_bounds = np.empty((8, 6), dtype=np.float64)
    for ix in (0, 1):
        for iy in (0, 1):
            for iz in (0, 1):
                cidx = ix + iy * 2 + iz * 4
                xl = xmin if ix == 0 else xmid
                xr = xmid if ix == 0 else xmax
                yl = ymin if iy == 0 else ymid
                yr = ymid if iy == 0 else ymax
                zl = zmin if iz == 0 else zmid
                zr = zmid if iz == 0 else zmax
                sub_bounds[cidx, 0] = xl
                sub_bounds[cidx, 1] = xr
                sub_bounds[cidx, 2] = yl
                sub_bounds[cidx, 3] = yr
                sub_bounds[cidx, 4] = zl
                sub_bounds[cidx, 5] = zr

    # Bucket particles into the eight octants. The temporary (8, n) buffer
    # is preallocated to avoid dynamic resizing inside @njit. Each row uses
    # only the first child_counts[c] entries.
    child_lists = np.full((8, n), -1, dtype=np.int32)
    child_counts = np.zeros(8, dtype=np.int32)
    for i in range(n):
        p_idx = indices[i]
        x = particles[p_idx, 0]
        y = particles[p_idx, 1]
        z = particles[p_idx, 2]
        ix = 0 if x < xmid else 1
        iy = 0 if y < ymid else 1
        iz = 0 if z < zmid else 1
        cidx = ix + iy * 2 + iz * 4
        child_lists[cidx, child_counts[cidx]] = p_idx
        child_counts[cidx] += 1

    # Recurse into every non-empty octant. The child root index is read
    # before the recursive call because next_node_ptr advances inside it.
    for c in range(8):
        cnt = child_counts[c]
        if cnt > 0:
            sub_indices = child_lists[c, :cnt]
            child_root = next_node_ptr[0]
            next_node_ptr = build_flat_tree(
                particles, sub_indices, sub_bounds[c].reshape(3, 2),
                masses, node_mass, node_com, node_children, node_bounds,
                node_leaf_idx, depth + 1, max_depth, next_node_ptr)
            node_children[node_idx, c] = child_root
        else:
            node_children[node_idx, c] = -1

    # Sentinel -1 marks this node as internal (not a leaf).
    node_leaf_idx[node_idx] = -1
    return next_node_ptr


# ---------------------------------------------------------------------------
# Per-particle force walk and parallel driver.
# ---------------------------------------------------------------------------

@njit
def acceleration_from_node(particles, target_idx, softening_sq,
                           node_mass, node_com, node_children, node_bounds,
                           node_leaf_idx, node_root, theta, G):
    """
    Walk the flat octree for one target particle and accumulate acceleration.

    Implemented as an explicit stack (length 2000 is more than enough for
    the depth-40 tree used here) rather than Python recursion, both because
    Numba prefers loops and because each prange thread needs a private stack.

    Args:
        particles: (N, 3) positions
        target_idx: Index of the particle receiving force
        softening_sq: Plummer softening squared
        node_mass, node_com, node_children, node_bounds, node_leaf_idx:
            Flat tree buffers
        node_root: Index of the root node (0 in the current driver)
        theta: Barnes-Hut opening angle
        G: Gravitational constant

    Returns:
        Tuple (ax, ay, az) acceleration components on target_idx.
    """
    ax = 0.0
    ay = 0.0
    az = 0.0
    # Per-thread DFS stack. Depth-40 trees fit very comfortably.
    stack = np.zeros(2000, dtype=np.int32)
    stack[0] = node_root
    top = 1
    while top > 0:
        top -= 1
        node = stack[top]

        mass_node = node_mass[node]
        if mass_node == 0.0:
            continue

        # Vector from target to centre of mass; softened distance for the
        # opening criterion and force.
        dx = node_com[node, 0] - particles[target_idx, 0]
        dy = node_com[node, 1] - particles[target_idx, 1]
        dz = node_com[node, 2] - particles[target_idx, 2]
        dist_sq = dx * dx + dy * dy + dz * dz + softening_sq
        dist = np.sqrt(dist_sq)

        leaf = node_leaf_idx[node]
        if leaf == -2:
            # Pseudo-leaf: max_depth reached with multiple particles. Use
            # the node centre of mass directly.
            factor = G * mass_node / (dist_sq * dist)
            ax += factor * dx
            ay += factor * dy
            az += factor * dz
        elif leaf != -1:
            # Single-particle leaf: skip self, otherwise compute the exact
            # pair force using the particle's actual position.
            if leaf != target_idx:
                px = particles[leaf, 0]
                py = particles[leaf, 1]
                pz = particles[leaf, 2]
                dx2 = px - particles[target_idx, 0]
                dy2 = py - particles[target_idx, 1]
                dz2 = pz - particles[target_idx, 2]
                dist_sq2 = dx2 * dx2 + dy2 * dy2 + dz2 * dz2 + softening_sq
                dist_cube = dist_sq2 * np.sqrt(dist_sq2)
                factor = G * mass_node / dist_cube
                ax += factor * dx2
                ay += factor * dy2
                az += factor * dz2
        else:
            # Internal node: opening criterion uses the cube side length
            # (x-axis span) versus the softened distance to the COM.
            size = node_bounds[node, 1] - node_bounds[node, 0]
            if size / dist < theta:
                # Far enough: approximate by the centre of mass.
                factor = G * mass_node / (dist_sq * dist)
                ax += factor * dx
                ay += factor * dy
                az += factor * dz
            else:
                # Too close: push all non-empty children for later visiting.
                for c in range(8):
                    child = node_children[node, c]
                    if child != -1:
                        stack[top] = child
                        top += 1
    return ax, ay, az


@njit(parallel=True)
def compute_all_accelerations_parallel(particles, node_root,
                                       node_mass, node_com, node_children,
                                       node_bounds, node_leaf_idx,
                                       softening_sq, theta, G):
    """
    Parallel per-particle force evaluation via Numba prange.

    Each iteration of the prange loop computes the acceleration on one
    particle independently, so the outer loop is embarrassingly parallel
    and scales with available CPU cores.

    Args:
        particles: (N, 3) positions
        node_root: Index of the root node
        node_mass, node_com, node_children, node_bounds, node_leaf_idx:
            Flat tree buffers
        softening_sq, theta, G: Force-kernel parameters

    Returns:
        (N, 3) array of accelerations.
    """
    N = particles.shape[0]
    acc = np.zeros((N, 3), dtype=np.float64)
    for i in prange(N):
        ax, ay, az = acceleration_from_node(
            particles, i, softening_sq,
            node_mass, node_com, node_children, node_bounds, node_leaf_idx,
            node_root, theta, G)
        acc[i, 0] = ax
        acc[i, 1] = ay
        acc[i, 2] = az
    return acc


# ---------------------------------------------------------------------------
# Main driver.
# ---------------------------------------------------------------------------

def main():
    """
    Load the initial state, integrate for n_steps, and stream snapshots to disk.

    The root-cube bounds are computed once from the initial extent and
    held fixed; the flat tree buffers are sized to 2 * N + 1000 nodes,
    which comfortably exceeds the 2N - 1 worst-case node count of an
    octree with N particles. The integrator is a symplectic leapfrog with
    an explicit initial half-kick.
    """
    os.makedirs(output_dir, exist_ok=True)
    data = np.load("data/galaxy_initial_state.npz")
    positions = data["pos"].astype(np.float64)
    velocities = data["vel"].astype(np.float64)
    masses = data["mass"].astype(np.float64)
    labels = data["label"]
    N = len(masses)
    print(f"Loaded {N} particles, total mass = {np.sum(masses):.2f}")

    # Root cube of the octree; padded by 5 units to absorb early motion.
    bounds_global = get_bounds(positions, margin=5.0)
    softening_sq = softening ** 2

    # Preallocate the flat tree buffers. 2 * N is an upper bound on the
    # number of octree nodes (N leaves + at most N - 1 internal); the +1000
    # is a comfort margin.
    max_nodes = N * 2 + 1000
    node_mass = np.zeros(max_nodes, dtype=np.float64)
    node_com = np.zeros((max_nodes, 3), dtype=np.float64)
    node_children = np.full((max_nodes, 8), -1, dtype=np.int32)
    node_bounds = np.zeros((max_nodes, 6), dtype=np.float64)
    node_leaf_idx = np.full(max_nodes, -1, dtype=np.int32)

    # Build the initial tree. next_node_ptr is a length-1 array used by
    # build_flat_tree as a mutable allocation pointer.
    indices_all = np.arange(N, dtype=np.int32)
    next_node_ptr = np.array([0], dtype=np.int32)
    build_flat_tree(positions, indices_all, bounds_global, masses,
                    node_mass, node_com, node_children, node_bounds,
                    node_leaf_idx, 0, 40, next_node_ptr)
    root = 0   # Root node index is always 0 after a fresh build.
    print(f"Initial tree built, used {next_node_ptr[0]} nodes.")

    # Initial half-kick so the per-step kernel can be uniform.
    acc0 = compute_all_accelerations_parallel(
        positions, root,
        node_mass, node_com, node_children, node_bounds, node_leaf_idx,
        softening_sq, theta, G)
    velocities_half = velocities + 0.5 * dt * acc0
    print("Initialization done.")

    # Step 0 snapshot uses the original (integer-time) velocities.
    np.savez(os.path.join(output_dir, "galaxy_step_0000.npz"),
             pos=positions, vel=velocities, mass=masses, label=labels)
    print(f"Step    0/{n_steps}  |  time =  0.00  |  elapsed =  0.0s")

    # Main loop: drift, rebuild tree, evaluate force, full kick on the
    # half-step velocities. This is the standard symplectic leapfrog.
    start_time = time.time()
    for step in range(1, n_steps + 1):
        # 1. Drift: x(t) -> x(t + dt) using the half-step velocity.
        positions = positions + velocities_half * dt

        # 2. Rebuild the tree at the new positions and recompute a(t + dt).
        next_node_ptr[0] = 0
        build_flat_tree(positions, indices_all, bounds_global, masses,
                        node_mass, node_com, node_children, node_bounds,
                        node_leaf_idx, 0, 40, next_node_ptr)
        acc = compute_all_accelerations_parallel(
            positions, root,
            node_mass, node_com, node_children, node_bounds, node_leaf_idx,
            softening_sq, theta, G)

        # 3. Kick: v(t + dt/2) -> v(t + 3dt/2) with a full-step kick.
        velocities_half = velocities_half + acc * dt

        # 4. Snapshot. v(t + dt) is recovered by undoing half a kick.
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
