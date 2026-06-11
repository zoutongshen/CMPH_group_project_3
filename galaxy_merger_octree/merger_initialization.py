"""
Initial conditions for a two-galaxy encounter.

Two disk galaxies (bulge + disk + dark halo) are placed on the z-axis
and given opposite bulk velocities so they approach each other in the
centre-of-mass frame. The downstream simulation (merger_simulation.py)
takes over and integrates the encounter; the specific outcome (fly-by,
ring, merger, ...) is set entirely by the tunable parameters in the
__main__ block below:

    n_bulge / n_disk / n_halo (per galaxy)  -- particle counts per species
    mass_bulge_total / mass_disk_total / mass_halo_total
                                            -- total mass per species (galaxy 1)
    mass_ratio                              -- m_galaxy_1 / m_galaxy_2; controls
                                               whether galaxy 2 is a smaller
                                               intruder (>1), equal-mass (=1),
                                               or larger (<1)
    separation                              -- half the initial galaxy-galaxy
                                               distance along z
    v_rel                                   -- relative approach speed in the
                                               COM frame; split between the two
                                               galaxies in inverse proportion to
                                               their masses
    tilt_1 / tilt_2                         -- rotation angle (rad) of each
                                               disk plane around x, i.e. the
                                               angle between the disk spin axis
                                               and the collision axis (z)
    b_impact                                -- non-zero offsets along x can be
                                               added to center_pos to give the
                                               encounter an orbital angular
                                               momentum

The output is a COLD merger IC: bulge and halo carry zero internal
velocity dispersion. The pre-equilibration that lets the bulge/halo
virialise around the rotating disk is done inside merger_simulation.py,
so this file only has to handle geometry.

Output: data/merger_initial_state.npz with the (pos, vel, mass, label)
layout consumed by the N-body driver. Labels use 0-5 to distinguish the
two galaxies as well as their species:

    0 = galaxy 1 bulge   3 = galaxy 2 bulge
    1 = galaxy 1 disk    4 = galaxy 2 disk
    2 = galaxy 1 halo    5 = galaxy 2 halo

Author: CMPH 2026 Project 3
Date: May 2026
"""

import os

import numpy as np
import matplotlib.pyplot as plt

from galaxy_initialization import (
    generate_disk_positions,
    generate_bulge_positions,
    generate_halo_positions,
)

G = 1.0


def assign_disk_rotation(pos_disk, all_pos_local, all_mass_local):
    """
    Assign tangential circular velocities to the disk particles.

    Each disk particle is given the velocity v_circ = sqrt(G M(<R) / R) in
    the direction tangent to the (x, y) plane, with the sign chosen so that
    the angular momentum points along +z (counter-clockwise as seen from
    +z). M(<R) is the total mass of all particles (bulge + disk + halo)
    enclosed within a spherical radius R, looked up via a sorted cumulative
    sum, so the whole assignment is O(N log N) rather than O(N^2).

    Args:
        pos_disk: (N_disk, 3) disk particle positions in the galaxy's local
            frame (disk in the x-y plane).
        all_pos_local: (N_total, 3) positions of every particle in this
            galaxy in the local frame.
        all_mass_local: (N_total,) per-particle masses for all_pos_local.

    Returns:
        (N_disk, 3) velocities in the local frame.
    """
    R = np.sqrt(pos_disk[:, 0] ** 2 + pos_disk[:, 1] ** 2)
    r_all = np.linalg.norm(all_pos_local, axis=1)

    # Cumulative enclosed mass as a function of sorted radius.
    order = np.argsort(r_all)
    r_sorted = r_all[order]
    m_cumsum = np.cumsum(all_mass_local[order])

    # M(<R) for each disk particle by inserting R into the sorted radii.
    idx = np.searchsorted(r_sorted, R, side="right")
    M_enc = np.zeros_like(R)
    valid = idx > 0
    M_enc[valid] = m_cumsum[idx[valid] - 1]

    # Floor R so v_circ does not diverge for particles very close to the centre.
    R_safe = np.maximum(R, 0.1)
    v_circ = np.sqrt(G * M_enc / R_safe)

    phi = np.arctan2(pos_disk[:, 1], pos_disk[:, 0])
    vx = -v_circ * np.sin(phi)
    vy = v_circ * np.cos(phi)
    vz = np.zeros_like(vx)
    return np.column_stack((vx, vy, vz))


def rotation_matrix_x(angle_rad):
    """3x3 rotation matrix around the x-axis by angle_rad."""
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    return np.array([[1.0, 0.0, 0.0],
                     [0.0,   c,  -s],
                     [0.0,   s,   c]])


def build_galaxy(n_bulge, n_disk, n_halo,
                 m_bulge, m_disk, m_halo,
                 center_pos, center_vel, tilt_angle_rad):
    """
    Build one rotated, translated galaxy and return its phase-space arrays.

    Steps:
        1. Sample bulge, disk and halo positions in the galaxy's local frame
           with the disk in the x-y plane.
        2. Assign circular velocities to the disk (rotation around +z) from
           the enclosed mass profile of the galaxy.
        3. Rotate every position and velocity by tilt_angle_rad around the
           x-axis to tilt the disk plane.
        4. Translate by center_pos / center_vel into the simulation frame.

    Args:
        n_bulge, n_disk, n_halo: Particle counts per species.
        m_bulge, m_disk, m_halo: Per-particle masses per species.
        center_pos: (3,) galaxy centre position in the simulation frame.
        center_vel: (3,) galaxy bulk velocity in the simulation frame.
        tilt_angle_rad: Rotation angle of the disk plane around x.

    Returns:
        Tuple (positions, velocities, masses) for this galaxy. Each array
        is ordered bulge, disk, halo so the caller can compose labels.
    """
    pos_bulge, _ = generate_bulge_positions(n_bulge)
    pos_disk, _ = generate_disk_positions(n_disk)
    pos_halo, _ = generate_halo_positions(n_halo)

    # Stack local arrays once for the enclosed-mass lookup.
    all_pos = np.vstack((pos_bulge, pos_disk, pos_halo))
    all_mass = np.concatenate((
        np.full(n_bulge, m_bulge),
        np.full(n_disk, m_disk),
        np.full(n_halo, m_halo),
    ))

    vel_bulge = np.zeros_like(pos_bulge)
    vel_disk = assign_disk_rotation(pos_disk, all_pos, all_mass)
    vel_halo = np.zeros_like(pos_halo)

    pos_local = np.vstack((pos_bulge, pos_disk, pos_halo))
    vel_local = np.vstack((vel_bulge, vel_disk, vel_halo))

    # Rotate around x to tilt the disk plane, then translate into place.
    R = rotation_matrix_x(tilt_angle_rad)
    pos_rot = pos_local @ R.T
    vel_rot = vel_local @ R.T

    positions = pos_rot + center_pos
    velocities = vel_rot + center_vel
    masses = all_mass

    return positions, velocities, masses


if __name__ == "__main__":
    np.random.seed(42)

    # Galaxy 1: particle counts and total mass per species.
    n_bulge_1 = 2000
    n_disk_1 = 6000
    n_halo_1 = 12000
    mass_bulge_total = 0.25
    mass_disk_total = 1.0
    mass_halo_total = 3.0
    m_bulge = mass_bulge_total / n_bulge_1
    m_disk = mass_disk_total / n_disk_1
    m_halo = mass_halo_total / n_halo_1

    # Galaxy 2: particle counts scale as 1/mass_ratio (per-particle masses
    # stay equal to galaxy 1, so total mass also scales as 1/mass_ratio).
    # mass_ratio = 1 -> equal-mass major merger; > 1 -> galaxy 2 is the
    # smaller intruder.
    mass_ratio = 1
    n_bulge_2 = n_bulge_1 // mass_ratio
    n_disk_2 = n_disk_1 // mass_ratio
    n_halo_2 = n_halo_1 // mass_ratio

    # Encounter geometry. The collision axis is x: galaxy 1 sits at -x
    # moving in +x, galaxy 2 sits at +x moving in -x. b_impact along y
    # or z can be added to center_pos to give the orbit angular momentum;
    # here it is zero (head-on along x). Velocities are split between
    # the two galaxies in inverse proportion to their masses so the
    # system has zero total momentum in the simulation frame.
    separation = 15.0
    v_rel = 0.5
    v_1 = v_rel / (mass_ratio + 1)
    v_2 = v_1 * mass_ratio

    # Disk-plane tilts around the x-axis. Disks are generated in the x-y
    # plane (spin axis along +z), so a rotation around x rotates the disk
    # within the plane perpendicular to the collision axis.
    tilt_1 = 0.0
    tilt_2 = np.deg2rad(30.0)

    print("Building galaxy 1...")
    pos1, vel1, mass1 = build_galaxy(
        n_bulge_1, n_disk_1, n_halo_1,
        m_bulge, m_disk, m_halo,
        center_pos=np.array([-separation, 0.0, 0.0]),
        center_vel=np.array([+v_1, 0.0, 0.0]),
        tilt_angle_rad=tilt_1,
    )

    print("Building galaxy 2...")
    pos2, vel2, mass2 = build_galaxy(
        n_bulge_2, n_disk_2, n_halo_2,
        m_bulge, m_disk, m_halo,
        center_pos=np.array([+separation, 0.0, 0.0]),
        center_vel=np.array([-v_2, 0.0, 0.0]),
        tilt_angle_rad=tilt_2,
    )

    # Labels 0-5 distinguish galaxy as well as species.
    labels1 = np.concatenate((
        np.full(n_bulge_1, 0, dtype=int),
        np.full(n_disk_1, 1, dtype=int),
        np.full(n_halo_1, 2, dtype=int),
    ))
    labels2 = np.concatenate((
        np.full(n_bulge_2, 3, dtype=int),
        np.full(n_disk_2, 4, dtype=int),
        np.full(n_halo_2, 5, dtype=int),
    ))

    positions = np.vstack((pos1, pos2))
    velocities = np.vstack((vel1, vel2))
    masses = np.concatenate((mass1, mass2))
    labels = np.concatenate((labels1, labels2))

    np.savez("data/merger_initial_state.npz",
             pos=positions,
             vel=velocities,
             mass=masses,
             label=labels)
    print("Initial state saved to data/merger_initial_state.npz")
    print(f"  Total particles: {len(positions)}")
    print(f"  Galaxy 1: bulge={n_bulge_1}, disk={n_disk_1}, halo={n_halo_1}, "
          f"tilt={np.rad2deg(tilt_1):.1f} deg")
    print(f"  Galaxy 2: bulge={n_bulge_2}, disk={n_disk_2}, halo={n_halo_2}, "
          f"tilt={np.rad2deg(tilt_2):.1f} deg")
    print(f"  Galaxy 1 mass = {np.sum(mass1):.2f}, galaxy 2 mass = {np.sum(mass2):.2f}"
          f"   (ratio {np.sum(mass1)/np.sum(mass2):.2f}:1)")
    print(f"  Initial separation along x: {2 * separation:.1f} (head-on)")
    print(f"  Relative approach velocity: {v_rel:.2f}"
          f" (v_1 = {v_1:.2f}, v_2 = {v_2:.2f})")

    # Face-on and edge-on snapshots of the initial configuration.
    plt.style.use("dark_background")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    species_style = {
        0: dict(s=1.5, c="red",    alpha=0.8, label="Primary bulge"),
        1: dict(s=0.8, c="white",  alpha=0.8, label="Primary disk"),
        2: dict(s=0.2, c="cyan",   alpha=0.6, label="Primary halo"),
        3: dict(s=1.5, c="orange", alpha=0.8, label="Intruder bulge"),
        4: dict(s=0.8, c="yellow", alpha=0.8, label="Intruder disk"),
        5: dict(s=0.2, c="magenta", alpha=0.6, label="Intruder halo"),
    }

    # Collision is along x, so both x-y and x-z views show the
    # encounter with the collision axis running horizontally.
    for lab in [2, 5, 1, 4, 0, 3]:
        mask = labels == lab
        ax1.scatter(positions[mask, 0], positions[mask, 1], **species_style[lab])
        ax2.scatter(positions[mask, 0], positions[mask, 2], **species_style[lab])

    ax1.set_title("Collision axis (x-y)")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_xlim(-25, 25)
    ax1.set_ylim(-20, 20)
    ax1.set_aspect("equal")
    ax1.legend(loc="upper right", fontsize=7)

    ax2.set_title("Collision axis (x-z)")
    ax2.set_xlabel("x")
    ax2.set_ylabel("z")
    ax2.set_xlim(-25, 25)
    ax2.set_ylim(-20, 20)
    ax2.set_aspect("equal")

    plt.tight_layout()
    os.makedirs("figure", exist_ok=True)
    plt.savefig("figure/merger_initial.png", dpi=200, bbox_inches="tight")
    print("Plot saved to figure/merger_initial.png")
