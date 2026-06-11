"""
Initial conditions for a three-component disk galaxy.

Generates positions and (zero) velocities for three particle species and
writes them to a single .npz file consumed by the N-body driver:

    Bulge  (label 0):  Hernquist-like profile,  rho(r) ~ 1 / [r(r+a)^3]
                       inverse-CDF sampling  r = a sqrt(u) / (1 - sqrt(u)).
    Disk   (label 1):  exponential radial profile, sech^2 vertical profile;
                       radial draw uses a Gamma(2, h) ansatz that mimics
                       R exp(-R/h).
    Halo   (label 2):  cored profile r^2 / (r^2 + gamma^2) with an
                       exponential cutoff exp(-(r/r_c)^2), sampled by
                       rejection.

Velocities are set to zero for all three species. The system therefore
starts cold and is expected to collapse and virialise during the run; the
intent here is to provide reproducible initial geometry rather than a
self-consistent equilibrium.

Author: CMPH 2026 Project 3
Date: May 2026
"""

import numpy as np
import matplotlib.pyplot as plt


def generate_disk_positions(
    N_particles: int, *, h: float = 1.0, z_0: float = 0.1
) -> tuple[np.ndarray, np.ndarray]:
    """
    Sample positions for an exponential disk with a sech^2 vertical profile.

    Radii are drawn from Gamma(shape=2, scale=h), whose density
    R exp(-R/h) matches the area element times an exponential surface
    density. Vertical positions follow z_0 * arctanh(2u-1) with u uniform,
    which is the inverse CDF of sech^2(z / z_0). Velocities are returned
    as zero so that the disk starts cold and develops rotation only
    through self-gravity during the simulation.

    Args:
        N_particles: Number of disk particles.
        h: Radial scale length.
        z_0: Vertical scale height.

    Returns:
        Tuple (pos, vel), each of shape (N_particles, 3).
    """
    phi = np.random.uniform(0, 2 * np.pi, N_particles)
    R = np.random.gamma(shape=2.0, scale=h, size=N_particles)

    u = np.random.uniform(0, 1, N_particles)
    eps = 1e-6
    # Clip away from 0 and 1 so that arctanh(2u-1) never diverges.
    u = np.clip(u, eps, 1 - eps)
    z = z_0 * np.arctanh(2 * u - 1)

    x = R * np.cos(phi)
    y = R * np.sin(phi)
    pos = np.column_stack((x, y, z))

    # Cold start: zero initial velocity in every direction.
    vx = np.zeros(N_particles)
    vy = np.zeros(N_particles)
    vz = np.zeros(N_particles)
    vel = np.column_stack((vx, vy, vz))

    return pos, vel


def generate_bulge_positions(
    N_particles: int, *, a: float = 0.1
) -> tuple[np.ndarray, np.ndarray]:
    """
    Sample positions for a Hernquist-like bulge with scale radius a.

    The radial coordinate is drawn from the inverse CDF
    r = a sqrt(u) / (1 - sqrt(u)) of the Hernquist enclosed-mass profile
    M(<r)/M = r^2 / (r + a)^2. Because that CDF has a heavy tail, u is
    clipped just below 1 to bound the maximum sampled radius; otherwise a
    handful of particles can land at unphysically large r. Angular
    coordinates are sampled isotropically on the sphere.

    Args:
        N_particles: Number of bulge particles.
        a: Hernquist scale radius.

    Returns:
        Tuple (pos, vel), each of shape (N_particles, 3); vel is zero.
    """
    u = np.random.uniform(0, 1, N_particles)
    # Clip u close to 1 to truncate the Hernquist tail; values too close
    # to 1 push r to extreme magnitudes (r ~ a / (1 - sqrt(u))).
    u = np.clip(u, 0, 0.9999)
    r = a * (np.sqrt(u) / (1.0 - np.sqrt(u)))

    phi = np.random.uniform(0, 2 * np.pi, N_particles)
    costheta = np.random.uniform(-1, 1, N_particles)
    sintheta = np.sqrt(1 - costheta**2)

    x = r * sintheta * np.cos(phi)
    y = r * sintheta * np.sin(phi)
    z = r * costheta
    pos = np.column_stack((x, y, z))

    # Cold start: zero initial velocity.
    vx = np.zeros(N_particles)
    vy = np.zeros(N_particles)
    vz = np.zeros(N_particles)
    vel = np.column_stack((vx, vy, vz))

    return pos, vel


def generate_halo_positions(
    N_particles: int, *, gamma: float = 1.0, r_c: float = 10.0
) -> tuple[np.ndarray, np.ndarray]:
    """
    Sample positions for a cored, exponentially-truncated dark matter halo.

    The radial density follows r^2 / (r^2 + gamma^2) * exp(-(r/r_c)^2):
    a soft core of size gamma and a Gaussian cutoff at r_c. The function
    samples r by rejection inside r in [0, 3 r_c]; the loop refills if
    fewer than N_particles candidates are accepted in one batch.

    Args:
        N_particles: Number of halo particles.
        gamma: Core radius regularising the r -> 0 behaviour.
        r_c: Cutoff radius for the Gaussian truncation.

    Returns:
        Tuple (pos, vel), each of shape (N_particles, 3); vel is zero.
    """
    r_list: list[float] = []
    r_max = 3.0 * r_c

    while len(r_list) < N_particles:
        r_candidate = np.random.uniform(0, r_max, N_particles)
        # Acceptance probability proportional to the radial profile.
        prob_r = (r_candidate**2 / (r_candidate**2 + gamma**2)) * np.exp(
            -(r_candidate**2) / (r_c**2)
        )
        y = np.random.uniform(0, 1, N_particles)
        accepted_r = r_candidate[y < prob_r]
        r_list.extend(accepted_r)

    r = np.array(r_list[:N_particles])

    phi = np.random.uniform(0, 2 * np.pi, N_particles)
    costheta = np.random.uniform(-1, 1, N_particles)
    sintheta = np.sqrt(1 - costheta**2)

    x = r * sintheta * np.cos(phi)
    y = r * sintheta * np.sin(phi)
    z = r * costheta
    pos = np.column_stack((x, y, z))

    # Cold start: zero initial velocity.
    vx = np.zeros(N_particles)
    vy = np.zeros(N_particles)
    vz = np.zeros(N_particles)
    vel = np.column_stack((vx, vy, vz))

    return pos, vel


if __name__ == "__main__":
    # Particle counts: halo dominates by number, then disk, then bulge.
    n_disk = 6000
    n_bulge = 2000
    n_halo = 12000

    # Mass ratios: dark matter halo dominates, then disk, then bulge.
    mass_disk_total = 1.0
    mass_bulge_total = 0.25
    mass_halo_total = 3.0

    # Equal-mass particles within each species; the species itself has a
    # different per-particle mass.
    m_disk = mass_disk_total / n_disk
    m_bulge = mass_bulge_total / n_bulge
    m_halo = mass_halo_total / n_halo

    pos_disk, vel_disk = generate_disk_positions(n_disk)
    pos_bulge, vel_bulge = generate_bulge_positions(n_bulge)
    pos_halo, vel_halo = generate_halo_positions(n_halo)

    # Concatenate in a fixed order: bulge, disk, halo. Downstream code
    # reads the label array, so this order is documentation, not load-bearing.
    positions = np.vstack((pos_bulge, pos_disk, pos_halo))
    velocities = np.vstack((vel_bulge, vel_disk, vel_halo))
    masses = np.concatenate((
        np.full(n_bulge, m_bulge),
        np.full(n_disk, m_disk),
        np.full(n_halo, m_halo),
    ))

    # Species labels used by the visualiser to colour points.
    labels = np.concatenate((
        np.full(n_bulge, 0, dtype=int),    # 0 = bulge
        np.full(n_disk, 1, dtype=int),     # 1 = disk
        np.full(n_halo, 2, dtype=int),     # 2 = halo
    ))

    np.savez("data/galaxy_initial_state.npz",
             pos=positions,
             vel=velocities,
             mass=masses,
             label=labels)
    print("Initial state saved to data/galaxy_initial_state.npz")
    print(f"  Total particles: {len(positions)}")
    print(f"  Bulge: {n_bulge}, Disk: {n_disk}, Halo: {n_halo}")
    print(f"  Total mass: {np.sum(masses):.2f}")

    # Plot initial configuration for visual verification.
    plt.style.use("dark_background")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    ax1.scatter(pos_halo[:, 0], pos_halo[:, 1], s=0.2, c="cyan", alpha=0.8, label="Halo")
    ax1.scatter(pos_disk[:, 0], pos_disk[:, 1], s=0.8, c="white", alpha=0.8, label="Disk")
    ax1.scatter(pos_bulge[:, 0], pos_bulge[:, 1], s=1.5, c="red", alpha=0.8, label="Bulge")
    ax1.set_title("Face-on View (x-y plane)")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_xlim(-8, 8)
    ax1.set_ylim(-8, 8)
    ax1.set_aspect("equal")

    ax2.scatter(pos_halo[:, 0], pos_halo[:, 2], s=0.2, c="cyan", alpha=0.8, label="Halo")
    ax2.scatter(pos_disk[:, 0], pos_disk[:, 2], s=0.8, c="white", alpha=0.8, label="Disk")
    ax2.scatter(pos_bulge[:, 0], pos_bulge[:, 2], s=1.5, c="red", alpha=0.8, label="Bulge")
    ax2.set_title("Edge-on View (x-z plane)")
    ax2.set_xlabel("x")
    ax2.set_ylabel("z")
    ax2.set_xlim(-8, 8)
    ax2.set_ylim(-4, 4)
    ax2.set_aspect("equal")

    plt.tight_layout()
    plt.savefig("figure/single_initial.png", dpi=300, bbox_inches="tight")
    print("Plot saved to figure/single_initial.png")
