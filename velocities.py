"""
Velocity sampling for the three-component galaxy IC.

Implements the velocity-initialisation prescription of Hernquist (1993a) as
described in Chapter 2.3.2 of the practicum manual (Naab 2006).

Disk -- anisotropic Gaussian in cylindrical coordinates:
    sigma_z^2(R)   = pi G Sigma_d(R) z_0                                (eq 2.40)
    sigma_R^2(R)   = sigma_R^2(R_crit) * Sigma_d(R) / Sigma_d(R_crit)   (eq 2.39)
    sigma_R(R_crit) = Q * 3.36 G Sigma_d(R_crit) / kappa(R_crit)        (eq 2.47-2.48)
    sigma_phi^2(R) = sigma_R^2 * kappa^2 / (4 Omega^2)                  (eq 2.53)
    <v_phi>^2(R)   = v_c^2 + sigma_R^2 [1 - kappa^2/(4 Omega^2) - 2R/h] (eq 2.56)

Bulge and halo -- isotropic Maxwellian, dispersion from the spherical Jeans
equation (eq 2.63) evaluated in the *total* potential.

The asymmetric-drift expression (eq 2.56) can go negative at intermediate
radii; following Hernquist's "softening" of v_R^2 at small radii (mentioned in
the manual), we clamp <v_phi>^2 at zero. Negative sigma_R^2 from the Toomre
prescription is similarly clamped.

CMPH Project 3 -- Zoutong Shen / Zhaoyang Chu, 2026.
"""

from typing import Optional, Tuple

import numpy as np

from initial_conditions import GalaxyParams, galaxy_initial_conditions
from potential import GalaxyPotential


_TOOMRE_Q = 1.5
_CRITICAL_RADIUS_FACTOR = 2.4


def disk_velocities(
    disk_position_array: np.ndarray,
    potential: GalaxyPotential,
    rng: np.random.Generator,
    toomre_q: float = _TOOMRE_Q,
    critical_radius_factor: float = _CRITICAL_RADIUS_FACTOR,
) -> np.ndarray:
    """
    Sample disk velocities at the given disk-particle positions.

    Args:
        disk_position_array:     shape (N_disk, 3).
        potential:               GalaxyPotential built with matching params.
        rng:                     numpy random generator.
        toomre_q:                Toomre Q at the critical radius (default 1.5).
        critical_radius_factor:  R_crit / h (default 2.4, manual page 16).

    Returns:
        Array of shape (N_disk, 3) -- velocities in Cartesian coordinates.
    """
    params = potential.params
    h = params.disk_scale_length
    z_0 = params.disk_scale_height
    critical_radius = critical_radius_factor * h

    cylindrical_R = np.sqrt(
        disk_position_array[:, 0] ** 2 + disk_position_array[:, 1] ** 2
    )
    azimuthal_angle = np.arctan2(
        disk_position_array[:, 1], disk_position_array[:, 0]
    )

    surface_density_at_particle = potential.disk_surface_density(cylindrical_R)
    surface_density_at_critical = potential.disk_surface_density(critical_radius)
    kappa_squared_at_critical = potential.epicyclic_frequency_squared(critical_radius)

    sigma_z_squared = np.pi * surface_density_at_particle * z_0
    sigma_R_at_critical_squared = (
        (toomre_q * 3.36 * surface_density_at_critical) ** 2
        / kappa_squared_at_critical
    )
    sigma_R_squared = (
        sigma_R_at_critical_squared
        * surface_density_at_particle
        / surface_density_at_critical
    )
    sigma_R_squared = np.maximum(sigma_R_squared, 0.0)

    v_circular_squared = potential.circular_speed_squared(cylindrical_R)
    omega_squared = v_circular_squared / np.maximum(cylindrical_R ** 2, 1.0e-30)
    kappa_squared = potential.epicyclic_frequency_squared(cylindrical_R)
    epicyclic_ratio = kappa_squared / (4.0 * np.maximum(omega_squared, 1.0e-30))

    sigma_phi_squared = sigma_R_squared * epicyclic_ratio

    asymmetric_drift_bracket = 1.0 - epicyclic_ratio - 2.0 * cylindrical_R / h
    mean_v_phi_squared = v_circular_squared + sigma_R_squared * asymmetric_drift_bracket
    mean_v_phi_squared = np.maximum(mean_v_phi_squared, 0.0)
    mean_v_phi = np.sqrt(mean_v_phi_squared)

    sigma_R = np.sqrt(sigma_R_squared)
    sigma_phi = np.sqrt(sigma_phi_squared)
    sigma_z = np.sqrt(sigma_z_squared)

    velocity_radial = rng.normal(0.0, sigma_R)
    velocity_vertical = rng.normal(0.0, sigma_z)
    velocity_azimuthal = rng.normal(mean_v_phi, sigma_phi)

    velocities = np.empty_like(disk_position_array)
    cos_phi = np.cos(azimuthal_angle)
    sin_phi = np.sin(azimuthal_angle)
    velocities[:, 0] = velocity_radial * cos_phi - velocity_azimuthal * sin_phi
    velocities[:, 1] = velocity_radial * sin_phi + velocity_azimuthal * cos_phi
    velocities[:, 2] = velocity_vertical
    return velocities


def _isotropic_velocities(
    positions: np.ndarray,
    dispersion_squared_at_particle: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Sample 3D Cartesian velocity components from N(0, sigma) i.i.d., with sigma
    set per particle. Yields a 3D isotropic Maxwell-Boltzmann distribution
    matching eq 2.65 of the manual.
    """
    sigma = np.sqrt(np.maximum(dispersion_squared_at_particle, 0.0))
    sigma_column = sigma[:, None]
    return rng.normal(0.0, 1.0, positions.shape) * sigma_column


def bulge_velocities(
    bulge_position_array: np.ndarray,
    potential: GalaxyPotential,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Sample bulge velocities via isotropic Jeans (eq 2.63 in the total potential).

    Args:
        bulge_position_array: shape (N_bulge, 3).
        potential:            GalaxyPotential with matching params.
        rng:                  numpy random generator.

    Returns:
        Array of shape (N_bulge, 3) -- isotropic Cartesian velocities.
    """
    radii = np.linalg.norm(bulge_position_array, axis=1)
    dispersion_squared = potential.bulge_velocity_dispersion_squared(radii)
    return _isotropic_velocities(bulge_position_array, dispersion_squared, rng)


def halo_velocities(
    halo_position_array: np.ndarray,
    potential: GalaxyPotential,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Sample halo velocities via isotropic Jeans (eq 2.63 in the total potential).

    Args:
        halo_position_array: shape (N_halo, 3).
        potential:           GalaxyPotential with matching params.
        rng:                 numpy random generator.

    Returns:
        Array of shape (N_halo, 3) -- isotropic Cartesian velocities.
    """
    radii = np.linalg.norm(halo_position_array, axis=1)
    dispersion_squared = potential.halo_velocity_dispersion_squared(radii)
    return _isotropic_velocities(halo_position_array, dispersion_squared, rng)


def galaxy_velocities(
    positions: np.ndarray,
    potential: GalaxyPotential,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Sample velocities for all particles, assuming the disk -> bulge -> halo
    ordering convention from `initial_conditions.galaxy_initial_conditions`.
    """
    params = potential.params
    num_disk = params.num_disk
    num_bulge = params.num_bulge
    end_bulge = num_disk + num_bulge

    velocities = np.empty_like(positions)
    velocities[:num_disk] = disk_velocities(positions[:num_disk], potential, rng)
    velocities[num_disk:end_bulge] = bulge_velocities(
        positions[num_disk:end_bulge], potential, rng
    )
    velocities[end_bulge:] = halo_velocities(positions[end_bulge:], potential, rng)
    return velocities


def make_galaxy_initial_conditions(
    params: Optional[GalaxyParams] = None,
    seed: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build complete IC (positions + velocities + masses) for one galaxy.

    Args:
        params: GalaxyParams instance; manual defaults if None.
        seed:   integer seed for reproducible sampling.

    Returns:
        positions:  shape (N_total, 3).
        velocities: shape (N_total, 3).
        masses:     shape (N_total,).
    """
    if params is None:
        params = GalaxyParams()

    positions, masses = galaxy_initial_conditions(params=params, seed=seed)
    potential = GalaxyPotential(params)

    rng = np.random.default_rng(
        None if seed is None else seed + 1_000_003
    )
    velocities = galaxy_velocities(positions, potential, rng)
    return positions, velocities, masses


def main() -> None:
    """Build the default-parameter IC and store the complete (x, v, m) to disk."""
    params = GalaxyParams()
    positions, velocities, masses = make_galaxy_initial_conditions(params=params, seed=0)

    print(f"Galaxy IC built -- {positions.shape[0]} particles")
    print(f"  net momentum magnitude / total mass = "
          f"{np.linalg.norm((velocities * masses[:, None]).sum(axis=0)) / masses.sum():.3e}")
    print(f"  net angular momentum L_z (disk-component dominant) = "
          f"{(masses * (positions[:, 0] * velocities[:, 1] - positions[:, 1] * velocities[:, 0])).sum():.4f}")
    print(f"  RMS disk in-plane speed = "
          f"{np.sqrt((velocities[:params.num_disk, :2] ** 2).mean()):.4f}")
    print(f"  RMS halo speed (all components) = "
          f"{np.sqrt((velocities[params.num_disk + params.num_bulge:] ** 2).mean()):.4f}")

    output_path = "data/spiral.npz"
    np.savez(output_path, positions=positions, velocities=velocities, masses=masses)
    print(f"\nSaved to {output_path}")


if __name__ == "__main__":
    main()
