"""
Equilibrium initial-condition spatial sampling for a Milky Way-like disk galaxy.

Follows Hernquist (1993a) as described in Chapter 2 of the practicum manual
(Naab 2006). Three components:

    1. Exponential disk × sech² vertical profile        (eq 2.38)
    2. Hernquist spheroidal bulge                        (eq 2.67)
    3. Truncated isothermal-sphere dark matter halo      (eq 2.57)

This module produces positions and per-particle masses only. Velocity
sampling (which requires the tabulated total potential, eq 2.56) is
implemented separately downstream.

Default parameters match the manual (page 18):
    h = 1, M_d = 1, z_0 = 0.2, M_b = 1/3, a = 0.1,
    M_h = 5.8, gamma = 1, r_c = 10.

Particles are ordered disk -> bulge -> halo, matching the spiral.ascii
convention from the manual (§3.0.3).

CMPH Project 3 -- Zoutong Shen / Zhaoyang Chu, 2026.
"""

import os
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


@dataclass(frozen=True)
class GalaxyParams:
    """Parameters of the three-component galaxy model (manual page 18)."""

    disk_scale_length: float = 1.0
    disk_scale_height: float = 0.2
    disk_mass: float = 1.0
    bulge_mass: float = 1.0 / 3.0
    bulge_scale: float = 0.1
    halo_mass: float = 5.8
    halo_core: float = 1.0
    halo_cutoff: float = 10.0
    num_disk: int = 6000
    num_bulge: int = 2000
    num_halo: int = 12000


def disk_positions(
    num_particles: int,
    scale_length: float,
    scale_height: float,
    rng: np.random.Generator,
    radial_truncation_factor: float = 15.0,
) -> np.ndarray:
    """
    Sample disk particle positions from rho ~ exp(-R/h) sech^2(z/z_0).

    The radial sampling inverts the analytic exponential-disk CDF
    P(<R) = 1 - (1 + R/h) exp(-R/h) numerically via a fine lookup table.
    The vertical sampling uses the analytic inverse of the sech^2 CDF:
    z = z_0 atanh(2u - 1).

    Args:
        num_particles:            number of disk particles to draw.
        scale_length:             radial scale length h.
        scale_height:             vertical scale height z_0.
        rng:                      numpy random generator.
        radial_truncation_factor: cap R at this multiple of h (default 15h
                                  encloses > 99.99% of disk mass).

    Returns:
        Array of shape (num_particles, 3) -- disk positions in Cartesian
        coordinates.
    """
    radius_max = radial_truncation_factor * scale_length
    radii_grid = np.linspace(0.0, radius_max, 10_000)
    cumulative = 1.0 - (1.0 + radii_grid / scale_length) * np.exp(
        -radii_grid / scale_length
    )

    uniform_radial = rng.uniform(0.0, cumulative[-1], num_particles)
    radii = np.interp(uniform_radial, cumulative, radii_grid)

    azimuths = rng.uniform(0.0, 2.0 * np.pi, num_particles)

    safe_uniform = rng.uniform(1.0e-12, 1.0 - 1.0e-12, num_particles)
    heights = scale_height * np.arctanh(2.0 * safe_uniform - 1.0)

    positions = np.empty((num_particles, 3))
    positions[:, 0] = radii * np.cos(azimuths)
    positions[:, 1] = radii * np.sin(azimuths)
    positions[:, 2] = heights
    return positions


def bulge_positions(
    num_particles: int,
    bulge_scale: float,
    rng: np.random.Generator,
    truncation_factor: float = 500.0,
) -> np.ndarray:
    """
    Sample bulge particle positions from the Hernquist (1990) profile.

    The Hernquist cumulative mass M(r) = M_b r^2 / (r + a)^2 admits an
    analytic inverse: r = a sqrt(u) / (1 - sqrt(u)) for u ~ U(0, u_max),
    where u_max corresponds to the truncation radius.

    Args:
        num_particles:     number of bulge particles to draw.
        bulge_scale:       Hernquist scale length a.
        rng:               numpy random generator.
        truncation_factor: cap r at this multiple of a (default 500 a
                           encloses > 99.6% of bulge mass).

    Returns:
        Array of shape (num_particles, 3) -- bulge positions in Cartesian
        coordinates, isotropically oriented.
    """
    radius_max = truncation_factor * bulge_scale
    u_max = (radius_max / (radius_max + bulge_scale)) ** 2

    uniform_mass = rng.uniform(0.0, u_max, num_particles)
    sqrt_mass = np.sqrt(uniform_mass)
    radii = bulge_scale * sqrt_mass / (1.0 - sqrt_mass)

    return _isotropic_positions(radii, rng)


def halo_positions(
    num_particles: int,
    halo_core: float,
    halo_cutoff: float,
    rng: np.random.Generator,
    radial_grid_extent_factor: float = 5.0,
    radial_grid_size: int = 20_000,
) -> np.ndarray:
    """
    Sample halo particle positions from the truncated isothermal-sphere
    variant rho_h(r) ~ exp(-r^2 / r_c^2) / (r^2 + gamma^2) (eq 2.57).

    The cumulative mass M(<r) has no closed form, so it is built by
    trapezoidal integration on a fine radial grid extending to
    `radial_grid_extent_factor` times r_c; sampling proceeds by inverse-CDF
    interpolation. The exponential factor exp(-r^2/r_c^2) makes the integral
    converge quickly, so 5 r_c is ample.

    Args:
        num_particles:             number of halo particles to draw.
        halo_core:                 core radius gamma.
        halo_cutoff:               cutoff radius r_c.
        rng:                       numpy random generator.
        radial_grid_extent_factor: outer extent of the CDF grid (units of r_c).
        radial_grid_size:          number of grid points (resolution of inversion).

    Returns:
        Array of shape (num_particles, 3) -- halo positions in Cartesian
        coordinates, isotropically oriented.
    """
    radius_max = radial_grid_extent_factor * halo_cutoff
    radii_grid = np.linspace(0.0, radius_max, radial_grid_size)

    density_unnormalised = np.exp(
        -((radii_grid / halo_cutoff) ** 2)
    ) / (radii_grid ** 2 + halo_core ** 2)
    integrand = radii_grid ** 2 * density_unnormalised

    spacing = np.diff(radii_grid)
    trapezoidal_increments = 0.5 * (integrand[1:] + integrand[:-1]) * spacing
    cumulative = np.concatenate(([0.0], np.cumsum(trapezoidal_increments)))
    cumulative /= cumulative[-1]

    uniform_mass = rng.uniform(0.0, 1.0, num_particles)
    radii = np.interp(uniform_mass, cumulative, radii_grid)

    return _isotropic_positions(radii, rng)


def _isotropic_positions(
    radii: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    """
    Place each particle at its sampled radius on a uniformly-random unit
    sphere direction.

    Args:
        radii: shape (num_particles,) -- radial distances from origin.
        rng:   numpy random generator.

    Returns:
        Array of shape (num_particles, 3).
    """
    num_particles = radii.shape[0]
    cos_polar = rng.uniform(-1.0, 1.0, num_particles)
    sin_polar = np.sqrt(1.0 - cos_polar ** 2)
    azimuths = rng.uniform(0.0, 2.0 * np.pi, num_particles)

    positions = np.empty((num_particles, 3))
    positions[:, 0] = radii * sin_polar * np.cos(azimuths)
    positions[:, 1] = radii * sin_polar * np.sin(azimuths)
    positions[:, 2] = radii * cos_polar
    return positions


def galaxy_initial_conditions(
    params: Optional[GalaxyParams] = None,
    seed: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build positions and per-particle masses for the full three-component
    galaxy.

    Particles are concatenated in the order disk -> bulge -> halo, matching
    the spiral.ascii convention (manual §3.0.3). All three components share
    the same coordinate origin (the galaxy's centre of mass).

    Args:
        params: GalaxyParams instance; if None, manual defaults are used.
        seed:   integer seed for reproducible sampling.

    Returns:
        positions: shape (N_total, 3).
        masses:    shape (N_total,).
    """
    if params is None:
        params = GalaxyParams()

    rng = np.random.default_rng(seed)

    positions_disk = disk_positions(
        params.num_disk,
        params.disk_scale_length,
        params.disk_scale_height,
        rng,
    )
    positions_bulge = bulge_positions(params.num_bulge, params.bulge_scale, rng)
    positions_halo = halo_positions(
        params.num_halo, params.halo_core, params.halo_cutoff, rng
    )

    positions = np.vstack([positions_disk, positions_bulge, positions_halo])

    mass_per_disk_particle = params.disk_mass / params.num_disk
    mass_per_bulge_particle = params.bulge_mass / params.num_bulge
    mass_per_halo_particle = params.halo_mass / params.num_halo

    masses = np.concatenate(
        [
            np.full(params.num_disk, mass_per_disk_particle),
            np.full(params.num_bulge, mass_per_bulge_particle),
            np.full(params.num_halo, mass_per_halo_particle),
        ]
    )
    return positions, masses


def main() -> None:
    """Build the default IC, print a summary, and store the result to disk."""
    params = GalaxyParams()
    positions, masses = galaxy_initial_conditions(params=params, seed=0)

    num_total = positions.shape[0]
    print(f"Galaxy IC -- {num_total} particles")
    print(
        f"  disk:  {params.num_disk:6d} particles, "
        f"M = {params.disk_mass:.4f}"
    )
    print(
        f"  bulge: {params.num_bulge:6d} particles, "
        f"M = {params.bulge_mass:.4f}"
    )
    print(
        f"  halo:  {params.num_halo:6d} particles, "
        f"M = {params.halo_mass:.4f}"
    )
    print(f"  total mass (sum of particle masses) = {masses.sum():.4f}")
    print(
        f"  expected total = {params.disk_mass + params.bulge_mass + params.halo_mass:.4f}"
    )
    print(f"  position bounding box (max |x_i|) = {np.abs(positions).max():.4f}")

    output_path = "data/galaxy_positions.npz"
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    np.savez(output_path, positions=positions, masses=masses)
    print(f"\nSaved to {output_path}")


if __name__ == "__main__":
    main()
