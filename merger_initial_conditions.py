"""
Two-galaxy merger initial condition for the section 3.0.5 exercise.

Builds two copies of the validated single-galaxy equilibrium model and
places them on a Keplerian two-body encounter orbit, reproducing the
manual's MERGER_000 setup: both disks spin in the same sense, one is
inclined by 30 degrees, and the pair is put on a bound (parabolic)
orbit with a close pericentre so that tidal tails form and the galaxies
coalesce into a spheroidal remnant.

Orbit construction (relative two-body problem, total point mass
M = M_1 + M_2, G = 1):

    parabolic speed at separation r_0 :  v_0 = sqrt(2 G M / r_0)
    specific ang. mom. for pericentre r_p (parabola, e = 1):
                                          l   = sqrt(2 G M r_p)
    tangential component at r_0        :  v_t = l / r_0
    radial (approaching) component     :  v_r = -sqrt(v_0^2 - v_t^2)

The two equal-mass galaxies are placed symmetrically about the centre
of mass, so the system carries zero net momentum.

References: Toomre & Toomre (1972); Barnes (1988); Hernquist (1992) --
standard disc-galaxy merger orbit setup.

CMPH Project 3 -- Zoutong Shen / Zhaoyang Chu, 2026.
"""

import argparse
from typing import Optional, Tuple

import numpy as np

from initial_conditions import GalaxyParams
from velocities import make_galaxy_initial_conditions


def rotation_matrix(angle_radians: float, axis: str) -> np.ndarray:
    """
    Return the 3x3 active rotation matrix about a Cartesian axis.

    Args:
        angle_radians: rotation angle.
        axis:          one of "x", "y", "z".

    Returns:
        Rotation matrix R such that ``(R @ v.T).T`` rotates row vectors.
    """
    cosine = np.cos(angle_radians)
    sine = np.sin(angle_radians)
    if axis == "x":
        return np.array([[1.0, 0.0, 0.0],
                         [0.0, cosine, -sine],
                         [0.0, sine, cosine]])
    if axis == "y":
        return np.array([[cosine, 0.0, sine],
                         [0.0, 1.0, 0.0],
                         [-sine, 0.0, cosine]])
    if axis == "z":
        return np.array([[cosine, -sine, 0.0],
                         [sine, cosine, 0.0],
                         [0.0, 0.0, 1.0]])
    raise ValueError(f"axis must be 'x', 'y' or 'z', got {axis!r}")


def two_body_encounter_velocity(
    total_mass: float, separation: float, pericentre: float
) -> Tuple[float, float]:
    """
    Radial and tangential relative-velocity components of a parabolic
    two-body orbit, evaluated at the initial separation.

    Returns:
        (radial_velocity, tangential_velocity); radial_velocity < 0
        (the galaxies are approaching).
    """
    parabolic_speed = np.sqrt(2.0 * total_mass / separation)
    specific_angular_momentum = np.sqrt(2.0 * total_mass * pericentre)
    tangential_velocity = specific_angular_momentum / separation
    radial_velocity = -np.sqrt(
        max(parabolic_speed ** 2 - tangential_velocity ** 2, 0.0)
    )
    return radial_velocity, tangential_velocity


def make_merger_initial_conditions(
    params: GalaxyParams,
    *,
    params_2: Optional[GalaxyParams] = None,
    seed: int = 0,
    separation: float = 30.0,
    pericentre: float = 5.0,
    inclination_degrees: float = 30.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Build the full two-galaxy merger IC.

    Galaxy 1 uses ``params`` and lies in the x-y (orbital) plane; galaxy 2
    uses ``params_2`` (defaulting to ``params`` for the equal-mass case)
    and is inclined by ``inclination_degrees`` about the y axis (same spin
    sense). The pair is placed on a parabolic two-body encounter orbit
    with the given initial ``separation`` and ``pericentre``. For unequal
    masses the centre-of-mass split of position and velocity is mass-ratio
    weighted, so the system carries exactly zero net momentum and the
    relative orbit matches the prescribed parabolic trajectory.

    Returns:
        positions, velocities, masses, galaxy_id  (galaxy_id is 0 for the
        first galaxy's particles, 1 for the second's).
    """
    if params_2 is None:
        params_2 = params

    positions_1, velocities_1, masses_1 = make_galaxy_initial_conditions(
        params=params, seed=seed
    )
    positions_2, velocities_2, masses_2 = make_galaxy_initial_conditions(
        params=params_2, seed=seed + 7919
    )

    # Incline the second galaxy (positions and velocities rotate together).
    tilt = rotation_matrix(np.radians(inclination_degrees), "y")
    positions_2 = positions_2 @ tilt.T
    velocities_2 = velocities_2 @ tilt.T

    mass_1 = float(masses_1.sum())
    mass_2 = float(masses_2.sum())
    total_mass = mass_1 + mass_2
    fraction_1 = mass_2 / total_mass  # galaxy 1 sits at +fraction_1 * separation
    fraction_2 = mass_1 / total_mass  # galaxy 2 sits at -fraction_2 * separation
    radial_velocity, tangential_velocity = two_body_encounter_velocity(
        total_mass, separation, pericentre
    )

    # Each galaxy is an independent finite-N realisation, so it carries a
    # small spurious bulk position/velocity offset (~1-2% of v at default
    # N). Remove each galaxy's own mass-weighted centroid and mean velocity
    # before placing it on the orbit, so the only motion is the prescribed
    # two-body orbit and the system has exactly zero net momentum.
    for positions, velocities, component_masses in (
        (positions_1, velocities_1, masses_1),
        (positions_2, velocities_2, masses_2),
    ):
        total = component_masses.sum()
        positions -= (
            positions * component_masses[:, None]
        ).sum(axis=0) / total
        velocities -= (
            velocities * component_masses[:, None]
        ).sum(axis=0) / total

    # Place each galaxy at its mass-weighted centre-of-mass offset and give
    # it the mass-weighted share of the relative velocity, so the net
    # momentum is exactly zero. Equal masses reduce to a 50/50 split.
    relative_velocity = np.array(
        [radial_velocity, tangential_velocity, 0.0]
    )
    positions_1 = positions_1 + np.array([fraction_1 * separation, 0.0, 0.0])
    positions_2 = positions_2 - np.array([fraction_2 * separation, 0.0, 0.0])
    velocities_1 = velocities_1 + fraction_1 * relative_velocity
    velocities_2 = velocities_2 - fraction_2 * relative_velocity

    positions = np.vstack([positions_1, positions_2])
    velocities = np.vstack([velocities_1, velocities_2])
    masses = np.concatenate([masses_1, masses_2])
    galaxy_id = np.concatenate([
        np.zeros(len(masses_1), dtype=np.int8),
        np.ones(len(masses_2), dtype=np.int8),
    ])
    return positions, velocities, masses, galaxy_id


def main() -> None:
    """Build and save a default merger IC, print orbit diagnostics."""
    parser = argparse.ArgumentParser(
        description="Build the two-galaxy merger initial condition."
    )
    parser.add_argument("--n-disk", type=int, default=6000)
    parser.add_argument("--n-bulge", type=int, default=2000)
    parser.add_argument("--n-halo", type=int, default=12000)
    parser.add_argument("--separation", type=float, default=30.0)
    parser.add_argument("--pericentre", type=float, default=5.0)
    parser.add_argument("--inclination", type=float, default=30.0)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", type=str, default="data/merger_ic.npz")
    args = parser.parse_args()

    params = GalaxyParams(
        num_disk=args.n_disk,
        num_bulge=args.n_bulge,
        num_halo=args.n_halo,
    )
    positions, velocities, masses, galaxy_id = make_merger_initial_conditions(
        params,
        seed=args.seed,
        separation=args.separation,
        pericentre=args.pericentre,
        inclination_degrees=args.inclination,
    )
    net_momentum = np.linalg.norm((velocities * masses[:, None]).sum(axis=0))
    print(f"Merger IC: {positions.shape[0]} particles "
          f"({(galaxy_id == 0).sum()} + {(galaxy_id == 1).sum()})")
    print(f"  separation = {args.separation}, pericentre = {args.pericentre}, "
          f"inclination = {args.inclination} deg")
    print(f"  |net momentum| / total mass = "
          f"{net_momentum / masses.sum():.3e}")
    np.savez_compressed(
        args.out, positions=positions, velocities=velocities,
        masses=masses, galaxy_id=galaxy_id,
    )
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
