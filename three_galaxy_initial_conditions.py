"""
Three-galaxy initial condition for the multi-galaxy merger extension.

Three identical galaxies are placed at the vertices of an equilateral
triangle in the x-y plane (one in the orbital plane, two inclined by
+/- inclination_degrees about the y axis for visual variety). Each
galaxy is given a tangential velocity equal to ``tangential_fraction``
times the circular-orbit velocity for the equilateral configuration, so
that with the default fraction = 0.5 they are sub-Keplerian and spiral
inward to a triple merger.

The equilateral-triangle circular-orbit velocity (G = 1):
    v_circ = sqrt( M_galaxy / (R * sqrt(3)) )

derived from the net gravitational force at one vertex from the other
two equal masses (force magnitude G M^2 / (R^2 sqrt(3)) toward the COM).

CMPH Project 3 -- extension/multi-galaxy.
"""

import os
import argparse
from typing import Tuple

import numpy as np

from initial_conditions import GalaxyParams
from merger_initial_conditions import rotation_matrix
from velocities import make_galaxy_initial_conditions


def make_three_galaxy_initial_conditions(
    params: GalaxyParams,
    *,
    seed: int = 0,
    ring_radius: float = 20.0,
    tangential_fraction: float = 0.5,
    inclination_degrees: float = 30.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Build the three-galaxy IC.

    Args:
        params:              shared galaxy structural parameters
                             (all three galaxies are identical).
        seed:                base RNG seed (each galaxy is offset by k*7919).
        ring_radius:         COM-to-vertex distance R for each galaxy.
        tangential_fraction: factor multiplying v_circ for the initial
                             tangential velocity (0 = rest, 1 = full
                             circular support, 0.5 = sub-Keplerian).
        inclination_degrees: tilt about y for galaxies 2 and 3
                             (galaxy 1 stays in the orbital plane).

    Returns:
        positions, velocities, masses, galaxy_id  (galaxy_id is 0, 1, 2
        for the three galaxies; same Cartesian layout as the two-body case).
    """
    galaxies = []
    for k in range(3):
        positions, velocities, masses = make_galaxy_initial_conditions(
            params=params, seed=seed + k * 7919
        )
        total_mass = masses.sum()
        positions -= (positions * masses[:, None]).sum(axis=0) / total_mass
        velocities -= (velocities * masses[:, None]).sum(axis=0) / total_mass
        galaxies.append([positions, velocities, masses])

    # Incline galaxies 2 and 3 in opposite senses for visual variety.
    tilt_positive = rotation_matrix(np.radians(inclination_degrees), "y")
    tilt_negative = rotation_matrix(np.radians(-inclination_degrees), "y")
    galaxies[1][0] = galaxies[1][0] @ tilt_positive.T
    galaxies[1][1] = galaxies[1][1] @ tilt_positive.T
    galaxies[2][0] = galaxies[2][0] @ tilt_negative.T
    galaxies[2][1] = galaxies[2][1] @ tilt_negative.T

    galaxy_mass = float(galaxies[0][2].sum())
    v_circular = np.sqrt(galaxy_mass / (ring_radius * np.sqrt(3.0)))
    v_tangential_magnitude = tangential_fraction * v_circular

    angles = np.array([0.0, 2.0 * np.pi / 3.0, 4.0 * np.pi / 3.0])
    for k, angle in enumerate(angles):
        position_offset = np.array([
            ring_radius * np.cos(angle),
            ring_radius * np.sin(angle),
            0.0,
        ])
        tangential_velocity = np.array([
            -v_tangential_magnitude * np.sin(angle),
            v_tangential_magnitude * np.cos(angle),
            0.0,
        ])
        galaxies[k][0] = galaxies[k][0] + position_offset
        galaxies[k][1] = galaxies[k][1] + tangential_velocity

    positions = np.vstack([g[0] for g in galaxies])
    velocities = np.vstack([g[1] for g in galaxies])
    masses = np.concatenate([g[2] for g in galaxies])
    galaxy_id = np.concatenate([
        np.full(len(galaxies[k][2]), k, dtype=np.int8) for k in range(3)
    ])
    return positions, velocities, masses, galaxy_id


def main() -> None:
    """Build and save a default three-galaxy IC, print diagnostics."""
    parser = argparse.ArgumentParser(
        description="Build the three-galaxy initial condition."
    )
    parser.add_argument("--n-disk", type=int, default=6000)
    parser.add_argument("--n-bulge", type=int, default=2000)
    parser.add_argument("--n-halo", type=int, default=12000)
    parser.add_argument("--ring-radius", type=float, default=20.0)
    parser.add_argument("--tangential-fraction", type=float, default=0.5)
    parser.add_argument("--inclination", type=float, default=30.0)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", type=str, default="data/three_galaxy_ic.npz")
    args = parser.parse_args()

    params = GalaxyParams(
        num_disk=args.n_disk,
        num_bulge=args.n_bulge,
        num_halo=args.n_halo,
    )
    positions, velocities, masses, galaxy_id = (
        make_three_galaxy_initial_conditions(
            params,
            seed=args.seed,
            ring_radius=args.ring_radius,
            tangential_fraction=args.tangential_fraction,
            inclination_degrees=args.inclination,
        )
    )
    net_momentum = np.linalg.norm(
        (velocities * masses[:, None]).sum(axis=0)
    )
    print(f"3-galaxy IC: {positions.shape[0]} particles "
          f"({(galaxy_id == 0).sum()} + {(galaxy_id == 1).sum()} + "
          f"{(galaxy_id == 2).sum()})")
    print(f"  ring_radius = {args.ring_radius}, "
          f"tangential_fraction = {args.tangential_fraction}, "
          f"inclination = {args.inclination} deg")
    print(f"  |net momentum| / total mass = "
          f"{net_momentum / masses.sum():.3e}")
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    np.savez_compressed(
        args.out, positions=positions, velocities=velocities,
        masses=masses, galaxy_id=galaxy_id,
    )
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
