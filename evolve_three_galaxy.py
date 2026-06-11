"""
Driver for three-galaxy merger runs (multi-galaxy extension).

Three identical galaxies are placed at the vertices of an equilateral
triangle and given sub-Keplerian tangential velocities, so they spiral
in under mutual gravity and eventually coalesce into a single remnant.

Run as a script:
    python evolve_three_galaxy.py
    python evolve_three_galaxy.py --tangential-fraction 0.7 --tstop 400 \\
        --out data/three_galaxy_tf07.npz

CMPH Project 3 -- extension/multi-galaxy.
"""

import os
import argparse
import time

import numpy as np

from initial_conditions import GalaxyParams
from integrator import run_simulation
from three_galaxy_initial_conditions import (
    make_three_galaxy_initial_conditions,
)


def run_three_galaxy(
    params: GalaxyParams,
    *,
    timestep: float,
    tstop: float,
    softening: float,
    seed: int,
    ring_radius: float,
    tangential_fraction: float,
    inclination_degrees: float,
    num_dumps: int,
) -> dict:
    """Build the three-galaxy IC and integrate."""
    positions, velocities, masses, galaxy_id = (
        make_three_galaxy_initial_conditions(
            params,
            seed=seed,
            ring_radius=ring_radius,
            tangential_fraction=tangential_fraction,
            inclination_degrees=inclination_degrees,
        )
    )

    num_steps = int(round(tstop / timestep))
    snapshot_interval = max(1, num_steps // num_dumps)
    print(f"Built three-galaxy IC: {positions.shape[0]} particles "
          f"(3 x {positions.shape[0] // 3})")
    print(f"Integrating: dt = {timestep}, tstop = {tstop}, "
          f"steps = {num_steps}, eps = {softening}")

    start = time.time()
    positions_history, velocities_history, snapshot_times = run_simulation(
        positions, velocities, masses,
        timestep=timestep, num_steps=num_steps, softening=softening,
        snapshot_interval=snapshot_interval, progress=True,
    )
    print(f"Finished in {time.time() - start:.1f} s")

    return {
        "positions_history": positions_history.astype(np.float32),
        "velocities_history": velocities_history.astype(np.float32),
        "masses": masses,
        "galaxy_id": galaxy_id,
        "snapshot_times": snapshot_times,
        "num_disk": params.num_disk,
        "num_bulge": params.num_bulge,
        "num_halo": params.num_halo,
        "softening": softening,
        "timestep": timestep,
        "ring_radius": ring_radius,
        "tangential_fraction": tangential_fraction,
        "inclination_degrees": inclination_degrees,
        "seed": seed,
    }


def parse_command_line() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description="Evolve a three-galaxy merger (equilateral triangle)."
    )
    parser.add_argument("--n-disk", type=int, default=6000)
    parser.add_argument("--n-bulge", type=int, default=2000)
    parser.add_argument("--n-halo", type=int, default=12000)
    parser.add_argument("--tstop", type=float, default=400.0)
    parser.add_argument("--dt", type=float, default=0.125)
    parser.add_argument("--eps", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--ring-radius", type=float, default=20.0)
    parser.add_argument("--tangential-fraction", type=float, default=0.5)
    parser.add_argument("--inclination", type=float, default=30.0)
    parser.add_argument("--num-dumps", type=int, default=200)
    parser.add_argument("--out", type=str,
                        default="data/three_galaxy.npz")
    return parser.parse_args()


def main() -> None:
    """Top-level entry point."""
    args = parse_command_line()
    params = GalaxyParams(
        num_disk=args.n_disk,
        num_bulge=args.n_bulge,
        num_halo=args.n_halo,
    )
    result = run_three_galaxy(
        params,
        timestep=args.dt,
        tstop=args.tstop,
        softening=args.eps,
        seed=args.seed,
        ring_radius=args.ring_radius,
        tangential_fraction=args.tangential_fraction,
        inclination_degrees=args.inclination,
        num_dumps=args.num_dumps,
    )
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    np.savez_compressed(args.out, **result)
    print(f"Saved trajectory to {args.out}")


if __name__ == "__main__":
    main()
