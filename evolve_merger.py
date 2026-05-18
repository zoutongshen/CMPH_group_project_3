"""
Driver for the section 3.0.5 exercise: collide two disk galaxies and
dump the trajectory for remnant analysis and an animation.

Builds the two-galaxy encounter initial condition
(`merger_initial_conditions.py`) and integrates it with the same fixed
leapfrog step and Plummer softening used for the isolated disk. The run
is carried to t ~ 200 (~2.6 Gyr) so that the galaxies complete the first
passage, raise tidal tails, sink by dynamical friction and coalesce into
a single relaxed remnant. The per-particle galaxy_id is stored so the
analysis can show which material ends up where.

Run as a script:
    python evolve_merger.py                       # 2 x 20000, eps = 0.1
    python evolve_merger.py --n-disk 24000 --n-bulge 8000 --n-halo 48000 \\
        --out data/merger_N160k_eps01.npz

CMPH Project 3 -- Zoutong Shen / Zhaoyang Chu, 2026.
"""

import argparse
import time

import numpy as np

from initial_conditions import GalaxyParams
from integrator import run_simulation
from merger_initial_conditions import make_merger_initial_conditions


def run_merger(
    params: GalaxyParams,
    *,
    timestep: float,
    tstop: float,
    softening: float,
    seed: int,
    separation: float,
    pericentre: float,
    inclination_degrees: float,
    num_dumps: int,
) -> dict:
    """
    Build the two-galaxy IC and integrate it to coalescence.

    Args:
        params:              per-galaxy structural / particle-count
                             parameters (each galaxy uses these).
        timestep:            leapfrog step dt, keyword-only.
        tstop:               end time in code units, keyword-only.
        softening:           Plummer softening length eps, keyword-only.
        seed:                RNG seed for galaxy 1 (galaxy 2 is offset),
                             keyword-only.
        separation:          initial galaxy separation, keyword-only.
        pericentre:          orbit pericentre, keyword-only.
        inclination_degrees: tilt of the second disk, keyword-only.
        num_dumps:           number of snapshots after t = 0.

    Returns:
        Dict ready to be saved (trajectory, masses, galaxy_id, metadata).
    """
    positions, velocities, masses, galaxy_id = (
        make_merger_initial_conditions(
            params,
            seed=seed,
            separation=separation,
            pericentre=pericentre,
            inclination_degrees=inclination_degrees,
        )
    )
    num_steps = int(round(tstop / timestep))
    snapshot_interval = max(1, num_steps // num_dumps)

    print(f"Built merger IC: {positions.shape[0]} particles "
          f"(2 x {positions.shape[0] // 2})")
    print(f"Integrating: dt = {timestep}, tstop = {tstop}, "
          f"steps = {num_steps}, eps = {softening}, "
          f"dumps = {num_steps // snapshot_interval + 1}")

    start = time.time()
    positions_history, velocities_history, snapshot_times = run_simulation(
        positions,
        velocities,
        masses,
        timestep=timestep,
        num_steps=num_steps,
        softening=softening,
        snapshot_interval=snapshot_interval,
        progress=True,
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
        "separation": separation,
        "pericentre": pericentre,
        "inclination_degrees": inclination_degrees,
        "seed": seed,
    }


def parse_command_line() -> argparse.Namespace:
    """Parse CLI arguments for the merger run."""
    parser = argparse.ArgumentParser(
        description="Evolve a two-galaxy merger (manual 3.0.5)."
    )
    parser.add_argument("--n-disk", type=int, default=6000)
    parser.add_argument("--n-bulge", type=int, default=2000)
    parser.add_argument("--n-halo", type=int, default=12000)
    parser.add_argument("--tstop", type=float, default=300.0)
    parser.add_argument("--dt", type=float, default=0.125)
    parser.add_argument("--eps", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--separation", type=float, default=30.0)
    parser.add_argument("--pericentre", type=float, default=5.0)
    parser.add_argument("--inclination", type=float, default=30.0)
    parser.add_argument("--num-dumps", type=int, default=200)
    parser.add_argument(
        "--out", type=str, default="data/merger_N40k_eps01.npz",
        help="output .npz path for the trajectory",
    )
    return parser.parse_args()


def main() -> None:
    """Top-level entry point."""
    args = parse_command_line()
    params = GalaxyParams(
        num_disk=args.n_disk,
        num_bulge=args.n_bulge,
        num_halo=args.n_halo,
    )
    result = run_merger(
        params,
        timestep=args.dt,
        tstop=args.tstop,
        softening=args.eps,
        seed=args.seed,
        separation=args.separation,
        pericentre=args.pericentre,
        inclination_degrees=args.inclination,
        num_dumps=args.num_dumps,
    )
    np.savez_compressed(args.out, **result)
    print(f"Saved trajectory to {args.out}")


if __name__ == "__main__":
    main()
