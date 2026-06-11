"""
Driver for the section 3.0.4 exercise: evolve a single disk galaxy in
isolation and dump the trajectory for later analysis.

This reproduces what the manual's VINE run produces: starting from the
equilibrium initial condition, integrate to tstop = 100 with the fixed
leapfrog step dt = 0.125 and Plummer softening eps, writing 100 dump
snapshots (plus the t = 0 state). Unlike ``stability_test.py`` -- which
only keeps scalar diagnostics -- this saves the full particle arrays so
that the circular-velocity, enclosed-mass and surface-density profiles
(Q4/Q5) and the face-on / edge-on images (Q1/Q2) can be computed
afterwards from the same trajectory.

Run as a script:
    python evolve_disk.py                       # manual N = 20000, eps = 0.1
    python evolve_disk.py --eps 1e-4 --out data/evolve_N20k_eps1em4.npz
    python evolve_disk.py --n-disk 24000 --n-bulge 8000 --n-halo 48000 \\
        --out data/evolve_N80k_eps01.npz

The output .npz holds positions_history, velocities_history (float32),
masses, snapshot_times, the disk/bulge/halo split and the run metadata.

CMPH Project 3 -- Zoutong Shen / Zhaoyang Chu, 2026.
"""

import os
import argparse
import time

import numpy as np

from initial_conditions import GalaxyParams
from integrator import run_simulation
from velocities import make_galaxy_initial_conditions


def run_isolated_disk(
    params: GalaxyParams,
    *,
    timestep: float,
    tstop: float,
    softening: float,
    seed: int,
    num_dumps: int,
) -> dict:
    """
    Build one galaxy IC and integrate it in isolation.

    Args:
        params:    galaxy structural / particle-count parameters.
        timestep:  leapfrog step dt (manual: 0.125), keyword-only.
        tstop:     end time in code units (manual: 100), keyword-only.
        softening: Plummer softening length eps, keyword-only.
        seed:      RNG seed for the IC sampling, keyword-only.
        num_dumps: number of snapshots after t = 0 (manual: 100).

    Returns:
        Dict with positions_history, velocities_history, masses,
        snapshot_times and the component split, ready to be saved.
    """
    positions, velocities, masses = make_galaxy_initial_conditions(
        params=params, seed=seed
    )
    num_steps = int(round(tstop / timestep))
    snapshot_interval = max(1, num_steps // num_dumps)

    print(f"Built IC: {positions.shape[0]} particles "
          f"(disk {params.num_disk}, bulge {params.num_bulge}, "
          f"halo {params.num_halo})")
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
        "snapshot_times": snapshot_times,
        "num_disk": params.num_disk,
        "num_bulge": params.num_bulge,
        "num_halo": params.num_halo,
        "softening": softening,
        "timestep": timestep,
        "seed": seed,
    }


def parse_command_line() -> argparse.Namespace:
    """Parse CLI arguments for the isolated-disk evolution run."""
    parser = argparse.ArgumentParser(
        description="Evolve a single isolated disk galaxy (manual 3.0.4)."
    )
    parser.add_argument("--n-disk", type=int, default=6000)
    parser.add_argument("--n-bulge", type=int, default=2000)
    parser.add_argument("--n-halo", type=int, default=12000)
    parser.add_argument("--tstop", type=float, default=100.0)
    parser.add_argument("--dt", type=float, default=0.125)
    parser.add_argument("--eps", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--num-dumps", type=int, default=100)
    parser.add_argument(
        "--out", type=str, default="data/evolve_N20k_eps01.npz",
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
    result = run_isolated_disk(
        params,
        timestep=args.dt,
        tstop=args.tstop,
        softening=args.eps,
        seed=args.seed,
        num_dumps=args.num_dumps,
    )
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    np.savez_compressed(args.out, **result)
    print(f"Saved trajectory to {args.out}")


if __name__ == "__main__":
    main()
