"""
Driver for unequal-mass merger runs (multi-galaxy extension).

Galaxy 1 uses the default GalaxyParams. Galaxy 2 has all three component
masses (disk, bulge, halo) scaled by ``mass_ratio``; particle counts are
scaled proportionally so the per-particle mass stays constant. With
``mass_ratio = 1.0`` this reduces to the equal-mass merger.

Run as a script:
    python evolve_unequal_merger.py --mass-ratio 0.5
    python evolve_unequal_merger.py --mass-ratio 0.25 \\
        --out data/merger_unequal_1to4.npz

CMPH Project 3 -- extension/multi-galaxy.
"""

import os
import argparse
import time

import numpy as np

from initial_conditions import GalaxyParams
from integrator import run_simulation
from merger_initial_conditions import make_merger_initial_conditions


def scaled_galaxy_params(base: GalaxyParams, mass_ratio: float) -> GalaxyParams:
    """
    Return a GalaxyParams with all masses (and particle counts) scaled by
    ``mass_ratio``. Structural scale lengths are unchanged, so the galaxy
    is geometrically identical but dynamically lighter.
    """
    return GalaxyParams(
        disk_scale_length=base.disk_scale_length,
        disk_scale_height=base.disk_scale_height,
        disk_mass=base.disk_mass * mass_ratio,
        bulge_mass=base.bulge_mass * mass_ratio,
        bulge_scale=base.bulge_scale,
        halo_mass=base.halo_mass * mass_ratio,
        halo_core=base.halo_core,
        halo_cutoff=base.halo_cutoff,
        num_disk=max(1, int(round(base.num_disk * mass_ratio))),
        num_bulge=max(1, int(round(base.num_bulge * mass_ratio))),
        num_halo=max(1, int(round(base.num_halo * mass_ratio))),
    )


def run_unequal_merger(
    params_1: GalaxyParams,
    *,
    mass_ratio: float,
    timestep: float,
    tstop: float,
    softening: float,
    seed: int,
    separation: float,
    pericentre: float,
    inclination_degrees: float,
    num_dumps: int,
) -> dict:
    """Build the unequal-mass IC and integrate."""
    params_2 = scaled_galaxy_params(params_1, mass_ratio)
    positions, velocities, masses, galaxy_id = make_merger_initial_conditions(
        params_1,
        params_2=params_2,
        seed=seed,
        separation=separation,
        pericentre=pericentre,
        inclination_degrees=inclination_degrees,
    )

    num_steps = int(round(tstop / timestep))
    snapshot_interval = max(1, num_steps // num_dumps)
    print(f"Built unequal-mass IC: galaxy 1 = {(galaxy_id == 0).sum()} "
          f"particles (M = {masses[galaxy_id == 0].sum():.3f}), "
          f"galaxy 2 = {(galaxy_id == 1).sum()} particles "
          f"(M = {masses[galaxy_id == 1].sum():.3f}, ratio = {mass_ratio})")
    print(f"Integrating: dt = {timestep}, tstop = {tstop}, "
          f"steps = {num_steps}, eps = {softening}")

    start = time.time()
    positions_history, velocities_history, snapshot_times = run_simulation(
        positions, velocities, masses,
        timestep=timestep, num_steps=num_steps, softening=softening,
        snapshot_interval=snapshot_interval, progress=True,
    )
    print(f"Finished in {time.time() - start:.1f} s")

    # Reuse evolve_merger's output schema so existing analysis/animation
    # scripts work unchanged; expose mass_ratio in metadata.
    return {
        "positions_history": positions_history.astype(np.float32),
        "velocities_history": velocities_history.astype(np.float32),
        "masses": masses,
        "galaxy_id": galaxy_id,
        "snapshot_times": snapshot_times,
        "num_disk": params_1.num_disk,
        "num_bulge": params_1.num_bulge,
        "num_halo": params_1.num_halo,
        "softening": softening,
        "timestep": timestep,
        "separation": separation,
        "pericentre": pericentre,
        "inclination_degrees": inclination_degrees,
        "mass_ratio": mass_ratio,
        "seed": seed,
    }


def parse_command_line() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description="Evolve an unequal-mass two-galaxy merger."
    )
    parser.add_argument("--n-disk", type=int, default=6000)
    parser.add_argument("--n-bulge", type=int, default=2000)
    parser.add_argument("--n-halo", type=int, default=12000)
    parser.add_argument("--mass-ratio", type=float, default=0.5,
                        help="m_galaxy2 / m_galaxy1")
    parser.add_argument("--tstop", type=float, default=300.0)
    parser.add_argument("--dt", type=float, default=0.125)
    parser.add_argument("--eps", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--separation", type=float, default=30.0)
    parser.add_argument("--pericentre", type=float, default=5.0)
    parser.add_argument("--inclination", type=float, default=30.0)
    parser.add_argument("--num-dumps", type=int, default=200)
    parser.add_argument("--out", type=str,
                        default="data/merger_unequal.npz")
    return parser.parse_args()


def main() -> None:
    """Top-level entry point."""
    args = parse_command_line()
    params_1 = GalaxyParams(
        num_disk=args.n_disk,
        num_bulge=args.n_bulge,
        num_halo=args.n_halo,
    )
    result = run_unequal_merger(
        params_1,
        mass_ratio=args.mass_ratio,
        timestep=args.dt,
        tstop=args.tstop,
        softening=args.eps,
        seed=args.seed,
        separation=args.separation,
        pericentre=args.pericentre,
        inclination_degrees=args.inclination,
        num_dumps=args.num_dumps,
    )
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    np.savez_compressed(args.out, **result)
    print(f"Saved trajectory to {args.out}")


if __name__ == "__main__":
    main()
