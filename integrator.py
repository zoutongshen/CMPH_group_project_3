"""
Velocity-Verlet (leapfrog) integrator for the N-body equations of motion.

Implements equation (2.12) of the practicum manual in the kick-drift-kick
form, which is time-reversible and symplectic — energy drift is bounded
over arbitrarily long runs at fixed dt, rather than growing secularly.

Provides:
    velocity_verlet_step — advance one step in place
    run_simulation       — main driver, returns trajectory snapshots

CMPH Project 3 — Zoutong Shen / Zhaoyang Chu, 2026.
"""

from typing import Optional, Tuple

import numpy as np
from numba import njit

from forces import gravitational_accelerations


@njit(cache=True)
def velocity_verlet_step(
    positions: np.ndarray,
    velocities: np.ndarray,
    accelerations: np.ndarray,
    masses: np.ndarray,
    timestep: float,
    softening: float,
) -> np.ndarray:
    """
    Advance the system by one timestep using kick-drift-kick velocity Verlet.

    Modifies `positions` and `velocities` in place. The new accelerations are
    returned so the caller can reuse them as the starting kick for the next
    step (avoids recomputing forces twice per step).

    Args:
        positions:     array of shape (num_particles, 3), mutated in place.
        velocities:    array of shape (num_particles, 3), mutated in place.
        accelerations: array of shape (num_particles, 3), the acceleration
                       evaluated at the *current* positions (start of step).
        masses:        array of shape (num_particles,).
        timestep:      integration timestep dt.
        softening:     Plummer softening length.

    Returns:
        Array of shape (num_particles, 3) — accelerations at the new positions.
    """
    half_dt = 0.5 * timestep

    velocities += half_dt * accelerations
    positions += timestep * velocities

    new_accelerations = gravitational_accelerations(positions, masses, softening)

    velocities += half_dt * new_accelerations
    return new_accelerations


def run_simulation(
    positions: np.ndarray,
    velocities: np.ndarray,
    masses: np.ndarray,
    *,
    timestep: float,
    num_steps: int,
    softening: float,
    snapshot_interval: Optional[int] = None,
    progress: bool = False,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Run the N-body simulation for `num_steps` steps.

    Snapshots of (positions, velocities) are recorded every `snapshot_interval`
    steps (and at t=0). If `snapshot_interval` is None, only the initial and
    final states are stored.

    Args:
        positions:         initial positions, shape (num_particles, 3).
        velocities:        initial velocities, shape (num_particles, 3).
        masses:            shape (num_particles,).
        timestep:          dt for the integrator (keyword-only).
        num_steps:         number of timesteps to run (keyword-only).
        softening:         Plummer softening length (keyword-only).
        snapshot_interval: store a snapshot every this many steps. Default:
                           initial + final only.
        progress:          show a tqdm progress bar.

    Returns:
        positions_history:  shape (num_snapshots, num_particles, 3).
        velocities_history: shape (num_snapshots, num_particles, 3).
        snapshot_times:     shape (num_snapshots,) — simulation times.
    """
    positions = positions.copy()
    velocities = velocities.copy()

    if snapshot_interval is None:
        snapshot_interval = num_steps

    num_snapshots = num_steps // snapshot_interval + 1
    positions_history = np.empty((num_snapshots, *positions.shape))
    velocities_history = np.empty((num_snapshots, *velocities.shape))
    snapshot_times = np.empty(num_snapshots)

    positions_history[0] = positions
    velocities_history[0] = velocities
    snapshot_times[0] = 0.0

    accelerations = gravitational_accelerations(positions, masses, softening)

    step_iterator = range(1, num_steps + 1)
    if progress:
        from tqdm import tqdm
        step_iterator = tqdm(step_iterator)

    snapshot_index = 0
    for step in step_iterator:
        accelerations = velocity_verlet_step(
            positions, velocities, accelerations, masses, timestep, softening
        )
        if step % snapshot_interval == 0:
            snapshot_index += 1
            positions_history[snapshot_index] = positions
            velocities_history[snapshot_index] = velocities
            snapshot_times[snapshot_index] = step * timestep

    return positions_history, velocities_history, snapshot_times
