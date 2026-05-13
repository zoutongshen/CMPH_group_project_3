"""
Conserved-quantity diagnostics for the N-body simulation.

Provides kinetic energy, gravitational potential energy (Plummer-softened),
total energy, and angular momentum. These are used both as numerical-quality
checks on the integrator (energy drift, angular-momentum drift) and as
physical observables for the merger remnant.

Units: G = 1, consistent with `forces.py` and the practicum manual.

CMPH Project 3 — Zoutong Shen / Zhaoyang Chu, 2026.
"""

import numpy as np
from numba import njit, prange


@njit(fastmath=True, cache=True)
def kinetic_energy(velocities: np.ndarray, masses: np.ndarray) -> float:
    """
    Return total kinetic energy of the particle system.

    KE = sum_i 0.5 m_i |v_i|^2

    Args:
        velocities: array of shape (num_particles, 3).
        masses:     array of shape (num_particles,).

    Returns:
        Scalar kinetic energy.
    """
    num_particles = velocities.shape[0]
    total_kinetic = 0.0
    for i in range(num_particles):
        speed_squared = (
            velocities[i, 0] ** 2
            + velocities[i, 1] ** 2
            + velocities[i, 2] ** 2
        )
        total_kinetic += 0.5 * masses[i] * speed_squared
    return total_kinetic


@njit(parallel=True, fastmath=True, cache=True)
def potential_energy(
    positions: np.ndarray, masses: np.ndarray, softening: float
) -> float:
    """
    Return total gravitational potential energy with Plummer softening.

    U = -sum_{i<j} m_i m_j / sqrt(|x_i - x_j|^2 + eps^2),  with G = 1.

    Each pair is counted once. The outer loop is parallelised via Numba's
    prange reduction on `total_potential`.

    Args:
        positions: array of shape (num_particles, 3).
        masses:    array of shape (num_particles,).
        softening: Plummer softening length.

    Returns:
        Scalar potential energy.
    """
    num_particles = positions.shape[0]
    softening_squared = softening * softening
    total_potential = 0.0
    for i in prange(num_particles):
        partial = 0.0
        xi = positions[i, 0]
        yi = positions[i, 1]
        zi = positions[i, 2]
        mi = masses[i]
        for j in range(i + 1, num_particles):
            dx = positions[j, 0] - xi
            dy = positions[j, 1] - yi
            dz = positions[j, 2] - zi
            distance_squared = (
                dx * dx + dy * dy + dz * dz + softening_squared
            )
            partial -= mi * masses[j] / np.sqrt(distance_squared)
        total_potential += partial
    return total_potential


def total_energy(
    positions: np.ndarray,
    velocities: np.ndarray,
    masses: np.ndarray,
    softening: float,
) -> float:
    """
    Return total mechanical energy E = KE + U.

    Wrapper combining `kinetic_energy` and `potential_energy`.

    Args:
        positions:  array of shape (num_particles, 3).
        velocities: array of shape (num_particles, 3).
        masses:     array of shape (num_particles,).
        softening:  Plummer softening length.

    Returns:
        Scalar total energy.
    """
    return kinetic_energy(velocities, masses) + potential_energy(
        positions, masses, softening
    )


@njit(fastmath=True, cache=True)
def angular_momentum(
    positions: np.ndarray, velocities: np.ndarray, masses: np.ndarray
) -> np.ndarray:
    """
    Return total angular momentum vector L = sum_i m_i (x_i x v_i).

    Args:
        positions:  array of shape (num_particles, 3).
        velocities: array of shape (num_particles, 3).
        masses:     array of shape (num_particles,).

    Returns:
        Array of shape (3,) — total angular momentum vector.
    """
    num_particles = positions.shape[0]
    total_angular = np.zeros(3)
    for i in range(num_particles):
        mass = masses[i]
        total_angular[0] += mass * (
            positions[i, 1] * velocities[i, 2]
            - positions[i, 2] * velocities[i, 1]
        )
        total_angular[1] += mass * (
            positions[i, 2] * velocities[i, 0]
            - positions[i, 0] * velocities[i, 2]
        )
        total_angular[2] += mass * (
            positions[i, 0] * velocities[i, 1]
            - positions[i, 1] * velocities[i, 0]
        )
    return total_angular
