"""
Gravitational accelerations for an N-body system with Plummer softening.

Implements equation (2.1) of the practicum manual (Naab 2006):

    F_i = -sum_{j != i} G m_i m_j (x_i - x_j) / (|x_i - x_j|^2 + eps^2)^(3/2)

We work in units where G = 1. The softening length suppresses the singularity
at zero separation and reduces spurious two-body relaxation from finite-N
particle noise.

CMPH Project 3 — Zoutong Shen / Zhaoyang Chu, 2026.
"""

import numpy as np
from numba import njit, prange


@njit(parallel=True, fastmath=True, cache=True)
def gravitational_accelerations(
    positions: np.ndarray, masses: np.ndarray, softening: float
) -> np.ndarray:
    """
    Return the gravitational acceleration on every particle (units with G = 1).

    The inner loop is parallelised over particles via Numba's prange. Plummer
    softening is folded into the squared distance so no branching is needed in
    the inner loop.

    Args:
        positions: array of shape (num_particles, 3) — particle positions.
        masses:    array of shape (num_particles,)   — particle masses.
        softening: Plummer softening length (same units as positions).

    Returns:
        Array of shape (num_particles, 3) — acceleration on each particle.
    """
    num_particles = positions.shape[0]
    accelerations = np.zeros((num_particles, 3))
    softening_squared = softening * softening

    for i in prange(num_particles):
        accel_x = 0.0
        accel_y = 0.0
        accel_z = 0.0
        xi = positions[i, 0]
        yi = positions[i, 1]
        zi = positions[i, 2]
        for j in range(num_particles):
            if j == i:
                continue
            dx = positions[j, 0] - xi
            dy = positions[j, 1] - yi
            dz = positions[j, 2] - zi
            distance_squared = dx * dx + dy * dy + dz * dz + softening_squared
            inverse_r_cubed = distance_squared ** (-1.5)
            mass_over_r3 = masses[j] * inverse_r_cubed
            accel_x += mass_over_r3 * dx
            accel_y += mass_over_r3 * dy
            accel_z += mass_over_r3 * dz
        accelerations[i, 0] = accel_x
        accelerations[i, 1] = accel_y
        accelerations[i, 2] = accel_z

    return accelerations
