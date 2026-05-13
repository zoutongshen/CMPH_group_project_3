"""
Two-body Kepler validation for the gravity + integrator stack.

Two equal masses on a circular orbit. With G = 1, equal masses m, and orbit
radius r about the centre of mass (separation d = 2r), the centripetal
balance gives

    v^2 = G m / (2 d)        and        T = 2*pi*sqrt(d^3 / (2 G m))

For m = 1 and d = 1, this is v = sqrt(1/2) and T = pi*sqrt(2) ~= 4.443.

Three checks are run:
    1. Periodic return — after k full periods the particles are close to their
       initial positions.
    2. Energy conservation — fractional drift below 1e-3 over 10 periods.
    3. Angular-momentum conservation — fractional drift below 1e-6.

Run as a script:
    python kepler_test.py

CMPH Project 3 — Zoutong Shen / Zhaoyang Chu, 2026.
"""

import numpy as np

from diagnostics import angular_momentum, total_energy
from integrator import run_simulation


def two_body_circular_initial_conditions() -> (
    "tuple[np.ndarray, np.ndarray, np.ndarray, float]"
):
    """
    Build the initial state for the two-body circular test.

    Returns:
        positions: shape (2, 3).
        velocities: shape (2, 3).
        masses: shape (2,).
        period: analytic orbital period.
    """
    positions = np.array(
        [[+0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]], dtype=np.float64
    )
    speed = np.sqrt(0.5)
    velocities = np.array(
        [[0.0, +speed, 0.0], [0.0, -speed, 0.0]], dtype=np.float64
    )
    masses = np.array([1.0, 1.0], dtype=np.float64)
    period = np.pi * np.sqrt(2.0)
    return positions, velocities, masses, period


def main() -> None:
    """Run the Kepler validation suite and print pass/fail per check."""
    positions, velocities, masses, period = two_body_circular_initial_conditions()
    softening = 1.0e-3
    num_periods = 10
    timestep = 1.0e-3 * period
    num_steps = int(np.round(num_periods * period / timestep))

    print(f"Two-body circular Kepler test")
    print(f"  Analytic period T = {period:.6f}")
    print(f"  dt / T = {timestep / period:.3e}, num_steps = {num_steps}")

    positions_history, velocities_history, snapshot_times = run_simulation(
        positions,
        velocities,
        masses,
        timestep=timestep,
        num_steps=num_steps,
        softening=softening,
        snapshot_interval=num_steps // 100,
        progress=False,
    )

    energies = np.array(
        [
            total_energy(
                positions_history[snap], velocities_history[snap], masses, softening
            )
            for snap in range(positions_history.shape[0])
        ]
    )
    angular_momenta = np.array(
        [
            angular_momentum(
                positions_history[snap], velocities_history[snap], masses
            )
            for snap in range(positions_history.shape[0])
        ]
    )

    energy_initial = energies[0]
    energy_drift = np.max(np.abs(energies - energy_initial) / np.abs(energy_initial))

    angular_initial_norm = np.linalg.norm(angular_momenta[0])
    angular_drift = np.max(
        np.linalg.norm(angular_momenta - angular_momenta[0], axis=1)
        / max(angular_initial_norm, 1.0e-30)
    )

    final_positions = positions_history[-1]
    initial_positions = positions_history[0]
    return_residual = np.max(np.abs(final_positions - initial_positions))

    energy_tolerance = 1.0e-3
    angular_tolerance = 1.0e-6
    return_tolerance = 1.0e-2

    energy_passed = energy_drift < energy_tolerance
    angular_passed = angular_drift < angular_tolerance
    return_passed = return_residual < return_tolerance

    print()
    print(f"  Energy drift           = {energy_drift:.3e}   "
          f"(tol {energy_tolerance:.0e})   "
          f"{'PASS' if energy_passed else 'FAIL'}")
    print(f"  Angular-momentum drift = {angular_drift:.3e}   "
          f"(tol {angular_tolerance:.0e})   "
          f"{'PASS' if angular_passed else 'FAIL'}")
    print(f"  Periodic-return residual after {num_periods} T = "
          f"{return_residual:.3e}   "
          f"(tol {return_tolerance:.0e})   "
          f"{'PASS' if return_passed else 'FAIL'}")
    print()

    if not (energy_passed and angular_passed and return_passed):
        raise SystemExit("Kepler test FAILED — see drifts above.")
    print("Kepler test PASSED.")


if __name__ == "__main__":
    main()
