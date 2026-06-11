"""
Validation tests for the tabulated GalaxyPotential.

Each test checks one analytic / asymptotic property of the manual's model.

Run as a script:
    python test_potential.py

CMPH Project 3 -- Zoutong Shen / Zhaoyang Chu, 2026.
"""

import numpy as np

from initial_conditions import GalaxyParams
from potential import GalaxyPotential


def check_total_enclosed_mass_asymptote() -> None:
    """Total M(<r) at the outer grid edge should equal M_d + M_b + M_h."""
    params = GalaxyParams()
    potential = GalaxyPotential(params)

    expected = params.disk_mass + params.bulge_mass + params.halo_mass
    obtained = potential.enclosed_mass(potential._radii[-1])
    relative_error = abs(obtained - expected) / expected

    tolerance = 0.01
    assert relative_error < tolerance, (
        f"asymptotic total mass {obtained:.4f} vs expected {expected:.4f}"
    )
    print(
        f"  total enclosed mass at r_max: {obtained:.4f} "
        f"(expected {expected:.4f}, rel.err = {relative_error:.3e}) -- PASS"
    )


def check_bulge_alone_against_analytic() -> None:
    """
    Build a potential with only the bulge active and verify the empirical
    M_bulge(<r) matches the analytic Hernquist formula M_b r^2 / (r+a)^2.
    """
    params = GalaxyParams(disk_mass=1.0e-12, halo_mass=1.0e-12)
    potential = GalaxyPotential(params)

    bulge_scale = params.bulge_scale
    test_radii = np.geomspace(1.0e-2, 10.0, 50)
    analytic = params.bulge_mass * test_radii ** 2 / (test_radii + bulge_scale) ** 2
    tabulated = potential.enclosed_mass(test_radii)

    relative_error = np.max(np.abs(tabulated - analytic) / (analytic + 1e-30))
    tolerance = 0.02
    assert relative_error < tolerance, (
        f"max relative error {relative_error:.3e} exceeds {tolerance}"
    )
    print(
        f"  Hernquist bulge M(<r) vs analytic: max rel.err = "
        f"{relative_error:.3e} (tol {tolerance}) -- PASS"
    )


def check_bulge_potential_derivative_against_analytic() -> None:
    """
    dPhi_bulge/dr should match the analytic Hernquist derivative
    d/dr [-G M_b / (r+a)] = G M_b / (r+a)^2.
    """
    params = GalaxyParams(disk_mass=1.0e-12, halo_mass=1.0e-12)
    potential = GalaxyPotential(params)

    bulge_scale = params.bulge_scale
    test_radii = np.geomspace(0.01, 10.0, 50)
    analytic = params.bulge_mass / (test_radii + bulge_scale) ** 2
    tabulated = potential.potential_derivative(test_radii)

    relative_error = np.max(np.abs(tabulated - analytic) / analytic)
    tolerance = 0.02
    assert relative_error < tolerance, (
        f"max relative error {relative_error:.3e} exceeds {tolerance}"
    )
    print(
        f"  Hernquist dPhi/dr vs analytic: max rel.err = "
        f"{relative_error:.3e} (tol {tolerance}) -- PASS"
    )


def check_positivity_of_kinematic_quantities() -> None:
    """v_circ^2 and kappa^2 must be non-negative everywhere on the grid."""
    params = GalaxyParams()
    potential = GalaxyPotential(params)

    test_radii = np.geomspace(0.05, 50.0, 200)
    v_circ_squared = potential.circular_speed_squared(test_radii)
    kappa_squared = potential.epicyclic_frequency_squared(test_radii)

    assert np.all(v_circ_squared >= 0.0), "v_circ^2 has negative values"
    assert np.all(kappa_squared > 0.0), "kappa^2 has non-positive values"
    print("  positivity of v_circ^2 and kappa^2 across the grid -- PASS")


def check_disk_surface_density_analytic() -> None:
    """Sigma_d(0) = M_d / (2 pi h^2); Sigma_d(R) decays as exp(-R/h)."""
    params = GalaxyParams()
    potential = GalaxyPotential(params)

    expected_central = params.disk_mass / (2.0 * np.pi * params.disk_scale_length ** 2)
    obtained_central = potential.disk_surface_density(0.0)
    assert np.isclose(obtained_central, expected_central), (
        f"Sigma_d(0): {obtained_central:.4f} vs {expected_central:.4f}"
    )

    obtained_at_h = potential.disk_surface_density(params.disk_scale_length)
    expected_at_h = expected_central / np.e
    assert np.isclose(obtained_at_h, expected_at_h), (
        f"Sigma_d(h): {obtained_at_h:.4f} vs {expected_at_h:.4f}"
    )
    print(
        f"  disk surface density: Sigma(0) = {obtained_central:.4f}, "
        f"Sigma(h)/Sigma(0) = 1/e -- PASS"
    )


def check_circular_speed_order_of_magnitude() -> None:
    """v_circ at R = h should be O(1) in code units (Milky-Way scaled to 262 km/s)."""
    params = GalaxyParams()
    potential = GalaxyPotential(params)
    v_circ_at_h = np.sqrt(potential.circular_speed_squared(1.0))
    assert 0.3 < v_circ_at_h < 1.5, (
        f"v_circ(R=h) = {v_circ_at_h:.4f} outside reasonable range"
    )
    print(
        f"  v_circ(R = h) = {v_circ_at_h:.4f} "
        f"(reasonable for Milky-Way-like model) -- PASS"
    )


def main() -> None:
    """Run all GalaxyPotential validation checks."""
    print("Galaxy potential table validation\n")
    check_total_enclosed_mass_asymptote()
    check_bulge_alone_against_analytic()
    check_bulge_potential_derivative_against_analytic()
    check_positivity_of_kinematic_quantities()
    check_disk_surface_density_analytic()
    check_circular_speed_order_of_magnitude()
    print("\nAll potential table checks PASSED.")


if __name__ == "__main__":
    main()
