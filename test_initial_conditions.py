"""
Validation tests for the IC spatial sampler.

For each component (disk, bulge, halo), the sampler's empirical cumulative
mass distribution is compared against the analytic profile from the manual.
Tolerances are set generously enough to pass on Poisson-noise-only samples
of the manual's particle counts.

Run as a script:
    python test_initial_conditions.py

CMPH Project 3 -- Zoutong Shen / Zhaoyang Chu, 2026.
"""

import numpy as np

from initial_conditions import (
    GalaxyParams,
    bulge_positions,
    disk_positions,
    galaxy_initial_conditions,
    halo_positions,
)


def _max_cdf_discrepancy(empirical: np.ndarray, analytic: np.ndarray) -> float:
    """Maximum absolute deviation between two CDF arrays evaluated on the same grid."""
    return float(np.max(np.abs(empirical - analytic)))


def check_total_mass_and_count() -> None:
    """Verify particle counts and total mass per component sum correctly."""
    params = GalaxyParams()
    positions, masses = galaxy_initial_conditions(params=params, seed=0)

    expected_total = params.num_disk + params.num_bulge + params.num_halo
    assert positions.shape == (expected_total, 3), (
        f"expected positions shape ({expected_total}, 3), got {positions.shape}"
    )
    assert masses.shape == (expected_total,), (
        f"expected masses shape ({expected_total},), got {masses.shape}"
    )

    disk_slice = slice(0, params.num_disk)
    bulge_slice = slice(params.num_disk, params.num_disk + params.num_bulge)
    halo_slice = slice(params.num_disk + params.num_bulge, expected_total)

    disk_total = masses[disk_slice].sum()
    bulge_total = masses[bulge_slice].sum()
    halo_total = masses[halo_slice].sum()

    assert np.isclose(disk_total, params.disk_mass), (
        f"disk mass: {disk_total} vs {params.disk_mass}"
    )
    assert np.isclose(bulge_total, params.bulge_mass), (
        f"bulge mass: {bulge_total} vs {params.bulge_mass}"
    )
    assert np.isclose(halo_total, params.halo_mass), (
        f"halo mass: {halo_total} vs {params.halo_mass}"
    )

    print("  total particle count + per-component masses: PASS")


def check_disk_radial_profile() -> None:
    """Empirical disk radial CDF vs analytic 1 - (1 + R/h) exp(-R/h)."""
    rng = np.random.default_rng(0)
    num_particles = 60_000
    scale_length = 1.0
    scale_height = 0.2
    positions = disk_positions(num_particles, scale_length, scale_height, rng)

    radii = np.sqrt(positions[:, 0] ** 2 + positions[:, 1] ** 2)
    sorted_radii = np.sort(radii)
    empirical_cdf = np.arange(1, num_particles + 1) / num_particles
    analytic_cdf = 1.0 - (1.0 + sorted_radii / scale_length) * np.exp(
        -sorted_radii / scale_length
    )

    discrepancy = _max_cdf_discrepancy(empirical_cdf, analytic_cdf)
    tolerance = 0.01
    assert discrepancy < tolerance, (
        f"disk radial CDF discrepancy {discrepancy:.4f} exceeds {tolerance}"
    )
    print(
        f"  disk radial profile: max CDF deviation = {discrepancy:.4f} "
        f"(tol {tolerance}) -- PASS"
    )


def check_disk_vertical_profile() -> None:
    """Empirical disk vertical CDF vs analytic 0.5(1 + tanh(z/z_0))."""
    rng = np.random.default_rng(1)
    num_particles = 60_000
    scale_length = 1.0
    scale_height = 0.2
    positions = disk_positions(num_particles, scale_length, scale_height, rng)

    sorted_heights = np.sort(positions[:, 2])
    empirical_cdf = np.arange(1, num_particles + 1) / num_particles
    analytic_cdf = 0.5 * (1.0 + np.tanh(sorted_heights / scale_height))

    discrepancy = _max_cdf_discrepancy(empirical_cdf, analytic_cdf)
    tolerance = 0.01
    assert discrepancy < tolerance, (
        f"disk vertical CDF discrepancy {discrepancy:.4f} exceeds {tolerance}"
    )
    print(
        f"  disk vertical profile: max CDF deviation = {discrepancy:.4f} "
        f"(tol {tolerance}) -- PASS"
    )


def check_bulge_radial_profile() -> None:
    """Empirical bulge radial CDF vs analytic Hernquist M(r) = M r^2 / (r+a)^2."""
    rng = np.random.default_rng(2)
    num_particles = 20_000
    bulge_scale = 0.1
    positions = bulge_positions(num_particles, bulge_scale, rng)

    radii = np.linalg.norm(positions, axis=1)
    sorted_radii = np.sort(radii)
    empirical_cdf = np.arange(1, num_particles + 1) / num_particles
    analytic_cdf = sorted_radii ** 2 / (sorted_radii + bulge_scale) ** 2

    # Truncation slightly suppresses the upper tail of the analytic CDF.
    truncation_factor = 500.0
    radius_max = truncation_factor * bulge_scale
    u_max = (radius_max / (radius_max + bulge_scale)) ** 2
    analytic_cdf_truncated = analytic_cdf / u_max

    discrepancy = _max_cdf_discrepancy(empirical_cdf, analytic_cdf_truncated)
    tolerance = 0.015
    assert discrepancy < tolerance, (
        f"bulge radial CDF discrepancy {discrepancy:.4f} exceeds {tolerance}"
    )
    print(
        f"  bulge radial profile: max CDF deviation = {discrepancy:.4f} "
        f"(tol {tolerance}) -- PASS"
    )


def check_halo_radial_profile() -> None:
    """Empirical halo radial CDF vs numerically-built reference CDF (eq 2.57)."""
    rng = np.random.default_rng(3)
    num_particles = 120_000
    halo_core = 1.0
    halo_cutoff = 10.0
    positions = halo_positions(num_particles, halo_core, halo_cutoff, rng)

    radii = np.linalg.norm(positions, axis=1)
    sorted_radii = np.sort(radii)
    empirical_cdf = np.arange(1, num_particles + 1) / num_particles

    reference_grid = np.linspace(0.0, 5.0 * halo_cutoff, 20_000)
    density_unnormalised = np.exp(
        -((reference_grid / halo_cutoff) ** 2)
    ) / (reference_grid ** 2 + halo_core ** 2)
    integrand = reference_grid ** 2 * density_unnormalised
    spacing = np.diff(reference_grid)
    cumulative = np.concatenate(
        ([0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * spacing))
    )
    cumulative /= cumulative[-1]
    analytic_cdf = np.interp(sorted_radii, reference_grid, cumulative)

    discrepancy = _max_cdf_discrepancy(empirical_cdf, analytic_cdf)
    tolerance = 0.01
    assert discrepancy < tolerance, (
        f"halo radial CDF discrepancy {discrepancy:.4f} exceeds {tolerance}"
    )
    print(
        f"  halo radial profile: max CDF deviation = {discrepancy:.4f} "
        f"(tol {tolerance}) -- PASS"
    )


def check_isotropy_of_spheroidal_components() -> None:
    """
    Bulge and halo should be statistically spherically symmetric.

    We test isotropy on the unit-vector directions of each particle (which
    factors out the heavy radial tail of the Hernquist profile that would
    dominate any moment-based test). For an isotropic distribution the
    mean of each Cartesian direction component has expected std 1/sqrt(3N).
    """
    rng = np.random.default_rng(4)
    bulge_pos = bulge_positions(50_000, 0.1, rng)
    halo_pos = halo_positions(50_000, 1.0, 10.0, rng)

    for name, positions in (("bulge", bulge_pos), ("halo", halo_pos)):
        radii = np.linalg.norm(positions, axis=1)
        nonzero = radii > 0.0
        unit_directions = positions[nonzero] / radii[nonzero, None]

        mean_direction = unit_directions.mean(axis=0)
        num_used = unit_directions.shape[0]
        expected_std = 1.0 / np.sqrt(3.0 * num_used)
        five_sigma_threshold = 5.0 * expected_std

        max_offset = float(np.max(np.abs(mean_direction)))
        assert max_offset < five_sigma_threshold, (
            f"{name} mean direction {mean_direction} exceeds 5-sigma "
            f"threshold {five_sigma_threshold:.4f}"
        )
        print(
            f"  {name} isotropy: max |<n_axis>| = {max_offset:.4f} "
            f"(5-sigma = {five_sigma_threshold:.4f}) -- PASS"
        )


def main() -> None:
    """Run all IC spatial-sampling validation checks."""
    print("IC spatial sampling validation\n")
    check_total_mass_and_count()
    check_disk_radial_profile()
    check_disk_vertical_profile()
    check_bulge_radial_profile()
    check_halo_radial_profile()
    check_isotropy_of_spheroidal_components()
    print("\nAll IC spatial-sampling checks PASSED.")


if __name__ == "__main__":
    main()
