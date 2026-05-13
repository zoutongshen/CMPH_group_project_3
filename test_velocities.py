"""
Validation tests for the velocity sampler.

Checks expected statistical and physical properties of each component's
sampled velocities. Tolerances are set to pass on N=20000 particle counts
with normal random-number variance.

Run as a script:
    python test_velocities.py

CMPH Project 3 -- Zoutong Shen / Zhaoyang Chu, 2026.
"""

import numpy as np

from initial_conditions import GalaxyParams, galaxy_initial_conditions
from potential import GalaxyPotential
from velocities import (
    bulge_velocities,
    disk_velocities,
    galaxy_velocities,
    halo_velocities,
)


def _build_test_galaxy(seed: int = 0):
    """Return (positions, masses, potential) for a default-params galaxy."""
    params = GalaxyParams()
    positions, masses = galaxy_initial_conditions(params=params, seed=seed)
    potential = GalaxyPotential(params)
    return positions, masses, potential


def check_bulge_isotropy_in_velocity_space() -> None:
    """Bulge velocities should have <v> ~ 0 and equal per-axis variance."""
    positions, masses, potential = _build_test_galaxy(seed=0)
    params = potential.params
    bulge_positions = positions[
        params.num_disk : params.num_disk + params.num_bulge
    ]
    rng = np.random.default_rng(42)
    velocities = bulge_velocities(bulge_positions, potential, rng)

    mean_velocity = velocities.mean(axis=0)
    variances = velocities.var(axis=0)

    typical_speed = np.sqrt(variances.mean())
    mean_fraction = float(np.max(np.abs(mean_velocity)) / typical_speed)
    variance_ratio_spread = float(variances.max() / variances.min() - 1.0)

    assert mean_fraction < 0.1, (
        f"bulge mean velocity fraction {mean_fraction:.3f} exceeds 0.1"
    )
    assert variance_ratio_spread < 0.1, (
        f"bulge per-axis variance spread {variance_ratio_spread:.3f} exceeds 0.1"
    )
    print(
        f"  bulge isotropy: |<v>|/sigma = {mean_fraction:.4f}, "
        f"variance ratio spread = {variance_ratio_spread:.4f} -- PASS"
    )


def check_halo_isotropy_in_velocity_space() -> None:
    """Halo velocities should have <v> ~ 0 and equal per-axis variance."""
    positions, masses, potential = _build_test_galaxy(seed=0)
    params = potential.params
    halo_positions = positions[params.num_disk + params.num_bulge :]
    rng = np.random.default_rng(43)
    velocities = halo_velocities(halo_positions, potential, rng)

    mean_velocity = velocities.mean(axis=0)
    variances = velocities.var(axis=0)

    typical_speed = np.sqrt(variances.mean())
    mean_fraction = float(np.max(np.abs(mean_velocity)) / typical_speed)
    variance_ratio_spread = float(variances.max() / variances.min() - 1.0)

    assert mean_fraction < 0.05, (
        f"halo mean velocity fraction {mean_fraction:.3f} exceeds 0.05"
    )
    assert variance_ratio_spread < 0.05, (
        f"halo per-axis variance spread {variance_ratio_spread:.3f} exceeds 0.05"
    )
    print(
        f"  halo isotropy: |<v>|/sigma = {mean_fraction:.4f}, "
        f"variance ratio spread = {variance_ratio_spread:.4f} -- PASS"
    )


def check_disk_is_rotating() -> None:
    """
    Disk should have substantial mean tangential velocity in +phi direction
    (or -phi, depending on Gaussian random sign of mean v_phi)
    and far smaller mean radial / vertical motion.
    """
    positions, masses, potential = _build_test_galaxy(seed=0)
    params = potential.params
    disk_positions = positions[: params.num_disk]
    rng = np.random.default_rng(44)
    velocities = disk_velocities(disk_positions, potential, rng)

    cylindrical_R = np.sqrt(
        disk_positions[:, 0] ** 2 + disk_positions[:, 1] ** 2
    )
    azimuthal_angle = np.arctan2(disk_positions[:, 1], disk_positions[:, 0])
    cos_phi = np.cos(azimuthal_angle)
    sin_phi = np.sin(azimuthal_angle)
    velocity_radial = velocities[:, 0] * cos_phi + velocities[:, 1] * sin_phi
    velocity_azimuthal = -velocities[:, 0] * sin_phi + velocities[:, 1] * cos_phi

    mean_radial = float(velocity_radial.mean())
    mean_azimuthal = float(velocity_azimuthal.mean())
    rms_radial = float(np.sqrt((velocity_radial ** 2).mean()))
    rms_azimuthal = float(np.sqrt((velocity_azimuthal ** 2).mean()))

    assert abs(mean_radial) < 0.05 * rms_radial, (
        f"disk mean radial drift {mean_radial:.3e} too large"
    )
    assert abs(mean_azimuthal) > 0.5 * rms_azimuthal, (
        f"disk mean rotation {mean_azimuthal:.3e} too small compared to RMS"
    )
    print(
        f"  disk rotation: <v_R> = {mean_radial:.4f} (rms {rms_radial:.4f}), "
        f"<v_phi> = {mean_azimuthal:.4f} (rms {rms_azimuthal:.4f}) -- PASS"
    )


def check_disk_dispersion_anisotropy() -> None:
    """Disk should be vertically colder than radially (sigma_z < sigma_R)."""
    positions, masses, potential = _build_test_galaxy(seed=0)
    params = potential.params
    disk_positions = positions[: params.num_disk]
    rng = np.random.default_rng(45)
    velocities = disk_velocities(disk_positions, potential, rng)

    cylindrical_R = np.sqrt(
        disk_positions[:, 0] ** 2 + disk_positions[:, 1] ** 2
    )
    azimuthal_angle = np.arctan2(disk_positions[:, 1], disk_positions[:, 0])
    cos_phi = np.cos(azimuthal_angle)
    sin_phi = np.sin(azimuthal_angle)
    velocity_radial = velocities[:, 0] * cos_phi + velocities[:, 1] * sin_phi
    velocity_vertical = velocities[:, 2]

    sigma_radial = float(np.std(velocity_radial))
    sigma_vertical = float(np.std(velocity_vertical))

    assert sigma_vertical < sigma_radial, (
        f"disk should be vertically colder: sigma_z {sigma_vertical:.4f} "
        f"vs sigma_R {sigma_radial:.4f}"
    )
    print(
        f"  disk dispersion anisotropy: sigma_z = {sigma_vertical:.4f} "
        f"< sigma_R = {sigma_radial:.4f} -- PASS"
    )


def check_no_nans_in_complete_ic() -> None:
    """Sanity: no NaN or inf values anywhere in the full IC velocities."""
    positions, masses, potential = _build_test_galaxy(seed=0)
    rng = np.random.default_rng(46)
    velocities = galaxy_velocities(positions, potential, rng)

    assert np.all(np.isfinite(velocities)), "non-finite velocity values present"
    print(f"  no NaN/inf in full IC velocities -- PASS")


def main() -> None:
    """Run all velocity-sampling validation checks."""
    print("Velocity sampling validation\n")
    check_bulge_isotropy_in_velocity_space()
    check_halo_isotropy_in_velocity_space()
    check_disk_is_rotating()
    check_disk_dispersion_anisotropy()
    check_no_nans_in_complete_ic()
    print("\nAll velocity-sampling checks PASSED.")


if __name__ == "__main__":
    main()
