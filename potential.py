"""
Tabulated spherical-approximation potential for the three-component galaxy.

Follows the Hernquist (1993a) prescription used by the practicum manual:
each component is reduced to its spherical-shell average, the total enclosed
mass M(<r) is built by numerical integration on a logarithmic radial grid,
and all kinematic quantities needed downstream are derived from M(<r) and
the local total density:

    v_circ^2(R) = G M(<R) / R                    (circular speed, eq 2.55)
    kappa^2(R)  = G M(<R) / R^3 + 4 pi G rho(R)   (epicyclic freq, eq 2.49)
    dPhi/dr     = G M(<r) / r^2                   (Newton's spherical thm)

The disk surface density Sigma_d(R) is kept *exact* (not spherically averaged)
because it enters the vertical Jeans relation (eq 2.40) and the Toomre Q
prescription (eq 2.47) directly.

The manual notes (figure 2.1 caption) that this spherical approximation
underestimates the disk's true circular speed by ~ 15 %; this is accepted in
exchange for analytic simplicity and matches Hernquist 1993a.

Units: G = 1.

CMPH Project 3 -- Zoutong Shen / Zhaoyang Chu, 2026.
"""

from math import erf
from typing import Union

import numpy as np

from initial_conditions import GalaxyParams


ArrayOrScalar = Union[float, np.ndarray]


class GalaxyPotential:
    """
    Spherical-approximation total potential of the three-component galaxy.

    All radial quantities are tabulated on a logarithmic grid at construction
    and interpolated linearly on access. Attributes are read-only by convention
    (leading underscore); use the public accessor methods.
    """

    def __init__(
        self,
        params: GalaxyParams,
        radial_grid_size: int = 2_000,
        inner_radius: float = 1.0e-3,
        outer_radius_factor: float = 10.0,
        polar_grid_size: int = 200,
    ) -> None:
        """
        Build the tabulated potential.

        Args:
            params:               GalaxyParams instance.
            radial_grid_size:     number of points on the log-radial grid.
            inner_radius:         smallest radius in the grid (avoid r=0).
            outer_radius_factor:  outermost radius = this * halo_cutoff.
            polar_grid_size:      number of theta points for the spherical
                                  average of the disk density.
        """
        self.params = params
        self._radii = np.geomspace(
            inner_radius, outer_radius_factor * params.halo_cutoff, radial_grid_size
        )
        self._build_density_tables(polar_grid_size)
        self._build_enclosed_mass_table()
        self._build_isotropic_jeans_tables()

    def _build_density_tables(self, polar_grid_size: int) -> None:
        """Populate spherical-shell-averaged density tables for each component."""
        polar_angles = np.linspace(0.0, np.pi, polar_grid_size)
        sin_polar = np.sin(polar_angles)
        cos_polar = np.cos(polar_angles)

        radii_column = self._radii[:, None]
        cylindrical_R = radii_column * sin_polar[None, :]
        height_z = radii_column * cos_polar[None, :]

        h = self.params.disk_scale_length
        z0 = self.params.disk_scale_height
        disk_prefactor = self.params.disk_mass / (
            4.0 * np.pi * h * h * z0
        )
        # sech^2(z/z0) -> 0 exponentially for |z| >> z0; clipping the cosh
        # argument at |x| = 350 avoids overflow of cosh(x)^2 while leaving the
        # density unchanged to many orders of magnitude.
        clipped_height_argument = np.clip(height_z / z0, -350.0, 350.0)
        disk_density_2d = (
            disk_prefactor
            * np.exp(-cylindrical_R / h)
            / np.cosh(clipped_height_argument) ** 2
        )
        self._rho_disk = 0.5 * np.trapezoid(
            disk_density_2d * sin_polar[None, :], polar_angles, axis=1
        )

        a = self.params.bulge_scale
        self._rho_bulge = (
            self.params.bulge_mass * a
            / (2.0 * np.pi * self._radii * (self._radii + a) ** 3)
        )

        r_c = self.params.halo_cutoff
        gamma = self.params.halo_core
        q = gamma / r_c
        alpha = 1.0 / (1.0 - np.sqrt(np.pi) * q * np.exp(q * q) * (1.0 - erf(q)))
        self._rho_halo = (
            self.params.halo_mass
            * alpha
            * np.exp(-((self._radii / r_c) ** 2))
            / (2.0 * np.pi ** 1.5 * r_c * (self._radii ** 2 + gamma ** 2))
        )

        self._rho_total = self._rho_disk + self._rho_bulge + self._rho_halo

    def _build_enclosed_mass_table(self) -> None:
        """Cumulative enclosed mass M(<r) via trapezoidal integration."""
        integrand = 4.0 * np.pi * self._radii ** 2 * self._rho_total
        spacing = np.diff(self._radii)
        increments = 0.5 * (integrand[1:] + integrand[:-1]) * spacing
        self._mass_enclosed = np.concatenate(([0.0], np.cumsum(increments)))

        self._v_circ_squared = self._mass_enclosed / self._radii
        self._epicyclic_squared = (
            self._mass_enclosed / self._radii ** 3
            + 4.0 * np.pi * self._rho_total
        )

    def _build_isotropic_jeans_tables(self) -> None:
        """
        Isotropic spherical Jeans equation (eq 2.63):

            sigma_r^2(r) = (1 / rho(r)) * integral_r^infty rho(r') dPhi/dr' dr'

        applied to the bulge and halo separately, using the *total* potential's
        radial derivative dPhi/dr = G M_total(<r) / r^2. The integral is built
        once for each component on the master radial grid.
        """
        potential_derivative = self._mass_enclosed / self._radii ** 2

        for component_name, component_density in (
            ("bulge", self._rho_bulge),
            ("halo", self._rho_halo),
        ):
            integrand = component_density * potential_derivative
            spacing = np.diff(self._radii)
            forward_increments = 0.5 * (integrand[1:] + integrand[:-1]) * spacing
            cumulative_forward = np.concatenate(
                ([0.0], np.cumsum(forward_increments))
            )
            cumulative_from_r = cumulative_forward[-1] - cumulative_forward

            safe_density = np.maximum(component_density, 1.0e-300)
            sigma_r_squared = cumulative_from_r / safe_density
            sigma_r_squared = np.maximum(sigma_r_squared, 0.0)

            setattr(self, f"_sigma_r_squared_{component_name}", sigma_r_squared)

    def bulge_velocity_dispersion_squared(self, r: ArrayOrScalar) -> ArrayOrScalar:
        """
        Isotropic radial velocity dispersion squared sigma_r^2(r) for the bulge,
        from spherical Jeans (eq 2.63) in the total potential. By isotropy,
        each Cartesian component shares this variance.
        """
        return np.interp(r, self._radii, self._sigma_r_squared_bulge)

    def halo_velocity_dispersion_squared(self, r: ArrayOrScalar) -> ArrayOrScalar:
        """
        Isotropic radial velocity dispersion squared sigma_r^2(r) for the halo,
        from spherical Jeans (eq 2.63) in the total potential.
        """
        return np.interp(r, self._radii, self._sigma_r_squared_halo)

    def enclosed_mass(self, r: ArrayOrScalar) -> ArrayOrScalar:
        """Total mass enclosed within spherical radius r."""
        return np.interp(r, self._radii, self._mass_enclosed)

    def total_density(self, r: ArrayOrScalar) -> ArrayOrScalar:
        """Spherical-shell-averaged total density at radius r."""
        return np.interp(r, self._radii, self._rho_total)

    def circular_speed_squared(self, cylindrical_radius: ArrayOrScalar) -> ArrayOrScalar:
        """
        Squared circular speed in the disk plane (spherical approximation):
        v_c^2(R) = G M(<R) / R.
        """
        return np.interp(
            cylindrical_radius, self._radii, self._v_circ_squared
        )

    def epicyclic_frequency_squared(
        self, cylindrical_radius: ArrayOrScalar
    ) -> ArrayOrScalar:
        """
        Squared epicyclic frequency in the disk plane (eq 2.49, spherical):
        kappa^2(R) = G M(<R) / R^3 + 4 pi G rho(R).
        """
        return np.interp(
            cylindrical_radius, self._radii, self._epicyclic_squared
        )

    def potential_derivative(self, r: ArrayOrScalar) -> ArrayOrScalar:
        """
        Radial derivative of the gravitational potential:
        dPhi/dr = G M(<r) / r^2 (spherical Newton's theorem).
        Used for the bulge/halo Jeans equation (eq 2.63).
        """
        mass = self.enclosed_mass(r)
        return mass / (np.asarray(r) ** 2)

    def disk_surface_density(
        self, cylindrical_radius: ArrayOrScalar
    ) -> ArrayOrScalar:
        """
        Disk surface density Sigma_d(R) = (M_d / 2 pi h^2) exp(-R/h).

        Analytic — independent of the spherical approximation. Used in the
        vertical Jeans relation (eq 2.40) and the Toomre stability criterion
        (eq 2.47).
        """
        h = self.params.disk_scale_length
        normalisation = self.params.disk_mass / (2.0 * np.pi * h * h)
        return normalisation * np.exp(-np.asarray(cylindrical_radius) / h)


def main() -> None:
    """Build the default-parameter potential and print a sanity summary."""
    params = GalaxyParams()
    potential = GalaxyPotential(params)

    expected_total_mass = params.disk_mass + params.bulge_mass + params.halo_mass
    outermost = potential._radii[-1]
    print("Galaxy potential table built.")
    print(
        f"  radial grid: {potential._radii[0]:.3e} ... "
        f"{outermost:.3e} (n = {potential._radii.size})"
    )
    print(
        f"  total enclosed mass at outer edge = "
        f"{potential.enclosed_mass(outermost):.4f}"
    )
    print(f"  expected (M_d + M_b + M_h) = {expected_total_mass:.4f}")
    print(
        f"  v_circ at R=1: {np.sqrt(potential.circular_speed_squared(1.0)):.4f}"
    )
    print(
        f"  Sigma_d at R=0: {potential.disk_surface_density(0.0):.4f}"
    )


if __name__ == "__main__":
    main()
