"""
Analytic answers to the section 3.0.6 questions (Q8, Q9, Q10).

These questions are pen-and-paper estimates about the *initial* galaxy
model and the merger remnant; this module evaluates the closed-form
expressions numerically so the printed numbers can go straight into the
answer write-up and the presentation.

Definitions used (all consistent with the manual):

  * Exponential disk    Sigma(R) = Sigma_0 exp(-R / r_d)
        total mass            M_d   = 2 pi Sigma_0 r_d^2
        enclosed (in plane)   M(<R) = M_d [1 - (1 + R/r_d) exp(-R/r_d)]
        half-mass radius      (1 + x) exp(-x) = 1/2,  x = R_half / r_d

  * Hernquist (1990) bulge   M(<r) = M_b r^2 / (r + a)^2
        half-mass radius      r_h = (1 + sqrt(2)) a

  * "Dynamical time" of a (mostly spherical) body of mean density rho:
    a test particle released from rest at radius r in a homogeneous
    sphere obeys  r'' = -(4 pi G rho / 3) r, i.e. simple harmonic motion
    with omega^2 = 4 pi G rho / 3. It reaches the centre in a quarter
    period, giving

        t_dyn = (pi / 2) / omega = sqrt( 3 pi / (16 G rho) )
              = (pi / 2) sqrt( r^3 / (G M(<r)) )                  (Q9)

  * For the disk the manual instead defines the dynamical time as the
    rotation period at the disk half-mass radius,
        t_dyn = 2 pi R_half / v_circ(R_half),  v_circ^2 = G M(<R)/R,
    which equals  2 pi sqrt( R^3 / (G M(<R)) )                    (Q8).

The simulation uses G = 1, fixed leapfrog step dt = 0.125 and Plummer
softening eps = 0.1 for the merger.

Run as a script:
    python analytic_questions.py
    python analytic_questions.py --merger data/merger_N80k_eps01.npz

CMPH Project 3 -- Zoutong Shen / Zhaoyang Chu, 2026.
"""

import argparse

import numpy as np

from initial_conditions import GalaxyParams
from merger_analysis import shrinking_sphere_centre, stellar_mask
from physical_units import milky_way_unit_system

GRAVITATIONAL_CONSTANT_CODE = 1.0
LEAPFROG_TIMESTEP_CODE = 0.125
MERGER_SOFTENING_CODE = 0.1


def exponential_disk_half_mass_radius(scale_length: float) -> float:
    """
    Half-mass radius of an exponential disk in the same length units as
    ``scale_length``.

    Solves ``(1 + x) exp(-x) = 1/2`` for ``x = R_half / r_d`` by bisection
    (the enclosed-mass fraction is monotone, so bisection is exact to
    machine tolerance and needs no external solver).
    """
    def enclosed_fraction(scaled_radius: float) -> float:
        """Mass fraction of an exponential disk inside ``x = R / r_d``."""
        return 1.0 - (1.0 + scaled_radius) * np.exp(-scaled_radius)

    low, high = 0.0, 20.0
    for _ in range(200):
        midpoint = 0.5 * (low + high)
        if enclosed_fraction(midpoint) < 0.5:
            low = midpoint
        else:
            high = midpoint
    return 0.5 * (low + high) * scale_length


def exponential_disk_enclosed_mass(
    radius: float, disk_mass: float, scale_length: float
) -> float:
    """In-plane mass of an exponential disk inside cylindrical ``radius``."""
    scaled_radius = radius / scale_length
    return disk_mass * (
        1.0 - (1.0 + scaled_radius) * np.exp(-scaled_radius)
    )


def hernquist_enclosed_mass(
    radius: float, total_mass: float, scale: float
) -> float:
    """Hernquist (1990) cumulative mass ``M r^2 / (r + a)^2``."""
    return total_mass * radius ** 2 / (radius + scale) ** 2


def hernquist_half_mass_radius(scale: float) -> float:
    """
    Half-mass radius of a Hernquist sphere.

    Setting ``r^2 / (r + a)^2 = 1/2`` gives ``r = (1 + sqrt(2)) a``.
    """
    return (1.0 + np.sqrt(2.0)) * scale


def truncated_isothermal_enclosed_mass(
    radius: float,
    halo_mass: float,
    halo_core: float,
    halo_cutoff: float,
    grid_extent_factor: float = 5.0,
    grid_size: int = 20_000,
) -> float:
    """
    Mass of the truncated-isothermal halo inside ``radius``.

    The density ``rho ~ exp(-r^2/r_c^2) / (r^2 + gamma^2)`` has no
    closed-form cumulative, so this integrates ``4 pi r^2 rho`` on the
    same fine grid the IC sampler uses and normalises the total to
    ``halo_mass`` over the truncation extent, reproducing the model the
    particles were actually drawn from.
    """
    outer_radius = grid_extent_factor * halo_cutoff
    radii_grid = np.linspace(0.0, outer_radius, grid_size)
    density = np.exp(
        -((radii_grid / halo_cutoff) ** 2)
    ) / (radii_grid ** 2 + halo_core ** 2)
    integrand = radii_grid ** 2 * density
    spacing = np.diff(radii_grid)
    increments = 0.5 * (integrand[1:] + integrand[:-1]) * spacing
    cumulative = np.concatenate(([0.0], np.cumsum(increments)))
    enclosed_fraction = np.interp(radius, radii_grid, cumulative) / cumulative[-1]
    return halo_mass * enclosed_fraction


def mean_density_dynamical_time(
    enclosed_mass: float, radius: float
) -> float:
    """
    Free-fall / dynamical time ``(pi/2) sqrt(r^3 / (G M(<r)))`` (Q9).

    This is the time for a particle released from rest at ``radius`` to
    reach the centre of a uniform sphere of the same mean density.
    """
    return 0.5 * np.pi * np.sqrt(
        radius ** 3 / (GRAVITATIONAL_CONSTANT_CODE * enclosed_mass)
    )


def circular_orbit_period(enclosed_mass: float, radius: float) -> float:
    """Period of a circular orbit, ``2 pi sqrt(r^3 / (G M(<r)))``."""
    return 2.0 * np.pi * np.sqrt(
        radius ** 3 / (GRAVITATIONAL_CONSTANT_CODE * enclosed_mass)
    )


def answer_q8(params: GalaxyParams) -> None:
    """Q8: disk total mass, half-mass radius and dynamical time."""
    units = milky_way_unit_system()
    scale_length = params.disk_scale_length
    disk_mass = params.disk_mass

    central_surface_density = disk_mass / (2.0 * np.pi * scale_length ** 2)
    half_mass_radius = exponential_disk_half_mass_radius(scale_length)

    # The disk material orbits in the *total* potential, so the rotation
    # period uses the mass of all three components inside R_half (the disk
    # mass there is exactly M_d/2 by definition of the half-mass radius).
    disk_inside = 0.5 * disk_mass
    bulge_inside = hernquist_enclosed_mass(
        half_mass_radius, params.bulge_mass, params.bulge_scale
    )
    halo_inside = truncated_isothermal_enclosed_mass(
        half_mass_radius, params.halo_mass,
        params.halo_core, params.halo_cutoff,
    )
    total_inside = disk_inside + bulge_inside + halo_inside

    dynamical_time = circular_orbit_period(total_inside, half_mass_radius)
    timestep_fraction = LEAPFROG_TIMESTEP_CODE / dynamical_time

    print("Q8 -- initial exponential disk")
    print(f"  total disk mass      M_d   = {disk_mass:.4g} code "
          f"= {disk_mass * units.mass_msun:.4g} M_sun")
    print(f"  central density    Sigma_0 = {central_surface_density:.5g} code "
          f"= {central_surface_density * units.surface_density_msun_pc2:.5g}"
          f" M_sun/pc^2")
    print(f"  half-mass radius   R_half  = {half_mass_radius:.4g} code "
          f"= {half_mass_radius * units.length_kpc:.4g} kpc "
          f"(= {half_mass_radius / scale_length:.4f} r_d)")
    print(f"  enclosed within R_half: disk {disk_inside:.4g}, "
          f"bulge {bulge_inside:.4g}, halo {halo_inside:.4g}, "
          f"total {total_inside:.4g}")
    print(f"  dynamical time (rotation period at R_half) "
          f"t_dyn = {dynamical_time:.4g} code "
          f"= {dynamical_time * units.time_year / 1.0e6:.4g} Myr")
    print(f"  leapfrog step dt = {LEAPFROG_TIMESTEP_CODE} code is "
          f"{timestep_fraction:.3e} of t_dyn "
          f"(~{1.0 / timestep_fraction:.0f} steps per orbit)\n")


def answer_q9(merger_path: str, params: GalaxyParams) -> None:
    """Q9: dynamical time of the merger remnant at fixed physical radii."""
    units = milky_way_unit_system()
    print("Q9 -- equation of motion and remnant dynamical time")
    print("  Homogeneous sphere: r'' = -(4 pi G rho / 3) r  (SHM),")
    print("  omega^2 = 4 pi G rho / 3; particle reaches r = 0 in a quarter")
    print("  period -> t_dyn = sqrt(3 pi / (16 G rho)) "
          "= (pi/2) sqrt(r^3 / (G M(<r))).")

    timestep_years = LEAPFROG_TIMESTEP_CODE * units.time_year
    radii_kpc = np.array([0.5, 1.0, 3.0, 5.0])

    try:
        trajectory = np.load(merger_path)
    except FileNotFoundError:
        print(f"\n  (merger file {merger_path!r} not found -- run "
              f"evolve_merger.py first to fill in the remnant table)\n")
        return

    positions = trajectory["positions_history"][-1].astype(np.float64)
    masses = trajectory["masses"].astype(np.float64)
    stars = stellar_mask(
        int(trajectory["num_disk"]),
        int(trajectory["num_bulge"]),
        int(trajectory["num_halo"]),
    )
    # Centre on the densest concentration so the enclosed mass is measured
    # about the remnant, robust to the diffuse tidal-debris envelope.
    centre = shrinking_sphere_centre(positions[stars], masses[stars])
    radii_code = np.linalg.norm(positions - centre, axis=1)
    radii_kpc_to_code = radii_kpc / units.length_kpc

    print(f"\n  total merger remnant ({len(masses)} particles, all "
          f"components), dt = {timestep_years / 1.0e6:.4g} Myr:")
    print("    r [kpc]   M(<r) [M_sun]   t_dyn [Myr]   dt / t_dyn")
    for radius_kpc, radius_code in zip(radii_kpc, radii_kpc_to_code):
        enclosed_mass = masses[radii_code < radius_code].sum()
        dynamical_time = mean_density_dynamical_time(
            enclosed_mass, radius_code
        )
        dynamical_time_years = dynamical_time * units.time_year
        print(f"    {radius_kpc:6.1f}   {enclosed_mass * units.mass_msun:12.4g}"
              f"   {dynamical_time_years / 1.0e6:10.4g}"
              f"   {timestep_years / dynamical_time_years:10.3e}")
    print()


def answer_q10(params: GalaxyParams) -> None:
    """Q10: bulge half-mass radius, its dynamical time, resolvability."""
    units = milky_way_unit_system()
    half_mass_radius = hernquist_half_mass_radius(params.bulge_scale)
    enclosed_mass = 0.5 * params.bulge_mass
    dynamical_time = mean_density_dynamical_time(
        enclosed_mass, half_mass_radius
    )
    timestep_fraction = LEAPFROG_TIMESTEP_CODE / dynamical_time

    print("Q10 -- bulge component")
    print(f"  Hernquist scale       a   = {params.bulge_scale:.4g} code "
          f"= {params.bulge_scale * units.length_kpc:.4g} kpc")
    print(f"  bulge half-mass radius r_h = {half_mass_radius:.4g} code "
          f"= {half_mass_radius * units.length_kpc:.4g} kpc")
    print(f"  dynamical time at r_h t_dyn = {dynamical_time:.4g} code "
          f"= {dynamical_time * units.time_year / 1.0e6:.4g} Myr")
    print(f"  leapfrog step dt = {LEAPFROG_TIMESTEP_CODE} code is "
          f"{timestep_fraction:.3e} of t_dyn "
          f"(~{1.0 / timestep_fraction:.0f} steps per orbit)")
    print(f"  Plummer softening eps = {MERGER_SOFTENING_CODE} code "
          f"= {MERGER_SOFTENING_CODE * units.length_kpc:.4g} kpc "
          f"= {MERGER_SOFTENING_CODE / params.bulge_scale:.2f} a")
    print("  => eps ~ a: the bulge's own scale radius is at the softening")
    print("     length, so its central cusp (r < eps) is smoothed away and")
    print("     not resolved; only the half-mass region and outwards, where")
    print("     there are still many steps per orbit, is followed faithfully.\n")


def main() -> None:
    """Print the worked answers to Q8, Q9 and Q10."""
    parser = argparse.ArgumentParser(
        description="Analytic answers to manual questions Q8-Q10."
    )
    parser.add_argument(
        "--merger", type=str, default="data/merger_N80k_eps01.npz",
        help="merger trajectory .npz used for the Q9 remnant table",
    )
    arguments = parser.parse_args()

    params = GalaxyParams()
    answer_q8(params)
    answer_q9(arguments.merger, params)
    answer_q10(params)


if __name__ == "__main__":
    main()
