"""
Conversion from internal (G = 1) code units to physical units.

This module answers manual question Q3 (derive the physical time unit in
years and the velocity unit in km/s) and additionally supplies the mass,
length and surface-density conversions needed for the Q4 physical-unit
plots.

The simulation runs with the gravitational constant set to G = 1 and two
chosen scales fixed by the manual (page 23):

    length unit  L = 1   ->   3.5  kpc
    mass unit    M = 1   ->   5.6 x 10^10  M_sun

With G dimensionless and equal to one, the time unit follows from
requiring Newton's constant to take its physical value in SI units:

    G_code = G_phys * [M] * [T]^2 / [L]^3 = 1
    =>  [T] = sqrt( [L]^3 / (G_phys * [M]) )
    =>  [V] = [L] / [T] = sqrt( G_phys * [M] / [L] )

Physical constants used (exactly the values quoted in the manual):

    G_phys  = 6.672 x 10^-11  m^3 kg^-1 s^-2
    1 kpc   = 3.0856 x 10^19  m
    1 M_sun = 1.989  x 10^30  kg

CMPH Project 3 -- Zoutong Shen / Zhaoyang Chu, 2026.
"""

from dataclasses import dataclass

# Physical constants (manual page 23, SI).
GRAVITATIONAL_CONSTANT_SI = 6.672e-11        # m^3 kg^-1 s^-2
METRES_PER_KILOPARSEC = 3.0856e19            # m
KILOGRAMS_PER_SOLAR_MASS = 1.989e30          # kg
SECONDS_PER_JULIAN_YEAR = 3.1557e7           # s (365.25 d)
METRES_PER_PARSEC = METRES_PER_KILOPARSEC / 1.0e3

# Manual-prescribed scaling of the code units to the Milky Way (page 23).
LENGTH_UNIT_KILOPARSEC = 3.5
MASS_UNIT_SOLAR_MASSES = 5.6e10


@dataclass(frozen=True)
class UnitSystem:
    """
    Multiplicative factors converting a code-unit quantity to physical units.

    A code-unit value is converted by multiplying with the matching factor,
    e.g. ``time_physical_years = time_code * units.time_year``.
    """

    length_kpc: float          # code length -> kiloparsec
    mass_msun: float           # code mass   -> solar mass
    time_year: float           # code time   -> year
    velocity_kms: float        # code speed  -> km/s
    surface_density_msun_pc2: float  # code Sigma -> M_sun / pc^2


def milky_way_unit_system() -> UnitSystem:
    """
    Build the unit system for the manual's Milky Way scaling.

    Returns:
        A :class:`UnitSystem` whose factors implement the derivation in this
        module's docstring.
    """
    length_metres = LENGTH_UNIT_KILOPARSEC * METRES_PER_KILOPARSEC
    mass_kilograms = MASS_UNIT_SOLAR_MASSES * KILOGRAMS_PER_SOLAR_MASS

    time_seconds = (
        length_metres ** 3
        / (GRAVITATIONAL_CONSTANT_SI * mass_kilograms)
    ) ** 0.5
    velocity_metres_per_second = length_metres / time_seconds

    length_parsecs = LENGTH_UNIT_KILOPARSEC * 1.0e3
    surface_density_msun_pc2 = MASS_UNIT_SOLAR_MASSES / length_parsecs ** 2

    return UnitSystem(
        length_kpc=LENGTH_UNIT_KILOPARSEC,
        mass_msun=MASS_UNIT_SOLAR_MASSES,
        time_year=time_seconds / SECONDS_PER_JULIAN_YEAR,
        velocity_kms=velocity_metres_per_second / 1.0e3,
        surface_density_msun_pc2=surface_density_msun_pc2,
    )


def main() -> None:
    """Print the derived units -- this is the worked answer to Q3."""
    units = milky_way_unit_system()
    print("Q3 -- physical units for the manual's Milky Way scaling\n")
    print(f"  length unit    1 -> {units.length_kpc:.4g} kpc")
    print(f"  mass unit      1 -> {units.mass_msun:.4g} M_sun")
    print(f"  time unit      1 -> {units.time_year:.5g} yr"
          f"  (= {units.time_year / 1.0e6:.4g} Myr)")
    print(f"  velocity unit  1 -> {units.velocity_kms:.5g} km/s")
    print(f"  surface dens.  1 -> {units.surface_density_msun_pc2:.5g}"
          f" M_sun/pc^2")
    print(
        "\n  tstop = 100 code units"
        f" -> {100.0 * units.time_year / 1.0e6:.3g} Myr"
    )
    print(
        f"  dt = 0.125 code units"
        f" -> {0.125 * units.time_year / 1.0e6:.4g} Myr"
    )


if __name__ == "__main__":
    main()
