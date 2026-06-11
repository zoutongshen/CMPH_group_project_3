"""
Pedagogical sketch of the parabolic two-body encounter orbit used to set
up the merger initial condition (``merger_initial_conditions.py``).

Numerically integrates the Kepler problem for equal-mass point particles
with the exact same initial separation, pericentre and parabolic-energy
setup as the simulation, and plots both galaxies' trajectories in the
centre-of-mass frame. A tiny Plummer softening is added so the head-on
(pericentre = 0) case stays numerically well-behaved at impact.

Not used by the production pipeline; just an illustration.

Run as a script:
    python orbit_illustration.py                # produces all four pericentres
    python orbit_illustration.py --pericentre 5 # one panel only
"""

import argparse
import os
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def integrate_kepler(
    separation: float,
    pericentre: float,
    total_mass: float,
    duration: float,
    timestep: float,
    softening: float = 0.3,
) -> np.ndarray:
    """Leapfrog-integrate two equal-mass point particles in COM frame."""
    # Velocities from two-body parabolic formulae (mirrors
    # ``two_body_encounter_velocity`` in merger_initial_conditions.py).
    parabolic_speed = np.sqrt(2.0 * total_mass / separation)
    specific_angular_momentum = np.sqrt(2.0 * total_mass * pericentre)
    tangential_velocity = specific_angular_momentum / separation
    radial_velocity = -np.sqrt(
        max(parabolic_speed ** 2 - tangential_velocity ** 2, 0.0)
    )

    # Two equal-mass bodies, COM at origin, motion in x-y plane.
    positions = np.array([
        [+0.5 * separation, 0.0],
        [-0.5 * separation, 0.0],
    ])
    rel_velocity = np.array([radial_velocity, tangential_velocity])
    velocities = np.array([+0.5 * rel_velocity, -0.5 * rel_velocity])

    num_steps = int(round(duration / timestep))
    history = np.zeros((num_steps + 1, 2, 2))
    history[0] = positions

    mass_half = 0.5 * total_mass  # each body holds half the total mass
    softening_sq = softening ** 2
    for step in range(num_steps):
        # Pairwise softened gravity (only one pair).
        delta = positions[1] - positions[0]
        distance_sq = float(delta @ delta) + softening_sq
        distance_cubed = distance_sq ** 1.5
        accel_on_0 = +mass_half * delta / distance_cubed
        accel_on_1 = -mass_half * delta / distance_cubed
        # Velocity-Verlet.
        velocities[0] += 0.5 * timestep * accel_on_0
        velocities[1] += 0.5 * timestep * accel_on_1
        positions += timestep * velocities
        delta = positions[1] - positions[0]
        distance_sq = float(delta @ delta) + softening_sq
        distance_cubed = distance_sq ** 1.5
        accel_on_0 = +mass_half * delta / distance_cubed
        accel_on_1 = -mass_half * delta / distance_cubed
        velocities[0] += 0.5 * timestep * accel_on_0
        velocities[1] += 0.5 * timestep * accel_on_1
        history[step + 1] = positions
    return history


def make_one_illustration(
    pericentre: float,
    output_path: str,
    *,
    separation: float = 30.0,
    total_mass: float = 2.0,
    duration: float = 180.0,
    timestep: float = 0.05,
    view_kpc: float = 22.0,
) -> None:
    """Render a single orbit illustration at the requested pericentre."""
    history = integrate_kepler(
        separation=separation,
        pericentre=pericentre,
        total_mass=total_mass,
        duration=duration,
        timestep=timestep,
    )
    pos_1 = history[:, 0, :]
    pos_2 = history[:, 1, :]
    separations = np.linalg.norm(pos_1 - pos_2, axis=1)
    pericentre_index = int(np.argmin(separations))
    actual_pericentre = float(separations[pericentre_index])
    print(f"  peri = {pericentre:>5.1f}  ->  numerical gap = "
          f"{actual_pericentre:.3f}  (saved {output_path})")

    figure, axis = plt.subplots(figsize=(8.5, 8.0))
    figure.patch.set_facecolor("#0a0a0a")
    axis.set_facecolor("#0a0a0a")

    # Trajectories.
    axis.plot(pos_1[:, 0], pos_1[:, 1], color="#5fa8ff", lw=1.6,
              label="Galaxy 1 (M = 1)")
    axis.plot(pos_2[:, 0], pos_2[:, 1], color="#ff7a5f", lw=1.6,
              label="Galaxy 2 (M = 1)")

    # Start positions.
    axis.scatter([pos_1[0, 0]], [pos_1[0, 1]], color="#5fa8ff",
                 s=80, zorder=5)
    axis.text(pos_1[0, 0] + 1.5, pos_1[0, 1] + 1.5,
              "start (1)\n(+15, 0)", color="#5fa8ff", fontsize=14)
    axis.scatter([pos_2[0, 0]], [pos_2[0, 1]], color="#ff7a5f",
                 s=80, zorder=5)
    axis.text(pos_2[0, 0] - 8.0, pos_2[0, 1] - 3.5,
              "start (2)\n(-15, 0)", color="#ff7a5f", fontsize=14)

    # Pericentre points.
    axis.scatter([pos_1[pericentre_index, 0]],
                 [pos_1[pericentre_index, 1]],
                 marker="*", color="#5fa8ff", s=220,
                 edgecolors="#ffffff", linewidths=0.6, zorder=6)
    axis.scatter([pos_2[pericentre_index, 0]],
                 [pos_2[pericentre_index, 1]],
                 marker="*", color="#ff7a5f", s=220,
                 edgecolors="#ffffff", linewidths=0.6, zorder=6)
    # Dashed line between them showing the pericentre gap.
    axis.plot(
        [pos_1[pericentre_index, 0], pos_2[pericentre_index, 0]],
        [pos_1[pericentre_index, 1], pos_2[pericentre_index, 1]],
        color="#ffd76f", ls="--", lw=1.0,
    )
    midx = 0.5 * (pos_1[pericentre_index, 0] + pos_2[pericentre_index, 0])
    midy = 0.5 * (pos_1[pericentre_index, 1] + pos_2[pericentre_index, 1])
    axis.text(midx + 0.4, midy + 0.4,
              f"pericentre gap\nr_p = {actual_pericentre:.1f}",
              color="#ffd76f", fontsize=14)

    # COM marker.
    axis.scatter([0], [0], marker="+", color="#ffffff", s=180, lw=2.0)
    axis.text(0.6, -1.0, "COM", color="#ffffff", fontsize=14)

    # Axes / labels.
    axis.set_xlim(-view_kpc, view_kpc)
    axis.set_ylim(-view_kpc, view_kpc)
    axis.set_aspect("equal")
    axis.set_xlabel("x  [code units; 1 = 3.5 kpc]", color="#cccccc")
    axis.set_ylabel("y  [code units]", color="#cccccc")
    axis.tick_params(colors="#888888")
    for spine in axis.spines.values():
        spine.set_color("#444444")
    axis.grid(True, color="#222222", lw=0.4)
    legend = axis.legend(loc="upper left", facecolor="#101010",
                         edgecolor="#444444", fontsize=14)
    for text in legend.get_texts():
        text.set_color("#cccccc")
    extra = ""
    if pericentre == 0.0:
        extra = "  (head-on: zero angular momentum)"
    axis.set_title(
        f"Parabolic two-body encounter: separation = {separation:.0f}, "
        f"pericentre = {pericentre:.0f}{extra}\n"
        "(point-mass limit; merger pipeline starts from this orbit)",
        color="#ffffff", fontsize=18,
    )

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    figure.savefig(output_path, dpi=130, facecolor=figure.get_facecolor())
    plt.close(figure)


def main() -> None:
    """Generate the four-pericentre illustration set (or a single one)."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pericentre", type=float, default=None,
        help="render a single pericentre instead of the full set",
    )
    parser.add_argument(
        "--out-dir", type=str, default="figures/orbit_illustrations",
    )
    args = parser.parse_args()

    pericentres = (
        [args.pericentre] if args.pericentre is not None
        else [0.0, 5.0, 10.0, 20.0]
    )
    print(f"Saving orbit illustrations to {args.out_dir}/")
    for pericentre in pericentres:
        output_path = os.path.join(
            args.out_dir, f"orbit_peri{int(round(pericentre)):02d}.png"
        )
        make_one_illustration(pericentre, output_path)


if __name__ == "__main__":
    main()
