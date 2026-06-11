"""
Compare the circular-velocity (rotation) curves of the two galaxy models:

  * Manual prescription (this project's reference run): M_d = 1, M_b = 1/3,
    M_h = 5.8 -- the masses on page 18 of the practicum manual.
  * Octree tree-code model: M_d = 1, M_b = 0.25, M_h = 3.0 -- the masses
    hard-coded in galaxy_initialization.py.

The rotation curve is fixed by the mass distribution alone, so it can be read
straight off the initial conditions -- no dynamical evolution is needed. For
each model we sort particles by spherical radius, form the cumulative enclosed
mass M(<r), and evaluate v_c(r) = sqrt(G M(<r) / r). Both curves use the SAME
estimator, so the only thing that moves between them is the mass model.

Code-unit -> physical conversion is the manual's Milky-Way scaling
(h = 3.5 kpc, v unit = 262 km/s), identical for both models since both fix
M_d = 1 and h = 1.

Run from the galaxy_merger_octree/ directory:
    python compare_rotation_curves.py
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import galaxy_initialization as gi

G = 1.0
LENGTH_KPC = 3.5      # code length unit -> kpc
VELOCITY_KMS = 262.0  # code velocity unit -> km/s
MW_VC_KMS = 220.0     # observed Milky-Way flat rotation speed (reference)


def circular_velocity(positions, masses, r_grid_code):
    """v_c(r) = sqrt(G M(<r) / r) from the spherical cumulative mass profile."""
    radius = np.linalg.norm(positions, axis=1)
    order = np.argsort(radius)
    radius_sorted = radius[order]
    mass_cumulative = np.cumsum(masses[order])
    enclosed_mass = np.interp(r_grid_code, radius_sorted, mass_cumulative)
    return np.sqrt(G * enclosed_mass / r_grid_code)


def build_octree_model_galaxy(seed=0):
    """Sample the tree-code model galaxy (his masses) and return pos, mass."""
    np.random.seed(seed)
    n_bulge, n_disk, n_halo = 2000, 6000, 12000
    mass_bulge, mass_disk, mass_halo = 0.25, 1.0, 3.0

    pos_bulge, _ = gi.generate_bulge_positions(n_bulge)
    pos_disk, _ = gi.generate_disk_positions(n_disk)
    pos_halo, _ = gi.generate_halo_positions(n_halo)

    positions = np.vstack((pos_bulge, pos_disk, pos_halo))
    masses = np.concatenate((
        np.full(n_bulge, mass_bulge / n_bulge),
        np.full(n_disk, mass_disk / n_disk),
        np.full(n_halo, mass_halo / n_halo),
    ))
    return positions, masses


def load_manual_model_galaxy(npz_path):
    """Load the reference run's t = 0 initial conditions (manual masses)."""
    data = np.load(npz_path)
    positions = data["positions_history"][0].astype(np.float64)
    masses = data["masses"].astype(np.float64)
    # Recentre on the mass-weighted centroid.
    centre = (positions * masses[:, None]).sum(axis=0) / masses.sum()
    return positions - centre, masses


def main():
    pos_octree, mass_octree = build_octree_model_galaxy()
    pos_manual, mass_manual = load_manual_model_galaxy(
        "../data/evolve_N80k_eps01.npz")

    print(f"Manual model  total mass = {mass_manual.sum():.2f} "
          "(M_d=1, M_b=1/3, M_h=5.8)")
    print(f"Octree model  total mass = {mass_octree.sum():.2f} "
          "(M_d=1, M_b=0.25, M_h=3.0)")

    r_grid = np.linspace(0.15, 15.0, 400)        # code units (~0.5 - 52 kpc)
    vc_manual = circular_velocity(pos_manual, mass_manual, r_grid)
    vc_octree = circular_velocity(pos_octree, mass_octree, r_grid)

    r_kpc = r_grid * LENGTH_KPC
    vc_manual_kms = vc_manual * VELOCITY_KMS
    vc_octree_kms = vc_octree * VELOCITY_KMS

    print(f"Manual model  peak v_c = {vc_manual_kms.max():.0f} km/s")
    print(f"Octree model  peak v_c = {vc_octree_kms.max():.0f} km/s")

    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(9, 6))

    ax.axhline(MW_VC_KMS, color="#888888", ls=":", lw=2,
               label=f"Milky Way (~{MW_VC_KMS:.0f} km/s)")
    ax.plot(r_kpc, vc_manual_kms, color="#4da6ff", lw=3,
            label=r"Manual model  ($M_h=5.8$)")
    ax.plot(r_kpc, vc_octree_kms, color="#ff7b4d", lw=3,
            label=r"Octree model  ($M_h=3.0$)")

    ax.set_xlabel("r  [kpc]", fontsize=18)
    ax.set_ylabel(r"$v_c(r)$  [km/s]", fontsize=18)
    ax.set_title("Circular-velocity curve: mass-model comparison", fontsize=19)
    ax.tick_params(labelsize=14)
    ax.set_xlim(0, r_kpc.max())
    ax.set_ylim(0, None)
    ax.legend(fontsize=15, frameon=False)

    fig.tight_layout()
    import os
    os.makedirs("../figures/comparison", exist_ok=True)
    out_path = "../figures/comparison/rotation_curve_compare.png"
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
