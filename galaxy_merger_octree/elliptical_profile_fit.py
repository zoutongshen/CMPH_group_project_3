"""
Fit a de Vaucouleurs R^(1/4) surface-density profile to the merger remnant.

Reads the last n_stack snapshots in data/snapshot_merger/, keeps only
the luminous component (bulge + disk; the dark halo is not visible in
real galaxies), bins particles by projected radius on the x-y plane,
and linearly fits

    log10(Sigma) = a + b * R^(1/4)

i.e. the straight-line form of the de Vaucouleurs law.

Tunable parameters at the top of the file:
    snapshot_dir    -- directory holding galaxy_step_*.npz
    output_path     -- where the figure is written
    n_stack         -- number of trailing snapshots to average over
    luminous_labels -- label values counted as luminous matter
    n_bins          -- number of radial bins
    fit_r_range     -- (R_min, R_max) for the linear fit window
"""

import os
import glob

import numpy as np
import matplotlib.pyplot as plt


snapshot_dir = "data/snapshot_merger"
output_path = "figure/elliptical_profile.png"
n_stack = 10
luminous_labels = (0, 1, 3, 4)
n_bins = 25
fit_r_range = (0.3, 15.0)


def main():
    files = sorted(glob.glob(os.path.join(snapshot_dir, "galaxy_step_*.npz")))
    if not files:
        raise FileNotFoundError(f"No snapshots found in {snapshot_dir}")
    files = files[-n_stack:]
    print(f"Averaging over the last {len(files)} snapshots in {snapshot_dir}")

    # Stack luminous particles from each snapshot, re-centring each on
    # its own COM so any residual drift does not smear the profile.
    R_all = []
    mass_all = []
    for f in files:
        snap = np.load(f)
        mask = np.isin(snap["label"], luminous_labels)
        pos = snap["pos"][mask].astype(np.float64)
        mass = snap["mass"][mask].astype(np.float64)
        com = np.average(pos, axis=0, weights=mass)
        pos -= com
        # Projected radius on x-y plane.
        R_all.append(np.sqrt(pos[:, 0] ** 2 + pos[:, 1] ** 2))
        mass_all.append(mass)
    R = np.concatenate(R_all)
    masses = np.concatenate(mass_all) / len(files)   # per-snapshot weight

    # Log-spaced annuli; surface density = mass in annulus / annulus area.
    edges = np.logspace(np.log10(0.1), np.log10(fit_r_range[1] * 2.0),
                        n_bins + 1)
    centres = np.sqrt(edges[:-1] * edges[1:])
    mass_in, _ = np.histogram(R, bins=edges, weights=masses)
    area = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    sigma = mass_in / area

    # Linear fit of log10(Sigma) vs R^(1/4) over the fit window.
    fit_mask = ((centres >= fit_r_range[0])
                & (centres <= fit_r_range[1])
                & (sigma > 0))
    x_fit = centres[fit_mask] ** 0.25
    y_fit = np.log10(sigma[fit_mask])
    slope, intercept = np.polyfit(x_fit, y_fit, 1)
    print(f"Fit: log10(Sigma) = {intercept:.3f} + ({slope:.3f}) * R^(1/4)")

    # Plot.
    fig, ax = plt.subplots(figsize=(7, 6))
    pos_mask = sigma > 0
    ax.plot(centres[pos_mask] ** 0.25, np.log10(sigma[pos_mask]),
            "o", color="k", label="N-body profile")
    x_line = np.linspace(x_fit.min(), x_fit.max(), 100)
    ax.plot(x_line, intercept + slope * x_line,
            "-", color="crimson", label=r"$R^{1/4}$ fit")
    ax.set_xlabel(r"$R^{1/4}$")
    ax.set_ylabel(r"$\log_{10}\,\Sigma$")
    ax.set_title(f"Merger remnant: averaged over last {len(files)} snapshots")
    ax.legend()

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    print(f"Figure saved to {output_path}")


if __name__ == "__main__":
    main()
