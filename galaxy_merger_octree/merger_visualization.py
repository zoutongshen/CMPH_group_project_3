"""
Render an animated GIF of a two-galaxy merger from saved snapshots.

Reads every galaxy_step_*.npz file under data/snapshot_merger, lays out two
matplotlib axes (face-on x-y and edge-on x-z), and animates the particle
positions over time using matplotlib.animation.FuncAnimation. Labels in the
merger snapshots use six classes so each galaxy keeps its own colours
through the encounter:

    0 = galaxy 1 bulge   3 = galaxy 2 bulge
    1 = galaxy 1 disk    4 = galaxy 2 disk
    2 = galaxy 1 halo    5 = galaxy 2 halo

Galaxy 1 is shown in warm tones (red bulge, white disk, cyan halo) and
galaxy 2 in cool/contrasting tones (orange bulge, yellow disk, magenta
halo). The species draw order is set via zorder so the small central
bulges remain visible above the diffuse halo cloud, and the axis limits
are wider than the bulk of the distribution so the galaxies never touch
the frame even at the widest separation.

Requires Pillow for the GIF writer: pip install Pillow.

Author: CMPH 2026 Project 3
Date: May 2026
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import glob

# Run-level parameters.
data_dir = "data/snapshot_merger"             # directory holding the snapshots
output_gif = "figure/merger_evolution.gif"    # output GIF filename
fps = 24                                      # animation frame rate
sample_every = 1                              # frame stride (raise to thin a dense run)
dpi = 100                                     # GIF resolution (higher dpi = larger file)

# Per-species colours, legend names, and point sizes for both galaxies.
# Galaxy 1 uses the same palette as the single-galaxy visualisation; galaxy
# 2 uses contrasting hues so the two systems remain distinguishable as they
# interpenetrate.
label_colors = {
    0: 'red',     1: 'white',  2: 'cyan',
    3: 'orange',  4: 'yellow', 5: 'magenta',
}
label_names = {
    0: 'Primary Bulge', 1: 'Primary Disk', 2: 'Primary Halo',
    3: 'Intruder Bulge', 4: 'Intruder Disk', 5: 'Intruder Halo',
}
point_sizes = {
    0: 1.5, 1: 0.8, 2: 0.2,
    3: 1.5, 4: 0.8, 5: 0.2,
}

# Per-species alpha. Halo (labels 2 and 5) is drawn faint so the diffuse
# dark-matter cloud does not wash out the disk's tidal features.
point_alphas = {
    0: 0.8, 1: 0.7, 2: 0.10,
    3: 0.8, 4: 0.7, 5: 0.10,
}

# Plot limits. The collision is along z (intruder travels from +z to -z
# at v ~ 1.2, so by t = 100 it can be at z ~ -100), but the ring of
# interest stays inside the primary's disk extent ~ +/-15 in x-y. Keep
# the x-y limits tight to show the ring clearly, and let z span a wider
# range so the intruder's trajectory stays in frame for longer.
xlim = (-25, 25)
ylim = (-25, 25)
zlim = (-60, 30)


# ---------------------------------------------------------------------------
# Collect snapshots and sort them by step index.
# ---------------------------------------------------------------------------

snapshot_files = sorted(glob.glob(os.path.join(data_dir, "galaxy_step_*.npz")))
steps = []
for f in snapshot_files:
    try:
        # File name pattern: galaxy_step_NNNN.npz; extract the integer step.
        step = int(os.path.basename(f).split('_')[2].split('.')[0])
        steps.append((step, f))
    except Exception:
        continue
steps.sort(key=lambda x: x[0])
steps = steps[::sample_every]

if len(steps) == 0:
    raise FileNotFoundError(
        f"No galaxy_step_*.npz files found in {data_dir}.")

print(f"Found {len(steps)} snapshot files (after sampling).")

# Labels are constant across the run: read them once from the first frame.
first_data = np.load(steps[0][1])
labels = first_data["label"]


# ---------------------------------------------------------------------------
# Build the figure: dark background and per-species scatter handles.
# ---------------------------------------------------------------------------

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
fig.patch.set_facecolor('black')
fig.suptitle("Galaxy Merger", fontsize=14, color='white')

# Draw order: halos first (largest, most diffuse), then disks, then bulges
# on top. Use the same relative ordering for both galaxies so neither
# bulge is hidden behind the other galaxy's halo.
scatters1 = {}   # face-on (x-y)
scatters2 = {}   # edge-on (x-z)
zorders = {2: 1, 5: 1, 1: 2, 4: 2, 0: 3, 3: 3}
draw_order = [2, 5, 1, 4, 0, 3]   # back-to-front
for label in draw_order:
    sc1 = ax1.scatter([], [], s=point_sizes[label], c=label_colors[label],
                      alpha=point_alphas[label], label=label_names[label],
                      zorder=zorders[label])
    sc2 = ax2.scatter([], [], s=point_sizes[label], c=label_colors[label],
                      alpha=point_alphas[label], zorder=zorders[label])
    scatters1[label] = sc1
    scatters2[label] = sc2


def _style_axis(ax):
    """Apply the dark-background colour scheme to a single axes object."""
    ax.set_facecolor('black')
    ax.tick_params(colors='white', which='both')
    for spine in ax.spines.values():
        spine.set_color('white')
    ax.title.set_color('white')
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')


ax1.set_title("Primary face-on (x-y)")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_xlim(xlim)
ax1.set_ylim(ylim)
ax1.set_aspect("equal")
_style_axis(ax1)
leg = ax1.legend(loc="upper right", fontsize=7, facecolor='black',
                 edgecolor='white', labelcolor='white', ncol=2)

ax2.set_title("Collision axis (x-z)")
ax2.set_xlabel("x")
ax2.set_ylabel("z")
ax2.set_xlim(xlim)
ax2.set_ylim(zlim)
ax2.set_aspect("equal")
_style_axis(ax2)

# Time stamp drawn in the top-left corner of the face-on panel.
time_text = ax1.text(0.02, 0.95, "", transform=ax1.transAxes, fontsize=10,
                     color='white')


# ---------------------------------------------------------------------------
# Per-frame update: read one snapshot and push the new offsets to scatters.
# ---------------------------------------------------------------------------

def update(frame_idx):
    """Refresh the scatter offsets and time label for one animation frame."""
    step, filepath = steps[frame_idx]
    data = np.load(filepath)
    pos = data["pos"]
    for label in range(6):
        mask = (labels == label)
        x = pos[mask, 0]
        y = pos[mask, 1]
        z = pos[mask, 2]
        scatters1[label].set_offsets(np.c_[x, y])   # face-on (x-y)
        scatters2[label].set_offsets(np.c_[x, z])   # edge-on (x-z)
    # Assumes dt = 0.01; matches the merger_simulation default.
    time_text.set_text(f"Time = {step * 0.01:.2f}")
    return list(scatters1.values()) + list(scatters2.values()) + [time_text]


# ---------------------------------------------------------------------------
# Render the animation and save to disk.
# ---------------------------------------------------------------------------

ani = animation.FuncAnimation(fig, update, frames=len(steps),
                              interval=1000 / fps, blit=True)

# PillowWriter defaults to a white frame background; pass the figure's
# facecolor through savefig_kwargs so the saved GIF stays black.
os.makedirs(os.path.dirname(output_gif), exist_ok=True)
writer = animation.PillowWriter(fps=fps)
ani.save(output_gif, writer=writer, dpi=dpi,
         savefig_kwargs={'facecolor': fig.get_facecolor()})
print(f"GIF animation saved to {output_gif}")

# Optional debug helper: show only the last frame as a static figure.
# update(len(steps) - 1)
# plt.show()
