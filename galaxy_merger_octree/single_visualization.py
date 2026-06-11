"""
Render an animated GIF of the galaxy evolution from saved snapshots.

Reads every galaxy_step_*.npz file under data/snapshot_single, lays out two
matplotlib axes (face-on x-y and edge-on x-z), and animates the particle
positions over time using matplotlib.animation.FuncAnimation. The figure
uses a dark background and per-species colouring (bulge red, disk white,
halo cyan) reminiscent of a real galaxy photograph.

The species draw order is set via zorder so that bulge points (small but
densely packed at the centre) are not occluded by the larger halo cloud.
The visualisation axis limits are wider than the bulk of the particle
distribution so the galaxy never touches the frame edge.

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
data_dir = "data/snapshot_single"      # directory holding the snapshots
output_gif = "figure/galaxy_evolution.gif"    # output GIF filename
fps = 10                               # animation frame rate
sample_every = 1                       # frame stride (raise to thin a dense run)
dpi = 100                              # GIF resolution (higher dpi = larger file)

# Per-species colours, legend names, and point sizes (labels: 0=bulge,
# 1=disk, 2=halo).
label_colors = {0: 'red', 1: 'white', 2: 'cyan'}
label_names = {0: 'Bulge', 1: 'Disk', 2: 'Halo'}
point_sizes = {0: 1.5, 1: 0.8, 2: 0.2}

# Plot limits, padded out so end-of-run disk/halo never touch the frame.
xlim = (-20, 20)
ylim = (-20, 20)
zlim = (-12, 12)


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
fig.suptitle("Galaxy Evolution", fontsize=14, color='white')

# One scatter handle per species so the per-frame update can just refresh
# their offsets. Draw halo first, disk on top, bulge highest, so the small
# central bulge stays visible above the diffuse halo cloud.
scatters1 = {}   # face-on (x-y)
scatters2 = {}   # edge-on (x-z)
zorders = {2: 1, 1: 2, 0: 3}
for label in [2, 1, 0]:
    sc1 = ax1.scatter([], [], s=point_sizes[label], c=label_colors[label],
                      alpha=0.6, label=label_names[label],
                      zorder=zorders[label])
    sc2 = ax2.scatter([], [], s=point_sizes[label], c=label_colors[label],
                      alpha=0.6, zorder=zorders[label])
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


ax1.set_title("Face-on (x-y)")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_xlim(xlim)
ax1.set_ylim(ylim)
ax1.set_aspect("equal")
_style_axis(ax1)
leg = ax1.legend(loc="upper right", fontsize=8, facecolor='black',
                 edgecolor='white', labelcolor='white')

ax2.set_title("Edge-on (x-z)")
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
    for label in [0, 1, 2]:
        mask = (labels == label)
        x = pos[mask, 0]
        y = pos[mask, 1]
        z = pos[mask, 2]
        scatters1[label].set_offsets(np.c_[x, y])   # face-on (x-y)
        scatters2[label].set_offsets(np.c_[x, z])   # edge-on (x-z)
    # Assumes dt = 0.01; matches the simulation default.
    time_text.set_text(f"Time = {step * 0.01:.2f}")
    return list(scatters1.values()) + list(scatters2.values()) + [time_text]


# ---------------------------------------------------------------------------
# Render the animation and save to disk.
# ---------------------------------------------------------------------------

ani = animation.FuncAnimation(fig, update, frames=len(steps),
                              interval=1000 / fps, blit=True)

# PillowWriter defaults to a white frame background; pass the figure's
# facecolor through savefig_kwargs so the saved GIF stays black.
writer = animation.PillowWriter(fps=fps)
ani.save(output_gif, writer=writer, dpi=dpi,
         savefig_kwargs={'facecolor': fig.get_facecolor()})
print(f"GIF animation saved to {output_gif}")

# Optional debug helper: show only the last frame as a static figure.
# update(len(steps) - 1)
# plt.show()
