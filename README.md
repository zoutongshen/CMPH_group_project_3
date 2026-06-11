# Project 3 — Galaxy–Galaxy Collision

Self-contained, from-scratch N-body simulation of disk-galaxy mergers,
following Thorsten Naab's *N-body Simulations and Galaxy Formation*
practicum manual (`../project_3_ideas/nbody_manual.pdf`).

The code is split into three decoupled layers:

1. **Simulation** — builds initial conditions, integrates gravity, and writes
   a plain `.npz` trajectory dump. Runs entirely headless; never imports
   matplotlib.
2. **Analysis / plotting** — separate scripts that *read* a trajectory dump
   and produce figures, profiles, or animations.
3. **Batch runners** — small shell scripts that launch a family of
   simulations (parameter scans) and then the matching animations.

This separation means every simulation can be run without producing a single
plot, and every figure can be regenerated from a saved `.npz` without re-running
the simulation.

---

## Requirements & setup

The pinned environment (numpy, numba, matplotlib, tqdm) lives in a virtualenv
one level up:

```bash
# From CMPH_code/project_3/
source ../.venv/bin/activate          # Python 3.9, numba 0.60, numpy 2.0, matplotlib
```

`numba` is required (the gravity kernel and integrator step are JIT-compiled);
`tqdm` is only used for the optional progress bar; `matplotlib` is only needed
for the analysis/plotting layer.

If you prefer your own interpreter, install the four packages with
`pip install numpy numba matplotlib tqdm`.

---

## Code files

### Foundation (physics core)

| File | Purpose |
|---|---|
| `forces.py` | Direct N² gravity with Plummer softening, Numba-JIT'd |
| `integrator.py` | Velocity-Verlet (leapfrog) step + `run_simulation` driver |
| `diagnostics.py` | Kinetic / potential / total energy and angular momentum |

### Initial conditions

| File | Purpose |
|---|---|
| `initial_conditions.py` | Hernquist (1993a) equilibrium spatial sampling (disk + bulge + halo); defines `GalaxyParams` |
| `potential.py` | Tabulated spherical-approximation potential: M(<r), v_c², κ², Σ_d, Jeans dispersions |
| `velocities.py` | Velocity sampling for all three components; full single-galaxy IC builder |
| `merger_initial_conditions.py` | Two-galaxy parabolic-encounter IC builder (§3.0.5) |
| `three_galaxy_initial_conditions.py` | Three-galaxy ring-configuration IC builder (extension) |

### Validation tests (run these first)

| File | Purpose |
|---|---|
| `kepler_test.py` | 2-body Kepler validation of the gravity + integrator stack |
| `test_initial_conditions.py` | Statistical validation of the IC spatial sampler |
| `test_potential.py` | Validation of the tabulated potential against analytic limits |
| `test_velocities.py` | Validation of the velocity sampler (isotropy, rotation, anisotropy) |

### Simulation drivers (write `.npz` trajectory dumps)

| File | Purpose |
|---|---|
| `evolve_disk.py` | §3.0.4 — evolve an isolated disk, dump full trajectory |
| `evolve_merger.py` | §3.0.5 — collide two equal-mass disks |
| `evolve_unequal_merger.py` | Extension — unequal-mass merger (`--mass-ratio`) |
| `evolve_three_galaxy.py` | Extension — three-galaxy ring encounter |
| `stability_test.py` | Evolve a galaxy in isolation, measure structural drift |

### Analysis & plotting (read a `.npz`, write figures)

| File | Purpose |
|---|---|
| `disk_profiles.py` | v_c / M(<r) / Σ profiles by component, physical units (Q4/Q5/Q6) |
| `disk_visualization.py` | Face-on / edge-on images + vertical-thickening of all 3 components (Q1/Q2) |
| `merger_analysis.py` | Remnant Σ(R) / R¹ᐟ⁴ / M(<r) vs progenitor (Q7) |
| `merger_animation.py` | Face-on + edge-on merger movie (gif/mp4) |
| `merger_snapshot_grid.py` | Four-panel time sequence of one merger |
| `scan_grid_plot.py` | Grid of static 4-panels across a parameter scan |
| `gif_grid.py` | Grid of synchronised animations across a parameter scan |
| `orbit_illustration.py` | Pedagogical sketch of the two-body encounter orbit |
| `ring_diagnostic.py` | Ring-formation diagnostic for the head-on z-axis collision |

### Analytic answers & unit conversion

| File | Purpose |
|---|---|
| `physical_units.py` | Code-unit → physical-unit conversion (answers Q3) |
| `analytic_questions.py` | Worked analytic answers Q8–Q10 (disk/bulge t_dyn) |

### Batch runners (shell)

| File | Purpose |
|---|---|
| `multi_runner.sh` | Run the unequal-mass + three-galaxy extension sims |
| `zaxis_runner.sh` | Run the z-axis (face-on) collision variants |
| `extended_runner.sh` | Re-run selected orbits out to tstop = 900 (long-term fate) |
| `scan_gif_runner.sh` | Build gif animations for every available scan trajectory |

### Documentation files

- `FORMULA_REFERENCES.md` — every physics formula, with manual eq # / literature / derivation
- `NOTES_DISK_HEATING.md` — investigation log of the IC vertical-heating finding (see "Findings")
- `EXTENSIONS.md` — motivation, setup, and results for every parameter-scan extension
- `ANSWERS_Q1_Q6.md` — written answers to the §3.0.4 manual questions Q1–Q6
- `ANSWERS_Q7.md` — written answer to the §3.0.5 merger question Q7
- `ANSWERS_Q8_Q10.md` — worked answers to the §3.0.6 analytic questions Q8–Q10
- `presentation/` — the reveal.js talk (`presentation/index.html`); see `presentation/README.md`

---

## Quick start: validate the stack

These run in seconds (the first call pays a one-off Numba compile) and print
PASS/FAIL lines:

```bash
python kepler_test.py              # energy + L conservation, periodic return
python test_potential.py           # M(<r), dΦ/dr, Σ_d vs analytic
python test_velocities.py          # isotropy, rotation, dispersion anisotropy
python test_initial_conditions.py  # radial/vertical CDFs, isotropy
```

Build a single galaxy's complete IC (positions + velocities + masses) and save
it to `data/spiral.npz`:

```bash
python velocities.py
```

---

## Running the simulations

Every driver builds its own IC, integrates, prints a short summary, and writes
a trajectory `.npz` to the path given by `--out`. **No plotting happens here.**

```bash
# Isolated disk (§3.0.4) — manual's literal N = 20k
python evolve_disk.py --out data/evolve_N20k_eps01.npz

# Isolated disk at N = 80k (cleaner physics — see Findings)
python evolve_disk.py --n-disk 24000 --n-bulge 8000 --n-halo 48000 \
    --out data/evolve_N80k_eps01.npz

# Equal-mass merger (§3.0.5)
python evolve_merger.py --out data/merger_N40k_eps01.npz

# Unequal-mass merger (extension)
python evolve_unequal_merger.py --mass-ratio 0.25 --out data/multi_unequal_1to4.npz

# Three-galaxy ring encounter (extension)
python evolve_three_galaxy.py --out data/multi_three_default.npz

# Head-on, face-on (z-axis) smash that forms a ring
python evolve_merger.py --collision-axis z --inclination 0 --pericentre 0 \
    --out data/multi_zaxis_headon.npz

# Isolated-galaxy stability / structural-drift test
python stability_test.py                          # default N = 20k
python stability_test.py --save-snapshots data/stability_N20k.npz
```

To reproduce a whole family of runs at once, use the shell runners (each is
idempotent — it skips outputs that already exist):

```bash
./multi_runner.sh        # unequal-mass scan (1:2, 1:4, 1:8) + three-galaxy
./zaxis_runner.sh        # z-axis collisions at i = 0, 90, 180 + head-on
./extended_runner.sh     # long-tstop (=900) reruns of selected orbits
```

### Simulation parameters — what to change

All drivers accept these flags (defaults in brackets):

| Flag | Meaning | Default |
|---|---|---|
| `--n-disk` / `--n-bulge` / `--n-halo` | particle counts per component (**more particles → less two-body noise, slower**) | 6000 / 2000 / 12000 |
| `--tstop` | total simulated time in code units (1 ≈ 13.05 Myr) | 100 (disk) / 300 (merger) |
| `--dt` | integrator timestep (**smaller → smaller energy drift, slower**) | 0.125 |
| `--eps` | Plummer softening length (**larger → smoother, suppresses finite-N heating**) | 0.1 |
| `--seed` | RNG seed for the IC sampling (reproducibility) | 0 |
| `--num-dumps` | number of trajectory snapshots written | 100 / 200 |
| `--out` | output `.npz` path | per script |

Merger-only flags (`evolve_merger.py`, `evolve_unequal_merger.py`):

| Flag | Meaning | Default |
|---|---|---|
| `--separation` | initial galaxy separation | 30 |
| `--pericentre` | closest approach of the parabolic orbit (`0` = head-on) | 5 |
| `--inclination` | disk tilt relative to the orbital plane, degrees | 30 |
| `--mass-ratio` | secondary/primary mass (unequal merger only) | 0.5 |
| `--collision-axis` | orbital plane: `x` (edge-on) or `z` (face-on) | `x` |

Three-galaxy flags (`evolve_three_galaxy.py`): `--ring-radius` (20),
`--tangential-fraction` (0.5), `--inclination` (30).

`stability_test.py` adds `--save-snapshots <path>` to dump the diagnostic
time-series. **Note:** its energy-drift tolerance assumes the production
particle counts and `--dt 0.125`; toy values (very small N, large dt) will
legitimately exceed the tolerance and report FAIL — that is the measurement
working, not a bug.

---

## Making the plots

The plotting scripts take a trajectory `.npz` (positional argument) and write a
figure to `--out`. They never run a simulation.

```bash
# Isolated-disk deliverables
python disk_visualization.py data/evolve_N20k_eps01.npz --out figures/disk_isolation/evolve_N20k_eps01_images.png
python disk_profiles.py      data/evolve_N20k_eps01.npz --out figures/disk_isolation/evolve_N20k_eps01_profiles.png

# Merger deliverables
python merger_analysis.py     data/merger_N40k_eps01.npz --out figures/baseline_merger/merger_N40k_eps01_remnant.png
python merger_snapshot_grid.py data/merger_N40k_eps01.npz --out figures/baseline_merger/merger_N40k_eps01_4panel.png
python merger_animation.py     data/merger_N40k_eps01.npz --stride 2 --view-kpc 90 --fps 15 --out figures/baseline_merger/merger_N40k_eps01.gif

# Parameter-scan summary grids (--kind picks the family)
python scan_grid_plot.py --kind pericentre  --out figures/pericentre_scan/scan_pericentre_grid.png
python gif_grid.py       --kind inclination --out figures/inclination_scan/scan_inclination_grid.gif

# Supporting illustrations / diagnostics
python orbit_illustration.py                                   # one panel per pericentre
python ring_diagnostic.py --trajectory data/multi_zaxis_headon.npz --out figures/z_axis/multi_zaxis_headon_ring_diagnostic.png

# Build all scan animations at once
./scan_gif_runner.sh
```

### Plotting parameters — what to change

| Flag (where it applies) | Meaning |
|---|---|
| `--out` | output figure / animation path |
| `--view-kpc` | half-width of the plotted field of view |
| `--max-points` | down-sample particles drawn (speed vs. detail) |
| `--stride` | use every Nth snapshot in an animation (speed vs. smoothness) |
| `--fps` | animation frame rate |
| `--frames` (`gif_grid.py`) | number of frames in a grid animation |
| `--times` (`*_grid`, `ring_diagnostic`) | which snapshot times to show |
| `--kind` (`scan_grid_plot`, `gif_grid`) | scan family: `pericentre`, `inclination`, `multi`, `zaxis` |
| `--pericentre` (`orbit_illustration`) | render a single pericentre instead of all four |
| `--merger` (`analytic_questions`) | print the merger-specific analytic answers |

---

## Units

Hernquist (1992/1993b) system (manual page 18), with **G = 1**:

- disk:  h = 1, M_d = 1, z₀ = 0.2
- bulge: M_b = 1/3, a = 0.1
- halo:  M_h = 5.8, γ = 1, r_c = 10

Scaled to the Milky Way (printed by `python physical_units.py`):
h = 3.5 kpc, M_d = 5.6×10¹⁰ M_☉, time unit ≈ 1.31×10⁷ yr (13.05 Myr),
velocity unit ≈ 262 km/s.

---

## Findings

**Vertical heating in the isolated-disk evolution is dominated by finite-N
two-body noise**, not by an IC bug. Quantitative scan:

| N (total) | z₀ drift over t = 0 → 100 |
|---|---|
| 20,000 (manual) | 67 % |
| 80,000 | 24 % |

The 24/67 ≈ 0.36 reduction on quadrupling N matches the `t_relax ∝ N / ln N`
scaling (predicted 0.28). The radial scale length h stays essentially fixed
(< 1 % drift) regardless, because angular momentum is conserved at machine
precision. Full log in `NOTES_DISK_HEATING.md`.

**Decision:** for the exercise runs (§3.0.4 and §3.0.5) we present *both* N = 20k
(manual's literal prescription) and N = 80k (clean physics), so the
N-dependence is shown as a deliberate science result rather than an artefact.

---

## Status

| Step | State |
|---|---|
| Foundation (forces, integrator, diagnostics) | ✓ Kepler test passes |
| IC spatial sampler | ✓ CDF-match tests pass |
| Tabulated potential | ✓ Analytic-limit tests pass |
| Velocity sampler | ✓ Isotropy / rotation / anisotropy tests pass |
| IC isolated-disk stability test | ✓ Finite-N heating characterised |
| §3.0.4 isolated disk + Q1–Q6 | ✓ See `ANSWERS_Q1_Q6.md` |
| §3.0.5 merger + Q7 | ✓ See `ANSWERS_Q7.md` |
| §3.0.6 analytic questions Q8–Q10 | ✓ See `ANSWERS_Q8_Q10.md` |
| Extensions (pericentre / inclination / mass-ratio / z-axis / three-galaxy) | ✓ See `EXTENSIONS.md` |
