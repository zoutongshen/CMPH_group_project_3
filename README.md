# Project 3 — Galaxy-Galaxy Collision

Self-contained from-scratch N-body simulation of a disk-galaxy merger,
following Thorsten Naab's *N-body Simulations and Galaxy Formation*
practicum manual (`../project_3_ideas/nbody_manual.pdf`).

## Code files

| File | Purpose |
|---|---|
| `forces.py` | Direct N² gravity with Plummer softening, Numba-JIT'd |
| `integrator.py` | Velocity-Verlet (leapfrog) integrator + `run_simulation` driver |
| `diagnostics.py` | Kinetic / potential / total energy and angular momentum |
| `initial_conditions.py` | Hernquist 1993a equilibrium IC spatial sampling (disk + bulge + halo) |
| `potential.py` | Tabulated spherical-approximation potential: M(<r), v_c², κ², Σ_d, Jeans dispersions |
| `velocities.py` | IC velocity sampling for all three components; full-IC builder |
| `stability_test.py` | Run isolated galaxy in isolation, measure structural drift |
| `physical_units.py` | Code-unit → physical-unit conversion (answers Q3) |
| `evolve_disk.py` | §3.0.4 driver: evolve isolated disk, dump full trajectory |
| `disk_profiles.py` | v_c / M(<r) / Σ profiles by component, physical units (Q4/Q5/Q6) |
| `disk_visualization.py` | Face-on / edge-on images of all 3 components (Q1/Q2) |
| `merger_initial_conditions.py` | Two-galaxy parabolic-encounter IC builder (§3.0.5) |
| `evolve_merger.py` | §3.0.5 driver: collide two disks, dump full trajectory |
| `merger_analysis.py` | Remnant Σ(R) / R¹ᐟ⁴ / M(<r) vs progenitor (Q7) |
| `merger_animation.py` | Face-on + edge-on merger movie (mp4/gif) |
| `analytic_questions.py` | Worked analytic answers Q8–Q10 (disk/bulge t_dyn) |
| `kepler_test.py` | 2-body Kepler validation of the gravity + integrator stack |
| `test_initial_conditions.py` | Statistical validation of the IC spatial sampler |
| `test_potential.py` | Validation of the tabulated potential against analytic limits |
| `test_velocities.py` | Validation of the velocity sampler (isotropy, rotation, anisotropy) |

## Documentation files

- `FORMULA_REFERENCES.md` — every physics formula, with manual eq # / original literature / derivation
- `NOTES_DISK_HEATING.md` — investigation log of the IC vertical heating finding (see "Findings" below)
- `ANSWERS_Q1_Q6.md` — written answers to the §3.0.4 manual questions Q1–Q6
- `ANSWERS_Q7.md` — written answer to the §3.0.5 merger question Q7
- `ANSWERS_Q8_Q10.md` — worked answers to the §3.0.6 analytic questions Q8–Q10

## Units

Hernquist 1992/1993b system (manual page 18):

- G = 1
- disk: h = 1, M_d = 1, z₀ = 0.2
- bulge: M_b = 1/3, a = 0.1
- halo: M_h = 5.8, γ = 1, r_c = 10

Scaled to the Milky Way: h = 3.5 kpc, M_d = 5.6×10¹⁰ M_☉,
time unit ≈ 1.31×10⁷ yr, velocity unit ≈ 262 km/s.

## How to run

```bash
# Activate the project venv (numba, numpy, matplotlib already installed)
source ../.venv/bin/activate

# Foundation validation
python kepler_test.py

# IC component validation
python test_initial_conditions.py
python test_potential.py
python test_velocities.py

# Build the full IC and save to data/spiral.npz
python velocities.py

# Isolated-disk stability test
python stability_test.py                                       # default: manual prescription (N=20k)
python stability_test.py --n-disk 24000 --n-bulge 8000 --n-halo 48000   # N=80k for cleaner physics
```

CLI flags for `stability_test.py`:
`--n-disk`, `--n-bulge`, `--n-halo`, `--tstop`, `--dt`, `--eps`, `--seed`,
`--save-snapshots <path>`.

## Current status

| Step | State |
|---|---|
| Foundation (forces, integrator, diagnostics) | ✓ Kepler test passes |
| IC spatial sampler | ✓ CDF-match tests pass |
| Tabulated potential | ✓ Analytic-limit tests pass |
| Velocity sampler | ✓ Isotropy / rotation / anisotropy tests pass |
| IC isolated-disk stability test | ✓ Done — finite-N heating characterised (see findings) |
| §3.0.4 isolated disk evolution + Q1–Q6 | ✓ Done — see `ANSWERS_Q1_Q6.md` |
| §3.0.5 merger + Q7 | ✓ Done — see `ANSWERS_Q7.md` |
| §3.0.6 additional analytic questions | ✓ Done — see `ANSWERS_Q8_Q10.md` |

## Findings so far

**Vertical heating in the isolated-disk evolution is dominated by finite-N
two-body noise**, not by any IC bug. Quantitative scan:

| N (total) | z₀ drift over t=0→100 |
|---|---|
| 20,000 (manual) | 67 % |
| 80,000 | 24 % |

The 24/67 ≈ 0.36 reduction quadrupling N matches the
`t_relax ∝ N / ln N` scaling (predicted 0.28). The radial scale length h
stays essentially fixed (< 1 % drift) regardless, because angular momentum
is conserved at machine precision. Full investigation log in
`NOTES_DISK_HEATING.md`.

**Decision (taken 2026-05-13):** for the actual exercise runs (§3.0.4 and
§3.0.5), run at *both* N=20k (manual's literal prescription) and N=80k
(clean physics) so we can present the N-dependence as a deliberate
science result rather than an artefact.

## Repository

Hosted at https://github.com/zoutongshen/CMPH_group_project_3 (private).
