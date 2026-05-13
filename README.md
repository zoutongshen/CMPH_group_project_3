# Project 3 — Galaxy-Galaxy Collision

Self-contained from-scratch N-body simulation of a disk-galaxy merger, following
Thorsten Naab's *N-body Simulations and Galaxy Formation* practicum manual
(`../project_3_ideas/nbody_manual.pdf`).

## Files

- `forces.py` — Direct N² gravity with Plummer softening (Numba-JIT'd)
- `integrator.py` — Velocity-Verlet (leapfrog) integrator + simulation driver
- `diagnostics.py` — Total / kinetic / potential energy and angular momentum
- `kepler_test.py` — 2-body Kepler validation (run as a script)
- `data/` — Output snapshots
- `figures/` — Plots

## Units

Hernquist 1992/1993b system, as in the manual:

- G = 1
- Exponential disk scale length h = 1
- Larger-galaxy disk mass M_d = 1
- z₀ = 0.2, M_h = 5.8, γ = 1, r_c = 10, M_b = 1/3, a = 0.1

Scaled to the Milky Way: h = 3.5 kpc, M_d = 5.6×10¹⁰ M_☉, time unit = 1.31×10⁷ yr,
velocity unit = 262 km/s.

## How to run

```bash
# Activate the venv (numba 0.60, numpy 2.0, matplotlib 3.9)
source ../.venv/bin/activate

# Validate the integrator
python kepler_test.py
```
