# Project 3 Extensions: parameter-space exploration of galaxy mergers

Everything beyond the manual's Q1-Q11. Each extension varies one parameter
of the standard two-galaxy merger setup; all share the same baseline so
the only knob that moves is the one being studied.

---

## Common baseline

The reference run that all extensions deviate from:

| parameter | value | physical |
|---|---|---|
| N | 40,000 (20k per galaxy = 6k disk + 2k bulge + 12k halo) | -- |
| eps (softening) | 0.1 | 0.35 kpc |
| dt | 0.125 | 1.63 Myr |
| tstop | 300 | 3.91 Gyr |
| seed | 0 | -- |
| separation | 30 | 105 kpc |
| pericentre | 5 | 17.5 kpc |
| inclination | 30 deg | -- |
| collision axis | x | -- |
| mass ratio | 1 : 1 | both galaxies use the same `GalaxyParams` |

Reference data file: `data/merger_N40k_eps01.npz`
Reference figures:   `figures/baseline_merger/merger_N40k_eps01_4panel.png`,
                     `figures/baseline_merger/merger_N40k_eps01.gif`,
                     `figures/baseline_merger/merger_N40k_eps01_remnant.png`

All extensions reproduce via runners checked into the repo
(`scan_runner.sh` on branch `extension/orbit-scans`,
`multi_runner.sh` + `zaxis_runner.sh` + `extended_runner.sh` on
branch `extension/multi-galaxy`).

---

## Extension 1: Pericentre scan

**Motivation.** The pericentre r_p controls how close the galaxies come
on first passage, and therefore (i) how much tidal disturbance they
inflict on each other and (ii) how much orbital energy dynamical
friction can extract. Sweeping r_p maps the boundary between **merger**
and **flyby**.

**Setup.** Same as baseline, vary only `--pericentre`.

**Grid.**

| label | r_p [code] | r_p [kpc] | tstop | file |
|---|---|---|---|---|
| scan_peri1   | 1   | 3.5  | 300 | `data/scan_peri1.npz` |
| baseline     | 5   | 17.5 | 300 | `data/merger_N40k_eps01.npz` |
| scan_peri10  | 10  | 35   | 300 | `data/scan_peri10.npz` |
| scan_peri10_long | 10 | 35 | 900 | `data/scan_peri10_long.npz` |
| scan_peri20  | 20  | 70   | 500 | `data/scan_peri20.npz` |
| scan_peri20_long | 20 | 70 | 900 | `data/scan_peri20_long.npz` |

Orbital geometry: `figures/orbit_illustrations/orbit_peri{00,05,10,20}.png`
(point-mass parabolic orbits, COM frame; for the eye to see the
relationship between r_p and orbit curvature).

**Figures.** `figures/pericentre_scan/`
- 4-panel grids: `scan_peri{1,10,20}_4panel.png`, `_long_4panel.png`
- Individual gifs: `scan_peri{1,10,20}.gif`, `_long.gif`
- Master comparison: `scan_pericentre_grid.{png,gif}`

**Results (outline -- flesh out tomorrow).**
- [ ] r_p = 1 (near head-on): rapid merger by t ~ 60, strongest tidal disruption, longest tails
- [ ] r_p = 5 (default): clean merger by t ~ 120, classic Toomre-tail morphology
- [ ] r_p = 10: extended interaction, eventual merger by t ~ 600-700 (visible in `_long_4panel.png` panel 4)
- [ ] r_p = 20: definitive flyby/escape -- galaxies leave the 90 kpc viewport before t = 300 and never return by t = 900
- [ ] **Threshold**: somewhere between r_p = 10 and r_p = 20 lies the merger / no-merger boundary at this mass / dynamical-friction regime
- [ ] Connection to orbital theory: friction efficiency scales with pericentre passage time; wider pericentre = faster pass = less energy extracted

**Presentation hook.** "How close do they need to come to merge?"
Three outcomes (rapid, delayed, escape) in one parameter sweep.

---

## Extension 2: Inclination scan

**Motivation.** The relative orientation of the two disks vs the orbital
plane controls the prograde/retrograde geometry of disk material at
pericentre, which Toomre & Toomre 1972 showed strongly affects tidal-tail
formation. Sweeping inclination maps the morphological zoo.

**Setup.** Same as baseline, vary only `--inclination` (rotation of disk 2
about the y-axis, so 0 deg = coplanar, 90 deg = polar, 180 deg = anti-aligned).

**Grid.**

| label | inclination [deg] | physical interpretation | file |
|---|---|---|---|
| scan_incl0   | 0   | coplanar prograde   | `data/scan_incl0.npz` |
| baseline     | 30  | classic Toomre setup | `data/merger_N40k_eps01.npz` |
| scan_incl60  | 60  | mildly tilted       | `data/scan_incl60.npz` |
| scan_incl90  | 90  | polar (disk 2 perpendicular to orbit plane) | `data/scan_incl90.npz` |
| scan_incl180 | 180 | anti-aligned (retrograde disk 2) | `data/scan_incl180.npz` |

**Figures.** `figures/inclination_scan/`
- 4-panel grids: `scan_incl{0,60,90,180}_4panel.png`
- Individual gifs: `scan_incl{0,60,90,180}.gif`
- Master comparison: `scan_inclination_grid.{png,gif}`

**Results (outline -- flesh out tomorrow).**
- [ ] incl 0 (coplanar): prograde-prograde, longest and most symmetric tidal tails
- [ ] incl 30 (default): asymmetric tails, classic Toomre morphology
- [ ] incl 60-90: tails progressively shorter / more diffuse as one disk's spin moves out of the orbit plane
- [ ] incl 180 (retrograde): minimal tidal response from disk 2 -- retrograde encounters do NOT raise prominent tails (Toomre's original result)
- [ ] **Punchline**: the spectacular tidal-tail systems we see in real galaxies (NGC 4676 "Mice", NGC 4038/9 "Antennae") require prograde geometry; retrograde mergers look much less dramatic

**Presentation hook.** "Why don't all mergers look like the Antennae?"
Answer: most don't have the right disk orientation.

---

## Extension 3: Mass-ratio scan

**Motivation.** Real galaxy mergers are rarely 1:1. Sweeping the mass
ratio maps from "major merger" (1:1, both disks destroyed, R^(1/4)
remnant per Q7) to "minor merger" (1:8, secondary shredded and
absorbed by an intact primary).

**Setup.** Same as baseline, but galaxy 2's `GalaxyParams` are scaled
by the mass ratio (all three components -- disk, bulge, halo masses
AND particle counts -- scale proportionally so the structural ratios
are preserved). The orbit's mass-weighted COM split kicks in too:
the lighter galaxy moves much more than the heavier one.

**Grid.**

| label | mass ratio (M_2 / M_1) | secondary N_particles | file |
|---|---|---|---|
| baseline (1:1) | 1.0   | 20,000 | `data/merger_N40k_eps01.npz` |
| 1to2 | 0.5   | 10,000 | `data/multi_unequal_1to2.npz` |
| 1to4 | 0.25  | 5,000  | `data/multi_unequal_1to4.npz` |
| 1to8 | 0.125 | 2,500  | `data/multi_unequal_1to8.npz` |
| 1to8_long | 0.125 | 2,500 | `data/multi_unequal_1to8_long.npz` (tstop = 900) |

Notation convention: "1to8" means M_2 = M_1 / 8 (galaxy 2 is the
secondary, galaxy 1 is the primary).

**Figures.** `figures/mass_ratio_scan/`
- 4-panel grids: `multi_unequal_1to{2,4,8}_4panel.png`, `_long_4panel.png`
- Individual gifs: `multi_unequal_1to{2,4,8}.gif`, `_long.gif`
- Master comparison: `multi_unequal_grid.{png,gif}`

**Results (outline -- flesh out tomorrow).**
- [ ] 1:1 (baseline): both progenitors destroyed, R^(1/4) elliptical remnant (Q7)
- [ ] 1:2: still classifiable as major; both disks heavily disrupted
- [ ] 1:4: borderline; primary survives mostly intact, secondary thoroughly disrupted
- [ ] 1:8 (minor merger): primary essentially untouched, secondary tidally shredded; debris spirals in. By t = 800 the secondary is fully absorbed.
- [ ] **Trend**: the "merger" smoothly becomes "accretion" as mass ratio drops; tidal structure in primary scales with secondary mass
- [ ] Connection to galaxy assembly: minor mergers dominate the cosmological merger rate, drive bulge growth without destroying the disk

**Presentation hook.** "Major vs minor mergers" --
the same physics smoothly interpolates between two qualitatively
different observational regimes.

---

## Extension 4: Three-galaxy merger

**Motivation.** Multi-body chaos. With two galaxies on a parabolic orbit
the dynamics are clean; adding a third immediately introduces
three-body chaos, no closed-form orbit, and very different morphological
outcomes depending on which two interact first.

**Setup.** Three equal-mass galaxies placed on an equilateral triangle
with sub-Keplerian tangential velocity. Built by
`three_galaxy_initial_conditions.py`; driven by `evolve_three_galaxy.py`.
Tangential speed v_circ = sqrt(M_galaxy / (R * sqrt(3))) -- slightly
under-bound so they spiral inward.

**Grid.** Only one run; the parameter is "three galaxies, default
geometry."

| label | description | file |
|---|---|---|
| three_default | equilateral triangle, sub-Keplerian | `data/multi_three_default.npz` |

**Figures.** `figures/three_galaxy/`
- 4-panel grid: `multi_three_default_4panel.png`
- Gif: `multi_three_default.gif`

**Results (outline -- flesh out tomorrow).**
- [ ] Initial geometry: triangle inscribed in r ~ 30 code units
- [ ] First interaction: two of the three meet first (chaos-dependent), form a pair while the third orbits
- [ ] Eventual outcome: full three-way merger into a single hotter, more disturbed remnant than any 1:1 merger
- [ ] Phase-mixing is more complete (more shuffling per relaxation time)
- [ ] **Punchline**: real compact galaxy groups (Stephan's Quintet, Seyfert's Sextet) undergo this kind of multi-body merging; the outcome is generally more violent than serial pairwise mergers

**Presentation hook.** "What if it's not just two?" Demonstrate the
qualitative jump in complexity going from 2-body to 3-body interactions.

---

## Extension 5: Collision axis (head-on vs in-plane)

**Motivation.** The standard setup has the orbit plane = disk plane
(x-axis collision -- disks approach edge-on). A z-axis collision puts
the orbit perpendicular to the disk planes, so the disks approach
face-on (one through the other's spin axis). This unlocks the
**Cartwheel ring** mechanism (Lynds & Toomre 1976) when the impact is
head-on.

**Setup.** Same as baseline, but `--collision-axis z` swaps the role
of position and velocity components so the orbital plane is x-z
instead of x-y. Disks themselves remain in x-y.

**Grid.**

| label | pericentre | inclination | physical interpretation | file |
|---|---|---|---|---|
| zaxis_i0      | 5 | 0   | face-on disks, swing-past | `data/multi_zaxis_i0.npz` |
| zaxis_i90     | 5 | 90  | one face-on / one edge-on | `data/multi_zaxis_i90.npz` |
| zaxis_i180    | 5 | 180 | anti-aligned face-on      | `data/multi_zaxis_i180.npz` |
| zaxis_headon  | **0** | 0 | true head-on smash (no tangential motion -- straight-line collision) | `data/multi_zaxis_headon.npz` |

**Figures.** `figures/z_axis/`
- 4-panel grids: `multi_zaxis_{i0,i90,i180,headon}_4panel.png`
- Individual gifs: `multi_zaxis_{i0,i90,i180,headon}.gif`
- Master comparison: `zaxis_grid.{png,gif}` (does NOT include headon by default)
- **Ring diagnostic**: `multi_zaxis_headon_ring_diagnostic.png` (t = 18 - 45),
  `multi_zaxis_headon_ring_zoom.png` (t = 21 - 40), each showing
  face-on snapshots of galaxy 1's disk + radial Sigma(R) profile with
  initial-disk reference overlay. Generated by `ring_diagnostic.py`.

**Results (outline -- flesh out tomorrow).**
- [ ] zaxis i0/i90/i180 (peri = 5): orbital swing-past with face-on disk impact; tidal response differs from in-plane case due to perpendicular geometry
- [ ] **HEAD-ON (peri = 0): Cartwheel ring forms!**
  - Compression t = 21 - 27 (disk falls inward toward impulsive potential)
  - Rebound t = 30: hollow centre + bright annulus visible
  - Peak ring t = 33: classic Cartwheel morphology, R_ring ~ 10-12 kpc
  - Dissipation t = 36+ (ring expands outward and dilutes)
  - Sigma(R) profiles confirm: profile becomes non-monotonic with a bump above the initial-disk baseline at R ~ 8-12 kpc, then moves outward and flattens
  - Ring expansion speed measured at ~50 km/s (Δ R = 5 kpc over Δ t = 6 code units = 0.08 Gyr) -- consistent with real Cartwheel-galaxy observations
- [ ] **Physics**: perpendicular impact gives every disk star an inward radial impulse; stars fall, overshoot, rebound outward in phase, forming a coherent outward density wave (Lynds & Toomre 1976; Theys & Spiegel 1977)
- [ ] **Why equal-mass differs from real Cartwheel**: real Cartwheel = small intruder through large target, so only the target rings. Our equal-mass setup rings both galaxies symmetrically.

**Presentation hook.** "Beyond the planar merger" -- demonstrate that
varying the orbit orientation unlocks a qualitatively different
morphological feature (the ring) that no in-plane geometry can produce.
This is the **flagship visual** of the project (see
`project_3_presentation_results.md` memory).

---

## Pipeline summary (for the methods slide)

**Runners** (all idempotent):
- `scan_runner.sh` (branch `extension/orbit-scans`): pericentre + inclination scans
- `multi_runner.sh` (branch `extension/multi-galaxy`): mass-ratio scan + three-galaxy
- `zaxis_runner.sh` (branch `extension/multi-galaxy`): four z-axis variants
- `extended_runner.sh`: re-runs peri10, peri20, 1to8 at tstop = 900 to clarify late-time outcomes
- `scan_gif_runner.sh`: generates per-trajectory gifs for everything above

**Visualization tools** (build from any .npz):
- `merger_snapshot_grid.py`: 4-panel face-on + edge-on snapshots
- `merger_animation.py`: per-trajectory gif (face-on)
- `scan_grid_plot.py`: master comparison png (one row per trajectory)
- `gif_grid.py`: master comparison gif (synchronised side-by-side)
- `ring_diagnostic.py`: face-on snapshots + radial Sigma(R) profile (Cartwheel diagnostic)
- `orbit_illustration.py`: parabolic orbit illustration in COM frame
  (the four-panel pericentre = 0 / 5 / 10 / 20 set)

**Code modifications made for extensions:**
- `evolve_merger.py`: added `--collision-axis x|z`
- `merger_initial_conditions.py`: added `params_2` (unequal mass) and
  `collision_axis` parameter; mass-weighted COM split for unequal mass
- `merger_analysis.py`: added `robust_stellar_mask` for N-galaxy / unequal-mass support
- `merger_animation.py`, `merger_snapshot_grid.py`: generalised from 2 to N galaxies
- New: `three_galaxy_initial_conditions.py`, `evolve_unequal_merger.py`, `evolve_three_galaxy.py`

---

## Suggested presentation framing (axes of evaluation)

Each extension answers a different "what-if" question about the merger:

1. **What pericentre is needed to merge?** (pericentre scan)
2. **How does disk orientation shape the tidal tails?** (inclination scan)
3. **When does a merger stop being a merger and become accretion?** (mass-ratio scan)
4. **What happens with three?** (three-galaxy)
5. **What's possible when the orbit isn't in the disk plane?** (z-axis -> Cartwheel ring)

Reuse the same baseline (`merger_N40k_eps01`) as the reference visual
across all five slides so the audience can see what each knob does
against a common backdrop.

**Lead with #5 (Cartwheel)** -- most visually striking, most unique
result, ties to a famous real galaxy.

---

## TODO for tomorrow

- [ ] Fill in "Results" bullets per extension with the actual late-time
      verdicts from the 4-panel grids and gifs
- [ ] Pick the 5-10 most striking figures for the slide deck
- [ ] Write 1-2 sentence captions per chosen figure
- [ ] Decide which master grids (pericentre, inclination, mass-ratio, z-axis)
      go in main slides vs backup
- [ ] Decide if the orbit_illustration set + ring_diagnostic deserve their
      own slides or appear as in-line evidence
- [ ] Final ordering: probably (1) brief recap of baseline merger (Q7),
      (2) Cartwheel ring as the headline, (3) the four parameter scans as
      "axes of evaluation"
