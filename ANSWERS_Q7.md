# Section 3.0.5 — Galaxy merger: answer to Q7

> **Q7.** Discuss the changes in cumulative mass and surface density
> profiles (of the merger remnant relative to the progenitor disks).

## Setup

Two copies of the validated equilibrium galaxy (`merger_initial_conditions.py`),
same spin sense, the second disk inclined 30 deg, placed on a parabolic
two-body encounter (initial separation 30, pericentre 5 code units),
eps = 0.1, dt = 0.125, integrated to t = 300 (~3.9 Gyr).

Runs analysed:

| label | N total (per galaxy) | file |
|---|---|---|
| test | 12 000 (6k stellar + halo) | `data/merger_test12k_rp5.npz` |
| production | 40 000 (20k/galaxy) | `data/merger_N40k_eps01.npz` |
| clean cross-check | 80 000 (40k/galaxy) | `data/merger_N80k_eps01.npz` |
| high-N fidelity | 160 000 (80k/galaxy) | `data/merger_N160k_eps01.npz` |

Figures: `figures/merger_*_remnant.png` (profiles, Q7),
`figures/merger_N*.gif` (animation). Reproduce with
`python evolve_merger.py …` then `python merger_analysis.py …`.

## What happens

The galaxies approach, are tidally distorted at first pericentre,
sink by dynamical friction on their dark haloes and coalesce into a
single concentrated body by t ~ 120 (~1.6 Gyr); by t ~ 180 it is a
relaxed, round, pressure-supported spheroid in which the two progenitor
star populations are completely phase-mixed. Tidal tails form but are
modest — expected, because these galaxies are halo-dominated
(M_halo = 5.8 vs M_disk = 1), which makes dynamical friction strong and
suppresses long tails relative to disk-dominated encounters
(Barnes 1988).

## Surface density (the headline result)

- **Progenitor:** an exponential disk — a straight line in
  log Sigma vs R (slope = −1/h), curved in the R^{1/4} projection.
- **Remnant:** **no longer exponential.** It is a straight line in
  log Sigma vs R^{1/4} over ~4–5 decades in Sigma — i.e. it obeys the
  **de Vaucouleurs R^{1/4} law**, the defining photometric signature of
  an elliptical galaxy. Measured fit slope:

  | run | de Vaucouleurs slope | remnant R_half [kpc] |
  |---|---|---|
  | 12 000 | −2.18 | 6.67 |
  | 40 000 | −2.22 | 5.71 |
  | 80 000 | −2.26 | 5.56 |
  | 160 000 | −2.25 | 5.58 |

  The slope is stable across a factor ~13 in N, so this is a physical
  result, not a finite-N artefact: **merging two cold exponential disks
  produces a hot R^{1/4} spheroid** (the Toomre 1977 "mergers make
  ellipticals" hypothesis, reproduced from scratch).

## Cumulative mass

- The stellar remnant contains ~2x the stellar mass of one progenitor
  (both galaxies' disks + bulges, conserved).
- M_star(<r) is **more centrally concentrated** in the inner few kpc
  (the violently-relaxed core is denser than either original disk
  centre) **and simultaneously more extended** at large r: a diffuse,
  loosely-bound stellar envelope of tidal debris is deposited out to
  tens of kpc. Mass is redistributed *both* inward and outward; the
  intermediate exponential disk structure is destroyed.
- The half-mass radius of the stellar remnant is ~5–7 kpc (run table
  above), comparable to a real intermediate-mass elliptical.

## Why: the physical mechanism

During the merger the gravitational potential varies violently on a
dynamical timescale. **Violent relaxation** (Lynden-Bell 1967)
scrambles particle energies independently of mass, erasing the disks'
cold ordered rotation and driving the system toward the
pressure-supported, centrally-cusped, R^{1/4} equilibrium that
characterises ellipticals. Phase mixing then smooths the remnant and
intermingles the two progenitor populations. Angular momentum that was
in ordered disk rotation is partly transported outward into the tidal
debris and the dark haloes.

## References

Toomre & Toomre (1972); Toomre (1977) — mergers form ellipticals;
Lynden-Bell (1967) — violent relaxation; Barnes (1988, 1992),
Hernquist (1992, 1993) — N-body disc-galaxy mergers; de Vaucouleurs
(1948) — the R^{1/4} law. All also cited in the manual's bibliography.

## Convergence / N-dependence

Unlike the isolated-disk vertical heating (which is finite-N dominated,
see `NOTES_DISK_HEATING.md`), the merger remnant's gross structure is
set by violent relaxation and is *not* finite-N sensitive: the
R^{1/4} slope moves −0.18 → −0.22 → −0.26 → −0.25 from N = 12k
to 40k to 80k to 160k, and the remnant half-mass radius settles at
6.7 → 5.7 → 5.6 → 5.6 kpc. Doubling the resolution from 80k to
160k shifts the slope by only +0.01 and the half-mass radius by
0.02 kpc (~0.4 %) — i.e. the structural conclusion was already
converged at 80k, and the 160k run (executed overnight, 4.3 h wall
on 10 cores) confirms it directly rather than merely extrapolating.
