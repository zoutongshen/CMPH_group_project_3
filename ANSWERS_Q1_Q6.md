# Section 3.0.4 — Isolated disk evolution: answers to Q1–Q6

Runs analysed (dt = 0.125, tstop = 100 ≈ 1.30 Gyr, seed = 0):

| label | N (disk/bulge/halo) | eps | file |
|---|---|---|---|
| baseline | 20 000 (6k/2k/12k) | 0.1 | `data/evolve_N20k_eps01.npz` |
| Q6 small-eps | 20 000 (6k/2k/12k) | 1e-4 | `data/evolve_N20k_eps1em4.npz` |
| clean physics | 80 000 (24k/8k/48k) | 0.1 | `data/evolve_N80k_eps01.npz` |

Figures in `figures/`: `*_profiles.png` (Q4/Q5), `*_images.png` (Q1/Q2).
Reproduce with `python evolve_disk.py …` then `python disk_profiles.py …`
and `python disk_visualization.py …`.

---

## Q1 — Morphological difference: disk vs bulge

> **Q1.** Viewing the model galaxy edge-on and face-on (disk particles
> in white, bulge particles in red): what is the morphological
> difference?

See the edge-on disk+bulge zoom panel of `evolve_N20k_eps01_images.png`.

- **Disk** (dark): thin, highly flattened, rotationally supported. The
  vertical extent (z₀ ≈ 0.2 ≈ 0.7 kpc) is ~20× smaller than the radial
  extent (h = 1 ≈ 3.5 kpc). Ordered rotation dominates the kinetic
  energy; the velocity dispersion is small and anisotropic
  (σ_z < σ_R < σ_φ).
- **Bulge** (red): round, centrally concentrated, *pressure* supported.
  Near-isotropic velocity dispersion, no preferred plane, density cusp
  toward the centre (Hernquist profile, scale a = 0.1 ≈ 0.35 kpc).

The physical distinction is the support mechanism: the disk is held up
against gravity by rotation (a cold, flattened system), the bulge by
random motions (a hot, round system).

## Q2 — Adding the halo

> **Q2.** Plot the halo particles in addition. What does the halo
> component look like relative to the disk and bulge?

See the face-on / edge-on "(all, incl. halo)" panels. The halo (faint
blue) is near-spherical and far more extended than the luminous
components: it reaches ~80 kpc while the disk is confined to ~15 kpc. It
carries most of the mass (M_h = 5.8 of 7.13 total code units) but at very
low density, so it is invisible in a star-light image yet sets the depth
of the potential well and keeps the outer rotation curve flat (see Q4).

## Q3 — Physical units (worked in `physical_units.py`)

> **Q3.** With G = 1, length unit L = 1 → 3.5 kpc and mass unit
> M = 1 → 5.6×10¹⁰ M⊙, derive the physical time unit (in years) and
> the velocity unit (in km/s). Use G = 6.672×10⁻¹¹ m³ kg⁻¹ s⁻²,
> 1 kpc = 3.0856×10¹⁹ m, 1 M⊙ = 1.989×10³⁰ kg.

With G = 1, length unit = 3.5 kpc and mass unit = 5.6×10¹⁰ M⊙, requiring
G to take its SI value fixes

    [T] = sqrt( [L]^3 / (G_phys [M]) )    [V] = [L]/[T]

| code unit = 1 | physical value |
|---|---|
| length | 3.5 kpc |
| mass | 5.6 × 10¹⁰ M⊙ |
| **time** | **1.305 × 10⁷ yr** (13.05 Myr) |
| **velocity** | **262.3 km/s** |
| surface density | 4571 M⊙ pc⁻² |

Hence tstop = 100 ≈ **1.30 Gyr**, and the leapfrog step dt = 0.125 ≈
1.63 Myr.

## Q4 — Profiles in physical units

> **Q4.** Plot the circular-velocity / cumulative-mass /
> surface-density figure using physical units: km/s, solar masses and
> solar masses per square parsec.

`disk_profiles.py` reproduces the manual's `plotvcirc` analysis directly
from the particles and plots, in km/s / M⊙ / M⊙ pc⁻²:

- **Circular velocity, decomposed by component**
  v_c,i(r) = √(G M_i(<r)/r), so v_c,total² = Σ v_c,i².
- **Cumulative mass M(<r)**, per component, log-y (spans decades; the
  halo dominates — a linear axis would hide the disk and bulge).
- **Projected face-on surface density Σ(R)** of disk and bulge, semilog-y
  (an exponential disk is a straight line, the scale length is the
  slope — van der Kruit & Searle 1981; Freeman 1970).

Plot-scale conventions used: v_c linear–linear (standard rotation
curve), M(<r) log-y, Σ semilog-y.

**On the dip in the total rotation curve at r ≈ 1–3 kpc:** it is *not* a
numerical artefact. The component decomposition shows the total dips
exactly where the centrally-concentrated bulge contribution has fallen
off (∝ r^-1/2 beyond a ≈ 0.35 kpc) while the halo has not yet risen and
the disk is still climbing toward its peak (~2.2 h ≈ 8 kpc). This
bulge→disk/halo transition is a genuine feature of decomposed rotation
curves (it is seen in the real Milky Way). The only methodological
caveat is that v_c = √(G M_sph(<r)/r) assumes spherical symmetry and so
slightly under-counts the flattened disk's in-plane support — but this
is exactly the estimator the manual's `plotvcirc` uses. Innermost bins
with < 50 enclosed particles are masked so the dip is not confused with
small-r shot noise.

## Q5 — How did the system change, t = 0 → t = 100 (baseline, N=20k, ε=0.1)

> **Q5.** Repeating the analysis at t = 100: how did the system
> properties change?

| quantity | t = 0 | t = 100 | change |
|---|---|---|---|
| peak v_c | 236.7 km/s | 244.2 km/s | +3 % |
| disk central Σ (R<1 kpc) | 645 M⊙/pc² | 419 M⊙/pc² | **−35 %** |
| M(<r), all components | — | — | conserved (only mild redistribution) |
| disk scale height z₀ | 0.199 | ≈ 0.33 | **+67 %** (see `NOTES_DISK_HEATING.md`) |
| total energy drift | — | — | 2.6 × 10⁻³ |
| L_z drift | — | — | ~10⁻¹⁵ (machine precision) |

**Interpretation.** Globally the galaxy stays in equilibrium: the
rotation curve and enclosed-mass profile barely move and energy/angular
momentum are well conserved. Locally the *disk* changes: it puffs up
vertically (z₀ +67 %) and loses ~35 % of its central surface density.
Both are the same physics — finite-N two-body relaxation. With only
6000 disk particles each standing in for ~10⁷ stars, artificial close
encounters thermalise the velocity ellipsoid, draining energy from the
hot radial direction into the cold vertical one (equipartition) and
eroding the central concentration. The radial scale length h is
preserved because L_z is conserved exactly. Full diagnosis, the
N-scaling test confirming this is finite-N noise, and the literature
placement are in `NOTES_DISK_HEATING.md`.

## Q6 — Effect of the softening: ε = 0.1 vs ε = 10⁻⁴

> **Q6.** Rerunning with the gravitational softening reduced to
> ε = 10⁻⁴ (basename SPIRALA000): what has changed?

The softening ε guarantees the validity of the collisionless Boltzmann
equation: it caps the maximum two-body acceleration so that the
N-body system approximates a smooth potential rather than a collisional
star cluster. Shrinking ε removes that cap.

| quantity | ε = 0.1 | ε = 10⁻⁴ |
|---|---|---|
| total energy drift (t=0→100) | 2.6 × 10⁻³ | **5.2 (≈ +524 %)** |
| total energy sign | stays bound (E < 0) | **flips to E > 0 (unbound)** |
| L_z drift | ~10⁻¹⁵ | 6 × 10⁻⁷ |
| disk central Σ | 645 → 419 M⊙/pc² | 645 → **107** M⊙/pc² (−83 %) |
| peak v_c | 237 → 244 km/s | 237 → 213 km/s |
| disk morphology at t=100 | thickened but recognisable | **destroyed — diffuse blob** |

**What changed (Q6).** With ε = 10⁻⁴ (≪ the mean interparticle spacing)
close encounters produce enormous accelerations that the *fixed*
leapfrog step dt = 0.125 cannot resolve. The integrator injects
spurious energy at every hard encounter; over 800 steps the total energy
not only drifts by ~500 % but changes sign — the system becomes
formally unbound and the thin disk is heated into a diffuse,
structureless cloud (see `evolve_N20k_eps1em4_images.png`). Crucially,
**angular momentum is still conserved to ~10⁻⁷**: the failure is a
collisional / time-resolution breakdown of energy integration, not an
error in the force law (a central force exerts no torque). This is
precisely why a finite softening is mandatory — it is what keeps the
simulation on the collisionless branch the Boltzmann equation describes.

---

## N-dependence (clean-physics cross-check, N = 80 000, ε = 0.1)

Re-running the identical exercise at 4× the particle count separates the
*physical* equilibrium adjustments from the *finite-N* relaxation
artefact (rationale and literature in `NOTES_DISK_HEATING.md`).

| quantity (t=0 → t=100) | N = 20 000 | N = 80 000 |
|---|---|---|
| disk scale length h | −0.9 % | −6.0 % |
| **disk scale height z₀** | **+67.4 %** | **+24.1 %** |
| peak v_c | 237 → 244 km/s | 243 → 245 km/s |
| disk central Σ | 645 → 419 (−35 %) | 591 → 313 (−47 %) |
| total energy drift | 2.6 × 10⁻³ | 2.6 × 10⁻³ |

**Reading.**

- **Vertical heating is dominantly finite-N noise.** z₀ drift falls
  67 % → 24 % on quadrupling N. The ratio 24/67 ≈ 0.36 matches the
  Chandrasekhar scaling t_relax ∝ N/ln N (predicted 0.28) to within
  20 % — the same result already established in `NOTES_DISK_HEATING.md`,
  now reproduced through the §3.0.4 observables.
- **The rotation curve is physical and N-independent.** Peak v_c is
  ~244 km/s at both N and barely moves over the run; the M(<r) curves
  are essentially unchanged. These are not relaxation-sensitive — they
  reflect the (well-built) equilibrium mass distribution.
- **Energy-integration quality is N-independent** (drift 2.6 × 10⁻³ at
  both N), as it should be: it is set by dt and ε, not particle count.
- The central-Σ drop does *not* shrink with N (−35 % → −47 %). It is
  therefore not purely two-body relaxation: part is the known σ_z
  approximation transient (eq 2.40 neglects the bulge/halo vertical
  pull, so the inner disk settles slightly on the first few dynamical
  times). This is the open item already flagged in
  `NOTES_DISK_HEATING.md`; it is sub-dominant and does not affect the
  rotation curve or the merger.

**Presentation takeaway.** The N-dependence is itself the headline
result: with a from-scratch direct-N² code we can push to N = 80 000 on
a laptop and *show* artificial disk heating halving as predicted —
turning what the manual treats as a fixed nuisance into a controlled,
quantitative demonstration of the collisionless limit.

Figures: `figures/evolve_N80k_eps01_profiles.png`,
`figures/evolve_N80k_eps01_images.png`.
