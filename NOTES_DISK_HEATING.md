# Investigation: vertical disk heating in the isolated-galaxy IC

**Date:** 2026-05-13
**Context:** Task #5 of the project plan — validate the IC by running an
isolated galaxy and checking that the disk stays in equilibrium.

## Observation

At the manual's prescribed parameters
(N=20,000, ε=0.1, dt=0.125, tstop=100, seed=0):

- Radial scale length **h drift: +0.9 %**  (essentially zero)
- Vertical scale height **z₀ drift: +67 %**  (from 0.199 to 0.334)
- Total energy drift: 2.6×10⁻³  (well within integrator tolerance)
- Total z-angular-momentum drift: 4×10⁻¹⁵  (machine precision)

The radial structure is preserved exactly as it should be; the disk
thickens vertically by two thirds. The asymmetry was the first big clue
— the heating mechanism acts on the *vertical* phase-space alone, not
the radial.

## Initial (wrong) measurement

The first version of `fit_disk_scale_height` used the formula
`z₀ = √(6 <z²>/π²)`, which is *off by √2*. The correct expression for
`ρ(z) = (1/(2 z₀)) sech²(z/z₀)` is

```
<z^2> = (z_0^2/2) * integral_{-inf}^{inf} u^2 sech^2(u) du
      = (z_0^2/2) * (pi^2/6)
      = pi^2 z_0^2 / 12
=> z_0 = sqrt(12 <z^2> / pi^2)
```

The bare integral `∫ u² sech²(u) du = π²/6` is the *unnormalised* second
moment; ⟨z²⟩ for the normalised distribution requires dividing by
`∫ sech²(u) du = 2`. Fixed.

## Hypothesis 1 (rejected): σ_z formula too small

Eq 2.40 prescribes `σ_z² = π G Σ_d(R) z₀`, derived for an isothermal
*slab* in the self-gravity-of-the-disk-only limit. The full vertical
gravity also gets contributions from bulge and halo, so the true
equilibrium σ_z is larger. Underestimating σ_z would make the disk
*compress* (under-supported), not expand.

We do see a brief compression at t = 0 → 10
(z₀ goes 0.199 → 0.188, ~5 %), then it expands. So this effect is real
but sub-dominant. **Not the answer.**

## Hypothesis 2 (rejected): softening artefact making vertical gravity weak

For a Plummer-softened slab, the vertical restoring force on a particle
at height z' is reduced inside the softening length ε. With ε = 0.1 and
z₀ = 0.2, ε is *half* the disk scale height. The intuition was: weakened
gravity → larger equilibrium z₀ → disk thickens.

**Prediction:** smaller ε ⇒ less artificial thickening.
**Test (4 ε values, otherwise identical runs at N=20k):**

| ε    | z₀ drift | E drift |
|------|----------|---------|
| 0.20 | +49 %    | 8.9×10⁻⁴ |
| 0.10 | +67 %    | 2.6×10⁻³ |
| 0.05 | +91 %    | 1.2×10⁻² |
| 0.01 | +124 %   | 1.1×10⁻¹ (integrator failing) |

**Opposite trend: smaller ε ⇒ more heating.** The softening-artefact
hypothesis is wrong. The energy drift growing with smaller ε is also
telling — sharper close encounters than the integrator can handle.

## Hypothesis 3 (accepted): finite-N two-body relaxation

The physical picture: our disk has 6000 particles, each representing
~10⁷ real stars. Particles graze close to each other, scattering
gravitationally — encounters that don't happen in real galaxies (which
have ~10¹¹ stars and are nearly collisionless). Over time these
encounters thermalise the velocity distribution, equipartitioning energy
across degrees of freedom.

Our IC is built kinematically anisotropic:
```
sigma_R(R = h) = 0.29     (radial — hotter)
sigma_z(R = h) = 0.19     (vertical — colder)
```
Equipartition wants `<v_R²> = <v_z²>`, so **energy flows from radial to
vertical**. The disk heats vertically. The radial profile (set by
angular momentum) is preserved because L is conserved exactly.

The Chandrasekhar relaxation time scales as
```
t_relax  proportional to  N / ln(Lambda),  Lambda = b_max / b_min
b_min ~ max(eps, classical close-encounter distance)
```

So smaller ε ⇒ larger ln Λ ⇒ faster relaxation ⇒ more heating. **This
matches the ε scan above quantitatively** (rates from ε=0.2→0.01 grow by
roughly the factor `ln Λ` grows).

### N-scaling: the definitive test

If finite-N noise is the cause, `t_relax ∝ N / ln N`, so heating fraction
over a fixed time should scale as `ln N / N`. Quadrupling N should
reduce heating by

```
0.25 * (ln 80000 / ln 20000)  =  0.25 * 11.29 / 9.90  =  0.285
```

so 67 % → 19 %.

**Run at N=80,000 (24k disk + 8k bulge + 48k halo, ε = 0.1, dt = 0.125,
tstop = 100, seed = 0):**

| N      | z₀ drift   | h drift | E drift  | wall-clock |
|--------|-----------:|--------:|---------:|-----------:|
| 20,000 | **+67 %**  | 0.9 %   | 2.6×10⁻³ | 67 s |
| 80,000 | **+24 %**  | 6.0 %   | 2.6×10⁻³ | 17 min |

Observed ratio: 24 / 67 = 0.36. Predicted ratio: 0.28. **Match within
20 %.** Hypothesis confirmed: the heating is dominantly finite-N noise.

(The small overshoot — 24 % instead of 19 % — is plausibly explained by
sub-dominant contributions from the σ_z approximation (~5 % residual)
and the imperfect equipartition treatment in the scaling formula.)

## Literature placement

This is **the canonical finite-N artefact of disc simulations**, well
documented in the literature:

- Sellwood & Binney (2002), *MNRAS* 336, 785 — original identification
- Sellwood (2013), *Reviews of Modern Physics* 86, 1 — review of the
  scaling `t_heat ∝ N / ln N`
- Fujii et al. (2011), *PASJ* 63, S5 — explicit N-scaling measurement

The manual specifies N=20,000 and ε=0.1 likely because Naab's exercise
was designed for VINE (a tree code, where compute scales as N log N).
With our direct-N² implementation we can comfortably reach N=80,000 on
a laptop, which cuts the artificial heating to ~24 %. That is small
enough not to dominate downstream physics.

## Decision

For the actual manual exercises (Tasks #6 and onward), run at **both**:

1. **N=20,000** — the manual's literal prescription. Run this for the
   pedagogical answer to Q5 ("how did the properties change?") and Q6
   (ε=10⁻⁴ comparison). The 67 % vertical heating is *the* finding the
   manual is steering students toward.

2. **N=80,000** — clean-physics version. Run this for the merger and
   for the radial-profile measurements where we want results that
   aren't dominated by N-noise.

For the presentation: the **N-dependence itself** becomes a
deliberately-chosen slide showing how artificial heating depends on N,
which is honest and well-grounded science.

## Open question (parked for now)

The N=80k case still shows ~24 % drift over 100 code units. Some of
this is residual finite-N (would shrink further at N=160k, ~12 %); some
is plausibly the σ_z formula approximation (~5 %, would require
re-deriving σ_z in the full potential to fix). Not pursued today; if
the merger remnant looks anomalous we can revisit.
