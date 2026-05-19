# Section 3.0.6 — Additional questions: answers to Q8, Q9, Q10

All numbers are produced by `analytic_questions.py` (run it to
reproduce). Code units have G = 1; the Milky-Way scaling is
[L] = 3.5 kpc, [M] = 5.6×10¹⁰ M_⊙, [T] = 1.305×10⁷ yr (Q3). The
integrator is leapfrog with a fixed step Δt = 0.125 code units
(= 1.63 Myr); the merger softening is ε = 0.1 code units (= 0.35 kpc).

---

## Q8 — Total mass, half-mass radius and dynamical time of the disk

> **Q8.** The surface density of the initial disk is
> Σ = Σ₀ exp(−r/r_d), with r_d the exponential scale length. What is
> the total mass of the disk? For disk systems the dynamical time is
> the rotation period at the half-mass radius of the disk — compute
> the half-mass radius and the dynamical time for the initial disk
> component. What fraction of the dynamical time is the integration
> time-step of the leapfrog integrator we have used?

**Total mass.** Integrating the exponential surface density over the
plane,

  M_d = ∫₀^∞ Σ₀ e^(−R/r_d) · 2πR dR = **2π Σ₀ r_d²**.

In the model the disk is normalised to M_d = 1 code unit = **5.6×10¹⁰
M_⊙**, so the central surface density is
Σ₀ = M_d /(2π r_d²) = 1/(2π) = 0.159 code units = **728 M_⊙ pc⁻²**
(r_d = 1 code = 3.5 kpc).

**Half-mass radius.** The in-plane enclosed mass of an exponential disk
is M(<R) = M_d [1 − (1 + R/r_d) e^(−R/r_d)]. Setting this to M_d/2 gives
the transcendental equation

  (1 + x) e^(−x) = ½,  x = R_half / r_d,

solved numerically: x = 1.678, so **R_half = 1.678 r_d = 5.87 kpc**.

**Dynamical time.** The manual defines the disk dynamical time as the
*rotation period at the half-mass radius*,
t_dyn = 2π R_half / v_circ(R_half) = 2π √(R_half³ / (G M(<R_half))).
The disk orbits in the **total** potential, so the enclosed mass at
R_half is the sum of all three components:

| component | M(<R_half) [code] |
|---|---|
| disk (½ M_d by definition) | 0.500 |
| Hernquist bulge | 0.297 |
| truncated-isothermal halo | 0.495 |
| **total** | **1.292** |

⇒ **t_dyn = 12.0 code units = 157 Myr.**

**Time-step fraction.** Δt / t_dyn = 0.125 / 12.0 = **1.0×10⁻²**, i.e.
about **96 leapfrog steps per orbit** at the disk half-mass radius — the
disk's bulk dynamics are very well time-resolved.

---

## Q9 — Equation of motion in a homogeneous sphere; remnant t_dyn

> **Q9.** For a homogeneous sphere of constant density ρ the enclosed
> mass is M(r) = (4/3)πr³ρ. What is the orbital period of a mass on a
> circular orbit? What is the equation of motion for a test particle
> released from rest at radius r in its gravitational field? How long
> does the particle need to reach r = 0? (This time-scale defines the
> dynamical time of a mostly-spherical system of mean density ρ.)
> Estimate the dynamical time (in physical units) of the total merger
> remnant at radii 0.5, 1, 3 and 5 kpc. How does it compare to the
> fixed integration time-step?

For a uniform sphere of density ρ the enclosed mass is
M(<r) = (4/3)π r³ ρ, so a test particle **released from rest** at radius
r feels

  r̈ = − G M(<r) / r² = − (4πGρ/3) r.

This is **simple harmonic motion** with angular frequency
ω² = 4πGρ/3 (independent of the starting radius — every shell reaches
the centre simultaneously). Starting from rest, the particle arrives at
r = 0 after a **quarter period**:

  t_dyn = (1/4)·(2π/ω) = (π/2)/ω = **√( 3π / (16 G ρ) )**
        = (π/2) √( r³ / (G M(<r)) ).

This is the standard dynamical (free-fall) time for a mostly-spherical
system of mean density ρ = M(<r) / [(4/3)π r³].

**Merger remnant** (total system, all components, from
`data/merger_N80k_eps01.npz`, centred with the shrinking-sphere method;
Δt = 1.63 Myr):

| r [kpc] | M(<r) [M_⊙] | t_dyn [Myr] | Δt / t_dyn |
|---|---|---|---|
| 0.5 | 4.1×10⁹  | 4.1  | 0.40 |
| 1.0 | 1.3×10¹⁰ | 6.5  | 0.25 |
| 3.0 | 5.3×10¹⁰ | 16.6 | 0.098 |
| 5.0 | 9.2×10¹⁰ | 27.4 | 0.060 |

The remnant is dense and so dynamically *fast* in the centre: at
r = 5 kpc the step is a comfortable 6 % of t_dyn (~17 steps per
free-fall time), but by r = 0.5 kpc Δt has grown to **40 %** of the
local dynamical time. The innermost ~1 kpc is therefore only marginally
time-resolved — which is acceptable here precisely because that region
is also smoothed by the ε = 0.35 kpc softening, so no real structure
exists there to integrate.

---

## Q10 — Bulge half-mass radius, dynamical time, resolvability

> **Q10.** Compute the half-mass radius of the bulge component of the
> disk model. What is the dynamical time of the bulge at this radius?
> Can we resolve bulge dynamics with the simulations performed here?

The bulge is a Hernquist sphere, M(<r) = M_b r²/(r+a)², with
a = 0.1 code = 0.35 kpc. Setting r²/(r+a)² = ½ gives
**r_h = (1+√2) a = 0.241 code = 0.845 kpc**.

Using the same mean-density dynamical time as Q9 with the bulge's own
enclosed mass (½ M_b):

  **t_dyn(r_h) = 0.456 code units = 5.95 Myr.**

(This is an *upper bound*: the halo and disk add mass inside r_h, which
would shorten the true orbital time and only sharpen the conclusion
below.)

**Can we resolve bulge dynamics? No, only marginally.** Two independent
limits both bite at the bulge scale:

1. **Spatial:** the Plummer softening ε = 0.1 code equals the bulge
   scale length a *exactly* (ε / a = 1.00). The bulge's defining
   structure — its central cusp at r ≲ a — is therefore washed out by
   the force softening; only the half-mass radius and outward is
   spatially meaningful.
2. **Temporal:** Δt / t_dyn(r_h) = 0.125 / 0.456 = **0.27**, i.e. only
   **~4 leapfrog steps per orbit** at the bulge half-mass radius (and
   far fewer deeper in). Compare the disk's ~96 steps/orbit (Q8): the
   bulge sits right at the edge of what the fixed step can follow.

So the global disk and the merger remnant on kpc scales are well
resolved, but the **bulge's internal structure is not** — it is
softened on its own scale radius and integrated with only a handful of
steps per dynamical time. Sub-kpc bulge dynamics from these runs should
be treated as qualitative, not quantitative.

---

## Q11

> **Q11.** Request to report errors, typos, or difficulties, and to
> suggest improvements to the exercise (personally or to
> naab@usm.lmu.de).

Q11 is the manual's request for feedback to the author
(naab@usm.lmu.de) — no scientific deliverable.
