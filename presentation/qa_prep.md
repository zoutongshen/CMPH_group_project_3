# Q&A prep — Project 3 (galaxy–galaxy mergers)

Anticipated questions from classmates and graders, grouped. Each answer is
kept to ≤3 sentences for live delivery; the `▸` line points at the slide or
backup that supports it. All numbers are from the code and `EXTENSIONS.md`,
kept in sync with the deck. The audience has the source, so it's fine to say
"it's in `evolve_merger.py`" — but lead with the physics.

---

## Fast cheat-sheet — the five axes

| axis | knob (baseline → range) | result | physics |
|---|---|---|---|
| Pericentre | `r_p` 5 → 1, 5, 10, 20 (3.5–70 kpc) | merge t≈60 / 120 / 600; r_p=20 escapes (no return by t=900) | drag ∝ ρ — wider pass barely interpenetrates → low density → little energy lost |
| Inclination | disk-2 tilt 30° → 0–180° | prograde = long symmetric tails; retrograde ≈ none | spin–orbit resonance lifts only co-rotating material |
| Mass ratio | 1:1 → 1:2, 1:4, 1:8 | 1:1 destroys both → elliptical; 1:8 secondary absorbed by t≈800 | merger→accretion is one continuum |
| Number | 2 → 3 galaxies | near-simultaneous collapse → hotter, more mixed remnant | formally chaotic, but this symmetric setup's order is seed-robust |
| Collision axis | x → z, head-on (r_p=0) | Cartwheel ring, peak t≈33, R≈10–12 kpc, ~50 km/s | inward impulse → overshoot → outward Σ(R) wave |

**If you remember four numbers:** ε = 0.1 = 0.35 kpc · v_c ≈ 240 km/s
(Milky-Way-like) · ring R ≈ 10–12 kpc expanding ~50 km/s · merge/flyby
threshold between r_p = 10 and 20.

---

## A. Method & numerics

**Q: Why direct N² instead of a tree code (Barnes–Hut)?**
At our N (≤10⁵) the JIT-parallel direct sum is competitive — 0.06 s vs 0.04 s
per force eval at N = 20k — and it's *exact* and conserves momentum to machine
precision because every pair force is symmetric. The tree's cell approximation
breaks that pair symmetry, so momentum drifts at the error level. We
cross-validated both on the same galaxy: 0.4% median force error, identical
profiles to t = 100; the tree only wins asymptotically (1.0 s → 0.15 s at 80k).
▸ slide 9 "Same gravity, two algorithms" + backup "Tree vs direct N²".

**Q: How did you pick the softening ε = 0.1?**
ε caps the force at close approach so a single near-collision can't dominate the
timestep; 0.1 code = 0.35 kpc, which is well below the disk scale length
(h = 3.5 kpc), so it suppresses spurious two-body scattering without smearing
real structure. It's the standard accuracy-vs-stability trade, and we checked
the residual disk heating scales with N and ε as theory predicts.
▸ slide 8 validation + backup "finite-N heating".

**Q: How well is energy conserved?**
We integrate with symplectic velocity-Verlet (kick–drift–kick leapfrog), so
energy is *bounded*, not secularly drifting. In the Kepler two-body test over
ten orbits the energy error is ~10⁻¹⁰ and angular momentum holds at machine
precision — that's what justifies the long 900-time-unit runs.
▸ slide 4 (integrator) + slide 8 (Kepler bullet).

**Q: Isn't the disk thickening a sign of a broken initial condition?**
No — we traced it to finite-N two-body relaxation, not an IC error: it's set by
N and ε, and quadrupling N cuts the z₀ drift from 67% to 24%, the predicted
factor. The radial structure stays stable (scale length within ~6%) and the
rotation curve stays flat at ~240 km/s, so the equilibrium is sound.
▸ slide 8 + backup "finite-N heating: N=20k vs 80k".

**Q: What are the code units, and how do you get physical numbers?**
We work in G = 1 units with disk mass = 1 and disk scale length = 1; fixing
h = 3.5 kpc and M_d = 5.6×10¹⁰ M_⊙ then sets the time unit to 13 Myr and the
velocity unit to 262 km/s. Every kpc / Gyr / km·s⁻¹ I quoted is just that
scaling applied.
▸ slide 5 caption.

**Q: Why a parabolic encounter orbit?**
It's the standard idealization of a first encounter falling together from a
cosmological turnaround — marginally bound, zero orbital energy — so the outcome
depends on the galaxies' structure and the friction, not on an arbitrary chosen
binding energy. It's also Toomre & Toomre's original choice, which keeps the
tidal-tail comparison clean.
▸ slide 13 baseline table.

**Q: On the pericentre slide — why does a *wider* pass lose less energy? Isn't a wide pass slower?**
The driver is interpenetration, not speed: a small pericentre drives the galaxies
deep into each other's high-density interior, where dynamical friction (drag ∝ ρ)
is strong, so lots of orbital energy is drained and they bind; a wide pericentre
only grazes the thin outskirts, so almost none is lost. Speed isn't an independent
knob here — the orbit is parabolic (E = 0), so a *closer* pass is actually faster
at pericentre; the `1/v²` term would weaken drag, but the density term dominates,
giving a net energy loss that scales ~1/r_p. ("Faster = less efficient" is only the
right intuition when you compare different orbital energies at *fixed* pericentre.)
▸ slide 14 (pericentre).

---

## B. Physics & interpretation

**Q: Why exactly do retrograde encounters not raise tails?**
A tail is a resonance: near-side stars that co-rotate with the orbit feel the
companion's pull in the same direction for an extended time and get flung out
into a long tail. In a retrograde disk the stars sweep past the perturbation the
opposite way, so the force reverses before it can do coherent work — no tail.
▸ slide 15.

**Q: How do you decide a "merger" actually happened vs a flyby?**
We track the separation of the two density centres (and the bound mass): a
merger is when the cores sink together and stay within a few softening lengths,
a flyby is when the separation keeps growing past the box. For r_p = 20 we
re-ran out to t = 900 to confirm it never turns around.
▸ slide 14 + backup "Pericentre — late-time flyby".

**Q: Is the three-galaxy result just one realization — would another seed differ?**
We checked: across four random seeds, galaxy 1 pairs with galaxy 2 first in all of
them, so the merger order is robust — it's set by the deterministic inclination
geometry (galaxy 1 flat, 2 at +30°, 3 at −30°), not by the random realization. The
three-body problem is chaotic in principle, but this symmetric setup is too
constrained to amplify the seed noise into a different order; every seed gives a
near-simultaneous three-way collapse to the same hot remnant. Seeing true
seed-sensitivity would need a deliberately asymmetric configuration.
▸ slide 17.

**Q: What does the Toomre Q = 1.5 actually buy you?**
Q is the disk's stability margin against its own self-gravity; Q ≳ 1 stops the
disk fragmenting into clumps on its own, so any structure that appears in a
merger is from the interaction, not a spontaneous instability. We set Q = 1.5 at
R = 2.4h.
▸ slide 7 (motion).

**Q: How robust is calling the remnant an elliptical?**
The signature is the straight line in the R¹ᐟ⁴ (de Vaucouleurs) plot — the
defining profile of an elliptical — and the fit gives an effective radius
R_e ≈ 4.8 kpc. As a model-independent cross-check the half-mass radius (5.6 kpc)
agrees, so the classification doesn't hinge on the fit itself.
▸ slide 12 (remnant).

---

## C. Cartwheel specifics

**Q: Is the ring a numerical artifact?**
No. The surface-density bump *moves outward* in radius over time — a static
artifact wouldn't propagate — and it's stable to our timestep and softening.
Independently, its expansion speed of ~50 km/s matches the observed Cartwheel.
▸ slide 19 + backup "Cartwheel — full ring diagnostic".

**Q: Why does your ring look different from the real Cartwheel?**
Because our collision is equal-mass, so both disks ring symmetrically. The real
Cartwheel is a small intruder punching through a much larger target, so only the
big target forms the ring while the small intruder flies on.
▸ slide 19 caveat.

**Q: Where does the Δv ≈ 2GM_p/(bV) impulse come from?**
It's the impulse approximation: integrate the intruder's transverse
gravitational acceleration along a straight-line pass, and the velocity kick
comes out proportional to the intruder mass over (impact parameter × speed). A
fast, central pass gives the cleanest, most coherent kick — which is exactly the
head-on, pericentre-0 case.
▸ slide 19.

**Q: How did you measure the ~50 km/s expansion?**
From the ring diagnostic: the Σ(R) bump moves ΔR = 5 kpc over Δt = 6 code units
(0.08 Gyr), giving ~50 km/s once converted to physical units.
▸ backup "Cartwheel — full ring diagnostic" (t = 18–45).

---

## D. Limits & next steps

**Q: What's the biggest limitation?**
It's collisionless — pure gravity, no gas — so no shocks, no star formation, no
dissipation. Real gas-rich ("wet") mergers funnel gas to the centre and can
trigger starbursts or rebuild a disk, which we can't capture; adding an SPH or
grid gas component is the obvious next step.
▸ slide 20 (limits).

**Q: Would using a tree code change any of your conclusions?**
No — we showed the tree and the direct sum agree to 0.4% on the same galaxy with
overlapping profiles. A tree would let us push to higher N (less disk heating,
sharper tails), but the physics of all five axes is already converged at the N
we used.
▸ slide 9 + backup "Tree vs direct N²".

**Q: Could you reproduce a *specific* observed system (the actual Antennae / Cartwheel)?**
In principle yes — that's a targeted-orbit-fitting problem: tune pericentre,
inclination, mass ratio and viewing angle to match one system, as Toomre &
Toomre did by hand. We did the complementary thing — map how each knob changes
the outcome — rather than fit one object.
▸ slide 20.

---

### Backup slides available (press ↓ from the talk)
- Cartwheel — full ring diagnostic (Σ(R), t = 18–45)
- Finite-N heating: N = 20k vs 80k
- Pericentre — orbit geometry & late-time flyby
- Inclination — per-run panels
- Tree vs direct N² — force accuracy & cost
- Mass ratio — 1:8 late-time absorption
