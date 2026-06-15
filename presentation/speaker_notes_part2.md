# Speaker notes — Part 2 (Zoutong: the five axes → Cartwheel → close)

Full spoken script for the slides from **"Beyond the baseline merger"** through
**"Key references"** (deck slides 13–21). Zoutong presents this half and takes
the clicker at the divider. Time cues match the deck and land the talk at 18:30.

This is the read-from sheet — the in-deck `<aside>` notes carry only the ⏱ cue
plus the **Knob / Result / Physics** bullets. Pace ≈ 140 words/min. The
`[pause]` on the Cartwheel slide is real screen time, not words to read.

Each section states explicitly **which parameter was changed** (vs the fixed
baseline) and **how to read the physical meaning** of the outcome.

---

## 13. Beyond the baseline merger · ⏱ 0:30 (cum 9:50)

**[Zoutong takes the clicker.]**
**Parameter changed:** none yet — this sets the rule for everything that follows.
**Physical interpretation:** the whole second half is one controlled experiment — fix the baseline, move exactly one knob at a time.

> Thanks, Zhaoyang. So far we've built a galaxy from scratch, checked that it's
> genuinely stable, and collided two of them into an elliptical. Everything from
> here is one controlled experiment. We fix this entire baseline — forty thousand
> particles, this softening, this timestep, this parabolic orbit — and we turn
> exactly one knob at a time. Five knobs, five questions. Let's start with the
> simplest one: how close do they actually have to pass?

---

## 14. Axis 1 — Pericentre: how close must they pass to merge? · ⏱ 1:10 (cum 11:00)

**Parameter changed:** pericentre distance `r_p` — baseline 5 (17.5 kpc), swept over 1, 5, 10, 20 code units = 3.5 to 70 kpc. Everything else held fixed.
**Physical interpretation:** what sets the threshold is how deeply the two galaxies interpenetrate. Small pericentre → deep overlap → high local density → strong dynamical friction (drag ∝ ρ) → lots of orbital energy lost → they bind. Wide pericentre → only grazes the low-density outskirts → little drag → flyby. (The orbit is parabolic, so speed isn't an independent knob: a closer pass is actually *faster* at pericentre — the density term dominates.)

> The one knob here is the pericentre — how close the two galaxies come on their
> first pass. Baseline is five code units, about seventeen kiloparsecs, and I sweep
> it from three and a half all the way out to seventy. Watch the difference across
> the panels. At the smallest pericentre — a near head-on grazing — they merge
> almost immediately, by about sixty time units, and throw out the strongest tails.
> At the baseline it takes roughly twice as long. Open it up to thirty-five
> kiloparsecs and the merger is delayed all the way to six or seven hundred time
> units. And at seventy kiloparsecs it's a clean flyby — they sail past each other
> and never come back, even when I run it out to nine hundred. So there's a sharp
> threshold somewhere between ten and twenty. The physics behind it is dynamical
> friction, and the key quantity is how deeply they interpenetrate. A close pass
> drives each galaxy deep into the other's dense interior, where the drag — which
> scales with the local density — is strong, so plenty of orbital energy gets bled
> off and they bind. A wide pass only skims the thin outskirts and loses almost
> nothing. One subtlety if it comes up: the orbit is parabolic, so the closer pass
> is actually the faster one at pericentre — it's really the density, not the speed,
> that's doing the work.

---

## 15. Axis 2 — Inclination: what shapes the tidal tails? · ⏱ 1:05 (cum 12:05)

**Parameter changed:** inclination of disk 2 relative to the orbital plane — baseline 30°, swept 0, 30, 60, 90, 180°. The orbit itself is unchanged.
**Physical interpretation:** tidal tails are a spin–orbit resonance — only disk material co-rotating with the encounter is lifted into long tails; counter-rotating material barely responds.

> Now I keep the orbit fixed and instead tilt the second disk. Baseline is thirty
> degrees; I go from zero — fully coplanar — all the way to a hundred and eighty,
> anti-aligned. The coplanar, prograde case throws the longest and most symmetric
> tails — these are the textbook Toomre tails. As I tilt the disk out of the orbit
> plane, the tails get shorter and more diffuse, and by a hundred and eighty degrees
> — retrograde — they almost completely vanish. The reason is a resonance. A tail
> forms when disk stars orbit in the same sense as the encounter, so the companion's
> tidal pull acts on them coherently for a long time and flings them out.
> Counter-rotating stars sweep past the perturbation in the opposite direction, the
> force reverses before it can do any real work, and you get almost no tail. So the
> spectacular tail systems — the Antennae, the Mice — aren't typical mergers at all.
> They're prograde survivors. Most real mergers have misaligned disks and look far
> less dramatic.

---

## 16. Axis 3 — Mass ratio: when does a merger become accretion? · ⏱ 1:05 (cum 13:10)

**Parameter changed:** mass ratio `M₂/M₁` — baseline 1:1, swept 1:2, 1:4, 1:8. The secondary's particle count scales with its mass too, so it is structurally the same galaxy, just lighter.
**Physical interpretation:** "merger" and "accretion" are not two categories but the two ends of one continuum, and the mass ratio is the dial.

> Here the knob is the mass ratio. The baseline is two equal galaxies; I shrink the
> second one to a half, a quarter, then an eighth of the primary — and crucially I
> scale its particle count down with its mass, so it's the same galaxy, just
> lighter. At one-to-one, both disks are destroyed and you get the elliptical we saw
> a moment ago. At one-to-eight it's a completely different story: the primary
> barely notices. The little secondary gets tidally shredded as it spirals in, and
> it's fully absorbed by about eight hundred time units, leaving the big disk
> essentially intact. So "merger" and "accretion" aren't two separate things —
> they're the two ends of one continuum, and the mass ratio is the dial that moves
> you along it. And this is exactly what matters for how galaxies actually grow:
> minor mergers like the one-to-eight dominate the cosmic merger rate. They build up
> a galaxy's central bulge without ever destroying its disk.

---

## 17. Axis 4 — Number: what if it's not just two? · ⏱ 1:00 (cum 14:10)

**Parameter changed:** number of galaxies — 2 → 3, equal mass, on an equilateral triangle with a slightly sub-Keplerian tangential speed so they spiral inward.
**Physical interpretation:** two bodies on a parabola are clean and deterministic; a third turns it into the chaotic three-body problem, so the path to the merger is sensitive to initial conditions.

> The fourth knob is just: what if there are three? I place three equal galaxies on
> an equilateral triangle, with a tangential speed a little below circular so they
> gradually spiral inward. Remember the two-body merger was clean and deterministic
> — a parabolic orbit you can write down on paper. Add a third body and you're in
> the three-body problem: there's no closed-form orbit, and which pair pairs up
> first is genuinely sensitive to the setup — change the random seed and you can
> change the order of events. They do all eventually merge, but into a remnant
> that's hotter and more thoroughly phase-mixed than any one-to-one collision —
> there's simply more shuffling per relaxation time. And real compact groups do
> exactly this: Stephan's Quintet is the famous example, and the end product is
> generally more violent than a tidy sequence of pairwise mergers.

---

## 18. Axis 5 — Collision axis: a head-on hit makes a ring · ⏱ 1:20 (cum 15:30) ★ flagship

**Parameter changed:** collision axis — baseline in-plane (along x); now perpendicular (along z) AND head-on (pericentre 0). The orbit is tipped out of the disk plane so the galaxies meet face-on, dead-centre.
**Physical interpretation:** the morphology is set not just by *what* collides but by the *geometry of how* — a perpendicular, central impact makes a ring that no in-plane geometry can.

> This is the one I want you to walk out remembering. Every axis so far kept the
> orbit in the disk plane. Now I tip it ninety degrees — the orbit is perpendicular
> to the disks, so they hit each other face-on — and I make it dead-centre,
> pericentre zero. A true head-on smash. Watch what happens.
>
> **[Pause — let the gif play, ~10 seconds. Let the ring open and expand on screen.]**
>
> The disk falls inward, rebounds, and a clean hollow ring opens up and marches
> outward. That's a Cartwheel galaxy — and nothing in any of the planar mergers we
> just saw could produce it. Look at what actually changed: it's the same code, the
> same two galaxies, the same gravity. The only thing I changed was the direction of
> approach. In the plane you get tails and an elliptical; perpendicular and central,
> you get a ring. So the morphology is set not just by what collides, but by the
> geometry of how. Let me show you why a ring forms at all.

---

## 19. Why a ring? — the impulse mechanism · ⏱ 1:15 (cum 16:45)

**Parameter changed:** none — this explains the flagship result.
**Physical interpretation:** the ring is a radial density wave, not a static structure; the outward-marching Σ(R) bump is the proof.

> So why a ring, physically? When the intruder passes straight through the centre,
> every star in the target disk gets a sudden inward tug toward the point of impact
> — an impulse that scales like two G times the intruder's mass, divided by the
> impact parameter and the passage speed. All those stars start falling inward
> together. But they overshoot the centre, and because they were all kicked
> coherently, they come back out in phase — and that synchronized rebound is a
> density wave that travels outward through the disk. The clinching evidence is on
> this plot: the surface-density bump doesn't sit still, it moves outward in radius
> over time. That's how we know it's a genuine wave and not a frozen feature. It
> peaks at a radius of about ten to twelve kiloparsecs and expands at roughly fifty
> kilometres per second — right in line with the real Cartwheel. One honest caveat:
> because our two galaxies are equal mass, both of them ring symmetrically. The real
> Cartwheel is a small intruder punching through a big target, so in nature only the
> big target rings.

---

## 20. Summary · ⏱ 1:00 (cum 17:45)

**Parameter changed:** all five axes, recapped.
**Physical interpretation:** one code + one baseline mapped a five-dimensional slice of merger outcomes.

> So, pulling it together: from one code and one baseline encounter, we mapped five
> axes of the merger parameter space. Pericentre gives a sharp merge-or-flyby
> threshold, set by dynamical friction. Inclination controls the tidal tails —
> prograde makes them, retrograde kills them. Mass ratio turns a merger smoothly
> into an accretion event. Going from two galaxies to three turns clean, predictable
> dynamics into three-body chaos. And the collision axis is the dramatic one — tip
> the orbit perpendicular and a head-on hit makes a Cartwheel ring. The honest
> limits: our model is collisionless, so there's no gas and no star formation; it's
> finite-N, which heats the disk a little; and the production runs use a direct sum
> rather than a live tree. Each of those is a clean next step — and a good place to
> start the conversation.

---

## 21. Key references · ⏱ 0:15 (cum 18:00)

> Two names worth saying out loud: Toomre and Toomre, who explained the tidal tails
> back in 1972, and Hernquist, whose method is how we built equilibrium galaxies in
> the first place. Both thread through this entire talk. Thank you — happy to take
> questions.

---

### Timing recap (Part 2)

| slide | budget | cum |
|---|---|---|
| 13 divider | 0:30 | 9:50 |
| 14 pericentre | 1:10 | 11:00 |
| 15 inclination | 1:05 | 12:05 |
| 16 mass ratio | 1:05 | 13:10 |
| 17 three-galaxy | 1:00 | 14:10 |
| 18 Cartwheel finale | 1:20 | 15:30 |
| 19 ring mechanism | 1:15 | 16:45 |
| 20 summary | 1:00 | 17:45 |
| 21 references | 0:15 | 18:00 |
| 22 thank you | 0:30 | 18:30 |

Part 2 content ≈ 8:40 (within the 8–10 min target). If you need more Q&A buffer,
the cheapest trims are three-galaxy 1:00 → 0:50 and reading references as a flash.
