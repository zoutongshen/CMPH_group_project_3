# Speaker notes — Part 1 (title → remnant)

Full draft notes for the slides **before** the "Beyond the baseline merger"
divider, archived here when the in-deck notes were reduced to skeleton bullets
for Zhaoyang to fill in his own words. Time cues stay in the deck.

---

## 1. Title

⏱ 0:20 (cum 0:20). One breath. "We built a galaxy from scratch, collided two
of them, and then explored how the outcome depends on the encounter." Go.

## 2. Motivation — Galaxies collide, and it shows

⏱ 0:50 (cum 1:10). Hook with the real systems. The point: one force (gravity)
produces a zoo of morphologies — what selects which? That's what we'll map.
The Cartwheel photo is the JWST poster hanging right outside the Sterrewacht
offices — "we'll make this exact object appear in our own simulation by the
end of the talk."

## 3. The workflow

⏱ 0:45 (cum 1:55, target ~2:00). The map for the talk: build → collide →
validate → explore. First half is how it's set up; second half is the five
axes, each a movie. Don't enumerate the axes here — they each get a slide.

## 4. Setup I — gravity & integration

⏱ 1:00 (cum 2:55). Two ideas only: (1) softened direct sum — ε is why a close
2-body kick can't blow up; (2) symplectic leapfrog so energy doesn't drift.
Flag ε now — it comes back in validation as one of the two knobs.

## 5. Setup II — building the galaxy: the mass

⏱ 0:50 (cum 3:45). Where the mass goes: pull positions from the three density
profiles — disk is a flat exponential, bulge & halo are spherical. Halo is a
truncated isothermal (flat rotation curve). Don't read the formulae — name each
component and move on. Next page: here's what that looks like.

## 6. Setup II — the built galaxy

⏱ 0:30 (cum 4:15). One beat — "this is the galaxy those profiles produce."
Face-on shows the exponential disk; edge-on shows the thin disk embedded in the
round bulge + halo. Then: positions alone aren't enough — we need velocities.

## 7. Setup III — building the galaxy: the motion

⏱ 0:55 (cum 5:10). The subtle half of the IC: each component needs velocities
matched to the SAME potential. Disk spins (Toomre Q ≳ 1 keeps it from
fragmenting); bulge & halo are pressure-supported via Jeans. This is the part
that's easy to get wrong — motivates the next slide.

## 8. Does it hold up?

⏱ 1:10 (cum 6:20). The "we turned a bug into a result" slide — brief but
important. The disk puffs up; we proved it's finite-N relaxation, controlled
by N and ε, NOT a broken IC. Quadrupling N cuts heating by the predicted
factor. Also point at the rotation curve — flat, peaking ~240 km/s: it
reproduces the Milky Way, an independent literature check, not just stability.
This buys trust for every result after. Don't dwell — converged numbers, move on.

## 9. Comparing the two galaxy models

⏱ 0:45 (cum 7:05). Compare the galaxy setup, not the method. Both build the same
three-component galaxy with the same density-profile SHAPES; the second is an
independent Barnes–Hut tree implementation — mention the method only if a
professor asks. The one input that differs: the tree model uses a lighter halo
(M_h = 3.0 vs the manual's 5.8) — same shape, ~half the mass — plus cold-collapse
ICs. Result: our curve is flat & Milky-Way-like (~240 km/s); the lighter halo
shifts the whole curve down a near-uniform ~30–40 km/s, just below the MW line
(it's a √M scaling of the dominant halo term, so it reads as a parallel offset).
Caveat if pressed: this figure proves the MASS-MODEL difference only — it is NOT
a tree-vs-direct-N² test, so don't claim the two methods are validated against
each other. (Partner's slot — hand off here if he presents it.)

## 10. The baseline merger

⏱ 1:15 (cum 8:20). Set it up in one line — two equal galaxies dropped on a
parabolic orbit — then LET IT PLAY. Stop talking ~10 s; name the phases as they
happen: first passage → tails throw out → cores sink → merge. Everything in the
second half deviates from THIS.

## 11. What's left behind — the remnant

⏱ 1:00 (cum 9:20). The classic result: a major merger of two disks makes an
elliptical. The R^1/4 straight line is the money plot — point at it. That
closes the "build & collide" half; now we turn the knobs.
R_e ≈ 4.8 kpc comes from the fit itself (R_e = (3.33/|slope|)^4) — say "the right
size for a 10^11 M_sun elliptical" (the size–mass relation). The model-independent
half-mass radius (5.6 kpc) agrees → a robustness cross-check; don't call IT R_e.
And remember: the STRAIGHT LINE is the proof it's an elliptical — R_e is just its
size, only meaningful because the R^1/4 fit holds.
