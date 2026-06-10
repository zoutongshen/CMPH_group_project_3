# Project 3 presentation — reveal.js deck

20-minute talk: *Galaxy–Galaxy Collisions, a from-scratch self-gravitating N-body study.*
Built as a single offline HTML deck so the gifs autoplay and the equations render
without any network connection.

## Present it

```bash
open index.html          # macOS — opens in your default browser
```

Then:

- **Arrow keys / space** — next/previous slide.
- **`S`** — open **speaker view**: current + next slide, my time-cue notes, and a
  running **timer**. Use this to rehearse to 20 min. Each slide's notes start with
  `⏱ <slide time> (cum <cumulative>)` so you can tell at a glance if you're on pace.
- **`F`** — fullscreen. **`Esc`** — slide overview. **`B`** — black the screen.
- **Backup slides** sit after "Thank you": press **↓** to descend into them during Q&A,
  **←** to climb back out. They never appear in the main left→right flow.

### Pacing checkpoints (say these times to yourself)

| By end of… | clock |
|---|---|
| Workflow roadmap | ~2:00 |
| Setup (gravity / positions / velocities) | ~4:45 |
| Validation (N &amp; ε) | ~5:55 |
| Baseline merger + remnant | ~8:10 |
| Axes 1–4 (pericentre → number) | ~12:30 |
| Cartwheel (collision axis, the finale) | ~15:05 |
| Thank you | ~16:25 |

Leaves ~3.5 min buffer/Q&A inside the 20.

**Structure:** Part 1 is setup-and-validation (build a galaxy → collide two →
check it's physics not numerics); Part 2 is the five encounter axes, each a gif you
narrate, building to the Cartwheel ring as the finale.

## What's still TODO (search the deck for `TODO` / `todo`)

- Confirm the **group number** for project 3 (may differ from projects 1/2) — title slide.
- Optional **observation image** (Antennae / Cartwheel) on the motivation slide, or delete it.
- One spoken sentence per **equation** (the maths is rendered; you supply the interpretation).
- A couple of small `TODO`s for numbers you may want to quote (remnant \(R_e\), tree-vs-N²).

The figures themselves are wired in and final — every `<img>` points at a copy under
`assets/` (verified to resolve).

## Disaster fallback (carry this too)

PDFs freeze gifs on their first frame, but a static PDF is a safe backup if the
presentation machine misbehaves. With reveal's built-in print mode:

```bash
# open this URL in Chrome, then Cmd-P → "Save as PDF", landscape, margins none:
open -a "Google Chrome" "index.html?print-pdf"
```

or, if you have Node:

```bash
npx decktape reveal index.html deck_backup.pdf
```

Last-resort backup for the **flagship Cartwheel clip**: keep
`assets/z_axis/multi_zaxis_headon.gif` open in Preview so you can play it by hand.

## Offline note

`reveal/revealjs/` (reveal.js 5.1) and `reveal/katex/` (KaTeX 0.16) are vendored, and
the 17 figures used in the deck are copied into `assets/` — nothing loads from the
internet and nothing lives outside this folder. The **`presentation/` folder is fully
self-contained** (~28 MB): zip it, copy it to a USB stick, or share it as-is and it just
works. (The originals stay in the gitignored top-level `figures/`; `assets/` is the
committed, presentation-pinned copy.)
