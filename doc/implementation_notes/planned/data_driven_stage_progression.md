# Data-driven stage progression in solve3D

Implementation note, 2026-10-19. **Proposal; nothing implemented.** Written
after PfCRT `latest9` (run at `8c8d8b42`), where the two progress-measuring
mechanisms added on 2026-10-18 overruled the planned stage schedule in every
restart. The question this note asks is not how to recalibrate the ladder, but
what the ladder should be *for*.

## The ladder is a prediction used as a schedule

`lpstages` builds the per-stage matching limits by stepping an FRC criterion
down from `frclims(1)` to `frclims(2)`, with `calc_lpinfo` converting each
criterion to a limit, and a linear `box_crop` ramp alongside it
(`simple_estimate_ssnr.f90`). Those FRCs are **2D class FRCs, measured before
any 3D refinement has happened**. The resulting rungs then carry a schedule's
authority: each decides a stage's band, its crop and, with `MAXITS`, how long
the stage runs.

So the schedule is a prediction, made from a different modality, about a
trajectory it cannot observe. Recalibrating its constants is working on the
wrong object; the mismatch is structural.

## What PfCRT measured

Planned limits, identical across all ten restarts:

| stage | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| lp (A) | 20.0 | 19.4 | 18.3 | 12.9 | 11.1 | 9.7 | 9.1 | 6.0 |

Three rungs sit inside the 2 A span 11.1-9.1, then a 3 A gap to the last.
Against that, the measured behaviour:

- **Stage 6 was skipped in 10 of 10 restarts** and stage 4 in 2 of 10; the
  effective schedule was 1, 2, 3, 4, 5 (briefly), 7, 8.
- **Stage 5 exited early in 10 of 10.** Once NU filtering is on the map
  crosses 10 -> 7.6 A in three iterations, passing all three bunched rungs at
  once.
- **The productive work is not where the iterations are.** Stages 1-4 spend
  72 iterations to reach about 10 A; stage 4's fifteen unfiltered iterations
  bought 0.33 A (10.350 -> 10.016); stage 5's *three* NU-filtered iterations
  bought 2.3 A (10.016 -> 7.762).
- Per-iteration partition cost rises only ~6x across the ladder (3.02 s at
  `box_crop` 88 to 18.10 s at 160, from `latest8`'s `BENCH`), so iteration
  count, not stage cost, is what the schedule is really allocating.

The schedule was therefore already being driven by measured resolution, with
the rungs acting as an upper bound that the mechanisms routinely overrode.
That is the awkward middle: a bound with a schedule's authority.

**The pre-NU stages have a job, and it is not resolution.** `nu_stage` 3 and
4 were tried on PfCRT and streptavidin and did not work; 5 is the earliest
the NU machinery can be inserted, and `nu_stage` now rejects anything lower.
The reason is that NU filtering needs an alignment good enough for its
local-resolution estimate: below that the estimate is meaningless, the
filtered reference reinforces the bad alignment instead of correcting it, and
the run goes into a hole.

That changes how the figures above should be read. Stages 1-4 are not 72
wasted iterations with a missing filter; they exist to earn an alignment the
filter can use, and they are the prerequisite for everything stage 5 onward
does. The surplus, if any, is *within* that job rather than a case against
it, and the question is whether those stages need that many iterations to
finish it -- measured on alignment quality, not on resolution. The evidence
that something is loose there is alignment-shaped: on the `latest6`
reference, stage 3's orientation overlap peaks at its fourth iteration
(0.50-0.52) and then declines monotonically for the remaining thirteen
(to 0.43-0.45) in every restart, and stage 1's overlap plateaus by its fourth
of twenty.

So a progression rule needs **two signals, not one**: alignment quality
(overlap) governs the pre-NU stages, measured resolution governs the stages
from 5 on. Keying the early stages to resolution would be the same category
error as the reverted `params%lp` exit -- steering on a quantity that is not
what the stage is there to improve.

## The proposal

Make the roles explicit. The ladder supplies:

- `lpstart` and `lpfinal`, as now;
- a **monotone ceiling** per stage, so the band can never run ahead of what
  the data supports;
- `box_crop` for whatever stage the run is actually in.

Progression between `lpstart` and `lpfinal` is then driven by measured
resolution rather than by a rung count: advance when the map has earned it,
stay while it is still gaining. Stages stop being fixed destinations and
become states the run passes through at its own rate. The existing
`solve3D_stage_behind_fsc05` (advance) and the `lpexit` gain guard (stay)
are already the two halves of that rule.

`f4d6accff` moved in an adjacent direction for the `vol1`-supplied path:
`lpstages_setlims` gained a `first_stage` argument and a linear lp schedule
from that stage onward. That is the explicit-limits variant, not the
FRC-driven `lpstages` the default path uses, but it is the nearer precedent
and this proposal should build on it rather than add a third scheme.

## What has to be true first

- **The gain guard has to be proven.** Progression keyed to measured
  resolution is only as good as the measurement, and `res05` is the signal
  that has already produced two regressions on this dataset: the FSC=0.5
  crossing is a ~40x noise amplifier on PfCRT's plateau, and the first
  `lpexit` fired on floating-point noise. The guard added on 2026-10-19
  (`res05` must have stopped improving, compared within the stage) needs two
  or three runs behind it before anything larger leans on it.
- **The ceiling must stay a hard bound.** Nothing about this proposal
  relaxes it: a fully data-driven band on a noisy signal is exactly the
  failure mode of `65bee3073` and `latest8`, in the other direction.
- **FSC=0.143 stays out of it.** solve3D has no gold-standard stage, so the
  0.143 crossing is inflated by alignment the halves share and must not steer
  a matching band. The progression signal is FSC=0.5, as the current
  promotion already uses.
- **Multi-state needs its own answer.** Its default schedule is five stages
  (`NSTAGES_INDEPENDENT`), and both current mechanisms exclude the last
  stage, so a progression rule that can consume stages would need an explicit
  floor on how few it may run.

## Not proposed

- Recalibrating the FRC criterion constants or the rung count. The rungs
  being badly spaced for PfCRT is a symptom; a hand-tuned ladder is still a
  prediction from 2D FRCs.
- Any change to `MAXITS`. A stage that cannot be cut while it is still
  gaining already makes its budget a ceiling rather than a quota, which is
  most of what a budget change would have bought.

## Caveats on the evidence

`latest9` ran at `8c8d8b42`, and HEAD is ten commits past it. Two of those
touch what is measured above: `f4d6accff` (frequency marching, and
`lpstages_setlims`) and `af13bc828` (reusing the search's crop in 3D
reconstruction, which bears on the `box_crop` and cost figures). The rung
values and per-iteration costs should be re-measured before they are used to
size anything. The behavioural findings -- stage 6 skipped in 10 of 10, the
NU stage exiting mid-climb in 10 of 10 -- are about the mechanisms and do not
depend on those commits.
