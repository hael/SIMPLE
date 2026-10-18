# The matching band: its contract, a broken clamp, and the test that was missing

Implementation note, 2026-10-18. **Proposal. Only the revert (section 2) is
implemented; sections 3 and 4 are not.** Written after a PfCRT regression at
commit `25d25688` that was caused by a change made two days earlier to the
matching low-pass limit, and found, as the two before it, by a maintainer
looking at density rather than by a test.

Three separate matters are kept apart below, because conflating them is what
made the diagnosis take three attempts:

1. who owns the matching band, and the contract that ownership implies;
2. the clamp that violated it (closed: reverted);
3. the FSC=0.5 promotion (open: needs a decision, not a patch);
4. the test that would have caught (2), and would bound (3).

## 1. Who owns the matching band in solve3D

The band is `params%kfromto(2)`. In solve3D it is owned entirely by the
controller:

- `stage_matching_lp` returns the static stage ladder `lpinfo(istage)%lp`
  (`simple_solve3D_controller.f90:706`);
- `promote_stage_lp_from_fsc05` may sharpen it at a stage boundary
  (`:718`, section 3);
- the result is set on the refine3D command line **unconditionally**
  (`:636`).

Because `l_lpset` is `cline%defined('lp')`
(`simple_parameters_phases.f90:597`) and no route strips `lp` for a solve3D
stage, `set_bp_range3D` always takes its `l_lpset` branch during solve3D. Its
FSC-derived routes -- the NU handoff, `lplim_crit`, the project `lp` field --
are unreachable there. They remain live for standalone refine3D.

**The contract: the band leads achieved resolution.** The ladder deliberately
runs ahead of the current map's FSC so the map has somewhere to climb into.
Clamping the band to the achieved FSC makes the two chase each other: the band
can only advance as the FSC advances, and the FSC only as fast as the band
permits.

This contract is currently expressed in exactly one place, implicitly:

```fortran
self%l_incrreslim = trim(self%incrreslim) == 'yes' .and. .not. self%l_lpset
```

(`simple_parameters_phases.f90:953`) -- the `+10`-shell headroom is switched
off whenever `lp` is set, because when the caller owns the band the FSC must
not second-guess it. Nothing states the rule in prose, and nothing tests it.

## 2. The clamp that violated it (closed)

`65bee3073` added a per-iteration clamp of `kfromto(2)` to the current
FSC=0.143 crossing in `set_bp_range3D`, intending to stop matching from
running past the resolution the data supports. Reverted in `fb4e653ac`,
byte-identical to the pre-clamp file. Three defects, the first fundamental:

- **FSC=0.143 must not decide any matching resolution in a run that is not
  gold-standard.** solve3D has no gold-standard stage by design
  (`simple_solve3D_controller.f90:416`: "gold-standard refinement belongs to
  `refine3D_auto`"), so its halves share alignment and the 0.143 crossing --
  the criterion that presumes independent halves -- is inflated. Feeding it
  back into the band closes a loop: matched noise raises the FSC, which
  sharpens the band, which matches more noise. The 0.143 crossing
  "improving" from 4.44 to 4.03 A while the maps visibly softened is that
  loop's signature. The invariant was already stated at
  `simple_commanders_refine3D.f90:66` -- "banded via lpstop (never lp:
  `l_lpset` would break gold standard)" -- and the clamp set exactly that
  band.

- **It clamped the ladder.** Per section 1 the band in solve3D *is* the
  controller's ladder, so the clamp did not restrain an FSC heuristic; it
  restrained the mechanism whose job is to lead. It fired on **14-29% of
  iterations across ten PfCRT restarts** (18-38 of ~133 each).
- **First-crossing on a non-monotonic FSC.** `get_find_at_crit`
  (`simple_math_ft.f90:190`) returns the first shell scanning outward below
  the criterion. PfCRT's FSC dips near 5.97 A and recovers by 5.17 A, so a
  single noisy shell collapses the band for that iteration.

Recorded in `doc/policies/NU/nonuniform_filtering_policy.md`. Do not
reintroduce it in `set_bp_range3D`.

The separate change in the same commit -- restoring the density window on the
Guinier B-factor fit -- stands, and is unrelated.

## 3. The FSC=0.5 promotion (open; a decision, not a patch)

This is independent of section 2 and predates it.

Past `FSC05_PROMOTE_MIN_STAGE = 2`, every stage boundary re-derives the band
from the project's FSC=0.5 crossing, taking the best-resolved populated state
(`simple_solve3D_controller.f90:718`, `:740`):

```fortran
res05  = project_best_fsc05_resolution(params)   ! minval over populated states
lp_new = max(min(lp, res05), lp_cap)
```

Promoting from 0.5 rather than 0.143 is **required**, not merely
conservative: solve3D is not gold-standard, so 0.143 is not available to it
(section 2). FSC=0.5 is the pre-gold-standard criterion and is legitimate on a
correlated pair.

The guard in that comment -- "only an FSC beyond the previous band is clean"
-- is implemented, through the direction of the test: promotion fires only
when `res05` is finer than the current band, so the crossing it trusts always
lies outside the matched region, where overfitting cannot have reached. That
is what makes promotion safe at all.

A second bound: `lp_cap = LPSTOP_BOUNDS(1) = 4.5 A`, and promotion only ever
sharpens, so it stops firing once the ladder itself is finer than 4.5 A.
PfCRT's final 3.98 A band is the ladder's own floor, not a promotion. The
instability below therefore moves the mid-stage trajectory -- search cost and
how fast the climb advances -- not the final band.

The problem is dataset-specific. PfCRT's FSC does not cross 0.5 steeply; it runs nearly
parallel to the threshold:

| res (A) | 11.9 | 9.7 | 7.96 | 6.9 | 5.97 | 5.17 | 4.85 | 4.70 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| FSC | .893 | .865 | .711 | .563 | .522 | .594 | .555 | .519 |

So the reported crossing is hypersensitive. Applying a uniform offset to the
Sep16 curve (box 300, 1.035 A):

| offset | FSC=0.5 crossing | FSC=0.143 crossing |
| --- | --- | --- |
| +0.00 | 4.63 A | 4.37 A |
| -0.02 | 4.70 A | 4.37 A |
| **-0.03** | **6.09 A** | 4.44 A |
| -0.05 | 6.21 A | 4.44 A |
| -0.08 | 7.06 A | 4.44 A |

A 0.03 FSC difference moves it 1.5 A; FSC=0.143 barely moves. **PfCRT's band
is therefore re-sampled at every stage boundary from the least stable
statistic in the pipeline, on its optimistic tail** (`minval`), for reasons
unrelated to the map improving. That is a plausible standing contributor to
PfCRT's run-to-run spread, and it is why PfCRT and not beta-galactosidase is
the dataset that keeps moving.

**Implemented 2026-10-18: `res05` gates the advance, the ladder sets its
size.** When the FSC=0.5 crossing is finer than the planned
`lpinfo(istage)%lp`, the stage runs at the next rung `lpinfo(istage+1)%lp`
instead of at `res05` itself. No new constants. It is less machinery than
what it replaces: the band is no longer set to an arbitrary measured value.

- A noisy crossing can now cost at most one rung early. The gate itself still
  flips with noise -- `res05` wandering 4.6-7.1 A against an `lp` of 6 A
  flips yes/no -- but a wrong yes is bounded.
- It does not accumulate: each stage starts from its own planned rung, so the
  lead is at most one rung at any boundary.
- The ladder is monotone decreasing by construction (`frc_crit - frc_stepsz`
  per stage in `lpstages`, and the linear fallback likewise), so the next
  rung is always finer.
- Everything in section 1 and the two guards above hold: steers on 0.5 only,
  the `res05 >= lp` return keeps the sharpen-only direction, monotone,
  floored at `lp_cap`.

**Stage shift, same day, and the real value of the approach.** Borrowing the
next rung's `lp` while staying at the stage's own crop is close to inert:
`set_bp_range3D` ends with
`kfromto(2) = min(kfromto(2), nyqcrop_ind)`, so the advanced limit is clipped
back to the crop Nyquist of `lpinfo(istage)%box_crop`. The stage banner then
prints a limit the search never used, and the stage still spends its full
budget at a sampling that cannot represent the resolution already achieved.

So where a stage can be skipped, skip it: `solve3D_stage_behind_fsc05` gates a
`cycle` in the frequency-marching loop. That drops the whole stage config --
`box_crop`, `smpd_crop`, `trslim`, `nspace`, `maxits` -- which is where the
wall-clock goes. Unlike the rung advance this compounds naturally: a map three
rungs ahead skips three stages. Excluded are the stages carrying
one-time work: `start_stage`, the last stage, and `SYMSRCH_STAGE` when an axis
is actually searched. `NOSYMSRCH_STAGE` is shiftable, because `l_symran` is
set only when `pgrp_start` outranks `pgrp` and so implies `l_srch4symaxis`:
with no axis search, `symmetrize` does nothing at all. The rung advance is
retained for the excluded stages.

**Mid-stage exit, same gate set.** The boundary test only fires between
stages, so a stage that becomes correct to leave halfway through still runs to
convergence. The internal `lpexit` key (emitted for the shiftable stages only)
adds the same test to `check_conv3D`, against `params%lp` -- the band in use
after the crop-Nyquist clamp, not the printed one. It reuses
`res05_state_avg`, which that function already computes, and only ever adds
convergence, so `minits` and the Cartesian rule are untouched.

This is within the policy's never-per-iteration rule, not an exception to it.
That rule exists because an FSC measured *inside* an opened band is
contaminated by noise the alignment just fitted; the test here asks only
whether the crossing lies *beyond* the band, which was never matched. It also
self-brakes: once the next stage opens the band past the crossing, the test
fails again until the map improves.

Compatibility was free: `emitted_lp` is initialised to `-1.` and the stage
report already does `if( emitted_lp(i) < 0. ) cycle`, while the add-on
manifest already prints "NOT RUN BY THE BASE RUN" for that sentinel. Unrun
stages were already modelled.

This is a real behaviour change -- no constant recovers the old output, so
there is no identity test for it. It needs a before/after robustness profile
on beta-galactosidase (steep FSC: the band should land close to where it
lands today) and PfCRT (plateau: the across-restart spread of the mid-stage
band should drop). Note it is also strictly more conservative than what it
replaces: a map far ahead of schedule no longer leaps, it gains one rung.

**Rejected: promoting from FSC=0.143** (with or without a backing-off
margin). Not a cost question -- 0.143 presumes gold-standard halves, which
solve3D does not have, so it must not decide any matching resolution here.
See section 2. The absence of a `res0143` project field is moot.

**Not folded in**, to keep the change single-purpose:
`project_best_fsc05_resolution`'s `minval` over states (the optimistic tail),
and logging the per-stage `res05`/rung/actual triple so the trajectory is
diffable. The second is worth doing, with the run-summary sidecar.

This changes results. It is judged by comparing robustness profiles before
and after on cases with known truth, per
`doc/implementation_notes/planned/test_suite_truth_and_robustness.md`, not by
a single PfCRT run. **A corollary applies either way: FSC=0.5 must not be a
pass/fail gate for PfCRT**, since on a plateau it is a ~40x amplifier of
noise. cFAR and the FSC=0.143 crossing are the stable metrics.

## 3a. This had been tried twice before, and written down

`doc/policies/3D/solve3D_policy.md` already recorded the experiment, in the
NU-stages paragraph: "The band therefore leads the FSC by 1.25-1.5x, which is
what carries the climb: with the band at the previous FSC=0.143 crossing the
alignment was starved of the shells where the merged reference still holds
signal and **PfCRT converged in 0/4 and 1/9 restarts (2026-09-16/18), with a
1.5x lead in 10/10 (2026-09-16)**, and a ceiling only pins the map."

And two band ceilings had already been tried and retired: the class-FRC
`lpfinal` (6.0 A on PfCRT) pinned the NU stages at 5.97 A on 2026-09-07, and
the ladder's hard 4.5 A bound pinned the 2026-09-08 run at exactly 4.50 A for
30 iterations while the handoff asked for 4.14 and 3.98 A.

The clamp of section 2 was the third attempt at the same thing. The knowledge
existed, in writing, in the policy file that owns the band. That is the case
for section 4: the mechanism has to be one that does not depend on someone
having read 900 lines of policy before touching `set_bp_range3D`.

## 4. The test (proposal)

The broader suite proposal's Progress section reads "Nothing yet", and no
tester asserts anything about the band: every `kfromto(2)` in a `*_tester`
module is a fixture allocation. So a change to the band contract has exactly
one detector -- a ~4,500 s run on a real dataset, read by eye, weeks later.
That is the whole of the mechanism today.

Following that proposal's rulings (every high-level test is a SIMPLE
commander, no outside driver, correctness on simulated data with known truth,
floors with a stated origin):

**Claim.** During a solve3D climb the matching band is at or ahead of the
achieved FSC=0.143, never behind it.

**0.143 as an observable, not a control input.** Section 2 forbids FSC=0.143
from *deciding* a matching resolution in a non-gold-standard run. The test
reads it to *detect* a regression, which is a different use: nothing in the
run's behaviour depends on it. The inflation works in the test's favour here
-- an inflated 0.143 reads finer than the truth, so "the band leads 0.143" is
a weak lower bound on the band, and a violation of it is unambiguous.

**Fixture.** Hermetic: an embedded model or phantom, no user-supplied
dataset. A short climb is enough -- the contract is per-iteration, so the
assertion does not need a converged map.

**Assertion.** At every iteration, `params%kfromto(2)` is greater than or
equal to the FSC=0.143 Fourier index of the current map. Reported per
iteration so a failure names where the band fell behind.

**What it catches.** `65bee3073` by construction: the clamp set
`kfromto(2) = fsc0143_ind` exactly, so the assertion fails on equality's wrong
side on 14-29% of iterations. It would also catch a recurrence of the
2026-09-18 static-bank revert, which tightened the handoff by the same
mechanism from the other direction.

**Floor.** The assertion is an inequality, not a threshold, so it needs no
calibrated floor -- which is why it can be a correctness gate rather than a
robustness profile. If a tolerance proves necessary (interpolation at the
crossing), its origin is stated in the code and it is set from a profile of at
least 20 runs, never from one.

**Bounding section 3.** A second, cheaper assertion: feed
`project_best_fsc05_resolution` a plateau-shaped synthetic FSC and fail if a
small perturbation moves the promoted band by more than a stated margin. This
turns the 40x amplification from a dataset surprise into a test failure, and
generalises past PfCRT. It is a test of the promotion's stability and is
independent of which option in section 3 is chosen.

**Registration.** Reachable through `simple_test_exec`, wired as the existing
`*_tester` modules and `simple_test_ui_class` entries are, and passing
`scripts/check_test_registry.py`.

## 5. Not proposed here

- No change to the band logic, the controller, or `set_bp_range3D`. Section 3
  is a decision to be taken on evidence; section 4 does not depend on it.
- The real-data manifest layer, cFAR and B-factor gates, and the rest of the
  dataset tier belong to the broader suite proposal, not to this note.
