# solve3D on small membrane proteins: diagnosis and experimental plan

Archived in `completed/` on 2026-11-17 by maintainer decision. The status and open items below are as
last written.

Status, 2026-11-16: in progress. Phase 1 read as Hans reported it (`lpstart=16.25` worse on PorF and PfCRT,
the export gate unchanged; the start cap is rejected) and the band-versus-sampling control run (Decision 13).
On Hans's rulings of 2026-11-16 the tree now holds the release fix (Decision 3, committed in `360fbeb4f`), the
late release through stage 4 (Decision 4, amended) and the best-first sticky draw (Decision 14), the latter two
uncompiled. Base: `603c44490`. Supersedes: the first analysis under this name and its
independent review, both written and left uncommitted on 2026-11-15 and merged here (section 2.8 lists what
the review corrected; the plan adopts its recommendations). Depends on: the logs under
`~/Processing/solve3D_regression/{PfCRT,PorF,exp_gate}/{current,last_working_*}/` on the Mac (read); the
`lpstart=16.25` restarts of PorF, the export gate and PfCRT that Hans started on 2026-11-15 (running, not
read when this was written). Related: [data-driven stage progression](../rejected/data_driven_stage_progression.md),
sections 3-4 of the [solve3D policy](../../policies/3D/solve3D_policy.md) and "Sticky sampling schedules"
in the [fractional-update policy](../../policies/importance_sampling_fractional_update_policy.md).

Terms. A *rung* is a stage's planned low-pass limit (`lpinfo(istage)%lp`); the *band* is the matching
limit an iteration used (`MATCHING LOW-PASS LIMIT`). *Overlap* is the fraction of the particles drawn in
an iteration whose projection direction stayed within the angular threshold of their previous one
(`mi_proj`); the *angular change* is the log's `DIST BTW BEST ORIS` (a running average per particle).
*Effective poses* is the iteration log's count of poses carrying a particle's posterior (1 = sharp). The
*sticky subset* is the class-balanced set of `nsample` particles that stages 1-2 refine (`sticky_sampling`);
its *release* is the first stage that samples freely, stage 3 (the add-on's *cohort* is a different set). A
*first-time* particle is drawn for the first time in the run, so its "previous" pose is its starting
assignment. A stage *consolidates* when its overlap rises and its posteriors sharpen. *One Fourier shell* is
the next coarser shell at the final reconstruction's box and pixel size.

## 1. Why

At `603c44490` one of five completed PorF restarts ends at 7.43 Å (FSC=0.143), against 4.11-4.22 Å in the
other four and in all ten restarts of `fa9b82b5`; PorF has to succeed in ten of ten. All three completed
export-gate restarts end at 7.17-7.40 Å, against 3.88-3.95 Å in the nine completed restarts of `3fb0b00b`
(on a different input, section 2); in those runs the helices were visible after stage 3, and now they
never appear. PfCRT is back at its reference (4.09-4.14 Å, against 4.14-4.25 Å at `db20cde7`) since the
in-plane posterior insertion was switched off by default (`603c44490`). Streptavidin, beta-galactosidase
and TRPV1 are perfect by Hans's reading. TRPV1 and the hard targets all have mask diameters of 150-170 Å,
so the size-planned ladder gives them the same rungs and cannot tell them apart.

## 2. What the logs show

Read from the `SOLVE3D_OUTPUT_RESTART*` logs, which print their commits, with a parser kept outside the
repository (Phase 2 puts a stage summary into the log itself).

| Target | Last working | Current (`603c4449`) | Active particles | Sampling |
| --- | --- | --- | --- | --- |
| PfCRT | `db20cde7`, 10 restarts | 5 restarts, 1 partial | 16,778 then, 16,832 now | full: `nsample` 25,000 is above 0.9 of the active particles |
| PorF | `fa9b82b5`, 10 | 5, 1 partial | 37,502 in 99-100 class units, both | fractional, `nsample` 10,000 |
| Export gate | `3fb0b00b`, 10 (9 completed) | 3, 1 partial | 129,650 in 94 classes then, 111,099 in 60 now | fractional, `nsample` 10,000 |

The stage profile, the baseline every phase compares against (bands and FSC=0.5 in Å; partial restarts
count where they reached the stage; "consolidated" is the readout of section 6: overlap at least 0.30 and
angular change at most 30 degrees at the end of stage 3 or 4):

| Target, build | Stage 2: band, FSC=0.5 at its end | Stage 3: band; at its end FSC=0.5, overlap | Stage 4: band; at its end FSC=0.5, overlap | Consolidated | Final FSC=0.143 |
| --- | --- | --- | --- | --- | --- |
| PfCRT `db20cde7` | 18.3, 18.3 | 18.3; 13.5, 0.05-0.06 | 13.5; 12.4-12.9, 0.43-0.45 | 10 of 10 | 4.14-4.25 |
| PfCRT `603c4449` | 18.3, 18.3 | 16.3; 12.9-13.5, 0.09-0.10 | 12.9-13.5; 12.4-12.9, 0.44-0.52 | 6 of 6 | 4.09-4.14 |
| PorF `fa9b82b5` | 19.5, 22.3-24.0 | 17.3; 12.5-14.2, 0.13-0.16 | 12.5; 9.5-10.8, 0.35-0.43 | 10 of 10 | 4.11-4.22 |
| PorF `603c4449` | 19.5, 22.3-24.0 | 17.3; 12.0-15.6, 0.13-0.17 | 12.0-14.9; 9.5-14.2, 0.25-0.47 | 5 of 6 | 4.11-4.22; 7.43 |
| Export gate `3fb0b00b` | 9.1, 13.2-13.9 | 8.5; 7.6-8.2, 0.70-0.82 | omitted in 8 of 10 | 10 of 10, in stage 3 | 3.88-3.95 |
| Export gate `603c4449` | 16.9, 18.2 | 15.8; 12.5, 0.01 | 12.5; 10.3-11.8, 0.01-0.02 | 0 of 4 | 7.17-7.40 |

### 2.1 PfCRT

Same data but for 54 particles, full sampling, no sticky subset. Both builds churn in stage 3 (overlap
0.05-0.10, angular change 44-60 degrees) and consolidate in stage 4 at 12.9-13.5 Å, the FSC=0.5 value. The
posteriors sharpen with the stage-3 band: 25-31 effective poses at 16.3 Å (`603c4449`), 49-55 at 18.3 Å
(`db20cde7`, whose rule could coarsen the stage). PfCRT needs no fix here; it is the control the ladder
change must not break.

### 2.2 PorF

Same input and ladder in both builds (mask 170 Å: 21.25, 19.07, 16.89, 14.71, 12.54, 10.36, 8.18,
6.00 Å). Every completed success consolidates in stage 4 at a band of 12.0-14.2 Å (overlap 0.33-0.47 and
3.8-6.2 effective poses at its end). The failed restart (restart 5 of `603c4449`) left stage 3 at
FSC=0.5 15.6 Å, the coarsest of the sixteen, with the worst FSC=0.143 (11.6 Å against 9.2-10.1 Å). Being
coarser than stage 4's rung, that FSC kept the rung (a 14.86 Å band): FSC=0.5 stayed at 14.18 Å in all 15
iterations, the overlap reached 0.25 (10.4 effective poses, 38 degrees), stage 5 reached 9.2 Å and stage 6
plateaued at 9.2/7.4 Å. PfCRT failed in the same way at `fa9b82b5` (restart 8, described in `db20cde75`'s
message): stage 4 at 14 Å, stalled at overlap 0.36 and 31 degrees.

The rule of 2026-11-13 (`set_stage_lp_from_fsc05`, the FSC=0.5 value when finer than the rung since
`a4881e6c4`) ran stage 4 coarser than the one-rung advance of `fa9b82b5` in three of the five other
restarts (13.6-14.2 Å against 12.5 Å). From the failed restart's stage-3 endpoint the advance would have
chosen the same 14.7 Å; whether `fa9b82b5` would have reached that endpoint is unknown, since no restart was
replayed. So the failure matches a mechanism seen before, and the current rule leaves several successful
restarts a smaller margin.

### 2.3 The export gate: the band

`3fb0b00b` planned its ladder from the 2D input's class FRCs (`lpstages`: the mean class FRC was 0.97 at
20 Å and still 0.86 at 9.1 Å), the inflated kind that gold-standard 2D removed on 2026-11-10: 20.0, 9.1,
8.45, 7.89, 7.40, 7.17, 6.96, 4.5 Å.
`603c4449` plans from the mask diameter (`set_lplims_from_size`, 155 Å): 19.38, 17.44, 15.50, 13.56, 11.63,
9.69, 7.75, 5.81 Å. The change was therefore not from FSC-driven to FSC-free. Honest class FRCs would not
give the old rungs (they moved PfCRT's stage 4 from 12.4 to 18.3 Å, the reason the size ladder replaced
them), and the 3D FSC=0.5 cannot lead either: after stages 1-2 it sits at or behind the band (18.2-18.3 Å
for PfCRT and the export gate, 22-24 Å for PorF at a 19.5 Å band). What worked was a band well ahead of the
FSC: the old stage 2 matched at 9.1 Å while its FSC=0.5 was 13-16 Å, and its stage 3 at 8.5 Å
took the overlap from 0.04 to 0.78 in ten iterations. The difference starts in stage 2, before the sticky subset
is released, so the release (section 2.4) cannot be the whole story.

### 2.4 The export gate: the sticky subset's release

Stages 1-2 refine one class-balanced sticky subset of 10,020 particles (9% of 111,099); sticky rescoring
(`sample4update_rescore`) advances its `updatecnt` every iteration, to 40 after 20 + 20 iterations. Stage 3
samples freely, each class's quota lowest `updatecnt` first (`sample4update_class`), and the sampling history
is reset only at the run's first iteration (`sample_ptcls4prob3D`). Three consequences:

- Never-searched particles take every draw until they catch up with the subset's count. With 60 classes at
  a quota of 167, the smallest class (574 particles) needs ceil(40 × 407 / 167) = 98 free iterations before
  its sticky particles are drawn again, larger classes longer; the runs have 48-49. No particle of the subset is
  searched again after stage 2.
- The first stage-3 iteration seeds the trailing chain from its own sample (`VOLASSEMBLE: SEEDED FULL-MASS
  TRAILING CHAIN ... THIS ITERATION USES THE CURRENT SAMPLE ONLY`, the absent-chain path of
  `blend_trailing_accumulators`). The stage-2 map guides that iteration's alignment, but from then on no
  image of the subset is in any map, while the chain's represented population (20,040) counts the subset and
  weights the first stage-3 sample twice.
- Stage 3 is fed by first-time particles aligned against an unconsolidated map, for long: the largest class
  (6,267 particles) is seen once only after 38 iterations (`Draws for one full sweep`), and every particle has
  been drawn only in stage 7.

Restart 1, by stage (each stage read at its last iteration):

| | Stage 3 (15.8 Å), iterations 41-48 | Stage 4 (12.5 Å), 49-55 | Stage 5 (NU), 56-67 | Stage 6, 68-76 |
| --- | --- | --- | --- | --- |
| First-time share of each draw | all to 47% | 40% to 19% | 15% to 8% | 8% to 2% |
| Particles drawn at least once, at the stage's end | 67.4% | 84.6% | 96.8% | 99.9% |
| Effective poses | 150-210 | 128-145 | 151-402 | 14-20 |
| Overlap | 0.000-0.005 | 0.007-0.015 | 0.003-0.013 | 0.050-0.093 |

PorF at a similar band: no first-time particle from the seventh free iteration on, 19-37 effective poses in
stage 3 and 5-16 in stage 4, overlap 0.33-0.47 at the end of stage 4. PorF's sweep is seven iterations, the
export gate's 38, which may explain why PorF tolerates the same release.

The statistics hide no stable population. At the end of stage 4 (iteration 55) 1,908 of the 10,020 drawn
particles are first-time and the overlap is 0.015, so even if all of it came from revisits their overlap is
at most 0.015 / (1 - 0.19) = 0.019. `check_conv3D` does average overlap and angular change over all drawn
particles (mask: sampled); separating the two groups diagnoses, it does not rescue.

The trailing reconstruction blends the current partial sums with the decayed accumulator of all earlier
iterations (`population_blend_weights`: current × u/f, chain × (1 - u)·N/M). A pose's contribution fades
by a factor (1 - u) per iteration; it is not replaced when the particle is drawn again. With u near 0.09 in
the late stages its half-life is about seven iterations.

Which of band and release matters more cannot be separated from these logs. Both differ from the run that
worked, the input differs too, and the fixed subset of `3fb0b00b` was each class's best-scoring particles by 2D
score (10,058 in all; the sampler ranked by score until `cf643a192`, 2026-11-11), consolidated in stage 3
before stage 5 drew new particles. Its release then cost nothing: FSC=0.5 stayed at 7.64 Å in the first
stage-5 iteration, whose chain also started from the new sample alone, and reached 7.17 Å in the next.

### 2.5 The frozen sticky subset on every fractional run

The release problem is not the export gate's alone. The count of 40 keeps much of PorF's sticky subset out of the
free sampling as well, and every active particle with `updatecnt` above 0 enters the final all-particle
reconstruction (`sample4rec` through `sample_rows` in `simple_rec3D_strategy`) with the pose of its last
search; the final missing-update pass searches only the particles with `updatecnt` 0.

| Target | Active | Sticky subset | Largest unit's sweep | Free iterations | Subset never searched again (estimate) |
| --- | --- | --- | --- | --- | --- |
| PorF (restarts 1 and 5) | 37,502 | 10,078 (26.9%) | 7 | 66-74 | about two thirds: 6,300-6,800 particles, 17-18% of all |
| Export gate | 111,099 | 10,020 (9.0%) | 38 | 48-49 | all: 9% of all |
| SLC (logs of 2026-11-05) | 429,600 | 25,000 (5.8%) | unknown | unknown | expected all (Phase 0 checks) |
| PfCRT | 16,832 | none (full sampling) | 1 | not applicable | not applicable |

Estimated from each log's unit table: a unit's sticky particles become eligible again after ceil(40 × rest / quota)
free iterations; units smaller than their quota (13 on PorF) are drawn whole every iteration. These
particles carry stage-2 poses (bands of 16.9-19.5 Å, overlap 0.09-0.12 at the end of stage 2) into the final
map. The fractional-update policy says that the sticky subset lets the shc search scan each selected particle in
full once and then refine it; nothing in it intends the subset to drop out of the rest of the run.

### 2.6 SLC

`solve3D pgrp=c1 mskdiam=130 nsample=25000 inpl_posterior=cone` (`nrestarts=3`, `nparts=10`, `nthr=8`) with
the code of 2026-11-15 gave the best SLC map so far. A 130 Å mask plans 16.25, 14.62, 13.00, 11.38, 9.75,
8.12, 6.50, 4.88 Å: stage 3 at 13.0 Å and stage 4 at 11.4 Å, the band where PfCRT and PorF consolidate. The
mask, `nsample` and the insertion differ from the defaults at once, and the run's particle count and stage
logs are not in hand, so it supports an experiment, not a constant. The export gate consolidated at 8.5 Å
when it worked.

### 2.7 Established and open

Established by the logs:

- No target consolidates at the size ladder's stage-3 band (15.8-17.3 Å). PfCRT and PorF consolidate in
  stage 4 at 12.0-14.2 Å, and the two stage-4 stalls on record ran at 14 and 14.9 Å. The export gate
  consolidated in stage 3 at 8.5 Å at `3fb0b00b` and never at `603c4449`.
- The consolidation readout of section 6 separates the successes from the failures in all 46 logs, with
  narrow margins.
- The sticky subset is not searched again after stage 2 on the export gate and mostly not on PorF; it enters the
  final map with stage-2 poses.
- On the export gate the stage-3 chain starts from first-time particles, and revisited particles are as
  unstable as new ones.

Open, each with the phase that answers it:

- Does a finer early band recover PorF in ten of ten, and the export gate? (Phase 1)
- Does a 16 Å start cost the global search of PfCRT, or of a large particle? (Phases 1 and 5)
- Does the frozen sticky subset cost final resolution? (Phase 2)
- Does the export gate also need its sticky subset consolidated before the release, or a larger one? (Phase 3)
- Does the input change contribute to the export gate's failure? (Phase 1, control)

### 2.8 What the review changed

The independent review corrected the first version of this analysis; each point was checked against the
logs and the code:

- PfCRT has 16,832 active particles (16,778 at `db20cde7`) under full sampling (`nsample` 25,000), not
  about 6,700.
- The first export-gate table read each stage one iteration early; the table above uses each stage's last
  iteration.
- The completed export-gate restarts end at 7.17-7.40 Å; 7.89 Å was a partial restart.
- The trailing reconstruction blends decayed accumulators; it does not keep a pose until its particle is
  drawn again.
- Revisited particles do not hide a stable population, so changing the exit statistics is a diagnostic, not
  a remedy.
- "The PorF failure predates the recent changes" was too strong (section 2.2).
- SLC supports an experiment, not a 13 Å constant.
- New from the review, quantified here: the sticky subset leaves stage 2 with `updatecnt` 40 and is not drawn again
  (sections 2.4-2.5).

## 3. Contract

- No candidate changes the last rung (`LPSTOP_BOUNDS`, mask/26.7). The crops follow the end stages
  (`lpstages_setlims`: stages 1 and 8 cropped to their limits, a linear box ramp between), so they change
  only where stage 1's crop lies above the 88-pixel floor. For PfCRT, PorF and the export gate it is at the
  floor under every ladder here, so `lpstart=16.25` changes their bands alone, not their crops or shift
  limits; the stage-3 cap changes no crop anywhere.
- Full-sampling runs (PfCRT) have no sticky subset; the release fix and the late release leave them untouched.
- Multi-state runs, the class-average route (`solve3D_cavgs`), add-ons (manifest ladder) and the
  input-volume and no-2D routes keep their ladders and their sampling.
- Shared-memory and distributed runs behave the same: the release fix acts in the solve3D commander between
  stages, before any part is scheduled.
- No new command-line key, except `nsample_ramp` (Decision 15, the thing compared). Experimental variants are
  patch builds on an export (`git archive`) of the phase's base, recorded with their diff.
- `inpl_posterior` stayed as each target's command had it (default `no`; SLC's `cone`) until `fd9af0add`
  (2026-11-16) removed the 3D in-plane posterior insertion and its key; later runs drop it from SLC's command.

## 4. Decisions and rulings

1. (2026-11-15) Scope: single-state particle-route solve3D (`set_lplims_from_size` with `l_cavgs=.false.`,
   `nstates=1`) and its fractional sampling.
2. (2026-11-15) Two ladder candidates. The *start cap*: the size rule's first limit is at most 16.25 Å
   (`LPSTART_SIZE_BOUNDS(2)` from 30 to 16.25 Å, single state only). It equals `lpstart=16.25` for every
   mask above 130 Å (PfCRT and PorF: 16.25, 14.79, 13.32, 11.86, 10.39, 8.93, 7.46, 6.00 Å; the export gate:
   16.25, 14.76, 13.27, 11.78, 10.29, 8.79, 7.30, 5.81 Å) and leaves masks of 130 Å and less as they are
   (streptavidin, SLC). The *stage-3 cap*, the fallback: when stage 3's rung is coarser than 13.0 Å, stages
   3-8 are replanned linearly from 13.0 Å to the final limit and stages 1-2 keep the size rule (PfCRT: 20.00,
   18.00, 13.00, 11.60, 10.20, 8.80, 7.40, 6.00 Å). Phase 1's verdict chooses between them; Hans rules before
   Phase 4.
3. (Proposed 2026-11-15; ruled 2026-11-16 by Hans: the frozen sticky subset is unintended, put in the fix;
   implemented the same day, uncompiled) The *release fix*: before the first freely sampled stage after the
   sticky stages, the solve3D commander sets the `updatecnt` of every particle to 0
   (`release_sticky_subset` in `simple_solve3D_utils`, through `set_nonzero_updatecnt(0)`, an existing
   routine that had no caller), the way it already sets the starting sample's counts
   (`set_updatecnt(1, pinds)`). It acts on the transition from a stage that carries `sticky_sampling` to one
   that does not, which only the single-state fractional schedule has. The free sampling then starts from a
   clean history, as at the run's first iteration: the subset competes with the never-searched particles in
   the first sweep, the trailing chain represents only particles whose images it holds, and a particle not
   searched after the release is assigned by the final missing-update pass instead of entering the final map
   with a stage-2 pose. Resetting to 1 (never-searched particles first) was considered and not chosen: the
   chain would still count a subset it does not contain, and on the export gate, whose largest class needs
   38 iterations per sweep, part of the subset would still not be searched. No other workflow holds a
   sticky subset into a free stage (audit of 2026-11-16, Progress), so the fix stays in solve3D.
4. (2026-11-15) The *late release* (experiment, Phase 3): the sticky subset is kept through stage 3
   (`sticky_sampling=keep`) and the release, with the fix of Decision 3, moves to stage 4. A patch build. It
   repeats what `3fb0b00b` did by accident: a fixed subset consolidates before new particles arrive. (Amended
   2026-11-16: the patch leaves `TRAILREC_STAGE_SINGLE` at 3, since while the subset is sticky only its
   particles have counts, so u = n/N = 1 and the trailing blend is a no-op.) (Amended 2026-11-16 by Hans's
   ruling, after the control of Decision 13: the subset is held through stage 4 and free sampling starts at
   stage 5, as at `3fb0b00b`, whose greedy draw (`greedy_sampling=yes`, `STOCH_SAMPL_STAGE` 5) took the same
   best-scoring particles every iteration through stage 4. Implemented in the tree, not as a patch:
   `STICKY_LAST_STAGE` = 4 in the controller; the sticky stages are those from a run's first stage, when that is
   an shc stage, through stage 4; a run that starts later (supplied volumes, stage 4) holds none. Hans asked
   the same day what trailing from stage 3 gained over a sticky stage 4: nothing, u = 1 there, so the sticky
   stages emit `trail_rec=no` and trailing starts with the free sampling; `TRAILREC_STAGE_SINGLE` stays 3 for
   runs without a sticky subset.) (Amended again 2026-11-16, Hans's ruling (b): one boundary for every
   single-state route. `TRAILREC_STAGE_SINGLE = STICKY_LAST_STAGE + 1`; every single-state fractional particle
   run, add-ons excepted, holds the sticky subset from its first stage through stage 4, whatever its entry
   stage, so a stage-4 entry (supplied volumes, `cavg_ini`, `cavg_ini_ext`) holds it for that stage and trails
   from 5, as `3fb0b00b`'s greedy draw did below `STOCH_SAMPL_STAGE` on every route. Multi-state rules are
   untouched: no sticky subset, no scored draw, trailing from stage 8.)
5. (2026-11-15) The *larger sample* (experiment, Phase 3): `nsample=25000`, SLC's setting. On the export gate
   the class quota becomes 417 and the largest class is seen once in 16 iterations instead of 38.
6. (2026-11-15) Replayable restarts (Phase 2): solve3D fixes its run seed (`ensure_run_seed`) before its
   first random draw, so the starting volume is covered, and prints it; `restarted_exec` gives restart i the
   seed base + i - 1 when `SIMPLE_SEED` is set and prints it, and otherwise each restart prints the seed it
   drew. Restart i of two arms then starts from the same seed, and any restart can be replayed alone.
7. (2026-11-15) Diagnostics (Phase 2; printed only, no decision changes): overlap and angular change split
   into first-time and revisited particles, with their counts; one summary line at the end of each stage
   (stage, iterations, band, FSC=0.5 and FSC=0.143 at its first and last iteration, overlap, angular change,
   effective poses, first-time share, exit).
8. (Needs a ruling before Phase 2) Explicit `lp`: with an explicit `lp` the controller still passes `lpstop`
   at the planned rung (`emit_refine3D_stage_cfg`) and `set_bp_range3D` matches to the coarser of the two, so
   `lp` can coarsen a stage but never sharpen it, while the banner prints the command-line value and the
   policy says "an explicit `lp` wins". Proposed: the code follows the policy (with an explicit `lp`,
   `lpstop` = `lp` in the non-NU stages). Alternative: keep the code, correct the banner and the policy.
9. (2026-11-15) The floors and verdict rules of section 6 and the verdict tables of Phases 1 and 3 were
   frozen on 2026-11-15, before any `lpstart=16.25` result was read, and are never moved.
10. (2026-11-15) Restart counts: a decisive verdict needs ten completed restarts (PorF and PfCRT; the export
    gate in the panel); a screening arm four (the export gate) or five (PorF); fewer is a lead, not a
    finding.
11. (2026-11-15) One executable per comparison: Phase 1 runs at `603c44490`, the build of the current logs;
    later arms record their commit and patch, and arms compared within a phase share the build and frozen
    copies of the input projects (PorF's `3_solve3D` input, the export gate's `4_solve3D` input with 111,099
    particles, PfCRT's `4_solve3D` input).
12. (2026-11-15) Logs go to `~/Processing/solve3D_regression/<target>/<arm>/`, one folder per arm
    (`lpstart16`, `b1_lpstart16`, `b1_lpstart16_nsample25k`, `b2_lpstart16`, `candidate`), with the arm's
    command line in a `cmd.txt` beside the logs.
13. (2026-11-16) The *band-versus-sampling control* on the export gate (Phase 3): a random 20,000-particle
    subset of its input (`selection nran=20000 oritype=ptcl2D` on a copy of the project), run on B1 with
    `lpstart=16.25` and `nsample=25000`, so every particle is searched in every iteration (full sampling: no
    sticky subset, no first-time particles, no trailing blend), the regime in which PfCRT's 16,832 particles
    consolidate at 13 Å; two restarts (folder `b1_lpstart16_sub20k_full`). Read on the consolidation
    readout of section 6 at the end of stage 3 or 4; the final is reported but not judged against 4.5 Å, since
    a fifth of the particles cannot reach it. Consolidated: the band suffices and the large set's fractional
    regime is the fault, which the late release and the larger sample address. Not consolidated: a 13 Å stage
    3 is still too coarse for this target, and the next arm is a finer early band (re-plan, section 9 stop 3).
14. (2026-11-16, Hans's ruling) The *best-first sticky draw*: the sticky subset is drawn by the particles' 2D
    objective values, best first, as `3fb0b00b` did, and from the first free stage the objective plays no part.
    `sticky_sampling=best` (a new value of the internal key; no command-line key) draws at the stage's first
    iteration like `yes`, with the same quotas, but takes the particles with the best `ptcl2D` `corr`, whatever
    their update counts (amended 2026-11-16 with ruling (b): a stage-4 entry's starting sample has counts of 1,
    which a tiered draw would skip), under every balance mode: `none` the best active particles, `class` the
    best of every class, `cavg` the best of every unit, each as many as its quota inside the nested group quota. The
    project owner reads the scores from the `ptcl2D` field it holds; without them (no 2D solution) the draw is
    random and the log says so. The rescoring (`keep`) and every free draw are unchanged (no score). Single
    state only: multi-state solve3D holds no sticky subset (a greedy sample would keep a minority state's start
    labels, the reason `cf643a192` removed score ranking).
15. (2026-11-16, Hans's ruling) The *sample ramp*, under the switch `nsample_ramp=yes|no` (default `no`; the
    key is the thing being compared). Single state, fractional particle runs that hold the sticky subset: stage
    k of the sticky stages draws `nsample*k/(STICKY_LAST_STAGE+1)`, the free stages draw `nsample`, and the
    default `nsample` under the switch is 25,000 (`NSAMPLE_RAMP_DEFAULT`): 5k, 10k, 15k, 20k, then 25k. Every
    sticky stage redraws `best` at its first iteration at its own size; the fixed 2D ranking and quotas that grow
    with the target make each subset contain the previous one (a ±1 change at a unit's margin is possible only
    for quotas below the stage number). The stage's `update_frac` is the run's fraction scaled to its sample; a
    stage-4 entry draws its starting sample at the stage-4 size. The full-sampling switch stays on `nsample` for
    the whole run (PfCRT keeps full sampling, untouched). Refused with several states and without a 2D solution
    (the random fallback of the best draw would not keep the subset). Switch off: the run is the one of
    `218f0f289`. Experiment: the ramp against `nsample=25000` without it on the same build (the larger-sample
    arm of Decision 5, which isolates the ramp), after the candidate's current restarts are read.

Defaults taken without a question: 16.25 Å for the start cap (the size rule's own start at a 130 Å mask)
and 13.0 Å for the stage-3 cap (SLC's stage 3); the consolidation thresholds of section 6; the seed stride
of one; the folder names; the wording of the new log lines.

## 5. Not in this plan

- Multi-state ladders and sampling, `solve3D_cavgs`, `solve3D_addon`, `refine3D_states` and
  `classify3D_refs`.
- The stage-4 FSC=0.5 rule (`set_stage_lp_from_fsc05`) stays; it is revisited only if a stage-4 stall
  happens at a band coarser than the next rung, where the one-rung advance would have gone further.
- Exit decisions driven by first-time statistics (Decision 7 prints them only).
- A data-driven release or band progression ([data-driven stage progression](../rejected/data_driven_stage_progression.md)).
- The in-plane posterior insertion, the shape of the class quotas and the 2D side.
- Analysis scripts in the repository: the stage summary line replaces the external parser.

## 6. Acceptance

Final-resolution floors, the gates (FSC=0.143 of the final reconstruction, the `POSTPROCESS` low-pass value):

| Target | Floor | Origin |
| --- | --- | --- |
| PfCRT | at most 4.31 Å | measured: the worst of 15 completed restarts (`db20cde7` 10, `603c4449` 5) is 4.25 Å; one shell coarser at box 300 and 1.035 Å per pixel |
| PorF | at most 4.27 Å | measured: the worst of 14 completed successes (`fa9b82b5` 10, `603c4449` 4) is 4.22 Å; one shell coarser at box 240 and 1.3 Å per pixel |
| Export gate | at most 4.5 Å | guessed: between the current input's failures (7.17-7.40 Å) and the old input's successes (3.88-3.95 Å); no run on the current input has succeeded |
| SLC, streptavidin, beta-galactosidase, TRPV1 | the worst completed restart of their `603c4449` profile, plus one shell | measured in Phase 0 and frozen when collected |

Consolidation readout, reported for every restart, not a gate: overlap at least 0.30 and angular change at
most 30 degrees at the end of stage 3 or stage 4. Measured on the 46 logs of section 2: the successes reach
0.33-0.82 and 2.6-27.5 degrees, the failures 0.01-0.25 and 38-68 degrees, and the PfCRT stall of
`fa9b82b5` (overlap 0.36 at 31 degrees) misses on the angle; the thresholds sit in these gaps and are
guessed within them. A restart that passes its floor without consolidating, or consolidates and then
fails, is flagged.

Verdicts. An arm passes on a target when every completed restart meets the floor and at least the count of
Decision 10 completed. With PorF's present rate of one failure in six, ten of ten happens by chance with
probability (5/6)^10 = 0.16, so a PorF pass is read with the consolidation readout, and the panel's ten
restarts are the count Hans set. The export gate's two outcomes lie 3 Å apart, so four restarts per
screening arm separate them.

Reported, not gated: the stage-3 and stage-4 maps (`_stage03`, `_stage04`) of one restart per arm and every
final map, for Hans to judge (the export gate's helices after stages 3-4, map quality on every target);
Progress records their paths.

## 7. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | no code | the data preconditions |
| 1 | runs on the base | the verdict table of Phase 1 |
| 2 | seed and diagnostics: printing, and a seeded stream where an unseeded one was; the release fix and the explicit-`lp` ruling: result-changing | build, fast gate, the unit suites of the exit, the Phase 2 runs |
| 3 | patch builds and runs | the verdict table of Phase 3 |
| 4 | result-changing defaults | build, fast gate, unit suites, the highlevel entries that run solve3D |
| 5 | runs | the floors of section 6 |
| 6 | documents | links resolve; skills only on approval |

## 8. Phases

Run times at `603c44490` for planning: PfCRT 46-51 min, PorF 38-53 min, the export gate 71-81 min per restart.

**Phase 0, freeze and collect.** The floors and verdicts are frozen in this note. Collected into their
folders: the logs of the best SLC run (three restarts) with its particle count; a `603c4449` profile of
streptavidin, beta-galactosidase and TRPV1 (at least three restarts each, their usual command lines in
`cmd.txt`); frozen copies of the input projects (Decision 11); and, if it still exists, the export gate's
input of `3fb0b00b` (129,650 particles in 94 classes).
Exit: the control floors computed and written into section 6; every missing item listed in Progress as
"not judged".

**Phase 1, the band (running).** `lpstart=16.25` on PfCRT, PorF and the export gate at `603c44490`,
otherwise as the current runs. Control, if the old export-gate input exists: two restarts of `603c44490`
with its default ladder on that input; a success there makes the input a contributor and moves the sampling
arms of Phase 3 forward, a failure clears the input. Read per restart: the bands and, at the end of stages
2-5, FSC=0.5, overlap, angular change and effective poses; the consolidation readout; the final.

| PfCRT | PorF | Export gate | Reading | Next |
| --- | --- | --- | --- | --- |
| pass | pass | pass | the band is the lever on all three | Phase 2, then Phase 4 with the start cap; Phase 3 only if the panel's export gate fails |
| pass | pass | fail | the band suffices for the smaller sets, not for the large one | Phase 2, then Phase 3 on the start cap |
| fail | any | any | the 16 Å start costs PfCRT's global search | the stage-3 cap as a patch build on PfCRT (10), PorF (10) and the export gate (4), then as above |
| pass | fail | any | a PorF failure at a band of 11.9 Å or finer, finer than every successful stage 4, points away from the band | read its stages 3-4; Phase 3 adds PorF (5) to each arm |

Exit: the verdict and the per-restart table in Progress; the maps of one restart per target collected.

**Phase 2, replayable restarts, diagnostics and the release fix (build B1).** Implement Decisions 6 and 7,
and 8 as ruled, on the base plus the release fix of Decision 3 (written 2026-11-16, ahead of this phase, at
Hans's ruling). The release fix's test is in `simple_solve3D_reseed_tester` (`check_sticky_release`, beside
the sticky-subset handoff test): stages 1-3 run through the owner's sampler with the commander's boundary call
before every stage on a project file; it releases nothing while the stages are sticky, at stage 3 every
count becomes 0, in memory and in the file, and the first free sweep
draws every particle, while without the release (the negative control) that sweep skips the subset and takes
every other particle. Still to test here: the restart seeds; the explicit `lp`, beside the FSC=0.5
stage-limit tests. No reconstruction reads `sample4rec` between the reset and the release stage's first draw:
the stage starts with the sigma2 set-up, which reuses the committed state, and the draw comes first in its
first iteration. Runs at the Phase-1 ladder (`lpstart=16.25`,
or the stage-3 patch if Phase 1 chose it): PorF 5, two of them with one seed (the determinism check); the
export gate 4; PfCRT 2, the negative control without a sticky subset. Expected: PfCRT within its floor; PorF's finals
at least as fine as Phase 1's, perhaps a shell finer without the stale subset; the first-time split visible on
the export gate.
Exit: build and fast gate pass; `simple_test_exec test=unit_reconstruction` (the solve3D stage rules), and
`test=unit_ori` and `test=unit_cart_align3D` where their testers changed, pass; the stage summary and the split
appear in the logs; the determinism check recorded as identical, equal readouts or different; the runs read
against the floors.

**Phase 3, the large set's release (only if the export gate is not recovered).** On B1 at the chosen ladder,
four export-gate restarts per arm: the larger sample (Decision 5) and the late release (patch build B2,
Decision 4), with the control of Decision 13 alongside. An arm recovers when every completed restart meets the
export gate's floor.

| Larger sample | Late release | Next |
| --- | --- | --- |
| recovers | recovers | the late release: no extra cost per iteration |
| fails | recovers | the late release |
| recovers | fails | the larger sample, by ruling: `nsample` 25,000 by default or a rule on the particle count |
| fails | fails | both together (4); if that fails too, stop (section 9) |

The chosen arm then runs PorF (5), and SLC (3) if it changes SLC's run (the late release does; the larger
sample does not, SLC already uses it).
Exit: the verdict and the chosen arm's patch in Progress.

**Phase 4, the candidate build.** Implement on B1 the ruled ladder (Decision 2) and, if Phase 3 chose one,
the late release or the `nsample` rule. Update `simple_solve3D_reseed_tester` (the ladder at masks of 80,
130, 155, 160, 170, 180 and 200 Å; the cap applied to single state only; the class-average route unchanged)
and policy section 4.
Exit: build, fast gate and the unit suites pass, and so do the highlevel entries that run single-state
solve3D on masks above 130 Å: `simulated_workflow_6vxx` and `solve3D_addon` (180 Å) and the consensus route of
`multistate_workflows` (200 Å); `simulated_workflow_1jxy` was removed in `fd9af0add`.

**Phase 5, the panel.** On the candidate build, without `lpstart` on the command line: PorF 10, PfCRT 10, the
export gate 10, SLC 3 (its best command), streptavidin, beta-galactosidase and TRPV1 3 each (their usual
commands). Add a large-particle target if one is in Hans's regression set, since the start cap moves stage 1
of masks above 200 Å furthest (from up to 30 Å to 16.25 Å) and no target here tests that.
Exit: every completed restart meets its floor; the consolidation readout reported; the maps collected; Hans
rules on making the candidate the default.

**Phase 6, close.** Policies: `solve3D_policy.md` (section 4, the ladder; section 3, the release and, if
chosen, the late release; the explicit `lp`) and `importance_sampling_fractional_update_policy.md` ("Sticky
sampling schedules", the release). Skills `simple-solve3d-importance-sampling` and
`simple-frac-update-trailing`: edits proposed, applied on approval. The report
`solve3D_small_membrane_regression_report.md`; this note moves to `completed/`; the README row updated.
Exit: links resolve; `scripts/check_guide_imports.py` passes if a skill changed.

## 9. Stop conditions

1. A control (SLC, streptavidin, beta-galactosidase, TRPV1 or a highlevel solve3D entry) misses its floor or
   gate in any restart of a candidate: no default changes; diagnose first.
2. All three targets fail in Phase 1: the band is not the lever; re-plan.
3. Phase 3 recovers no arm, nor the two together: a 13 Å stage 3 does not suffice for the export gate;
   re-plan (a finer early band needs Decision 8's ruling or a patch).
4. A file outside the table needs an edit: stop, add its row with the reason, then edit.
5. A data precondition of Phase 0 is missing: that target is not judged, and the panel waits for it.

Not a stop: the determinism check finding that one seed gives two different runs; pairing is then a matched
start only, and Progress says so.

Gates no ruling may waive: the floors of section 6 as frozen; PorF ten of ten in the panel; the fast gate;
the unit tests added in Phases 2 and 4.

## 10. File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/solve/simple_solve3D_controller.f90` | the ladder constants and the stage-3 cap; the sticky and trailing stages of the late release (`STICKY_LAST_STAGE`, `holds_sticky_subset`, `solve3D_first_free_stage`); `sticky_sampling=best` at the first sticky stage; the sample ramp (`solve3D_stage_nsample`, `NSAMPLE_RAMP_DEFAULT`, `best` at every ramped sticky stage, Decision 15); the explicit-`lp` ceiling | 2 (Decision 8), 3, 4 |
| `src/main/solve/simple_solve3D_utils.f90` | `set_lplims_from_size` (the cap, single state); the release helper; the interfaces of `solve3D_first_free_stage`, `solve3D_stage_nsample` and `solve3D_nsample_default(params)`; `nsample_ramp` stripped from the child command lines | 2, 3, 4 |
| `src/main/commanders/simple/simple_commanders_solve3D.f90` | the release at the stage boundary; the run seed; the stage summary, if printed there; the unit table's planned draws counted from the first free stage; under `nsample_ramp` the default, the refusals, the stage-4 entry's starting sample and the per-stage log line | 2, 3 |
| `src/main/solve/simple_solve3D_reseed_tester.f90` | ladder, release, sticky-stage and best-first tests, under `none`, `class` and `cavg`; the ramp's emitted sizes and its growth (`check_ramp_growth`); explicit-`lp` tests | 2, 3, 4 |
| `src/main/params/simple_parameters.f90`, `src/main/params/simple_parameters_phases.f90` | `sticky_sampling=best` and `l_sticky_best` (Decision 14; added 2026-11-16 on Hans's ruling); `nsample_ramp` and `l_nsample_ramp` (Decision 15) | 3 |
| `src/main/params/simple_parameters_parse.f90`, `src/main/params/simple_parameters_tester.f90`, `src/main/ui/simple/simple_ui_solve3D.f90` | the `nsample_ramp` key: registry, binding test, solve3D input (Decision 15; added 2026-11-16) | 3 |
| `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90` | the owner's draw passes the `ptcl2D` scores under `best` (Decision 14; added 2026-11-16) | 3 |
| `src/main/simple_external_reference_pose_initialization.f90` | the pose initialization removes the stage's `sticky_sampling` from its own command line, since a stage-4 entry now carries `best` (Decision 4, ruling (b); added 2026-11-16) | 3 |
| `src/main/ori/simple_oris.f90` | the interfaces of the samplers' optional scores (Decision 14; added 2026-11-16); the comment of `class_sample_quotas` (square-root quota restored) | 3 |
| `src/main/simple_convergence.f90`, `src/main/simple_convergence_tester.f90` | the first-time split; the stage summary, if printed there | 2 |
| `src/main/exec/simple_exec_helpers.f90` | per-restart seeds in `restarted_exec` | 2 |
| `src/utils/math/simple_rnd.f90` | a getter of the run seed, if needed | 2 |
| `src/main/ori/simple_oris_sampling.f90`, `src/main/ori/simple_oris_tester.f90` | `set_nonzero_updatecnt` and its test, if needed; the best-first draw of `sample4update_cnt` and `sample4update_class` and its test; the square-root quota (`sqrt_capped_shares`, `alloc_unit_quotas`, `class_sample_quotas`) and its tests restored from `15d6280a9` (2026-11-16, Hans's request) | 2, 3 |
| `src/main/ori/simple_ori.f90`, `src/main/ori/simple_oris.f90`, `src/main/ori/simple_oris_getters.f90`, `src/main/ori/simple_ori_tester.f90`, `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90` | `has_search_score` (the search score alone) for the shc stages' dense first pass, so it reaches never-searched particles; `has_been_searched` unchanged; its test; the caller (added 2026-11-17 on Hans's ruling) | 3 |
| `doc/policies/3D/solve3D_policy.md` | sections 3-4 | 2, 3, 4, 6 |
| `doc/policies/importance_sampling_fractional_update_policy.md` | sticky sampling schedules; the sampling units (the one scored draw, the square-root quotas) | 2, 3, 6 |
| `doc/policies/3D/classify3D_refs_policy.md`, `doc/policies/heterogeneity/refine3D_states_policy.md`, `doc/algorithms/sampling_and_fractional_updates.md`, `doc/algorithms/solve3d.md` | the sticky schedule; the square-root quotas (added 2026-11-16 on Hans's request) | 3 |
| `.github/skills/simple-solve3d-importance-sampling/`, `.github/skills/simple-frac-update-trailing/` | on approval | 2 (the sampling skill, approved 2026-11-16), 3 (both, requested 2026-11-16), 6 |
| `doc/implementation_notes/README.md` | the topic row | 6 |
| this note | Progress | all |

## Progress

(One entry per phase: date, machine, build and patch, the per-restart table, the verdict, the map paths.)

**2026-11-16, Decision 3 ahead of Phase 2 (Mac session, base `603c44490`, uncompiled: no compiler in the
session).** Hans ruled the frozen sticky subset unintended and asked for the fix. Files:
`simple_solve3D_utils.f90` (`release_sticky_subset`: read `ptcl3D`, `set_nonzero_updatecnt(0)`, write the
segment, one log line), `simple_commanders_solve3D.f90` (the stage loop calls it at every stage boundary, before
the stage starts; it acts at the transition from a stage with `sticky_sampling` to one without), `simple_solve3D_reseed_tester.f90`
(`check_sticky_release`, which calls it before every stage, with the negative control without it,
from `test_gold_standard_stages`),
`solve3D_policy.md` section 3 and `importance_sampling_fractional_update_policy.md` "Sticky sampling
schedules". Owed: build, fast gate and `simple_test_exec test=unit_reconstruction`. The audit for the same
fault elsewhere found none: multi-state solve3D runs `shc` in stages 1-2 and samples freely from stage 1, with
no sticky subset; `refine3D_states` holds a sticky subset in every frequency block for at most three
iterations (`FREQUENCY_BLOCK_NITS`) and draws the next block's subset lowest count first, so the counts grow
evenly and nothing is frozen out; `refine3D_states_cavgs` updates fully (it refuses `nsample`, `update_frac`
and `sticky_sampling`); `classify3D_refs` sets no sticky subset; solve2D's
sticky stage reproduces its subset without advancing `updatecnt` (`sample4update_reprod`), so the subset leaves
it at a count of 1; `solve3D_addon` runs only the base run's last stage. The fix covers every caller of the
solve3D stage loop: the first chunk of `solve3D_chunks`, the pose student's and SINGLE's solve3D runs. Terms:
the first version of this note called the sticky subset the cohort, the name the library gave up in
`0fb4a7494`; the note now uses the library's name. On Hans's approval the same day the skill
`simple-solve3d-importance-sampling` describes the release (`SKILL.md` working rules and focus list,
`references/coupling-map.md` sampling by stage, bookkeeping and the owner's draw); `simple-frac-update-trailing`
was left as it is, its statements still true; `scripts/check_guide_imports.py` passes.

**2026-11-16, Phase 1, export gate (partial).** `lpstart=16.25` did not recover the export gate (Hans; the
logs, and the PorF and PfCRT arms, are not yet in `~/Processing/solve3D_regression`, so the per-restart table
and the ladder verdict are still owed). By the verdict table this opens Phase 3 for the export gate. Prepared,
not applied: the late-release patch for B2, `late_release_experiment.diff` (`simple_solve3D_controller.f90`
only, on top of B1; `patch --dry-run -p1` applies cleanly), which keeps the sticky subset through stage 3 for a
single-state run that began with the sticky stages (`STICKY_LAST_STAGE_LATE`). On a B2 build two assertions of
`simple_solve3D_reseed_tester` fail by design (stage 3 is sticky). Added Decision 13, the band-versus-sampling
control, frozen before its run. Expected for B1 alone on the export gate (the Phase 2 arm): no recovery, since
the release fix changes which particles the free stages draw and what the final map holds, not whether stage 3
consolidates.

**2026-11-16, Phase 1 verdict (as Hans reported it; the logs are on another machine, the per-restart table is
owed).** `lpstart=16.25` at `603c44490`: half of the PorF restarts fail (pooled baseline 1 of 16), PfCRT fails
3 of 10 (baseline 0 of 15), the export gate is unchanged (every restart fails). Verdict-table row 3 (PfCRT
fails) and stop condition 2 (all three fail) both apply; the start cap is rejected. The stage-3 cap was not run:
from stage 3 on it differs from `lpstart=16.25` by at most 0.3 Å.

**2026-11-16, Decision 13 control (as Hans reported it; readout owed).** The random 20,000-particle subset of
the export gate at `nsample=25000` (full sampling) and `lpstart=16.25`: the map recovers, and the stage-3
volume shows helices, though the run gets there later than `3fb0b00b` did (Hans's reading). Read: with every
particle refined every iteration the export gate consolidates, so the full set's fractional regime is the
fault. Hans ruled the late release through stage 4 (Decision 4, amended) and the best-first sticky draw
(Decision 14).

**2026-11-16, Decisions 4 (amended) and 14 implemented (session on the Mac's folders; uncompiled, no compiler
in the session).** The controller holds the sticky subset from a run's first stage, when that is an shc stage,
through `STICKY_LAST_STAGE` = 4 (`holds_sticky_subset`, single state, particles, fractional updates) and emits
`sticky_sampling=best` at the first sticky stage, `keep` after; `release_sticky_subset` therefore acts at stage
5 unchanged. `TRAILREC_STAGE_SINGLE` stays 3 (u = 1 while the subset is held); its comment says so. The
samplers `sample4update_cnt` and `sample4update_class` take optional per-row scores and then order a tier best
first (`order_best_first`) instead of shuffling it; the quotas are untouched. `sample_ptcls4update3D` passes
`ptcl2D` `corr` under `best` (`sticky_draw_scores`, with a log line either way). The commander's unit table counts
its planned draws from the first free stage (`solve3D_first_free_stage`). Tests: `test_best_first_draw`
(`simple_oris_tester`: `none` with a partial and a whole lower tier, `class`, `cavg`); in
`simple_solve3D_reseed_tester` the stage table (stage 1 `best`, 2-4 `keep`, free from 5), a stage-2 start held
to stage 4, a stage-4 start and the class-average route with none, `check_early_stage_sticky` through stage 4
(random fallback without `ptcl2D`), and `check_sticky_release` over stages 1-5 under `none`, `class` and `cavg`,
asserting the best rows of every unit's quota, no release before stage 5, the release at 5 and the first free
sweep (the negative control unchanged). Policies: `solve3D_policy.md` section 3, the fractional-update policy
(the `best` rule, the sampling units' one scored draw, sticky sampling schedules). Owed: build, fast gate,
`simple_test_exec test=unit_ori` and `test=unit_reconstruction`; the skill `simple-solve3d-importance-sampling`
(its working rule says the sampler never uses objective values) on approval. Proposed runs, at the default
ladder: the export gate 4, PorF 5 (it holds a sticky subset), PfCRT 2 (full sampling, no sticky subset: the
negative control); read the log line `STICKY SUBSET DRAWN BY THE 2D OBJECTIVE`, the release line before stage
5, and the consolidation readout at the end of stages 3 and 4.

**2026-11-16, follow-ups (same session, uncompiled).** Hans's style rulings applied: declaration blocks aligned in
columns where `scores` was added (interface and implementation), `present(scores)` read once into `l_scores`
before the loops. The sticky stages emit `trail_rec=no` (Decision 4's amendment); the reseed tester's stage table
expects trailing from stage 5, and a stage-4 start trailing from its first stage. Skills updated on Hans's
request: `simple-solve3d-importance-sampling` (`SKILL.md`, `references/coupling-map.md`) and
`simple-frac-update-trailing` (`SKILL.md`, `references/frac-update-contract.md`). Policies made consistent with
the code: besides the sticky schedule, the class-unit quotas are equal shares capped at the populations
(`group_quotas`, `alloc_unit_quotas` since `a027fe3bd`), where `solve3D_policy.md`,
`importance_sampling_fractional_update_policy.md`, `classify3D_refs_policy.md`, `refine3D_states_policy.md`
and `doc/algorithms/sampling_and_fractional_updates.md` still said square-root shares;
`doc/algorithms/solve3d.md` updated for the schedule. Open question from Hans: release earlier to trail sooner
(more signal) or keep the hold through stage 4; the runs' readout at the end of stage 3 decides whether
`STICKY_LAST_STAGE` = 3 is safe.

**2026-11-16, ruling (b) implemented (uncompiled).** `TRAILREC_STAGE_SINGLE = STICKY_LAST_STAGE + 1`, the
`trail_rec=no` override in `emit_refine3D_stage_cfg` removed. `holds_sticky_subset` tests the stage alone (one
state, particles, fractional, not an add-on, first stage through 4; `l_addon` passed from `set_cline_refine3D`
through `emit_refine3D_stage_cfg`, and to `solve3D_first_free_stage` from the commander). The `best` draw ranks
by score alone: `sample4update_cnt` takes the best-scoring active rows, `sample4update_class` each unit's
best-scoring particles, quotas unchanged; the tiered keyed draw without scores is unchanged. The external-reference
pose initialization deletes `sticky_sampling` from its command line, which it copies from the stage-4 line. Checked
for copies of the stage line that would inherit `keep`/`best`: the missing-update pass already deletes it, the
sigma2 commit and the multi-state final assignment do not sample by it (and multi-state lines never carry it).
Multi-state audit: `holds_sticky_subset` returns false for `nstates /= 1`, `TRAILREC_STAGE_MULTI` stays 8,
`refine3D_states` keeps `sticky_sampling=yes` (no score). Tests: `test_best_first_draw` expects the best rows
whatever their counts (an updated row taken, never-updated low scorers not preferred); the reseed tester expects a
`cavg_ini` and a stage-4 entry to draw `best` and not trail in stage 4, to trail from 5, a stage-5 entry to trail
from its first stage, an add-on stage to hold no subset, full sampling to hold none and not trail;
`check_early_stage_sticky(params, 4, 4)` added. Policies (`solve3D_policy.md` section 3, the fractional-update
policy), `doc/algorithms/{sampling_and_fractional_updates,solve3d}.md` and both skills updated.

**2026-11-16, square-root quotas restored (after `8363af240`; uncompiled).** The follow-up entry above was
wrong about the quotas: I took the code for the rule and rewrote the documents to equal shares. Hans recalled
pushing square-root shares. The history: `15d6280a9` (2026-11-05 17:15) introduced them (`sqrt_capped_shares`,
the floors of the real quotas plus the remainder to the least-updated units with a fractional quota, tests
asserting 2, 3, 4 and 1, 3, 3, 3, 9); `a027fe3bd` (2026-11-06 01:01, the agent-friendly run, whose plan was
written 2026-11-04) replaced them with equal rounds in `simple_oris_sampling.f90`, the tester, the comment in
`simple_oris.f90` and both skill references, without a word in its message or report and without touching the
policies; the rest of `15d6280a9` (the class default, the cavg group floor, the single-state missing pass)
survived. The run evidently worked from a base older than `15d6280a9`. Every build since, the last-working
`3fb0b00b`, `fa9b82b5` and `db20cde7` included, drew with equal shares, so the revert does not separate the
working from the failing builds; it does shape the best-first sticky subset (equal shares give a small odd class
the same quota as a large good one). Restored from `15d6280a9`: `alloc_unit_quotas`, `sqrt_capped_shares` and
`class_sample_quotas` (replacing `group_quotas`), the `sample4update_class` and `class_sample_quotas` comments,
the assertions of `test_nested_class_quota` and `test_class_sample_sweep`. Today's tests recomputed (checked with
a small model of the rule that reproduces `15d6280a9`'s expected values): under `cavg` with units of 4 | 4, 8 and
4 of 16, the quotas are 1.46 | 1.05, 1.49, so 1, 1 and 2 particles (the remainder to the larger of the tied
least-updated units), and the best rows are 3, 6, 10, 13 in `test_best_first_draw` and `check_sticky_release`;
the `class` cases (equal populations) and the sweep counts are unchanged. The documents rewritten in the entry
above are back to the square-root rule, with the revert and the restore recorded in `solve3D_policy.md` and the
fractional-update policy.

**2026-11-16, first reading of the candidate (Hans; restarts still running, logs not read here).** With the
best-first sticky subset held through stage 4 (`8363af240`, square-root quotas restored in `218f0f289`), "everything
seems to be working": the targets look promising part-way through. Owed when the restarts complete: the build
each arm ran (`8363af240` alone drew with equal shares, `218f0f289` with square-root shares), the per-restart table
against the floors of section 6, the consolidation readout at the end of stages 3-4, the stage-3/4 maps of one
restart per target, and the controls of section 9 (SLC, streptavidin, beta-galactosidase, TRPV1, the highlevel
solve3D entries) before any default is declared.

**2026-11-16, Decision 15 implemented under its switch (uncompiled; base `218f0f289`).** `nsample_ramp=yes|no`
(default `no`) in `parameters`, the registry and the solve3D UI (advanced, next to `nsample`), stripped from every
child command line with the sampling keys. The controller's `solve3D_stage_nsample` is the one rule (nsample, or
`nsample*istage/5` in the stages that hold the sticky subset); `emit_refine3D_stage_cfg` emits it as the stage's
`nsample` and scales the run's `update_frac` to it, and under the switch emits `best` at every sticky stage instead
of `keep`. `solve3D_nsample_default(params)` returns 25,000 under the switch. The commander refuses the switch with
several states and without a 2D solution, logs the ramp (or that full sampling leaves none), draws a stage-4
entry's starting sample at its stage's size, and logs `STAGE n DRAWS m PARTICLES PER ITERATION`. Switch off, the
emitted lines are those of `218f0f289`. Tests: the stage table under the switch (200, 400, 600, 800 of 1000 and
the scaled fractions, `best` at every sticky stage, 1000 and the run's fraction from stage 5, a stage-4 entry at 800,
several states and full sampling unramped, switch off `keep`); `check_ramp_growth` under `none`, `class` and
`cavg` (40 rows, nsample 20: sizes 4, 8, 12, 16, each stage's subset containing the previous through every iteration;
the `cavg` layout gives a unit nothing at stage 1, checked with the quota model to stay nested);
`test_nsample_ramp_binding`. Policies, the sampling skill (`SKILL.md`, `coupling-map.md`) updated; the solve3D
policy's stale "`cavg` by default" corrected to `class`. Owed: build, fast gate, `unit_reconstruction`,
`unit_ori`, the parameters suite; the experiment of Decision 15 after the current restarts are read.

**2026-11-17, the dense first pass reaches never-searched particles (uncompiled; base `76110d3f5`).** Asked what a
particle the ramp adds (updatecnt 0) meets, I traced it: the shc neighbourhood search forces its dense first pass
(every reference scored into the table, the draw over the posterior of all of them) only when
`has_been_searched` is false, and that function counted Euler angles as evidence of a search. solve3D gives every
particle a random start orientation (`delete_3Dalignment`, then `rnd_oris`), so no particle ever took the full
scan: an added particle's first pose was drawn from the handful of references that first beat its random pose,
in effect a random pose inserted into the map. The policies' rationale for the shc stages ("scans a never-searched
particle in full", written with `cf643a192`) never held on this route. Hans ruled the fix unconditional. The dense
pass now keys on a new predicate, `has_search_score` (the search score `corr` alone, which only a search writes).
Reach: stage 1's first iteration on every single-state run (the dense draw over 500 references against the noise
start instead of shc from the random pose), the particles the ramp adds in stage 2, and a
`prob_neigh_mode=shc` refine3D on a never-searched or imported `ptcl3D` without `corr`; `prob` stages, multi-state
(`shc`) and `refine3D_states` are untouched. Test `test_has_search_score` (`unit_ori`'s ori suite): a fresh record
and a random-start record with kept shifts carry no search score, a search writes one, the solve3D restart clears
it, and `has_been_searched` still counts the random angles as a pose. Policies and the sampling skill reference say
what now happens. The current candidate's runs predate the fix: stage 1's first iteration differs on any rerun.

**2026-11-17, correction (Hans's fast gate).** The first version changed `has_been_searched` itself, which I had
wrongly taken for the dense pass's only caller. It is also read by `is_virgin_field` (and through it the 3D
matcher's strategy choice, the reference preparation, refine3D's random start and the RELION export), by the 2D
search and by the 2D matcher. `unit_project` caught it: `test_relion_export_poses` sets Euler angles without `corr`,
the field became virgin and the export wrote no 3D poses (7 failed checks). `has_been_searched` is restored; the
dense pass reads `has_search_score`, so nothing else changes.

**2026-11-17, Decision 15 withdrawn (Hans's ruling).** The sample ramp did not speed the runs up; `nsample_ramp`
and its code are removed (the `src` side of `76110d3f5` reversed: the key, `l_nsample_ramp`,
`solve3D_stage_nsample`, `NSAMPLE_RAMP_DEFAULT`, the ramp's refusals and log lines, its tests). The sticky stages
are those of `218f0f289` again: stage 1 draws `best`, stages 2-4 `keep`, every stage at `nsample`. The
`has_search_score` dense pass of the entry above stays; it serves stage 1's first iteration.

**2026-11-18, square-root quotas and the dense first pass reverted (Hans's ruling after the gate runs; base
`a8a59331`).** On the same particles and command, today's build against the Nov16 build (`8363af24`), 10 restarts
each: PorF broke its pseudo-symmetry in stage 5 later (never in 4 restarts of 10, in 0 of 10 before; the sample's
view-imbalance correction 5.0x against 6.8x), FlhB likewise (p = 0.001 on the break iteration), and PfCRT, full
sampling with the quotas inactive, was slower in stage 1 (11.6 against 7.9 iterations a shell behind, p = 0.03; one
restart never recovered). The net source diff between the builds leaves two behavioural changes on solve3D's
default path: `218f0f289` (square-root quotas, the fractional sample) and `e94681175` (the dense first pass, stage
1 iteration 1); everything else is net-zero there (the in-plane insertion's removal restores the `no` branch,
`keep_centring` and the refine3D_auto defaults are the workflow's own, the shift solver's default is unchanged).
Both reversed: equal shares again (the code, tests and documents of `218f0f289`); the dense pass keys on
`has_been_searched` again and `has_search_score` is gone. The pass is not needed: the sticky subset is searched at
every shc iteration, so its particles climb from their random start and refine from there. A sampling policy that
adds particles during the shc stages must bring its own first-pass rule (the solve3D policy, the fractional-update
policy and the sampling skill's coupling map say so).
