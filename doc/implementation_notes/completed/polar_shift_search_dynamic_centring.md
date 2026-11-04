# Dynamic centring: the disk per particle, across the stages

Proposal, 2026-11-02 (Claude, for Hans's ruling before implementation). Companion
to `polar_shift_search.md`, whose rulings 15-18 and Phase 4 this builds on.
Completed with the plan on 2026-11-03. Ruled and built the same day (ruling 19; the as-built record is the note's
progress entry of 2026-11-02): the two metrics of section 2.6 with the rule of
section 2.1, the box's mask per disk particle (B), the subspace pool in the
state mode and the geometric neighbourhood in the geom mode, `trs_window` and
the drain exits retired, the 2D particle disk cache removed on the way.

## 1. The oversight

The centring is a stage: solve2D's stage 2 and solve3D's stage 3 search the
disk of 3 sigma for every sampled particle, drain the particles into the next
stage's window, and end when few still need the disk (rulings 16, S1, S3). Every
stage after searches the window from the coarse window map (Phase 4). The
design assumed that the particles the centring stage sees are the particles the
workflow has, and the beta-gal bench, which samples all 4463 particles every
iteration (`NSAMPLE_DEFAULT_2D` = 200000, `nsample` above the active count),
never contradicted it.

The workflows sample. solve3D draws `update_frac` of the particles per
iteration from the lowest `updatecnt` tiers (`sample4update_cnt`, class or
class-average balanced under `balance`), fills in every fifth iteration
(`sample_ptcls4fillin`), and its default `nsample` is 10000 against
`UPDATE_FRAC_MAX` = 0.9 and the full-sampling switch at 0.9 of the active
count. solve2D holds a sticky random subset through stage 1 and resamples by
`updatecnt` from stage 2 (`sample_ptcls4update2D`). The stream's pool receives
new particles every iteration and clears their `updatecnt` and `sampled` on
injection (`simple_stream_pool2D`, decision 5). So:

- a particle first sampled after the centring stage ended never sees the disk.
  With a centring stage of k iterations at update fraction f, the share of the
  particles it never sampled is about (1 - f)^k: at f = 0.3 and k = 4, 24 %;
  at f = 0.1, 66 %. Those particles carry their picking error into the window
  stages, which cannot recover a 20-30 A offset (the finding that motivated the
  work);
- a particle whose later move exceeds its window is never centred again: the
  drain's second rule (`shincr` above the window) is read only in the centring
  stage;
- the stream has no centring stage at all: the pool's refine2D runs the window
  from the first iteration, under the default key, with the window map.

The centring belongs to the particle, not to the stage.

## 2. The proposal

Replace the centring stage by a per-particle state that every table fill reads
and writes, in 2D and 3D, in the dense, sparse and neighbourhood tables alike.

### 2.1 The state (ruled 2026-11-02: both metrics, section 2.6)

Each searched particle carries two numbers, written by the consumer (the 2D and
3D probabilistic strategies) after the polish:

- `shincr`, the last move: the magnitude of the polished increment, in pixels;
- `shsig`, the particle's prior width for its next search, in pixels: the
  posterior width of section 2.6, inflated to the radius searched when the move
  reached that radius's bound. Below zero is the sentinel, *no history in this
  workflow*.

A sampled particle's next box is the smallest rung of the ladder of radii
(section 2.6) that covers three of its prior widths, floored at the stage's
`trs` and capped at the disk's 3 sigma; the sentinel takes the disk. With the
two rungs the fills have today, the window and the disk, the rule reads: the
disk when `shsig` is the sentinel or exceeds a third of the stage's `trs`, the
window otherwise. The fill writes nothing.

The sentinel is set once per workflow, at the start of solve2D and solve3D (the
commanders, in place of the refine2D/refine3D masters' per-run reset: ruling 1,
per workflow; the 2D shifts stay the 3D's start, the 3D disk is searched about
them under the pick prior), and on every injection: `import_particles` for a
new project, the stream's pool on the chunk's arrival (beside its
`clean_entry('updatecnt', 'sampled')`); the fill-in and missing-update passes
need nothing since their particles carry the state from the workflow's start. A
standalone refine2D or refine3D sets it at `startit = 1` only.

The alternative sentinel, `updatecnt == 0` at the fill, needs no new write but
fails in the stream (the pool clears `updatecnt` on every injection, for the
new and the old alike) and in a continued refine3D (the counts persist while
the centring state should not). The particle's own record is the cleaner
carrier.

### 2.2 The disk in every fill

The dense fills have both paths since Phase 4 (`process_particle_with_map` with
`l_disk`). Two changes make the disk independent of the stage:

- the disk's radius is 3 sigma in the stage's pixels whatever the stage's
  `trs`, bounded by the box (today it is `min(3 sigma, trs)`, which is why the
  controllers raise the centring stage's `trs`). The polish of a disk particle
  runs in a second search object per thread with limits of 3 sigma (the pair
  the consolidation had before Phase 4 and for the same reason: the box polish
  clamps to the object's limits); the window object keeps `trs`;
- the drained particles' window is the stage's `trs`, so the window map's radius
  is `trs` in every stage and the second map geometry of the centring stage
  (`trs_window`) goes.

The sparse 2D fill (`prob_neigh`, solve2D's stages 3-6 and the stream's pool)
scores a bounded, extremal-decayed set of classes at the seed and refines the
best of them. A disk particle there is a particle without a trustworthy class:
it scores every active class at the origin (one dense row, 90 posterior values
at 5.6 us each), screens the pool holding 99 % of that posterior with the disk
map, and polishes the best `smpl_ncls` of the screen; the window particles keep
the sparse path. The fill's candidate buffer takes the extra rows unchanged.

The neighbourhood table (solve3D's stages 6-8, `prob_neigh` in its modes)
scores the neighbourhood of the particle's stored direction. A disk particle
there has a direction as stale as its shift: an injected particle's direction is
the random start's, a mover's may be a neighbour's. Its disk is the subspace's:
the coarse scoring over the `nspace_sub` representatives at the origin that the
`state` and `geom` modes already do (`find_peak_subspaces`,
`build_pooled_neighborhood`), the disk map over the pooled neighbourhood of the
best representatives, and the polish of the best of the screen; in the `shc`
and `snhc` modes (solve3D's stages 1-2, no shifts) nothing changes. This is the
one place the proposal widens a search beyond the shift: a particle that needs
the disk also gets the subspace's direction search, which is the right thing for
an injected particle and a small cost for a mover.

### 2.3 What goes

- the centring stage as a stage: `CENTRING_STAGE` and `solve2D_centring_stage`
  in 2D, the stage-3 raise of `trs` and `cfg%trs_window` in 3D, the `trs_window`
  parameter and its UI entries;
- the drain as a stage exit: the `shift_stable_its` blocks of `check_conv2D` and
  `check_conv3D`, `DISK_NEED_FRAC`, the tester cases `test_shift_map_exit` and
  `test_shift_map_exit3D`; the stages end on their own rules again. The master
  keeps printing the disk set it is about to search, split into the never
  centred and the movers, so a run shows its injection rate;
- the masters' per-run reset of `shincr` (the commanders set the sentinel once).

### 2.4 The mask of a disk particle (ruled 2026-11-02: B)

The centring stage runs under the box's mask because a particle 20-30 A off its
pick is cut by the molecular mask at the pick. Ruled: a particle that needs the
disk is prepared under the box's mask whatever the stage. The matcher's particle
preparation (`prepimg4align` in `simple_matcher_2Dprep` and the 3D batch) takes
a per-particle mask radius: the box's for a particle whose state says the disk
(the sentinel or a width above the window's third), the stage's for the others.
The read-side branch is the only place the state is read outside the fills; the
fill and the preparation apply the same rule to the same field, so they agree.
stage 2 (2D) keeps the box's mask as a stage policy for the mask estimate at
its exit (`mskdiam=0`), which is a mask question, not a centring one.

### 2.5 Costs

From the beta-gal runs of 2026-11-02 (4 parts of 1116 particles, 6 threads each,
the fills' wall time per part and iteration):

| path | 2D (90 classes) | 3D (solve3D stage 3) |
|---|---|---|
| the disk (screen + polish) | 5.1 ms per particle | 14 ms per particle |
| the window map | 3.6 ms (dense), 0.5 ms (sparse) | 11 ms |
| the window search without the key | 1.2 ms (dense), 0.4 ms (sparse) | 8 ms |

The disk is 1.4 times the window map in 2D and 1.3 times in 3D, so moving a
particle's disk from the centring stage to the iteration it is first sampled
changes the total work little: every particle is centred once, the movers again,
and the work follows the sampling instead of the stage. Under fractional
sampling the never-centred pool drains geometrically with the sampling, so the
first stage with shifts costs a disk iteration spread over its iterations rather
than one expensive first iteration. The sparse 2D stages are the one place the
disk is far above the stage's own cost (5 ms against 0.5 ms a particle): a
stream pool with 10 % new particles per iteration would spend half its fill on
them, which is the price of centring them at all.

### 2.6 The two metrics and one rule (ruled 2026-11-02)

Ruled: record the last move and the posterior width, and combine them. The
combination, and why it is the sensible one:

**The width.** The coarse map's draw already evaluates the posterior over every
cell of the box searched (the shift cells of the map at every rotation). The
posterior width `s` is the standard deviation of the shift under that posterior
at the physical temperature (`prob_beta` = 1: a measurement of uncertainty, not
the exploration-tempered draw), in pixels: the square root of the trace of the
covariance of the cell shifts in the particle's axes, from the first and second
moments the draw's loop accumulates (the mean needs each cell's shift rotated
into the particle's axes; the second moment is rotation-invariant). About the
mean, not the mode, so a two-peaked posterior is as wide as it is uncertain.
Under the disk the weights carry the pick prior with the data term, as the draw
does; under the window the data term alone. Measured on the coarse band, `s` is
an upper bound of the full-band width: conservative, more disks rather than
fewer.

**The move.** The width is the posterior's within the box that was searched. A
polish whose move `d` reached that box's bound (within half a map step of the
radius `R` searched) says the posterior was truncated there: the mass beyond
the edge was never seen, and `s` understates the uncertainty by an unknown
amount. The one honest statement is that the uncertainty is at least the radius
searched. So the particle's prior width for its next search is

    sigma_next = s                 if d < R - step/2
    sigma_next = max(s, R)         if d >= R - step/2

and `shsig` stores `sigma_next`. The consumer computes it, since it alone knows
`s`, `d` and `R` together; the fill reads one number.

**The box.** The next search covers three prior widths, the confidence the disk
itself uses (3 sigma of the pick prior): the smallest rung of the ladder of
radii at or above `3 sigma_next`, floored at the stage's `trs` (no box smaller
than the window: a settled particle searches the window as today, and no
particle can lock itself into a box its posterior cannot leave) and capped at
the disk. With two rungs, the window and the disk, a particle takes the disk
when it was never searched, or when its last move reached its box's bound (the
rule in force today, which this reproduces exactly), or when its posterior is
wider than a third of the window although it did not move; and the window
otherwise. A third rung (the geometric mean of the two, say) is one more map
geometry when the widths on data show particles between.

**The pick prior at the first centring only (as built, 2026-11-02).** The
first fixture run re-centred correctly placed particles under the pick prior
and lost 56 points: the draw tempers the data term by `prob_beta` (0.3) and
not the prior, so at a re-centring the prior weighed more than three times its
physical weight and pulled the particle back towards its pick. A re-centring
is now the data's alone; the pick prior regularises a particle's first
centring, where the references are poorest and the pick is the only
information. The bound test is a tenth of the radius, or half a map step where
the map is finer, since a window map of two rings has a step equal to its
radius.

**The width as measured, and the Laplace width (2026-11-02, after Hans's beta-gal run).**
On data the coarse map's width reads 0.00 for nearly every particle from the
second iteration on, with 90th percentiles of 0.27 px early and 0.00 later. That
is a resolution limit of the measurement: the coarse band ends at shell
`kscreen`, so the map's radial step is `box / (4 kscreen)`, 3.6 px at beta-gal's
stage 4 against a window of 4.4 px (the window map has two rings, the origin
and the rim), and at the physical temperature the loss differences between
cells are tens of nats, so all the mass sits on one cell. A read of 0.00 means
"surer than one cell", which is the truth for a well aligned particle, and the
rule "three widths beyond the window" needs widths resolved below 1.5 px, which
the map cannot give. As it stands the width is a two-state flag, multi-modal at
the coarse band or not, and the re-centring is carried by the move at the
bound, which drains to a few movers an iteration on data. Two changes within
the map do not help: the draw's temperature (0.3 times the physical) inflates
the widths with the exploration's flatness, not the particle's uncertainty;
a floor at the cell's size makes every particle wide against a window of one
cell.

The informative version, the **Laplace width** at the polished optimum on the
full band (the loss's curvature at the committed pose), with the review's
qualifications, is its own planned note since 2026-11-03:
`../planned/polar_shift_laplace_width.md`. Not built; Hans's call.

**What the width does not do.** It never enters the likelihood as a prior. The
posterior re-used as a prior on the same particle image would count the data
twice; the pick prior stays the only prior in the disk's draw and in the 3D
insertion's loss, and the particle's width chooses the domain of the search,
not the weight of its cells. This is what keeps the rule mathematically clean:
coverage of the posterior at three standard deviations, the truncation
correction from the move, and the prior untouched.

**The record.** `shincr` keeps the move for the diagnostics; `shsig` holds
`sigma_next` (a new slot beside `shincr`). The master prints, per iteration, the
disk set by cause (never searched, bound reached, wide) and the median and 90th
percentile of the widths, the retained-mass measurement the review asked for.
A class change (2D) or a direction change beyond the neighbourhood (3D) needs no
guard of its own: the polish at the new class or direction measures a new width
and move, and a wrong-class centring shows as a bound hit or a wide posterior at
the next search.

## 3. Interfaces that change

- `shincr` and `shsig` (a new slot): the sentinel's writers (the solve2D and
  solve3D commanders, `import_particles`, the stream's pool injection,
  refine2D/refine3D at `startit = 1`); the width's accumulation in
  `shiftmap_draw` (the moments at the physical temperature, returned beside the
  draw); the consumer's `sigma_next`; the readers (the four fills, the particle
  preparation's mask, the prob_align drivers' print).
- the particle preparation (`prepimg4align`, 2D and 3D): the box's mask for a
  particle whose state says the disk.
- `simple_eul_prob_tab2D`: the sparse fill's disk path; the dense fill's second
  search object and the disk radius.
- `simple_eul_prob_tab`: the second search object and the disk radius.
- `simple_eul_prob_tab_neigh`: the disk path over the subspace pool.
- `simple_polarft_shiftmap`: `setup_shiftmaps` with the disk's radius 3 sigma
  bounded by the box and the window's `trs` (no `trs_window`).
- `simple_solve2D_controller`, `simple_solve3D_controller`: the centring stage
  retired; stage 2's box mask kept in 2D for the mask estimate.
- `simple_convergence`, `_tester`: the drain exits retired.
- `simple_parameters`, parse, phases, the UIs: `trs_window` retired.
- `simple_commanders_solve2D`: the adoption of the mask estimate at "the
  centring stage's exit" becomes stage 2's exit.
- the policies (`solve2D_policy.md`, `refine3D_policy.md`, `solve3D_policy.md`,
  the stream's pool policy) and the note.

## 4. Implementation plan

1. The state and the dense fills: the sentinel's writers, the second search
   object and the disk radius in both dense fills, the print; `trs_window` kept
   for one build so the fixtures compare the stage and the dynamic forms.
2. The sparse 2D fill's disk path; the neighbourhood table's disk path over the
   subspace pool.
3. The width in the draw and the consumer's `sigma_next`; the preparation's
   per-particle mask; the controllers: the centring stage retired, the drain
   exits retired, `trs_window` retired.
4. Tests: a fill test with a mixed field (sentinel, mover, drained) counting
   the disk set and checking the window particle's row is the window's; the
   convergence tester's drain cases removed; the fixtures of section 5.
5. Documentation.

Each step builds and passes the gate on its own; steps 1-2 are the science,
3-5 the plumbing.

## 5. Validation

- The fixtures with sampling: solve2D on the 400-particle random-start fixture
  with `nsample = 100` (a quarter of the particles per iteration, so three
  quarters are injected after the first iteration) and solve3D with
  `nsample = 100`; the recovered shifts of the particles first sampled after
  what was the centring stage, against the truth, under the present tree and
  the dynamic one. This is the test the oversight failed.
- The beta-gal bench with `nsample` at solve3D's default (10000 is full
  sampling for 4463 particles; use 1000 to exercise the drain) and the
  class/map quality against the present tree.
- The stream's pool test (`simple_stream_pool2D_tester`) with the sentinel on
  injection; a pool run with the key on, the new particles' shifts after their
  first iteration.
- The cost table of section 2.5 re-measured per path.

## 6. Rulings (Hans, 2026-11-02)

1. The sentinel resets per workflow. The 2D shifts stay the 3D's start.
2. A particle that needs the disk is prepared under the box's mask, whatever
   the stage (B).
3. The re-centring rule is the move that reached the box's bound (the rule in
   force), now the truncation correction of the width (section 2.6).
4. The neighbourhood table's disk candidates come from the subspace pool.
5. `trs_window` is retired with the stage.
6. Both metrics are recorded, the last move and the posterior width, combined
   as section 2.6: the width chooses the next box at three standard
   deviations, the move at the bound inflates it to the radius searched, the
   prior in the likelihood stays the pick's.

Left open for the build: the ladder's rungs (two to start), the bound's
tolerance (half a map step), the width's calibration between the coarse and the
full band (an upper bound until measured).
