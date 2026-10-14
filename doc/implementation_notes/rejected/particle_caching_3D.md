# Caching for 3D refinement: what the code allows today, and what would pay

Written 2026-10-13 for the discussion with Susan and Joe. These are thoughts
for a decision, read against the code as it stands; no code was written.
The companion is `particle_io_few_particle_stacks.md`, which covers how the
bytes are fetched; this note covers what is done with them and whether any of
it can be kept. The last two sections were added the same day from the
ribosome benchmark files and change the conclusion.

## What the ribosome benchmark files say

The files are those of the consolidated base run of the ribosome
(`solve3D`, 29,101 particles, 1,666 sampled per state and iteration for 3
states, 5 parts of 16 threads, box 336, one spinning disk, 7,103 s in all;
`/Users/elmlundho/Processing/ribo/bench`, tally in its `tally/`
subdirectory). The 78 refine3D iterations took 6,572 s by the strategy
files; summed per stage, with the matcher figures taken from the slowest
part of each iteration because the iteration waits for it:

| stage (box_crop, mode) | iterations | total | probabilistic pre-step | reprojection model | other setup | scheduling overhead | matcher, slowest part | of which particle preparation | of which partial reconstruction | assembly | of which volassemble |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 88, shc | 12 | 403 | 8 | 60 | 31 | 78 | 79 | 10 | 11 | 92 | 49 |
| 112, shc | 20 | 925 | 13 | 157 | 56 | 128 | 266 | 26 | 22 | 213 | 133 |
| 132, prob | 17 | 1,357 | 627 | 121 | 51 | 117 | 103 | 12 | 24 | 258 | 188 |
| 160, prob | 15 | 1,429 | 743 | 124 | 41 | 101 | 101 | 12 | 27 | 249 | 195 |
| 182, prob | 12 | 2,259 | 1,642 | 199 | 34 | 80 | 101 | 11 | 27 | 149 | 101 |
| 182, greedy; 336, sigma | 2 | 199 | 0 | 74 | 6 | 13 | 95 | 11 | 0 | 3 | 41 |
| all | 78 | 6,572 | 3,033 | 734 | 217 | 517 | 745 | 81 | 111 | 964 | 707 |
| share | | 100% | 46% | 11% | 3% | 8% | 11% | 1.2% | 1.7% | 15% | 11% |

Seconds; 362 s (5.5%) are unaccounted in the strategy files. "Scheduling
overhead" is the matcher/scheduler phase minus the slowest part's own time,
so it holds script generation, process launch, the parts' startup before the
matcher timer starts, polling (1 s interval) and the merge of the alignment
documents.

Three things follow.

First, the particle pipeline is small here. Particle preparation (read,
normalise, transform, crop, shift, mask, polarize) is 1.2% of the run and
partial reconstruction (read again, prepare, grid, write) 1.7%, taken from
the slowest part; the probability-table parts prepare the same particles once
more and have no timer, so add about another percent. Everything a particle
cache can touch is three to four percent of this run. In the original layout
(19,520 stacks) the same run took 8,011 s in the stages, so the scattered
reads cost about 1,440 s, 18% of that run (with the caveat that the two runs
reached different three-state solutions, so not all of the difference is
reading); that is the per-stack open and seek cost, which the read-path work
and a cache attack from different sides, and it is the ceiling for both.

Second, the probabilistic pre-step is 46% of the run and 73% of a stage-5
iteration (138 of 190 s). In `refine=prob` the exhaustive scoring of every
sampled particle against every reference moved out of the matcher into the
`prob_tab` parts (`fill_tab_range` in `simple_eul_prob_tab.f90`: all 3,000
references at the seed shift, then shift refinement of the best ones), which
is why the matcher's own search falls to under a second per part in the prob
stages. The pre-step has no benchmark file of its own, so the split between
the parts' table fill, the shift refinement, the table write and read, the
master's assignment (`ref_assign`) and two process launches is unknown. It
needs the matcher's timer pattern in `exec_prob_tab` and `exec_prob_align`
before anything is decided about it.

Third, the per-iteration work that does not scale with particles is 37%:
the master's reprojection model (volumes read, masked, filtered, projected
into 1,000 × 3 × 2 polar references and written for the parts to read; 4 to
17 s per iteration), volume assembly (4 to 13 s), the scheduling overhead
(6.7 s per iteration, twice in the prob stages since `prob_tab` launches
parts as well) and the remaining setup. At 5,000 particles per iteration on
80 threads the matcher's own work is 6 to 14 s, so the iteration is too small
to amortise its orchestration. This is the regime of every stream chunk.

So for this run, and for the stream chunks it emulates, a particle cache
cannot change the picture; the levers are the probabilistic pre-step and the
fixed per-iteration cost. The sections below stand as the analysis of what a
3D cache would be and how to build one without repeating the August result;
they are worth acting on only where the particle pipeline dominates, which is
the whole-dataset run on scattered stacks from a cold disk, not the chunk.

## The short version of the cache analysis

The wish is to shift, phase-flip, mask, crop and polarize each particle once,
and afterwards only update the shifts on the cached polar transforms. The
mechanics for that exist. The polar calculator stores every particle's polar
transform for all shells up to the Nyquist of the cropped box
(`pfts_ptcls(pftsz, kfromto(1):interpklim, :)` in `simple_polarft_core.f90`,
with `interpklim = box_crop/2`), so a stage's low-pass limit is a read limit
and nothing else. It shifts a stored particle transform by a phase ramp in
one call (`shift_ptcl`, `simple_polarft_geom.f90`); the shift is in original
pixels and the ramp is `2*pi*k*s/box` at shell `k`, which is the same number
on every crop. And the reconstruction already applies the shift in Fourier
space after its edge taper (`gen_rec_plane`, `simple_matcher_3Drec.f90`), so
on the reconstruction side there is no shift-before-mask problem at all.

What made the August 3D cache slower was arithmetic, not the mask.
Fractional sampling touches each particle once in every N/n iterations
(N active particles, n sampled per iteration). A cache built up front over
all N particles has to be read back many times before it has paid for its
build, and the per-stage crop of solve3D invalidated it before it could. The
fix in `666e95960`, one crop for all stages, made the early stages pay for the
final crop instead. Both have a cheap remedy: fill the cache on first touch
instead of building it, and store the record at the run's largest crop and
crop it further on read, which is exact.

The shift-before-mask problem of the search has an answer that costs one
phase ramp per particle: cache the masked polar transform together with the
shift it was built with, and apply the difference between the current shift
and the cached one. The mask is then off the particle's centre by that
difference. The search accepts exactly this misregistration already: during
the shift search the trial shift is applied to the reference (`pft_ref *
shmat` in `simple_polarft_corr.f90`) while the particle's mask stays put, so
the two masks are misregistered by up to the shift bound on every
evaluation. A per-particle tolerance and a re-preparation beyond it keep the
error at the level the search tolerates anyway.

The number that decides whether any of this is worth doing is not in the
code but in the benchmark files. Every refine3D part writes
`REFINE3D_BENCH_ITER###_PART##.txt` (benchmarking is on unconditionally,
`L_BENCH_GLOB = .true.` in `simple_defs.f90`), with the iteration split into
particle preparation, reference preparation, orientation search, project
metadata I/O and partial reconstruction. A cache can shrink the first of
these and the read-and-prepare part of the last; it cannot touch the search
or the reference projection. On the ribosome chunk run those files put the
particle pipeline at three to four percent of the run (the section above),
so there a cache is a side show whatever its design.

## What a sampled particle costs today

Per iteration, in probabilistic distributed mode, a sampled particle is read
three times in two freshly launched processes: for the probability table, for
the search, and for the reconstruction (the I/O note, section 2.2). Each read
is of the full-box image, 451 KB for the ribosome at box 336.

For the search (`prepimg4align`, `simple_matcher_2Dprep.f90`) the image is
noise-normalised against the full-box mask, transformed at the full box,
Fourier-cropped to `box_crop`, shifted by minus its stored shift (a phase
ramp) and phase-flipped when `ctfflag` is `yes`; then inverse-transformed at
the cropped box, soft-masked at `msk_crop`, padded and transformed again at
the padded cropped box, and interpolated onto the polar grid
(`polarize_ptcl_pft`). The `|CTF|` matrices on the polar grid are evaluated per
batch (`create_polar_absctfmats`) and the rotation transforms of particle
times CTF are memoised (`memoize_ptcls`). The probability table runs the same
preparation (`build_batch_particles3D` from `simple_commanders_prob.f90`).

For the reconstruction (`calc_3Drec` then `prep_imgs4rec`) the image is
noise-normalised, transformed at the full box, cropped and inverse-transformed
(`prep_rec_observation`, cropped branch, `simple_matcher_ptcl_io.f90`), then
edge-tapered, padded and transformed at the padded cropped box; the shift,
CTF and whitening are applied in Fourier space (`gen_rec_plane`) before
gridding.

So the full-box transform, the most expensive single step per particle, is
paid three times per sampled particle per iteration, and nothing derived from
a particle survives from one pass to the next, nor from one iteration to the
next. The matcher kills the polar calculator between the search and the
reconstruction to bound peak memory (`simple_strategy3D_matcher.f90`), which is
why the reconstruction reads again.

## One record serves both passes, exactly

The 2D cache record is the noise-normalised, Fourier-cropped particle, stored
in real space at `box_crop` (`doc/policies/2D/particle_cache_policy.md`,
section 2). In today's 3D code this record is bit-exact for both consumers.
For the search, `prepimg4align_cached` exists and performs the remaining
steps in the same order as `prepimg4align`. For the reconstruction, the
cropped branch of `prep_rec_observation` is `norm_noise_fft_clip_shift(lmsk,
obs, [0,0])` followed by an inverse transform, which is the record itself;
the taper, padding and transform come after it, from the cropped image, in
the uncached path too. The August policy text warned that cached
reconstruction changed the numerics (the taper acting at `box_crop`); the
reconstruction has since been reorganised to crop first in every case, so
that caveat is gone.

Nested cropping is exact as well: a Fourier crop of a cropped record equals
the crop from the original, because the normalisation precedes the crop. One
record at the largest `box_crop` of the run therefore serves every stage of
solve3D, read back at the stage's crop, and the early stages keep their small
search grids.

## The arithmetic under fractional sampling

A particle is touched about once every N/n iterations. A cache built up front
costs N preparations and N writes before the first iteration, and earns back
three full-box preparations and the read bytes per touch. On the ribosome,
N = 436k and n = 10k (solve3D's default `nsample`) give N/n of about 44
iterations before every particle has been touched once; a solve3D stage is
shorter than that, which is why a per-stage cache could never pay. On a
stream chunk, N = 29k and n = 5k, the ratio is 6.

Filling on first touch changes the picture. The first touch costs today's
preparation plus one write of the cropped record; every later touch of that
particle is served from the record. There is no build pass, so the cache is
never slower than today by more than the write bandwidth, and it reaches full
coverage after N/n iterations without anybody waiting for it. Over a run of
T iterations the share of touches served from the cache is about 1 - N/(nT).
For a long refine3D on the whole set it approaches one; for a short run on a
large set it stays small, and that is the honest limit of any cache under
fractional sampling. Susan's observation that copying particles to local
flash pays only in the late stages, when the set is small and the iterations
many, is the same arithmetic seen from the other side.

Sizes for the ribosome at `box_crop` 160: the cropped record is 102 KB per
particle, 44 GB for all 436k; the raw image is 451 KB, 197 GB in all. A polar
record is about 1.4 times the cropped record (`pftsz` times the shells to the
crop Nyquist, 8 bytes each; `pftsz` is about `pi*msk_crop`). The per-iteration
working set at n = 10k falls from 4.5 GB read three times to about 1 GB per
pass. On the production server (251 GB of memory, one spinning disk) the
whole cropped cache stays in the page cache once written; the records of a
part can also live on a node-local disk.

## The shift and the mask

The search prepares the particle centred, shifted by minus its stored shift,
and masks it at the box centre. If the polar transform is cached at the shift
it was built with, s_c, and the current shift is s, applying the ramp for
s - s_c to the cached transform gives the transform of the centred particle
under a mask displaced by s - s_c. The signal is untouched as long as the
displacement stays inside the margin between `msk_crop` and the particle's
radius; what changes is which peripheral noise pixels fall inside the mask.
This is the misregistration the shift search produces on every evaluation,
because the search shifts the reference, whose own mask moves with it,
against the particle whose mask stays put.

A design that stays within what the search tolerates: each cached polar
record carries the shift it was built with; the per-iteration preparation
applies the difference (`shift_ptcl`), evaluates the `|CTF|` matrices and
memoises as today (the memoised rotation transforms depend on the ramp, the
shell power sums do not); when the difference exceeds a tolerance, the
particle is re-prepared from its cropped record and the polar record
replaced. The shift bound `trs`, in original pixels, is the natural tolerance,
since it is the misregistration the search already accepts; a tighter value
is a ruling. In the first solve3D stages, where shifts move by many pixels,
this re-prepares often, at today's cost plus a write; from the middle stages
on the shifts drift by fractions of a pixel and the records stand. Being
"loose with the mask" is then not a separate measure: the user's mask margin
is the tolerance budget, and the tolerance is explicit.

The phase flip belongs in the polar record (flip and shift are both Fourier
multiplications and commute). The `|CTF|` matrices are one evaluation per
polar sample and cheap; keep computing them per batch until a measurement
says otherwise.

One caveat that applies to the polar record and not to the cropped one: the
number of rotations, `2*pftsz`, follows `msk_crop` and so the stage's crop. A
polar record built at the final crop has more rotations than an early stage's
calculator, so either the early stages run with the final `pftsz` (what the
August fix did, and what made them slower) or the polar record is keyed by
stage crop and refilled per stage. With first-touch filling from the cropped
record, a refill costs three cropped-box transforms and the polar
interpolation per particle, not a full-box read, so per-stage polar records
are affordable.

## Where the records can live

In the distributed route every part is a new process per iteration, so the
records have to be in files, one record file and index per part, the part
owning its particle range. The persistent-worker overlay does not change
this: the workers execute generated scripts (`persistent_worker_policy.md`),
so a worker's memory does not survive an iteration either. In the
shared-memory route (`refine3D_inmem_strategy`, one process for all
iterations) the records can stay in memory for the run, as a separate object
next to the polar calculator, which is killed between search and
reconstruction.

For the lifetime, the 2D cache's rules carry over: owned by the workflow
master (solve3D, not its refine3D children) and deleted at exit, one record
set for the whole run. A stream can use the same mechanism with per-row
validity (stack path and modification time per row) instead of the global
fingerprint, since a growing particle set only appends; the ruling of
2026-10-12 against rewriting the stacks into one file stands, because
nothing here rewrites or replaces the stacks, and no full pass is made.

## What this does not fix

The per-stack cost of the first touch (two opens, two header reads, a seek
on a spinning disk) is paid at least once per particle regardless; that is
the subject of the I/O note. The probability table stays a separate process
with its own read. The search and the reference projection are untouched.
The add-on passes of the stream are their own track (ruled).

## The cache steps, if a run shows the particle pipeline to matter

The first step is a cache with the lifetime of one pass and no files: keep
the cropped records of the sample in memory from the search into the
reconstruction, so the third read and its full-box transform disappear. This
is the "pass-sized particle buffer" already on the I/O note's candidate list;
it is bit-exact and small (n/nparts records of 102 KB).

The second is the cropped record on disk, filled on first touch, at the run's
largest crop, read by both passes in every process of the run. This is the
August design with the two repairs. It is exact for both consumers and gives
the full gain of the full-box transform and the read bytes on every touch
after the first.

The third is the polar record on top of it, for the search pass only, with
the snapshot-and-difference rule and per-stage keying. It removes the three
cropped-box transforms and the polar interpolation per particle per search
pass. Whether it earns its tolerance rule depends on what the second step
leaves in "particle preparation".

## What the measurement puts first instead

The global-sigma refactoring is finished and validated first, as planned.
Then the probabilistic pre-step, which has its own plan:
`doc/refactoring_notes/planned/prob_search_speedup.md`. Its section 2
counts what a sampled particle costs in the pre-step from the code (the
exhaustive table against 3,000 references, then a 25-point coarse shift
grid and an L-BFGS-B refinement for each of 17 to 55 best references, the
refinement being the larger part by the counts), and its Phase 0 gives
`exec_prob_tab` and `exec_prob_align` benchmark files of their own so that
the 46% is split into named parts before anything is built.

The per-iteration fixed cost (37%) is a separate question with its own
candidates: the reprojection model is projected once by the master and read
by every part twice per iteration in the prob stages; volume assembly runs
at 4 to 13 s per iteration; each distributed step pays about 6.7 s of launch,
startup and merge. For a chunk of this size on one node, the shared-memory
route (`refine3D_inmem_strategy`) pays none of the launches and none of the
model writes and reads, at the price of running the matcher on all threads of
one process; whether it is faster is a measurement, not a guess.

The whole-dataset run on scattered stacks is the case the cache analysis was
written for. Its benchmark files are the ones to tally next; if the particle
pipeline is a large share there, the three cache steps apply, in that order,
after the read-path steps of the I/O note that remove the first-touch cost.
