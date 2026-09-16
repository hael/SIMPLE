# abinitio2D seeded restart: initialize from the previous 2D clustering

Implementation note, 2026-09-16, revised the same day after review (§8
records the decisions). Implemented the same day on top of master
`6a41a40f0` (uncommitted; §9 lists the files); compilation and the test
runs of §6 are Hans's. Companion pages:
`doc/policies/2D/abinitio2D_policy.md`,
`doc/policies/importance_sampling_fractional_update_policy.md`,
`claude/refine3D_auto_registration_pass.md` (the 3D precedent for the
full-particle `refine=prob` pass).

## 0. Observation and hypothesis

Repeated rounds of `abinitio2D` → class selection → `abinitio2D` on the
cleaned set keep producing junk classes, even when the input set is clean.

Hypothesis (Hans): the run is initialized randomly, so every round
re-creates junk classes from scratch regardless of how clean the data are.
The fix has two coupled parts: (i) seed the references and the particle
partition from the previous 2D clustering, balanced in population and
representative of the previous class distribution even when the requested
`ncls` differs; (ii) since the seed is already a good solution, skip the
stochastic exploration stages (Gaussian-reference stage, sticky random
subset, SNHC) and go straight to the probabilistic stages, opened by one
full-particle pass, the way `refine3D_auto` now starts with a registration
pass.

Constraints set at review: the seed is built from project metadata only
(class labels, `corr`, populations); no image is read, registered or
split before the search starts. The low-pass schedule and the stage
schedule from the first probabilistic stage onward stay exactly as they
are today.

## 1. What the code does today

Everything below is read from `simple_commanders_abinitio2D.f90`,
`simple_abinitio2D_controller.f90`, `simple_cluster2D_strategy.f90`,
`simple_strategy2D_matcher.f90`, `simple_strategy2D_alloc.f90`,
`simple_strategy2D_srch.f90`, `simple_matcher_smpl_and_lplims.f90` and
`simple_eul_prob_tab2D.f90`.

### 1.1 Initialization is random, and the previous clustering is erased

- `exec_abinitio2D_workflow` defaults `cls_init=rand`. On a fresh run
  (`start_stage == 1`) it calls `spproj_field%delete_2Dclustering` and
  `clean_entry('updatecnt','sampled')` before stage 1: class, in-plane and
  shift of every particle are wiped. Nothing of the previous round survives
  except `state` (the selection).
- References: `init_cluster2D_refs` (cluster2D strategy) →
  `init_standard_refs`: `rand` = `noise_imgfile` (pure noise images),
  `ptcl` = `random_selection_from_imgfile` (random raw particles),
  `randcls` = random class labels + `make_cavgs`. `refs=` on the command
  line bypasses generation (`inirefs` rescales the stack to `box_crop`) but
  does not change the schedule below, and the project partition is still
  erased.
- Iteration 1 of stage 1 runs with `extr_iter=1`, which forces
  `ctrl%l_greedy=.true.` (`init_ctrl`): a greedy assignment of every
  sampled particle against noise references, i.e. an essentially random
  partition. `is_fresh_2D_start` is true, so `prep_strategy2D_glob` zeroes
  all shifts and `prep4srch` draws a random "previous" class per particle.

### 1.2 The stage schedule

`NSTAGES_CLS=6`, `PHASES=[4,6]`, `EXTR_LIM_LOCAL=20`, `ITS_INCR=5`,
`PROBREFINE_STAGE=3`.

| stage | lp | refine | extr | refs / regularization | sampling (only if `nptcls_eff > nsample`) |
|---|---|---|---|---|---|
| 1 | lpstart | `snhc_smpl` | extremal, greedy at it 1 | noise/`refs=`, `gauref=yes`, `ml_reg=no`, `trs=0`, `center=no` | sticky random subset (`sample4update_rnd` once, `sample4update_reprod` after), no carry-over |
| 2 | march | `snhc_smpl` | extremal | `cavgs_iterNNN`, `ml_reg` | `sample4update_cnt`, fractional restore |
| 3-4 | march | `prob_snhc` (sparse) or `prob` | extremal (`extr_iter` tracks `which_iter`) | same | same |
| 5 | lpstop | `prob_snhc`/`prob` | `extr_iter = extr_lim+1` (off) | same | same, `fillin=yes` when sampled |
| 6 | e/o (FRC band) | `prob` | off | same | same, `fillin=yes` when sampled |
| terminal | – | `greedy`, all particles, 1 it | – | – | only if any stage sampled |

Phase-1 stage lengths are `istage*maxits/PHASES(1)` (`phase1_stage_endit`,
guarded by `sample_coverage_nits`): with the defaults, stage 1 ends at
iteration 5, stage 2 at 10, stage 3 at 15, stage 4 at 20; phase-2 stages
run `nits_per_stage`. Note that `NSAMPLE_DEFAULT_2D=200000`: for sets
below 200k active particles `update_frac=1`, no sampling happens anywhere
and there is no terminal pass. "Sticky sampling" therefore only exists on
large sets; what is always present in stages 1-2 is the noise start, the
Gaussian reference regularization, `trs=0`, and the extremal SNHC
neighbourhood.

The commander already supports starting at a later stage: the stream
checkpoint path (`execute_abinitio2D_staged` with `start_stage > 1`) sets
`endit` on the cluster2D command line, skips `inirefs` and
`delete_2Dclustering`, and validates/rebuilds the canonical sigma2 state
through `ensure_resume_sigma_state`. The seeded restart reuses exactly this
entry.

### 1.3 Two mechanisms that can manufacture or preserve junk classes

Neither is proven; both are testable with the runs in §6.

1. Random start. After the greedy iteration 1, every class average is the
   mean of a random subset: a blurred "everything" image. Differentiation
   relies on the SNHC neighbourhood shrinking over 20 iterations. Classes
   that never acquire a coherent seed (or lock onto a mis-centred subset,
   `trs=0`) end as sinks for particles nothing else wants. This is the
   mechanism the hypothesis names.
2. Dense probabilistic assignment is class-driven, not particle-driven.
   `ref_assign_likelihood` sorts every class column, then each round every
   class offers its best remaining particle and one class is drawn from the
   top-K softmax of that frontier. `class_exists` is forced `.true.` for all
   classes ("classes must be able to recover"), so a junk class keeps being
   fed; it cannot starve. This keeps populations balanced by construction
   but also means that once a junk reference exists, nothing in the
   assignment removes it. If seeding alone does not cure the symptom, this
   is the next suspect (see §6, T1 outcome b).

### 1.4 Existing pieces the design reuses

- `oris%get_class_sample_stats` → `class_sample` (per class: `pinds`
  sorted best-to-worst by `corr`, `pop`); `sample_balanced_1/2` (the
  "increment every class that still has particles" allocation);
  `split_class`, `expand_classes`, `remap_cls` (`simple_oris_reshape.f90`,
  `simple_oris_sampling.f90`). All metadata-only.
- `make_cavgs` from labels (the `cls_init=randcls` path), writing
  merged/even/odd stacks and `FRCS_FILE`; `inirefs` (rescale a `refs=`
  stack to `box_crop`).
- The checkpoint-resume entry described in §1.2.
- The terminal greedy pass (`execute_terminal_pass`): the pattern for a
  one-iteration, all-particle cluster2D child command line.
- The `refine3D_auto` registration pass (`run_registration_pass`): the
  before/after reassignment diagnostic.

`cls_split` (diffusion-map/kPCA + k-medoids) and `transform_ptcls` exist
and could split a parent on image content; they are deliberately not used
here (review constraint: no image processing before the search).

## 2. Contract

New value `cls_init=prev` for `abinitio2D` (advanced UI). `cluster2D`
rejects it in v1 (its own non-virgin branch of `init_cluster2D_refs`
already regenerates references from existing labels; not touched).
`abinitio2D_chunks` and the stream keep `cls_init=rand` explicitly and are
out of scope.

The only hard error: no previous clustering, i.e. `ptcl2D` is virgin or no
`state>0` particle carries `class >= 1`. Everything else is repaired with a
warning that states the count:

- `cls2D` absent or without `state`: every class index present in `ptcl2D`
  counts as accepted.
- Accepted parents = classes with `cls2D%state>0` (when available) whose
  active population is `>= MINCLSPOPLIM` (5). `state>0` particles whose
  class is rejected, absent from `cls2D`, or below the floor become
  unassigned (§3.3).
- No even/odd partition at all (`get_nevenodd() == 0`): `partition_eo`,
  before the sigma2 state is validated or rebuilt, since that state is
  built on the halves. Otherwise `eo` is kept untouched. There is no
  per-particle repair: for particle fields `isthere('eo')` is false for
  `eo=0` (`oriparam_isthere` treats zero-valued parameters as absent), so a
  per-particle test cannot distinguish "even" from "missing"; the first
  implementation had one and relabelled every even particle after
  `calc_pspec`, which the pass then rejected as an empty sigma2 half.
  Missing `e3`/shift: zero.
- `ncls` (required UI input) is the target `K`. `M` = number of accepted
  parents. `K == M`, `K > M` and `K < M` are all supported.
- No previous class-average stack is needed: seed references are rebuilt
  from the partition.

Outputs in addition to the normal ones: `seed_lineage.txt`
(`new_class parent_class seed_pop`, extended with the final population at
the end of the run) and the log lines of §4.3.

## 3. Seed partition: balanced and representative, metadata only

Let `N` be the number of active particles, `N_i` the population of accepted
parent `i` (`i = 1..M`), and `n* = N/K` the target population per seed
class. The whole procedure reads `class`, `state` and `corr` from `ptcl2D`
and `state` from `cls2D`; it never opens a stack.

### 3.1 Allocation of seed classes to parents

Each parent receives `c_i` seed classes with `Σ c_i = K`, by largest
remainder (Hamilton) on `K·N_i/N`:

1. `c_i = floor(K·N_i/N)`; if `K >= M`, raise every `c_i` to at least 1.
2. Distribute the remaining `K − Σ c_i` one at a time to the parents with
   the largest fractional remainder (ties by larger `N_i`).
3. If `K < M`, the `M − K` parents with the smallest `N_i` (after step 2
   those are the ones with `c_i = 0`) are dropped from the seed (§3.3).

Properties: expected seed population is `n*` for every seed class
(balanced); the number of seed classes per parent is proportional to its
population, so the seeded reference set represents the previous view
distribution rather than the previous class count (representative). A
parent with `N_i < n*/2` may still get one class when `K >= M`; that is
intentional (no accepted parent is silently lost when the user asked for
at least as many classes).

### 3.2 Splitting a parent into `c_i > 1` children: rank interleaving

Sort the parent's particles best-to-worst by `corr` (this is what
`get_class_sample_stats` already returns) and deal them out round-robin:
rank 1 → child 1, rank 2 → child 2, …, rank `c_i` → child `c_i`, rank
`c_i+1` → child 1, and so on. Children keep the parent's `e3`/shift (same
frame, nothing to transform).

What this gives: children of equal size (to within one particle) and of
identical `corr` distribution, each a faithful lower-count copy of the
parent, so every seed reference is a clean average of a clean class. The
children start indistinguishable and the pass distributes the parent's
population between them by likelihood; the following stages differentiate
them on whatever structure the parent contained (the situation an
over-populated view is in at the end of any normal run, where it occupies
several classes). That is redundancy, not junk, and the population
balancing of §1.3 item 2 works for it rather than against it.

Rejected alternatives, both metadata-only:

- Contiguous `corr` quantile blocks (best block → child 1, …, worst block
  → child `c_i`): children are balanced in size but deliberately unequal
  in fit; the worst block of a clean parent is a low-SNR average of the
  same view, i.e. exactly the kind of reference the seed is meant to
  avoid creating, and under the frontier assignment it cannot starve.
- `split_class` (random halves): statistically the same as interleaving
  but non-deterministic and not exactly balanced in `corr`. Interleaving
  is the deterministic version of it.

### 3.3 `K < M` and unassigned particles

Particles of dropped parents, of rejected classes and of parents below the
population floor are set to `class=0` and keep their `e3`/shift. They do
not contribute to the seed references (`cavger_transf_oridat` skips
`icls<1 .or. icls>ncls`) and are assigned during the pass (§4.1) like any
particle whose previous class is invalid: `fill_tab_range` skips the shift
seed (`cxy=0`) and evaluates every class; `prep4srch` draws a random
populated "previous" class as the fallback. This path exists today; §7
lists what to confirm about it.

Decision at review: drop, do not merge. Merging the dropped parents into a
similar surviving parent would need class-average registration and pose
composition, and dropping the least populous parents is the honest reading
of `K < M`.

### 3.4 Writing the seed

1. Write `class` (1..K) for every active particle, `class=0` for the
   unassigned set; `clean_entry('updatecnt','sampled')`; write `ptcl2D`.
2. `make_cavgs` on the seed labels (`ml_reg=yes`) → `start2Drefs.mrc`
   (+ `_even`/`_odd`) at the working `box_crop`/`smpd_crop`, the same
   names `inirefs` and `init_cluster2D_refs` use today, then `refs=` on the
   pass command line. `FRCS_FILE` from this `make_cavgs` is what the pass
   would read if it were not `l_lpset`; it is (§4.1), so it is only
   informative.
3. `seed_lineage.txt` (seed class, parent, seed population; the final
   population is filled in after `gen_final_cavgs`). `cls2D` is left as the
   previous run wrote it until the final `make_cavgs` replaces it: writing
   K entries into it mid-run would break `map_cavgs_selection` consistency.

## 4. Search schedule under `cls_init=prev`

Decision at review: no new low-pass or stage logic. The run is the
existing run entered at `PROBREFINE_STAGE` (3), preceded by one
full-particle pass. Everything from stage 3 to the terminal greedy pass is
untouched: same `lp` per stage, same `refine` policy, same iteration
counts, same extremal schedule, same sampling and fractional restore.

### 4.1 Seed pass

`it_pass = phase1_stage_endit(maxits, FRAC_UPDATE_STAGE, update_frac)`,
the iteration at which stage 2 would have ended (10 with the defaults), so
that stage 3 starts from `endit = it_pass` and runs exactly the iterations
it runs today (11..15). One `cluster2D` invocation, built the way
`execute_terminal_pass` builds its command line:

- `refine=prob` (dense: every class is evaluated for every particle;
  `prob_snhc` would only evaluate a sparse neighbourhood around the seed
  label, which defeats a global reassignment);
- every active particle: `update_frac`, `nsample` and `fillin` deleted;
- `maxits=minits=1`, `startit=which_iter=it_pass`, `refs=start2Drefs.mrc`,
  output `cavgs_iter<it_pass>` which stage 3 picks up through the existing
  `CAVGS_ITER_FBODY//int2str_pad(iter-1,3)` naming;
- `lp = stage_parms(PROBREFINE_STAGE)%lp` (`l_lpset`), i.e. the band stage
  3 opens at today;
- `extr_iter=extr_lim+1` (`neigh_frac=0`: no extremal class neighbourhood,
  in-plane sampled from the top-2 peaks in `fill_tab_range`, the near-greedy
  regime the terminal and phase-2 passes already run in);
- `ml_reg`, `center`, `trs`, `objfun`, `box_crop`, `smpd_crop`,
  `restore_cavgs=yes` as for any stage-3 iteration.

`startit=it_pass > 1` matters: `is_fresh_2D_start` must be false, otherwise
`prep_strategy2D_glob` zeroes the shifts after `prob_tab2D` built its table
against the stored ones, and `prep4srch` discards the seed class. Because
the cluster2D strategies only `clean_entry('updatecnt','sampled')` when
`startit==1`, abinitio2D does it itself (§3.4 step 1). The pass sets
`updatecnt=1` on every active particle (`sample4update_all` with
`incr_sampled=.true.`), so the count-biased sampling of stages 3+ on large
sets starts uniform; stage 3 keeps `sample4update_cnt` with fractional
restore, and the terminal greedy pass runs whenever any stage sampled, as
today.

Sigma2: the resume entry runs `ensure_resume_sigma_state` before stage 3
(reuse the project's committed canonical state if its identity validates,
else `calc_pspec`); the seed pass runs after it, so it sees the same state.

### 4.2 Stages after the pass

Stages `PROBREFINE_STAGE..nstages` and the terminal pass, unchanged. In
code this is `execute_abinitio2D_staged(cline, PROBREFINE_STAGE, 0,
it_pass, ...)` semantics with the seed and the pass inserted before the
stage loop; the checkpoint-only validation (`l_checkpoint`) is not
engaged.

### 4.3 Diagnostics

- `>>> ABINITIO2D SEED: M accepted parents -> K seed classes; dropped d
  parents (p particles), unassigned u particles; seed pops min/median/max`.
- `>>> ABINITIO2D SEED PASS: refine=prob, all particles, iteration
  it_pass, lp=x A (stage 3 limit)`.
- `>>> ABINITIO2D SEED PASS REASSIGNED: c% changed class, s% moved shift >
  1 px, r% of seed classes retained >= 50% of their members` — poses read
  from the project before/after the pass, same construction as
  `REGISTRATION PASS REASSIGNED`.
- `seed_lineage.txt` extended after the final `make_cavgs` with the final
  population per seed class, so a class's history (parent → seed pop →
  final pop) can be read in one place when judging which classes went bad.

## 5. Code touch points

- `simple_parameters.f90` (`cls_init` doc string), `simple_parameters_phases.f90`
  (accept `prev`), `simple_ui_params_common.f90` (`cls_init` choices),
  `simple_ui_cluster2D.f90` (help text for abinitio2D only).
- `oris%reseed_classes` (`simple_oris_reshape.f90`, interface in
  `simple_oris.f90`, next to `split_class`/`expand_classes`): Hamilton
  allocation on `class_sample` populations, rank interleaving, label
  writing, lineage arrays. No image types; the repo convention of adding
  behaviour to existing types was preferred over a detached module. Unit
  test `test_reseed_classes` in `simple_oris_tester.f90` (K = M identity,
  K > M proportional and balanced with interleaved children sharing the
  `corr` distribution, K < M dropping the smallest parent, rejected and
  state=0 particles ignored).
- `simple_abinitio2D_controller.f90` (stage policy owner): nothing in the
  stage tables; `abinitio2D_seed_pass_iter` (= `phase1_stage_endit` of
  stage `PROBREFINE_STAGE-1`, never below 2) and
  `set_cline_cluster2D_seed_pass`, mirroring the terminal-pass cline
  construction.
- `simple_commanders_abinitio2D.f90`: when `cls_init=prev`,
  `start_stage=PROBREFINE_STAGE` with `endit=it_pass` through the existing
  resume path (so `ensure_resume_sigma_state` runs, `inirefs`,
  `delete_2Dclustering` and `partition_eo` do not), then
  `seed_from_previous_clustering` (validation/repair per §2, `make_cavgs`
  → `start2Drefs*`), `execute_seed_pass` with the reassignment diagnostic,
  and `write_seed_lineage` after `gen_final_cavgs`. Rejects the mode with
  stream checkpointing or `nstages < 3`.
- `simple_cluster2D_strategy.f90`: `init_standard_refs` names the mode in
  its rejection message.
- `doc/policies/2D/abinitio2D_policy.md`: §2 (seeded entry), §3 (ownership
  of the seed module), §4 (no sticky stage when seeded), §7 checklist item
  "does `cls_init=prev` preserve existing `eo` and never call
  `delete_2Dclustering`?".

Nothing in `simple_strategy2D_matcher.f90`, `simple_eul_prob_tab2D.f90`,
`simple_cluster2D_strategy.f90`, the classaverager or the streaming code
changes. The default `cls_init=rand` schedule must be byte-identical in
its `>>> STAGE` summary lines before and after.

## 6. Test protocol (runs are Hans's)

T1, the hypothesis. One project where the loop was observed: round r
output → selection → prune → round r+1. Run round r+1 three ways on the
identical pruned project with the same `ncls`: A `cls_init=rand` (today),
B `cls_init=randcls` (random labels → blurred averages, not noise: tests
whether any non-noise start suffices), C `cls_init=prev`. Count junk classes
with the same criterion used for the selection. Outcomes: (a) C ≪ A:
hypothesis confirmed, ship; B ≈ C says the noise/greedy start is the
culprit, B ≈ A says the seed content is; (b) C ≈ A: junk is not an
initialization artefact; next suspect is the class-driven frontier
assignment (§1.3, item 2), to be tested by logging per-class population
and mean assignment distance per iteration for the classes that end up
junk.

T2, the pass. From `SEED PASS REASSIGNED`: near-zero reassignment at the
stage-3 band means the pass could be dropped and stage 3 entered directly
(same reading as the 3D registration-pass test); a few percent with the
final class quality responding means it is doing its job.

T3, `K ≠ M`. `K = 2M` and `K = M/2` on the same project: seed populations
(should be ≈ `n*`), whether interleaved children of one parent have
differentiated by the e/o stage or still show the same average (read from
`seed_lineage.txt` plus the ranked final stack), and the junk count.

T4, large set (> 200k active particles): sampled stages after the pass,
fractional restore, terminal pass; compare wall time against the unseeded
run (expected cheaper: the pass costs one full iteration, stages 1-2 cost
ten).

Regression: an unseeded default run before/after the change, `>>> STAGE`
lines identical.

## 7. Things to verify during implementation (not design decisions)

- `class=0` particles on the search side: shift semantics in
  `fill_tab_range` (seed `cxy=0` with the stored shift kept) match what an
  out-of-range `icls_prev` does today; `record_sparse_eval` shifts are
  increments in the same frame; `prep4srch` fallback for `prev_class=0` on
  a non-fresh start.
- The resume path with `start_stage=PROBREFINE_STAGE` and `endit=it_pass`
  without `l_checkpoint`: `set_cline_cluster2D_stage` for stage 3 yields
  `startit=it_pass+1`, `maxits=15` (defaults), `refs=cavgs_iter010`.
- `sigma2_group_iter`/canonical sigma2 naming with the pass at
  `which_iter=it_pass` after a reused or rebuilt state.
- `cleanup_distributed_iteration_artifacts` at `startit=it_pass` without
  `update_frac` still deletes stale partial sums (it does when
  `.not. l_update_frac`).
- `make_cavgs` at `box_crop`/`smpd_crop` from the abinitio2D command line
  (which carries them) writes `start2Drefs*` at the working scale, so
  `inirefs`-style rescaling is not needed; if it does not, route through
  `inirefs`.
- `partition_eo` on a field where some particles already have `eo` only
  fills the missing ones (else fill them explicitly).

## 8. Decisions (Hans, 2026-09-16)

1. Split engine: no image processing before the search. Rank interleaving
   by `corr` (§3.2). PC1 bisection and `cls_split` are out.
2. `K < M`: drop the smallest parents, let the pass reassign their
   particles.
3. Low-pass: unchanged; the pass uses the stage-3 limit, no FRC-derived
   band.
4. Stages after the pass: unchanged from stage 3 on; no stage selection
   logic.
5. Switch: `cls_init=prev`.
6. Refuse only when no previous clustering exists; everything else is
   repaired with a warning.

The seed pass is reported by its own `>>> ABINITIO2D SEED PASS` line; the
`>>> STAGE` summary block is unchanged.

## 9. Files touched (uncommitted)

`src/main/ori/simple_oris.f90`, `src/main/ori/simple_oris_reshape.f90`,
`src/main/ori/simple_oris_tester.f90`, `src/main/simple_abinitio2D_controller.f90`,
`src/main/commanders/simple/simple_commanders_abinitio2D.f90`,
`src/main/strategies/parallelization/simple_cluster2D_strategy.f90`,
`src/main/params/simple_parameters.f90`, `src/main/params/simple_parameters_phases.f90`,
`src/main/ui/simple_ui_params_common.f90`, `doc/policies/2D/abinitio2D_policy.md`,
this note. Not compiled; the §7 items are still to be confirmed at the
first run.
