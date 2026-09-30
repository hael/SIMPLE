# Continuous Cartesian pose refinement (`refine=cont`): refactoring plan

- Status: contract accepted by the maintainer; one ruling outstanding (O3); documentation only, no source change made
- Scope: five-parameter Cartesian pose refinement and its integration in `refine3D` and `refine3D_auto`
- Source baseline: `master` at `712b4d86`
- Last updated: 2026-09-30

This document is self-contained. It holds the contract, the current-state
inventory, the target architecture, the phased refactoring, and the test plan.
Per repository policy the maintainer compiles and runs; nothing here has been
compiled or run.

## 1. Goal

Give continuous Cartesian pose refinement the same boundaries as the rest of
SIMPLE:

- `refine=cont` is a refinement mode of `refine3D` that runs with no polar
  object alive, and `refine3D_auto` accepts it;
- the matcher references exactly one continuous-pose symbol, `strategy3D_cont`;
- the numerics sit in two encapsulated classes, a calculator owned by
  `builder` and an optimizer owned by the strategy, as in the polar path;
- both objectives follow `objfun` with the polar definitions, and the
  reference follows the polar branch;
- a pose polish after discrete search is available as a scheduled Cartesian
  pass;
- there is no separate `refine3D_pose_cont` program.

Two rules hold for the whole sequence. Ownership moves and scientific changes
are separate commits; Phase 3 is the only phase that changes numerics.
Existing tests are reused and moved with the code they test; every new
behavior lands with its test in the same phase (section 10).

## 2. Current state

Production code (3,147 lines in eight modules):

| Module | Lines | Holds |
| --- | --- | --- |
| `volume/simple_cartesian_pose_refiner.f90` | 1,182 | Cartesian reference, particle preparation, objective and gradient, shift and joint Levenberg-Marquardt (LM) |
| `interp/simple_cartesian_fourier.f90` | 138 | Stateless embedding, extraction and KB gathers; also used by `reconstructor_pcg` and flex PCA |
| `strategies/search/simple_pose_cont_refine3D_adapter.f90` | 745 | Policy records, reference and particle workspaces, reference file I/O, observation preparation, coordinate conversion, transaction control |
| `strategies/search/simple_strategy3D_pose_cont.f90` | 276 | Standalone `refine=pose_cont` strategy |
| `strategies/search/simple_pose_cont_run_stats.f90` | 377 | Per-thread accounting and per-iteration report files |
| `commanders/simple/simple_commander_refine3D_pose_cont.f90` | 82 | Program `refine3D_pose_cont`, extending `commander_refine3D_auto` |
| `strategies/parallelization/simple_refine3D_pose_cont_workflow.f90` | 196 | Stage sequencing for that program |
| `ui/simple/simple_ui_refine3D_pose_cont.f90` | 151 | Its UI |

Shared files touched: `simple_strategy3D_matcher`, `simple_matcher_refvol_utils`,
`simple_matcher_ptcl_batch`, `simple_refine3D_strategy`, `simple_strategy3D_prob`,
`simple_euclid_sigma2`, `simple_parameters*`, `simple_ui_refine3D`,
`simple_exec_refine3D`, `simple_commanders_refine3D`, `simple_convergence`,
`simple_defs_ori`, `simple_binoris`, `simple_refine3D_fnames`.

What is wrong with it:

1. The matcher hosts two continuous-pose pipelines (`pose_cont=yes` after a
   discrete winner, `refine=pose_cont` instead of the search). 124 of its 717
   lines mention continuous-pose state; it holds the only three `select type`
   downcasts in non-test strategy code and a 55-line numerical routine
   (section 8.1).
2. The adapter exports eleven public types and nine free procedures and has no
   type named after the module. Most procedures forward to the refiner or wrap
   one or two `image`/`ori` calls (section 8.2).
3. `strategy3D_pose_cont` cannot satisfy the `strategy3D` contract: the base
   embeds the PFTC search object `s`, which it never initialises, and it needs
   a seven-argument `bind_context` and a `get_result` that are not in the
   abstract interface.
4. The objective ignores `objfun`: both `euclid` and `cc` runs optimize a
   normalized correlation on sigma-whitened data. `corr` is never written.
5. Shared helpers import the adapter and take `optional` feature flags
   (`cartesian_only`, `pose_cont_particles`); `euclid_sigma2` gained a second
   constructor and a setter that lets an outside module write `sigma2_part`.
6. Parameters: four switches, string comparisons in four modules, validation
   keyed on a program name, `cc_emit_sigma` silently forced on, LM policy
   hard-coded in adapter type defaults.
7. The particle record grew from 50 to 52 reals for two telemetry flags.
8. Reporting exists three times: per-particle fields, a stats module with its
   own per-partition files, and `system_clock` timers in the particle loop.
9. `commander_refine3D_auto` was given four override hooks so that a second
   commander could inherit from it.
10. 32 `error stop` statements in library code, several reachable from the
    OpenMP particle loop.

What is sound and is kept: rollback to the input pose after a failed solve;
explicit shift and rotation bounds anchored at the seed; explicit even/odd
reference selection; the KB derivative in `simple_kbinterpol`; the stateless
`simple_cartesian_fourier` layer; the numerical tests.

## 3. Contract

Accepted by the maintainer on 2026-09-30.

| Id | Ruling |
| --- | --- |
| C1 | The matcher uses one strategy object, `strategy3D_cont`, and no other continuous-pose symbol. |
| C2 | `simple_pose_cont_refine3D_adapter` is dissolved; each piece moves to its owner. |
| C3 | The feature is restructured as encapsulated classes: private components, `new`/`kill`, type-bound behavior, one type per module named after it. |
| C4 | Both objectives are supported and follow `objfun`, mirroring the polar implementation. |
| C5 | Cartesian `cc` is unweighted by sigma2, as in the polar branch. |
| C6 | Cartesian `euclid` is the polar loss, normalized by the particle's weighted power, with `exp(-L)` as the stored score. |
| C7 | `refine=cont` is a refinement mode of `refine3D`, testable on its own and independent of the polar machinery. |
| C8 | Pose polishing after discrete search is supported. It is scheduled per pass: a `refine=cont` pass follows a discrete pass. A per-particle second stage inside one pass is not part of this plan; it can be added later on the same objects. |
| C9 | The Cartesian and polar objectives are not directly comparable. No decision compares a value from one with a value from the other. |
| C10 | The reference follows the polar branch of `refine3D`, including how deapodization is handled. |
| C11 | `refine3D_auto` supports `refine=cont`. In that mode it bypasses all polar machinery, including the discrete greedy initialization, and runs a purely Cartesian refinement. It is therefore a continuation mode: the project must already hold 3D poses from a discrete workflow. |
| C12 | `refine=cont` in `refine3D_auto` follows the `nsample` policy. |
| C13 | `refine3D_pose_cont` is not a separate program. |
| C14 | The refined Cartesian score is persisted in `corr`. |
| C15 | The particle record returns to 50 reals. The improved flag uses the unused slot 42; "attempted" equals "sampled in a Cartesian pass" and needs no field. |
| C16 | One LM route is kept, chosen from the Phase 0 baseline; `pose_cont_route` is removed and LM policy becomes private constants of the optimizer. |
| C17 | Names: folder `src/main/cftc`, classes `cartesianft_calc` and `cartesianft_pose_opt`, builder field `cftc`, and an abstract `strategy3D_pftc` parent for the polar strategies. |
| C18 | The end-to-end gate gets its own high-level CTest entry, `cont_refine3D_1jxy`; `SIMPLE_CTEST_BUDGET` goes from 29 to 30. |
| C19 | The polish switch is `pose_cont=yes\|no`. When on, a Cartesian pass follows every discrete iteration, in `refine3D` and in `refine3D_auto`. |
| C20 | `refine3D_auto` keeps `objfun=euclid`; `objfun=cc` with `refine=cont` is available through `refine3D` only. |
| C21 | `objfun_den` with `refine=cont` is rejected in parameter validation for now. |

Facts verified in the source that the contract relies on:

- C10. The polar branch prepares the reference in `read_mask_filter_refvols`,
  pads it by 2 in real space (`pad_fft`), expands it (`expand_cmat`) and
  gathers with normalized KB weights scaled by pf^3 (`interp_fcomp_oversamp`).
  It applies no inverse KB envelope at the gather; deapodization happens when
  the volume is reconstructed (`reconstructor%invenv1d`). The Cartesian
  production path (`new_physical_reference`) already does the same on the same
  prepared volume, and `test_matched_projector_boundary` pins its gather
  against `interp_fcomp_oversamp`. The unused `new_inverse_envelope_reference`
  is removed.
- C14. The discrete search recalculates the previous score each iteration
  (`strategy3D_srch%prep4srch` calls `gen_objfun_vals` at the previous
  reference), so a Cartesian value in `corr` does not enter a later polar
  search. A Cartesian pass writes the Cartesian score for every particle it
  processes, accepted or not. Under `nsample` or `update_frac` it processes a
  sample, so the project then holds Cartesian scores for that sample and
  earlier scores for the rest. One consumer ranks the stored `ptcl3D` score
  across particles: `make_projdir_class_samples` in `refine3D_states`
  (`balance=yes` with fractional updates) sorts `corr` within projection
  bins. That is a comparison C9 forbids, and O3 settles it. Every other
  reader of the stored score is a report; the class-biased sampling of the
  other workflows ranks `ptcl2D` scores.
- C15. `simple_binoris%read_particle_record` reads records narrower than the
  current width and is pinned by `test_legacy_narrow_particle_records`. It
  rejects wider ones. Projects written between 2026-09-21 and the revert carry
  52-real records, so the reader is kept and extended to read the first
  `N_PTCL_ORIPARAMS` values of a wider record.
- C8. The reference files written for the Cartesian path mirror the
  reprojection-model handoff (the master materializes, the matcher only
  reads), so they are a real handoff and stay; their ownership moves
  (section 6.4).

## 4. Open decisions

| Id | Decision | Recommendation | Blocks |
| --- | --- | --- | --- |
| O2 | The convergence rule for `refine=cont` (section 7.2). | Until a rule has been chosen and tested, early convergence is disabled for `refine=cont`: a run goes to `maxits`. The candidate statistics are logged every iteration so the rule can be chosen from data. | Enabling an early stop; no phase |
| O3 | Ranking of stored scores of mixed origin (section 3, C14). Options: (a) the consumer rescores: before it ranks, `refine3D_states` runs one `refine=eval` pass over every active particle, so the ranked values come from one polar evaluation; (b) store the Cartesian score in a field of its own and leave `corr` polar, which needs a record slot and reverses C14; (c) forbid projection-balanced sampling on a project that has seen a Cartesian pass, which needs a record of that history. | (a). It needs no new field, follows the rule that a score is recalculated by whoever needs it, and also removes the existing mix of scores from different iterations in that ranking. Cost: one full evaluation pass at the start of `refine3D_states`, only when `balance=yes` with fractional updates. | Phase 7 |

## 5. Conventions to mirror from the polar branch

| Aspect | Polar branch | Cartesian today | Target |
| --- | --- | --- | --- |
| Objective selection | `params%cc_objfun` in `gen_objfun_vals` | Fixed in `pose_cont_config`, independent of `objfun` | Dispatch on `params%cc_objfun` |
| `objfun=cc` | Normalized correlation, no sigma2 weighting (`gen_corrs`) | `1 - cc` on sigma-whitened data | Unweighted (C5) |
| `objfun=euclid` | `L = sum_k (k/sigma2_k) sum_p abs(X - CTF R)^2 / wsqsum(X)`; score `exp(-L)` | `0.5 sum abs(r)^2 / sigma2`; not normalized; no score | Same loss and score (C6) |
| Shell weight | `k/sigma2_k`: polar Jacobian times inverse variance | `1/sigma2` per Cartesian pixel | Equivalent measure; pinned by test N14 |
| Stored score | `corr` holds cc or `exp(-L)` | Not written | Written by every Cartesian pass (C14) |
| sigma2 | `esig%calc_sigma2(pftc, ...)` after assignment | Own residual routine plus an outside setter | Sigma owner, from the Cartesian calculator |
| CTF | Particle phase-flipped in `prepimg4align` for `CTFFLAG_YES`; reference times abs(CTF) | Unflipped particle, signed CTF (abs for `CTFFLAG_FLIP`) | Equivalent for both objectives; pinned by test N5 |
| Reference | Prepared volume, pad by 2, KB gather times pf^3, no inverse envelope at the gather | The same | Unchanged (C10) |
| Band limit | `kfromto` adopted from the reprojection-model header | Standalone route skips the adoption | Travels with the Cartesian model |

## 6. Target architecture

### 6.1 Dependency direction

```text
refine3D / refine3D_auto commander
    -> refine3D iteration strategy (shared-memory or distributed)
        -> matcher
            -> strategy3D_cont                       (one particle)
                -> builder%cftc   : cartesianft_calc (references, particles, objective)
                -> cartesianft_pose_opt              (LM, bounds, accept or reject)
                    -> simple_kbinterpol, simple_cartesian_fourier
```

The calculator and optimizer import no commander, matcher, UI, project I/O or
parallelization module.

### 6.2 Classes

`src/main/cftc/simple_cartesianft_calc.f90`, type `cartesianft_calc`, owned by
`builder` as `cftc` (the role `pftc` has):

- state-by-half Cartesian references and their `read`/`write`;
- the prepared particles of the current batch and per-thread scratch;
- objective and gradient for `cc` and `euclid`, dispatched on `cc_objfun`;
- the residual used by the sigma owner.

Built from these routines of the refiner: `load_pose_reference`,
`prepare_pose_particle`, `sample_fourier_with_grad`, `predict_unweighted_pose`,
`shift_normal_terms`, `pose_normal_terms`, `finalize_cc_normal_terms`,
`pose_objective_gradient`, `prepared_pose_sigma_contribution`,
`pose_cont_ctf_value`; and from the adapter: `pose_cont_reference_workspace`,
`load_reference_slot`, `write_pose_cont_reference_artifact`.

`src/main/cftc/simple_cartesianft_pose_opt.f90`, type `cartesianft_pose_opt`
(the role `pftc_shsrch_grad` has):

- the LM solve, damping, cumulative bounds, accept or reject;
- takes and returns `ori`; rotation matrices, cropped-grid shifts, normal
  equations and stage records are private.

Built from `refine_shift_lm`, `refine_pose_lm`, `refine_prepared_pose_lm`,
`right_increment_rotation`, `build_pose_lm_system`, `solve_pose_cholesky`,
`increase_lm_damping`, `rejection_terminal_status`,
`apply_pose_parameter_mask`; and from the adapter:
`run_pose_cont_transaction` with its helpers, the pose, config, limits and
result records, and the `ori` conversions.

`simple_cartesian_fourier` stays in `src/main/interp`. Sources under
`src/main` are globbed recursively, so the new folder needs only `cftc.inf`.

Ownership, lifetime and thread-safety:

| Object | Owner | Lifetime | Written when | Read when |
| --- | --- | --- | --- | --- |
| Calculator references (state by half) | `builder%cftc` | One matcher pass | `read`, before the particle loop, serial | Particle loop, concurrently |
| Calculator prepared particles | `builder%cftc` | One batch | `prepare_batch`, in its own parallel loop, one particle per iteration | Particle loop, concurrently |
| Calculator work images and residual buffers | `builder%cftc`, indexed by thread | One matcher pass | By the thread that owns the index | By that thread only |
| Optimizer | A component of `strategy3D_cont` | One particle: built in `new`, released in `kill` | By the one thread running that particle's `srch` | By that thread only |

- The objective, gradient and prediction methods take the calculator as
  `intent(in)`, so concurrent evaluation on different particles cannot write
  shared state. The refiner's evaluation routines already have this form.
- The residual method writes only the thread-indexed buffers. Today the
  sigma routine allocates four per-shell arrays on every call; those become
  the thread-indexed buffers.
- The optimizer holds fixed-size state only (the 5x5 system, bounds, status);
  it allocates nothing during a solve and nothing is shared between
  optimizers. This mirrors `strategy3D_srch` owning its `pftc_shsrch_grad`
  objects.

### 6.3 Strategy hierarchy

```text
strategy3D                 abstract; spec and the four deferred methods only
  strategy3D_pftc          abstract; owns s (strategy3D_srch)
    greedy, greedy_inpl, greedy_smpl, greedy_sub, shc, shc_smpl,
    snhc_smpl, eval, prob              (nine files change their parent)
  strategy3D_cont          Cartesian transaction for one particle
```

The matcher dereferences `ptr%s` in exactly two places (lines 595 and 602).
Both move: the shift increment is taken from the project orientation before
and after `srch` (`s%prev_shvec` is that stored shift), and the `cont_inpl_*`
fields are written by the search object that produced them.

`strategy3D_cont` exposes `new(params, spec, build)`, `srch(os, ithr)`,
`oris_assign()` and `kill()` and nothing else.

### 6.4 Shared owners

| Concern | Owner after the refactor |
| --- | --- |
| Reference handoff | `materialize_reprojection_model_from_volumes`, `read_reprojection_model`, `remove_ref_section_files` produce, read and remove the Cartesian model as they do the polar one; the file contract is section 6.6 |
| Particle preparation | A Cartesian sibling of `prepimg4align` in `simple_matcher_2Dprep`; `build_batch_particles3D` fills the active calculator from the single read |
| Per-thread storage | `prep_strategy3D` / `clean_strategy3D` |
| sigma2 | `euclid_sigma2` gains a Cartesian `calc_sigma2` variant taking `cartesianft_calc`; `set_particle_contribution` is removed |
| Representation switch | One derived logical in `parameters` (`l_cartesian_refine`), set once from `refine`; no string comparisons downstream |
| Reporting | Per-particle fields read by `simple_convergence`; no stats module, no stats files, no timers in the particle loop |

### 6.5 User-facing surface

- `refine=cont` in `refine3D` and `refine3D_auto` (replaces `refine=pose_cont`).
- `pose_cont=yes|no` schedules a Cartesian pass after every discrete
  iteration (C19).
- Removed: `pose_cont_mode`, `pose_cont_route`, the program
  `refine3D_pose_cont`, the forced `cc_emit_sigma=yes`.

### 6.6 Cartesian model file contract

Modelled on the polar reprojection model (`reprojection_model_even.bin` and
`_odd.bin`: a four-integer header, then the payload, validated by
`reprojection_model_header_compatible`). To be confirmed at the Phase 4
checkpoint.

| Item | Contract |
| --- | --- |
| Files | `cartesian_model_even.bin`, `cartesian_model_odd.bin`, named by `refine3D_cart_model_fname(half)` in `simple_refine3D_fnames`. They replace the per-state `pose_cont_reference_stateNN_*.mrc` files. |
| Layout | Stream binary. Header: five default integers `[format_version, box_crop, nstates, kfrom, kto]` and one default real `smpd_crop`. Payload: `nstates` real-space volumes of `box_crop^3` default reals, state 1 first, each the prepared reference that the polar branch pads and projects. |
| Version | `format_version = 1`. A reader rejects any other value. |
| Compatibility on read | `box_crop` and `nstates` equal the run's; `smpd_crop` equal within single-precision tolerance; `1 <= kfrom <= kto <= fdim(box_crop) - 1`; the even and odd headers are identical. The checks are a logical function, `cartesian_model_header_compatible`, which the reader calls; on a mismatch the reader stops the run and names the file and the field. |
| Band limit | The reader sets `params%kfromto` and `params%lp` from the header, as `adopt_reprojection_model_range` does for the polar model. No other source of `kfromto` exists in a Cartesian pass. |
| Producer | `materialize_reprojection_model_from_volumes`, through `cartesianft_calc%write`, from the volumes it has just prepared in memory. It writes after preparation has succeeded, replacing the previous files, as it does for the polar model. With `pose_cont=yes` the same call writes both models from the same prepared volumes. |
| Reader | `read_reprojection_model`, through `cartesianft_calc%read`, in the matcher. Shared-memory and distributed runs read the same files. |
| Removal | `remove_ref_section_files`. Nothing else deletes them. |
| Size | `nstates * box_crop^3 * 4` bytes per half (64 MB per state at box 256). |

The alternative is to keep one MRC volume per state and half, which can be
opened in a viewer, with a small header file beside them. It needs a
completeness check across several files where the single file needs none.

## 7. `refine3D_auto` with `refine=cont`

### 7.1 Stages

| Stage of `exec_refine3D_auto` | Today | With `refine=cont` |
| --- | --- | --- |
| Hard default `refine=prob_neigh` | Always set | `cont` is accepted and kept |
| Sampling (`nsample`, `update_frac`, `trail_rec`) | Set from the active particle count | Unchanged (C12) |
| `ref_pose_init=cc` | Optional polar CC pose initialization | Bypassed; the combination is rejected |
| Sigma bootstrap (`calc_pspec`) and startup `reconstruct3D` | Run | Run; neither is polar |
| NU low-pass seeding | Run | Run |
| Registration pass (`refine=greedy`, polar) | On by default (`regpass`) | Bypassed |
| Main loop | `refine=prob_neigh`, `nspace`, `nspace_sub`, `prob_align`, polar reprojection model | `refine=cont`; no projection grid, no probability tables, no polar model |
| Final reconstruction (`calc_final_rec`) | Run | Run |

This is a continuation mode (C11). It does not produce poses; it refines the
ones a discrete workflow left in the project (`abinitio3D`, or a polar
`refine3D_auto` or `refine3D` run). The commander checks this at entry,
before the sigma bootstrap, with the existing seed check (positive state,
half-set assignment, finite Euler angles and shifts for every active
particle), and stops with a message that names the requirement. The UI text
of `refine` on `refine3D_auto` says the same.

`set_refine3D_auto_sampling` stays as it is. Because `refine=cont` is not a
probabilistic mode, the matcher draws the sample itself through
`sample_ptcls4update3D` (update-count sampling under `balance=no`).

### 7.2 Convergence (open, O2)

- The single-state rule in `check_conv3D` is
  `frac_srch%avg > 99 .and. mi_proj > overlap`, gated by `minits`.
  `refine3D_auto` sets `overlap=0.99`.
- Per particle, `mi_proj` is 1 when the symmetry-aware angular distance
  between the previous and the new orientation is at most `angthres_mi_proj`
  (default 2 degrees). It is not tied to the discrete grid.
- A local optimizer reports `frac=100`, as deterministic neighbourhood
  refinement already does, so the `frac` gate is always met.
- Under `refine=cont` the rule therefore reduces to "99% of the updated
  particles moved less than 2 degrees", which a pass that only removes grid
  error could meet at once, and it ignores shifts.
- Precedent for a motion-based rule: 2D `refine=inpl` converges on
  `dist_inpl%avg < 0.5`.

Interim rule until O2 is decided: `check_conv3D` does not declare convergence
under `refine=cont`. A run goes to `maxits`, which `refine3D_auto` derives
from the target number of updates per particle. `minits` would only delay the
faulty decision. With `pose_cont=yes` in a discrete mode the rule is
unchanged, because the polish leaves the discrete search's convergence fields
in place (Phase 9).

Candidates to evaluate on the logged statistics: an angular threshold derived
from the current resolution and mask diameter (the registration pass already
computes `res / (mskdiam/2)`); a rule on mean orientation and shift motion
(`dist`, `shincarg`); the fraction of improved particles falling below a
floor; an FSC plateau.

## 8. Code inventory

### 8.1 Continuous-pose code in `simple_strategy3D_matcher.f90`

Line numbers at `712b4d86`.

| Lines | What | Disposition |
| --- | --- | --- |
| 23-24, 31-38 | Imports from strategy, stats and adapter modules | Keep only `strategy3D_cont` |
| 55-56, 713-714 | `do_pose_cont_polish`, `do_pose_cont_strategy` | Remove |
| 83-88 | Six declarations of adapter and stats types | Remove |
| 115-116, 389-412 | Skip model adoption; seed validation | Model read; strategy or parameters |
| 136, 139-140 | Empty-partition sigma and stats | Generic sigma prep; stats removed |
| 164, 170, 242 | Skip `memoize_refs`, `prep_strategy3D`, `clean_strategy3D` | Behind `prep_strategy3D` and builder |
| 182-183, 226-231, 240 | Per-thread stats allocate, merge, report | Remove |
| 188-191 | Clear `pose_cont_*` fields on every run | Remove with the fields (C15) |
| 211 | Skip PFTC `calc_sigma2` | Sigma owner |
| 221 | `frac_greedy` guard | Strategy-neutral counting |
| 246-251, 270 | Teardown variants | Generic teardown |
| 345-371 | Mode validation in `init_ctrl` | `parameters` validation |
| 447-468 | Reference and image allocation branch | `read_reprojection_model`, `alloc_ptcl_imgs` |
| 477-486 | Batch-build branch | `build_batch_particles3D` |
| 494-501, 555-593 | Downcasts, `bind_context`, `get_result`, timers, field writes | Remove; the strategy writes its own fields |
| 595, 602 | `ptr%s` dereferences | Section 6.3 |
| 610-669 | `run_pose_cont_after_pftc` | Remove; the polish is a scheduled pass (C8) |
| 671-678 | `wall_time_seconds` | Remove |

### 8.2 Exports of `simple_pose_cont_refine3D_adapter.f90`

"Users" are production modules outside the adapter.

| Export | What it is | Users | New owner |
| --- | --- | --- | --- |
| `pose_cont_reference_workspace` | State-by-half array of refiners; methods choose even or odd and forward | matcher, strategy | `cartesianft_calc` |
| `cartesian_pose_data` | Re-export of the refiner's prepared-particle type | matcher, strategy | `cartesianft_calc`, private |
| `pose_cont_particle_workspace` | Raw-particle copies kept across PFTC preparation | matcher, `ptcl_batch` | Removed; a pass has one representation |
| `prepare_pose_cont_observation` | Noise-normalise, clip, mask, expand | strategy | `simple_matcher_2Dprep` |
| `pose_cont_observation_spec`, `pose_cont_observation`, `pose_cont_particle_spec` | Argument bundles | matcher, strategy | Removed with the wrappers |
| `pose_cont_pose` | Rotation matrix and shift | matcher, strategy | `cartesianft_pose_opt`, private |
| `pose_cont_seed_from_orientation`, `pose_cont_pose_to_orientation` | `ori` to matrix and crop-scaled shift, and back | matcher, strategy | `cartesianft_pose_opt`, private |
| `shift_native_to_crop`, `shift_crop_to_native` | Multiply by the crop factor | none | Folded into the conversion |
| `pose_cont_config`, `pose_cont_limits`, their constructors, `POSE_CONT_ROUTE_*` | LM policy | matcher, strategy, stats, `refine3D_strategy` | Private constants of the optimizer (C16) |
| `POSE_CONT_OBJECTIVE_*` | Objective selector independent of `objfun` | stats | Removed (C4) |
| `pose_cont_transaction_result`, `pose_cont_stage_result`, status codes | Result and accounting of one solve | matcher, strategy, stats | `cartesianft_pose_opt`, private |
| `run_pose_cont_transaction` (private) | Staged route, cumulative bounds, accept or reject | through the workspace | `cartesianft_pose_opt` |
| `write_pose_cont_reference_artifact` | Validated `image%write` | `refvol_utils` | `cartesianft_calc%write`, called by the model materializer |
| `remove_pose_cont_reference_artifacts` | Delete the files | none | `remove_ref_section_files` |
| `load_reference_slot` (private) | Read, check dimensions and sampling, build refiner | through the workspace | `cartesianft_calc%read` |

### 8.3 File disposition

Paths are under `src/main/` unless they start with `src/` or `doc/`.

| File | Disposition | Phase |
| --- | --- | --- |
| `commanders/simple/simple_commander_refine3D_pose_cont.f90` | Delete | 1 |
| `strategies/parallelization/simple_refine3D_pose_cont_workflow.f90` and tester | Delete | 1 |
| `ui/simple/simple_ui_refine3D_pose_cont.f90` | Delete | 1 |
| `doc/policies/3D/refine3D_pose_cont_policy.md` | Delete | 1 |
| `commanders/simple/simple_commanders_refine3D.f90` | Remove the four hooks (1); accept `refine=cont` (8) | 1, 8 |
| `volume/simple_cartesian_pose_refiner.f90` and tester | Split into the two `cftc` classes; delete after parity | 2, 7 |
| new `cftc/simple_cartesianft_calc.f90`, `cftc/simple_cartesianft_pose_opt.f90`, testers, `cftc.inf` | Create | 2 |
| `simple_builder.f90` | Add `cftc` | 4 |
| `strategies/search/simple_matcher_refvol_utils.f90` | Own the Cartesian model beside the polar one; drop adapter import and `cartesian_only` | 4 |
| `strategies/search/simple_strategy3D.f90` | Representation-neutral | 5 |
| new `strategies/search/simple_strategy3D_pftc.f90` | Owns `s` | 5 |
| nine polar strategy files | `extends(strategy3D_pftc)` | 5 |
| `strategies/search/simple_matcher_ptcl_batch.f90`, `simple_matcher_2Dprep.f90` | Cartesian preparation sibling; drop adapter import and `optional` flags | 6 |
| `sigma2/simple_euclid_sigma2.f90` | Cartesian `calc_sigma2`; remove the setter | 6 |
| `strategies/search/simple_strategy3D_alloc.f90` | Skip PFTC allocations when `l_cartesian_refine` | 6 |
| `strategies/search/simple_strategy3D_pose_cont.f90` and tester | Replace with `simple_strategy3D_cont.f90` and tester | 7 |
| `strategies/search/simple_strategy3D_matcher.f90` | Remove every entry of section 8.1 | 5, 7 |
| `strategies/search/simple_pose_cont_refine3D_adapter.f90` and tester | Delete | 7 |
| `strategies/search/simple_pose_cont_run_stats.f90` and tester | Delete | 7 |
| `strategies/parallelization/simple_refine3D_strategy.f90` | Remove stats aggregation (7); bypass `prob_align` and the polar model for `refine=cont` (8); schedule the polish (9) | 7, 8, 9 |
| `src/defs/simple_defs_ori.f90`, `src/fileio/simple_binoris.f90`, `simple_convergence.f90` | Record width 50, slot 42, tolerant reader, statistics (C15) | 7 |
| `params/simple_parameters*.f90` | Objective override removed (3); derived flag (4); surface (8) | 3, 4, 8 |
| `ui/simple/simple_ui_refine3D.f90`, `exec/simple_exec_refine3D.f90` | Remove the program (1); `refine=cont`, `pose_cont`, and `refine` on `refine3D_auto` (8) | 1, 8 |
| `commanders/test/simple_commanders_test_class.f90`, `ui/simple_test/simple_test_ui_class.f90` | Suite table and `suite=` list follow every tester move | each |
| `strategies/search/simple_pose_cont_1jyx_tester.f90` | Port to the new classes, then to `strategy3D_cont` | 2, 7 |

## 9. Phases

Each phase is one reviewable change that compiles and passes the fast gate.
Stop for maintainer review at every phase boundary. The tests named in each
phase are defined in section 10.

### Phase 0. Baseline

- Maintainer records, at `712b4d86`: results and per-sub-suite times of
  `unit_cart_align3D` and `lib_cart_align3D`; one `refine3D refine=pose_cont`
  iteration and one `pose_cont=yes` iteration on the beta-gal set (FSC, mean
  motion, accept fraction, peak RSS, wall time) for both LM routes; one
  ordinary `refine3D_auto` run with continuous refinement off; the nightly
  `simulated_workflow_1jxy` and `abinitio3D_addon` metrics.
- From the `lib_cart_align3D` run (E25), record the metrics that later serve
  as pre-refactor floors: rotation and shift error before and after, accept
  fraction, mean objective before and after, and the truth-map correlation of
  the perturbed and refined reconstructions. Where the tester computes a
  metric without printing it, a print line is added; no assertion changes.
- Record the files each run leaves in its directory.

Exit: numbers and file inventories pasted into section 14.

### Phase 1. Remove the separate program (C13)

- Delete the commander, the workflow module and its tester, the UI file, the
  exec case and the policy file.
- Remove `pose_cont_mode` from `parameters`, the parser and validation,
  including the program-name checks in `parameters_phases` (lines 853-863)
  and the `prg` case at line 821.
- Remove the four hooks from `commander_refine3D_auto` and restore the inline
  defaults, registration-pass setup and main-stage call.
- Tests: retire E22-E24 and E26 (section 10.2).

Exit: `exec_refine3D_auto` differs from `e86286f2^` only by unrelated later
commits; the ordinary `refine3D_auto` baseline reproduces; registry check
clean.

### Phase 2. Numerical owners (C3, C17)

- Create `cartesianft_calc` and `cartesianft_pose_opt` by moving code without
  changing formulas, defaults or the current objective. Both LM routes stay.
- Private components, `new`/`kill`, type-bound methods. Replace `error stop`
  with `THROW_HARD` for caller-contract violations and typed statuses for
  per-particle failures.
- Keep `simple_cartesian_pose_refiner` and the adapter compiling as thin
  forwarders so production is untouched; they gain no behavior.
- Tests: the refiner and adapter testers keep running unchanged against the
  forwarders, which is the parity proof. Their checks (E4-E19) are copied to
  the calculator and optimizer testers and run against the new interfaces;
  port E25; add N1 and N2.

Exit: every Phase 0 numerical result reproduces within the existing
tolerances through both the old and the new testers; neither class imports
matcher, commander, UI or project modules.

### Phase 3. Objectives follow `objfun` (C4, C5, C6)

The scientific change, reviewed on its own.

- Dispatch on `params%cc_objfun`; implement the `cc` and `euclid` definitions
  of section 5; remove `POSE_CONT_OBJECTIVE_*` and the `cc_emit_sigma`
  override in `parameters_phases`.
- Tests: rewrite E8 and E15 for the new definitions; retire E20; add N3-N7;
  run E25 under each objective, with floors from the injected error and the
  Phase 0 values (section 10.1), not from its own first run.

Exit: both objectives pass formula, gradient and recovery tests; N6 shows the
Cartesian and polar objectives select the same grid pose.

### Phase 4. Builder and reference handoff (C10)

- Add `builder%cftc`; derive `l_cartesian_refine` once in `parameters`.
- Implement the file contract of section 6.6: `write`/`read` on the
  calculator, a logical header-compatibility function, the calls from the
  model materializer and reader, the filename function, and removal in
  `remove_ref_section_files`.
- Delete `new_inverse_envelope_reference`.
- Tests: retire E6; keep E10 as the guard on the gather; add N8 and N9.

Exit: one owner produces, validates, reads and removes the references; a
distributed worker builds the calculator from declared files only; no file
I/O in the optimizer.

### Phase 5. Representation-neutral strategy base (C17)

- Move `s` from `strategy3D` to a new abstract `strategy3D_pftc`; change the
  parent of the nine polar strategies.
- Matcher: shift increment from the project orientation before and after
  `srch`; `cont_inpl_*` written by `strategy3D_srch`.
- Tests: `unit_pftc_align2D3D` unchanged; add N10; rerun the Phase 0 polar
  guards (section 10.5).

Exit: ordinary `refine3D` and `refine3D_auto` baselines reproduce; the matcher
contains no `ptr%s`.

### Phase 6. Shared helpers select the representation

- `build_batch_particles3D` prepares the active calculator's particles from
  the single read; remove `build_batch_particles3D_cartesian`, the particle
  workspace argument and the `optional` pad-image arguments.
- `prep_sigmas_objfun` and `euclid_sigma2`: Cartesian `calc_sigma2`; remove
  `cartesian_only`.
- `prep_strategy3D` / `clean_strategy3D` and the builder toolbox allocate only
  what the active representation needs.
- Tests: move E13 (observation part) and E18 to their new owners; add
  N11-N14 and N26.

Exit: no shared helper imports a continuous-pose module or takes a feature
flag; one sigma path; N14 passes.

### Phase 7. `strategy3D_cont` and matcher cut-over (C1, C2, C14, C15)

- Implement `strategy3D_cont` on the Phase 2-6 objects: seed from the project
  field, one solve against the particle's state and half reference, commit
  pose and Cartesian score together, write the Cartesian score at the seed
  when the solve is rejected, sigma through the sigma owner, convergence
  fields (`dist`, `dist_inpl`, `shincarg`, `mi_proj`, `frac`).
- Matcher: the strategy case and nothing else; remove every entry of
  section 8.1. The in-matcher polish is unavailable from here until Phase 9.
- Delete the adapter, the old strategy, the stats module, the forwarders of
  Phase 2, the master-side stats aggregation, and their testers.
- Record width back to 50; improved flag in slot 42; `read_particle_record`
  reads the first 50 values of a wider record; `simple_convergence` derives
  attempted from `sampled`.
- Apply the ruling on O3 before the first Cartesian score is written to
  `corr`. Under the recommended option, `prepare_refine3D_states_class_sampling`
  runs a `refine=eval` pass before it ranks.
- Tests: delete the old refiner and adapter testers with the forwarders
  (their checks already run in the calculator and optimizer testers); retire
  E21; move the seed-validity check of E19 to the `strategy3D_cont` tester;
  port E25 to `strategy3D_cont`; add N15-N19, N24 and N27.

Exit: the matcher's only continuous-pose symbol is `strategy3D_cont`; no
`select type`; ordinary refinement baselines reproduce; no consumer ranks
stored scores of mixed origin.

### Phase 8. `refine=cont` surface and `refine3D_auto` (C7, C11, C12, C16, C18, C20, C21)

- UI and parameters: the refine value becomes `cont`; `refine` is exposed on
  `refine3D_auto` with `prob_neigh|cont`; `pose_cont_route` is removed and the
  kept LM route becomes private; validation in
  `validate_parameter_consistency` (`oritype=ptcl3D`, initialized poses,
  `inpl_cont` interaction, C20, C21).
- `simple_refine3D_strategy`: for `refine=cont`, no `prob_align` child and no
  polar model.
- `exec_refine3D_auto`: section 7.1, including the entry check that makes it
  a continuation mode.
- Convergence: no early stop under `refine=cont` (section 7.2); mean
  orientation motion, mean shift motion, improved fraction and FSC are logged
  every iteration.
- Tests: reduce the route tests E16 and E17 to the kept route; add N20
  (stages a-e and h), N21 and N25.

Exit: `refine3D_auto refine=cont` completes from a project whose poses came
from a polar run, runs to `maxits`, updates `nsample` particles per
iteration, and creates no polar model, probability table or registration
pass; shared-memory and distributed runs agree.

### Phase 9. Polish as a scheduled pass (C8, C9, C19)

- `pose_cont=yes` makes the iteration strategy follow every discrete pass
  with a `refine=cont` pass over the same particle sample
  (`sample4update_reprod`), in `refine3D` and in the main loop of
  `refine3D_auto`.
- One reconstruction per iteration, from the polished poses: the discrete
  pass writes no partial reconstructions when a polish follows.
- The polish pass writes pose, `corr`, the improved flag and sigma. It leaves
  the convergence fields of the discrete search (`dist`, `dist_inpl`,
  `mi_proj`, `frac`) in place, so the main-loop convergence rule keeps
  measuring the discrete search and is not met trivially by the small polish
  motion.
- Tests: add N22 and N23; extend N20 with stages f and g.

Exit: with `pose_cont=no` every baseline reproduces; with `pose_cont=yes` the
pose error against the simulation truth is no worse than without it, and a
discrete iteration that follows a polish runs unaffected.

### Phase 10. Documentation and closing checks

- Add the continuous-pose section to `doc/policies/3D/refine3D_policy.md` and
  `refine3D_auto_policy.md`.
- Update `doc/policies/test_environment_policy.md`: the sub-suite table of
  section 1.1, the long-running table of 1.2, the process budget, and the
  "Where is my test now?" table.
- Regenerate `doc/code_overview/fortran-indexes` once, here.
- `SIMPLE_UNIT_ORDER=reverse simple_test_exec test=unit_cart_align3D`.

Exit: stale-symbol scan empty for `pose_cont_refine3D_adapter`,
`pose_cont_run_stats`, `strategy3D_pose_cont`, `cartesian_pose_refiner`,
`pose_cont_mode`, `pose_cont_route`, `refine3D_pose_cont`; this file moves to
`doc/refactoring_notes/completed/`.

## 10. Test plan

### 10.1 Principles

From `doc/policies/test_environment_policy.md`:

- A test asserts through `simple_test_utils`, derives its expected value
  independently of the code under test, names its tolerances, fixes its seeds,
  and pins sign, handedness and index conventions.
- A unit test is a routine in `simple_<thing>_tester.f90` beside the code it
  tests, exporting one `run_all_<thing>_tests`, registered as a sub-suite in
  `simple_commanders_test_class.f90` and listed in the `suite=` help of
  `simple_test_ui_class.f90`.
- Fast tier: hermetic, in-process, one thread, each sub-suite well under a
  second, whole gate under 30 s. Library tier: the same rules at realistic
  size. High-level tier: a pipeline on simulated data, gated against the
  simulation truth with declared floors and a `metrics.tsv`.
- A new CTest entry is an owner decision; the one new entry of this plan is
  accepted (C18).

Specific to this refactoring:

- An existing test is moved, not rewritten, when its guarantee still holds.
  It moves in the same commit as the code it tests.
- While a forwarder exists (Phases 2-6), the old tester keeps running against
  it. That is the parity proof for the move.
- A test is retired only when its subject is deleted or its guarantee is
  replaced by a decision in section 3. The commit names the reason and the
  replacing test.
- A new behavior is not complete until its test is in the same phase. No
  phase exits on manual inspection alone.
- Cross-representation tests compare optima and physical quantities, never
  raw objective values (C9).
- A floor is never taken from the first run of the code being accepted. It
  comes from the simulation truth (the error must fall by a stated factor
  relative to the injected error), from an analytic bound (the rotation error
  must end below the basin width at the band limit, `lp / (mskdiam/2)`
  radians), or from the pre-refactor implementation measured at Phase 0 on
  the same fixture, with a stated margin.
- A guarantee that ends in `THROW_HARD` cannot be asserted in-process. Such
  checks are written as logical functions that the caller turns into a stop,
  and the function is what the test asserts.
- Fixtures stay at the sizes already in use (box 16 and 24 in the fast tier)
  and any file is named `tmp_<tester>_*` and removed, also on failure.

### 10.2 Existing tests and their disposition

Fast tier, `unit_cart_align3D`, unless noted.

| Id | Tester and routine | Pins | Disposition | Phase |
| --- | --- | --- | --- | --- |
| E1 | `cartesian_fourier`: KB derivative (fast polynomial, normalized stencil, stencil switch) | KB kernel and its derivative | Keep in place | - |
| E2 | `cartesian_fourier`: packed gather derivative | Gather gradient against finite differences | Keep in place | - |
| E3 | `cartesian_fourier`: neutral extract (embed/crop, crop envelope, packed gathers, plane extraction) | Parity with the pre-extraction operations | Keep in place | - |
| E4 | `pose_refiner`: `test_particle_preparation_contract` | Prepared-particle validity, shell capping | Move to the calculator tester | 2 |
| E5 | `pose_refiner`: `test_sigma_shell_contract` | Shell whitening inputs | Move to the calculator tester | 2 |
| E6 | `pose_refiner`: `test_reference_envelope_contract` | Inverse-envelope constructor | Retire with the constructor; E10 guards the gather | 4 |
| E7 | `pose_refiner`: `test_shift_phase_sign` | Shift sign on the native pixel scale | Move to the calculator tester | 2 |
| E8 | `pose_refiner`: `test_ncc_objective_formula` | 1-NCC formula, gain invariance | Move (2); rewrite for unweighted cc as N3 (3) | 2, 3 |
| E9 | `pose_refiner`: `test_five_parameter_gradient` | Gradients of both objectives | Move (2); extend as N5 (3) | 2, 3 |
| E10 | `pose_refiner`: `test_matched_projector_boundary` | Gather equals `interp_fcomp_oversamp` | Move to the calculator tester; permanent guard for C10 | 2 |
| E11 | `pose_refiner`: `test_rotation_increment` | Right increment keeps orthogonality | Move to the optimizer tester | 2 |
| E12 | `pose_refiner`: `test_shift_solver`, `test_joint_solver`, `test_tiny_accepted_reduction_stop` | Recovery, masks, guards, stopping | Move to the optimizer tester | 2 |
| E13 | `pose_adapter`: `test_observation_and_coordinate_adapters` | Observation equals the established particle path; crop scaling | Observation part to the calculator tester (6); scaling part to the optimizer tester (2) | 2, 6 |
| E14 | `pose_refiner`: `test_invalid_and_unobservable_inputs` | Invalid input leaves the pose untouched | Move to the optimizer tester | 2 |
| E15 | `pose_refiner`: `test_ncc_solver` | NCC solve on a gain-scaled particle | Move (2); rerun under the new cc (3) | 2, 3 |
| E16 | `pose_adapter`: `test_route_transaction_contracts` | Stage accounting of both routes | Move to the optimizer tester (2); reduce to the kept route (8) | 2, 8 |
| E17 | `pose_adapter`: `test_rollback_transaction_contracts` | Bound rejection, no-improvement, invalid preparation preserve the pose | Move to the optimizer tester | 2 |
| E18 | `pose_adapter`: `test_sigma_endpoint_contracts` | Sigma contribution on accept and rollback | Calculator tester (2); sigma-owner test N13 (6) | 2, 6 |
| E19 | `pose_adapter`: `test_reference_workspace_lifecycle`, `test_inpl_pose_cont_handoff`; `pose_strategy`: `test_seed_contract` | Even/odd reference lifecycle; seed round trip with metadata intact; seed validity | Lifecycle to the calculator tester (2); round trip to the optimizer tester (2); seed validity to the `strategy3D_cont` tester (7) | 2, 7 |
| E20 | `pose_strategy`: `test_outer_objective_sigma_policy` | The predicate behind the forced `cc_emit_sigma` | Retire with the predicate; replaced by N3 and N13 | 3 |
| E21 | `pose_statistics`: `test_partition_aggregation` | Stats-file aggregation | Retire with the module; replaced by N19 | 7 |
| E22-E24 | `pose_workflow`: `test_registration_stage_policy`, `test_post_matcher_stage_policy`, `test_standalone_final_stage_policy` | Child command lines of the removed program | Retire with the module; replaced by N20 and N21 | 1 |
| E25 | `lib_cart_align3D`: `pose_cont_1jyx` (library tier) | 5,000 simulated 1JYX particles perturbed by 15 degrees and 2 pixels: objective, rotation error and shift error fall; refined map beats the perturbed one | Keep; port to the new classes (2), both objectives with declared floors (3), through `strategy3D_cont` (7) | 2, 3, 7 |
| E26 | `unit_ui`, UI visibility: `test_refine3D_pose_cont_policy` | Surface of the removed program | Retire; replaced by N21 | 1 |
| E27 | `unit_project`, binoris: `test_legacy_narrow_particle_records`, `test_particle_segment_roundtrip` | Narrow legacy records; record is `N_PTCL_ORIPARAMS` reals | Keep; joined by N17 | - |
| E28 | `unit_numerics`, KB kernel: `test_apod_mat_3d_fast_grad` | The stencil derivative | Keep in place | - |

### 10.3 New tests

Each is implemented in the phase shown, in the tester named, with the
expected value from the source named.

| Id | Phase | Tier and tester | Pins | Expected value from |
| --- | --- | --- | --- | --- |
| N1 | 2 | fast, `cartesianft_calc` | `new`/`kill` repeated; even and odd references built from different volumes give different predictions for the same pose | The two volumes, gathered in the test |
| N2 | 2 | fast, `cartesianft_pose_opt` | `ori` in, `ori` out: state, half, `proj` and every non-pose field unchanged by a solve | The input record |
| N3 | 3 | fast, calculator | `cc` formula; invariant to particle gain; unchanged when sigma2 is altered (negative control for C5) | Closed form over the test arrays |
| N4 | 3 | fast, calculator | `euclid` loss and `exp(-L)` score; invariant to a uniform scaling of sigma2; changes under a shell-dependent scaling | Brute-force sum in the test |
| N5 | 3 | fast, calculator | Analytic gradients of both objectives, without CTF, with CTF and with phase flip | Central differences |
| N6 | 3 | fast, calculator with a small `polarft_calc` | For a noise-free particle placed on a grid pose, the Cartesian objective evaluated over the grid poses peaks at the pose where `gen_objfun_vals` peaks, for both objectives | The placed pose |
| N7 | 3 | fast, optimizer | Recovery of injected shift and rotation under each objective | The injected pose |
| N8 | 4 | fast, calculator | Model file of section 6.6: write, read, identical prediction for every state and half; header carries version, box, state count, `kfrom`, `kto` and sampling; the reader adopts `kfromto`; the compatibility function is false for a wrong version, box, state count, sampling, an out-of-range band, and even and odd headers that differ | The written arrays; hand-built headers |
| N9 | 4 | fast, calculator | `remove_ref_section_files` removes the Cartesian model; no file survives the test | Directory listing |
| N10 | 5 | fast, `unit_pftc_align2D3D`, refine3D in-plane state | The search object's previous shift equals the shift stored in the project before the search, so the matcher's increment from the project orientation is the one formerly read from `s%prev_shvec` | The stored shift |
| N11 | 6 | fast, calculator | Batch preparation equals the established path (`norm_noise_fft_clip_shift`, `ifft_mask_fft`, `expand_ft`) particle by particle | That path, called in the test |
| N12 | 6 | fast, calculator | Threaded batch preparation (team of three) equals the serial result | The serial result |
| N13 | 6 | fast, calculator with `euclid_sigma2` | Cartesian `calc_sigma2` at a committed pose equals the per-shell residual; an exact match gives zero | Brute-force per-shell sum |
| N14 | 6 | fast, calculator with a small `polarft_calc` | Per-shell sigma2 from the Cartesian residual and from the polar residual agree for a noise-only particle at a grid pose | Agreement within the standard error of a shell |
| N15 | 7 | fast, `strategy3D_cont` | Lifecycle through a `class(strategy3D)` pointer: `new`, `srch`, `kill`, repeated; identity seed accepted; state-zero particle rejected; runs with `build%pftc` never constructed | Construction of the fixture |
| N16 | 7 | fast, `strategy3D_cont` | Accepted solve commits pose, `corr`, `dist`, `shincarg`, `mi_proj`; rejected solve leaves the pose bit-identical and writes the Cartesian score at the seed | Seed record; objective evaluated in the test |
| N17 | 7 | fast, `unit_project`, binoris | A 52-real record file reads back with its first 50 values; narrow records still read | Hand-written payload |
| N18 | 7 | fast, `unit_ori` | The improved flag lives in slot 42: set, get, record round trip, flag name | The value set |
| N19 | 7 | fast, `pose statistics` sub-suite on `simple_convergence` | Attempts equal the sampled count; improved percentage from slot 42 | Hand count over a small `oris` |
| N20 | 8, 9 | high-level, new entry `cont_refine3D_1jxy` (C18) | Stages below | Simulation truth |
| N21 | 8 | fast, `unit_ui`, UI visibility | `refine3D` offers `refine=cont` and `pose_cont`; `refine3D_auto` offers `refine` with `prob_neigh` and `cont`; `pose_cont_route`, `pose_cont_mode` and the program `refine3D_pose_cont` are absent | The contract of section 6.5 |
| N22 | 9 | fast, `strategy3D_cont` | A particle not in the pass sample is untouched: pose, `corr` and flags unchanged | The input record |
| N23 | 9 | fast, `strategy3D_cont` | In a polish pass the convergence fields written by the discrete search are unchanged; in a `refine=cont` pass they are written from the seed-to-result motion | The input record; motion computed in the test |
| N24 | 7 | fast, `unit_pftc_align2D3D`, refine3D in-plane state | After an `eval` pass over particles whose stored scores were overwritten with arbitrary values, every active particle's `corr` is the polar objective at its stored pose (the rescoring that O3 option (a) relies on) | `gen_objfun_vals` evaluated in the test |
| N25 | 8 | fast, `pose statistics` sub-suite on `simple_convergence` | Under `refine=cont`, convergence is not declared for any overlap and search fraction, and the motion statistics are still computed; in a discrete mode the rule is unchanged | The rule of section 7.2 |
| N26 | 6 | fast, calculator | Objective, gradient and residual evaluated concurrently on different particles by a team of three equal the serial results bit for bit | The serial results |
| N27 | 7 | fast, `strategy3D_cont` | A batch run through `srch` by a team of three gives the same poses, scores and flags as the serial run | The serial run |

Stages of the high-level gate N20. Fixture: the recipe of E25 (the embedded
1JYX model, the same box, sampling, CTF spread, noise level and seeds), so
that the Phase 0 numbers of the pre-refactor implementation apply to it; the
project is seeded with the truth poses perturbed by the same known rotation
and shift. Metrics go through `test_gate` to `metrics.tsv`. Floors follow
section 10.1: reduction relative to the injected error, the analytic bound,
and no worse than the Phase 0 value within a stated margin.

| Stage | Run | Gate |
| --- | --- | --- |
| a | `refine3D refine=cont objfun=euclid`, shared-memory | Pose error against the truth (`pair_pose_error`) below its floor and below the perturbed input; map against the truth (`validate_reconstructed_volume`) |
| b | The same with `objfun=cc` | The same |
| c | Stage a with `nparts=2` | Poses and map agree with stage a within a declared tolerance |
| d | All of a-c | No polar reprojection model and no probability table in the run directory; only declared files remain |
| e | `refine3D_auto refine=cont` from the perturbed project | Runs to `maxits`; no pose-initialization or registration pass ran; `nsample` particles updated per iteration; pose error below its floor |
| f (Phase 9) | `refine3D` discrete mode with `pose_cont=yes` and with `pose_cont=no` | Pose error with the polish no worse than without; the discrete iteration after a polish runs and its scores are finite |
| g (Phase 9) | `refine3D_auto pose_cont=yes` from the perturbed project | Completes; a Cartesian pass follows every main-loop iteration over the same sample; one reconstruction per iteration; pose error no worse than stage f without the polish |
| h | The normal entry path: a polar `refine3D_auto` run from the perturbed project, then `refine3D_auto refine=cont` on its output project | The continuation accepts the polar run's poses without any seeding step; pose error and map against the truth are no worse after the continuation than after the polar run |

### 10.4 Suite layout after the refactoring

| Entry | Sub-suites, in table order | Tester modules |
| --- | --- | --- |
| `unit_cart_align3D` | `Cartesian Fourier`, `Cartesian calculator`, `pose optimizer`, `pose strategy`, `pose statistics` | `simple_cartesian_fourier_tester`, `simple_cartesianft_calc_tester`, `simple_cartesianft_pose_opt_tester`, `simple_strategy3D_cont_tester`, `simple_convergence_tester` |
| `lib_cart_align3D` | `pose 1JYX recovery` | `simple_strategy3D_cont_1jyx_tester` |
| `cont_refine3D_1jxy` (C18) | high-level gate N20 | commander in `simple_commanders_test_highlevel.f90` |

Removed sub-suites: `pose refiner` and `pose adapter` (Phase 7, split into the
calculator and optimizer sub-suites), `pose workflow` (Phase 1). The
`suite=` selectors become `cartesian_fourier`, `cartesian_calculator`,
`pose_optimizer`, `pose_strategy`, `pose_statistics`.

### 10.5 Guards on the polar code

The refactoring touches code the polar path shares (Phases 5-7). These
existing suites are rerun and compared with Phase 0 at each of those phases:

- fast: `unit_pftc_align2D3D`, `unit_ori`, `unit_project`, `unit_ui`,
  `unit_reconstruction`;
- high-level: `simulated_workflow_1jxy`, `simulated_workflow_6vxx` and
  `abinitio3D_addon`, which drive the polar strategies through the matcher;
- by hand: the ordinary `refine3D_auto` run of Phase 0.

### 10.6 Budget and records

- The fast gate's time for `unit_cart_align3D` is read from
  `build/test_runs/ctest_fast.log.timing.txt` at Phase 0 and after every
  phase. A sub-suite that approaches a second is reduced or moved to the
  library tier before the phase exits.
- During Phases 2-6 some checks run twice (old tester on the forwarder, new
  tester on the class). That duplication ends in Phase 7.
- `SIMPLE_CTEST_BUDGET` goes from 29 to 30 in Phase 8, with the new entry (C18).
- Every retirement in section 10.2 is stated in its commit message with the
  id of the replacing test.

## 11. Checks at every phase

- `git diff --check`
- `python3 scripts/check_test_registry.py . --verbose`
- `python3 scripts/check_descr.py .`
- no multiline `THROW_HARD`/`THROW_WARN`
- stale-symbol scan for whatever the phase moved or deleted
- `./compile_debug.sh` (build and fast gate), run by the maintainer

## 12. Risks and controls

| Risk | Control |
| --- | --- |
| Numerical drift during moves | Phase 2 moves formulas unchanged; old and new testers both pass before Phase 3 changes anything |
| Regression in the nine polar strategies | Phase 5 is mechanical and stands alone; section 10.5 guards rerun before the Cartesian child is added |
| Two reference authorities | One materializer and filename family; validity check on read; N8 and N9 |
| Shared-memory and distributed divergence | Same child command and file contract; N20 stage c |
| sigma2 scale differs between representations | N14, before any pass mixes them (Phase 9) |
| Memory and I/O | One representation per pass; with `pose_cont=yes` every iteration reads its particle sample twice; peak RSS, stack reads and wall time per iteration measured against Phase 0 |
| Project compatibility | Record width 50; tolerant reader; N17 |
| Fast gate grows | Section 10.6 |
| Scores of mixed origin are ranked | O3; N24; Phase 7 exit |
| A Cartesian run stops before it has stabilized | No early stop until O2 is decided; N25 |
| A floor canonizes a regression | Floors from truth, an analytic bound and the Phase 0 pre-refactor values (section 10.1) |
| Data race in the particle loop | Ownership table of section 6.2; N12, N26, N27; the nightly four-thread run of the fast tier |

## 13. Review checkpoints

1. Calculator and optimizer public interfaces, before production callers move (Phase 2).
2. Objective definitions and their tests (Phase 3).
3. Reference and particle lifecycle, and the file contract of section 6.6 (end of Phase 4).
4. Polar baselines after the hierarchy split (end of Phase 5).
5. O3, before Phase 7; matcher contents (end of Phase 7).
6. User-facing surface and the `refine3D_auto` stage bypass (end of Phase 8).
7. Polish scheduling, reconstruction and convergence-field handling (Phase 9).

## 14. Progress

| Date | Phase | Commit | Evidence | Status |
| --- | --- | --- | --- | --- |
| 2026-09-30 | Plan | - | Contract C1-C21 accepted; test plan added | Awaiting Phase 0 |
| 2026-09-30 | Plan review | - | Six review findings incorporated: mixed-origin scores (O3), no early stop under `refine=cont` (O2), continuation-only `refine3D_auto` mode, model file contract (6.6), ownership and thread-safety (6.2), floors (10.1) | O3 awaiting ruling |
