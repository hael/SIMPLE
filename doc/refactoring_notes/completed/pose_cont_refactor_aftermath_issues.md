# `pose_cont` refactoring aftermath

Status: active, reviewed against the source at `108c8c9da` (2026-10-07). This is the single living
record of what the completed [`pose_cont` refactoring](../completed/pose_cont_refactoring.md) and its
[review report](../completed/pose_cont_refactoring_report.md) left to do. Line numbers are at
`108c8c9da` and will drift; the routine names are the stable reference.

## Work list

State on 2026-10-07: everything below except V1 is implemented in the working tree on top of
`108c8c9da`, uncommitted and not yet compiled or run.

| # | Change | Kind | State |
| --- | --- | --- | --- |
| C1 | Publish the particle field before file-backed assembly | Correctness fix | Implemented |
| C2 | Make `params_polish` an allocatable component | Compile-time policy | Implemented |
| C3 | Shorten five file headers | Comment rule | Implemented |
| C4 | Convergence rule for `refine=cont` (motion of most particles) | Feature | Implemented |
| T1 | Pin the trailing counts of a shared-memory fractional run (tests C1) | Test | Implemented |
| T2 | Pin the `inpl_cont` phase derivation and the polar in-plane route before a polish | Test | Implemented |
| D1 | Correct the route sentence; add the publication invariant and the C4 rule | Policy | Implemented |
| D2 | Add `refine=cont` and the publication invariant to three skills | Skills | Implemented |
| V1 | Compile, run the gates, record the results | Validation | Open |

The sigma2 fallback and the LM route are closed (end of this note).

## C1. Publish the particle field before file-backed assembly

**Where.** `inmem_execute_iteration` in
`src/main/strategies/parallelization/simple_refine3D_strategy.f90`, the line after the main
`refine3D_exec` call (line 760):

```fortran
call refine3D_exec(params, build, cline, params%which_iter, converged, l_write_partial_recs .and. .not. l_polish)
! a Cartesian pass draws its particle sample in the matcher (no prob_align child writes it),
! and the assembly reads the sample of the iteration (trailing update fractions) from the
! project file, so the matcher's field goes to disk before the assembly
if( params%l_cart_refine .and. l_write_partial_recs ) call build%spproj%write_segment_inside(params%oritype)
```

**What is wrong.** In shared memory the matcher changes `sampled`, `updatecnt`, pose and state in
`build%spproj_field` and writes only `algndoc.simple`. For a discrete pass without a polish the
project file is written by `inmem_finalize_iteration`, after the assembly. The gridding
`volassemble` builds its own toolbox from the project file, and `determine_trailing_update_fraction`
(`simple_commanders_rec_distr.f90:982`) rereads `ptcl3D` and takes `N`, `n`, the realized
fraction and the first-time rows from `get_group_update_counts` and `get_state_update_fracs`. It
therefore sees the previous iteration's sample. Read from the source, not observed in a run:

- with `trail_rec=yes`, the first iteration of a fresh run stops: `inmem_initialize` clears
  `sampled` and `updatecnt` at `startit=1` and writes the project, so `get_group_update_counts`
  finds no sample and calls `THROW_HARD('requires previous sampling')`;
- in later iterations the trailing blend uses the previous iteration's `N` and `n`;
- with or without trailing, `refresh_state_populations` counts the previous iteration's states;
  those counts select the NU matching low-pass in multi-state runs
  (`update_project_nu_alignment_lowpass`).

Affected: shared memory, `rec_backend=gridding`, `volrec=yes`, no polish, every discrete mode except
`sigma`. In the probabilistic modes the `prob_align` child publishes the sample, but the poses and
states the parent changes still reach disk late. Not affected: the distributed path
(`merge_algndocs` writes the project before `volassemble`), the polish branch (it writes before and
after the polish pass), and PCG assembly (`assemble_refine3D_pcg` reads the in-memory builder).

**Change.** Replace the line and the three comment lines above it with:

```fortran
if( publish_before_file_assembly(params, l_write_partial_recs, l_polish) ) &
    &call build%spproj%write_segment_inside(params%oritype)
```

and add, beside `polish_follows`:

```fortran
!> gridding volassemble rereads the project file, so the field of this iteration goes to disk first;
!! the polish branch publishes its own; a sigma pass writes no orientations
logical function publish_before_file_assembly( params, l_write_partial_recs, l_polish ) result( l_publish )
    type(parameters), intent(in) :: params
    logical,          intent(in) :: l_write_partial_recs, l_polish
    l_publish = l_write_partial_recs .and. .not. l_polish .and. trim(params%refine) /= 'sigma' .and. &
        &(params%l_cart_refine .or. trim(params%rec_backend) == 'gridding')
end function publish_before_file_assembly
```

- The `l_cart_refine` term keeps today's write for Cartesian passes on PCG. It is redundant there,
  since PCG assembly uses the in-memory builder, but removing it is a separate change.
- The `sigma` exclusion matches the matcher (`ctrl%do_write_oris` is false in sigma mode) and
  `inmem_finalize_iteration`.
- Leave both writes in the polish branch (lines 765 and 778). Do not move a bare
  `if( l_write_partial_recs )` write above the polish branch: it would duplicate the pre-polish
  write.

**Done when** T1 passes and the fractional-update policy states the invariant (D1).

## C2. `params_polish` as an allocatable component

**Where.** `refine3D_inmem_strategy`, `simple_refine3D_strategy.f90:67`:
`type(parameters) :: params_polish`. `doc/policies/compile_time_policy.md` ("Do not embed heavy
types inline") requires such a component to be allocatable.

**Change.**

- Declare `type(parameters), allocatable :: params_polish`, as the stream stages do.
- In the polish branch of `inmem_execute_iteration`, add
  `if( .not. allocated(self%params_polish) ) allocate(self%params_polish)` before
  `init_params_and_build_strategy3D_tbox(cline_polish, self%params_polish)`.
- In `inmem_cleanup`, now a no-op, add
  `if( allocated(self%params_polish) ) deallocate(self%params_polish)`.

Keep it a component, not an iteration-local variable. The toolbox built with it by
`init_params_and_build_strategy3D_tbox` lives until the next iteration rebuilds it, so the object
must too.

## C3. File headers

The rule (`.github/skills/simple-modern-fortran/SKILL.md`): one `!@descr:` line, then at most five
comment lines. Current header lengths after the `!@descr:` line:

| File | Lines |
| --- | --- |
| `src/main/cftc/simple_cartft_calc.f90` | 31 |
| `src/main/cftc/simple_cartft_pose_opt.f90` | 27 |
| `src/main/strategies/search/simple_strategy3D_cont_1jyx_tester.f90` | 27 |
| `src/main/strategies/search/simple_strategy3D_cont_tester.f90` | 13 |
| `src/main/strategies/search/simple_strategy3D_cont.f90` | 11 |

Cut each to five lines. Move contract text that is still true to `refine3D_policy.md` (the
Cartesian section) and drop plan-section and ruling references. The two `cftc` testers already
comply.

## C4. Convergence rule for `refine=cont`

**Before.** `check_conv3D` never declared convergence for a pure `refine=cont` pass, so the run
went to `maxits`. A polish iteration (`pose_cont=yes`) uses the discrete pass's rule and is not part
of this change.

**Rule.** A refine=cont run has converged when most particles moved by at most half a degree and at
most 1 A between consecutive poses. Made exact:

- **Particle.** A particle of this pass's sample (`sampled` at its maximum, `state > 0`) is stable
  when all three hold:
  - its transaction ended `CARTFT_ACCEPTED` or `CARTFT_NO_IMPROVEMENT`;
  - `dist + dist_inpl <= 0.5` degrees;
  - `shincarg * smpd <= 1` A.
- **Iteration.** The iteration passes when at least 90% of the sampled particles of every populated
  state are stable.
- **Run.** The run converges when consecutive passing iterations have together sampled at least 90%
  of the active particles: every active particle with `sampled >= g0`, where `g0` is the sample
  generation at which the streak began. Under full sampling one passing iteration is enough. A
  failing iteration resets the streak. `minits` still applies, and the rule never extends `maxits`
  or the `refine3D_auto` update budget.

Notes on the rule:

- **Rotation measure.** `dist + dist_inpl` (projection-direction plus in-plane angle, both
  symmetry-aware from `sym_dists`) bounds the rotation angle from above for small angles. The rule
  is therefore conservative and needs no new rotation field.
- **Shift measure.** `shincarg` is the shift change in native pixels, so the 1 A bound converts
  through `params%smpd`, not `smpd_crop`.
- **Why the status is needed.** A rejected transaction keeps the seed pose, so invalid
  preparation, an unreliable update, bound rejection, invalid numerics and the iteration limit all
  show zero motion. Without the status, a pass that failed on most particles would look converged.
  Counting those outcomes as unstable also caps failures at 10%.
- **Why per state.** States do not change under `refine=cont`, but a small state could hide behind
  a large one in a global fraction.
- **For comparison.** 0.5 degrees moves a point at a 100 A radius by 0.87 A; the two bounds are
  matched for a particle of about 230 A diameter.

**Code (as implemented).**

1. **Typed outcome per particle.** Name the spare slot 42 `pose_cont_status`:
   - in `src/defs/simple_defs_ori.f90`, add `I_POSE_CONT_STATUS = 42`, a branch in
     `get_oriparam_ind` and `get_oriparam_flag`, drop `I_SPARE_ORIPARAM42`, and leave only slots
     above `I_LAST_NAMED_ORIPARAM` in `oriparam_is_spare`;
   - in `src/main/ori/simple_ori.f90:1059-1060` and `:1081-1082`, print slots `1..I_LAST_NAMED_ORIPARAM`;
   - in `src/main/project/simple_binoris_tester.f90:195`, slot 42 is now named.

   The record width does not change and projects written earlier read 0 there
   (`CARTFT_NOT_ATTEMPTED`). The distributed path merges the slot with the rest of the record.
2. **Write it.** In `oris_assign_cont` (`simple_strategy3D_cont.f90`), set
   `pose_cont_status = real(self%status)` for every searched particle, in pure and polish passes
   alike.
3. **Decide.**
   - Constants in `simple_defs_conv.f90`: `CONT_CONV_ROT_DEG = 0.5`, `CONT_CONV_SHIFT_A = 1.0`,
     `CONT_CONV_FRAC = 0.9`.
   - In `simple_convergence.f90`: the elemental `cont_particle_stable(status, rot_deg, shift_a)`
     (public, for the tester), the component `cont_streak_g0` (0 when no streak runs), and the
     method `check_cont_conv`, called from `check_conv3D` under `l_cart_refine` in place of the
     interim block.
   - It logs the stable percentage per state (multi-state), of the least stable state, and the
     streak coverage. It writes `CONT_STABLE_PCT` and `CONT_STREAK_COVERAGE_PCT` to the stats
     file, and `get` returns both.

   The shared and distributed strategies both call `check_conv3D` on the merged field, so the rule
   is the same. The distributed strategy keeps its existing extra condition
   `iter >= startit + 2`; the shared one has none.
4. **Tests.**
   - `simple_convergence_tester.f90` (N25): `test_cont_particle_stable` (inclusive bounds, every
     `CARTFT_*` outcome) and `test_cont_rule_full_sample`, on a 10-particle fixture:
     - all stable and 9 of 10 converge, 8 of 10 do not;
     - all no-improvement at zero motion converges, all bound-rejected does not;
     - 3 px fails at `smpd` 1 and passes at 0.3;
     - a state at 60% vetoes;
     - `minits` holds.
   - `test_cont_rule_streak`: half samples complete the streak in two passing iterations, and a
     failing iteration resets it.
   - `test_cartesian_pass_never_converges` became `test_cartesian_pass_statistics`. A field without
     a status does not converge, and the discrete and polish checks are unchanged.
   - `simple_strategy3D_cont_tester`: the accepted solve records `CARTFT_ACCEPTED`, the rejected
     one a rejection code.
   - `simple_binoris_tester`: slot 42 is named and round-trips.
5. **N20 stage e.** With `maxits` on the command line, `refine3D_auto` sets `minits = maxits`, so
   `e_ran_to_maxits` still holds. The stage now reports `e_cont_stable_pct` from the stats file.
6. **Policy.** The rule replaces the interim text in `refine3D_policy.md`. Section 7 of
   `refine3D_auto_policy.md` states that an automatic `maxits` can stop before the update budget.

**Calibration.** The bounds are absolute and stated in physical units, so nothing has to be fitted
before landing. Read the logged stable fraction per iteration in N20 and in a beta-gal
`pose_cont=only` run. If runs stop early with visibly worse FSC, raise `CONT_CONV_FRAC` before
changing the bounds.

## T1. Trailing counts of a shared-memory fractional run

**Owner.** N20 stage f, `run_refine3D_discrete` in
`src/main/commanders/test/simple_commanders_test_highlevel_cont.f90`. It already runs the affected
path: shared memory, gridding, `refine=neigh`, `update_frac=0.3`, two iterations, with
`pose_cont=no` and `pose_cont=yes`. No new CTest entry.

**Change.**

1. Set `trail_rec=yes` in `run_refine3D_discrete`.
2. After each of the two runs, still in the stage directory (before `simple_chdir(gate_root)`), read
   `refine3D_trail_manifest_fname(1)` with `trail_chain_manifest%read`, and take `nrep` from the final
   project with `out%os_ptcl3D%get_group_update_counts('state', 1, nrep, nsmp)`.
3. Check `nint(manifest%get_mrep()) == nrep(1)` through the existing `gate%check`, named
   `f_<tag>_trailing_counts_from_current_sample`.

Why this is exact: with `ufrac_trec` unset the applied fraction equals the realized one, and
`population_blend_weights` then returns `mnew = N` (`simple_oris_sampling.f90:81-102`). The
manifest of the last assembly must therefore carry the `N` of the final project. Iteration 2 of a
0.3 sample has a larger `N` than iteration 1, so a stale count differs.

On the current ordering the `pose_cont=no` run stops at the first assembly (C1, first consequence),
and the `pose_cont=yes` run passes. After C1 both pass.

## T2. `inpl_cont` derivation and the polar in-plane route before a polish

What is already pinned: the internal default and its absence from every UI
(`test_inpl_cont_policy`), the `strategy3D_cont` refusal of `inpl_cont=yes`, and the seed read
from the stored pose (`test_seed_contract`). The polish seeds itself from the project file the
pre-polish write leaves, so the handoff needs no further assertion.

**Change.**

- **Phase derivation.** In `simple_strategy3D_inplane_tester.f90`, add `test_inpl_cont_derivation`:
  build `parameters` through `params%new(cline)` from a minimal command line (`smpd`, `box`,
  `mskdiam`, `ctf=no`, `objfun=cc`, `nthr=1`, as `simple_cavg_registration_tester` does) and assert:
  - `refine=cont`: `inpl_cont == 'no'`, `l_cart_refine`, not `l_cont_polish`;
  - `refine=cont pose_cont=yes`: `l_cont_polish`;
  - `refine=shc`: `inpl_cont == 'yes'`, not `l_cart_refine`.
- **Route before a polish.** In N20, after stage f with `pose_cont=yes`, check that some particle
  of the last sample has `cont_inpl_attempted > 0`. The polish pass runs with `inpl_cont=no` and
  does not clear that field, so a non-zero value comes from the discrete pass. After stage a (pure
  `refine=cont`), check that no particle has `cont_inpl_attempted > 0`.

## D1. Policy

- `doc/policies/3D/refine3D_policy.md:499` says the optimizer runs a "shift stage then joint stage;
  the only route, C16". Replace it with: one bounded LM transaction per particle, the joint
  five-parameter stage; `cont_route=shift_then_joint` is an internal test seam. The paragraph at
  lines 518-523 is already correct.
- With C1, add to `doc/policies/importance_sampling_fractional_update_policy.md`: before an assembly
  that rereads the project file, the iteration strategy writes the particle field the matcher left;
  pose, state, `sampled`, `updatecnt` and the partial reconstructions describe one iteration.

## D2. Skills

None of these mention `refine=cont`, `pose_cont`, `athres_cont` or `cont_route` today.

- `.github/skills/simple-refine3d/SKILL.md` and `.github/skills/simple-main-strategies/SKILL.md`:
  `refine=cont` (Cartesian pass, `strategy3D_cont`, `src/main/cftc`), the polish
  (`pose_cont=yes`, scheduled by the iteration strategy over the discrete sample), `athres_cont`
  as its own rotation bound, and `cont_route` as a test-only seam.
- `.github/skills/simple-frac-update-trailing/SKILL.md`: the C1 invariant.

## V1. Runs to record

For each, record the command, the commit and the result in this note. First compile the working
tree (`./compile_debug.sh`) and run the fast gate.

- Fast gate (`units`), including `unit_cart_align3D`. The last observed pass was at `1aaa42ba`,
  before the review.
- E25 (`lib_cart_align3D`).
- N20 (`cont_refine3D_1jxy`), after T1 and T2.
- The polar guards (1JYX, 6VXX). 6VXX failed once from de-novo variability and passed on rerun;
  keep the `solve3D_addon` guard open until a current pass is recorded.

## Closed

- **Sigma2 fallback during a polish.** The polish opens no second sigma2 transaction. It overlays
  the residuals of valid Cartesian slots on the discrete pass's pending range, so a particle the
  polish cannot evaluate keeps its discrete residual. Pinned by N30 (`test_polish_sigma_fallback`).
- **LM route.** Production runs the joint five-parameter stage through `cartft_pose_opt%refine`.
  `refine_pose`, `refine_shift`, `refine_joint` and the step and iteration overrides stay public
  for the testers only. They have no production caller at `108c8c9da`; do not add one. Narrow
  them when the testers move to a deliberate test interface.

The note moves to `completed/` when V1 is recorded.
