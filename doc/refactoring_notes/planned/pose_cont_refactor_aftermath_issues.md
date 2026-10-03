# pose_cont refactoring: aftermath issues

Status: open, 2026-10-02. Collects the maintainer's comments on the delegated run and the open
questions of [the review report](../completed/pose_cont_refactoring_report.md), with the code
facts each one rests on. The completed plan is
[pose_cont_refactoring.md](../completed/pose_cont_refactoring.md); its contract ids (C*, O*) and
test ids (E*, N*) are used here. Line references are to the working tree of 2026-10-02 and
will drift.

| # | Issue | State |
| --- | --- | --- |
| 1 | `inpl_cont` before the Cartesian polish | Behaviour already as ruled; pin it |
| 2 | Sigma2 preparation before the polish | Resolved at review |
| 3 | Assembly reads the iteration's realized sample (Cartesian fix, polar lag) | Open: one focused correction |
| 4 | Convergence rule under `refine=cont` (O2) | Open: design below |
| 5 | Open review questions of the report | Open: listed |

## 1. `inpl_cont` before the Cartesian polish

Ruling (maintainer): no polar in-plane refinement in a pure `refine=cont` pass. In a polar
neighbourhood refinement with the Cartesian polish (`pose_cont=yes`), `inpl_cont` follows the
polar refinement: when it is on, the polar continuous in-plane step runs in the discrete pass
and precedes the polish.

The code already does this. A Cartesian pass overrides `inpl_cont` to `no`
(`simple_parameters_phases.f90`, the `l_cart_refine` block), which applies to `refine=cont`
and to the polish pass, and nothing else touches it: `check_polish_request` and
`polish_follows` do not read it, `refine3D_auto` never sets it (its `pose_cont=yes` main loop is
`prob_neigh` with the user's `inpl_cont`, default `yes`; `pose_cont=only` is `refine=cont` and
gets `no`). The discrete pass of a polish iteration keeps the discrete parameters (shared
memory: `cline_build`; distributed: the job description with only `volrec=no` changed). The
in-plane step runs in `refine_selected_continuously` (search-loop modes) or
`refine_assignment_continuously` (prob modes) and `assign_ori` stores the refined `e3` and
shift; the project goes to disk (shared memory) or is merged (distributed) before the polish
reads it as its seed (`strategy3D_cont`, seed from the stored pose).

Work:

- State the ruling in `refine3D_policy.md` (polish paragraph) and `refine3D_auto_policy.md`.
- Pin it with a test: a polish iteration with `inpl_cont=yes` reports continuous in-plane
  attempts in the discrete pass and the polish seed equals the stored, in-plane-refined pose;
  a `refine=cont` pass reports none. N20 stage f or g is the natural carrier (convergence
  statistics already print `CONTINUOUS IN-PLANE RUNS IMPROVED`).

## 2. Sigma2 preparation before the polish

Review comment: "Prepare sigma state again and make the polish contributions authoritative.
Do not accept as currently implemented. It loses the freshly computed discrete fallback."

Resolved on 2026-10-02 (plan section 15, Review row; N30). The polish no longer opens a second
transaction: it starts from the discrete pass's range of its part
(`euclid_sigma2%read_pending_range`), overlays the contributions of particles whose Cartesian
slot is valid, and replaces the range atomically, so a particle the polish cannot evaluate keeps
the residual of its discrete pose. Left here for the record and to close with the next build
(N30 passed on the first run).

## 3. Assembly reads the iteration's realized sample

Two review items are one problem. The run made a shared-memory Cartesian pass write its sample
to the project before assembly, because trailing reconstruction reads it from the file; the
report also noted that polar non-probabilistic shared-memory modes assemble with the previous
iteration's sample. The ownership rule they both violate: assembly must consume the current
iteration's realized sample, whatever the representation and execution mode
(`importance_sampling_fractional_update_policy.md`: counts "from the merged project";
`sampled` and `updatecnt` consistent before trailing consumes them; shared memory and
distributed follow the same workflow).

Where the sample comes from and where it goes:

- Drawn by `sample_particles_for_update` (`simple_strategy3D_matcher.f90`): reproduced for prob
  modes and the polish (`sample4update_reprod`), drawn fresh otherwise, which updates `sampled`
  and `updatecnt` in memory. The matcher writes `algndoc`, not the project.
- Gridding assembly (`exec_volassemble`, `simple_commanders_rec_distr.f90`) builds its own
  builder from disk; `determine_trailing_update_fraction` reads N, n, the realized fraction
  f = n/N and the new-row count from the project file and weights the trailing blend with them.
  The partial reconstructions themselves come from the current in-memory sample, so in a lagging
  path the partials and their weights disagree.

Which paths write the sample before assembly:

| Path | Before assembly |
| --- | --- |
| Distributed, any mode | Yes: `merge_algndocs` writes the merged project |
| Shared memory, prob modes | Yes: the `prob_align` child samples and writes; the parent reads it back |
| Shared memory, `refine=cont` | Yes: the run's fix |
| Shared memory, any polish iteration | Yes: the writes around the polish |
| Shared memory, non-prob polar, no polish | No: written only in `inmem_finalize_iteration`, after assembly |

Affected: shared memory, gridding backend, trailing reconstruction on (`update_frac <= 0.99`
with `trail_rec=yes`), no polish, a matcher-drawn sample: `shc`, `shc_smpl`, `snhc_smpl`,
`neigh`, `greedy`, `greedy_inpl`, including the fill-in and `update_missing` routes. Not
affected: distributed runs, the PCG backend (it assembles from the in-memory builder),
`refine3D_auto` (its main loop is `prob_neigh`), polish iterations. The error is largest while
`updatecnt > 0` coverage is still growing.

Correction (one focused change in `simple_refine3D_strategy.f90`, shared memory): after
`refine3D_exec` and before assembly, write the particle segment whenever partial
reconstructions are written, for every mode, replacing the Cartesian special case
(`if( params%l_cart_refine .and. l_write_partial_recs )`) with the general rule
(`if( l_write_partial_recs )`); the writes the polish needs stay. Distributed needs nothing.
State the contract in the fractional-update policy and the `simple-frac-update-trailing` skill:
the iteration's realized sample is on disk before assembly in both execution modes.

Test: shared-memory `refine3D refine=shc update_frac=0.5 trail_rec=yes` for two iterations on a
small simulated set; after each assembly the f, N, n and new-row count `volassemble` logs equal
the in-memory sample of that iteration (they equal the previous iteration's today). Compare with
the same run distributed, which is correct now.

## 4. Convergence rule under `refine=cont` (O2)

Interim rule (in place): no convergence is declared under `refine=cont`; the run goes to
`maxits`. This is the safe choice. The discrete rule (`check_conv3D`: search-space fraction
above 99 % and projection overlap above 0.99, `mi_proj` from `angthres_mi_proj`, default 2
degrees; shifts logged, not used) would stop a continuous run as soon as the small grid error
is removed and ignores shifts. A polish iteration is judged by the discrete rule on the
discrete pass's fields, which is correct and stays.

Design (maintainer):

- Statistics over the iteration's sample: symmetry-aware angular motion of the committed pose
  (seed to result), shift motion, and the improvement fraction (the improved flag), all of
  which a Cartesian pass already records.
- Converged when, for at least two consecutive iterations, the mean (or a high quantile of the)
  angular motion is below an angular tolerance tied to the current resolution and mask radius,
  the shift motion below a shift tolerance, and the improvement fraction below a floor.
- Angular tolerance from the current resolution and mask radius: the registration pass already
  uses `res / (mskdiam/2)` (radians) as its basin width; a fraction of it is the natural
  tolerance. The shift tolerance follows the same resolution (a fraction of `res/smpd` pixels).
- Evaluated separately for full and fractional updates: with `update_frac < 1` the statistics
  cover only the sampled particles, so the consecutive-iteration condition must count
  iterations per coverage (or require the condition over enough iterations to visit every
  particle once), not raw iterations.

Open points for the design: which quantile; how the tolerances scale (fractions to fix on the
logged statistics of N20 stage e and the beta-gal `pose_cont=only` runs before any floor is
declared); `minits`; whether `refine3D_auto pose_cont=only` keeps its update budget as the upper
bound. The decision belongs in the plan's O2 and the policy once settled.

## 5. Open review questions of the report

- Optimizer API: whether `refine_shift`, `refine_joint`, `refine_pose` and the optional
  arguments of `cartft_pose_opt%new` (step and iteration caps, `shift_first`) stay public; they
  are tester seams today, production passes only the bounds.
- The joint LM route (C16 ruling): measure its beta-gal timing and objective in production
  (`pose_cont=only` and the polish) against the Phase 0 numbers of both routes.
- Rerun the polar guards (`simulated_workflow_1jxy`, `_6vxx`, `solve3D_addon`) on the merged
  tree: two commits after `714a15a99` changed the simulated-workflow validation the run used.
- Pending runs of the review changes: `unit_cart_align3D` (N31 model-power gate, N33 group sigma2, `test_rotation_bound`,
  E16 both routes), `unit_ui` (N21), `lib_cart_align3D` (E25), `cont_refine3D_1jxy` (N20; stages
  e-h are the first workflow runs of the joint route).
- Comment pass over the run's new code (`cftc/`, `strategy3D_cont`, testers) under the comment
  rule of `simple-modern-fortran`; several file headers exceed it.
- Regenerate `doc/code_overview/fortran-indexes`; update the skills `simple-refine3d` and
  `simple-main-strategies` for `refine=cont`, the polish, `athres_cont` and `cont_route`.
- Every adopted recommendation (O4-O8, the 6.6 file layout) remains open to reversal; O6 and C16
  were revised by ruling on 2026-10-02.
