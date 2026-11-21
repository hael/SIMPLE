# automsk=order: the order weight as a reference mask on both backends

Status, 2026-11-21: implemented locally; real-data validation outstanding. Rulings: section 3.
Base: `f9fd7ca21` (tree clean). Supersedes:
[the PCG order prior](../rejected/pcg_order_prior.md), moved to `rejected/` during Phase 4. Depends on:
the NU evidence margin and label field (validated: yes, in production since `7e81790da`); the density
envelope (validated: yes, production); the order weight and core guard (validated: offline only, six pairs,
2026-11-07/08, under the order prior; maps approved by Hans); a benefit of attenuation to alignment
(validated: no; Phase 7 is that test).

## Terms

- **Half maps, base pair.** The even and odd reconstructions from the two halves of the particle set; the
  base pair is the pair no prior has shaped. It gives the FSC (Fourier shell correlation).
- **Gridding, PCG.** The two reconstruction backends. PCG (preconditioned conjugate gradients) solves the
  reconstruction as a linear system and can solve on a *support*, a real-space region outside which the
  solution is held at zero.
- **NU filtering.** Nonuniform filtering: a bank of low-pass candidates (20 Å down to 4 Å, capped at
  FSC=0.143/1.5) competes voxel by voxel on cross-half prediction. The winning *label field* gives the
  `_nu_filt` references (each half filtered label by label) and `_nu_locres` (the label as a resolution in
  Å). It runs per state in `nonuniform_filter_state` (`simple_nu_state_filter`), which gridding
  `volassemble` and the PCG master both call.
- **Evidence margin.** Per voxel, the gain of the best finer candidate over the 20 Å candidate in cross-half
  prediction, smoothed at `amsklp` (8 Å). About zero on solvent and on a detergent micelle, positive on
  ordered protein.
- **Density envelope `E`.** `automask3D` of the pair average: low-pass at `envmsklp` (20 Å), Otsu, largest
  component, dilation of at least 7.5 Å, cosine skirt. It keeps every density present at 20 Å, the micelle
  included. File `automask3D_stateNN.mrc`. Its *core* is the binary before dilation.
- **Order weight `w_f`.** A logistic of the evidence margin whose threshold is the `(1-f)` quantile of the
  margin over the dense set and whose width is the spread of the margin below the threshold: 1 on ordered
  density, 0 on solvent and disordered density. `f` is the *ordered fraction* of the dense set.
- **Core, density around the core.** The core is the dense voxels (inside the envelope's core) whose NU label
  is at or finer than 1.5 × the pair's FSC=0.5 resolution; the density around the core is the dense voxels
  outside it, on a membrane protein mostly the micelle.
- **Kept signal `K`.** A region's cross-half signal `sum y_e y_o`, in which the halves' independent noise
  cancels, after the mask over before it.
- **Lag one.** Made in iteration N-1, used in iteration N.

## 1. Why

Hans retired the PCG order prior (`pcg_order=yes`) on 2026-11-21: in the loop it did not do its job. Its
record (`pcg_order_prior.md`, sections 11-13) says where it failed: the ridge in the PCG operator is
nonlocal (the operator removed 2 to 4.5 times the core signal the closed form predicted on the membrane
pairs), it acts as a local low-pass where the data precision falls with resolution, it lags one iteration
(off at the first iteration and at every new box), and it costs a second cold PCG solve per half, state and
iteration. Those are exactly the parts that tied it to PCG. Its order weight and its core guard need no
backend, and `automsk=nu` already has the plumbing a backend-free mode needs: it works on the live
competition inside the one state filter both backends call, and applies its mask to the current references.
`automsk=order` keeps that plumbing and changes only the mask. The comparison of the two workflows at the
base is in the appendix.

## 2. Contract

- **The mask.** `A(r) = E(r) · M(r)`, `M(r) = 1 / (1 + λ (1 - w_f(r)))` with `λ = 1`: a voxel without evidence
  of order keeps half its amplitude, an ordered voxel all of it. `A` multiplies the composed `_nu_filt` pair
  in `nonuniform_filter_state`, after the skirt median of `E` has been removed. Nothing else in the loop
  sees `M`.
- **Order is density plus one multiplication.** Under `automsk=order` the NU competition, its noise level,
  the label field, the handoff, `_nu_locres`, the background clamp (density), the PCG solve support
  (density envelope of the lag-one prior-free pair), the FSC (gridding: density, phase-randomized; PCG:
  on its support) and every published map other than the `_nu_filt` products are those of
  `automsk=density` on the same build.
- **Off-state invariant.** With `automsk=no|density|nu` no line of the order-mask producer runs, and every
  output is byte-identical to the base at one thread (with `yes` spelled `density`).
- **Backends and modes.** Gridding and PCG, shared-memory and distributed, run the same code (the state
  filter) and must behave the same. Both backends define the solvent region the same way (decision 14).
- **No default changes.** `automsk=order` is opt-in everywhere. The only default value that changes is the
  spelling `yes` → `density` where `yes` is the default today (`exec_refine3D_workflow`, which serves
  `refine3D_auto` and `refine3D_local`, and `estimate_diam`).

## 3. Decisions and rulings

1. (2026-11-21, Hans) The PCG order prior is retired: the keys `pcg_order` and `pcg_order_lambda`, the program
   `order_prior` and the code behind them are deleted, with no alias.
2. (2026-11-21, Hans) `automsk=no|density|nu|order`. `yes` is renamed `density` everywhere, the standalone
   mask programs included (`no|density|tight` there); `yes` is rejected, with no alias. The `automsk` field
   (`character(len=5)` at the base) is widened, since `density` would be cut to `densi`.
3. (2026-11-21, Hans) Depth `λ = 1`: `M` runs from 1/2 to 1. Under the Euclidean objective a micelle kept at
   1/2 leaves 25 % of the image/reference mismatch that full removal causes.
4. (2026-11-21, Hans; amended the same day by decision 14) The quantile domain is the dense set: the
   density envelope's core, restricted to the NU domain, on both backends. The first reading of the
   ruling was the `mskdiam` sphere ("as it is done now"). On PCG the sphere outside the solve support holds
   exact zero/zero voxels, which the NU code treats as boundary conditions and keeps out of every statistic
   (`simple_nu_filter_evidence`: "a degenerate cluster at gap ~0 that would otherwise pin the lower
   quartile"); on gridding the same voxels hold solvent noise. A quantile and a width taken over the sphere
   would therefore draw the line between ordered and disordered density in different places on the two
   backends. The dense set lies at least the dilation width (7.5 Å or more) inside the PCG support, so both
   backends reconstruct every voxel of it.
5. (2026-11-21, Hans: "a much safer core") The core's resolution limit is 1.5 × the pair's FSC=0.5
   resolution (the order prior used FSC=0.143). `nonuniform_filter_state` therefore takes the pair's FSC=0.5
   resolution as a new argument; both callers hold it (`res05s(state)`) before they call it.
6. (2026-11-21, Hans) `solve3D` accepts any `automsk` on its command line and forwards it, with
   `envfsc=yes`, to every stage from stage 7 on; stages 1-6 keep `automsk=no envfsc=no`. Without `automsk`
   on the command line solve3D is unchanged.
7. (2026-11-21, Hans) solve3D's final reconstruction ships with the `automsk` the run used.
8. (2026-11-21, Hans) `nu_filt3D` gains `automsk=order` and replaces `prg=order_prior` as the offline check.
9. (2026-11-21) The ordered fraction: the smallest `f` on the grid 0.01, 0.02, …, 0.99 whose threshold is
   positive and whose `K_core(f) = sum_core M_f² y_e y_o / sum_core y_e y_o >= 0.95`, with `y_e`, `y_o` the
   composed halves after the skirt median of `E` is removed. The order prior's guard (0.95, ruled 2026-11-08)
   is kept; the strength is not scanned (decision 3 fixes it). The scan reads the core voxels once per `f`
   and uses no transform.
10. (2026-11-21) Off this iteration: fewer than 1000 core voxels, `sum_core y_e y_o <= 0`, or no `f`
    qualifies. Then `A = E` (the references equal those of `automsk=density`), the log says why, and the
    stale order artifact is deleted.
11. (2026-11-21) The background level always comes from `E`'s skirt, in the state filter and in the matcher
    fallback: `zero_env_background` with `E`, then a multiplication by `A`. `zero_env_background` takes the
    median over the voxels strictly between 0 and 1, which on `A` would be the attenuated interior.
12. (2026-11-21) The consumers:

    | consumer | `density` | `nu` | `order` |
    | --- | --- | --- | --- |
    | PCG solve support | density envelope, lag-one prior-free pair | same | same |
    | NU background clamp | density | NU envelope when valid | density |
    | `_nu_filt` references | × density | × NU envelope when valid | × `A`, current iteration |
    | gridding FSC mask | density, phase-randomized | lag-one NU, phase-randomized | density, phase-randomized |
    | PCG FSC | on the support, no post-hoc mask | same | same |
    | matcher mask while NU references are missing | (non-NU filtering: density) | lag-one NU, then density | lag-one `A` after the density background, then density |
    | final reconstruction, registered mask | density | NU, then density | `A`, then density |
    | final reconstruction, `filt_mode` | `none` | NU mode kept | NU mode kept |
    | refine3D provenance list | density artifact | + NU artifact | + order artifact |

13. (2026-11-21) One log line per state and iteration, `>>> NU ORDER MASK: STATE n, ...`: `f`, the threshold
    and width, `K_core`, the share of the dense signal around the core removed (at most 75 % at `λ = 1`),
    the kept signal of the dense voxels outside the core whose label is finer than the coarsest bank member,
    and the core's share of the dense set; or `no attenuation this iteration` with the reason.
14. (2026-11-21, Hans: "We should be consistent in how we define the solvent region between PCG and
    gridding. That's my only constraint.") Under `automsk=order` the solvent region is the density envelope
    (`automask3D` at `envmsklp`, `binwidth`, `edge`) in every consumer on both backends, and the order mask
    reads nothing that one backend defines differently from the other:
    - `E`, the dense set and the core come from the current base pair in `nonuniform_filter_state`, the same
      code on both backends;
    - the quantile, the width, `K_core` and the share removed are taken over the dense set only
      (decision 4), so PCG's zero/zero voxels never enter them;
    - the evidence margin enters only through quantile ranks and the spread below the threshold, so the
      global noise level by which each backend normalises it (`nu_objective_noise_scale`, over observed
      voxels) cancels;
    - the NU-evidence envelope, whose null differs by backend (the dilation ring on a constrained PCG pair,
      the solvent majority of the sphere on gridding), is written as a diagnostic and used by nothing under
      `order`.

    One difference stays, inherited from the standing `envfsc` contract and not changed here: PCG must hold
    its support before it solves, so its FSC (and the FSC=0.5 resolution of decision 5) is taken on the
    lag-one density envelope inside the estimator, while gridding applies the current density envelope post
    hoc with the phase-randomized correction. Both use the same envelope construction. Phase 0 measures how
    far apart the two backends' envelopes and FSC=0.5 resolutions are on the same particles and
    orientations, and Phase 2 reports the same for the order masks (section 6).
15. (2026-11-21, Hans, after review: "take the final reconstruction's mask settings from the last stage's
    command line only"; refines decision 7) The final reconstruction takes `automsk`, `envfsc` (given, or
    implied by an active `automsk`), the NU `filt_mode` under `automsk=nu|order`, `nu_msk_sig` and
    `amsklp` from the command line of the refinement it continues from (`final_rec_mask_settings` in
    `simple_final_rec`), for both the child command line and the registered mask, never from the calling
    workflow's parameters. Under solve3D that is the last stage's command line, so a run whose last stage is
    before stage 7 (multi-state at its default of five stages, any `nstages` up to 6) ships without the
    requested mask and with `envfsc=no`; solve3D warns at entry and runs `automsk=no` (9c17b8ca8). Before the ruling the final
    reconstruction of such a run mixed solve3D's parameters (`envfsc=yes`, NU `filt_mode`) with the last
    stage's `automsk=no`. `envmsklp` and `binwidth`, the envelope's construction parameters, still come from
    the caller's parameters; the stage command line is a copy of the caller's, so the values are the same.

Defaults taken without a question:

- The artifact is `order_envmask3D_stateNN.mrc` (`ORDER_ENVMASK_FBODY = 'order_envmask3D_state'` in
  `simple_defs_automask`), overwritten every iteration. The dedicated module follows the repository's
  rule against adding workflow constants to definitions re-exported by the core API.
- The constants live beside `NU_ENVMASK_*` in `simple_nu_filter.f90`: `NU_ORDER_LAMBDA = 1.0`,
  `NU_ORDER_CORE_RES_FAC = 1.5`, `NU_ORDER_CORE_KEEP = 0.95`, `NU_ORDER_CORE_MIN_VOX = 1000`,
  `NU_ORDER_FRAC_STEP = 0.01`.
- "Dense" is the density envelope's core (the binary before dilation that `automask3D` already returns),
  restricted to the NU domain and to observed voxels (`nu_observed_mask`). On gridding every voxel is
  observed, so the restriction changes nothing; on PCG it keeps out a dense voxel outside the solve support,
  which should not occur, and the log counts any it removes. The order prior used separate density weights
  (smoothed `|x|` at max(8 Å, 2 × FSC), Otsu, logistic); they go.
- The producer is a new submodule `simple_nu_filter_ordermask` of `simple_nu_filter`, in two levels: a
  routine on arrays (margin, label resolution, dense set, composed halves, FSC=0.5 resolution → `M`, `f`,
  statistics) that the tests call directly, and a routine on the live NU state that fetches the margin
  (`get_nu_evidence_map`) and the label resolution (`get_nu_local_resolution_map`) and calls the first.
  `nonuniform_filter_state` and `nu_filt3D` both call the second. It takes from `simple_pcg_order_prior`
  the quantile thresholds and widths (`evidence_thresholds`) and the logistic (`build_order_weights`).
- The tester is `simple_nu_filter_ordermask_tester`, sub-suite `order mask` of `unit_reconstruction`; it
  replaces the sub-suites `order prior` and `order prior offline`.
- solve3D's start stage is `AUTOMSK_START_STAGE` in `simple_solve3D_controller.f90` (read through
  `solve3D_automsk_stage()`), tied to `PCG_START_STAGE` as Hans's 9c17b8ca8 set it for `automsk=yes`; this
  plan first had it as its own constant.
- solve3D refuses `automsk=nu|order` at entry when NU filtering starts after stage 7 (`nu_stage > 7`), since
  those modes need an NU `filt_mode` in every stage they reach.
- Ruling 7 is applied to the shared final reconstruction (`simple_final_rec`), so refine3D's final
  reconstruction follows the same rule. It already inherits `automsk` on both routes; only the branches
  that name `nu` extend to `order`.
- `nu_filt3D` takes `automsk=no|order`. Under `order` its commander sets `filt_mode=nonuniform` before
  `params%new`, computes the pair's FSC as gridding assembly does (`evaluate_halfmap_pair`, backend
  `gridding`, density mask, phase-randomized), caps the bank at FSC=0.143/1.5 as assembly does, builds `E`
  with its core, and writes the order mask, the masked `_nu_filt` pair and average, the evidence map and the
  label resolution. `nu_envmsk` is unchanged. Without `automsk` its outputs are byte-identical to the base.
- `reconstructor_pcg`'s `shell_diag` and `get_shell_diagonal` go with the order prior, their only consumer
  (the compiler probe of Phase 4 proves it); `get_realspace_diagonal` stays (atomicity prior).
- `solve3D_addon` replays `automsk` from the base run's manifest (it is already a replay key): no change
  beyond dropping `pcg_order` and `pcg_order_lambda` from its refusal list and from the manifest key list.
- The solve3D stage-boundary reconstruction keeps stripping `automsk` (`strip_refine3D_planning_keys`); the
  masks are built by stage 7's first assembly.
- Phase 7 runs three datasets × two arms (`density`, `order`) × three restarts on PCG, plus one gridding
  restart per arm and dataset: 24 solve3D runs.

## 4. Not in this plan

- `automsk=nu`'s threshold and the [cross-validated threshold note](nu_evidence_mask_cross_validated_threshold.md)
  (only its status line changes).
- The PCG solve-support contract, the NU competition, the NU-evidence envelope.
- The atomicity prior beyond deleting its exclusion with `pcg_order`.
- `solve3D_cavgs`, `solve3D_chunks`, `refine3D_states_cavgs` and `reconstruct3D_wave`: they keep
  `automsk=no`.
- FLEX's PCG (`simple_flex_pca_pcg`).
- The NICE layer: it names no `automsk` value (checked at the base).
- Any default other than the spelling `yes` → `density`.

## 5. Change types and gates

| Phase | Kind | Gate |
| --- | --- | --- |
| 0 | baseline | builds, fast gate, the identity cases twice, the data present |
| 1 | identical | the Phase 0 identity cases byte-identical at one thread (`yes` spelled `density`) |
| 2 | identical for `no`, `density`, `nu`; result-changing under `order` | off-state identity; `order` against `density` on one assembly differs only in the declared files; sub-suite `order mask` |
| 3 | identical without `automsk`; offline report under `order` | `nu_filt3D` byte identity without `automsk`; the six pairs reported, not gated |
| 4 | identical (no run used `pcg_order` in the identity cases) | the Phase 0 identity cases; `pcg_recon` and `pcg_frac_update` pass; `rg` finds no retired symbol |
| 5 | identical without `automsk`; result-changing with it | `solve3D_cavgs stage rules` sub-suite with the new assertions; one solve3D smoke run per backend |
| 6 | docs and tests only | `check_guide_imports.py`, `check_test_registry.py`, `check_cli_chain.py` |
| 7 | science runs | reported, not gated |
| 8 | close | the report |

## 6. Acceptance

Numbers and their origin:

- `K_core >= 0.95`: guessed (2026-11-08), kept by Hans for the order prior and carried by decision 9.
- 1000 core voxels: guessed (2026-11-08), carried.
- `λ = 1`, 1.5 × FSC=0.5: ruled (2026-11-21). The dense set as the quantile domain: ruled and amended the
  same day (decisions 4 and 14).
- Unit test, backend consistency (decision 14): a fixture and its constrained twin, in which every voxel
  outside a support around the dilated envelope carries margin 0 and zero halves (as a density-constrained
  PCG pair does), give the same `f` and byte-identical `M` on the dense set. Origin: construction (the
  dense set is the same in both, and nothing outside it enters). Negative control, run once in Phase 2 and
  recorded: the same comparison with the quantile taken over the whole sphere must differ, which shows the
  test can fail.
- Unit test, exactness: the reported `K_core` equals the ratio measured on the masked output within 1e-4
  (relative); measured, carried from the order prior's core-criterion test.
- Unit test, the membrane fixture: the mean of `M` over the micelle is at most 0.6. Guessed from the
  fixture's construction (the micelle carries no evidence, so `M` tends to 1/2); Phase 0 dry-runs the
  fixture in numpy and freezes the bound before Phase 2 runs the test. If the dry run puts the mean above
  0.6, the fixture changes in Phase 0, not the bound in Phase 2.
- Unit test, threads: the same `f` at 1 thread and at 11, 16 and 24 threads, each repeated five times, and
  the same mask bytes for the same `f` (the thread rule for new OpenMP code).

Identity cases, frozen in Phase 0 (byte comparisons with `cmp` at one thread and a fixed seed, of the
published maps, half maps, `_nu_filt` products, masks, FSC files and the project's `ptcl3D` and `out`
segments; excluded: logs, `BENCH/` files, and lines holding absolute paths, run identifiers or times):

- `refine3D maxits=1 nthr=1 vol1=<truth volume> filt_mode=nonuniform ml_reg=yes` with
  `rec_backend=gridding|pcg` × `automsk=no|yes|nu`, six cases, on the simulated 1jxy project that the
  `cont_refine3D_1jxy` entry builds (Phase 0 runs the entry once and keeps the project and its truth
  volume). A case enters the gate only if its two base runs are byte-identical; the others are recorded and
  left out. (The registry entry itself declares no identity outputs and runs at eight threads, so
  `test=identity` on it would compare nothing.)
- Phase 4 adds the PCG operator's own gates: `pcg_recon` and `pcg_frac_update` pass as at the base (they are
  not registry entries, so they are gated on their pass, not on bytes).

The order-against-density comparison (Phase 2): the two PCG cases and the two gridding cases above with
`automsk=density` and `automsk=order`, on the Phase 2 build. They may differ only in the `_nu_filt` files,
`order_envmask3D_stateNN.mrc`, the refine3D provenance record and the digests it holds, and the logs.

Reported, not gated (Hans reads them):

- Phases 0 and 2, the two backends on the same particles: the gridding and PCG cases of the identity set
  start from the same `vol1`, so their first iteration should assign the same orientations (Phase 0 checks
  the `ptcl3D` segments and says so if they differ). On the base (`automsk=yes`): the Dice overlap of the
  two density envelopes and of their cores, and the two FSC=0.5 resolutions. On the Phase 2 build
  (`automsk=order`): the Dice overlap of the two dense sets, the two chosen `f`, the mean of
  `|M_gridding - M_PCG|` over the common dense set, and the two FSC=0.5 resolutions, beside the base
  numbers. This is the measurement of decision 14.
- Phase 3: per pair, `f`, `K_core`, the share removed around the core, and the maps
  `<pair>_order_nu_filt.mrc` beside `<pair>_nu_filt.mrc`, under `~/Processing/automsk_order/offline/<pair>/`.
- Phase 7: per dataset and arm, the final FSC=0.143 and FSC=0.5 of the final map, the stage 7-8
  references, every `NU ORDER MASK` line, and the run time, under
  `~/Processing/automsk_order/<target>/<arm>_<backend>_<restart>/`. Three restarts per arm are a lead, not a
  finding; the report says so.

## 7. Common rules and stop conditions

- No new key beyond those ruled (`automsk=order`, the `density` spelling, `automsk` on `nu_filt3D` and on
  solve3D's UI). No floor, tolerance or assertion relaxed unless this plan says so. Fix as you go, with a
  file-table row added before the edit. New OpenMP code is tested at 11 to 24 threads, repeated. The
  `THROW_HARD`/`THROW_WARN` arguments are never `//&`-continued (an `errmsg` variable instead). No
  short-circuit on allocatable arrays in one `if`. No `goto`.
- Policies change in the phase that changes the behaviour they describe: `automasking_policy.md` (the
  living description of the mode) in Phases 1 and 2, `nonuniform_filtering_policy.md` in Phases 2 and 4,
  `reconstruct3D_pcg_policy.md` in Phase 4, `solve3D_policy.md` in Phases 4 and 5.

Stop when:

- an off-state identity case differs and the difference is not a defect of the base path (a base defect is
  fixed in the phase and recorded);
- the order-against-density comparison differs in an undeclared file;
- a floor of section 6 would have to move;
- a file outside the table needs an edit that changes behaviour;
- the code contradicts a fixed decision (for example, a caller of `nonuniform_filter_state` that has no
  FSC=0.5 resolution, or a path on which the envelope's core is not available);
- an available validation exposes a regression against a non-data gate. Missing real-data inputs are recorded
  as outstanding and do not block implementation (user ruling, 2026-11-21).

Gates no ruling waives: the off-state identity cases at one thread; the sub-suites `order mask` and
`nonuniform filtering`; `scripts/run_fast_gate.sh`; `check_test_registry.py`; `check_cli_chain.py`.

## 8. Phases

**Phase 0, baseline (base `f9fd7ca21`).** Build Release and Debug; run `scripts/run_fast_gate.sh build`,
`simple_test_exec test=unit_reconstruction`, `check_test_registry.py` and `check_cli_chain.py` (both pass at
the base, checked 2026-11-21). Run `cont_refine3D_1jxy` once and keep its project; build the six refine3D
identity cases on it, each run twice; record which are deterministic. Run `pcg_recon` and `pcg_frac_update`
and record their outcome. Record the two backends' agreement on the base (section 6: envelope and core
Dice, FSC=0.5 resolutions, whether the orientations agree). Verify every command line of this plan against
the UI's accepted keys. Dry-run the membrane and soluble unit fixtures in numpy (`λ = 1`, the core at 1.5 ×
the fixture's FSC=0.5, the quantile over the dense set) and freeze the micelle bound. Real-data inputs,
outstanding when absent under the 2026-11-21 user ruling: the six `_unfil`
pairs in `~/Processing/solvent_prior/<pair>` (streptavidin, PfCRT, PorF, embb, exp_gate, tat; `mskdiam` 80,
160, 170, 150, 155, 150 Å as listed in the order prior's section 9; Phase 0 confirms each pair's value from
its folder); the solve3D inputs of PfCRT and PorF (the frozen copies of the regression study,
`~/Processing/solve3D_regression/{PfCRT,PorF}/`), a streptavidin project, and for each the solve3D command
line Hans gives, recorded in `cmd.txt`.
Exit: both builds and the gate pass (base failures recorded, not judged); the identity table lists each case
as deterministic or excluded; the bound is frozen in the plan's Progress. Missing real-data items and their
dependent runs remain outstanding without blocking implementation (user ruling, 2026-11-21).

**Phase 1, the rename.** `yes` → `density` in the parameter validation, the derived `envfsc` rule, the
`exec_refine3D_workflow` and `estimate_diam` defaults, the matcher, the UI choices and descriptors, the tests and the
comments; widen `automsk`; `automasking_policy.md` and the policies that spell `automsk=yes`.
Exit: the identity cases byte-identical with `automsk=density` in place of `yes`; `automsk=yes` is rejected
by parameter validation; `rg "automsk[=: ]*'?yes"` over `src/`, `doc/policies/` and `.github/skills/`
finds nothing.

**Phase 2, the order mask.** The producer and its tester (submodule pair); `ORDER_ENVMASK_FBODY`; the new
FSC=0.5 argument of `nonuniform_filter_state` and both callers; the order branch in the state filter
(request the envelope's core in this mode too, remove `E`'s skirt median, multiply by `A`, write or delete
the artifact, the log line); the matcher fallback (decision 11); `simple_final_rec` and
`prepare_bootstrap_rec_cline` (`order` where `nu` is named); the refine3D provenance list; parameter
validation (`order` needs an NU `filt_mode`, implies `envfsc=yes`); the UI choices; `automasking_policy.md`
and `nonuniform_filtering_policy.md`.
Exit: sub-suite `order mask` passes (fraction of the logistic; membrane fixture: `K_core >= 0.95`, the
micelle bound, `f` the smallest qualifying; soluble fixture: the rule keeps the core; exactness; backend
consistency against the constrained twin, with its negative control recorded; the off path equals
`automsk=density` byte for byte; threads); `nonuniform filtering` passes; the off-state identity cases pass;
the order-against-density comparison differs only in the declared files; the cross-backend numbers of
section 6 recorded beside the base's; one distributed PCG case (`nparts=2`) and its shared-memory twin choose
the same `f` (the largest difference of their order artifacts is recorded, not gated).

**Phase 3, the offline check.** `automsk` on `nu_filt3D` (UI, commander) as in the defaults; run it on the
six pairs. Exit: `nu_filt3D` without `automsk` byte-identical to the base on one pair; the six reports and
maps written to `~/Processing/automsk_order/offline/`, with the per-pair table in Progress.

**Phase 4, the order prior removed.** Every row of the file table tagged 4: the PCG strategy's re-solve
branch, order pair, weight, statistics and their uses; the operator's order terms, `solve_regularized`'s
order arguments, `shell_diag`; the modules, program, commander and exec case; the parameters, parser,
validation and UI inputs; solve3D's forwarding, utilities, manifest lists and add-on refusal; the
refine3D_states_cavgs requirement and the wave refusal; `write_support_provenance`'s argument; the
`_nu_evidence` write gate and the apply pair of the state filter; the file-name constants; the sub-suite
registration. `pcg_order_prior.md` moves to `rejected/` with the banner "negative result: in the loop it
did not do its job (Hans, 2026-11-21); its order weight and core guard carry into `automsk=order`", a
decision-log entry, `reconstruct3D_pcg_policy.md`, `nonuniform_filtering_policy.md`, `solve3D_policy.md`
and the other policies that name the keys.
Exit: `rg "pcg_order|order_prior|l_pcg_order|set_order_prior|get_shell_diagonal|PCG_ORDER_"` over `src/`
finds nothing; the build passes with no new warning; the identity cases byte-identical; `pcg_recon` and
`pcg_frac_update` pass; `check_test_registry.py` and `check_cli_chain.py` pass.

**Phase 5, solve3D.** The commander keeps the user's `automsk` (default `no`), refuses `nu|order` with
`nu_stage > 7`; the controller's mask policy sets `automsk` and `envfsc=yes` from `AUTOMSK_START_STAGE`;
the final reconstruction inherits the mode (decision 7); the solve3D UI offers `automsk`
(`no|density|nu|order`); `solve3D_policy.md` sections 2, 8 and 9. The reseed tester's assertion "every stage
keeps spherical support" becomes: stages up to 6 carry `automsk=no envfsc=no`, stages 7 and 8 the given
mode with `envfsc=yes`, and without the key every stage carries `no`.
Exit: sub-suite `solve3D_cavgs stage rules` passes with the new assertions; one solve3D run per backend
with `automsk=order` on the `cont_refine3D_1jxy` project reaches its final reconstruction with an order
artifact and `vol_msk` registered from it; a solve3D run without `automsk` gives stage command lines
identical to the base's.

**Phase 6, docs and tests.** The five skill files that name the order prior, `code_base_map.md` and the
Fortran indexes, the implementation-notes README rows (order prior, NU masking), the status line of the
cross-validated threshold note, the test policy's `unit_reconstruction` row, `simple_test_ui_class`'s suite
list. Exit: `check_guide_imports.py`, `check_test_registry.py`, `check_cli_chain.py` and the fast gate pass.

**Phase 7, the loop test.** On each dataset: one common `solve3D_cavgs` start (kept), then solve3D with
`cavg_ini_ext=yes` on copies of it, the dataset's command line plus `automsk=density` or `automsk=order`;
`rec_backend=pcg` three restarts per arm, `rec_backend=gridding` one per arm. Exit: the 24 run folders, each
with `cmd.txt` and its log; the per-dataset table in Progress; the maps left for Hans.

**Phase 8, close.** `automsk_order_report.md` beside the plan (section 8 of
`doc/how2s/how2write_agent_plans.md`), the plan to `completed/`, the commit message.

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/params/simple_parameters.f90` | widen `automsk`, its doc string; delete `pcg_order*`, `l_pcg_order*` | 1, 2, 4 |
| `src/main/params/simple_parameters_parse.f90` | delete the `pcg_order*` registrations | 4 |
| `src/main/params/simple_parameters_phases.f90` | values, `envfsc` rule, `order` validation; delete the `pcg_order` block and the `atm_prior` exclusion | 1, 2, 4 |
| `src/main/ui/simple_ui_params_common.f90` | `automsk` and `automsk_refine3D` choices and descriptors, the `envfsc` text | 1, 2 |
| `src/main/ui/simple/simple_ui_refine3D.f90` | delete `pcg_order*` (two programs) | 4 |
| `src/main/ui/simple/simple_ui_reconstruct3D.f90` | delete `pcg_order*` | 4 |
| `src/main/ui/simple/simple_ui_solve3D.f90` | delete `pcg_order*`; add `automsk` | 4, 5 |
| `src/main/ui/simple/simple_ui_heterogeneity.f90` | `automsk` descriptors | 1, 2 |
| `src/main/ui/simple/simple_ui_filter.f90` | `automsk` on `nu_filt3D`; delete `order_prior` | 3, 4 |
| `src/main/ui/simple/simple_ui_mask.f90`, `simple_ui_res.f90` | if needed (they use the common `automsk`) | 1 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | sub-suite list | 2, 4 |
| `src/main/commanders/simple/simple_commanders_refine3D.f90` | default `density`; bootstrap route `nu|order`; delete the `pcg_order` requirement | 1, 2, 4 |
| `src/main/commanders/simple/simple_commanders_imgproc.f90` | `estimate_diam` default `density` | 1 |
| `src/main/commanders/simple/simple_commanders_volops.f90`, `simple_commanders_mask.f90` | message; if needed | 1 |
| `src/main/commanders/simple/simple_commanders_resolest.f90` | `nu_filt3D` under `automsk=order` | 3 |
| `src/main/commanders/simple/simple_commanders_solve3D.f90` | keep the user's `automsk`, entry check; delete the `pcg_order` refusal keys | 4, 5 |
| `src/main/commanders/simple/simple_commanders_order_prior.f90` | delete | 4 |
| `src/main/exec/simple_exec_filter.f90` | delete the `order_prior` case | 4 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | sub-suite registration | 2, 4 |
| `src/main/commanders/test/simple_commanders_test_highlevel_cont.f90` | `yes` → `density` | 1 |
| `src/main/strategies/search/simple_matcher_refvol_utils.f90` | `density`; the order fallback | 1, 2 |
| `src/main/strategies/parallelization/simple_volassemble_strategy.f90` | pass FSC=0.5 | 2 |
| `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90` | pass FSC=0.5; delete the order prior | 2, 4 |
| `src/main/strategies/parallelization/simple_rec3D_wave_strategy.f90` | delete the `pcg_order` refusal | 4 |
| `src/main/volume/simple_nu_state_filter.f90` | FSC=0.5 argument, order branch; delete the apply pair and the evidence write | 2, 4 |
| `src/main/nu_filt/simple_nu_filter.f90` | interfaces, `NU_ORDER_*` constants, comment of `nu_filter_vols` | 2, 4 |
| `src/main/nu_filt/simple_nu_filter_ordermask.f90` | new: the producer | 2 |
| `src/main/nu_filt/simple_nu_filter_ordermask_tester.f90` | new: sub-suite `order mask` | 2 |
| `src/main/nu_filt/simple_nu_filter_envmask.f90`, `simple_nu_filter_vars.f90` | comments `yes` → `density` | 1 |
| `src/main/volume/simple_pcg_order_prior.f90`, `simple_pcg_order_prior_tester.f90` | delete (routines move to the producer) | 2, 4 |
| `src/main/volume/simple_order_prior_offline.f90`, `simple_order_prior_offline_tester.f90` | delete | 4 |
| `src/main/volume/simple_reconstructor_pcg.f90` | delete the order terms, `shell_diag`, `get_shell_diagonal` | 4 |
| `src/main/volume/simple_halfmap_diagnostics.f90` | drop `order_prior` from `write_support_provenance` | 4 |
| `src/main/volume/simple_refine3D_provenance.f90` | the order artifact | 2 |
| `src/main/simple_final_rec.f90` | `order` where `nu` is named; registered mask; comment | 2, 4 |
| `src/main/simple_convergence.f90` | if needed | 2 |
| `src/defs/simple_defs_fname.f90` | delete `PCG_ORDER_*` | 4 |
| `src/defs/simple_defs_automask.f90` | own the order artifact name without extending the core API | 4 |
| `src/main/solve/simple_solve3D_controller.f90` | `AUTOMSK_START_STAGE`, mask policy; delete the `pcg_order` forwarding | 4, 5 |
| `src/main/solve/simple_solve3D_utils.f90` | delete the `pcg_order` lines | 4 |
| `src/main/solve/simple_solve3D_manifest.f90` | drop `pcg_order*` from the key list | 4 |
| `src/main/solve/simple_solve3D_reseed_tester.f90` | `density`; delete the order-prior asserts; the mask asserts | 1, 4, 5 |
| `src/main/image/simple_image_msk.f90` | if needed | 2 |
| `doc/policies/3D/automasking_policy.md` | the living description of the mode | 1, 2 |
| `doc/policies/NU/nonuniform_filtering_policy.md` | the order step; the apply-pair paragraph goes | 2, 4 |
| `doc/policies/3D/reconstruct3D_pcg_policy.md`, `refine3D_policy.md`, `refine3D_auto_policy.md` | `density`; the order prior goes | 1, 4 |
| `doc/policies/3D/solve3D_policy.md` | the order prior goes; the mask forwarding | 4, 5 |
| `doc/policies/3D/wave_model_policy.md`, `solve3D_chunks_policy.md`, `solve3D_addon_policy.md`, `doc/policies/heterogeneity/refine3D_states_cavgs_policy.md` | the `pcg_order` refusals go | 4 |
| `doc/policies/test_environment_policy.md` | the `unit_reconstruction` row | 6 |
| `doc/implementation_notes/completed/pcg_order_prior.md` | to `rejected/` | 4 |
| `doc/implementation_notes/completed/pcg_decision_log.md` | the 2026-11-21 entry | 4 |
| `doc/implementation_notes/README.md`, `planned/nu_evidence_mask_cross_validated_threshold.md` | rows, status line | 6 |
| `.github/skills/simple-nonuniform-regularization/`, `simple-main-nu-filt/`, `simple-main-root/`, `simple-main-volume/` | the order prior goes; `automsk=order`; `density` | 1, 6 |
| `doc/code_overview/code_base_map.md`, `doc/code_overview/fortran-indexes/` | the modules | 6 |
| `doc/implementation_notes/planned/automsk_order.md` | Progress | all |

## Progress

(one entry per phase, written by the run)

### Phase 0 — blocked at data preconditions (2026-11-21)

The checkout is at the planned base `f9fd7ca21`. On this machine
`~/Processing/solvent_prior/` is absent, so none of the six specified frozen
`_unfil` pairs can be verified. The PfCRT and PorF regression directories exist,
but no `cmd.txt` or `.simple` project was found under `solve3D_regression/`.
Other `_unfil` pairs exist elsewhere, but their identity with the frozen inputs
has not been established. The user subsequently ruled that implementation and available validation proceed;
the six real-data pairs, three command lines, offline comparisons and 24 science
runs remain outstanding until their inputs are available.

The Release build and its complete UI validation passed at `f9fd7ca21`.
`pcg_recon` and `pcg_frac_update` passed after that build. The focused
`cont_refine3D_1jxy suite=nu_handoff` passed, and its seeded simulated project,
particles and truth map were captured for the identity cases. Earlier pre-build
executables were stale, so their results are not baseline evidence.
The Release fast gate passed 15/15 (113.3 s, warning-only budget exceeded under
concurrent builds). All six one-thread baseline cases (`gridding|pcg` ×
`no|yes|nu`) ran twice: their published MRC and FSC files and their `ptcl3D`
and `out` segments were identical after normalizing run paths. The two `yes`
backends assigned identical Euler angles and shifts. Their density envelope
Dice at mask > 0.5 was 0.9656; FSC=0.5 resolutions were 3.671 Å (gridding)
and 3.744 Å (PCG). Reconstructing the binary cores with the standalone masker
(zero dilation, retaining only mask values of one) gave core Dice 0.9253.
The isolated Debug baseline build and complete UI validation passed; the
Release baseline was installed after the fast gate.
A NumPy dry run with 4096 core margins uniformly spanning 0.9–1.1 and 4096
peripheral margins spanning 0–0.05, unit half-map products and the ruled
quantile/width formula selected `f=0.50`, kept 1.000 of core signal, and gave
mean peripheral `M=0.5563`; the 0.6 fixture bound is frozen.
For the backend-consistency negative control, adding 24,576 solvent voxels of
margin 0.37 versus 0 to the quantile domain selected thresholds 0.37 and
0.0480 respectively; mean absolute mask difference on the dense set was
0.0217. Using the dense domain alone gave the same result for both twins.

### Phase 1 — density spelling (2026-11-21)

The Release rebuild and complete UI validation passed. All six one-thread
identity cases matched the baseline MRCs, FSC files, `ptcl3D` and `out`
segments (run paths normalized); `yes` was spelled `density` in the two relevant
cases. `automsk=yes` was rejected by parameter validation. The static CLI chain,
test registry and `git diff --check` passed.

### Phase 2 — order mask (2026-11-21)

The Release build, UI validation, installation and 15/15 fast gate passed.
The original `no|density|nu`
cases matched the Phase 0 MRCs, FSC files and normalized project segments on
both backends. Density versus order matched outside the declared reference
artifacts; PCG diagnostic differences were timing fields. The expanded order
mask sub-suite passed 15 assertions, including the soluble, off-state and
thread-count cases; the nonuniform-filtering sub-suite passed. On the simulated pair, gridding selected `f=0.89`,
`K_core=0.9500`; PCG selected `f=0.90`, `K_core=0.9520`. The observed dense
sets had Dice 0.9334, and mean absolute order-factor difference on their
intersection was 0.0166 (the density core has envelope value one). Both FSC=0.5
resolutions remained 3.671 and 3.744 Å. A two-part distributed PCG run selected
the same `f=0.90` as shared memory; the largest mask difference was 8.11e-6.
A gridding run with `mskdiam=12` had fewer than 1000 core voxels and disabled
attenuation; its MRCs, FSC and project segments matched density mode exactly.
A pre-existing stale order artifact was also removed on this off path. A later
review verified the matcher fallback for compatible order, missing order, and
incompatible density artifacts; regional diagnostics now use the actual
coarsest bank cutoff and mark nonpositive denominators unavailable.

### Phase 3 — offline check (2026-11-21)

The Release build and UI validation passed. `nu_filt3D` without `automsk`
reproduced the baseline filtered halves, average and local-resolution map byte
for byte on the simulated pair. With `automsk=order` it produced the mask,
evidence, local resolution and filtered maps (`f=0.90`, `K_core=0.9558`,
removed around core 0.7423). The six named real-data offline runs remain
outstanding because their frozen pairs are unavailable here.

### Phase 4 — retire the PCG order prior (2026-11-21)

The Release build and complete UI validation passed without compiler warnings.
The six one-thread baseline cases matched the original MRCs, FSC files and
normalized project segments. `pcg_recon`, `pcg_frac_update` and the 33-assertion
atomicity-prior sub-suite passed. The retired source symbols are absent;
`pcg_order=yes` is rejected as an unknown key. The wider mode field also
rejects the invalid value `densityx`. The former prior's record was moved to
`rejected/` with its existing history preserved.

### Phase 5 — solve3D (2026-11-21)

The stage-rules sub-suite passed 1133 assertions, including all four modes in
stages 1–8; `automsk=order nu_stage=8` was rejected at entry. The final six
off-state identity comparisons matched the base. On the captured simulated
1JYX project, both gridding and PCG `solve3D automsk=order` completed through
the final reconstruction, produced `order_envmask3D_state01.mrc`, and registered
it as `vol_msk`. These are smoke tests, not the missing real-data science runs.
They preceded the final matcher-fallback safeguard and diagnostic reporting
fixes; the affected fallback cases, order-mask suite and backend comparisons
were rerun afterward.
After the fallback and diagnostic review fixes, all six off-state comparisons
and the order-versus-density comparisons passed again. The final fast gate
passed 15/15 in 52.57 s and the build was installed.

### Phase 6 — documentation and checks (2026-11-21)

The code map and Fortran indexes were regenerated. The guide-import, test
registry, CLI-chain, FLEX-DAG, descriptor and submodule-linkage checks passed.
The three new source files remain untracked and must be included when the
maintainer chooses what to commit. No Linux or BOX run was performed.

### Phase 7 — real-data validation outstanding

The six frozen half-map pairs, the three dataset command lines, offline
comparisons and 24 solve3D science runs remain outstanding. The simulated
smokes above establish execution and artifact handoff, not alignment benefit.

### Review fixes (2026-11-21, interactive session, not built)

Written in the desktop session on the working tree after the review; nothing was compiled
or run there. Static checks run after the edits: `check_test_registry.py`,
`check_cli_chain.py`, `check_guide_imports.py`, `check_submodule_linkage.py` and
`check_descr.py` pass.

- **Decision 15 (one defect).** `prep_final_rec_cline` took `envfsc` and the NU
  `filt_mode` from the caller's parameters and `automsk` from the refinement command
  line. A solve3D run whose last stage is before stage 7 (multi-state at its default of
  five stages) with `automsk` given therefore shipped an NU-filtered final map with a
  density-masked FSC although no stage had masked. Now `final_rec_mask_settings`
  derives all mask settings from the refinement command line, for the child command line
  and for `register_state_map`, which no longer takes `params`; solve3D logged when the
  requested mask could not act (after the rebase: 9c17b8ca8's entry warning instead). New assertions in `test_gold_standard_stages` (sub-suite
  `solve3D_cavgs stage rules`): stage 5 → `no`, no envelope FSC, no NU mask filter mode;
  stage 8 → `order` with both; `density` at stage 8 → envelope FSC, classical map.
- **Thread test.** `select_nu_order_mask` has no parallel region, so the 11/16/24-thread
  loop could not fail. Removed; the producer is documented as serial on purpose (the
  choice of `f` must not depend on a reduction order). The Phase 2 entry's "thread-count
  cases" refers to that loop.
- **Matcher fallback.** Under `order` the fallback refused a density or order artifact at
  another box (`state_mask_is_compatible`); the `density` and `nu` fallbacks rescale with
  `read_and_crop`. One rule now for all three: an artifact is used when it exists. Under
  `order`, either artifact missing falls back to the density artifact, and that one missing
  to the sphere.
- **Order artifact closure.** `wait_for_closure` after writing `order_envmask3D_stateNN.mrc`
  in the state filter (distributed workers read it through the fallback) and in
  `nu_filt3D`.
- **`nu_filt3D` output name.** The order mask is now `<outvol>_order_envmask.mrc`
  (`ORDER_ENVMASK_SUFFIX` in `simple_defs_automask`), beside the program's other outputs.
- **Text and style.** `automasking_policy.md`: the stale `no|yes|nu` validation line, the
  order artifact, the order fallback rule, decision 14 (one solvent region on both
  backends) and decision 15; `solve3D_policy.md`: decision 15. The solve3D UI help for
  `automsk` says it acts from stage 7. The background comment of `simple_nu_filter_vars`
  (it had named the NU-evidence background `automsk=density`); the misplaced interface
  comment of `get_nu_evidence_map` and comments for the three new interfaces; the
  `NU_ORDER_*` block aligned like `NU_ENVMASK_*`; no `;`-joined statements in the producer
  or tester; the blank-line pairs left in four UI blocks.

Recorded, not changed: the negative control for the backend-consistency test was run in
NumPy in Phase 0, not in Fortran in Phase 2. The standalone refusal of `automsk=nu|order`
is checked by program name, for `mask` and `automask` only. No Progress entry reports the
unobserved dense-voxel count of a PCG run. `solve3D_addon` from an `automsk=order` base
run has not been exercised. Under `order` the final reconstruction registers the soft order
mask (0.5 to 1 inside the envelope) as `vol_msk`; nothing in SIMPLE reads it back.

Still to run on the Mac build: `unit_reconstruction suite=order_mask` and
`suite=solve3d_cavgs_stage_rules`, the six off-state identity cases, the
order-against-density comparison, `nu_filt3D automsk=order` once, and the fast gate.

### Rebase onto `5481d42c5` (2026-11-21, interactive session, not built)

Master had gained Hans's 9c17b8ca8, which gave solve3D `automsk=yes|no` from stage 7 on its own:
`set_refine3D_mask_policy`, `solve3D_automsk_stage()` with `AUTOMSK_START_STAGE` tied to
`PCG_START_STAGE`, a warning and `automsk=no` (recorded in the manifest) for a run that ends before
stage 7, a binary UI input, tester assertions and policy text. The conflicts were resolved on that
structure: the mask policy and the entry check now take `no|density|nu|order` (the entry check keeps
this plan's `nu_stage` refusal), this plan's duplicate `solve3D_automsk_start_stage()`, its inline
policy and its late log line are gone, the UI input is the four-value one, the tester asserts
`density` where 9c17b8ca8 asserted `yes` and drops the `pcg_order` asserts, and `solve3D_policy.md`
merges both texts without `pcg_order`. The static checks above pass on the merged tree.

## Appendix: the two workflows at the base

| | `automsk=nu` | `pcg_order=yes` |
| --- | --- | --- |
| backends | gridding and PCG | PCG only |
| statistic | evidence margin at `amsklp` | the same margin, plus the label field |
| threshold | null median + `nu_msk_sig` (3) MADs, over the solvent majority of the sphere or the density envelope's dilation ring | `(1-f)` quantile, `f` chosen by the core guard |
| shape | binary: MRF, components, 1 Å growth, 6 Å edge | soft logistic, no morphology |
| acts on | the `_nu_filt` references (multiply) and the background clamp | the operator (a ridge `λ h (1-w) x²`), then a second cold solve; the references through the ML replay of that pair |
| micelle | zeroed | attenuated by up to `1/(1+λ)` where the operator is flat, less where it is not |
| timing | current evidence for the current references; lag one only for the gridding FSC mask and the matcher fallback | lag one throughout; off at the first iteration and at a new box |
| validity | null size, solvent majority | core size, core signal, a qualifying setting |
| prediction | none: the mask is what is applied | the operator lost 2-4.5 × the predicted core signal on membrane pairs |
| cost | one evidence pass and morphology | a second cold solve per half, band transforms |
| final map | NU artifact as mask, NU `filt_mode` | prior-free |

Both start from the same statistic. A fixed MAD multiple does not transfer between noise levels (in the
numpy check of the cross-validated threshold note the right multiple ran from about 7 to about 30); the
core guard adapts because the core is defined from the pair's own FSC. The order prior's strength cannot be
chosen from the half maps either: cross-half prediction prefers no removal at every noise level (same
note, section 4), which is why `λ` is ruled rather than scanned.
