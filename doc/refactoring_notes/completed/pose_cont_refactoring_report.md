# Continuous Cartesian pose refactoring (`refine=cont`): review report

Report of the delegated run that carried out the plan
[pose_cont_refactoring.md](pose_cont_refactoring.md), Phases 0-10, on the Dell on
2026-10-01. Every phase ended DONE; Phase 3 stopped once, on the 4 A E25 floor, and resumed
under ruling R1. The plan holds the contract and one section 15 row per phase (files,
decisions, evidence, log paths).

The run worked against `714a15a99` and committed nothing. Its change was brought into the
main checkout as uncommitted working-tree changes on top of `adaaea4c7` (three textual
conflicts with the later commits resolved by hand). The per-phase diffs (`phase_N.diff`),
the cumulative ones (`full_after_N.diff`) and the complete diff against `714a15a99`
(`full.diff`) are in `~/agent_runs/runs/pose_cont/result/review/` on the Mac and in
`~/pose_cont_autorun/review/` on the Dell. `full.diff` also records a deletion of
`CMakeFiles/cmake.check_cache`, an artifact of the copy to the Dell and not part of the
change. Paths of the form `scratch/...` below are under `~/pose_cont_autorun/` on the Dell;
the bulky run outputs there were deleted during the run (ruling R6), the logs are kept.
Changes made in the main checkout after the transfer, at review, are in "Changes at review"
below; the phase table and the evidence describe the run itself.

## What changed, phase by phase

| Phase | Diff | Change |
| --- | --- | --- |
| 0 | `phase_0.diff` | Baseline only: E25 print lines; 3 runs each of the polar guards; beta-gal runs of both LM routes, the polish and `refine3D_auto`; inventories. E25 re-measured at 1 000 particles (R2) |
| 1 | `phase_1.diff` | Program `refine3D_pose_cont`, its workflow, UI, policy and the four `refine3D_auto` hooks removed; `exec_refine3D_auto` identical to `e86286f2^` (C13) |
| 2 | `phase_2.diff` | New `src/main/cftc`: `cartft_calc` (references, particle slots, objective/gradient/normal terms) and `cartft_pose_opt` (LM transaction); old refiner and adapter become forwarders (C3, C17). Interfaces of checkpoint 1 in `review/phase_2_interfaces/` |
| 3 | `phase_3.diff` | Objectives follow `objfun`: `cc` uniform without sigma2, `euclid` normalized polar loss with `exp(-L)` (C4-C6); O4 preparation, O5 (a) taper, O6 bounds `trs`/`athres` (the rotation bound became `athres_cont`, default 10 degrees, by ruling of 2026-10-02); forced `cc_emit_sigma` gone |
| 4 | `phase_4.diff` | `build%cftc`; prepared reference volumes as `cart_refvols_{even,odd}.bin` (6.6, C22); `l_cart_refine`; inverse-envelope constructor deleted |
| 5 | `phase_5.diff` | `strategy3D` representation-neutral; abstract `strategy3D_pftc` owns `s`; the nine polar strategies extend it |
| 6 | `phase_6.diff` | Batch preparation of Cartesian slots (`prep_cart_batch`, `prepimg4align_cart`); generic `euclid_sigma2` `new`/`calc_sigma2`; no polar state under `l_cart_refine` |
| 7 | `phase_7.diff` | `strategy3D_cont`, the matcher's only continuous-pose symbol (C1); adapter, old strategy, refiner, stats module deleted (C2); record slot 51 `corr_cart`, 52 the improved flag (O8, C14, C15) |
| 8 | `phase_8.diff` | `refine=cont`; `refine3D_auto pose_cont=no|yes|only` (O7 b) with entry seed check and stage bypass (C11, C12); one LM route (C16); no early stop under `refine=cont` (O2 interim); CTest `cont_refine3D_1jxy`, budget 30 (C18) |
| 9 | `phase_9.diff` | Polish as a scheduled `refine=cont` pass after every discrete pass over the same sample, shared memory and distributed, one reconstruction per iteration (C8, C19) |
| 10 | `phase_10.diff` (written by the driver) | Policies (`refine3D_policy.md`, `refine3D_auto_policy.md`, `test_environment_policy.md`), superseded implementation notes marked, plan moved to `completed/` |

## Decisions taken under delegation (where recorded)

- O4, O5 (a), O6, O7 (b), O8 and the 6.6 layout adopted as recommended: Phase 3, 4, 7, 8 rows;
  sections 5, 6.2, 6.6, 6.5 "As built". Checkpoints 1-7 recorded, not waited for (section 15).
- C16 route: `shift_then_joint` kept (the polar shift-first staging, safe in the capture
  experiment, E25 floors measured on it); `joint` was 30-40% faster on beta-gal at equal
  objective. Phase 8 row, checkpoint 6. The optimizer's step/iteration optionals stay as
  tester seams; production uses the defaults. Revised at review (2026-10-02): production runs
  `joint`; `shift_then_joint` stays for the single-pass capture tests (`cont_route`, internal).
- Shell membership `nint(r)` inside the Nyquist circle; a zero-power observation is an
  invalid preparation (Phase 3).
- `refine=cont`: `inpl_cont` overridden to `no` rather than rejected; `oritype=ptcl3D` checked
  on `refine3D` only (its children inherit `refine`); initialized poses checked by
  `check_cont_seeds` in both iteration strategies and in `refine3D_auto`, before the
  random-orientation fallback (Phase 8, 6.5).
- `refine3D_auto pose_cont=only`: final reconstruction unchanged (its `bootstrap_rec3D` sigma
  pass, when needed, stays polar); `pose_cont=yes` polishes the main loop only (7.1, Phase 8/9).
- Polish pass named `refine=cont pose_cont=yes` (no new key for the workers), refused on the
  command line; sigma2 prepared again before it (Phase 9, checkpoint 7). Withdrawn at review
  (2026-10-02): the polish writes into the discrete pass's sigma2 transaction, starting from its
  residual rows (plan section 15, Review row; N30).
- A shared-memory Cartesian pass writes its sample to the project before the assembly
  (trailing reconstruction reads it from the file) (Phase 8).
- Implementation notes describing the old design marked superseded, not rewritten (Phase 10).

## Deviations from the plan, and why

- E25 at 1 000 particles (R2); its 4 A case starts from the 8 A poses (R1, after the Phase 3
  stop; finding: from 15 deg directly at 4 A, `euclid` leaves 20% of the particles in local
  minima, `cc` 0.4%).
- Polar guards once, library and high-level tests on Release only (R3); FSC no criterion at
  2 000 particles or fewer (R5 as amended; Phases 8-9 applied the first wording).
- N20 stage h: the polar `refine3D_auto` with its defaults scrambles the 15-degree seeds, and
  1JYX is D2 in a c1 gate; stage h uses `ref_pose_init=cc` and a symmetry-aware error. Recorded
  in 10.3 before each rerun; floors unchanged.
- N20 stage g compared with stage h's polar run rather than stage f (f is a `refine3D`
  against the truth map, not like-for-like) (10.3).
- `doc/code_overview/fortran-indexes` not regenerated (your instruction); they still list the
  removed modules.
- Files added to 8.3 when a phase needed them (each row says "added in Phase N").

## Evidence against the Phase 0 baseline

- E25, `lib_cart_align3D` (1 000 particles, ground truth): Phase 0 cc 8 A median rotation
  0.194 deg, basin 0.9870; now cc 0.188 deg / 0.9910, euclid 0.199 / 0.9710, 4 A 0.127 / 0.126
  deg; 47/47 checks. Not bit-reproducible across threaded runs (the fixture inputs differ at
  1e-8; serial runs are identical; Phase 7 finding).
- `cont_refine3D_1jxy` (new): `refine=cont` euclid/cc 0.91/0.23 deg from 15-degree seeds;
  distributed equal to shared memory particle by particle; `refine3D_auto pose_cont=only`
  0.94 deg; the polish improves `refine3D` 0.90 -> 0.63 deg and `refine3D_auto` 1.52 -> 0.22 deg
  (symmetry-aware) (`scratch/logs/cont_refine3D_1jxy_p9.log`). N29 on the nanoparticle within
  its margin.
- Polar guards (Phases 5-7, once each): all PASS, at or near the Phase 0 spread, with one
  confirmation run where outside (Phase 5-7 rows).
- Beta-gal (5 513 particles; half-map FSC usable under R5):

| Run | Wall | Peak RSS | FSC 0.5/0.143 (A) |
| --- | --- | --- | --- |
| Phase 0 polar `refine3D_auto` | 18:25 | 10.74 GB | 4.08 / 3.55 |
| `refine3D_auto pose_cont=only` (Phase 8) | 10:47 | 2.58 GB | 4.03 / 3.51 |
| `refine3D_auto pose_cont=yes` (Phase 9) | 27:44 | 10.84 GB | 4.08 / 3.51 |
| Phase 0 single `refine=pose_cont` pass | 3:34 | 2.53 GB | 4.35 / 3.75 |
| single `refine=pose_cont` pass, Phase 7 code | - | - | 4.35 / 3.80 |

- Fast gate 13/13, 2.7 s; `unit_cart_align3D` 634 checks, 0.78 s, also in reverse order.

## Changes at review (2026-10-02, main checkout, after the transfer)

Made on top of the transferred tree, uncommitted; plan section 15 has one row each (Review,
Ruling O6, Ruling C16), sections 4 and 10.3 the rulings and the new tests.

- Transfer: 57 files byte-identical to the Dell result, 13 merged three ways with the six
  commits after `714a15a99`. Three conflicts resolved by hand: the `refine` doc comment in
  `simple_parameters.f90` (`pose_cont` replaced by `cont`), a comment the comment refactor
  reworded inside the block the run re-indented in `simple_strategy3D_matcher.f90`, and the
  N29 block before the new PASS line in `simple_commanders_test_single.f90`. The plan took
  the run's version (a superset of the plan commits); the deletion of
  `CMakeFiles/cmake.check_cache` in `full.diff` was not applied.
- Polish sigma2 fallback (review finding P1). The polish opened a second sigma2 transaction,
  which discarded the discrete pass's fresh residuals, so a particle with an invalid Cartesian
  slot kept the previous generation's sigma2 under its new discrete pose. The polish now
  writes into the discrete pass's transaction: `euclid_sigma2%read_pending_range` (called by
  `prep_sigmas_objfun` when `l_cont_polish`) starts from the discrete range of its part, valid
  Cartesian contributions overlay it, and `write_sigma2` replaces the range through a staged
  file and an atomic rename (range files are otherwise created once). Shared memory and
  distributed (`simple_refine3D_strategy.f90`); N30.
- Polar/Cartesian sigma2 interchangeability (finding P2). N31 compares the two for a
  signal-bearing particle (nonzero reference, astigmatic CTF, phase flip, pose off the polar
  grid). A ring at radius k and the pixels with `nint(r) = k` sample a structured residual
  differently, so single shells differ by up to about 25 %; the gate is the pixel-count
  weighted band total within 4 % and every shell within [0.6, 1.5], taken from an independent
  model of the fixture. The first run failed (band total 0.72) on a fixture artifact: the
  preparation's phase flip moved background noise above the first CTF zero (k 11.3) into the
  mask; the defocus is now 0.15/0.18 um. N32 is the canonical round trip of a generation a
  Cartesian pass writes and a polar sigma owner reloads.
- Exact-zero policy (maintainer): the reducer's rule stands. A record with any non-positive
  or non-finite shell stays in the particle rows and is excluded from the groups; its particle
  reloads its half's group (N32). A numerically exact Cartesian match is zero only up to
  rounding (N13).
- Rotation bound (ruling on O6): `athres_cont`, default 10 degrees, the halfwidth of the total
  rotation from the seed, wide enough to capture wrongly assigned orientations and independent
  of `athres` (the polar neighbourhood radius) and `prob_athres` (which narrows to the mean
  angular motion). Parameter, validation in (0, 180], `strategy3D_cont`, UI of `refine3D` and
  `refine3D_auto`; E25 and N20 set 15, the Phase 0 bound; `test_rotation_bound` and N21.
- LM route (ruling on C16): production runs the joint stage alone (30-40 % faster; repeated
  passes find what one misses). `shift_then_joint` stays, through `cartft_pose_opt%new
  (shift_first=)` and the internal parameter `cont_route` (default `joint`, not in the UI);
  E25 and N20 stages a-c keep it and their floors, N20 e-h and the polish run joint. E16
  checks both routes.
- Policy text (`refine3D_policy.md`) updated for the polish's sigma2, `athres_cont` and the
  route.

Test status: N30 and N32 passed on the first build. N31 (revised fixture), the athres_cont and
route changes have not been run yet: `unit_cart_align3D`, `unit_ui`, `lib_cart_align3D` and
`cont_refine3D_1jxy`. Stages e-h of N20 are the first workflow runs of the joint route.

## Open for your review

Carried forward, with the maintainer's comments and the code facts, in
[pose_cont_refactor_aftermath_issues.md](../planned/pose_cont_refactor_aftermath_issues.md).


- O2: the convergence rule under `refine=cont` (interim: none, runs to `maxits`).
- Whether the optimizer's stage solvers, `refine_pose` and its optional arguments stay public.
- Polar non-probabilistic shared-memory modes with trailing reconstruction read the previous
  iteration's sample at assembly (pre-existing; fixed only for Cartesian passes).
- The polish costs about 2 min per beta-gal iteration on 2 parts; `joint` would be faster (C16).
  Done at review: production runs `joint` (2026-10-02); its beta-gal timing is not measured yet.
- Rerun the polar guards on the merged tree: two commits after `714a15a99` changed the
  simulated-workflow validation the run's guards used.
- The run's new code predates the comment rule of `simple-modern-fortran`; a comment pass over
  `cftc/` and `strategy3D_cont` is outstanding.
- Regenerate `doc/code_overview/fortran-indexes`; consider skill updates (`simple-refine3d`,
  `simple-main-strategies`) for `refine=cont` and the polish.
- Every adopted recommendation (O4-O8, 6.6) remains open to reversal.
