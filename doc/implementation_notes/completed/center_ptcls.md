# center_ptcls: centring particles at scale

Status, 2026-11-14: delivered (Phases 1 to 3, 2026-11-13), uncommitted for Hans's review; closed with the 2D student.
The algorithm page is `doc/algorithms/center_ptcls.md`. Approved by Hans 2026-11-13, carried out interactively.
Base: `f2e907c24`. Supersedes: none.
Depends on: refine2D's probabilistic dense table with the polar shift map (validated: yes, in production in
solve2D) and solve2D's reference generation (validated: yes). It pauses `completed/class_student_2D.md`,
whose student needs centred particles (that plan's decision 9 and its probes after Phase 0).

## 1. Why

On the real sets the picking error is 8 to 10 Å at the median and 16 to 22 Å at the 90th percentile
(`class_student_2D.md`, Phase 0). Both students tolerate about 5 Å: the 3D student's cone holds its mode for
96% of particles at 6.4 Å of shift, and for 28% at 12.8 Å. Today only a full solve2D centres particles, at
about 30 iterations of reading every particle. The probes measured the alternatives on trpv1's 35,645
picks against the full run's converged shifts:
- classical centring on the rotational average: 19% within 3 Å;
- one search pass against 100 full-resolution class averages: 30% (median residual 4.8 Å).

Two converged runs agree for 69%.

## 2. The workflow (Hans, 2026-11-13)

1. Sample about 10,000 particles.
2. Generate a few references from the sample: refine2D downsampled to the library's crop at a low
   resolution, with very few classes.
3. Search every particle against those references in one pass: exhaustive in-plane rotation and a wide
   shift prior. Write the shifts.

## 3. Decisions

1. **Name and layering.** The program is `center_ptcls` (`center` exists and centres volumes; `center` is
   also a key). It runs whole programs (refine2D twice), so it is a controller: `commander_center_ptcls`
   with its `exec_center_ptcls(cline)` in `src/main/commanders/simple/`, as `exec_refine3D_workflow` is.
   Other applications call the commander on a child command line, as controllers do for whole programs.
2. **Defaults** (revised after Phase 1, Hans, 2026-11-13): a 10,000-particle sample (all particles when
   fewer), 10 classes, references at `lp` 20 Å, the pass at `lpstop` 12 Å, `shift_sigma` 10 Å (a 30 Å
   disk). bgal is too small a set to optimise on.
3. **The crop** comes from the library's rule at lp (`lpstages_setlims`): at 20 Å it gives
   `smpd_target = 6.7` Å, so the 88-pixel floor on bgal and trpv1. It is not hard-coded.
4. **References:** the classes' own averages, with `center=yes`. A rotational average is not used (it
   centred 19% within 3 Å).
5. **The pass:** refine2D's probabilistic dense table under the polar shift map. Every particle starts as
   not yet centred, so the map searches the 3-sigma disk coarse-to-fine and measures each particle's
   posterior shift width.
6. **Outputs:** only the shift (`x`, `y`), its width (`shsig`) and the move (`shincr`), in `ptcl2D` and
   `ptcl3D`. A later run then treats the particles as centred: its dynamic centring searches the window,
   not the disk. Every other particle parameter is left as it was. A particle whose drawn class carried no measured width
   (an unrefined candidate) is not moved and keeps the sentinel, so a later run searches its disk (3.6% of
   trpv1's picks, 1.25% of the test's). The pass runs on a work copy of the
   project, and its shifts are copied back by the project owner in one write.
7. **No new key.** The sample size, class count, the two bands and the shift prior are existing keys
   (`nsample`, `ncls`, `lp`, `lpstop`, `shift_sigma`) with this program's defaults.
8. **Its own CTest entry** (2026-11-13, raising `SIMPLE_CTEST_BUDGET` from 44 to 45, for Hans to confirm):
   `center_ptcls` is a program of its own, called by others, and its workflow test runs its own route on the
   multi-state fixture, independent of every other entry's programs, as `pose_student`'s does.

## 4. Acceptance (measured, frozen before Phase 2's runs)

On trpv1 (35,645 picks) and bgal (5,513), against each set's converged full solve2D:
- the residual's median and 90th percentile;
- the shares within 3 Å and within 5 Å;
- the wall time beside the full solve2D's.

Reported beside classical centring (19% within 3 Å) and the single pass against 100 classes (30%). The
expected outcome (guessed) is a median residual near 5 Å, the students' tolerance. Not gated: the downstream
check that solve2D started from the centred particles reaches the same class resolution, with fewer disk
searches.

## 5. Phases

- **Phase 1, prototype with existing programs.** `selection`, refine2D on the sample, one refine2D pass over
  all particles with the references. Measure section 4's numbers for lp 20 and 15 Å, and `shift_sigma` 10
  and 13 Å. Exit: the numbers in Progress and the defaults confirmed or changed.
- **Phase 2, the program.** `commander_center_ptcls` and `exec_center_ptcls`, the UI and exec chains, a
  workflow test in miniature, the outputs of decision 6. Exit: the fast gate, the unit and UI suites, the
  workflow test, Phase 1's numbers reproduced by the program on trpv1 and bgal.
- **Phase 3, close.** An algorithm page (`doc/algorithms/center_ptcls.md`), the test policy's row, this note
  completed; the 2D student plan resumed.

### File table

| File | What | Phases |
| --- | --- | --- |
| `src/main/commanders/simple/simple_commanders_center_ptcls.f90` (new) | the commander and its controller | 2 |
| `src/main/ui/simple/simple_ui_*.f90` (the 2D group) | the program's registration | 2 |
| `src/main/exec/simple_exec_*.f90`, `src/main/apis/simple_exec_api.f90` | routing | 2 |
| `src/main/commanders/test/` (a workflow test), `simple_test_registry.f90`, `production/CMakeLists.txt`, `src/main/ui/simple_test/*.f90` | the test | 2 |
| `doc/algorithms/center_ptcls.md` (new), `doc/policies/test_environment_policy.md`, `doc/implementation_notes/README.md` | docs | 3 |
| `doc/implementation_notes/completed/center_ptcls.md`, `completed/class_student_2D.md` | this note; the resumed plan | 1-3 |

## Progress

### Phase 1, the prototype with existing programs (2026-11-13, the Dell, base `f2e907c24`)

**Commands** (`/data2/agent_runs/class_student/runs/cp_trpv1`, `cp_bgal`, Release build of the base,
24 threads):
- `selection nran=10000` (all of bgal's 5,513);
- `refine2D ncls=<n> cls_init=rand refine=prob objfun=euclid lp=20 box_crop=88 center=yes shift_sigma=10`
  on the sample;
- one `refine2D refs=<the sample's last class averages> refine=prob maxits=1 lp=<pass> box_crop=88` over
  all particles.

The references for each set are a converged full solve2D (trpv1: 100 classes, 468 s; bgal: 50 classes,
58 s).

| set, method | median residual | within 3 Å | within 5 Å | within 6.4 Å | time |
| --- | --- | --- | --- | --- | --- |
| trpv1, raw picks | 8.0 Å | 0.09 | 0.24 | 0.36 | |
| trpv1, classical centring | 5.7 Å | 0.18 | 0.42 | 0.57 | |
| trpv1, 5 references at lp 20, pass at lp 20 | 4.3 Å | 0.33 | 0.57 | 0.67 | 23 + 9 s |
| trpv1, 5 references at lp 20, pass at lp 15 | 3.9 Å | 0.37 | | | 23 + 9 s |
| trpv1, 10 references at lp 20, pass at lp 12 | 3.3 Å | 0.45 | 0.67 | 0.76 | 30 + 8 s |
| trpv1, 20 references at lp 20, pass at lp 15 | 3.3 Å | 0.45 | | | 34 + 8 s |
| bgal, raw picks | 10.1 Å | 0.07 | 0.16 | 0.25 | |
| bgal, 5 references, pass at lp 15 | 6.7 Å | 0.13 | | | 10 + 2 s |
| bgal, 10 references, pass at lp 12 | 5.6 Å | 0.25 | 0.45 | 0.55 | 10 + 2 s |
| bgal, 20 references, pass at lp 15 | 6.3 Å | 0.31 | 0.44 | 0.50 | 9 + 2 s |
| bgal, ceiling: one pass against the full run's 50 classes | 2.1 Å | 0.57 | 0.66 | 0.70 | |

Further results:
- **Not worth changing:** the wider shift prior (13 Å) changed nothing. References at lp 15 were worse than
  at 20. Forcing 15 reference iterations instead of refine2D's own stop (7 to 13) changed nothing.
- **The floor:** two converged solve2D runs agree within 3 Å for 69% (trpv1) and 76% (bgal) of the same
  particles.
- **bgal's tail:** its 90th percentile (18 to 22 Å) does not move under any method, the ceiling included. The
  worse-matching half of its particles carries that tail (its 90th percentile is 25 Å, the better half's
  13 Å). These are particles whose centre no run agrees on.

**What it means, in short.**
- On trpv1 the workflow centres two thirds of the picks within 5 Å and three quarters within 6.4 Å (raw:
  24% and 36%), in 38 s against solve2D's 468 s.
- On bgal it roughly triples the share within 3 to 5 Å, but stays well below the ceiling of converged
  classes. Its low-resolution references, not the pass, are the limit there.
- Ten references help on both sets. The pass gains from a finer band than the references (12 against 20 Å).
- The residual against one converged run partly measures that run's class centres. The test that decides
  is downstream: the students' reads on centred particles, and solve2D started from them.

### Phase 2, the program (2026-11-13, the Dell, on `f2e907c24`)

**Code.**
- `simple_commanders_center_ptcls.f90` (new): `commander_center_ptcls` with its controller
  `exec_center_ptcls`.
- `center_ptcls` registered in `simple_ui_refine2D.f90` (`nparts` optional) and routed in
  `simple_exec_refine2D.f90`.
- The workflow test `center_ptcls`: `simple_commanders_test_highlevel_center.f90` (new), the parent module's
  commander, the router, the test UI, the registry, CTest. The CTest budget was raised to 45 (decision 8).

**Defect found and fixed.** A pass numbered after the references' iterations searched only the window,
because a run marks particles not yet centred only at its first iteration. Numbered so, the pass reached 38%
within 3 Å on trpv1, against 45% for the same references as a first iteration. The pass now runs as the first
iteration of its own run in `center_pass/`.

**Particles left not centred.** A particle whose drawn class carried no measured width is not moved and
keeps the sentinel (decision 6): 3.6% of trpv1's picks, 1.25% of the test's.

**Measurements** (Release copy of the tree, 24 threads, `SIMPLE_SEED=1`, the program's defaults):

| set | median residual | within 3 Å | particles not centred | time |
| --- | --- | --- | --- | --- |
| trpv1 | 3.5 Å | 0.43 | 1,288 (3.6%) | 40 s |
| bgal | 7.7 Å | 0.22 | | 11 s |

These match the prototype's 3.3 to 3.5 Å and 0.43 to 0.45 on trpv1, and its 0.25 on bgal, within the
variation between random samples.

**The workflow test**, 2,000 beta-galactosidase particles with planted shifts within 6 px: all four checks
pass. Reported: planted median 10.4 Å, residual median 6.4 Å, 18% within 3 Å, 36% within 5 Å, 1.25% not
centred, 5 s with the fixture cached (Release); 10 s on the Debug build, which checks array bounds, also
passing.

**Gates.**
- the fast gate, 15 of 15;
- `unit_ui`, 1046 of 1046;
- `check_test_registry.py` and `check_cli_chain.py` consistent.

### Phase 3, close (2026-11-13)

- `doc/algorithms/center_ptcls.md` (new, entry 18 of the algorithm index), and the test policy's row and
  entry.
- This note stays under `planned/` until Hans has reviewed the diff and confirmed decision 8. The 2D student
  plan resumes from its Phase 1 after that.

