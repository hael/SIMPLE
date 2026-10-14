# Stream tests in the new test environment: handover to Ruben

Hans, 2026-09-24. Updated 14 October 2026 against the current source: the delivered items
are marked as such below, and the sections that still ask for work say what is left.

Your seven stream test programs (`simple_test_exec test=...` in
`src/main/commanders/test/simple_commanders_test_stream.f90`) have been moved into the
two-tier test environment (plan: `doc/refactoring_notes/completed/uniform_test_environment_refactoring.md`,
section 9.7, "stream"). Unlike the SINGLE tests, these already checked real things, and two of
them (optics assignment and the sieve) check exact truth. The move kept every check: each
`THROW_HARD` became an assertion with the same condition and message. This note says where each
test lives now, what changed, and what each one should check beyond counts and files.

## Where the tests are now

| was | now | runs |
|---|---|---|
| `sieve_cavgs` | `test_collect_and_reject_hard_gates` in `src/main/sieve/simple_ptcl_sieve_tester.f90`, sub-suite `particle sieve` of `unit_project` | fast gate, every build with tests |
| `assign_optics` | `src/main/stream/simple_stream_tester.f90`, sub-suite `optics assignment` of `lib_stream` | nightly (`ctest -L library`) |
| `gen_pickrefs` | same module, sub-suite `picking references` of `lib_stream` | nightly |
| `pick_extract` | same module, sub-suite `pick and extract` of `lib_stream` | nightly |
| `master` | `run_stream_heartbeat_tests` in `src/main/stream/master/simple_stream_master_tester.f90` (`test_fork_gui_status`), sub-suite `stream heartbeat` of `forked_process`; the heartbeat assembly from status records is `test_stream_heartbeat` in `src/utils/gui/simple_gui_assembler_tester.f90` | platform entry (`ctest -L platform`); the assembler test runs in `unit_ui` |
| `preproc` | unchanged commander, CTest entry `stream_preproc` | nightly (`ctest -L highlevel`) |
| `abinitio2D_stream` | retired (below) | - |

Why they went where they did:

- The sieve test runs `ptcl_sieve` directly, in-process, on two 64² class averages. The other
  `ptcl_sieve` unit tests were already in `unit_project`, so it joined them in the fast gate.
- Optics assignment, picking references and pick and extract run in-process on one thread and
  start no workers. That makes them a library suite. Optics assignment takes a minute, because
  the stage imports a project only once it is `LONGTIME` = 60 s old.
- `master` never started the stream master. It tested `gui_assembler%assemble_stream_heartbeat`
  with seven real forked children. Forking is excluded from the fast gate, so it runs beside the
  forked-process lifecycle tests in the `forked_process` platform entry. Its name is now
  `stream heartbeat`. Since the master refactor it forks one child and reads it through
  `fork_gui_status` (see the heartbeat section below); the seven-child form is gone.
- `preproc` submits its jobs to the local queue, so it stays a workflow entry of its own.

Running them:

```text
simple_test_exec test=unit_project suite=particle_sieve       # the sieve tests, with collect-and-reject
simple_test_exec test=lib_stream                              # the three in-process stages
simple_test_exec test=lib_stream suite=pick_and_extract       # one of them
simple_test_exec test=forked_process suite=stream_heartbeat   # the heartbeat
simple_test_exec test=preproc                                 # the preprocessing workflow
```

Each suite runs in its own dated directory and writes its report there
(`SIMPLE_TEST_<suite>_<date>/simple_test_<suite>_report.txt`). A `lib_stream` test removes its
fixture directory when all its checks pass. When one fails, it keeps the directory and prints its
path.

## What changed during the move

- **Assertions.** A failed `THROW_HARD` ended the process, so one broken check hid all later ones.
  Now every check is counted and reported, and the suite fails at the end. A check that reads a
  file is skipped when an earlier assertion found the file missing.
- **Fixture directories.** Every run used to leave a `test_<name>_<pid>` directory behind. Now
  the `lib_stream` tests delete theirs when they pass, and the sieve test uses the sieve tester's
  workspace, which is always removed.
- **Seeds.** `parameters%new` calls `seed_rnd`, which read `/dev/urandom`. So every test that runs a
  commander drew unseeded numbers from its first commander on: the movie simulator's noise and
  positions, and class initialisation. `seed_rnd` now reads the environment variable
  `SIMPLE_SEED`. When it is set, the seed is that integer advanced by 7919 per call, so successive
  commanders draw different but reproducible numbers, and distributed workers inherit it. When it
  is unset, production behaves as before. `simple_test_exec` sets `SIMPLE_SEED=20260923` for its process. The
  suite runner also reseeds before every sub-suite, so `suite=<name>` draws the same numbers as the
  full run, and a run by hand draws what a CTest run draws.
- **`abinitio2D_stream` is retired.** It ran one iteration of `solve2D` on 24 noise-free
  particles. Its command line was hand-written and differs from the one the stream's chunk code
  builds (`cls_init` ptcl vs rand, `rank_cavgs` no vs yes, no `chunk=yes`, `objfun` default vs
  euclid, `refine` snhc_smpl vs prob_snhc). It checked counts and files, but never that the two
  particle families end up in different classes. `solve2D` itself runs nightly in both
  `simulated_workflow` systems. If you want a test of the chunk path, see the last section.

## What a test has to do here

The rules are the same as for SINGLE (`single_area_tests_handover.md`, "What a test has to do
here"):

- Assert through `simple_test_utils`, one assertion per guarantee, with a message that states it.
- Keep `THROW_HARD` for broken fixtures.
- Derive expected values independently of the code under test, and name each tolerance.
- Seed every draw a test makes itself with `set_fixed_seed(<n>)`.
- Never assert on timings.
- Keep fast sub-suites well under a second.

## What each test should check

### `particle sieve`, collect and reject (fast)

This test is already strict: one class kept, one blank class rejected, exact counts, sentinels,
previews and the latest-product metadata. Two additions are still open:

1. **The rejection reason (partly done).** The sieve now stores the reason without blanks
   (`reason_token` in `simple_ptcl_sieve.f90`, for example `coarse_reject:` followed by the reason
   text with underscores), so the full reason round-trips through the orientation reader. The test
   still asserts only that the value contains `coarse_reject:` and that something follows it. Open:
   assert the whole reason for the blank class (NO_COMPONENT).
2. **The other hard gates.** Add one class for each other hard gate that
   `evaluate_cavg_quality_hard_reject` applies in the sieve context
   (`CAVG_QUALITY_CONTEXT_SIEVE`). Assert that each is rejected, with its own reason. The test
   fixture still has only the one blank class.

Keep it small. It runs in the fast gate.

### `optics assignment` (lib_stream)

This one checks exact truth. Missing:

1. **A one-cluster control.** Five micrographs with shifts within `tilt_thres` of each other
   should give one optics group of five, with the centroid at their mean.
2. **The threshold.** Two clusters just closer than `tilt_thres` should give one group, and just
   farther apart should give two.
3. **Restart.** Run the stage a second time into the same `outdir`. `map_count` continues from
   the latest optics map (`get_latest_optics_map_id`), so the second run writes `optics_map_2`,
   and the groups should come out the same as in the first run.
4. **The minute it costs.** `LONGTIME` (60 s) is hard-coded in the stage's watcher. If Joseph
   agrees to make the watcher's age a parameter, this test drops to a few seconds and could move
   to the fast gate.

### `picking references` (lib_stream)

The test's description promised "normalized rotation and mirror outputs", but it checks only
counts and metadata. The source images are synthetic, so the truth is known:

1. **Rotations and mirrors.** Output reference (i, r, m) must equal output (i, 0, 0) rotated by
   r·360°/`nrots` and mirrored when m = 1. Test with a normalised correlation of at least 0.99
   against your own rotation and mirror of (i, 0, 0). Read the stack order from `make_pickrefs`
   and write it into the test's comment.
2. **Normalisation.** State what `make_pickrefs` guarantees (for example, mean 0 and sd 1 inside
   the mask), and assert it for every reference.
3. **Diameter.** Source i is `square(5+i)` plus `square(2+i)` shifted by (7+2i, −5+i) pixels, so
   its extent follows from how `image%square` is defined. `diam_max` (in Å) should be within two
   pixels of that extent for the largest source. `box_for_pick` should be at least
   `diam_max`/`smpd`, if that is the rule; write the rule into the test.

### `pick and extract` (lib_stream)

`nboxes_max=3` caps the picks at exactly the number the test asserts, so over-picking cannot fail
it, and the positions are never compared with `PARTICLE_COORDS`. Four changes:

1. **Remove the cap,** or set it well above three (10). Assert exactly three picks.
2. **Positions.** Read the box file and match each box to one placed particle within 2 pixels,
   one to one. Check whether `add_window` takes the window's corner or its centre, and which one
   the box file stores. Compare in one convention and write it into the test's comment.
3. **Content.** Every extracted particle should correlate with the reference at 0.9 or more
   (normalised correlation, inside a mask).
4. **A negative control.** A noise-only micrograph gives no picks, or the micrograph is rejected
   (state 0). Assert whichever the stage guarantees.

### `stream heartbeat` (forked_process)

The forked-process test (`test_fork_gui_status`, `simple_stream_master_tester.f90`) forks one
child and covers only the running state and then the finished state of `fork_gui_status`. The
heartbeat assembly itself (`test_stream_heartbeat`, `simple_gui_assembler_tester.f90`, in `unit_ui`) is
fed status records, not live processes: it covers all-running, all-finished, and one failed
record beside running ones (stage `failed`, master `error`). The master's aggregate status has more
branches (`assemble_stream_heartbeat`), and `forked_process` can produce each of them. Open:

1. **Failed and error from a live child.** `kill` one child (SIGKILL) while the others run: that
   stage reports `failed` and the master reports `error`. No heartbeat test kills a child; the failed
   branch is reached only through a hand-made record. When every child has failed or stopped,
   with at least one failed and none running, the master reports `failed`; that all-failed case
   is not tested either.
2. **Skipped in the heartbeat.** `test_skipped_stays_skipped` (`stream master` suite) checks that a
   skipped stage stays skipped and forks nothing, but the heartbeat tests never feed a skipped
   status: assert that it reports `skipped` and does not count as running.
3. **Shared process.** `initial_picking` and `opening2D` report the same pid, because they share
   one process in the master. No test compares the two.

### `preproc` (workflow entry `stream_preproc`)

**Delivered.** The test now keeps `simulate_movie_params.txt` and `optimal_movie_average.mrc`
beside each movie (numbered per movie), matches micrographs to movies by file name, and fails on:
recovered defocus more than 0.10 µm from the simulator's recorded value (`DEFOCUS_TOL`; the two
astigmatic axes may be exchanged), recovered frame motion more than 0.50 pixels from the simulated
shifts (`MOTION_TOL`), and an integrated micrograph whose correlation with its own optimal
average is below a floor (`CORR_FLOOR`, 0.50) or does not exceed its correlation with another
movie's average by `CORR_MARGIN` (the wrong-movie control). Whether those checks pass on a given machine is for the
nightly run to say.

Workflow entries end on the first failed `THROW_HARD`, which CTest reports as a failure, so the
test can keep its style. If the checks grow, convert them to assertions like the others.

### If you want the chunk path tested (`abinitio2D_stream`, retired)

A useful replacement tests what the stream actually sends to `solve2D`:

1. Build the chunk command line with the chunk code. `simple_stream_chunk2D_utils` no longer
   exists (the chunk code moved to `src/main/solve/simple_solve2D_chunk.f90`, which takes its
   command line from its caller, and the 2D pool is in `src/main/stream/pool2D`); whichever
   routine builds the stream's chunk command line would have to be factored into a function a
   test can call (Joseph's code). This replacement is optional.
2. Use two particle families with noise (seeded) and enough iterations to converge.
3. Assert class purity: at least 90 % of each family in one class, and the two classes different.

At that size it belongs in `lib_stream`.

## Chained tests (5 October 2026)

Decision 19 of `doc/refactoring_notes/completed/stream_fix_plan_2026-10-05.md` adds two chained
tests, in `src/main/stream/simple_stream_chain_tester.f90`, sub-suites of `lib_stream`:

| sub-suite | stages | fixture |
|---|---|---|
| `sieve to 3D` | p05, p06, p07 | simulated particle stacks (a three-blob truth volume, CTF, noise) in completed reference-picking sets, `moldiam.txt`, and reference picking's finished marker |
| `movies to 3D` | p01, p02, p04, p05, p06, p07 | simulated movies of reprojections of the truth, picking references reprojected from it; preprocessing is finalised once every movie is processed |

Both drive the production stage types in one process, without the master: each stage is made
with its own `new()` on its commander's command-line defaults (the commanders' `set_*_cline`
routines are public for this), in a folder of its own, and the driver gives the stages a pass in
turn until multistate 3D completes its first run or three hours pass. The stages submit their
jobs to the local queue, so they need the installed executables (`SIMPLE_PATH`), like
`stream_preproc`, and take tens of minutes. They check that each stage made its product: the
sieve accepts particles and ends its intake, the pool takes the final set and publishes, 3D
completes a run with a resolution. They were first written without a build. Fixture fixes since
then landed in `c03b029fc` and `75cba3c19` (2026-10-06/07) and `eed4f9609` (2026-10-13, stage 4
after class-average initialization); the Phase 5 note
(`phase5_workflow_gates_and_nightly_runner_handover.md`) records them. The sub-suites carry
assertions and a three-hour limit (`TIME_LIMIT_S`) and run inside `lib_stream` (timeout 7200 s
since 2026-10-13). This note does not record a passing run of either; that evidence belongs to
the nightly run.

```text
simple_test_exec test=lib_stream suite=sieve_to_3d
simple_test_exec test=lib_stream suite=movies_to_3d
```

