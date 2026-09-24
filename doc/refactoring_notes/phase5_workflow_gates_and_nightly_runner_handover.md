# Phase 5 of the test environment: handover to Ruben

Since 2026-09-25 this is the live document of the test environment: the plan is archived
(`doc/refactoring_notes/completed/`), Phase 5 and the checks below are what is left of it, and
the day-to-day rules are `doc/policies/test_environment_policy.md`.

Hans, 2026-09-24. Phase 5 of the test-environment plan is yours
(`doc/refactoring_notes/completed/uniform_test_environment_refactoring.md`: section 12, Phase 5;
section 5.2.2 on the workflow gates; section 5.2.3 on the nightly run). It has two parts:

- **Part A:** turn the simulated workflows into gates on the truth they were simulated from.
- **Part B:** design and write the nightly runner that runs the extensive tier unattended on
  the dedicated machine and archives what it found.

Phases 0 to 3 are done and the build-time gate is in place. This note describes what exists,
what each workflow should compare with its simulation, and what the runner has to do. It
builds on your two earlier handovers (`single_area_tests_handover.md`,
`stream_area_tests_handover.md`); their workflow items belong to Part A.

## What exists today

Every test runs through `simple_test_exec`; there is no other test executable. CTest registers
25 entries, each in its own working directory `build/test_runs/<entry>`. Every entry gets
`SIMPLE_SEED=20260923`, which makes `seed_rnd`, and so every commander, draw reproducible
numbers.

| label | entries | run by |
|---|---|---|
| `fast` | 13 area suites `unit_<area>`, about 5 s together | every `compile_*.sh --compile-tests` build (`scripts/run_fast_gate.sh`) |
| `library` | `lib_reconstruction`, `lib_cart_align3D`, `lib_heterogeneity`, `lib_single`, `lib_stream` | nobody yet: the nightly run |
| `workflow` | `simulated_workflow_6vxx`, `simulated_workflow_1jxy`, `single_workflow`, `pcg_recon`, `simulate_particles`, `stream_preproc` | nobody yet: the nightly run |
| `platform` | `forked_process`, plus `coarrays`, `flex_gpu` and `openmp_offload` when the build has the capability | by hand, or by the nightly run where the machine has the capability |

None of the library or workflow entries has run in its current form, so none has a recorded
runtime. The first nightly run records them.

The number of CTest entries is fixed in `SIMPLE_CTEST_BUDGET` (`production/CMakeLists.txt`);
configuration fails when it doesn't match. Adding an entry, for example the second picker of
`simulated_workflow`, needs Hans's agreement and a raised budget with the reason recorded in
the plan. Checks grow inside existing entries.

## Part A: gates on the simulation truth

What the workflows check today is that the pipeline completes: files exist, counts match,
classes are produced. None of them compares its result with the model that generated the data.
For each workflow, the truth is already available at no cost, because the data comes from
embedded coordinates (6VXX, 1JYX, the Pt nanoparticle) with known orientations, defocus, shifts
and positions.

The rules (plan, section 5.2.2):

- **Declared floors.** Every floor is declared in the test, next to where its value comes from:
  the first measured runs plus a stated margin, or a physical argument. It is loosened only with
  a written justification in the commit. There are no blessed baselines, no platform keys and
  no report registry.
- **Failing.** A failed floor makes the entry exit non-zero. A workflow entry may keep failing
  through `THROW_HARD`, as today, or collect its checks with `simple_test_utils` and fail at the
  end. Either way, the message must name the metric, its value and the floor.
- **Metrics file.** Each workflow writes its metrics to a small machine-readable file in its run
  directory, for example `metrics.tsv` with one `name value floor pass` line per metric. The
  nightly runner collects these (Part B). Define the format once and use it in every workflow.
- **Truth independent of the code under test.** For example, compare the final map with the map
  `pdb2mrc` makes from the same coordinates, not with a map from an earlier stage of the same
  run.

### `simulated_workflow_6vxx`, `simulated_workflow_1jxy`

`exec_test_simulated_workflow` (`simple_commanders_test_highlevel.f90`) runs these stages:
- `pdb2mrc`: the truth map;
- `reproject`: 100 projections with known orientations;
- `simulate_movie`: 10 movies of 10 particles, with a defocus per movie and a frame drift per
  movie in `simulate_movie_params.txt`;
- import, motion correction, CTF estimation;
- picking, `segdiam` by default;
- extraction, `abinitio2D`, `abinitio3D`.

What to add:

1. **CTF.** The estimated `dfx`/`dfy` of each micrograph against the defocus it was simulated
   with, matched by movie name.
2. **Motion.** The estimated frame shifts against the simulated drift, after removing the
   reference frame's offset.
3. **Picking.** Recall and precision of the picked coordinates against the positions the
   simulator placed the particles at, within a stated radius. Write the positions out if
   `simulate_movie` doesn't keep them yet.
4. **Map.** The ab initio map is neither docked onto the truth map nor necessarily of the
   right hand, so the comparison is: same grid, dock, choose the hand, then FSC. Every step
   exists in SIMPLE; this is how to use them.
   - **Same grid.** Read the box and sampling of the final `abinitio3D` volume from its header
     (`find_ldim_nptcls`, `find_img_smpd`); it may be downscaled relative to the particles. Make
     the truth map on that grid with `atoms%pdb2mrc(smpd=<that smpd>, vol_dim=[box,box,box],
     mol=<the system's molecule_data>)`; it moves the atomic centre to the box centre when it is
     off-centre, and stops if the box is smaller than the molecule. Don't reuse the
     `0_pdb2mrc` map of the same run: it is on the particle grid.
   - **Dock.** Use `dock_vols` (`simple_dock_vols`, the engine behind `simple_exec
     prg=dock_volpair`). Call `new(truth, target, smpd, hp, lp, mskdiam)`, then `srch()`,
     `get_dock_info(eul, shift, cc)` and `rotate_target(target, docked)`. Search at low
     resolution: `lp` around 15–20 Å and `hp` well below the particle size. Use the workflow's
     `mskdiam`. The rotation search covers the whole sphere in C1, and the shift search a sixth
     of the (clipped) box, so a symmetry-equivalent or off-centre ab initio map docks too. The
     two volumes must have identical dimensions (`dock_vols` stops otherwise).
   - **Hand.** Write the mirror of the target with `image%mirror('x')`, the same flip that
     `volops mirr=x` and `postprocess` apply. Dock both the map and its mirror. Keep the one
     with the higher docking `cc`, and record which one it was as a metric (`hand_flipped` 0
     or 1). Choose by the docking score, not by the FSC under test.
   - **FSC.** Apply the same soft spherical mask (`mask3D_soft` with the workflow's `mskdiam`)
     to the truth map and the docked map. Fourier transform both, compute `truth%fsc(docked,
     fsc)`, then `get_resolution(fsc, res, res_fsc05, res_fsc0143)`. These are the steps the
     `fsc` commander (`simple_commanders_resolest`) takes with spherical masking.
   - **Floor.** The resolution at FSC 0.143 must be at or better than a floor per system. Also
     record the resolution at FSC 0.5 and the docking `cc`.
5. **Poses, if the bookkeeping allows it.** Once picks are matched to simulated particles,
   compare the recovered orientations with the projections they came from:
   - first compose each recovered orientation with the docking rotation;
   - apply the mirror too if the hand was flipped;
   - then take the minimum distance over the point group (`sym%sym_dists`);
   - report the fraction within a stated angle.

   If that matching is too fragile, say so and leave it out; the map FSC is the primary gate.
6. **A defect to fix while you're there.** Every stage runs with `NTHR = 4`, a constant in the
   routine, whatever the entry's 8 threads. `single_workflow` had the same defect (`nthr=40`),
   and it now takes `params%nthr`.
7. **The second picker.** The plan expects both pickers to run. Propose it as an extra entry or
   as a second pass inside the same entry (no new process), and ask Hans.

### `single_workflow`

This is listed in `single_area_tests_handover.md`: stage sizes, map and orientation accuracy
against the simulated nanoparticle, and an atomic-model check. Please also confirm the
`smpd=0.358` the CTest entry passes.

### `stream_preproc`

This is listed in `stream_area_tests_handover.md`, section "`preproc`":
- keep `simulate_movie_params.txt` and the optimal average instead of deleting them;
- defocus of movie i against 1.5 + 0.25·(i−1) µm;
- the estimated frame shifts against the simulated ones;
- the integrated micrograph against the optimal average.

### `simulate_particles`

It checks stack presence, image count, box, sampling and one orientation record per image.
Add:
- the CTF parameters recorded per particle equal the ones applied;
- the noise level matches the requested SNR within a stated tolerance, measured on a copy
  simulated without signal or from the known signal variance;
- the recorded orientations equal the input orientations.

### `pcg_recon`

This one already gates 14 stages of the PCG operator with fixed seeds. Record its runtime. If
one thread runs it in seconds, tell Hans: it could move to the fast gate. That would need its
checks rewritten on `simple_test_utils`.

### Not in scope

`mini_stream` is a manual test by review verdict: it needs a user's movies, gain reference and
cluster. It becomes a gate only if a simulated movie set replaces the user data. That is a
separate proposal if you want it.

## Part B: the nightly runner

Requirements, from plan section 5.2.3 and what the entries need:

1. **Unattended and scheduled** (cron or launchd) on the dedicated machine, with a lock so that
   two runs never overlap, and with a clear failure when the machine or checkout isn't in the
   expected state.
2. **A clean build of a known commit.**
   - Fetch, check out the branch to test, and run `./compile_clean.sh --compile-tests`, which
     builds and then runs the fast gate. A failed gate stops the night and is the first line of
     the report.
   - Record the commit, host, compiler version (`gfortran --version`) and CMake options.
3. **The extensive tier.**
   - Run `ctest -L library` first; the library suites can run side by side.
   - Then run `ctest -L workflow`; the workflow entries are `RUN_SERIAL`, so each owns the
     machine.
   - Add `ctest -L platform` where the machine has the capability.
   - Repeat the fast tier threaded, outside CTest: its `fast` entries pin `OMP_NUM_THREADS=1`
     in their CTest environment, which a variable set in the shell does not override. For each
     `unit_<area>` (`ctest -N -L fast` lists them), in a scratch directory with the environment
     the entries get (`SIMPLE_PATH`, `SIMPLE_SEED=20260923`):
     `OMP_NUM_THREADS=4 SIMPLE_UNIT_ORDER=reverse simple_test_exec test=unit_<area>`; a nonzero
     exit fails the night. The unit suites take their thread count from the environment, so
     this is the only run that exercises the threaded paths of the fast suites (seconds).
   - Optional, weekly rather than nightly: an instrumented build (`--coverage`) of the fast and
     library tiers with a `gcovr` summary in the archive, as a report and not a gate; it shows
     which production modules no test reaches.
   - Keep ctest's per-entry status and time: `--output-junit` in recent CMake, or parse the log
     as `scripts/ctest_budget.py` does for the fast gate.
4. **The summary.**
   - Collect every entry's status and time, and every `metrics.tsv` from
     `build/test_runs/<entry>/`.
   - Write a dated directory **outside `build/`** (the compile scripts delete `build/`), for
     example `~/simple_nightly/<date>_<commit>/`.
   - The directory holds `summary.md` (or `.txt`): status, times, metrics against floors, and the
     failing entries with the tail of their logs. It also holds the raw ctest log and the
     metrics files.
   - Append one line per entry and metric to a cumulative history file (CSV or TSV), so that a
     regression can be dated and a growing runtime seen.
5. **Notify where the team looks** (mail, a chat channel or a shared directory; agree it with
   Hans). The message is short: pass or fail, what failed, the link or path to the summary.
6. **Fit the night.**
   - Each entry keeps its CTest `TIMEOUT`, and the runner has an overall limit.
   - The summary flags an entry that grew past its share compared with the history.
   - Trimming such an entry is a reviewed change, as for the fast gate.
7. **Deliverables.**
   - The runner script(s) in `scripts/`, for example `nightly_run.sh` plus a summary writer.
   - A short install note: machine, schedule, where the archive lives, how to run a night by
     hand.
   - What was built, recorded in this handover (the plan is archived and no longer updated).
   - After the first complete night, the measured runtimes of every library and workflow entry,
     recorded here; the generated test inventory shows run times when the dossier script is given
     the night's timing file (`--timing`).

## What is still weak in the nightly tier (from the earlier handovers)

- `lib_single`: `nanoparticle atoms` and `C-alpha molecules` assert nothing yet. The report shows
  them as "completed" and they cannot fail.
- `lib_stream`:
  - picking references don't check rotations, mirrors or normalisation;
  - pick and extract caps the pick count at the number it asserts;
  - the stream heartbeat lacks the failed, error, skipped and restarting states.

These are not Phase 5 proper. The runner will run them every night from now on, which is the
reason to close them early.

## Outstanding checks, carried over from the plan

The plan closed with these observations still to make (its section 16, criterion 10: what has not
been observed is listed, not claimed). Record each here when it is made.

- A green `./compile_debug.sh --compile-tests` build of the tree with the in-module self-test
  batch and the build-time code map and test inventory (Hans).
- One `simple_test_exec test=pcg_recon`: its seed changed on 2026-09-25, and the run also gives
  its runtime (Hans).
- One build without `--compile-tests` (`BUILD_TESTS=OFF`), to confirm it links without the
  test-only sources (Hans).
- The offload branch of `simple_openmp_offload_tester` in an offload build (Cyril).
- The first night of the runner, with the runtimes of the library and workflow entries (Ruben;
  the exit of Part B).

## Who to ask

Hans reviews floors and any new CTest entry. Joseph owns the stream production code that the
stream gates exercise.
