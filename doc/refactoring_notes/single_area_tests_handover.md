# SINGLE tests in the new test environment: handover to Ruben

Hans, 2026-09-23. Your six SINGLE test programs (`simple_test_exec test=...` in
`src/main/commanders/test/simple_commanders_test_single.f90`) have been moved into the
two-tier test environment (plan: `doc/refactoring_notes/uniform_test_environment_refactoring.md`,
section 9.7, "single"). They were moved as they were. Only one of them checked anything
beyond "it did not crash", so most of them pass whatever the code does. This note says
where each one lives now and what it has to check before it counts as a test.

## Where the tests are now

| was | now | runs |
|---|---|---|
| `detect_calpha` | `src/main/nano/simple_calpha_finder_tester.f90`, sub-suite `C-alpha finder` of `unit_single` | fast gate, every `--compile-tests` build |
| `simulate_nanoparticle`, `detect_atoms`, `atoms_stats` | one pipeline, `exec_test_atoms_stats`, run as sub-suite `nanoparticle atoms` of `lib_single` | nightly (`ctest -L library`) |
| `detect_calpha_molecules` | unchanged commander, run as sub-suite `C-alpha molecules` of `lib_single` | nightly |
| `single_workflow` | unchanged commander, CTest entry `single_workflow` | nightly (`ctest -L workflow`) |

`unit_single` also holds the `atoms` sub-suite (`test_atoms` in `simple_atoms`), which used to sit in
`unit_project`. `simulate_nanoparticle` and `detect_atoms` were the first stages of `atoms_stats`,
so they are no longer separate cases. Their stages still run, inside `atoms_stats`.

Running them:

```text
simple_test_exec test=unit_single                          # both fast sub-suites
simple_test_exec test=unit_single suite=c_alpha_finder     # one sub-suite
simple_test_exec test=lib_single suite=nanoparticle_atoms  # the nanoparticle pipeline
simple_test_exec test=atoms_stats smpd=0.358 element=Pt    # the same pipeline by hand
simple_test_exec test=detect_calpha_molecules [smpd= angstep= thres=]
simple_test_exec test=single_workflow smpd=0.358 element=Pt
```

Each suite runs in its own dated directory and writes its report there
(`SIMPLE_TEST_<suite>_<date>/simple_test_<suite>_report.txt`).

## Two defects fixed during the move

- The CTest entry `single_workflow` passed only `element=Pt`. `smpd` is a required key of the
  program, so `simple_test_exec` printed the usage line and stopped with status 0. The entry
  "passed" without running anything. It now passes `smpd=0.358`. **Please confirm 0.358 Å, or
  give the sampling you want.** `lib_single` uses the same value for the nanoparticle pipeline.
- Every stage ran with `nthr=40` (module constant `NTHR`), whatever the machine or the
  CTest entry's 8 threads. The stages now take `params%nthr`, which defaults to `OMP_NUM_THREADS`.

## What a test has to do here

- Assert through `simple_test_utils` (`assert_true`, `assert_int`, `assert_real`, ...), one
  assertion per guarantee, with a message that states the guarantee ("every simulated atom is
  detected within 0.5 A"). The suite runner counts the assertions, reports them and fails the run.
  A sub-suite with no assertions shows up in the report as "completed" and can never fail. That
  is the state of `lib_single` today.
- Do not use `THROW_HARD`, `stop 'message'` (which exits with status 0) or printed "PASS"/"FAIL"
  lines as checks. Keep `THROW_HARD` for broken fixtures only.
- Derive every expected value independently of the code under test: closed forms, the
  generating model, or a Python/numpy emulation. Name the tolerances and say where each one
  comes from.
- Make it reproducible. Seed every random draw with `set_fixed_seed(<n>)` from
  `simple_test_utils`. `seed_rnd`, which every commander calls through `parameters%new`, reads
  `/dev/urandom` unless the environment variable `SIMPLE_SEED` is set; CTest sets it for every entry
  (since the stream review), so export `SIMPLE_SEED=20260923` to reproduce a CTest run by hand.
  Never assert on timings.
- Size: fast sub-suites (`unit_single`) must stay well under a second each. The whole fast gate
  is budgeted at 30 s on one thread. Library and workflow suites can take minutes.
- Write files in the current directory (the suite's directory) and delete the ones you don't
  need to keep.

## What each test should check

### `C-alpha finder` (fast, `simple_calpha_finder_tester.f90`)

Today it asserts only that some candidate exists, and that some candidate lies within 1.5 Å of
some residue centre. One good peak out of three residues passes. The map is synthetic and the
truth is known exactly, so this test can be strict:

1. Every residue's Cα (the three `centers`) has its own candidate within 1.5 Å. Match candidates
   to Cα one-to-one by nearest distance; don't count one candidate twice.
2. No candidate is closer to an N or C site than to that residue's Cα.
3. The number of candidates is at most `npeaks` (10). Every candidate's score is at or above
   `score_threshold`, if the PDB or CSV carries the score.
4. A negative control: an empty map, or one Gaussian blob without the backbone geometry,
   gives no candidate above the threshold.
5. Optionally, orientation independence: apply a different `rotation` to the residues and
   require the same recovery. It is cheap at 32³.

### `nanoparticle atoms` (nightly, `exec_test_atoms_stats`)

Today it simulates a Pt nanoparticle (box 160, diameter 20 Å), runs `detect_atoms` and
`atoms_stats`, and checks nothing. The simulation writes its own ground truth
(`simatms.pdb`) next to the map (`outvol.mrc`), so everything can be pinned:

1. Simulation: `simatms.pdb` and `outvol.mrc` exist, and the map has the requested box and
   sampling.
2. Detection on the noise-free map (`outvol_ATMS.pdb` against `simatms.pdb`): the detected atom
   count equals the simulated count, or is within an agreed 1%. After one-to-one
   nearest-neighbour matching, every detected atom lies within 0.5 Å of its simulated atom. Report
   and bound the RMS position error.
3. Statistics (`atoms_stats`, CSV files from `simple_nanoparticle%write_csv_files`): the
   nearest-neighbour distance distribution peaks at the Pt value. For fcc Pt with a = 3.92 Å that
   is a/√2 = 2.77 Å; assert the peak within a stated tolerance. Also check any other statistic
   whose true value the simulation fixes (atom count, diameter).
4. A noisy variant (add Gaussian noise at a fixed seed and a stated SNR), with looser
   recall and position floors. That is the case that shows whether detection works on data.

### `C-alpha molecules` (nightly, `exec_test_detect_calpha_molecules`)

A benchmark on the built-in 6VXX and 1JYX models. It prints top-N and top-2N recall and precision,
and stops (via `THROW_HARD`) only if the built-in Cα counts (2916, 4044) change.

1. Turn the Cα counts into `assert_int`.
2. Set floors for top-N recall and precision on both molecules at the default settings
   (`smpd=1.3`, `angstep=45`, `thres=0.25`). Take them from one measured run, minus a stated
   margin, and assert them. The floors are what protects the finder from regressions. Without
   them the benchmark reports numbers that nobody checks.
3. If the finder is still being tuned, keep the floors conservative and raise them as it improves.
   Record each change of floor in the commit message.

### `single_workflow` (nightly workflow)

It simulates a Pt nanoparticle and 5000 reprojections along a GLC-like trajectory. It adds noise,
denoises the trajectory, imports the particles, and runs `analysis2D_nano` and
`autorefine3D_nano`. The only check is that `analysis2D_nano` wrote `startvol.mrc`. The generating
model and the orientations (`glc_trajectory_oris.txt`) are known, so the pipeline can be held
to them:

1. Every stage produced its output with the right size. There are 5000 reprojections in the stack,
   in the trajectory and in the denoised trajectory, and 5000 particles in the project after the
   import.
2. The final 3D map against the simulated volume (`1_simulate_nanoparticle/outvol.mrc`): an FSC
   resolution, or a real-space correlation after alignment, better than a floor.
3. The refined orientations against the generating ones: median angular error below a floor,
   modulo the point-group ambiguity.
4. The atomic model from the refined map against `simatms.pdb`, as in `nanoparticle atoms`, with
   floors appropriate to the resolution.

These are the "simulation-truth gates" of Phase 5 of the plan. `single_workflow` can be the first
workflow to have them.
