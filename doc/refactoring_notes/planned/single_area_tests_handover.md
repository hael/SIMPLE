# SINGLE tests in the new test environment: handover to Ruben

Hans, 2026-09-23. Updated 14 October 2026 against the current source.

Your six SINGLE test programs (`simple_test_exec test=...` in
`src/main/commanders/test/simple_commanders_test_single.f90`) were moved into the
two-tier test environment (plan: `doc/refactoring_notes/completed/uniform_test_environment_refactoring.md`,
section 9.7, "single"). They were moved as they were. Only one of them checked anything
beyond "it did not crash". Most of what this note asked for has since been delivered; this
version records that briefly and keeps only the items that are still open. Whether any of the
checks passes on a given platform is a matter for the nightly run
(`phase5_workflow_gates_and_nightly_runner_handover.md`), not this note.

## Where the tests are now

| was | now | runs |
|---|---|---|
| `detect_calpha` | `src/main/nano/simple_calpha_finder_tester.f90`, sub-suite `C-alpha finder` of `unit_single` | fast gate, every build with tests |
| `detect_calpha_molecules` | deleted (`dff901ec7`, 2026-09-25); the policy table records it as deleted | - |
| `simulate_nanoparticle`, `detect_atoms`, `atoms_stats` | one pipeline, `exec_test_single_atoms_stats`, CTest entries `single_atoms_stats` (Pt) and `single_atoms_stats_cdsew`, `_cdsez`, `_cdser` (CdSe wurtzite, zincblende, rocksalt) | nightly (`ctest -L highlevel`) |
| `single_workflow` | `exec_test_single_workflow`, split into CTest entries `single_workflow_fcc` and `single_workflow_wurtzite` | nightly (`ctest -L highlevel`) |

`unit_single` also holds the `atoms` sub-suite (`simple_atoms_tester` since 2026-09-25; about 90
assertions), `cif2mrc` and `species`. `lib_single` now holds only `pdb2mrc`. The old
`nanoparticle atoms` sub-suite of `lib_single` no longer exists, and neither does the `workflow`
label.

Running them:

```text
simple_test_exec test=unit_single                          # the fast sub-suites
simple_test_exec test=unit_single suite=c_alpha_finder     # one sub-suite
simple_test_exec test=single_atoms_stats smpd=0.358 element=Pt   # element=CdSeW, CdSeZ, CdSeR for the others
simple_test_exec test=single_workflow suite=fcc smpd=0.358       # suite=wurtzite for the CdSe fixture
```

Each suite runs in its own dated directory and writes its report there
(`SIMPLE_TEST_<suite>_<date>/simple_test_<suite>_report.txt`).

## Delivered since the handover

- **The two defects of the move are fixed.** The CTest entries pass `smpd=0.358` (they used to omit
  a required key and "passed" without running); the stages take `params%nthr` instead of the
  module constant `NTHR=40`. **Please still confirm 0.358 A, or give the sampling you want.**
- **`atoms` sub-suite** in `unit_single`: assertions accumulate, the former run-only parts assert,
  PDB files are removed (see "Changed after this handover" below).
- **`atoms_stats` pipeline** (`single_atoms_stats`, four entries). It now fails on: simulated and
  detected-atom volumes (dimensions, sampling, finite non-constant density), simulated/detected
  density correlation at least 0.80, detection recall and precision at least 0.90, atom RMS
  position error at most 0.75 A, the statistics CSV files (atom count against the detected PDB,
  diameter, coordination numbers, per-coordination groups), and for the binary crystals the
  interior coordination and the fitted lattice constant within 2 % of the table. It uses a local
  `fail()` and `THROW_HARD` rather than `assert_*`. The position tolerance is **1.0 A**, not the
  0.5 A asked for here.
- **`single_workflow`** (`fcc`, `wurtzite`). It compares the final map with the simulated volume
  (soft-masked correlation at least 0.90, FSC 0.143 at 5 A or better), and checks the pair pose
  error against the generating orientations: the `refine=cont` run and the polish may not be
  worse than their predecessor by more than 0.5 degrees. The orientation check asked for as item 3
  below is therefore present in this relative form only. It ends with `THROW_HARD` on failure.
- **The `C-alpha molecules` benchmark is gone.** The request for top-N recall and precision floors on the
  built-in 6VXX and 1JYX models is moot, and the floors were never set: the surviving C-alpha finder
  suite has no recall or precision floors.

## What a test has to do here

- Assert through `simple_test_utils` (`assert_true`, `assert_int`, `assert_real`, ...), one
  assertion per guarantee, with a message that states the guarantee ("every simulated atom is
  detected within 0.5 A"). The suite runner counts the assertions, reports them and fails the run.
  A sub-suite with no assertions shows up in the report as "completed" and can never fail.
- Do not use `THROW_HARD`, `stop 'message'` (which exits with status 0) or printed "PASS"/"FAIL"
  lines as checks. Keep `THROW_HARD` for broken fixtures only. (The two nanoparticle workflows above
  still end on `fail()` / `THROW_HARD`; a workflow entry ends on the first failure, which CTest
  reports as a failure. Converting them is optional.)
- Derive every expected value independently of the code under test: closed forms, the
  generating model, or a Python/numpy emulation. Name the tolerances and say where each one
  comes from.
- Make it reproducible. Seed every random draw with `set_fixed_seed(<n>)` from
  `simple_test_utils`. `seed_rnd`, which every commander calls through `parameters%new`, reads
  `/dev/urandom` unless `SIMPLE_SEED` is set, and `simple_test_exec` sets it to 20260923 for its
  process, so a run by hand draws what a CTest run draws.
  Never assert on timings.
- Size: fast sub-suites (`unit_single`) must stay well under a second each. The whole fast gate
  is budgeted at 60 s of real time (`scripts/run_fast_gate.sh`; it was 30 s until `43808e946`).
  Library and workflow suites can take minutes.
- Write files in the current directory (the suite's directory) and delete the ones you don't
  need to keep.

## Still open

### `C-alpha finder` (fast, `simple_calpha_finder_tester.f90`)

It still has two assertions: that some candidate exists, and that some candidate lies within
1.5 A of some residue centre. One good peak out of three residues passes. The map is synthetic and
the truth is known exactly, so this test can be strict:

1. Every residue's Calpha (the three `centers`) has its own candidate within 1.5 A. Match candidates
   to Calpha one-to-one by nearest distance; don't count one candidate twice.
2. No candidate is closer to an N or C site than to that residue's Calpha.
3. The number of candidates is at most `npeaks` (10). Every candidate's score is at or above
   `score_threshold`, if the PDB or CSV carries the score.
4. A negative control: an empty map, or one Gaussian blob without the backbone geometry,
   gives no candidate above the threshold.
5. Optionally, orientation independence: apply a different `rotation` to the residues and
   require the same recovery. It is cheap at 32^3.

Recall and precision floors on the real 6VXX and 1JYX models cannot be asked for any more (the
benchmark is deleted); if you want them back, they need a new, cheaper fixture.

### `single_atoms_stats` (nightly)

Detection recall, precision, RMS position and density correlation are covered (above). Missing:

1. **The nearest-neighbour peak.** The nearest-neighbour distance distribution of the detected atoms
   should peak at the generating value. For fcc Pt with a = 3.92 A that is a/sqrt(2) = 2.77 A;
   assert the peak within a stated tolerance. Today only the binary crystals have a lattice
   check (fitted constant within 2 %); the Pt case has none.
2. **A noisy variant**: Gaussian noise at a fixed seed and a stated SNR, with looser recall and
   position floors. That is the case that shows whether detection works on data. (The
   `species_discovery` entry renders noisy maps for its own purpose; it is not this check.)
3. Optionally, tighten the position tolerance from 1.0 A toward the 0.5 A asked for here, and say
   which of the two the detector's resolution supports.

### `single_workflow` (nightly)

The final-map comparison and the pair pose error are covered (above). Missing:

1. **Per-stage counts.** The fixture uses `NREPROJS = 1000` reprojections (not the 5000 of the
   original workflow). Assert that the reprojection stack, the trajectory, the denoised trajectory
   and the project after import each hold 1000 images/particles. Today the only check before the
   final map is that `analysis2D_nano` wrote `startvol.mrc`.
2. **An absolute orientation floor.** The pair pose error is compared between stages, not against a
   floor; a median angular error floor against the generating orientations (modulo the point-group
   ambiguity) would make it an absolute gate.
3. **The atomic model.** The atomic model from the refined map against `simatms.pdb`, as in
   `single_atoms_stats`, with floors appropriate to the resolution. This is the same open item as
   "SINGLE final atomic-model comparison" in the Phase 5 note; a separate atom-statistics test is
   not that end-to-end check.

These are the "simulation-truth gates" of Phase 5 of the plan.

## Changed after this handover (2026-09-25)

- `test_atoms` left `simple_atoms` for `src/main/nano/simple_atoms_tester.f90` (sub-suite `atoms`
  of `unit_single`). Its checks accumulate instead of stopping at the first failure, the parts
  that only ran (PDB I/O, ANISOU, density simulation, validation) now assert, and its PDB files
  are removed.
- `atom_validate` and `map_validate` cut their per-atom window one voxel off the atom: `ang2vox`
  is 1-based and `window_slim` adds one to the corner, so the atom sat at voxel `atom_box/2`, not
  at `atom_box/2+1` where `center_inbox` puts the simulated atom. The corner is now
  `ang2vox - 1 - atom_box/2`. In a numpy emulation the correlation of an atom with its own
  simulated density rises from 0.35-0.44 to 0.99, so the per-atom scores (and the beta column
  they are written to) of both routines go up; `map_validate` compares two maps through the same
  window, so its scores change only through the mask now being centred on the atom.
- Routines only `test_atoms` called are removed: `cc_res` (its sum was never initialised),
  `find_masscen` (a copy of `get_geom_center`), `rotate`, `geometry_analysis_pdb`, `get_num`,
  `does_exist`, `print_atom` and `get_atom_corr`.
