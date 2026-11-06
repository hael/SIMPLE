# Writing a test in SIMPLE: the rules in short

The rule set is `doc/policies/test_environment_policy.md`; section numbers below are its.

## What a test is (2.1-2.3)

- It can fail, and its failure reaches the exit status: assertions of `simple_test_utils`, one per
  guarantee, the message stating the guarantee ("an integer shift is the circular shift
  out(x) = in(x + s)", not "shift test 3").
- The expected value is derived independently of the code: a closed form, a brute force in the test,
  an exact array operation, or an emulation (numpy) whose numbers are pasted with a comment saying how
  they were made. Never the routine's own output: that pins the bug with the behaviour.
- Tolerances are named and derived from the arithmetic (precision, interpolation order, a standard
  error). A tolerance widened until the test passed is a finding.
- It pins what the routine guarantees, not a random start, an unconverged iteration or a timing.
- A negative control where cheap (an unrelated image does not correlate).
- Conventions are asserted: sign, handedness, origin, index base.
- Every draw from a fixed seed.
- Not a test: a program that only runs, a wrapper of commanders, a benchmark, a demonstration.
  Those are developer programs or stay in a local worktree (3.4).
- Four kinds, never mixed: identity (byte-identical outputs of two builds at one thread count,
  `test=identity`), correctness (truth metrics on simulated data against floors), robustness (a
  distribution over seeds, `test=profile`, never a single-run gate), performance (timings, never
  pass/fail).

## Where it goes (3.1-3.2)

- The area of the code it tests (`unit_<area>`), not of its author or caller.
- Fast tier (the default; every build) only if all hold: it can fail; hermetic (no network, no
  user file, no dependence on core count, load, clock or another test's files); in-process on one
  OpenMP thread (a threaded path opens its own team of three with `num_threads`); restores the
  working directory and releases module memos (`unmemoize_mask_coords`); well under a second per
  sub-suite with the smallest fixture that exercises the path; no sleeps, polls or timers.
- Otherwise `lib_<area>` (nightly, minutes). Never shrink a test until it fits the gate.

## A tester (4.1, 4.3)

- `simple_<thing>_tester.f90` beside the module, one-line `!@descr:` naming what it tests.
- A module, `private`, using `simple_test_utils`, exporting only `run_all_<thing>_tests`, which
  prints `**** running all <thing> tests ****`; each `test_<what>` prints its name first.
- When it needs private entities of the module it tests, it is a submodule
  (`submodule (simple_<thing>) simple_<thing>_tester`): the module declares
  `module subroutine run_all_<thing>_tests()` in an interface block and makes it public, the suite
  table imports it from the module. The submodule does not re-`use` a name the module has (also
  through `simple_core_module_api`), declares no procedure or constant with a module's name, and
  cannot restore an intrinsic the module shadows (`TINY` hides `tiny`). Never make an entity public
  only for a tester, except a plain module procedure the submodule calls: gfortran before 16 gives a
  private one internal linkage and the Linux build does not link (`check_tester_linkage.py`); a
  private binding, constant, variable or separate module procedure links.
- A fixture or check several tests share is a module `simple_test_<x>_fixture` in
  `src/main/commanders/test/`, never a public procedure of a tester.
- No `stop`, `error stop`, `simple_end` or `report_summary`; `THROW_HARD` only for a broken fixture.
- Register with `add_suite` in the area's `suites_<area>` table
  (`src/main/commanders/test/simple_commanders_test_class.f90`) and in the area's `suite=` list
  (`src/main/ui/simple_test/simple_test_ui_class.f90`), in table order; names at most 32 characters.
- Large objects `class(T), allocatable` and `allocate`d (compile-time policy).

## Fixtures (4.2)

- Generate them: `gauran`, `gauimg`, `square`; molecules from the two committed PDB entries
  (`simple_molecule_data`, then `pdb2mrc`); particles and movies from the simulation commanders.
- `set_fixed_seed(<n>)` after `parameters%new`, which reseeds.
- In memory where possible. Files only when I/O is the subject, inside `enter_fixture(tag, ...)` /
  `leave_fixture(...)` (kept when a check failed); files outside a fixture directory are named
  `tmp_<tester>_*` and deleted, also after a failure.
- A loop that prints through `logfhandle` diverts it to a scratch unit and restores it.

## Before handing over (4.5)

- The gate's checks: `scripts/check_descr.py`, `check_test_registry.py`, `check_flex_dag.py`,
  `check_guide_imports.py`, `check_cli_chain.py`, `check_tester_linkage.py`; then the build and the fast gate.
  The gate runs them with the build's Python (`SIMPLE_PYTHON3`, floor 3.7 in the top-level
  `CMakeLists.txt`), not the `python3` on `PATH`: standard library only, nothing newer than the
  floor (policy, section 4.5).
- `SIMPLE_UNIT_ORDER=reverse` on the area suite: sub-suites do not depend on each other.
- A threaded routine gets a check that forces a team and compares with the serial result; nothing
  else runs the fast tier on more threads.
- A reader or writer gets round trips at the sizes where its layout changes (small, odd,
  non-square; where a header or record length changes).
- A test that finds a defect: the fix lands in the same change and the test pins it. A deleted
  test: the commit says why.

## Workflow gates (6.2)

High-level gates check truth, not completion; floors come from a 20-seed profile by the floor rule
(worst minus s, median minus 3s, control plus 3s), every gate has a negative control, and every
value it decides on is a `metrics.tsv` row (`name value floor pass sense`, sense `min`, `max`,
`check`, `report` or `success`). Read section 6.2 before adding or moving a floor.
