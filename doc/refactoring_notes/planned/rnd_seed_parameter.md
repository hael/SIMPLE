# The random seed through `parameters`, not the environment

Plan, 2026-10-11, for review. Re-verified against the code 14 October 2026:
not started (`SIMPLE_SEED`, `put_simple_seed` and `simple_seed_base` are all
still in the tree, there is no `rnd_seed` key), awaiting review. The file
table and the callers below were completed on that date.

## Problem

`seed_rnd` (`simple_rnd`), which every commander calls from `parameters%new`,
reads the environment variable `SIMPLE_SEED`. When it holds an integer, the
seed is that integer advanced by 7919 for every earlier call in the process;
otherwise it is drawn from `/dev/urandom`. The tests use it to be
reproducible, and the worker processes a test launches inherit it. It is a
configuration path outside `parameters`: not parsed, not in the job
descriptions, not in the project, and honoured by every executable, so a
value left in a shell fixes the draws of any production run without a
trace.

## Decision

The seed becomes an ordinary integer key, `rnd_seed` (`seed` alone would read
like the `trail_seed` and `frozen_seed` volume keys). Unset, production draws
from `/dev/urandom` as now. The environment variable goes.

1. `simple_rnd` keeps the per-call advance, but its base is a module variable
   set by a new `set_rnd_seed(base)` (which also restarts the call count),
   not the environment. `seed_rnd_fixed` is unchanged.
2. `parameters` gets `integer :: rnd_seed = 0`, registered with the other
   integer keys. `parameters%new` parses first, then calls
   `set_rnd_seed(rnd_seed)` when the key is given and differs from the
   current base, then `seed_rnd`. A command line without the key, or with
   the current base (the copies a commander makes for its sub-commanders),
   leaves the base and the count alone.
3. `cmdline%gen_job_descr`, through which every worker's job description is
   built, adds `rnd_seed` from the module base when one is set and the
   command line has none. A worker then seeds from its own command line, and
   the seed is in the job description and the project.
4. `simple_test_exec` sets the base to 20260923 first (`fix_test_seed` in
   `simple_test_utils`, now a call of `set_rnd_seed`).
   `set_fixed_seed(base, propagate=.true.)` sets the base the same way.
   `put_simple_seed`, the `setenv` binding and `SIMPLE_SEED` in
   `simple_rnd` are deleted.

## The draws do not change

A test process today draws from base 20260923, advanced per call and
restarted by `set_fixed_seed`; a worker starts its own count at 0. With the
base held in the module and handed to the worker on its command line, both
sequences are the same, so every test draws what it draws today and no floor
moves. A user who gives `rnd_seed=` gets reproducible draws, recorded with
the job.

`parameters%new` is not the only caller of `seed_rnd` (checked by `rg`
14 October 2026), and the argument has to hold for each of them, since all
of them advance the same per-process count:

- Reseeding in the middle of a workflow or a stage: `simple_micproc.f90`
  (~284, before the k-means repeats), `simple_cluster_seed.f90` (~27,
  `gen_labelling`), `simple_relion.f90` (~508), `simple_stream_pool2D.f90`
  (~198), `simple_nanoparticle_utils.f90` (~62), `simple_ctf_estimate_fit.f90`
  (~220), `simple_commanders_solve2D.f90` (~1085) and the probabilistic
  tables (`simple_eul_prob_tab.f90` ~223, `simple_eul_prob_tab2D.f90` ~220 and
  ~342, `simple_eul_prob_tab_neigh.f90` ~157 and ~356).
- In testers: `simple_commanders_test_highlevel_workflow.f90` (~172, ~198).
- `seed_rnd_fixed` called directly, which this plan leaves unchanged:
  `simple_flex_pca_tester.f90`, `simple_flex_pca_pcg_tester.f90`,
  `simple_cls_expansion_strategy.f90` (production code, seeding per class
  from `FLEX_CLS_SEED_BASE`) and `simple_test_utils.f90`.

The plan keeps the per-call advance and the count (only the source of the
base moves), so each call draws what it draws today as long as the order of
the calls in the process is unchanged; the change to `parameters%new`
(below) must not add or drop a call. What has to be checked, not assumed, is
that every process that calls `seed_rnd` gets the base: a worker through
`gen_job_descr`, a test process through `set_rnd_seed`, and any other
process launched without a job description (to be found by reading the
launchers) would otherwise fall back to `/dev/urandom` where the
environment variable used to reach it.

Today `parameters_phases.f90` (~28) calls `seed_rnd` before `parse_inputs`
(~30); the plan needs it after parsing, so `rnd_seed` is known when the seed
is made, as decision 2 says. That is a move of one call, in the same
function.

## Files

| file | change |
|---|---|
| `src/utils/math/simple_rnd.f90` | module base and `set_rnd_seed`; `simple_seed_base` reads the module, not the environment; comments |
| `src/main/params/simple_parameters.f90`, `simple_parameters_parse.f90`, `simple_parameters_phases.f90` | the `rnd_seed` key (the args list is generated from the declaration); seeding after parsing |
| `src/main/simple_cmdline.f90` | `gen_job_descr` adds `rnd_seed` when a base is set |
| `src/main/simple_cmdline_tester.f90` | `gen_job_descr` carries `rnd_seed` from the base and keeps a given one |
| `src/utils/simple_test_utils.f90`, `production/simple_test_exec.f90` | `fix_test_seed` and `propagate` through `set_rnd_seed`; no `setenv` |
| `src/main/commanders/test/simple_commanders_test_class.f90` (~852) | the comment on `set_fixed_seed` and "the `SIMPLE_SEED` count" |
| `production/CMakeLists.txt` (~190-191) | the comment on `SIMPLE_TEST_ENV` that says `simple_test_exec` fixes `SIMPLE_SEED` |
| `doc/policies/test_environment_policy.md` (lines ~154, 484, 485, 548) | the seed described as above: `fix_test_seed`, the `propagate` row of `set_fixed_seed`, the `fix_test_seed` row, and the sentence that `SIMPLE_SEED` fixes every draw |
| `doc/refactoring_notes/planned/stream_area_tests_handover.md` (~62-64), `single_area_tests_handover.md` (~60), `phase5_workflow_gates_and_nightly_runner_handover.md` (~31) | the seed described as above |
| `doc/implementation_notes/planned/species_discovery_status.md` (~216-218), `species_discovery_report.md` (~381), `species_discovery_completed.md` (~2294, 2313, 2338, 2543) | the status note says that production `seed_rnd` still honours `SIMPLE_SEED` and that this plan awaits review: reword when this lands; the report and the note quote `SIMPLE_SEED=20260923` as the setting of their runs, which stays true as a record of those runs unless the maintainer prefers otherwise |

Not touched: the completed notes
`doc/refactoring_notes/completed/view_balanced_fractional_sampling_refactoring.md`
(and its report) and `uniform_test_environment_refactoring.md` (and
`testing_standards_assessment.md`) describe the environment-variable
mechanism as it was when they were written. They stay as they are.

The generated indexes under `doc/code_overview/fortran-indexes/`
(`module_index.md`, `api_index.md`, `symbol_index.csv`) list
`simple_seed_base`, `put_simple_seed` (and its `c_setenv` binding) and carry
line numbers of `fix_test_seed`; they are regenerated after the change, not
edited by hand.

## Checks

- `scripts/check_test_registry.py` and `scripts/check_descr.py`.
- Build, then `ctest -L fast`, `ctest -L library` and `ctest -L highlevel`
  on the Mac: the reports and `metrics.tsv` values equal those of the build
  before the change.
- `grep -rn SIMPLE_SEED src production scripts .github` finds nothing (today
  it finds `simple_rnd.f90`, `simple_test_utils.f90`,
  `simple_commanders_test_class.f90` and `production/CMakeLists.txt`).

## Out of scope

Workers of one job draw the same sequence today (each starts at the same
base) and keep doing so; offsetting the base by the part number would change
every distributed test's draws.
