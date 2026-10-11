# The random seed through `parameters`, not the environment

Plan, 2026-10-11, for review.

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

## Files

| file | change |
|---|---|
| `src/utils/math/simple_rnd.f90` | module base and `set_rnd_seed`; `simple_seed_base` reads the module, not the environment; comments |
| `src/main/params/simple_parameters.f90`, `simple_parameters_parse.f90`, `simple_parameters_phases.f90` | the `rnd_seed` key (the args list is generated from the declaration); seeding after parsing |
| `src/main/simple_cmdline.f90` | `gen_job_descr` adds `rnd_seed` when a base is set |
| `src/main/simple_cmdline_tester.f90` | `gen_job_descr` carries `rnd_seed` from the base and keeps a given one |
| `src/utils/simple_test_utils.f90`, `production/simple_test_exec.f90` | `fix_test_seed` and `propagate` through `set_rnd_seed`; no `setenv` |
| `doc/policies/test_environment_policy.md`, the stream, SINGLE and phase-5 handover notes, `production/CMakeLists.txt` comment | the seed described as above |

## Checks

- `scripts/check_test_registry.py` and `scripts/check_descr.py`.
- Build, then `ctest -L fast`, `ctest -L library` and `ctest -L highlevel`
  on the Mac: the reports and `metrics.tsv` values equal those of the build
  before the change.
- `grep -rn SIMPLE_SEED src production scripts .github` finds nothing.

## Out of scope

Workers of one job draw the same sequence today (each starts at the same
base) and keep doing so; offsetting the base by the part number would change
every distributed test's draws.
