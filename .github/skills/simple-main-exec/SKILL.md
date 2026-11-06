---
name: simple-main-exec
description: Use when working in SIMPLE's src/main/exec subsystem, including the execution-router modules (exec_<topic>_commander) that dispatch program names to commander objects for simple_exec, single_exec and simple_test_exec, the shared execution helpers (restarted, scripted and asynchronous execution, the job description written to the project, master thread counts), and how they connect to the private worker dispatch of simple_private_exec and the simple_stream dispatch.
---

# `src/main/exec`: the routers between the front doors and the commanders

Owns the routers `exec_<topic>_commander(which, cline, l_silent, l_did_execute)` of `simple_exec`,
`single_exec` and `simple_test_exec`, and `simple_exec_helpers.f90`. A router dispatches a program
name to a commander and nothing else (layering policy section 2). Not here: registration and keys
(`ui/`, `params/`), the workers (`dispatch_private_prg` in `production/simple_private_exec_driver.f90`)
and the stream programs (the `select case` in `production/simple_stream.f90`).

## Rules

Silent ones first: nothing reports these.

- The front door (`production/<exe>.f90`) calls every router in turn; the first case that matches
  runs, and the later routers return at once. A name in two routers runs only in the first, and
  nothing checks the production routers for that. A router's topic need not be the UI category
  (`cls_expansion` is a refine2D program routed by `simple_exec_denoise`).
- After the commander returns, the front door appends the job description to the `jobproc` segment
  of the project `cline` names (`update_job_descriptions_in_project`: reads and writes the non-data
  segments), except under the test and private executables and for program names containing
  `print`, `info` or `report`. A naming rule: a new program that must leave no job description has
  one of these in its name, and only such a program.
- Beyond dispatch (deviations from policy section 2, listed in its section 8: do not copy):
  `single_exec` sets `mkdir=yes` when it is not given, and some SINGLE cases set `mkdir` over what
  the user gave (`atoms_register`, `crys_score`, `detect_atomic_species`, `simulate_nanoparticle`,
  `trajectory_graphene_subtr`, `trajectory_denoise`: no; `atoms_stats`: yes).
  `nrestarts=` on the cases that test it (`solve2D`, `solve3D`, `solve3D_cavgs`, the `refine3D`
  family, `classify3D_refs`, `autorefine3D_nano`) runs the program that many times as child
  processes (`restarted_exec`); a child's non-zero exit is only reported.
- `script=yes` builds a `parameters` with `mkdir=no`, writes a queue script and exits before routing
  (`script_exec`). `solve3D_cavgs` with `nrestarts_collapse` routes to a second commander, and
  `solve3D` with `screen` to `exec_screen` (unreachable: `screen` is no accepted key); the front
  doors rewrite the project's job descriptions after the commander (policy section 8).
- `simple_exec_helpers.f90` is library code: commanders get `set_shmem_flag` (deletes `nparts` when
  it is 1) and `set_master_num_threads` (sets the OpenMP team from `SLURM_CPUS_PER_TASK`, else at
  most `NTHR_SHMEM_MAX`) through `simple_commanders_api`, and strategies call them too. A module-level
  import added there reaches every importer of that umbrella (skill `simple-main-apis`).
- Loud: a name no router matches stops after the chain ("not recognized"); a case for a name the UI
  does not register never runs, because `cmdline%parse` stops first: add no alias case (the stream's
  `initial_analysis` and `opening2D` cases were dropped 2026-11-05).
  `l_silent` only drops the git version and timer lines.
- Test routers (`simple_test_exec_class|highlevel|single|stream.f90`): each test UI program has
  exactly one case, checked by `scripts/check_test_registry.py` in the fast gate.

## Wiring

A new program is a case in the router of its executable (its UI `executable`). A new router module
goes into that executable's API module (`src/main/apis/<exe>_api.f90`) and its call chain in
`production/<exe>.f90`. A worker is a case in `dispatch_private_prg`, its commander in
`simple_private_exec_api` and its keys in `ui/simple_private_prgs.f90` unless it has a public name.
The whole chain: `.github/skills/simple-main-ui/references/cli_change_chain.md`.

## Tests

- No unit suite covers the routers or `simple_exec_helpers`. The test routers are checked by
  `scripts/check_test_registry.py` (`scripts/run_fast_gate.sh build` runs it).
- A routing change: run the program, or the high-level entry that runs it
  (`doc/policies/test_environment_policy.md` section 1.2). `simple_exec prg=<name>` alone prints the
  usage from the UI and stops before routing, so it shows registration, not the route.

## See also

`doc/policies/layering_and_project_ownership_policy.md` (section 2), `doc/policies/ui_layer_policy.md`
(private programs); skills `simple-main-ui`, `simple-main-apis`, `simple-main-commanders`.
