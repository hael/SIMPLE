---
name: simple-main-exec
description: Use when working in SIMPLE's src/main/exec subsystem, including the execution-router modules (exec_<topic>_commander) that dispatch program names to commander objects for simple_exec, single_exec and simple_test_exec, and how they connect to the private worker dispatch of simple_private_exec.
---

# SIMPLE `src/main/exec`

`exec/` is the routing layer between the executable front door and commanders.

## Read First

- `simple_exec_project.f90`
- `simple_exec_preproc.f90`
- `simple_exec_refine2D.f90`
- `simple_exec_refine3D.f90`
- `simple_test_exec_*.f90`

## Role

- Groups related programs by topic
- Holds `select case` dispatch to commander instances
- Keeps `production/*.f90` entrypoints simple

## Shape

- A router is `exec_<topic>_commander(which, cline, l_silent, l_did_execute)`:
  it returns at once when `l_did_execute` is already set, otherwise sets it and
  dispatches, resetting it in `case default`.
- `production/simple_exec.f90` (likewise `single_exec.f90`, `simple_test_exec.f90`)
  calls every router in turn and stops if none executed. A new router module goes
  into the executable's API module (`src/main/apis/*_exec_api.f90`) and that call
  chain.
- Workers and internal steps are not routed here: `simple_private_exec` dispatches
  them in `dispatch_private_prg` of `production/simple_private_exec_driver.f90`
  (`prg=refine2D`/`refine3D` are the partition workers, `refine2D_distr` the 2D
  master).
- Test routers (`simple_test_exec_class|highlevel|single|stream.f90`): each
  program of the test UI has exactly one case; `scripts/check_test_registry.py`
  fails the build otherwise.

## Working Rule

If a new UI program exists but never runs, the missing wire is often here.
