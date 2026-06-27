---
name: simple-main-exec
description: Use when working in SIMPLE's src/main/exec subsystem, including the execution-router modules that dispatch program names to commander objects for SIMPLE, SINGLE, and test workflows.
---

# SIMPLE `src/main/exec`

`exec/` is the routing layer between the executable front door and commanders.

## Read First

- `simple_exec_project.f90`
- `simple_exec_preproc.f90`
- `simple_exec_cluster2D.f90`
- `simple_exec_refine3D.f90`
- `simple_test_exec_*.f90`

## Role

- Groups related programs by topic
- Holds `select case` dispatch to commander instances
- Keeps `production/*.f90` entrypoints simple

## Working Rule

If a new UI program exists but never runs, the missing wire is often here.
