---
name: simple-main-project
description: Use when working in SIMPLE's src/main/project subsystem, including the sp_project type, segment ownership for mic/stk/ptcl/cls/out/optics data, binary project I/O and segment writes, merging distributed alignment documents, the out-segment artifact registry, reading projects of earlier releases, and project-level helpers (even/odd partition, state weight sets, particle subsets, solve3D_addon supersets).
---

# SIMPLE `src/main/project`

This is the project-state backbone of SIMPLE.

## Read First

- `simple_sp_project.f90`
- `simple_sp_project_core.f90`
- `simple_sp_project_io.f90`
- `simple_sp_project_ptcl.f90`
- `simple_sp_project_cls.f90`
- `simple_sp_project_stk.f90`
- `simple_sp_project_out.f90`

## Structure

- `simple_sp_project.f90` declares the public type and interfaces
- Behaviors live in submodules by data segment responsibility (`_core`, `_io`, `_mic`, `_stk`, `_ptcl`, `_cls`, `_out`, `_optics`)
- Segment ownership is important: mic, stk, ptcl2D, cls2D, cls3D, ptcl3D, out, optics, projinfo, jobproc, compenv
- Companion modules: `simple_eo_partition` (project-owned even/odd halves), `simple_state_weight_set`, `simple_ptcl_subset`, `simple_mic_subset`, `simple_project_superset`, `simple_ptcl_layout`

## Persistence

- Write only the segments you changed (`write_segment_inside(oritype)`); a full
  `write` is for command-line steps and run ends. A workflow that owns its
  project in memory writes it once per iteration (the PCG master and the refine2D
  master write more: `doc/policies/layering_and_project_ownership_policy.md`
  section 8).
- `merge_algndocs(..., write=.false.)` merges distributed workers' documents
  into memory without touching the file; the default also writes the segment.
- There is no file locking: one process writes a project file at a time. Workers
  read the file and write part files, not the project. State that must agree
  across processes has one maker, the master, and verifying consumers
  (`ensure_field_halves`). Rules: the layering policy, section 5.
- `update_projinfo` names the run's own project file; writes go to the file in
  `projinfo`.
- Register run artifacts in the `out` segment through its getters and adders
  (`add_vol2os_out`, `add_fsc2os_out`, `get_vol`, ...), even when nothing reads
  them yet.
- Project files of earlier releases must still read (the only backwards
  compatibility release 4 keeps): see the earlier-release cases of the
  `eo_partition` and `state_weight_set` testers. Contract gaps are repaired
  once by `fix_projfile` (`doc/policies/project_stack_indexing_policy.md`),
  not by run-time fallbacks.

## Connections

- Used by almost every workflow
- Common neighbors: `ori/`, `fileio/`, `stream/`, `exec`, `commanders`

## Working Rule

When changing project behavior, think in terms of segment contracts and persistence, not just in-memory fields.
