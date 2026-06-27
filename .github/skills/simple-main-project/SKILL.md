---
name: simple-main-project
description: Use when working in SIMPLE's src/main/project subsystem, including the sp_project type, segment ownership for mic/stk/ptcl/cls/out/optics data, binary project I/O, and project-level metadata and artifact management.
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

## Structure

- `simple_sp_project.f90` declares the public type and interfaces
- Many behaviors live in submodules by data segment responsibility
- Segment ownership is important: mic, stk, ptcl2D, cls2D, cls3D, ptcl3D, out, optics, projinfo, jobproc, compenv

## Connections

- Used by almost every workflow
- Common neighbors: `ori/`, `fileio/`, `stream/`, `exec`, `commanders`

## Working Rule

When changing project behavior, think in terms of segment contracts and persistence, not just in-memory fields.
