---
name: simple-main-star
description: Use when working in SIMPLE's src/main/star subsystem, including STAR file reading and writing, the Fortran wrappers of the C++ starfile library, Relion interoperability, starproject import/export utilities, and the stream's STAR export.
---

# SIMPLE `src/main/star`

This folder owns STAR and Relion interoperability.

## Read First

- `simple_starfile.f90`: optics, micrograph and 2D-particle tables written for RELION-compatible STAR files
- `simple_starfile_wrapper.f90`: module `simple_starfile_wrappers`, Fortran bindings to the C++ library in `src/extlibs/starfile`
- `simple_starproject.f90`, `simple_starproject_utils.f90`: STAR project import/export
- `simple_starproject_stream.f90`: STAR export of the stream's stages and 2D pool
- `simple_relion.f90`

## Connections

- Reached from `commanders/simple/simple_commanders_starproject.f90`, `simple_commanders_relion.f90` and project export (`simple_sp_project`), from the refine2D class-average output (`cls2D` star file), from the stream stages and 2D pool, and through the wrappers from `motion/simple_motion_model` and `ctf/simple_ctf_estimate_fit`
- Tests: `unit_project` ("STAR file", "STAR project", "STAR stream export")

## Working Rule

Keep external-format adaptation concerns here rather than leaking STAR/Relion assumptions across the rest of the codebase.
