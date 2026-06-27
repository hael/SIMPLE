---
name: simple-main-star
description: Use when working in SIMPLE's src/main/star subsystem, including STAR file parsing/wrapping, Relion interoperability, starproject import/export utilities, and stream-aware STAR project handling.
---

# SIMPLE `src/main/star`

This folder owns STAR and Relion interoperability.

## Read First

- `simple_starfile.f90`
- `simple_starfile_wrapper.f90`
- `simple_starproject.f90`
- `simple_starproject_utils.f90`
- `simple_relion.f90`

## Connections

- Commonly reached from project import/export commanders
- Also touches stream workflows and external tool interoperability

## Working Rule

Keep external-format adaptation concerns here rather than leaking STAR/Relion assumptions across the rest of the codebase.
