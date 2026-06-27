---
name: simple-main-class
description: Use when working in SIMPLE's src/main/class subsystem, including class averaging, class FRCs, restore paths, and 2D class-product generation used by cluster2D and related workflows.
---

# SIMPLE `src/main/class`

This folder owns class-average products and related statistics.

## Read First

- `simple_classaverager.f90`
- `simple_classaverager_core.f90`
- `simple_classaverager_restore.f90`
- `simple_class_frcs.f90`

## Connections

- Called heavily from `cluster2D`, `mkcavgs`, and related commanders
- Depends on `image/`, `project/`, and `ori/`
- Often participates in restore/restart flows, not just fresh averaging

## Working Rule

When a bug looks like “class products are wrong,” check both the averaging core and the project/orientation state feeding it.
