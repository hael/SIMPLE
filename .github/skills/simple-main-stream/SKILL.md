---
name: simple-main-stream
description: Use when working in SIMPLE's src/main/stream subsystem, including online/mini-stream orchestration, chunk and microchunk processing, staged pipeline steps, watchers, queue/message definitions, and streaming variants of preprocessing and cluster2D workflows.
---

# SIMPLE `src/main/stream`

This folder owns the streaming pipeline.

## Read First

- `simple_stream_p00_master.f90`
- `simple_stream_p01_preprocess.f90`
- `simple_stream_p02_assign_optics.f90`
- `simple_stream_p03_opening2D.f90`
- `simple_stream_p04_refpick_extract.f90`
- `simple_stream_p05_sieve_cavgs.f90`
- `simple_stream_p06_pool2D.f90`

## Structure

- Pipeline stages are explicit and numbered
- Many `_new` variants exist; compare carefully before changing shared behavior
- Chunking, watchers, microchunking, and cluster2D-specific helpers are important supporting layers

## Working Rule

Do not assume batch semantics apply unchanged to stream code; stage boundaries and artifact cadence matter here.
