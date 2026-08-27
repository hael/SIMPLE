# Online Streaming Pipeline

## Problem

Process an acquisition while movies are still arriving, yet preserve the same
scientific stage boundaries used by batch processing. Streaming is an artifact
flow with restartable stages, not one continuously mutable global estimator.

## Stage graph

The GUI master launches six long-lived stage workers and exchanges framed
progress/control messages with them. Scientific data move between stages as
completed project artifacts watched on disk:

```text
movies
  -> P01 preprocessing
  -> P02 optics assignment
  -> P03 initial 2D analysis
  -> P04 reference picking and extraction
  -> P05 particle sieving
  -> P06 global pooled 2D classification.
```

Each watcher records history so a restart does not import the same completed
project twice. A stage publishes a project only after its local work is
complete; downstream stages never consume half-written in-memory state.

## P01: preprocessing

Incoming movies are grouped into bounded sets. Each set runs the production
motion-correction, CTF-estimation, and segmentation-picking algorithms. The
stage publishes corrected micrographs, CTF/quality metadata, initial
coordinates, thumbnails, histograms, and time trends. Acquisition metadata and
fitting policy remain explicit inputs to the same batch-capable commands.

## P02: optics assignment

Completed preprocessing projects are imported into a growing micrograph
project. Micrographs are assigned to optics groups, current optics and
micrograph STAR products are exported, and a bounded rolling set of optics maps
is retained. This stage changes grouping metadata, not corrected pixels.

## P03: initial analysis

Segmentation-picked particles feed a first 2D classification. Its selected
class averages provide two products: early quality feedback and the reference
bank required for stronger reference-based picking. This is the bootstrap that
turns a reference-free opening into a reference-guided stream.

## P04: reference picking and extraction

The stage watches for eligible micrographs and current picking references,
runs the exhaustive Pearson reference picker, extracts particle boxes, and
publishes per-set particle projects. Reference revisions and project cadence
are explicit stage artifacts rather than hidden shared objects.

## P05: sieving

New particle projects enter a staged `ptcl_sieve` cycle:

1. collect completed chunks and reject low-quality class averages;
2. generate coarse classification chunks;
3. promote accepted results to fine chunks unless coarse-only mode is active;
4. submit pending chunks to the worker queue.

The sieve therefore bounds contamination and computational commitment before
particles enter the global pool. Final-ingestion state is signaled only after
upstream input has remained idle for the configured interval.

## P06: global pooled 2D

Accepted sieve sets are imported into a growing project and classified with
sampled stochastic-neighborhood `cluster2D`. Dynamic resolution limits and
particle thresholds let the pool evolve as data arrive. Stepwise mode imports
only enough pending sets to reach the next threshold, deferring the rest until
the pool is ready. Snapshot requests copy stable class products without
changing the live estimator.

## Scientific invariants

- Batch and stream stages call the same motion, CTF, picker, extraction, and 2D
  algorithms wherever their semantics overlap.
- Stage cadence may differ from batch execution; an artifact becomes visible
  downstream only after the producing stage completes it.
- Project/watcher history is part of estimator identity: duplicate import would
  double-count data.
- GUI messages report or alter declared stage controls; they are not a second
  scientific data path.
- Streaming class rejection and global pooling are distinct. Rejected chunk
  classes do not silently re-enter the pool.

## Implementation

- Master: `src/main/stream/simple_stream_p00_master.f90`.
- Stages: `src/main/stream/simple_stream_p01_preprocess_new.f90` through
  `src/main/stream/simple_stream_p06_pool2D_new.f90`, with
  `simple_stream_p03_initial_analysis.f90` at stage 3.
- Watchers/state: `src/main/stream/simple_stream_watcher.f90` and
  `src/main/stream/simple_stream_state.f90`.
- Chunk/pool helpers: `src/main/stream/simple_stream_chunk*.f90` and
  `src/main/stream/simple_stream_pool2D_utils.f90`.

