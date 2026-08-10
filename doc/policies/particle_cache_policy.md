# Particle Cache Policy

## 1. Purpose and Scope

This document defines the policy for the downscaled particle disk cache
(`cache=yes`, `cache_dir=<dir>`), implemented in
`src/main/strategies/search/simple_ptcl_cache.f90` and consumed by the 2D and
3D matcher workflows.

The samplers draw a fresh ~nsample subset of the particles every iteration and
each selected particle is read at full box and Fourier-cropped to `box_crop`.
In probabilistic mode the same particle is read up to three times per
iteration: probabilistic scoring, search, and reconstruction. The cache writes
the iteration-independent prefix of that preprocessing to disk once and serves
all subsequent reads from it, trading disk space at `cache_dir` for a roughly
`(box/box_crop)^2` reduction in read volume.

Because the samplers sweep near-disjoint subsets before revisiting anything, a
partial cache would never hit; the cache always covers all active particles.

## 2. What Is Cached

An entry is the noise-normalized (against the full-box `lmsk`),
Fourier-cropped particle, stored as a real-space `box_crop` image. Fourier
crop -> inverse FFT -> forward FFT is an exact round trip, so reading an entry
back reproduces the same coefficients `prepimg4align` would have computed.
The shift is not cacheable (re-read from the project every iteration) and the
CTF is deliberately left out so the cache does not depend on CTF parameters.

Cache files (stack, index, key) live in `cache_dir` (else
`SIMPLE_PTCL_CACHE_DIR`, else the execution directory). The basename carries
the project name, `box_crop`, and a hash of the absolute execution directory,
so concurrent runs sharing a fast disk cannot collide and every distributed
rank recomputes the same name locally (qsys workers cd into the master's
directory).

## 3. Consumers

- 2D alignment (`prob_tab2D`, `cluster2D_exec` search): exact substitution.
- 2D class-average restoration (`cavger_update_sums` with
  `cropped_ptcls=.true.`): deliberate numerics change — edge taper and
  gridding source grid live at `box_crop`, no second noise normalization,
  shift and CTF pixel size scaled to the cropped grid.
- 3D alignment (`refine3D` search, `prob_tab`, `prob_tab_neigh` via
  `build_batch_particles3D`): exact substitution through
  `prepimg4align_cached`.
- 3D matcher reconstruction (`calc_3Drec`, `calc_projdir3Drec` via
  `prep_imgs4rec(cached=.true.)`): mirrors the 2D restoration treatment
  exactly. Sound because `box*smpd == box_crop*smpd_crop` — the padded grids
  cover the same physical extent, so a given Fourier index denotes the same
  spatial frequency on either grid, making CTF and shift handling exact; the
  only numerics delta versus uncached execution is the edge taper acting at
  `box_crop`.

Not consumers, by design: the one-off starting-volume `reconstruct3D`
(may run before the cache exists; reads each particle once), `volassemble`
(no particle reads), streaming 2D variants (long-lived processes over
changing particle sets would thrash the validity fingerprint), and the flex /
offload reconstruction paths. Cached reconstruction is opt-in by explicit
argument at every level (`calc_3Drec`/`calc_projdir3Drec`, `init_rec`,
`prep_imgs4rec` all take `cached`, defaulting to full-size behavior), and
only the refine3D matcher passes a value — one that went through
`ptcl_cache_assert_ready`, so an opportunistically visible leftover cache can
never be picked up by standalone reconstruction or mix cached and uncached
partials across ranks.

## 4. Validity Contract

The key file is the commit record; no reader accepts cache files without a
matching key. It records:

- the geometry line: `box`, `box_crop`, `smpd`, `smpd_crop`, `msk`,
  `oritype`, particle and stack counts, and the execution directory verbatim
  (so a name-hash collision from a different directory can never validate);
- the source fingerprint: every stack in project order with its particle
  range, physical size and mtime, plus every particle's `stkind`/`indstk`
  mapping (the range-derived fallback is covered by the per-stack
  `fromp`/`top`). A reordered or remapped project over unchanged stacks
  therefore cannot validate a stale cache.

The check is conservative: a touched file forces a rebuild. Only the rank
that decides whether to rebuild pays the full fingerprint; consumers check
the geometry line and inherit the verdict.

## 5. Lifecycle and Ownership

The process that builds — or adopts a valid leftover, e.g. after a killed
predecessor in the same execution directory — owns the cache files and
removes them on normal exit and on hard exception, via the
`cache_cleanup_glob` hook in `simple_defs` (called from `simple_exception`
and the tails of `simple_exec`/`single_exec`). Workers never take ownership,
so a dying worker cannot delete the cache under the other ranks or a
resubmitted part. Deletion is key-file-first, so a partially completed
cleanup can never leave a cache that still validates.

Cache-enabled `abinitio3D` uses the final active downscaling-ladder
`box_crop` for every refine3D stage. The stable key lets the owner fast path
reuse one cache throughout all eligible stages. Stage low-pass limits still
follow the ladder, while `cache=no` retains the stage-specific crop schedule.
The user's `cache=yes` request is re-stamped onto every stage command line, so
a stage-local fallback does not permanently disable later stages. A stage that
declines the cache releases its files; a later eligible stage may rebuild the
same final-crop cache.

Per-iteration `prob_align` calls hit an ownership fast path in
`ptcl_cache_ensure` (same owned key name, key file exists) and skip the full
revalidation.

## 6. Space Budget and Uniform Fallback

Before building, the exact predicted size (`nsel * box_crop^2 * 4` bytes plus
headers) is checked against the free space at the destination (statvfs); the
cache may claim at most 25% of it. Over budget — or with no active
particles, a denoised primary source, or a non-particle oritype — the whole
run falls back to uncached execution *uniformly*: `disable_cache` flips
`cache=no` on both params and the command line before any worker command
line is generated.

Uniformity is mandatory, not best-effort: restoring class averages or
reconstructing from cropped particles is deliberately not the same
preprocessing as from full-size ones, so ranks must never mix modes.
`ptcl_cache_assert_ready` hard-stops a worker that expected a cache and
cannot find one (the classic case: node-local `cache_dir` not visible to
every rank).

## 7. Eligibility

The cache is refused, uniformly and at every decision level (`in_use`,
`assert_ready`, `ensure`), when:

- `box_crop >= box` (nothing to gain);
- the primary particle source is the denoised stack (`ptcl_src=den`) — the
  entries derive from the raw stacks and would be the wrong pixels. A
  denoised *objective* (`objfun_den=yes`) is eligible: those images are read
  separately into a matcher-lifetime full-size workspace;
- `oritype` is not `ptcl2D`/`ptcl3D` — cls3D "particles" are class averages
  in `os_out`, which the stack fingerprint cannot see.

## 8. Known Limitations

- **Signals and Fortran runtime errors bypass the cleanup hook.** SIGKILL,
  `scancel`, node failure, or a runtime abort leaves the files behind.
  Correctness is preserved (a rerun validates or rebuilds), and a rerun in
  the same execution directory reclaims the orphan, but a run in a fresh
  directory strands the old files in `cache_dir` until manually removed.
- **Fixed-crop cost in cache-enabled `abinitio3D`.** Early stages use the
  final ladder crop instead of their smaller stage-local crops. This increases
  early-stage cache size and matching work, but avoids a sequential full-size
  cache rewrite at every crop transition. Uncached runs keep stage-local crops.
- **Node-local `cache_dir` on multi-node jobs.** The master builds on its own
  node; workers on other nodes then hard-stop with the uniformity error.
  Intentional, but there is no replicate-to-node-local mechanism.
- **Mid-run activation of previously inactive particles** (after the
  ownership fast path has skipped revalidation) fails loudly at read time
  rather than being served silently; this does not occur in the current
  workflows.
- The cached reconstruction and restoration paths are accepted on
  statistical parity (FSC trajectory, final maps) with uncached execution,
  not bit equality; the alignment paths are bit-exact.

## 9. Future Extensions

- **Orphan scavenger**: sweep stale `ptcl_cache_*` file sets in `cache_dir`
  at `ptcl_cache_ensure` time (age- or dead-key-based) to reclaim
  signal-killed leftovers automatically.
- **Benefit predicate**: `ensure` knows `nsel`; with the planned iterations,
  sampling schedule, and fixed crop ratio it can compute predicted bytes saved
  versus the one-time build and larger early-stage matching cost, then decline
  uneconomical caching through the existing uniform fallback.
- **Stage-local reads from the final-crop cache**: Fourier cropping is nested
  (normalization precedes any crop), so readers could truncate final-crop
  records to each stage's planned crop. This would preserve one cache build
  while recovering the smaller early-stage matching grids.
- **UI promotion**: `cache`/`cache_dir` are `UI_VIS_DEVELOPER`; promote after
  validation, and consider defaulting `cache=yes` for `abinitio3D`.
- **Denoised-source entries**: per-source cache entries would lift the
  `ptcl_src=den` exclusion if that path becomes I/O-bound in practice.
