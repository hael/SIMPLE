# Particle Cache Policy

## 1. Purpose and Scope

This document defines the policy for the downscaled particle disk cache
(`cache=yes`, `cache_dir=<dir>`), implemented in
`src/main/strategies/search/simple_ptcl_cache.f90` and consumed by the 2D
matcher workflows only. The 3D workflows (`refine3D` family, `solve3D`,
probabilistic 3D preparation, matcher reconstruction) do not consume the cache
and reject `cache=yes` with a hard error; caching never worked reliably for 3D
and was removed from those paths.

The samplers draw a fresh ~nsample subset of the particles every iteration and
each selected particle is read at full box and Fourier-cropped to `box_crop`.
In probabilistic mode the same particle is read twice per iteration:
probabilistic scoring (`prob_tab2D`) and search, whose batch also feeds
class-average restoration. The cache writes
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

- 2D alignment (`prob_tab2D`, `refine2D_exec` search): exact substitution.
- 2D class-average restoration (`cavger_update_sums` with
  `cropped_ptcls=.true.`): deliberate numerics change — edge taper and
  gridding source grid live at `box_crop`, no second noise normalization,
  shift and CTF pixel size scaled to the cropped grid.
Not consumers, by design: every 3D workflow (alignment, probabilistic
preparation, and reconstruction all read the full-size originals; the
`refine3D` execution strategies throw on `cache=yes`), the one-off
starting-volume `reconstruct3D`, `volassemble` (no particle reads), streaming
2D variants (long-lived processes over changing particle sets would thrash
the validity fingerprint), and the flex / offload reconstruction paths.
`flex_pca` keeps a cache of its own that shares this module's machinery
(Section 9).

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
`cache_cleanup_glob` hook in `simple_defs` (called from `simple_error`
and the tails of `simple_exec`/`single_exec`). Workers never take ownership,
so a dying worker cannot delete the cache under the other ranks or a
resubmitted part. Deletion is key-file-first, so a partially completed
cleanup can never leave a cache that still validates.

`solve2D` holds `box_crop` fixed across its stages, so the stable key lets
the owner fast path reuse one cache for every stage; there is no per-stage
invalidation.

Per-iteration `prob_align2D` calls hit an ownership fast path in
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

Uniformity is mandatory, not best-effort: restoring class averages from
cropped particles is deliberately not the same preprocessing as from
full-size ones, so ranks must never mix modes.
`ptcl_cache_assert_ready` hard-stops a worker that expected a cache and
cannot find one (the classic case: node-local `cache_dir` not visible to
every rank).

## 7. Eligibility

The cache is refused, uniformly and at every decision level (`in_use`,
`assert_ready`, `ensure`), when:

- `box_crop >= box` (nothing to gain);
- `oritype` is not `ptcl2D`/`ptcl3D` — cls3D "particles" are class averages
  in `os_out`, which the stack fingerprint cannot see.

## 8. Known Limitations

- **Signals and Fortran runtime errors bypass the cleanup hook.** SIGKILL,
  `scancel`, node failure, or a runtime abort leaves the files behind.
  Correctness is preserved (a rerun validates or rebuilds), and a rerun in
  the same execution directory reclaims the orphan, but a run in a fresh
  directory strands the old files in `cache_dir` until manually removed.
- **Node-local `cache_dir` on multi-node jobs.** The master builds on its own
  node; workers on other nodes then hard-stop with the uniformity error.
  Intentional, but there is no replicate-to-node-local mechanism.
- **Mid-run activation of previously inactive particles** (after the
  ownership fast path has skipped revalidation) fails loudly at read time
  rather than being served silently; this does not occur in the current
  workflows.
- The cached class-average restoration path is accepted on statistical
  parity with uncached execution, not bit equality; the alignment paths are
  bit-exact.

## 9. Shared Machinery

`simple_ptcl_cache` exports three helpers for other disk caches:

- `ptcl_cache_dir`: the location rule of Section 2;
- `ptcl_cache_run_token`: the hash of the execution directory used in the
  cache names;
- `ptcl_cache_space_ok`: the budget of Section 6, at most a quarter of the free
  space at the destination; an unknown free space does not block the build.

Their only other user is the `flex_pca` plane cache
(`src/main/flex/fit/simple_flex_pca_plane_cache.f90`). It holds each selected
particle's prepared, padded Fourier transform on the `box_crop` grid, and is
used with `cache=yes` when `box_crop < box`.

Both caches follow the same rules:

- they live where `ptcl_cache_dir` says, and the directory is created when
  missing;
- a stale cache is deleted before the budget check, so the check measures the
  space the rebuild can use;
- over budget, the run reads its particles without the cache;
- the cache is built under a temporary name that contains the run token and
  `_part`, and its commit record is written last. The particle cache renames
  its stack and index into place and then writes the key. The plane cache
  writes its header record into the finished file and then renames the file. A
  reader therefore finds either no cache or a complete one.

The plane cache differs in its lifetime and validation. Its published name
carries `box_crop` and a hash of the absolute project path, not the run token.
The run that built it deletes it once the master has its embedding, because
nothing after the embedding reads planes (and with `mkdir=yes` every run works
on its own project copy, so no later run could adopt it). The header records the project path and modification
time, the geometry and a hash of the master's particle selection. The master
accepts only the exact selection. A worker accepts the cache when its partition
lies within it, and otherwise reads its particles natively.

## 10. Future Extensions

- **Orphan scavenger**: sweep stale `ptcl_cache_*` file sets in `cache_dir`
  at `ptcl_cache_ensure` time (age- or dead-key-based) to reclaim
  signal-killed leftovers automatically.
- **Benefit predicate**: `ensure` knows `nsel`; with the planned iterations,
  sampling schedule, and crop ratio it can compute predicted bytes saved
  versus the one-time build, then decline uneconomical caching through the
  existing uniform fallback.
- **UI promotion**: `cache`/`cache_dir` are `UI_VIS_DEVELOPER`; promote after
  validation, and consider defaulting `cache=yes` for `solve2D`.
