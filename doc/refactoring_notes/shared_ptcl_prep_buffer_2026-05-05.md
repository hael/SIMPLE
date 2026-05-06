# Shared particle-prep buffer for `prepimg4align` / `prep_imgs4rec`

Date: 2026-05-05
Status: implementation note; Tier 0 default, Tier 1 stats-passing path is the default contract
Scope: `refine3D` and the shared 2D matcher/classaverager path

## Context

In the current `refine3D_exec` batch loop, when `do_write_partial_recs` is on, every batch
particle is touched by two distinct preprocessing pipelines:

1. **`prepimg4align`** (in [`simple_matcher_2Dprep.f90`](../../src/main/strategies/search/simple_matcher_2Dprep.f90))
   called from `build_batch_particles3D` / `build_batch_particles2D` in
   [`simple_matcher_ptcl_batch.f90`](../../src/main/strategies/search/simple_matcher_ptcl_batch.f90).
2. **`prep_imgs4rec`** (in [`simple_matcher_3Drec.f90`](../../src/main/strategies/search/simple_matcher_3Drec.f90))
   called from `maybe_restore_batch` in
   [`simple_strategy3D_matcher.f90`](../../src/main/strategies/search/simple_strategy3D_matcher.f90)
   *after* the search has run for the same batch.

The 2D matcher has a related pattern: `build_batch_particles2D` prepares particles for
alignment, keeps a full-box copy in `ptcl_imgs`, and `cavger_update_sums` later prepares
that copy for Cartesian class averaging.

To survive the search step in 3D, the real-space full-box image is duplicated into a
batch-sized scratch array `ptcl_rec_imgs` (`copy_fast` in `build_batch_particles3D`).
In 2D, the corresponding retained full-box copies are `ptcl_imgs`.

## Verified Normalization Order

Before Tier 0, the reconstruction and class-average prep route was not using the intended
order.

`image%norm_noise_taper_edge_pad_fft` previously:

1. tapers edges in-place on `self`;
2. computes the mask-based noise mean and variance on the already tapered `self`;
3. applies that normalization while copying into the padded output;
4. FFTs the padded output.

This was corrected in Tier 0. The current route computes the mask-based noise stats on
the untapered full-box particle and then passes those stats into the taper/pad/FFT path.
This matches the alignment path contract, where `image%norm_noise_fft_clip_shift`
computes mask-based noise stats before the FFT-shift/normalization pass.

The optimization plan below therefore assumes the intended route is already in place:
compute noise stats on the untapered full-box particle, then perform any tapering,
padding, FFT, or downstream CTF/gridding work needed by each consumer.

## Per-Particle Comparison After Correcting Normalization Order

| Step | `prepimg4align` (alignment) | `prep_imgs4rec` / `cavger_update_sums` |
| --- | --- | --- |
| Input | full-box real-space `img` | full-box real-space image copy |
| Noise normalization | full-box mask reduction over `build%lmsk` before FFT | same intended full-box mask reduction over `build%lmsk` before taper/pad/FFT |
| Real-space envelope | alignment mask after clipped FFT/IFFT path | edge taper before padded FFT |
| FFT | full-box image, then clip | zero-padded `boxpd` image |
| Crop in FT | clip to `box_crop` | none |
| Shift | applied in FT (`shvec = -(x,y) * crop_factor`) | passed to `gen_fplane4rec` |
| CTF | phase-flip only when `CTFFLAG_YES` | full CTF / CTF^2 (+ sigma2 when `l_ml_reg`) via `gen_fplane4rec` |
| Output | `img_out` (FT, `box_crop`), `img_out_pd` (FT, `box_croppd`) | `fplane_type` or class-average splat input |

The two pipelines still run on **different Fourier lattices** (`box_crop` vs `boxpd`),
apply **different real-space envelopes** after the shared normalization point, and use
**different CTF treatments**. The FFTs themselves are still not directly reusable.

The shared work, once the normalization order is fixed, is:

* the `norm_noise_*` step (full-box mask reduction over `build%lmsk`, identical mean / sigma);
* the extra pass over the same particle in a different OMP loop, with colder caches;
* the full-box real-space copy that exists only to bridge the search and restore/class-average steps.

## Hard Parallelization Constraint

Do **not** move `grid_ptcl` / `update_rec` into the particle-prep OpenMP loop.

`update_rec` currently prepares f-planes in parallel but grids them serially into
`build%eorecvols`. That serial gridding order must be preserved unless the reconstructor
is explicitly redesigned around thread-local reconstruction accumulators or another
proven race-free mechanism.

This note therefore keeps the current parallelization model:

* particle image prep and f-plane generation may stay parallel;
* reconstruction gridding remains serial;
* class-average accumulation keeps its existing synchronization/loop structure.

## Suggested Implementation Tiers

### Tier 0 — Restore Intended Normalization Order

`image%norm_noise_taper_edge_pad_fft` now bases noise normalization on the untapered
full-box particle, matching the intended alignment/reconstruction contract.

Validation should compare old/new outputs and timings for:

* `refine3D` with partial reconstruction on;
* the shared 2D matcher path with class averaging on;
* at least one CTF-on case;
* an `l_ml_reg` case if available.

This was the correctness prerequisite for sharing normalization work. From the theoretical
standpoint, normalization/tapering order no longer blocks the shared-prep design.

### Tier 1 — Hoist the Noise Normalization

The stats-passing route is now the default contract. The shared carrier is
`type(noise_stats)` in `simple_type_defs`, and the image prep methods consume precomputed
stats directly.

Current API shape:

* compute `(mean, invstd)` once from the full-box image and `build%lmsk` via
  `image%calc_noise_stats`;
* pass those stats into:
  * `image%norm_noise_fft_clip_shift`
  * `image%norm_noise_fft_clip_shift_ctf_flip`
  * `image%norm_noise_taper_edge_pad_fft`

The stats-passing route is safer because the current alignment prep mutates
`build%imgbatch(iptcl_batch)`, while reconstruction/class averaging needs an independent
full-box real-space copy.

Batch integration:

* `build_batch_particles2D` can store per-particle stats for `cavger_update_sums`;
* `build_batch_particles3D` can store per-particle stats for the post-search
  `prep_imgs4rec` route;
* `ptcl_noise_stats` is a required argument for the batch/restore interfaces;
* callers without shared overlap still compute stats explicitly before image prep.

Expected wins:

* saves one full-box mask reduction per particle when both paths need the same particle;
* applies to `refine3D`;
* applies to the shared 2D matcher/classaverager path.

### Tier 2 — Generate Reconstruction F-Planes After Search, Not During Batch Build

Do not fuse reconstruction prep into `build_batch_particles3D` before the search. The
search can update particle orientation, shift, state, weight, and sigma2. Reconstruction
prep must consume the post-search metadata.

A safer experiment is to keep the restore step after the search, but reduce the amount of
state carried into it:

* keep `maybe_restore_batch` after the search;
* keep `prep_imgs4rec` parallel;
* keep `update_rec` serial;
* pass precomputed normalization stats or already-normalized full-box copies into
  `prep_imgs4rec`;
* benchmark whether this removes enough repeated work to justify the extra interface.

This preserves the current metadata timing and the current gridding parallelization.

### Tier 3 — Reduce the Per-Batch Real-Space Copy

`build_batch_particles3D` currently does
`imgs4rec(iptcl_batch)%copy_fast(build%imgbatch(iptcl_batch))` because `prepimg4align`
overwrites `imgbatch(i)`.

Once Tier 1 is validated, investigate whether the retained copy can be slimmer:

* normalized full-box copy rather than raw full-box copy;
* per-batch copy only when restoration/class averaging is enabled;
* thread scratch only if the restore step can still run after search without losing the
  needed post-search metadata.

Dropping the copy entirely is not part of the first safe implementation because restore
still runs after search.

## What Does Not Help

* Sharing the FFT itself — the two pipelines target different lattices (`box_crop` vs
  `boxpd`), and the reconstruction/class-average lattice is preceded by its own real-space
  envelope and padding. Skip.
* Sharing the shift — different conventions; alignment bakes it into the cropped FT,
  reconstruction passes it to `gen_fplane4rec`. Skip.
* Sharing the CTF — phase-flip versus full CTF / CTF^2 / sigma2. Skip.
* Parallelizing `grid_ptcl` inside the particle-prep loop. Skip unless reconstruction
  accumulation is redesigned to be race-free.

## 2D Applicability

`build_batch_particles2D` prepares particles for alignment and keeps full-box copies in
`ptcl_imgs`. `cavger_update_sums` later calls `norm_noise_taper_edge_pad_fft` on those
copies before class-average splatting.

Therefore:

* Tier 0 is now the default 2D class averaging behavior;
* Tier 1 applies directly to 2D when alignment and class averaging are both active;
* Tiers involving `prep_imgs4rec`, `update_rec`, and `grid_ptcl` are 3D-only.

## Recommended Order

1. Compare current timings for 2D benchmarks and a representative `refine3D` run with
   Tier 0 as default behavior.
2. Benchmark the default **Tier 1** stats-passing route for 2D and representative
   `refine3D` runs.
3. For 3D, keep reconstruction prep after search and keep gridding serial; experiment only
   with sharing normalization stats or normalized full-box copies.
4. Keep `prep_imgs4rec` as a public entry point for `calc_3Drec` (the offline
   reconstruction path that does not run alignment) because that caller has no overlap to
   exploit.
