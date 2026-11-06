---
name: simple-main-motion
description: Use when working in SIMPLE's src/main/motion or src/main/preprocess subsystems, including reference-free movie motion correction (global hybrid alignment, patched alignment and the 18-term polynomial, frame and dose weighting, EER fractions, defect curation), the binary motion model (.mmodel) and its STAR export, gain references and gain flipping (flip_auto analysis), the intg/forctf products, particle extraction from micrographs or motion-corrected movie frames (ptcl_extractor), power-spectrum thumbnails, and the nano time-series aligner.
---

# `src/main/motion` and `src/main/preprocess`: motion correction and particle extraction

`motion/` corrects one movie without references: `motion_correct_iter` (one movie, its products and
micrograph row), `simple_motion_correct` (frames, gain, defects, scaling, global then patch
alignment, sums), `motion_align_hybrid`, `motion_patched`, `motion_model` (the `.mmodel`),
`motion_correct_utils` (gain, EER fractions, polynomial basis), the gain-flip analysis and
`motion_align_nano` (SINGLE time series). `preprocess/` holds `ptcl_extractor` (particles from the
micrograph, or from the movie through the `.mmodel`) and `pspec_thumb_iter`. `motion/` also holds
`simple_dose_series_manifest` (the dose-series half maps of `refine_motion_model`, `dose_series.txt`;
its tests run in `unit_reconstruction suite=dose_series_manifest`). Scheduling, gain
preparation and project writes belong to the strategies in `strategies/parallelization` and the
stream's preprocessing stage; CTF to `ctf/`; box files to `pick/`.

## Rules

Silent ones first: nothing reports these.

- Apply the gain once. The scheduler flips: `flip_gain` writes `./<gain>_flip<MODE>.mrc` and
  repoints `gainref`, then the caller sets `params%flipgain` and the command line's `flipgain` to
  `no` (the strategies, the stream's worker command). Frames are multiplied by the gain as given; the
  `.mmodel` stores that file and the `flipgain` in force, and the extractor reapplies exactly that.
  A path that flips without the reset flips twice. X reverses the first index, Y the second; a
  non-EER gain must match the frames (`THROW_HARD`). The stream resolves `flip_auto` and `generate`
  before workers start.
- Motion runs on the scaled grid, `smpd/scale_movies` (after EER upsampling). `motion_correct_iso`
  overwrites the caller's `ctfparams%smpd` with it and the row gets `smpd`, `xdim`, `ydim`: the CTF
  fit of `forctf` uses that same `ctfparams`. `trs` is in movie pixels. Offsets, coefficients and
  patch geometry in the `.mmodel` are scaled-grid pixels; STAR export converts copies to movie
  pixels with frame 1 as reference. Never rescale model arrays in place; consumers use the saved
  model geometry, not the project's sampling.
- Offsets are frame displacements: frames are shifted by their negative, the extractor windows at
  `corner - offset`. `mcconvention=simple` refers them to the central frame; `first`/`relion`
  re-reference to frame 1 before the frames are shifted.
- `nframesgrp` is for time series: `motion_correct_iter%iterate` stops on it otherwise (the temporary
  grouped movie would be recorded as the movie).
- `forctf` is the plain (1/nframes) sum without dose weights, made before dose weighting touches the
  frame buffers, and transient (deleted after the CTF fit); `intg` uses the alignment frame weights
  and, with `total_dose`, dose weights. `apply_dose_weighing(..., mode='GG')` is called in both sums
  and in the extractor (`prepare_frames`): change them together.
- Patches are accepted only when both fit RMSDs are below 4 scaled pixels (5 for `first`/`relion`),
  after one retry on a halved grid; else the sums are global-only. Consumers honour
  `patch_accepted`: rejected patch arrays may still be in the model.
- `simple_motion_correct` keeps the movie in module variables and `simple_ft_expanded`'s transfer
  matrix is module-level: one movie at a time per process, never movies on parallel threads; every
  exit path runs the kill calls.
- The `.mmodel` is the only motion metadata SIMPLE reads (`mc_starfile` is for other software).
  Readers reject other `FILE_VERSION`s: a layout change bumps it and the tester's binary fixture.
- Extractor: `init_mic`; `init_mov` falls back to `init_mic` with only a warning when `mcmodel` is
  missing; `init_model` (`refine_motion_model`) stops. Coordinates are box-file corners on the
  `intg` grid; centre = corner + box/2, +1 for `pix2polycoords`; motion is referred to the model's
  `fixed_frame` (central), so `extractfrommov=yes` needs `mcconvention=simple`: motion correction records
  the convention on the micrograph row and `init_mov` stops on `first`/`relion` (which register `intg` to
  frame 1); a row of an earlier run without the key passes.
  `local_polynomial = .false.` in the extractor is an experiment to keep. Per-thread buffers are
  sized `nthr_glob`: call `extract_particles*` outside a parallel region. `post_process` negates for
  `pcontrast=black`, removes the ramp, normalises noise.

## Where to start

- `doc/policies/motion_correction_policy.md` (the convention tables), then
  `simple_motion_correct_iter.f90` (one movie; with `tseries=yes` no `forctf`, model or STAR).
- Callers: `strategies/parallelization/simple_preprocess_strategy.f90`, `simple_extract_strategy.f90`,
  `simple_refine_motion_model_strategy.f90`; `stream/stages/simple_stream_stage_preprocess.f90` (gain).

## Tests

- `build/bin/simple_test_exec test=unit_motion` (under 1 s): `suite=motion_gain`, `motion_model`,
  `particle_extractor`, `fraction_dose_metadata`. Also `test=unit_numerics suite=shift_search`,
  `test=unit_image suite=image` (flips), `test=unit_stream suite=preprocessing`.
- `test=stream_preproc` (five simulated movies through the stage and its worker jobs),
  `test=lib_stream suite=pick_and_extract`, `test=lib_reconstruction suite=fraction_dose_export`.
- Alignment, sums or extraction: `test=simulated_workflow suite=6vxx` (or `1jxy`). Its map
  correlation is a robustness profile: one run shows breakage, not improvement.

## See also

`doc/policies/motion_gain_analysis_policy.md`, `doc/policies/refine_motion_model_policy.md`,
`doc/algorithms/motion_correction.md`; skills `simple-main-ctf`, `simple-main-pick`, `simple-main-stream`.
