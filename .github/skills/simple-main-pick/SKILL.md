---
name: simple-main-pick
description: Use when working in SIMPLE's src/main/pick subsystem, including segmentation-based (diameter) picking, batch picking with one diameter-bin and box decision, reference-based picking and its correlation kernel, Gaussian (moldiam) picking, the picker iterator, box files and their coordinate convention, and the picking used by pick, pick_extract, preprocessing and the stream's initial analysis.
---

# `src/main/pick`: particle picking

Owns finding particles on integrated micrographs and writing them as box files: the reference picker
(`pickref`, its BLAS correlation `pickref_corr_batch`, `multiref_merge` for several Gaussian diameters),
the segmentation picker (`picksegdiam`), its batch driver with one diameter-bin and box decision
(`segdiam_bin_picker`, the stream's picker in p03), and the per-micrograph dispatch (`picker_iter`,
`simple_picker_utils`). It does not make references (`make_pickrefs` in
`strategies/parallelization/simple_pick_strategy.f90`), loop micrographs or write part documents (that
strategy; `exec_pick_extract` in `commanders/simple/simple_commanders_pick.f90`, the p04 job), extract
(`strategies/parallelization/simple_extract_strategy.f90`, `preprocess/`), or pick-preprocess the
stream's first micrographs (`segdiampick_preprocess` in `stream/simple_mini_stream_utils.f90`).

## Rules

Silent ones first: nothing reports these.

- Box files: one line per particle, `x y box box` then `-3` (reference picker) or the diameter in A
  (segmentation). `x y` is the base-0 corner of a `box`-pixel window in pixels of the integrated
  micrograph at its own sampling; the centre is corner + box/2 (`preprocess/simple_particle_extractor`).
  A new box size keeps the centre: corner - (new - old)/2, as in `report_boxfile`, `segdiam_bin_picker`
  and extraction. Extraction takes the box from the first box file unless `box=` is given, and stops
  on mixed sizes. A column change reaches every `boxfile` reader (`rg boxfile src/main`).
- No picks, no file: with nptcls = 0 nothing is written and the returned name is empty;
  `sp_project%set_boxfile(i, name, nptcls=0)` sets `nptcls` 0 and no `boxfile` key. Consumers drop a
  micrograph without `nptcls` > 0 and a `boxfile` (`stream/shared/simple_mic_selection.f90`).
- `simple_pickref` holds the micrograph in module variables (`mic_raw`, `smpd_raw`, `ldim_raw`):
  `read_mic_raw_pickref` must precede `pickref%new` for each micrograph, and `new` silently uses the
  last one read. Reference picking takes one micrograph at a time per process, never several on
  threads; its threading is inside `new` and the distance filter.
- References must be at the micrograph's pixel size: `pickref%new` ignores their header and takes
  maxdiam = `smpd_raw` x reference box. `make_pickrefs` rescales them (`trust_header=yes`). The written
  box is `find_larger_magic_box` of 2 maxdiam / smpd (`BOX_EXP_FAC` = 1).
- Units: diameters, `moldiam`, `moldiam_max` and `thres` (peak separation) are in A, positions and
  boxes in pixels (`picksegdiam`'s default box `BOX_EXP_FAC*moldiam_max/smpd_raw`, the Gaussian
  references' `(1+BOX_EXP_FAC)*moldiam/smpd + 2`). Segmentation runs at 4 A/pixel; reference matching
  at 4 then 2 A/pixel, or the micrograph's if coarser. `picksegdiam` stores and writes diameters
  without the 8 A of its 2-pixel erosion: `segdiam_bin_picker` adds them back both to decide the bins
  and to assign each pick to one (its box files keep the stored diameter).
- `segdiam_bin_picker`: the first `pick` decides the accepted bins and the box, later calls reuse
  them; `new(accepted_bins, box)` supplies both (`THROW_HARD` if one comes alone). `pick` writes box
  files and thumbnails in the working directory, deletes from `os_mic` the micrographs with state 0,
  no picks or no box file, and writes the caller's project (a deviation listed in
  `doc/policies/layering_and_project_ownership_policy.md` section 8: add no more project writes here).
  Micrograph indices held across the call are stale.
- The `*_diams` text file (`picksegdiam%write_pos_and_diams`: `x y box box diam`) is the hand-off from
  pick-preprocessing (mic keys `mic_den`, `mic_topo`, `mic_bin`, `mic_diam`) to `segdiam_bin_picker`,
  which reads columns 1-3 and 5 and deletes the four keys. Change it in `picksegdiam`,
  `segdiam_bin_picker` and `stream/simple_mini_stream_utils.f90` together.

## Where to start

- `simple_picker_iter.f90`: `picker=new` with `pickrefs` (reference) or `moldiam`/`multi_moldiams`
  (Gaussian; returns no `thumb_den`), `picker=segdiam` (the thumbnail of its denoised micrograph, as
  the reference picker); `simple_picker_utils.f90` names the outputs and sets the box. The pick strategy
  stores `thumb_den` only when a picker returned one.
- `simple_segdiam_bin_picker.f90` (header has the lifecycle); its caller
  `stream/stages/simple_stream_stage_initial_analysis.f90`.
- `doc/algorithms/particle_picking.md` (both pickers), `doc/policies/stream/reference_generation_policy.md`,
  `doc/policies/segdiampick_mics_policy.md` (`mini_stream`'s estimate).

## Tests

- No unit sub-suite runs a picker. `build/bin/simple_test_exec test=lib_stream suite=pick_and_extract`
  (`pick_extract` on three planted particles: every centre within tolerance, box file, stack) and
  `suite=picking_references` (`make_pickrefs`).
- `test=unit_stream suite=reference_picking` and `suite=initial_analysis` cover the stages'
  bookkeeping only.
- Behaviour: `test=simulated_workflow suite=6vxx` (picks with `segdiam`; `picker=new` uses the
  reprojections as references), `test=mini_stream_6vxx`, `test=stream_chains` (the movie route runs
  the initial analysis's `segdiam_bin_picker`). `simulated_workflow` is a robustness profile: one run
  shows breakage, not improvement.

## See also

Skills `simple-main-stream`, `simple-main-image` (`image_bin`, segmentation), `simple-main-project`.
