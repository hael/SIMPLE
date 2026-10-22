# segdiampick_mics Policy

## 1. Purpose and Scope

`segdiampick_mics` in `src/main/stream/simple_mini_stream_utils.f90` estimates the
particle diameter from the first `mic_to` micrographs of a project, returns a box
size in pixels and a mask diameter in Angstroms, writes one box file per
micrograph, and writes the project's mic segment.

Its only caller is the `mini_stream` validation program
(`commanders/simple/simple_commanders_validate.f90`). The stream does not use it:
the stream picks with `segdiam_bin_picker` (`src/main/pick/simple_segdiam_bin_picker.f90`),
described in [stream/reference_generation_policy.md](stream/reference_generation_policy.md).
The `moldiam_max` argument is not read.

This document covers the diameter collection, the selection of diameters, the box
and mask derivation, the re-picking, the constants and the known defects.

## 2. Algorithm

### Stage 1: diameter collection

The micrographs are the project's active micrographs that have an `intg`
(`get_mics_table`); `mic_to` is cut to their number with a warning. For each of the
first `mic_to` micrographs:

1. If the mic row has `mic_den`, `mic_topo` and `mic_bin` (written by
   `segdiampick_preprocess` when preprocessing pre-picks), the three names are
   taken from the row and the diameters from the file named by `mic_diam` when it
   exists. The `mic_den`, `mic_topo`, `mic_bin` entries, and `mic_diam` when its
   file exists, are deleted from the row after the import.
2. Otherwise `picksegdiam%pick` segments and picks the raw micrograph, with the
   diameter bins `MOLDIAMS_PICK`; the density, topographic and binary images are
   named by `append2basename` of the micrograph. A micrograph without boxes is
   skipped.

All diameters go into one array. `2 * SMPD_SHRINK1` Angstroms are then added once to
every diameter, for the 2-pixel erosion of the binarization.

### Stage 2: statistics

`calc_stats` gives the mean, median, standard deviation, minimum and maximum of the
corrected diameters. `mad = mad_gau(diams, median)` is `1.4826 * MAD`.

### Stage 3: selection

`sortmeans` partitions the diameters into `NQ_DIAMS = 10` clusters. For each
cluster the absolute z-score of its mean is `|mean_i - median| / mad`. Clusters 1
and 10 are never accepted. Clusters 2 to 9 are accepted when their z-score is below
`SIGMA_CRIT`. The diameters of the accepted clusters form the selected set, and
`calc_stats` is applied to it again.

### Stage 4: box size

```
box_in_pix = find_magic_box( BOXFAC * nint( median_selected / smpd ) )
```

with the median of the selected set and `BOXFAC = 3`.

### Stage 5: mask diameter

```
mskdiam = min( (box_in_pix - COSMSKHALFWIDTH) * smpd, max_selected )
mskdiam = min( mean_selected + SIGMA_CRIT_MSK * sdev_selected, mskdiam )
```

### Stage 6: re-picking and output

Every micrograph of `1 .. mic_to` is picked again from its binary image with the
diameter range `[min_selected, max_selected]`. A micrograph with boxes gets a box
file (named after the micrograph, with `box_in_pix` as box) registered with
`set_boxfile`. When the density image exists, its thumbnail (`_den.jpg`) is stored
in the row as `thumb_den`. The mic segment is written with
`write_segment_inside('mic')` to the project's own file.

## 3. Constants

| Name | Value | Role |
|---|---|---|
| `SMPD_SHRINK1` | 4.0 A | sampling of the segmentation; twice it is the erosion correction |
| `NQ_DIAMS` | 10 | number of `sortmeans` clusters |
| `SIGMA_CRIT` | 2.0 | accepted-cluster threshold on the z-score of the cluster mean |
| `BOXFAC` | 3 | box factor on the median selected diameter, in pixels |
| `SIGMA_CRIT_MSK` | 2.5 | mask cap, in standard deviations above the mean selected diameter |

All five are active. `COSMSKHALFWIDTH` is the global soft-edge half-width.

## 4. Invariants

- The erosion correction is applied exactly once, before any statistic.
- The box, the mask diameter and the re-picking range come from the selected set,
  not from the full array.
- Pre-picked `mic_den`, `mic_topo` and `mic_bin` entries are deleted when imported
  and are not used by a later call.

## 5. Ownership

`segdiampick_mics` owns the estimate, the selection, the box and mask derivation, the
re-picking and the write of the mic segment. `picksegdiam` owns the segmentation,
the picking and the per-micrograph diameter files (`get_diameters`,
`write_pos_and_diams`). `sortmeans` (`simple_math.f90`) owns the quantization,
`calc_stats` and `mad_gau` (`simple_stat.f90`) the statistics, `find_magic_box` the
FFT-friendly rounding.

## 6. Known defects

These are the current behaviour.

- The import of a `mic_diam` file reads only its first line. The file has one box per
  line (`x y box box diameter`, `write_pos_and_diams`), so the five values of the
  first box enter the statistics as five diameters and the other boxes are ignored.
- There is no guard on an empty diameter array. When no micrograph yields a box the
  array is not allocated at the erosion correction.

## 7. Do Not

- Do not change the erosion correction without the segmentation picker's
  binarization.
- Do not use the full-array statistics for the box, the mask diameter or the
  re-picking range.
