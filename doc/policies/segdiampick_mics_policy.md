# segdiampick_mics Policy

## 1. Purpose and Scope

This document defines the algorithmic and implementation policy for
`segdiampick_mics` in `simple_mini_stream_utils.f90`.

`segdiampick_mics` estimates particle diameter from a mini-batch of
micrographs, returns a box size in pixels and a mask diameter in Ångströms
for downstream picking and refinement, and writes particle box positions to
the project. It is the diameter-estimation entry point for the mini-batch
stream pipeline.

This document covers:

- the diameter collection stage
- the statistical filtering stage (clustering and outlier rejection)
- the box-size and mask-diameter derivation rules
- the re-picking stage
- parameters and their roles
- the two clustering backends (k-means and HAC)

## 2. Algorithm

### Stage 1 — Diameter collection

For each micrograph in `1 .. mic_to`:

1. If the project already records `mic_den`, `mic_topo`, and `mic_bin` for a
   micrograph (written by `segdiampick_preprocess`), import diameters from the
   corresponding `mic_diam` text file. Delete the `mic_den`, `mic_topo`,
   `mic_bin`, and `mic_diam` entries from the project after import so they are
   not reused on a subsequent call.
2. Otherwise, run `picksegdiam%pick` against the raw micrograph. Retrieve
   diameters via `picksegdiam%get_diameters`. Micrographs that yield zero
   boxes are skipped silently.

All per-micrograph diameter vectors are concatenated into a single flat array
`diams_arr`.

A correction of `2 × SMPD_SHRINK1` Ångströms is added to every collected
diameter to compensate for the two-pixel erosion step performed during binary
mask generation in the segmentation picker.

### Stage 2 — Global diameter statistics

`calc_stats` is called on `diams_arr` to populate a `stats_struct` containing
mean, median, standard deviation, minimum, and maximum. The global
robust-scale estimate `mad` is computed via `mad_gau(diams_arr, median)`,
which returns `1.4826 × MAD`.

### Stage 3 — Clustering and outlier rejection

Diameters are clustered to identify and discard outlier size populations.
Two backends are available, selected by the optional `l_hac` argument.

#### 3a. k-means path (default, `l_hac` absent or `.false.`)

`sortmeans` partitions `diams_arr` into exactly `NQ_DIAMS = 10` clusters.
The first and last quantile (indices 1 and `NQ_DIAMS`) are always excluded
from the accepted set regardless of their Z-scores; they represent the
extreme-size tail populations.

#### 3b. HAC path (`l_hac = .true.`)

`hac_1d_fast` performs threshold-based hierarchical agglomerative clustering
on `diams_arr` with a linkage threshold equal to the global standard deviation
`diam_stats%sdev`. The number of clusters is adaptive and determined by the
data. All clusters are eligible for acceptance; there is no mandatory exclusion
of edge clusters.

#### 3c. Cluster filtering (both paths)

Two independent Z-score criteria are evaluated per cluster, both must pass for
a cluster to be accepted:

**Mean Z-score** — the absolute Z-score of the cluster centroid with respect
to the global median and MAD:

$$|z_i| = \frac{|\mu_i - \text{med}(\text{diams})|}{1.4826 \cdot \text{MAD}(\text{diams})}$$

Clusters with $|z_i| \geq \sigma_{\text{crit}}$ are rejected. These represent
populations whose mean diameter is a statistical outlier from the bulk.

**Population MAD Z-score** — the Z-score of the cluster population count with
respect to the median and MAD of the population-count vector across all
clusters:

$$z_{\text{pop},i} = \frac{\text{med}(\text{pops}) - \text{pop}_i}{1.4826 \cdot \text{MAD}(\text{pops})}$$

The sign convention is positive = below-median population = sparse cluster.
Clusters with $z_{\text{pop},i} \geq \sigma_{\text{crit}}$ are rejected. This
discards anomalously small clusters that are likely artefacts or contaminants
rather than the target particle population.

If the population MAD is zero (all cluster populations are equal), the
population Z-score criterion is suppressed (all clusters pass the population
test). The mean Z-score criterion still applies.

The threshold for both criteria is `SIGMA_CRIT = 2.0`.

Diameters belonging to accepted clusters are collected into `diams_arr_ts`
(threshold-selected). Post-filtering statistics are printed for diagnostic
purposes.

### Stage 4 — Box size derivation

The box size in pixels is:

```
box_in_pix = find_magic_box( nint( 1.5 × diam_stats%maxv / smpd ) )
```

where `diam_stats%maxv` is the **maximum** diameter in the threshold-selected
set and `smpd` is the pixel size in Ångströms. `find_magic_box` rounds up to
the nearest FFT-friendly size. The 1.5× factor provides clearance around the
largest accepted particle.

### Stage 5 — Mask diameter derivation

```
mskdiam = min( (box_in_pix - COSMSKHALFWIDTH) × smpd,
               diam_stats%maxv × 1.2 )
```

`COSMSKHALFWIDTH` is a global constant that reserves a cosine-edge half-width
at the box boundary. The mask diameter is capped at 1.2× the maximum accepted
diameter. The stricter of the two bounds is taken.

### Stage 6 — Re-picking with diameter constraints

After diameter estimation, all micrographs in `1 .. mic_to` are re-picked
with `picksegdiam%pick` using the diameter range
`[diam_stats%minv, diam_stats%maxv]` from the threshold-selected set as an
acceptance gate. Box positions and diameters are written to per-micrograph
`.box` files and registered in the project via `spproj%set_boxfile`.

Density thumbnails (`.jpg`) are generated for micrographs whose density image
(`*_DEN.mrc`) exists on disk and are stored in the project mic segment.

The mic segment is written to disk at the end of the subroutine.

## 3. Parameters

| Symbol | Value | Role |
|---|---|---|
| `SMPD_SHRINK1` | 4.0 Å | Shrinkage factor used in segmentation picking; applied as erosion correction |
| `NQ_DIAMS` | 10 | Number of k-means clusters in the default path |
| `SIGMA_CRIT` | 2.0 | Z-score rejection threshold for both mean and population criteria |
| `SIGMA_CRIT_MSK` | 2.5 | Reserved for alternative mask-diameter clipping rule; currently inactive |
| `BOXFAC` | 3 | Reserved for alternative box-size rule; currently inactive |

## 4. Invariants

- `diams_arr` must be non-empty before the clustering stage is reached. No
  guard is currently applied; callers must ensure `mic_to` covers at least one
  micrograph with detectable particles.
- The erosion correction (`+2 × SMPD_SHRINK1`) is applied exactly once, before
  any statistics are computed.
- The clustering stage operates on the corrected diameters. Box-size and
  mask-diameter derivation operate on the **post-filtering** statistics
  (`diams_arr_ts`), not the pre-filtering full array.
- Re-picking uses the diameter bounds from the threshold-selected statistics,
  not the pre-filtering bounds.
- Pre-picked `mic_den`/`mic_topo`/`mic_bin`/`mic_diam` project entries are
  deleted after import and must not be reused.

## 5. Ownership

`segdiampick_mics` owns the full estimation-to-picking pipeline for a
mini-batch: diameter collection, filtering, box/mask derivation, re-picking,
and project output.

`picksegdiam` owns the segmentation-based picking and the per-micrograph
diameter file I/O.

`hac_1d_fast` (in `simple_math.f90`) owns the 1D HAC clustering algorithm.
`sortmeans` (in `simple_math.f90`) owns the k-means quantization. Neither
clustering routine owns the filtering logic; the Z-score thresholding lives
entirely in `segdiampick_mics`.

`find_magic_box` owns FFT-friendly box-size rounding.

`calc_stats` and `mad_gau` own the statistical summaries. `median_nocopy`
owns non-destructive median computation for the population-count vector.

## 6. Do Not

- Do not change the erosion correction factor without updating the segmentation
  picker's binarization step correspondingly.
- Do not apply the diameter bounds from the pre-filtering `diam_stats` to the
  re-picking stage; always use post-filtering bounds.
- Do not add a mandatory edge-cluster exclusion to the HAC path; the adaptive
  cluster count makes fixed-index exclusion meaningless.
- Do not suppress the population Z-score when MAD is nonzero; only suppress it
  when MAD is exactly zero (all populations equal).
- Do not collapse the two clustering paths into a single code path that
  post-processes `sortmeans` output with HAC; they are independent policies.
