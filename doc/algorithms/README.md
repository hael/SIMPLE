# Algorithms

Mathematical and operational descriptions of SIMPLE's critical production
algorithms, written for a computer-science-literate reader. Each chapter states
what is computed, which objects own the result, and which guards make the
workflow well-defined. Source locations appear at the end of each chapter;
detailed lifecycle and refactoring constraints remain in `doc/policies` and
`doc/implementation_notes`.

## Scope

This series includes a workflow or operator when it materially determines a
corrected image, particle coordinate, pose, class, reconstructed map,
resolution estimate, or heterogeneity state. It does not attempt to document
GUI plumbing, schedulers, file formats, performance-only rewrites, isolated bug
fixes, or research paths without production callers.

The word *production* means that the path has a public or workflow-owned caller
in the current tree. Opt-in methods are labeled. A design note or test harness
alone is not enough.

## End-to-end pipeline

1. [Motion correction](motion_correction.md) aligns dose-fractionated movie
   frames globally and locally.
2. [CTF modeling and estimation](ctf_estimation.md) fits defocus,
   astigmatism, and optional phase shift from micrograph power spectra.
3. [Particle picking](particle_picking.md) finds candidate coordinates by
   segmentation or exhaustive reference matching.
4. [Cluster2D and class averaging](cluster2d_class_averaging.md) assigns
   particles to aligned 2D classes and restores CTF-aware class averages.
5. [Ab initio 2D](abinitio2d.md) stages 2D classification from initial
   references to dense final assignments.
6. [Sampling and fractional updates](sampling_and_fractional_updates.md)
   separates outer particle sampling from inner candidate sampling and defines
   class/volume restoration.
7. [Ab initio 3D](abinitio3d.md) frequency-marches particle orientations and
   state assignments from starting models to inspectable 3D volumes.
8. [Refine3D](refine3d.md) iterates reference preparation, probabilistic
   pre-alignment, hard pose assignment, and volume assembly.
9. [3D reconstruction](reconstruction.md) describes the Cartesian gridding and
   opt-in PCG backends behind the dense half-map seam.
10. [Nonuniform filtering](nonuniform_filtering.md) selects a spatially varying
    local low-pass field and hands its frontier to matching.
11. [NU-evidence envelope masking](nu_evidence_envelope_mask.md) derives a
    reproducibility-based molecular envelope from the NU candidate bank.
12. [Heterogeneous refinement](heterogeneous_refinement.md) describes the
    production multi-state wrappers and their all-particle final maps.
13. [Flex PCA](flex_pca.md) estimates a projection-aware low-rank covariance
    embedding and state initializers.
14. [Streaming](streaming_pipeline.md) composes preprocessing, optics,
    picking, sieving, and global 2D classification as an online artifact flow.

[Continuous in-plane refinement](continuous_inplane_refinement_abinitio2D.md)
is a cross-cutting pose-polishing chapter used by 2D and 3D matching rather
than a separate pipeline stage.

## Shared conventions

- Particle images are compared with references in Fourier space. A pose
  contains a projection direction, an in-plane angle, and a 2D shift; a
  multi-state pose also contains a state label.
- Even and odd halfsets remain separate wherever half-map independence or FRC/
  FSC estimation is required. Merged products are derived from the halves.
- `sampled` identifies the current outer particle subset and `updatecnt`
  records cumulative update history. Downstream restoration uses realized
  state, not just the requested `update_frac`.
- Probabilistic pre-alignment samples a bounded discrete candidate set. The
  matcher then commits one hard assignment. This is not a full soft-assignment
  EM reconstruction.
- Low-pass values are resolutions in Angstrom; a smaller value admits a higher
  spatial frequency. Crop changes preserve the physical field of view:
  `box * smpd = box_crop * smpd_crop`.
- Reconstruction, filtering, masking, and FSC estimation have separate
  ownership. A mask or bandwidth selected for one role is not automatically
  valid for another.
- Shared-memory and distributed routes are intended to preserve the same
  scientific workflow. Process lifetime and reduction mechanics may differ;
  mathematical inputs and durable artifacts may not.

## Production boundaries

The continuous five-parameter 3D pose solver, generative flex reconstruction,
and experimental PCG prior research are intentionally absent as production
chapters. Their design and validation records remain in
`doc/implementation_notes` until a production caller and scientific acceptance
gate exist.

