# SIMPLE: the algorithms

SIMPLE is a single-particle cryo-EM reconstruction platform. Its input is a
set of electron micrographs (or the movies they are integrated from) in which
many copies of one macromolecule lie in random orientations; its output is a
3D density map of that molecule, together with the per-particle parameters
that produced it. Everything between input and output is estimation: of
motion, of the microscope's transfer function, of particle positions, of
poses, of the map itself, and of how reliable each of those estimates is.

This series describes those estimators. It is written for a reader who is
comfortable with Fourier transforms, least squares, and Monte Carlo methods
but who has never opened the source tree, and it deliberately says nothing
about how the code is organized. Each chapter states what quantity is
estimated, from which data, under which model, by which update rule, and why
the rule has the shape it has.

## The problem in one equation

A particle image `y_i` is a noisy 2D projection of an unknown 3D density `x`,
seen through the microscope's contrast transfer function:

```text
y_i = C_i S_i P(R_i) x + n_i,        n_i ~ N(0, Sigma_i).
```

`P(R)` projects along the direction given by rotation `R`, `S` translates by an
unknown in-plane shift, `C` applies the CTF, and `n` is noise whose power
`sigma2_i(k)` varies with spatial frequency `k`. The unknowns are `x` and, for
every one of `10^5` to `10^7` particles, its pose `(R_i, shift_i)`. Given the
poses, `x` is a linear least-squares problem; given `x`, each pose is a small
nonlinear registration problem. Nobody knows either at the start.

Almost everything in the library follows from four decisions about how to
attack that equation.

**1. Work in Fourier space.** By the projection-slice theorem, the Fourier
transform of a projection is a central section of the 3D transform. A
projection is therefore a *gather* from the 3D Fourier grid, reconstruction is
the adjoint *scatter*, a shift is a phase ramp, the CTF is a multiplication,
and, with references stored in polar coordinates, an in-plane rotation is an
index shift and all rotations at once are one FFT. The noise is diagonal in
this basis, so whitening is a per-shell weight. Chapters
[Cluster2D](cluster2d_class_averaging.md), [refine3D](refine3d.md), and
[reconstruction](reconstruction.md) are the three faces of this choice.

**2. Alternate.** Poses and maps are estimated by coordinate ascent: fix the
map, assign each particle one pose; fix the poses, reconstruct the map. The
same two-step structure appears in 2D, where the "map" is a set of class
averages and the pose has no out-of-plane component, and in the continuous
heterogeneity model, where the pose is fixed and the alternation is between a
latent coordinate and a basis of volumes ([flex PCA](flex_pca.md)).

**3. Continue from coarse to fine.** The alternation converges to the
nearest local optimum, and from a random start that optimum is bad. Every
workflow therefore controls three quantities over time: the admitted
bandwidth (low-pass limit), the stochasticity of the pose search (random
visiting order, first-improvement acceptance, sampled rather than argmax
choices, annealed toward exhaustive search), and the fraction of particles
updated per iteration. [Ab initio 2D](abinitio2d.md) and
[ab initio 3D](abinitio3d.md) are those schedules; the stochastic machinery
they share is in [sampling and fractional updates](sampling_and_fractional_updates.md).

**4. Split the data in half.** Every reconstruction is made twice, from
disjoint even and odd particle sets, and every claim about resolution, every
filter, every regularizer, and every mask is derived from the agreement
between the halves. This is what turns a fitted map into a measured one. The
half-map comparison is the Fourier shell correlation in
[refine3D](refine3d.md); the same idea made local gives
[nonuniform filtering](nonuniform_filtering.md) and the
[NU-evidence envelope](nu_evidence_envelope_mask.md).

## Two objectives

The per-particle comparison used throughout is one of

```text
Euclidean:    L_i(pose)  = sum_k (k/sigma2_i(k)) |y_i(k) - C_i(k) [S P x](k)|^2
                           / sum_k (k/sigma2_i(k)) |y_i(k)|^2
Correlation:  cc_i(pose) = <y_i, C_i [S P x]> / (||y_i|| ||C_i [S P x]||)
```

summed over the Fourier band currently admitted. The Euclidean form is the
whitened negative log-likelihood, normalized by the particle's own energy so
that particles are on one scale; that common scale is what allows the
importance weights `exp(-L)` of the probabilistic modes and the signal priors
of ML regularization to use the same numbers. The correlation form is
scale-free and is used before a noise model exists.

Assignment is always hard. Probabilistic modes draw one candidate from a
bounded importance distribution and commit it; no chapter describes a
marginalization over poses. What the probabilistic modes add is a *global*
step between scoring and committing, a balanced stochastic assignment that
stops a few references from absorbing the data.

## Reading order

The chapters follow the data through a project. Each is self-contained and
uses the outline Problem, Model, Algorithm, Rationale, Implementation, with
source pointers confined to the last section.

**Preprocessing: from movie to particle.**

1. [Motion correction](motion_correction.md). A translation per frame from
   leave-one-out correlation, then an 18-term space-time polynomial for local
   deformation, accepted only if it reproduces the measured patch motion.
2. [CTF estimation](ctf_estimation.md). Defocus, astigmatism, and phase
   shift fitted to the Thon rings by multi-start differential evolution and
   gradient refinement, with a resolution and ice diagnostic.
3. [Particle picking](particle_picking.md). Reference-free segmentation to
   bootstrap, then exhaustive Pearson correlation against class averages with
   non-maximum suppression under a physical exclusion radius.

**2D: classes without a model.**

4. [Cluster2D and class averaging](cluster2d_class_averaging.md). The basic
   alternating estimator: stochastic hill climbing over classes and in-plane
   poses, then Wiener-type restoration of CTF-corrected class means.
5. [Ab initio 2D](abinitio2d.md). The coarse-to-fine schedule that takes
   Cluster2D from random noise references to a dense all-particle assignment.

**The stochastic machinery.**

6. [Sampling and fractional updates](sampling_and_fractional_updates.md).
   Which particles are updated per iteration, how candidate poses are drawn
   from importance weights and assigned in a balanced way, and how a partial
   update is blended into the running estimate without losing normalization.

**3D: from noise to map.**

7. [Ab initio 3D](abinitio3d.md). Eight stages of bandwidth, angular
   sampling, search mode, and symmetry, starting from a noise volume and
   random orientations.
8. [Refine3D](refine3d.md). One round of 3D alternation: reproject, assign,
   reconstruct, measure the FSC, regularize, filter.
9. [3D reconstruction](reconstruction.md). The fixed-pose inverse problem by
   Kaiser-Bessel gridding, or by preconditioned conjugate gradients on the
   exact normal equations.

**Resolution as a local quantity.**

10. [Nonuniform filtering](nonuniform_filtering.md). A spatially varying
    cutoff chosen by cross-half prediction error under a robust loss and an
    ordered smoothness prior.
11. [NU-evidence envelope masking](nu_evidence_envelope_mask.md). The same
    evidence turned into a molecular envelope by a two-label Markov random
    field.

**Heterogeneity.**

12. [Heterogeneous refinement](heterogeneous_refinement.md). The 3D
    estimator with a discrete state label per particle, and how states are
    initialized from external references without importing their biases.
13. [Flex PCA](flex_pca.md). Continuous variability as a low-rank volume
    covariance fitted by EM through the projection operator, and its latent
    coordinates turned into state maps.

**Online processing.**

14. [Streaming](streaming_pipeline.md). Steps 1 to 5 recast as an online
    algorithm: bounded chunks, a learned class-quality sieve, and a growing
    pool whose coarse-to-fine schedule is written in iterations rather than
    stages.

Two cross-cutting chapters stand outside the pipeline order.
[Continuous in-plane refinement](continuous_inplane_refinement_abinitio2D.md)
polishes a committed pose to sub-grid angular precision using the fact that
the polar correlation is a trigonometric series in the angle, while leaving
the search that selected the pose untouched.
[Automatic gain-flip detection](automatic_gain_flip.md) is a small
calibration estimator that decides the orientation of the detector gain
reference from the first movies of a session.

## Conventions

- A pose is a projection direction, an in-plane angle, and a 2D shift; in
  multi-state problems it also carries a state label.
- Low-pass limits are quoted in Angstrom; a smaller value admits higher
  spatial frequency. Cropping preserves the physical field of view,
  `box * smpd = box_crop * smpd_crop`, so a crop changes sampling, not content.
- Merged products (maps, class averages, filters) are derived from the even
  and odd halves, never the reverse.
- Fractional updates act on a subset of particles per iteration; restoration
  blends the subset's contribution with the previous estimate using the
  realized fraction per class or per state, not the requested one.
- Named constants and defaults are those of the current tree. Where a
  behavior is opt-in it is labeled as such; research paths without a
  production caller are not described.

## Where the rest lives

The series says nothing about lifecycle, ownership, execution routing, or
file formats. Those are in `doc/policies` (the rules a subsystem must
preserve), `doc/implementation_notes` (design and validation records for
individual features), `doc/refactoring_notes` (records of structural
changes), and `doc/code_overview` (a map of the source tree).
