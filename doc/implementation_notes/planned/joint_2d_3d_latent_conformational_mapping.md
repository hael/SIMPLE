# Joint 2D/3D Latent Affinity Propagation

**Date:** July 21, 2026
**Status:** Concise implementation proposal

## Goal

Combine two particle-level descriptions of conformational variation:

- local 2D diffusion coordinates from groups of nearby projection directions;
- global 3D coordinates and uncertainties from `flex_analysis`.

The result should be one conformational label per particle and one validated 3D
map per state. SIMPLE's existing affinity-propagation (AP) implementation
should perform the clustering.

## Core Idea

Use AP twice:

```text
2D diffusion coordinates + 3D flex coordinates
                    |
          local AP in each angular patch
                    |
          view-specific particle exemplars
                    |
       global AP over the local exemplars
                    |
       states containing exemplars from many views
                    |
             state-specific 3D maps
```

The first pass detects conformational structure where projections are directly
comparable. The second pass identifies which local clusters represent the same
3D state.

## Terminology

Three different objects must not be confused:

- **Orientation bin:** an atomic, disjoint projection-direction assignment.
  Each particle has one orientation bin.
- **Angular patch:** a temporary neighborhood made from several nearby bins.
  Patches overlap so neighboring analyses share particles.
- **Conformational state:** the final global AP cluster. Each particle has one
  hard state label.

For example:

```text
orientation bins:      A   B   C   D   E

patch centered at B:  [A   B   C]
patch centered at C:      [B   C   D]
patch centered at D:          [C   D   E]
```

A particle in bin `C` is analyzed in three local diffusion charts. It may have
three local coordinate vectors and three temporary local labels, but it still
has one orientation and one final conformational state.

The project stores the particle once. Patch membership is represented by
particle-index lists, so images need not be duplicated. When evidence from
overlapping patches is combined, per-particle patch weights satisfy

```text
sum[g in G_i] omega_ig = 1
```

and the particle is counted once in reconstruction.

## Inputs

For particle `i`:

- `m_i`: canonical `flex_analysis` coordinate;
- `S_i`: posterior covariance of `m_i`;
- `u_ig`: diffusion coordinate in angular patch `g`;
- `pind_i`: stable project particle index.

The `flex_analysis` basis volumes `B_q` and consensus map `V_mean` are also
required to map state coordinates back to volumes.

All joins use project particle index, never row order.

## Algorithm

### 1. Build Overlapping Angular Patches

Use the existing projection grid as patch centers. Each patch contains its
center bin and nearby orientation bins.

Within a patch:

1. align/transport particles to the patch center;
2. subtract the projected consensus map;
3. use the same CTF/noise convention as `flex_analysis`;
4. construct a diffusion map from the residual images.

Orientation should select candidate neighbors, but residual-image similarity
should determine graph weights. The current `ptcl3D` default `graph=ori`
weights edges using orientation distance alone and therefore describes viewing
geometry rather than conformation. Until a hybrid graph exists, `graph=euc`
on properly transformed residual images is the better starting point.

Persist each chart's particle indices, coordinates, eigenvalues, graph
diagnostics, and patch weights. Do not immediately discard the coordinates
after local clustering.

### 2. Local Joint Affinity Propagation

For particles `i` and `j` in patch `g`, combine normalized 2D and 3D distances:

```text
D_g(i,j) = alpha_2 D_2g(i,j) + alpha_3 D_3(i,j)
```

- `D_2g` is diffusion distance inside patch `g`.
- `D_3` is a distance between `flex_analysis` coordinates, preferably using
  `S_i + S_j` so uncertain particles contribute less strongly.

Run the existing AP wrapper:

```fortran
call cluster_dmat(D_g, 'aprop', nlocal, centers, labels, ap_pref=pref)
```

Use the median off-diagonal similarity as the initial preference. Retain:

- the exemplar particle index;
- local member particle indices;
- patch identifier;
- preference and similarity sum;
- 3D mean/covariance of the local members.

Each local AP cluster is now a meta-item for the global pass.

### 3. Global Affinity Propagation

Compare local exemplars `a` and `b` using:

```text
D_meta(a,b) = beta_3 D_3cluster(a,b)
            + beta_o H_ab [1 - overlap_agreement(a,b)]
```

- `D_3cluster` compares their member distributions in the global 3D latent
  space.
- `overlap_agreement` measures whether particles shared by the two patches
  retain the same local membership.
- `H_ab` is zero when the patches have insufficient overlap, making missing
  overlap evidence neutral.

Run AP on `D_meta`. Its output maps every local exemplar to a global
conformational state.

A global state is therefore an angular atlas of view-specific exemplars. The
single meta-exemplar selected by AP identifies the cluster, but final particle
assignment uses the state's complete exemplar set rather than only that one
view.

### 4. Assign Particles and Create 3D Maps

Map each particle's local AP labels through the global exemplar labels. If
overlapping patches disagree, resolve the assignment using the closest
available state exemplars plus the particle's 3D coordinate.

The baseline output is one hard AP state per particle. An optional soft weight
may be derived from the final state-distance margin for reconstruction, but it
must be described as a similarity weight, not a Bayesian posterior.

For state `k`, calculate its mean 3D latent coordinate `mu_k` and form the fast
eigenvolume prediction:

```text
V_k = V_mean + sum_q mu_k(q) B_q
```

Then reconstruct even/odd maps directly from particles assigned to the state.
The particle reconstructions are the physical validation target; the
eigenvolume maps are model predictions.

## Why Two AP Passes?

A single AP pass over all particles is unsuitable because:

1. diffusion coordinates are comparable only inside their own angular patch;
2. one particle exemplar represents only one projection direction;
3. the existing AP implementation uses dense `N x N` matrices.

The two-pass design keeps local comparisons valid, represents each state with
many viewing directions, and confines dense AP matrices to patches and the much
smaller exemplar set.

The current `aff_prop` object stores seven dense real matrices. Through
`cluster_dmat`, the input distance and converted similarity matrices are also
live. A conservative memory estimate is therefore at least

```text
36 * N * N bytes
```

per AP call. Check this before allocation and fail clearly when a patch or
meta-exemplar set is too large.

## Selecting the Number of States

AP determines cluster count through its diagonal preference:

- start with the existing median-similarity preference;
- scan a small range of global preference quantiles;
- record preference, exemplar count, convergence warning, and similarity sum;
- select the simplest solution stable across half-sets and angular holdouts.

Do not normally use `nclust_max`: the current wrapper enforces that limit with
a k-medoids merge, so the final result would no longer be a pure AP solution.

## SIMPLE Changes

### Reuse

- `simple_aff_prop.f90`: existing AP message passing and deterministic
  tie-breaking;
- `simple_clustering_utils.f90`: `cluster_dmat(..., 'aprop', ...)` and AP
  preference handling;
- `simple_projected_latent_result.f90`: in-memory particle indices, `z`,
  posterior covariance, residuals, and mode variances.

### Add or Refactor

- `simple_cls_split_strategy.f90`
  - separate diffusion embedding from k-medoids;
  - return non-destructive, particle-indexed local charts;
  - generate overlapping angular patches.
- `simple_diff_map_graphs.f90`
  - add orientation-gated, residual-image-weighted affinities.
- `simple_projected_latent_result.f90`
  - serialize posterior covariance as well as coordinates.
- new `simple_joint_latent_aff_prop.f90`
  - joint distances;
  - local and global AP;
  - preference scans;
  - final assignments and diagnostics.
- new high-level `flex_map3D` workflow
  - orchestrate features, AP, mapping, reconstruction, and validation.

The existing AP API may need small additions to expose convergence status,
memory estimates, or per-item preferences. The clustering algorithm itself
does not need to be replaced.

## Outputs

```text
flex_map3D_particles.txt
    particle, z..., state, assignment_margin, optional_weights..., residual

flex_map3D_exemplars.txt
    exemplar_particle, patch, local_population, global_state

flex_map3D_states.txt
    state, population, z_centroid..., validation_metrics

flex_map3D_patches.txt
    patch, population, rank, AP_preference, exemplar_count, overlap_consistency

flex_map3D_stateNN_eigenvol.mrc
flex_map3D_stateNN_even.mrc
flex_map3D_stateNN_odd.mrc
```

Preserve existing project `class` fields by default, and never overload the
record-activation `state` field. Updating project assignments must be explicit.

## Minimal Implementation Sequence

1. Preserve/export local diffusion charts and `flex_analysis` covariance.
2. Implement local joint AP and write exemplars without changing the project.
3. Implement global exemplar AP and final particle labels.
4. Write eigenvolume-predicted maps.
5. Reconstruct and validate even/odd state maps.
6. Add preference scanning and optional similarity-weighted reconstruction.

## Required Validation

- known-state simulation with realistic orientations, CTFs, and noise;
- homogeneous null data, which must not yield reproducible state differences;
- even/odd reproducibility of state assignments and maps;
- angular holdout and patch-dropout stability;
- 2D-only, 3D-only, and joint-distance ablation;
- repeated AP runs with identical exemplars and convergence;
- comparison of eigenvolume-predicted maps with particle reconstructions.

## Main Safeguards

- Use orientation only for locality; use image residuals for conformational
  affinity.
- Normalize overlapping patch weights so particles are not double counted.
- Keep feature extraction fixed during the first AP implementation.
- Represent a global state with all its view-specific exemplars.
- Bound dense AP memory before allocation.
- Export results before allowing destructive project updates.

## Summary

The method is cluster tracking with a global 3D anchor:

1. overlapping angular patches supply local diffusion features;
2. local AP selects particle exemplars using joint 2D/3D distance;
3. global AP groups those exemplars into conformational states;
4. the `flex_analysis` basis maps each state into 3D;
5. particle reconstructions validate the maps.

Particles may appear in several angular analyses, but each particle has one
orientation, one final conformational state, and one contribution to each
reconstruction.
