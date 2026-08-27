# SIMPLE Algorithm Descriptions

## Status

Completed on 2026-08-27. This is the single living design, review, and
validation record for the mathematical algorithm series in `doc/algorithms`.

## Contract

The series follows the editorial model in `/Users/elmlundho/src/X/docs/algorithms`:

- write for a computer-science-literate reader;
- state the problem, mathematical objects, transformations, selection rules,
  guards, and handoffs before discussing implementation;
- distinguish production behavior from optional, experimental, and proposed
  behavior;
- distinguish scientific semantics from execution and performance mechanics;
- keep implementation pointers in a short final section;
- link to policy or implementation notes for detailed lifecycle constraints.

A chapter is included when a production workflow or operator materially
determines corrected images, particle coordinates, poses, classes, maps,
resolution estimates, or heterogeneity states. GUI plumbing, schedulers,
file-format support, performance-only rewrites, isolated bug fixes, and
unimplemented research designs are outside the series.

## Coverage

The collection retains and indexes the existing chapters on motion correction,
continuous in-plane refinement, and the NU-evidence envelope. New chapters
cover:

1. CTF modeling and estimation.
2. Segmentation and reference-based particle picking.
3. Cluster2D and CTF-aware class averaging.
4. Staged ab initio 2D classification.
5. Probabilistic pre-alignment, outer sampling, and fractional restoration.
6. Staged ab initio 3D initialization.
7. Base refine3D projection matching and volume assembly.
8. Cartesian gridding and PCG reconstruction.
9. Nonuniform filtering and matching-bandwidth handoff.
10. Multi-state and competitive heterogeneous refinement.
11. Projection-aware flex PCA initialization.
12. The streaming preprocessing-to-global-2D pipeline.

The collection does not describe continuous five-parameter 3D pose polishing,
generative flex reconstruction, or proposed PCG priors as production features.
Those remain implementation/design notes until their production callers and
scientific acceptance gates exist.

## Review Findings

- SIMPLE is not one linear algorithm. The useful organizing spine is the
  end-to-end SPA workflow, with cross-cutting chapters for probabilistic
  sampling, reconstruction, and heterogeneity.
- Probabilistic modes implement bounded candidate sampling followed by a hard
  particle update. They are not a volume-integrated soft-assignment EM method.
- `update_frac` selects an outer particle subset. Candidate importance sampling
  happens inside that fixed subset, and restoration consumes realized
  `sampled`/`updatecnt` state.
- Reference-picker SGEMM batching changes evaluation cost but preserves the
  exhaustive Pearson search and coordinate-selection semantics.
- The PCG backend solves a fixed-pose weighted least-squares system. It does not
  optimize particle poses, and its base half maps remain the FSC authority.
- NU filtering, NU-evidence masking, density FSC masking, and spherical NU
  support are different objects and must not be collapsed in prose.
- Multi-state wrappers own stage and sampling policy; base `refine3D` continues
  to own probability tables, hard assignment, partial reconstruction, and
  volume assembly.
- `flex_pca` is production heterogeneity initialization. The removed
  `flex_analysis` pipeline and proposed generative-flex workflows are not.

## Implementation Plan

1. Add an indexed `doc/algorithms/README.md` with scope, pipeline, shared
   conventions, and production/experimental labels.
2. Add the twelve missing chapters listed above, reusing existing chapters by
   link rather than duplicating them.
3. Cross-check every implementation pointer against the current tree.
4. Check relative links, Markdown structure, whitespace, and Git diff.
5. Record completion and remaining scientific review below.

## Validation Criteria

- Every linked local Markdown file exists.
- Every named implementation file or directory exists in the current checkout.
- No chapter presents a design-only path as production behavior.
- No chapter describes probabilistic pre-alignment as full soft EM.
- No chapter conflates NU support, NU evidence, reference masking, or FSC masks.
- Existing unrelated worktree changes remain untouched.
- `git diff --check` passes for the documentation changes.
- Compilation and runtime tests are not required for documentation-only edits
  and remain with the user under repository policy.

## Completion Record

- Added `doc/algorithms/README.md` and the twelve planned chapters.
- Preserved and indexed the existing motion-correction, continuous in-plane,
  and NU-evidence chapters.
- Updated the existing NU-evidence chapter from the retired global noise scale
  to the production radial MAD whitening profile.
- Verified every relative Markdown link and every named implementation or
  policy target against the current tree.
- Verified balanced code fences, blank lines after headings, and absence of
  trailing whitespace across all 15 indexed algorithm chapters.
- `git diff --check` passed for tracked changes; an explicit scan covered the
  new untracked Markdown files, which `git diff --check` does not inspect.
- Compilation and runtime tests were not run, in accordance with repository
  policy for documentation-only work.

Remaining review is scientific/editorial: domain owners should confirm that
the selected set is the desired public definition of “critical” and that the
chapter depth is appropriate for the intended audience. No known broken links,
missing source anchors, or production/design boundary errors remain.
