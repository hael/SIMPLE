# Obsfield Refnode Reconstruction Design Note

## Status

This is a deferred design note, not current implementation policy.

The `polar=obsfield` refnode idea is attractive because it could avoid the dense
Cartesian observation-field gather while still avoiding the old all-particle by
all-reference common-line restore. The recent prototype showed that the basic
pieces can be made to run, but it also exposed two hard requirements:

1. The refnode path must be fed reconstruction-space dimensions, not
   matching-space dimensions.
2. The interpolation/data-point growth with `nspace_next` and `pftsz_next` may
   be the limiting feasibility problem.

Until both issues are explicitly designed around, this path should remain
paused. The current Cartesian `polar=obsfield` path and the `polar=no` path are
the safer baselines.

## Core Requirement

The refnode graph is not a matching object. It is a reconstruction-output object.
Therefore, when refine3D promotes the final assembly grid, every refnode
reconstruction component must be built with the promoted reconstruction
dimensions:

```text
nspace_rec = assembly output reference space
pftsz_rec  = assembly output polar section size
```

In practice this means:

- use `nspace_next` whenever assembly is promoted to `nspace_next`
- use `pftsz_next` whenever reconstruction output is promoted to `pftsz_next`
- build the refnode orientation grid from the reconstruction eulspace, not the
  matching eulspace
- write worker partial sums with headers matching reconstruction dimensions
- assemble only partial sums with exactly matching `pftsz`, shell range, and
  `nrefs`

A design that inserts particles into refnodes built for current matching
`nspace/pftsz`, then asks assembly to read them as `nspace_next/pftsz_next`, is
mode mixing. It will either fail with header mismatches or silently produce
references for the wrong discrete model.

The existing parameter helper only makes the `nspace` promotion explicit. A
production refnode design would also need a matching helper for `pftsz`, for
example:

```text
assembly_ref_nspace()
assembly_ref_pftsz()
```

or one combined dimension object used consistently by worker, assembly, and
partial-sum I/O code.

## Matching vs Reconstruction Split

The implementation should keep two concepts separate:

```text
matching polar object:
  nspace, pftsz, search k-range
  used for scoring and pose/state updates

reconstruction polar object:
  nspace_rec, pftsz_rec, interpolation k-range
  used for particle insertion, refnode partial sums, FSC, and restoration
```

Trying to make one `polarft_calc` instance serve both purposes is fragile when
`nspace_next` or `pftsz_next` differs from the matching setup. The clean design
is either two explicit objects or one object with explicit reconstruction
substate and headers that cannot be confused with matching arrays.

The reconstruction particle input also has to follow the `polar=no` semantics:
raw reconstruction-prepped particle Fourier planes, final shifts, CTF factors,
and ML/noise weighting. It should not consume images that were already modified
for matching-only alignment convenience.

## Refnode Estimator

The exact-node version should define nodes directly from the discrete polar
reference model:

```text
node(k, projection, irot)
  = direction( R_projection * [polar_x(irot,k), polar_y(irot,k), 0] )
```

Particle Fourier samples are then splatted into nearby nodes on the same shell:

```text
N(node) += W(sample, node) * particle_weight * CTF * sample
D(node) += W(sample, node) * particle_weight * CTF^2
```

The numerator and denominator must receive the same interpolation weights.
Regularization, FSC estimation, low-resolution handoff, and even/odd semantics
should then follow the existing reconstruction path as closely as possible.

A nearest-node insert is not a credible final estimator. It is useful only as a
debugging lower bound. The viable version needs a compact windowed splat whose
kernel and support can be compared against the Cartesian gridding/gather
operator.

## Scaling Risk

The main open question is feasibility.

For an exact refnode graph, the number of shell nodes scales roughly like:

```text
nk * pftsz_rec * nspace_rec
```

or half that for the internal no-mirror convention, depending on where the
mapping is stored. Either way, promotion to both `nspace_next` and `pftsz_next`
increases the node count directly.

The particle insertion cost scales like:

```text
particle_samples * neighbor_nodes_per_sample * symmetry_factor
```

A dense final reference grid therefore creates two pressures:

- more nodes and larger numerator/denominator arrays
- more candidate nodes per sample unless the support is aggressively reduced

The prototype's diagnostic scale already reached millions of shell nodes for
moderate `nspace` and `pftsz`. That is manageable for a graph direction cache,
but the true cost includes even/odd numerator and denominator arrays, worker
partial I/O, assembly reduction, and neighbor lookup overhead. This must be
benchmarked before treating the approach as production viable.

## Required Performance Design

A future implementation should not scan shell nodes during particle insertion.
It needs a bounded neighbor lookup per shell, such as:

- angular hash bins with shell-specific support
- a prebuilt spherical cell list
- compressed reference nodes with an explicit node-to-reference map
- adaptive shell budgets that cap candidates while preserving reconstruction
  quality

The exact refnode graph is the best first correctness target, but it may not be
the best production representation. A compressed graph should be considered only
after exact-node tests establish the scientific target.

Useful feasibility metrics:

- per-particle insertion time at final `nspace_next/pftsz_next`
- memory per state and per half
- partial-sum file size and assembly reduction time
- refnode vs Cartesian obsfield `corr`, `rel_l2`, and `scaled_rel_l2`
- final FSC/map quality against `polar=no`

## Refactoring Tasks Before Revival

1. Add explicit reconstruction polar dimension helpers.

   `nspace_next` and `pftsz_next` need to be promoted together through one
   reconstruction-dimension contract.

2. Split matching and reconstruction polar setup.

   Matching should keep its current dimensions and search range. Refnode
   reconstruction should allocate and write only reconstruction-sized partial
   sums.

3. Make refnode partial-sum headers strict.

   Assembly must reject stale or mixed files when `nrefs`, `pftsz`, or shell
   bounds differ.

4. Build exact operator tests.

   Start with synthetic single-particle/single-shell cases and compare against
   Cartesian obsfield and `polar=no` gridding behavior before running full
   refinement.

5. Benchmark exact-node insertion.

   Measure final promoted dimensions before optimizing. If exact-node insertion
   is already too slow, move directly to a compressed shell graph design.

6. Preserve the reconstruction statistical model.

   FSC, ML regularization, low-resolution merged insertion, trailing
   reconstruction, CTF weighting, symmetry, and final shifts must match the
   `polar=no` reconstruction semantics.

## Open Questions

- Is exact-node insertion feasible at realistic `nspace_next/pftsz_next`, or is
  compression mandatory?
- Should the splat kernel be a Cartesian KB-equivalent footprint, a shell-local
  angular KB kernel, or a calibrated approximation?
- How should support shrink with shell radius and reference density without
  creating holes?
- Can worker partial I/O remain polar-reference-shaped, or does the compressed
  graph require a new reduction format?
- Does the quality gain over Cartesian `polar=obsfield` justify the added
  reconstruction-specific machinery?

## Current Recommendation

Do not revive the refnode path as a small patch. It needs to be treated as a
separate reconstruction-dimension refactor with an explicit performance gate.

The minimum credible prototype should use `nspace_next` and `pftsz_next` for the
reconstruction graph and partial sums, keep matching dimensions separate, and
report both quality and insertion/assembly cost before any production wiring.
