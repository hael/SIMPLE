# Observation-backed polar restore design note

## Status
This is a design note, not an implementation policy. A first testable prototype now
exists in `src/main/pftc/simple_fgrid_obsfield.f90` and is wired as
`polar=obsfield`.

The existing `polar=direct` class-average path should remain unchanged. The class-average
case is currently fast and robust, and it exercises a different regime: few input
observations, no particle CTF variation, and much less pressure from particle-count
scaling.

This note is about a possible particle-path replacement for the direct common-line
restore formulation. The current `polar=obsfield` prototype now uses a dense
volume-like Fourier grid. The strategy is deliberately simple:

```text
particle Fourier samples
  -> nearest-cell accumulation on an expanded 3D Fourier grid
  -> KB interpolation of polar central sections from those raw grid sums
```

There is no insertion-side interpolation, no common-line redistribution, and no
intermediate reconstructed volume. The accumulated grid values are the object that is
sampled to make polar references.

Symmetry follows the Cartesian reconstructor convention: each particle plane is
inserted once for each symmetry-related orientation. Extraction does not apply
additional symmetry operators; it interpolates from the already symmetry-expanded
observation field.

## Motivation
The current `polar=direct` particle path avoids reprojection from a Cartesian volume, but
it still scales poorly with particle batch size. In the common-line formulation, each
particle contributes to every projection direction except the one it is assigned to.
That makes the dominant restore cost scale roughly like:

```text
nparticles_in_part * nrefs * radial_samples * symmetry_factor
```

By contrast, the `polar=no` Cartesian path has a different split:

```text
particle observations
  -> splat/interpolate once into a Cartesian Fourier volume
  -> interpolate all requested polar central sections out of that volume
```

The Cartesian route pays for a volume interpolation step, but it separates particle
accumulation from reference extraction. Once the volume has been built, extracting more
reference sections does not revisit every particle-reference pair.

The direct polar path removes the Cartesian volume, but at the cost of coupling
particle count and reference count directly. That is likely why it behaves beautifully
for around 100 class averages but becomes much less attractive for particle batches
around 1024 or larger.

## Hypothesis
We may want an intermediate representation that keeps the useful part of the Cartesian
path's scaling without averaging observations into a dense volume too early.

The idea is to accumulate particle observations once into a geometry-aware Fourier
grid field, then extract polar central sections directly from that field. This would
replace:

```text
particle -> all common-line target references
```

with:

```text
particle -> observation field -> requested polar reference sections
```

This is conceptually similar to the Cartesian volume route. The current prototype is
deliberately dense so that its inner loops look like ordinary volume work rather than
hash-table work. It is a numerator/density field over the same expanded Fourier grid
points that the Cartesian gridding path would have touched, preserving the useful
accumulate-once/query-many scaling without carrying every raw observation and without
building a real-space reconstruction.

## Representation
The observation field should be a gridded numerator/density representation, not a raw
observation list.

The current fast prototype inserts each mapped particle Fourier sample into its nearest
native 3D grid cell. The single-cell update is weighted by a center-normalized KB
factor computed from the sample's sub-cell offset, but it does not splat into
neighboring cells. Polar central sections are then gathered from the field. Each grid
cell preserves:

- The grid address.
- The accumulated complex numerator.
- The accumulated CTF2/density denominator.
- The accumulated nearest-cell KB insertion weight, which also acts as occupancy.
- The even/odd ownership through separate even and odd fields.

A useful mental model is:

```text
dense field:
  grid cell -> raw accumulated numerator and raw accumulated density
```

That is the core of the method. We first assemble the experimental Fourier evidence on
the expanded 3D grid with one-cell weighted nearest-cell deposition. Afterwards, polar
references are obtained by interpolating from those accumulated grid-point values.

This gives up the sparse memory target for the first serious performance test. The
point is to answer whether the observation-field numerical model is viable when the
implementation does not spend most of its time on sparse hash lookups.

## Extraction Model
For a requested polar Fourier sample at coordinate `q`, extraction inspects nearby
dense grid cells and applies the same KB gather logic used when sampling a gridded
Fourier object:

```text
num(q) = sum_g K(q - g) * num(g)
den(q) = sum_g K(q - g) * den(g)

value(q) = num(q) / regularized_den(q)
```

The gather must be evaluated in the same coordinate system as the deposited grid
cells. In the current prototype, insertion assigns samples to integer native Fourier
grid cells, so `q`, `g`, the KB stencil origin, and the `apod_mat_3d_fast` location
argument are all in native grid units. The padded plane is only the source from which
particle samples are read; it is not the coordinate system of the observation field.

Only occupied grid cells participate in the sums. Unoccupied cells are unknown
Fourier evidence, not valid zero-valued measurements, so they are masked out during
the KB gather. This is especially important for the nearest-cell prototype because
most dense grid cells remain untouched by a particle part. A cell is considered
occupied when its accumulated insertion weight `grid_w` is positive.

The exact form of `regularized_den(q)` must be compatible with the existing ML
regularization model. The important distinction is that regularization is applied to
the extracted query-local denominator, not to a shell-averaged common-line count that
has already lost the local observation geometry.

This deliberately has two stages, like `polar=no`: nearest-cell accumulation during
particle insertion and KB gather during polar extraction.

The insertion loop is parallelized with two native-`h` colors. All Fourier samples for
one `h` value are handled by one thread, and h values in the same color are separated
by two native grid units before and after rotation, so nearest-cell destination writes
cannot collide. This removes OpenMP atomics from the hot insertion loop while avoiding
the wider h-coloring needed by full KB splats. Observation/rejection counters are
reduced across threads. There is no finalization or secondary grid build step;
extraction reads directly from the accumulated dense arrays.

The nearest-cell KB factor is applied once, during insertion, to both the complex
numerator and CTF2 density. Extraction uses `grid_w` only to decide whether a cell has
observations; multiplying by `grid_w` again during gather would double-apply the
insertion confidence. `grid_w` is still retained as a useful path toward later
centroid or normalization variants if those become worth testing.

Insertion also follows the reconstructor's nonredundant first-axis convention: samples
whose nearest 3D cell lies fully on the redundant negative first-axis side are skipped
during insertion and recovered by Friedel conjugation during polar extraction.

Symmetry is handled at the same stage. `insert_plane_oversamp` builds one rotation
matrix for the assigned orientation and one for each `sym%apply` result, then deposits
the particle's Fourier plane through every symmetry-related matrix. Thus the dense
field contains virtual symmetry-related observations before any polar section is
queried.

## Reconstructor Coordinate Convention
The Cartesian reconstructor keeps its established insertion-side KB splat path, but
`reconstructor%insert_plane_oversamp` now evaluates the destination window and KB
weights in native destination-grid coordinates. This matches the expanded Fourier
volume lattice that receives the update. The padded Fourier plane remains only the
source from which samples are read with `hp=h*pf` and `kp=k*pf`.

This convention matches the coordinate system used by `fgrid_obsfield` extraction and
avoids weighting native destination cells with a padded-coordinate stencil.

## Distributed Execution Model
The design should fit the existing distributed part model. No worker needs to hold all
particles.

Each worker part can do:

```text
1. Read or receive the particles assigned to this part.
2. Build a local observation field from those updated particles only.
3. Extract numerator and denominator polar sections directly from that field.
4. Write the same style of even/odd partial polar sums used by current restore code.
```

The master can then do:

```text
1. Read all worker partial polar numerator and denominator arrays.
2. Sum them additively.
3. Run the existing polar restoration/normalization stage on the summed arrays.
```

This keeps the distributed reduction contract close to what already works. The new
data structure is local to a worker part and does not need to be serialized globally
unless future profiling shows a benefit.

## Critical Invariants
The representation must not collapse into either of the two paths we are trying to
avoid.

1. Do not use `calc_comlin_contrib`.

   That would return us to the common-line formulation and its unresolved self and
   in-plane weighting problem.

2. Keep the field as a numerator/density observation field.

   The current prototype is dense for performance testing, but it should not become a
   complete Cartesian reconstruction path with extra real-space reconstruction steps.

3. Do not special-case the assigned projection direction through a `SELFW`-style
   in-plane contribution.

   The goal is to query an observation field. If an observation is geometrically near a
   requested sample, it contributes through the same kernel and weighting model as all
   other observations.

4. Keep class averages on the current direct path by default.

   The class-average path is the case where `polar=direct` already has the right
   performance and convergence behavior. `polar=obsfield` is exposed for controlled
   testing, but it should not displace the direct default for class averages without
   benchmark evidence.

5. Preserve additive even/odd numerator and denominator semantics.

   The output of the new particle path should be reducible with the same distributed
   partial-sum logic as the current polar restore arrays.

## Expected Advantages
- Particle observations are inserted once per worker part, instead of once per
  particle-reference pair.
- Reference extraction is a regular interpolation query over nearby accumulated 3D
  grid values.
- The method removes one interpolation stage compared with the Cartesian path while
  keeping the favorable "accumulate once, query many" scaling.
- The method avoids common-line modulo-180 and self/in-plane contribution bookkeeping.
- CTF and per-particle weights enter the density model during insertion, then
  propagate through the same gather used for the numerator.

## Main Risks
- Memory pressure is higher than the sparse prototype because even and odd fields are
  dense expanded-grid numerator, denominator, and cell-weight arrays.
- If the dense prototype works, a later block-sparse/tiled representation may be needed
  to recover memory without returning to scalar hash-table overhead.
- The density and regularization model must be made geometry-aware without reproducing
  the failed geometric regularizer behavior.
- Symmetry expansion increases insertion work by `nsym`, because virtual
  symmetry-related observations are deposited into the field up front.
- Parallel insertion currently relies on two-color h scheduling for nearest-cell write
  safety. If insertion-side splatting returns, the wider `insert_plane_oversamp`
  race-avoidance logic would need to return with it.

## Possible Implementation Shape
The numerical backend likely belongs near the polar restore/pftc boundary, while the
workflow integration belongs in the 3D matcher strategy.

Prototype building block:

- `src/main/pftc/simple_fgrid_obsfield.f90` defines a part-local dense observation
  field and an `insert_plane_oversamp` fill method that assigns particle Fourier
  samples to their nearest 3D grid cells with one-cell KB weighting.
- Symmetry is handled during insertion by applying every symmetry operator to the
  particle orientation and depositing through the resulting rotation matrices.
- The current prototype stores accumulated grid-cell numerator, denominator, and
  insertion-weight values, not raw observation records.
- `polar_cavger_insert_ptcls_obsfield` is wired as the `polar=obsfield` restore mode.
  It fills the dense field, immediately extracts even/odd polar partial sums from the
  accumulated grid arrays, and adds them to the existing `pfts_even/odd` and
  `ctf2_even/odd` payloads.

Candidate ownership:

- `src/main/pftc`: observation field type, insertion, query extraction, and
  interpolation kernel use.
- `src/main/strategies/search`: particle-batch orchestration and distributed partial
  output integration.
- Existing restore I/O: reuse the current polar partial-sum payload shape if possible.

The first validation should be intentionally limited:

1. C1 symmetry only.
2. Weighted nearest-cell dense-grid accumulation without multi-cell insertion splatting.
3. Existing interpolation kernel and CTF/noise conventions.
4. Existing even/odd partial-sum output format.
5. Class-average code path untouched.

After that, add memory-control strategies only once numerical equivalence has been
established.

## Validation Plan
A useful validation sequence would be:

1. Tiny C1 no-CTF synthetic case.

   Compare extracted polar sections from the observation field against polar
   sections extracted from the Cartesian `polar=no` volume built from the same
   observations.

2. Add CTF/noise weighting.

   Compare numerator and denominator shell statistics against the Cartesian reference
   and the current direct-polar path.

3. Distributed additivity check.

   Split the same particles into several parts, extract partial polar arrays, sum them,
   and compare with a single-part extraction.

4. Performance scaling check.

   Measure separately:

   - observation insertion time,
   - extraction time per reference,
   - memory per worker part,
   - cost as `nparticles_in_part` grows,
   - cost as `nrefs` grows.

5. Convergence tests.

   Only after the above checks pass should this be tested in full `abinitio3D` particle
   runs.

## Decision Point
This direction is worth considering because it targets the specific scaling failure:
the current direct particle path couples particle count and reference count too tightly.
The observation field would preserve the experimental observation geometry while
recovering the `polar=no` path's useful separation between accumulation and reference
extraction.

The key design test is whether it can remain observation-backed rather than becoming a
second Cartesian volume implementation or a renamed common-line implementation.
