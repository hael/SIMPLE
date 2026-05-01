# Reference-Node Shell Polar Reference Assembly Draft

## Status

This is a future design note, not current implementation policy.

The immediate priority remains to bring the current `polar=obsfield` path up to
par with `polar=no` in reconstruction quality and workflow behavior. This note
records a possible successor design to revisit only after the current obsfield
path has been made scientifically comparable and benchmarked.

## Motivation

The useful long-term target is a polar pathway that reuses registered particle
polar representatives to generate the next iteration's polar reference central
sections without:

- reconstructing a Cartesian volume at every iteration
- scaling like the current `polar=yes` particle-reference common-line route
- depending on ad hoc line interpolation between target references
- diverging from the ML regularization and shift semantics of `polar=no`

The current `polar=obsfield` prototype is an important step because it separates
particle accumulation from reference extraction. However, it still uses a dense
Cartesian Fourier observation field as the intermediate representation. The
design below asks whether the natural intermediate object for a polar refinement
path is instead a shell-local graph whose nodes are the actual polar reference
samples required by the discrete search space.

## Core Idea

Each registered particle polar image samples great circles on 3D Fourier shells.
For radial frequency `k`, a polar particle sample

```text
P_i(k, alpha)
```

with assigned orientation `R_i` corresponds to the 3D Fourier point:

```text
q = R_i * [k cos(alpha), k sin(alpha), 0]
```

The direction of that point,

```text
qhat = q / |q|
```

is an observation of a spherical function on the shell with radius `k`.

The proposed representation is therefore not a generic latitude/longitude
spherical image. For the first implementation, the shell nodes should be exactly
the directions already requested by SIMPLE's discrete reprojection model:

```text
node_s,k,j,irot =
    direction( eulspace%get_mat(j) * [polar_x(irot,k), polar_y(irot,k), 0] )
```

where `j` is the discrete projection direction and `irot` is the in-plane polar
sample index. Each state, half, and shell owns numerator and denominator values
on these reference nodes:

```text
for each state, half, shell k, projection j, and in-plane sample irot:
    reference-node numerator   N_s,half,k(j,irot)
    reference-node denominator D_s,half,k(j,irot)
```

Particle observations are accumulated once into nearby shell nodes using a
compact angular distance kernel. This should be a windowed scatter, not a
nearest-cell assignment: each particle Fourier component may contribute to every
reference node inside the chosen shell-local support radius, with numerator and
denominator weights formed from the same kernel. The next iteration's polar
references are generated directly from the assembled reference-node shell graph.
This keeps the polar reprojection geometry discrete and explicit, instead of
creating a second spherical geometry that later has to be interpolated back onto
the reference space.

With exact reference nodes, normalization is the reprojection-model generation
step. Once the node numerators and denominators have been reduced and
regularized, the normalized node array is already the polar central-section
model consumed by matching. There is no separate shell-to-reference projection
stage:

```text
normalized node_s,half,k(j,irot) == POLAR_REFS_s,half(j,irot,k)
```

The design problem is therefore shifted to fast insertion: given an observed
particle-shell direction, identify all in-window reference nodes quickly enough
that the method remains local in the particle loop.

This graph is a 3D coupling object, not a set of independent 2D projection
averages. A particle sample contributes to nearby shell nodes across the
assembled central-section geometry. A design that only accumulates a particle
into its assigned projection bin would be a 2D class-average path and is not the
target described here.

## Geometry Choice

The observation-field shell geometry is a modeling and implementation choice.
The existing KB machinery is used here as a compact distance-weighting kernel
over nearby cells or nodes, not as a strict requirement that the intermediate
object be a Cartesian Fourier grid. The geometry should therefore be selected
for stable, efficient generation of the discrete polar reprojection model.

Prototype order:

1. **Reference-node shell graph.** Use the exact `(proj, irot)` directions for
   each shell as graph nodes. This is the least clever representation, matches
   the emitted polar references exactly, and gives the cleanest correctness
   comparisons against the Cartesian obsfield path.
2. **Compressed reference-node shell graph.** Once the exact-node graph is
   scientifically acceptable, merge duplicate or near-duplicate shell directions
   within a controlled tolerance and keep an explicit node-to-reference mapping.
   This should be treated as a performance/memory optimization, not as the first
   scientific prototype.

Generic latitude/longitude shell maps are not preferred for this path. They are
easy to implement but impose an unrelated angular grid, suffer from uneven area
and pole behavior, and require an extra interpolation back to the actual
reference samples. Equal-area shell grids such as HEALPix/Fibonacci nodes may be
useful later, but they are secondary to proving the exact reference-node graph.

Auxiliary search geometries are allowed only for neighbor identification. For
example, an equal-area hash, Fibonacci binning, or another coarse angular index
may be used to find candidate reference nodes on a shell. That auxiliary index
must not become the obsfield representation and must not define the emitted
polar model. The source of truth remains the exact `(proj, irot, k)` node array.

## Estimator

For each accepted particle sample, find the reference-node shell directions
inside the same radial shell's compact angular window and accumulate:

```text
N_s,half,k(node) += W(qhat_i, node)
                    * w_i * conj(CTF_i(k, alpha)) * P_i(k, alpha)
                    / sigma2_noise(k)

D_s,half,k(node) += W(qhat_i, node)
                    * w_i * |CTF_i(k, alpha)|^2
                    / sigma2_noise(k)
```

where:

- `s` is the assigned state
- `half` is even or odd ownership
- `w_i` is the particle/sample weight used by the matching policy
- `sigma2_noise(k)` follows the same shell range and noise policy as `polar=no`
- CTF factors and incremental shift phases are evaluated for the same registered
  sample
- `W(qhat_i, node)` is a compact shell-local distance weight, using the existing
  KB kernel shape or an explicitly documented compatible replacement
- nodes outside the support window receive zero contribution

The support radius is a tunable estimator parameter. A nearest-cell prototype is
useful only as a debugging lower bound; it is not the intended estimator because
it becomes brittle as the reference-node set becomes dense. With many projection
directions, assigning each noisy sample to only one cell can underuse the
available local angular information and make occupancy artifacts look like
signal. Windowed scatter instead lets the node graph behave like a genuine
local kernel estimator on the shell while still emitting the exact requested
polar reference geometry after normalization.

The denominator must receive exactly the same window weights as the numerator.
Any post-normalization occupancy correction should be treated as an explicit
regularization choice, not as a hidden compensation for nearest-neighbor
dilution.

After distributed or shared-memory reduction, normalize with the same
ML-regularized logic as the Cartesian branch:

```text
F_s,half,k(node) = N_s,half,k(node) / (D_s,half,k(node) + prior_s,k)
```

The prior must be derived from the compatible even/odd estimates with the same
FSC/SSNR policy used by the current reconstruction path. The reference-node
shell representation changes where observations are stored, not the statistical
model.

## Direct Polar Model Handoff

Once the node fields have been assembled and regularized, reference generation
is an array handoff:

```text
ref_s,j,half(k,irot) = F_s,half,k(node(j,irot))
```

This decouples particle count from reference count:

```text
particle samples -> reference-node shell graph -> requested polar references
```

The expensive particle loop is paid once per iteration or stage. The exact-node
prototype does not support "free" additional reference directions after
accumulation; the node graph must be built for the assembly output reference
space before insertion. If a compressed graph is used, reference generation is a
small gather through the explicit node-to-reference map.

There is no separate field-to-reprojection-model operation in the exact-node
prototype. Writing `POLAR_REFS.bin`, `POLAR_REFS_even.bin`, and
`POLAR_REFS_odd.bin` should be a layout conversion from normalized node arrays
to the established polar reference file contract.

## Fast Windowed Neighbor Identification

The graph construction cost must be amortized. For each shell, implementation
should precompute the lookup structure needed to scatter an observed particle
direction into nearby reference nodes. The particle loop should never scale as
all particles times all reference sections.

A concrete first design is:

```text
for each shell k:
    build node_dir(3,node) from all (proj,irot) samples
    build a coarse angular hash/bin index -> candidate node ids

for each particle sample qhat_i on shell k:
    bin = hash(qhat_i)
    candidates = node ids in bin and neighboring bins
    accept nodes with dot(qhat_i,node_dir) >= cos(theta_support)
    weight accepted nodes by KB-shaped distance W
```

The coarse hash is only an acceleration structure. Correctness depends on the
final dot-product or distance test against the exact reference-node directions.
The support radius and neighboring-bin expansion must be conservative enough
that no node inside the chosen KB support is missed. The desired output of this
lookup is a compact stencil of accepted node ids and weights, not a single best
node.

Implementation should avoid `acos` in the hot loop when possible. A dot-product
cutoff and a monotonic chord-distance weight,

```text
d_chord = sqrt(max(0, 2 - 2 * dot(qhat_i,node_dir)))
```

should be sufficient unless validation shows the angular metric itself is
needed.

The window size should be tied to shell radius and reference-node density. A
fixed angular window is simple but may over-smooth low-frequency shells or leave
high-frequency shells too sparse. A practical first policy is to expose a small
set of options:

- a conservative default support based on local node spacing per shell
- an integer or real support multiplier for tests
- a nearest-node/debug mode only for validation and performance bracketing

Diagnostics should report both candidate counts and accepted neighbor counts.
If accepted-neighbor counts grow with `nspace` until insertion becomes
reference-count dominated, the support policy is too broad. If too many nodes
remain weakly populated or empty at high frequency, the support policy is too
narrow or the lookup is missing in-window nodes.

## Shift Semantics

The reusable polar particle representative must be prepared from a centered
particle image. The old origin shift is therefore intentionally baked in before
mask application and polarization. This is critical because the mask and polar
Fourier representation should be built from the best currently centered particle,
not from an uncentered raw image.

The shell accumulator should then apply only the incremental shift discovered by
the current optimization step. That incremental shift is applied as a Fourier
phase factor during shell accumulation:

```text
prepared particle = image shifted by the old origin estimate, then masked and polarized
shell insertion   = prepared polar sample * exp(i * phase(incremental_shift))
```

If the orientation store records cumulative shifts, the accumulation phase must
be formed from the delta relative to the shift already baked into the prepared
particle. It must not apply the full absolute shift again.

This is a critical invariant. A shell estimator that silently double-applies
or drops incremental shifts would be scientifically worse than the Cartesian path
even if its scaling were attractive.

## Symmetry

Symmetry should be handled during observation accumulation, matching the
Cartesian reconstructor convention:

```text
for each assigned particle orientation:
    insert the observed great circle for R_i
    insert the same observed great circle for each sym * R_i
```

The normalized node field should already include the symmetry-expanded
observations. Writing polar references must not apply a second independent
symmetry expansion.

## Relationship To Existing Paths

`polar=no`:

```text
particle Fourier planes -> Cartesian Fourier volume -> polar reference sections
```

Current `polar=yes`:

```text
particle polar representatives -> common-line contributions to target refs
```

Current `polar=obsfield`:

```text
particle Fourier planes -> dense Cartesian observation field -> polar refs
```

Proposed reference-node shell assembly:

```text
particle polar representatives -> reference-node shell graph -> normalized polar refs
```

The proposed path keeps the attractive accumulate-once structure of `polar=no`
and `polar=obsfield`, but stores the object in the exact geometry where polar
reference sections naturally live. After normalization, the object is no longer
an intermediate representation; it is the emitted polar reference model.

## Bridge Design: Cartesian Accumulation With Reference-Node Cache

An intermediate design worth testing before direct particle-to-node
accumulation is:

```text
particle Fourier planes
  -> dense Cartesian observation field
  -> reference-node shell cache
  -> polar refs
```

This keeps the current `polar=obsfield` insertion hot loop unchanged. After all
particle insertions and distributed reductions are complete, assembly resamples
the restored Cartesian observation field directly onto the exact
`(proj, irot, k)` reference nodes. Polar reference writing then becomes a copy
from the node cache instead of a repeated Cartesian KB gather for every
requested central-section sample.

The scientific status of this bridge is different from the direct node
estimator above. It should first be treated as a cache of the Cartesian
obsfield estimator, not as a new estimator. That means the first validation
question is whether the cached node array reproduces current obsfield polar
references closely enough to be useful. If it does, the cache can isolate the
performance question: does paying one shell projection per assembly beat
paying the Cartesian gather for every reference section?

The cost model is:

```text
current obsfield extraction:
    nrefs * polar_samples_per_ref * Cartesian_KB_stencil

node-cache extraction:
    node_cache_samples * Cartesian_KB_stencil
    + layout conversion to POLAR_REFS*
```

This bridge is attractive when the number of requested references is large
enough that the one-time shell projection amortizes well. Its main drawbacks
are the transient memory footprint of holding both representations and the
extra sampling step introduced by Cartesian-to-node resampling.

## Distributed Execution Model

Workers should not need global particle state.

Each worker part can:

```text
1. prepare or reuse registered particle polar representatives
2. accumulate local reference-node shell numerator and denominator fields
3. write partition-local node fields, or write normalized partial polar refs
```

The assembly step can:

```text
1. reduce node numerator and denominator fields across parts
2. apply ML regularization and optional symmetry validation
3. write normalized node arrays as POLAR_REFS.bin, POLAR_REFS_even.bin, and POLAR_REFS_odd.bin
4. optionally materialize Cartesian MRC representatives at stage boundaries
```

The cleanest conceptual model is to reduce node fields first, normalize once
centrally, and then write the normalized nodes in the established polar reference
layout. A more memory-conservative variant could normalize or emit partial
reference sections per worker and reduce those, but that gives up part of the
value of a shared shell representation and should be treated as a fallback.

## Cartesian Representatives

The node shell field is not itself a Cartesian reconstruction. For diagnostics and
abinitio stage outputs, a Cartesian representative can be generated at stage
boundaries by gridding the regularized node field into a Cartesian Fourier
volume and transforming/postprocessing through the usual volume-domain policy.

This diagnostic volume should not become the source of truth for the next polar
matching iteration. Its purpose is inspection, FSC reporting, stage bookkeeping,
and comparison with the current `polar=no` workflow.

## Expected Advantages

- Particle accumulation scales with particle polar samples, not
  particle-reference pairs.
- Polar references are mutually consistent because they are sections through one
  assembled shell graph per state and half.
- Normalization directly produces the polar reference model; no extra
  shell-to-reference interpolation or reprojection step is required in the
  exact-node prototype.
- The method avoids direct common-line redistribution between target references.
- It creates a natural home for reusing prepared polar particle representatives.
- It may remove the poor scaling behavior of `polar=yes` while preserving a
  genuinely polar intermediate representation.
- The first prototype emits references from the exact same discrete geometry
  used by matching, making compare-mode failures easier to localize.

## Risks And Open Questions

- Reference-node shell sampling must be accurate enough at high frequency and
  near sparse angular coverage regions.
- The first graph may be large: `nspace * pftsz * nshells` nodes per state/half,
  plus neighbor stencils. Memory layout and compression must be measured before
  production use.
- A particle sample must scatter to a bounded local stencil of shell nodes, not
  just to a single nearest cell except in explicit debug mode.
  Accidentally falling back to all-particles-times-all-references would recreate
  the poor scaling this design is meant to avoid.
- The neighbor acceleration structure must be conservative. Missed in-support
  nodes would create hard-to-see angular holes even if the normalized array shape
  looks correct.
- Candidate lists must stay small enough for the particle loop. Per-shell
  candidate and accepted-neighbor statistics are required before benchmarking
  reconstruction quality.
- ML regularization must remain comparable to `polar=no`; shell occupancy must
  not become an unprincipled extra normalization.
- Shift handling must be audited before any benchmark is meaningful.
- CTF zeros and uneven angular coverage may require careful denominator floors.
- Memory layout and cache behavior may determine whether the theoretical scaling
  advantage survives in practice.
- Symmetry expansion must avoid double counting and preserve even/odd separation.
- Compression of near-duplicate nodes must preserve an explicit mapping back to
  every requested reference sample and must be introduced only after the exact
  reference-node graph validates.

## Validation Plan

A serious prototype should be judged against `polar=no` and current
`polar=obsfield` on the same datasets:

1. fixed-orientation synthetic data with known shifts and CTFs
2. asymmetric small membrane-protein-like targets
3. larger symmetric targets where current paths already behave well
4. abinitio3D stage transitions where Cartesian representatives are inspected
5. shared-memory and distributed runs with identical scientific outputs

The acceptance bar should be:

```text
same or better reconstruction quality than polar=no
and
better scaling than current polar=yes
and
no new shared-memory/distributed workflow divergence
```

If it cannot meet those conditions, it should remain a research note rather
than another production polar path.

Prototype diagnostics should include:

- whole-array max/rms/relative-rms against Cartesian obsfield polar extraction
- per-shell max/rms/relative-rms
- per-shell node count, populated-node count, and populated fraction
- neighbor-stencil min/mean/max by shell
- coarse-bin candidate count min/mean/max by shell
- accepted-neighbor count min/mean/max by shell
- missed-support sentinel tests for the neighbor index
- denominator statistics by shell and half
- number of emitted reference samples that map to missing or under-supported
  nodes
- for the compressed graph, compression ratio and max angular merge distance

## Non-Goals

- Do not use this note to delay fixing the current `polar=obsfield` correctness
  and parity issues.
- Do not introduce another shared-memory-only reference preparation path.
- Do not move assembled-reference construction into the matcher.
- Do not use Cartesian diagnostic volumes as the hidden source of truth for
  polar matching.
- Do not weaken the existing `POLAR_REFS*` handoff contract without a separate
  explicit policy change.
