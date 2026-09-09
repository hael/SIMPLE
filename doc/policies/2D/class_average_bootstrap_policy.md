# Class-Average Bootstrap Policy

This document defines the first implementation policy for stochastic
class-average expansion before `abinitio3D_cavgs`. It should be read alongside
[abinitio2D_policy.md](abinitio2D_policy.md) and
[abinitio3D_cavgs_policy.md](abinitio3D_cavgs_policy.md).

The workflow creates one `abinitio3D_cavgs`-compatible bootstrap package from
an existing class-average project. It is an offline class-average assembly
stage, not a new 2D classification, not a particle-realignment step, and not a
source of final particle orientations.

In plain language, the workflow:

1. copies the accepted original class averages into the start of the output
   stack, preserving their original class-index order;
2. derives one anchor-particle count from the median parent-class population
   and public `frac_best`;
3. selects the best-ranked particles within each parent class as that parent's
   shared anchor;
4. uses the public oversampling factor to create the same number of stochastic
   children for every accepted parent;
5. assembles each child from the same parent anchor plus a different random
   subset of the parent's remaining particles;
6. writes merged, even, and odd bootstrap stacks plus lineage files.

## 1. Public Contract

Public controls:

```text
osmpl_fac        integer oversampling factor; default 2
frac_best        default 0.5; applied to the median parent-class population
```

If `K` is the number of accepted original classes, the final stack contains
`K * osmpl_fac` class averages.

Internal first-version policy values:

```text
original block      fixed yes
inactive parents    excluded
ensembles per run   1
```

The command generates one bootstrap ensemble per execution. Users who want
multiple stochastic ensembles should run the workflow multiple times. The first
public interface should not expose `seed`, `nensembles`, `anchor_target`,
`ncls_target`, per-half minima, copy-fill behavior, child-specific scale
factors, or cross-parent combination.

## 2. Terms

Parent class

: One original 2D class average from the input project.

Accepted parent

: A parent class with `cls2D%state > 0`. Only accepted parents are copied or
  expanded.

Original block

: The leading part of the output stacks. It contains copied accepted parent
  class averages in original parent class-index order. This makes mapping back
  to the original class-average set direct.

Bootstrap child

: A synthetic class average assembled from particles belonging to exactly one
  parent class.

Anchor subset

: The best-ranked particles in a parent class that are copied into every child
  of that parent. The anchor keeps all children tied to the same parent
  projection signal.

Variable subset

: A random subset of the remaining non-anchor particles assigned to one child.
  This is what makes children from the same parent different.

Objective-rank ordering

: The existing best-to-worst ordering of particles inside a 2D class according
  to the 2D objective-function value. The bootstrap command consumes this
  ordering; it does not recompute alignments or scores.

Parity-complete anchor

: An anchor subset containing both even and odd contributors. Expanded parents
  must have a parity-complete anchor so every child can produce merged, even,
  and odd averages.

Child capacity

: The maximum number of children a parent can support. In the first policy this
  is the number of non-anchor particles, because each child must receive at
  least one variable particle.

Oversampling factor

: The requested multiplier for the accepted parent class count. With
  `osmpl_fac=2`, the output contains the original block plus one stochastic
  child per accepted parent, for a total of `2*K` class averages.

## 3. Required Inputs

The command requires:

- a project with `cls2D`;
- merged, even, and odd class-average stacks matching `cls2D`;
- source-particle membership for each parent class;
- frozen 2D class assignment, in-plane, and shift metadata;
- objective-function values, or an equivalent best-to-worst source-particle
  ordering, within each parent class.

Only parent classes with `cls2D%state > 0` are accepted. Missing merged/even/odd
stacks, missing membership, missing objective ranking, or stack-size mismatch
are hard errors.

## 4. Ownership

This workflow belongs near explicit class-average assembly:
`simple_commanders_mkcavgs.f90`, the classaverager modules, or a new commander
with the same ownership style.

It must not live in `simple_strategy2D_matcher.f90`,
`simple_cluster2D_strategy.f90`, `simple_oris_sampling.f90`, or `volassemble`.
It may perform explicit offline particle-stack reads. The online matcher
single-read restoration contract is unaffected.

## 5. Formula Symbols

For accepted parent class `C_i`, the formulas below use:

```text
P_i      source particles assigned to C_i
N_i      number of source particles in P_i
A_i      anchor subset copied into every stochastic child of C_i
V_i.j    variable subset assigned to child j
Y_i      number of stochastic children allocated to C_i
Y_child  fixed number of stochastic children per accepted parent
Y_cap_i  maximum allowed Y_i
```

Each bootstrap child has exactly one parent. Particles or images from different
parents are never combined.

## 6. Output Ordering

The output package contains:

```text
cavgs_bootstrap_001.mrc
cavgs_bootstrap_001_even.mrc
cavgs_bootstrap_001_odd.mrc
cavgs_bootstrap_001_manifest.txt
cavgs_bootstrap_001_membership.txt
```

Merged, even, and odd stacks must have identical ordering:

1. copied parent class averages in original parent class-index order;
2. stochastic children, grouped by parent class in original parent class-index
   order and then by child number.

The leading original block makes stack entries `1:K` map directly to the
accepted original parent classes, where `K` is the number of accepted parents.
Bootstrap children are appended after that block and must never shift copied
parent entries.

The original project and original class-average stacks are not modified.

## 7. Oversampling Factor

The requested `osmpl_fac` directly controls the final class count and preserves
the original class projection-direction distribution by expanding every
accepted parent equally.

```text
K = number of accepted parent classes, equivalent to ncls_orig for this workflow
Y_child = osmpl_fac - 1
N_out = K * osmpl_fac
```

`osmpl_fac` must be an integer greater than or equal to 1. If `osmpl_fac=1`,
the output is the original block only, with the manifest recording that no
stochastic children were requested.

## 8. Anchor Selection

The run-level anchor target is derived from the median accepted parent
population:

```text
N_med           = median_i(N_i)
N_anchor_target = max(2, nint(frac_best * N_med))
```

Using the median avoids letting one very large or very small class determine
the anchor size for the whole run.

Each parent receives a realized anchor size:

```text
N_anchor_i = min(N_anchor_target, max(0, N_i - Y_child))
frac_best_i = real(N_anchor_i) / real(N_i)
```

The `max(0, N_i - Y_child)` bound leaves at least one non-anchor source
particle for each stochastic child. Large classes therefore get a smaller
realized anchor fraction and can be diversified more. Small classes get a
larger realized anchor fraction, which protects their signal while preserving
the minimum variable support needed for `osmpl_fac`.

`A_i` is selected from the best-ranked source particles in `C_i` according to
the existing 2D objective-function ordering. The workflow consumes this
ordering as fixed input; it does not recompute alignments or objective values.

For a parent to receive stochastic children, `A_i` must contain both even and
odd contributors when both are present in the parent. If a parity-complete
anchor cannot be formed, that parent has zero stochastic child capacity and is
copied only in the original block.

The manifest must record `N_med`, `N_anchor_target`, `N_anchor_i`, and
`frac_best_i`.

## 9. Child Count Feasibility

After anchor selection:

```text
R_i     = N_i - |A_i|
Y_cap_i = R_i if A_i is parity-complete, otherwise 0
```

This means each allocated child receives at least one non-anchor source
particle, and every child has even/odd support through the shared anchor.

Every accepted parent must support the same child count:

```text
Y_i = Y_child
```

If any accepted parent has `Y_cap_i < Y_child`, the command stops and reports
the limiting parent classes. Exact-copy fill-in is not part of the first
workflow. This keeps the output distribution proportional to the original
accepted class-average distribution instead of preferentially expanding rare or
high-population parents.

## 10. Child Assembly

For each expanded parent:

1. copy `A_i` into every child;
2. randomly partition the non-anchor particles into `Y_i` balanced variable
   subsets without replacement;
3. assemble child `C_i.j` from `A_i union V_i.j` using existing 2D
   assignment/alignment metadata;
4. assemble merged, even, and odd child averages independently but in identical
   stack order.

Under the default partition rule, every non-anchor source particle appears in
exactly one child for that parent, and every anchor source particle appears in
every child.

The command must not perform new particle alignment, class reassignment,
child-to-parent orientation consensus, or mapping of synthetic child
orientations back to original particles.

## 11. Metadata

The manifest records one row per output class-average entry, including copied
parents and stochastic children:

```text
output_class_id
entry_kind
original_block_index
child_class_id
parent_class_id
parent_state
parent_count
child_count
child_even_count
child_odd_count
osmpl_fac
children_per_parent
children_cap_for_parent
median_parent_count
anchor_target
frac_best_input
frac_best_realized
anchor_count
anchor_even_count
anchor_odd_count
variable_count
variable_even_count
variable_odd_count
run_seed
parent_seed
status
```

The membership file records source membership:

```text
child_class_id
parent_class_id
source_particle_id
source_stack_id
source_stack_index
parity
role
```

where `role` is `anchor` or `variable`.

All stochastic choices must be traceable. The first public interface does not
expose a seed; it records an internal run seed and deterministic parent seeds.
Exact reproduction should use the recorded manifest and membership file.

## 12. Handoff

The bootstrap package is intended for `abinitio3D_cavgs` seed generation.
Synthetic child orientations remain attached only to the bootstrap package. If
`abinitio3D_cavgs` needs project mapping internally, the bootstrap route should
use a temporary bootstrap project and avoid writing synthetic orientation state
into the original project.

Downstream continuation should use the original data, not synthetic child
assignments.

## 13. Failure Conditions

The command stops when:

- `cls2D` is absent or has no accepted parent classes;
- `osmpl_fac` is missing or less than 1;
- merged, even, or odd stacks are missing;
- stack sizes, dimensions, or sampling are inconsistent;
- source membership or objective ordering is unavailable;
- a child would have an empty even or odd half;
- any accepted parent cannot support `osmpl_fac - 1` stochastic children;
- manifest or membership output cannot be written.

## 14. Review Checks

Before implementation, confirm:

- public controls are only `osmpl_fac` and `frac_best`;
- one command execution produces exactly one ensemble;
- original parent copies occupy the first `K` stack entries in original index
  order;
- every accepted parent receives exactly `osmpl_fac - 1` stochastic children;
- every child has exactly one parent and no cross-parent combination occurs;
- `N_anchor_target` is derived from median population and public `frac_best`;
- anchors are selected from objective-rank best particles;
- expanded parents have parity-complete anchors;
- every expanded parent satisfies `Y_cap_i >= osmpl_fac - 1`;
- no new particle alignment, reassignment, orientation consensus, or synthetic
  mapping back to the original project is introduced;
- merged, even, and odd outputs are complete and identically ordered;
- lineage is sufficient to audit every stochastic choice.

Useful validation diagnostics are child-count distributions, even/odd support,
child-vs-parent cavg correlation, cavg-vs-reprojection correlation after
`abinitio3D_cavgs`, and sensitivity to `osmpl_fac` and `frac_best`.
