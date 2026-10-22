# Class-Average Bootstrap Policy

This document defines the stochastic class-average expansion before
`solve3D_cavgs` (`bootstrap_cavgs`). It should be read alongside
[solve2D_policy.md](solve2D_policy.md) and
[solve3D_cavgs_policy.md](../3D/solve3D_cavgs_policy.md).

The workflow creates one `solve3D_cavgs`-compatible bootstrap package from
an existing class-average project. It builds a bootstrap project from the
particle memberships of the usable classes and runs `make_cavgs` on it. It is
an offline class-average assembly stage, not a new 2D classification, not a
particle-realignment step, and not a source of final particle orientations.

In plain language, the workflow:

1. keeps every usable original class as the leading block of the output
   stacks, in original class-index order; `make_cavgs` assembles the block from
   the class's particles;
2. derives one anchor-particle count from the median parent-class population
   and `frac_best`;
3. selects the best-ranked particles within each parent class as that parent's
   shared anchor;
4. uses the oversampling factor to create the same number of stochastic
   children for every usable parent;
5. assembles each child from the same parent anchor plus a different random
   subset of the parent's remaining particles;
6. writes merged, even, and odd bootstrap stacks, the bootstrap project, and
   the manifest and membership files.

## 1. Public Contract

Controls:

```text
osmpl_fac        integer oversampling factor, at least 1; default 2
frac_best        in (0,1]; default 0.5; applied to the median parent-class population
refs             name of the merged output stack; default cavgs_bootstrap_001.mrc
```

If `K` is the number of usable original classes (Sections 2 and 9), the final
stack contains `K * osmpl_fac` class averages.

Internal policy values:

```text
original block      fixed yes
inactive parents    excluded
ensembles per run   1
```

The command generates one bootstrap ensemble per execution. Users who want
multiple stochastic ensembles run the workflow multiple times. The command has
no `seed`, `nensembles`, `anchor_target`, `ncls_target`, per-half minima,
copy-fill behavior, child-specific scale factors, or cross-parent combination.

## 2. Terms

Parent class

: One original 2D class average from the input project.

Accepted parent

: A parent class whose `cls2D` state is not 0 (a class without a state is
  accepted) and whose population is not zero.

Usable parent

: An accepted parent that passes the checks of Sections 8 and 9. Only usable
  parents are copied or expanded; the other accepted parents are skipped and
  recorded.

Original block

: The leading part of the output stacks. It contains one average per usable
  parent, in original parent class-index order, assembled from all of the
  parent's particles. This makes mapping back to the original class-average
  set direct.

Bootstrap child

: A synthetic class average assembled from particles belonging to exactly one
  parent class.

Anchor subset

: The best-ranked particles in a parent class that are copied into every child
  of that parent. The anchor keeps all children tied to the same parent
  projection signal.

Variable subset

: A random subset of the remaining non-anchor particles assigned to one child.
  This is what makes children from the same parent different. Its particles
  have the role `member` in the files.

Objective-rank ordering

: The best-to-worst ordering of particles inside a 2D class by the particle
  `corr`. The bootstrap command consumes this ordering; it does not recompute
  alignments or scores.

Parity-complete anchor

: An anchor subset containing both even and odd contributors. Expanded parents
  have a parity-complete anchor so every child can produce merged, even, and
  odd averages.

Child capacity

: The maximum number of children a parent can support. It is the number of
  non-anchor particles, because each child must receive at least one variable
  particle.

Oversampling factor

: The requested multiplier for the usable parent class count. With
  `osmpl_fac=2`, the output contains the original block plus one stochastic
  child per usable parent, for a total of `2*K` class averages.

## 3. Required Inputs

The command requires a project with non-empty `ptcl2D`, `cls2D` and `stk`
segments. For each parent it reads:

- the source-particle membership (the `class` of the `ptcl2D` rows);
- the objective value (`corr`) of each member, which gives the ranking;
- the even/odd label (0 or 1) of each member;
- the in-plane and shift parameters of each member, which are carried into the
  bootstrap project for `make_cavgs`.

It reads no class-average stack. A parent that fails these requirements is
skipped (Section 9). The command stops when a segment is empty or no parent is
usable (Section 13).

## 4. Ownership

This workflow belongs near explicit class-average assembly:
`commander_bootstrap_cavgs` in `simple_commanders_mkcavgs.f90`.

It must not live in `simple_strategy2D_matcher.f90`,
`simple_refine2D_strategy.f90`, `simple_oris_sampling.f90`, or `volassemble`.
Its particle reads are the explicit offline reads of `make_cavgs`. The online
matcher single-read restoration contract is unaffected.

## 5. Formula Symbols

For usable parent class `C_i`, the formulas below use:

```text
P_i      source particles assigned to C_i
N_i      number of source particles in P_i
A_i      anchor subset copied into every stochastic child of C_i
V_i.j    variable subset assigned to child j
Y_i      number of stochastic children allocated to C_i
Y_child  fixed number of stochastic children per usable parent
Y_cap_i  maximum allowed Y_i
```

Each bootstrap child has exactly one parent. Particles or images from different
parents are never combined.

## 6. Output Ordering

The output package contains, for the default `refs`:

```text
cavgs_bootstrap_001.mrc
cavgs_bootstrap_001_even.mrc
cavgs_bootstrap_001_odd.mrc
cavgs_bootstrap_001_manifest.txt
cavgs_bootstrap_001_membership.txt
```

and the bootstrap project, written to the command's project file. The stack
and file names follow `refs`.

Merged, even, and odd stacks have identical ordering:

1. the original block, one entry per usable parent in original parent
   class-index order;
2. stochastic children, grouped by parent class in original parent class-index
   order and then by child number.

The leading original block makes stack entries `1:K` map directly to the
usable parents, where `K` is the number of usable parents; the manifest lists
the parent class of every entry. Bootstrap children are appended after that
block and never shift original entries.

The bootstrap project goes to the project file of the execution: with the
default `mkdir=yes` that is the copy in the execution directory, and the
original project and original class-average stacks are not modified; with
`mkdir=no` the input project file is overwritten.

## 7. Oversampling Factor

The requested `osmpl_fac` directly controls the final class count and preserves
the original class projection-direction distribution by expanding every
usable parent equally.

```text
K = number of usable parent classes
Y_child = osmpl_fac - 1
N_out = K * osmpl_fac
```

`osmpl_fac` must be an integer greater than or equal to 1. If `osmpl_fac=1`,
the output is the original block only, the manifest has no `child` rows, and
no anchor checks apply.

## 8. Anchor Selection

The run-level anchor target is derived from the median population of the
accepted parents:

```text
N_med           = median_i(N_i), rounded to an integer, at least 1
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
the minimum variable support needed for `osmpl_fac`. A parent whose realized
anchor is smaller than 2 is skipped (`too_few_particles_for_osmpl_fac`).

`A_i` is selected from the best-ranked source particles in `C_i` according to
the objective-rank ordering. The workflow consumes this ordering as fixed
input; it does not recompute alignments or objective values.

For a parent to receive stochastic children, `A_i` must contain both even and
odd contributors. When the best-ranked `N_anchor_i` particles hold one half
only, the last anchor slot is given to the best-ranked particle of the other
half. When the parent has no particle of the other half, it is skipped
(`cannot_make_parity_complete_anchor`).

The manifest records `anchor_target` and the realized anchor size `n_anchor` of
every entry.

## 9. Child Count Feasibility and Skipped Parents

After anchor selection:

```text
R_i     = N_i - |A_i|
Y_cap_i = R_i if A_i is parity-complete, otherwise 0
```

This means each allocated child receives at least one non-anchor source
particle, and every child has even/odd support through the shared anchor.

Every usable parent supports the same child count:

```text
Y_i = Y_child
```

A parent with `Y_cap_i < Y_child` is skipped (`too_few_non_anchor_particles`).
Exact-copy fill-in does not exist. Expanding every usable parent equally keeps
the output distribution proportional to the original usable class-average
distribution instead of preferentially expanding rare or high-population
parents.

A skipped parent is absent from the whole output, the original block included.
The reasons are:

```text
inactive_cls2D_state                 cls2D state is 0
zero_population                      no particle has the class
missing_particle_membership          the class has no member list
invalid_eo_assignment                a member has an even/odd label other than 0 or 1
too_few_particles_for_osmpl_fac      realized anchor smaller than 2
anchor_leaves_no_child_particles     the anchor takes every member
cannot_make_parity_complete_anchor   no member of the missing half
too_few_non_anchor_particles         fewer non-anchor members than children
```

The manifest lists every skipped parent with its class, population, `cls2D`
state and reason, and the log reports
`>>> BOOTSTRAP_CAVGS CLASSES INPUT/USED/SKIPPED/OUTPUT`. The output has
`K * osmpl_fac` classes for the `K` usable parents.

## 10. Child Assembly

For each expanded parent:

1. copy `A_i` into every child;
2. randomly partition the non-anchor particles into `Y_i` balanced variable
   subsets without replacement;
3. write the bootstrap project: one row per member of each original class
   (role `original`) and one row per member of each child (roles `anchor` and
   `member`). Each row is a copy of the source particle row with `class` set to
   the output class, `stkind` renumbered to the stacks used, `indstk` kept (the
   physical image does not change), `state` set to 1 and `pind_prev` set to the
   source row. The rows are ordered by source stack, so the stack ranges stay
   contiguous. `ptcl3D` is a copy of `ptcl2D`; `cls2D` has one row per output
   class with `bootstrap_parent`, `bootstrap_child` (0 for an original class)
   and `bootstrap_anchor`; `cls3D` is a copy of `cls2D`;
4. run `make_cavgs` on the bootstrap project with `ncls = K * osmpl_fac`, which
   assembles merged, even, and odd averages of every output class
   independently but in identical stack order, from the stored class, in-plane
   and shift parameters.

Under the default partition rule, every non-anchor source particle appears in
exactly one child for that parent, and every anchor source particle appears in
every child.

The command does not perform new particle alignment, class reassignment,
child-to-parent orientation consensus, or mapping of synthetic child
orientations back to original particles. The return path is
`unbootstrap_cavgs` (Section 12).

## 11. Metadata

The manifest (`<refs body>_manifest.txt`) opens with comment lines and has one
row per output class-average entry, copied parents and stochastic children:

```text
# SIMPLE bootstrap_cavgs manifest
# refs: <refs>
# ncls_input: <classes in cls2D>
# ncls_used: <usable parents>
# ncls_skipped: <skipped parents>
# ncls_expected_if_all_input_valid: <ncls_input * osmpl_fac>
# ncls_out: <ncls_used * osmpl_fac>
# osmpl_fac: <osmpl_fac>
# frac_best: <frac_best>
# anchor_target: <N_anchor_target>
# skipped_classes: parent_cls parent_pop cls_state reason      (when parents are skipped)
# skipped <parent_cls> <parent_pop> <cls_state> <reason>       (one line per skipped parent)
out_cls role parent_cls child_id parent_pop n_anchor n_members
```

`role` is `original` (`child_id` 0, `n_members` the parent population) or
`child` (`child_id` 1 to `osmpl_fac - 1`, `n_members` the anchor plus the
variable particles). `n_anchor` is the parent's realized anchor size.

The membership file (`<refs body>_membership.txt`) has one row per row of the
bootstrap project's particle field:

```text
synthetic_pind out_cls role parent_cls child_id source_pind source_rank source_corr eo stkind indstk
```

`synthetic_pind` is the row in the bootstrap project and `source_pind` the row
in the input project. `source_rank` is the rank within the parent in the
objective-rank ordering, `source_corr` the particle's `corr`, `eo` its half,
and `stkind` and `indstk` the stack row of the input project and the physical
image index. `role` is `original`, `anchor` or `member`.

The command takes no seed and records none; the random partition uses the
process random generator. All stochastic choices are traceable through the
membership file, and an exact reproduction uses the manifest and membership
file.

## 12. Handoff and Return Path

The bootstrap package is intended for `solve3D_cavgs` seed generation.
Synthetic child orientations remain attached to the bootstrap project and are
never written into the original project. Downstream continuation uses the
original data, not synthetic child assignments.

`unbootstrap_cavgs` maps a refined bootstrap project back. `projfile` is the
bootstrap project and `projfile_orig` the original project it was made from.
The command clears the 3D alignment of the original `cls3D`, sets all its
states to 0, and copies `corr`, `proj`, the Euler angles, the shift and the
`state` of the `cls3D` entries of the original-block classes
(`bootstrap_child=0`) to their `bootstrap_parent` classes, then maps the
classes to the particles (`map2ptcls`). Only the original block returns;
children are never mapped back. When the bootstrap `cls2D` lacks
`bootstrap_parent` and `bootstrap_child`, the mapping is recovered by majority
vote of `ptcl2D:pind_prev`, with a warning. The result is written to the copy
of `projfile_orig` in the execution directory (`mkdir=yes`, the default) or to
`projfile_orig` itself (`mkdir=no`).

## 13. Failure Conditions

The command stops when:

- the `ptcl2D`, `cls2D` or `stk` segment is empty;
- `osmpl_fac` is less than 1, or `frac_best` is not in (0,1];
- no `cls2D` class is both active and populated;
- no parent is usable after the checks of Sections 8 and 9;
- `make_cavgs` stops on the bootstrap project (for example a particle stack
  that cannot be read).

## 14. Review Checks

A change to the command keeps:

- the controls `osmpl_fac`, `frac_best` and `refs`;
- one command execution producing exactly one ensemble;
- the usable parents' copies in the first `K` stack entries in original index
  order;
- exactly `osmpl_fac - 1` stochastic children for every usable parent;
- exactly one parent per child and no cross-parent combination;
- `N_anchor_target` derived from the median population and `frac_best`;
- anchors selected from the objective-rank best particles;
- parity-complete anchors on expanded parents;
- a reason in the manifest and a count in the log for every skipped parent;
- no new particle alignment, reassignment, orientation consensus, or mapping
  back to the original project inside `bootstrap_cavgs`;
- merged, even, and odd outputs that are complete and identically ordered;
- lineage sufficient to audit every stochastic choice.

Useful validation diagnostics are child-count distributions, even/odd support,
child-vs-parent cavg correlation, cavg-vs-reprojection correlation after
`solve3D_cavgs`, and sensitivity to `osmpl_fac` and `frac_best`.
