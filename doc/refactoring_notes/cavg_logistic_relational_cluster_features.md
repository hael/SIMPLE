# Clustering Evidence in Class-Average Logistic Rejection

Date: 2026-07-14  
Status: proposal

## Opportunity

The current class-average logistic model treats each class as one feature row.
It does not directly use relationships between class averages. A pairwise
similarity matrix, for example the maximum correlation over in-plane rotation,
contains additional evidence about whether a class belongs to a coherent family
or is isolated from the rest of the dataset.

The cluster identifier itself should not be used as a logistic feature. Cluster
numbers are arbitrary, change between datasets, and have no transferable
meaning. The useful information is instead contained in continuous quantities
such as local density, cluster cohesion, and assignment confidence.

## Proposed Model

Let `x_i` be the existing normalized feature vector for class average `i`, and
let `g_i` contain relational features derived from the similarity matrix. The
logistic model becomes:

```text
logit P(accept_i) = beta_0 + beta * x_i + gamma * g_i
```

This preserves the current classifier and model artifact: relational evidence
appears as additional scalar features rather than as a separate clustering
decision.

## Candidate Relational Features

For pairwise similarity

```text
s_ij = max correlation between class i and class j over in-plane rotation
```

the most useful first features are:

1. `corr_local_density`

   Mean or median similarity to the nearest `k` other class averages, excluding
   self-correlation. This measures whether the class has nearby support without
   requiring a hard cluster assignment.

2. `corr_to_medoid`

   Similarity to the assigned cluster medoid. For a medoid, use a leave-one-out
   within-cluster statistic rather than its self-correlation.

3. `cluster_assignment_margin`

   Similarity to the assigned medoid minus similarity to the next-best medoid.
   This measures whether the cluster assignment is clear or ambiguous.

4. `cluster_cohesion`

   Median pairwise similarity within the assigned cluster. This is shared by
   cluster members and measures the compactness of the group.

5. `log_cluster_support`

   Log class count or summed particle population in the cluster. Particle
   population is likely more meaningful, but both definitions should be tested
   because cluster splitting can distort either quantity.

All features should be dataset-normalized with the existing robust procedure.
Singleton clusters need explicit, deterministic handling; they must not inherit
a perfect self-similarity value.

## Prefer Soft Neighborhoods First

The lowest-risk first experiment is `corr_local_density`, calculated from a
fixed nearest-neighbor rule. It avoids choosing a cluster count and is less
sensitive to unstable cluster boundaries.

Hard-cluster features can then be added if medoid similarity, assignment margin,
or cluster support provide information beyond local density. The classifier
should remain free to reject a member of an otherwise coherent cluster or accept
an unusual but individually convincing class.

## Possible Consensus Extension

A later model could use the baseline logistic predictions of similar class
averages:

```text
neighbor_quality_i = weighted mean of baseline logits for neighbors j != i
```

The final logistic model would combine the class's own evidence with this
neighbor consensus. This is a compact, interpretable graph model and may help
when several related classes independently show moderate evidence.

Training requires out-of-fold baseline predictions. Using fitted predictions
from the same rows would leak training fit into the consensus feature. At apply
time, baseline predictions are computed first, relational features second, and
the final logistic probability last.

Consensus should remain soft evidence. Coherent fuzzy balls, ice, or other junk
can reinforce one another just as valid molecular classes can.

## Existing SIMPLE Machinery

`calc_cluster_cavgs_dmat` already constructs signal, correlation, resolution,
and hybrid distance matrices, including pairwise correlation through full
in-plane search. `align_and_score_cavg_clusters` already calculates medoid
assignments, population, homogeneity, cluster compactness, and related scores.

The initial implementation should reuse these calculations or extract their
common distance-matrix and summary logic. It should not duplicate the expensive
pairwise image comparisons inside `simple_cavg_quality`.

## Evaluation

Compare the following on identical hard-gate survivors:

1. existing logistic baseline;
2. baseline plus `corr_local_density`;
3. baseline plus the compact hard-cluster feature set;
4. baseline plus cross-fitted neighbor consensus.

Hold out complete projects or specimens rather than random class-average rows.
Evaluate chunk and pool separately, and report manual-good recall, manual-bad
specificity, worst-dataset behavior, and visually inspected error stacks.

Important checks are:

- whether relational features add information beyond population and resolution;
- whether smooth junk forms misleading high-density clusters;
- sensitivity to the number of classes and the chosen neighbor or cluster rule;
- stability when classes are added to or removed from the dataset;
- runtime of the all-pairs in-plane rotational comparison;
- coefficient stability under grouped resampling.

## Recommended First Step

Add only `corr_local_density` and, if a hard clustering is already available,
`corr_to_medoid` and `cluster_assignment_margin` as diagnostic columns with zero
coefficients in current models. If these show held-out conditional value, fit a
small relational feature family before attempting neighbor-consensus scoring.
