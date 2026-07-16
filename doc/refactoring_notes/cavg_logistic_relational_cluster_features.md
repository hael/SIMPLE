# Pairwise-Relationship Neighbourhood Features for Class-Average Logistic Rejection

Date: 2026-07-16

Status: analyze-only CC/FRC and signal-channel foundation implemented;
train/validate/apply model path pending

## Implemented First Slice

The initial analyze-only foundation is now implemented behind
`relational_features=corr_knn_v1`. It does not alter class states or the
existing 14-feature model path. The implementation currently provides:

- named parameters for the correlation policy and neighbourhood definition;
- a shared `calc_cavg_pairwise_algninfo` provider reused by the existing
  CC/FRC distance calculation, so the full-search transform is not discarded;
- upper-triangular CC, FRC, rotation, shift, and mirror primitives;
- explicit neighbour identities, correlations, FRC distances, and normalized
  CC weights;
- the eight recommended wide `corr_knn_v1` summaries;
- a long-form candidate bank over several `k` values for CC and CC-anchored
  FRC reducers;
- a self-contained reducer regression test.

The opt-in command is:

```text
simple_exec prg=model_cavgs_rejection quality_mode=analyze \
  relational_features=corr_knn_v1 projfile=<project> mskdiam=<diameter>
```

The default search policy is `hp=100 Å`, `lp=15 Å`, `trs=10 px`, matching
Joe's coarse-sieve defaults. `relational_knn=5` and
`relational_weight_tau=0.05` remain explicit analysis parameters.

## Implemented Second Slice: Signal-Statistics Channels

The pairwise layer now also exposes every component produced by the existing
`calc_sigstats_dmats` pass:

- raw radial power-spectrum Euclidean distance;
- native histogram TVD, JSD, and Hellinger distance;
- the equal-weight, min-max-normalized histogram composite;
- the existing equal-weight power/histogram signal-statistics composite.

The legacy clustering calculation now consumes the new component result and
moves out the same final composite, preserving its numerical path. Analyze
schema version 2 persists all component values in the upper-triangular pair
table and on each CC-selected neighbour row. The long-form candidate bank
applies `anchor_nearest`, mean, median, quartiles, and spread reducers to FRC
and every signal distance channel while keeping rotational-max CC as the sole
neighbourhood anchor.

Metadata records the channel semantics, provider version, raw versus composite
scales, histogram bin count, power band, mask radius, and histogram intensity
range. This is sufficient to interpret or recompute reducers without reading
the images again.

The implementation boundary now stops before relational artifact I/O,
learning, grouped validation, or apply. It also stops before weighted
aggregation of every base quality feature and cross-channel rank-agreement
statistics; those remain cheap additions if the first empirical runs justify
them.

## Goal

Extend the existing class-average logistic model with continuous per-class
neighbourhood evidence derived from pairwise relationship matrices. Maximum
in-plane rotational correlation is the first and primary relationship, but the
architecture must accept additional similarity or distance matrices without
changing the logistic lifecycle or duplicating neighbourhood code.

There is no clustering step in this design:

- no affinity propagation or k-medoids;
- no cluster count or cluster identifier;
- no assignment of classes to groups;
- no cluster-wide accept, reject, or rescue decision;
- no recursive invocation of `cluster_cavgs` and no intermediate project.

Instead, each class average is a node in a similarity graph. Its relational
features summarize the nearby nodes in that graph and are ordinary inputs to a
regularized logistic model:

```text
logit P(accept_i) = beta_0 + beta_x . x_i + gamma . g_i
```

Here `x_i` is the current 14-feature vector and `g_i` is a small label-free
neighbourhood vector constructed from one or more named pairwise channels. The
model remains a per-class classifier: related images contribute evidence, but
never impose a collective decision.

The relational feature schema is opt-in. Existing commands, analysis files,
model artifacts, hard-gate behaviour, and 14-feature predictions remain on
their current code paths unless the new schema is explicitly selected.

## What to Retain from Joe's POC

Joe's coarse-sieve proof of concept provides two useful observations:

1. class averages should be compared with the existing full in-plane search,
   including rotation, bounded shifts, mirror choice, masking, and the coarse
   low-pass limit;
2. a low `bp_center_edge_var` value is more informative when similar class
   averages also have low values.

The POC turns the second observation into a hard cluster vote: reject all
members when at least 70% of a cluster have `bp_center_edge_var < 4.0`. The
relational model retains the signal but removes both the cluster and the vote.
It supplies the logistic learner with a similarity-weighted neighbourhood
failure fraction and lets training determine its importance.

The threshold of `4.0` remains only as the definition of a diagnostic feature.
It is equivalent to comparing the stored logarithmic feature with `log(4.0)`.
It is not an acceptance threshold, and no neighbour can directly change
another class's state.

## Rotational-Max Correlation Graph

For class averages `i` and `j`, define:

```text
s_ij = max CC(i, transform(j; theta, dx, dy, mirror))
```

over the same full in-plane search used by Joe's `cluster_cavgs` call. The
matrix is symmetric, has an excluded diagonal, and is calculated only among
rows eligible for model scoring after the context's validity gates.

For each class `i`, let `N_k(i)` be the `k` highest-correlation neighbours,
excluding `i`. When fewer than `k` eligible classes exist, use all available
non-self neighbours and record the actual neighbour count.

For weighted neighbourhood quantities use:

```text
w_ij = exp((s_ij - max_l s_il) / tau),  j in N_k(i)
```

where `tau` is a feature-schema parameter recorded in every analysis and model
artifact. Subtracting the row maximum stabilizes the exponential without
changing relative weights.

This is a local graph construction, not a clustering algorithm. Neighbourhoods
may overlap, need not be reciprocal, and do not partition the dataset.

## Extensible Pairwise-Relationship Layer

The rotational-max correlation matrix must enter the quality model through a
general pairwise-relation contract rather than a CC-specific array passed
directly into feature code. A relationship channel consists of:

```text
name                 stable channel name, for example rotmax_cc
provider             implementation and provider version
semantics            similarity | distance
direction            higher_is_closer | lower_is_closer
values               symmetric nclass x nclass matrix
valid_pairs          symmetric validity mask; diagonal is always invalid
parameters           canonical calculation settings
transform            optional distance-to-affinity transform
```

The channel validator requires the same class ordering as the base feature
rows, symmetric dimensions and validity, finite valid entries, and an excluded
diagonal. Missing pair values are represented only by `valid_pairs`; magic
numeric sentinels are not allowed.

Initial and plausible future channels include:

| Channel | Semantics | Existing or likely source |
| --- | --- | --- |
| `rotmax_cc` | similarity | raw aligned CC matrix from `calc_cc_and_res_dmats` |
| `rotmax_frc_distance` | distance | pairwise FRC distance already calculated beside CC |
| `signal_stats_distance` | distance | signal/power/histogram distances used by `calc_sigstats_dmats` |
| `radial_power_distance` | distance | current `pspecs%calc_distmat` output used by `calc_sigstats_dmats` |
| `ctf_power_defocus_distance` | distance | future integrated CTF-corrected power/defocus comparison |
| `feature_space_distance` | distance | future distance over selected per-class measurements or embeddings |

The feature schema, not the model scorer, declares which channels it consumes
and how. Adding a channel therefore requires a provider plus schema entries;
it does not require another classifier implementation.

### Anchor and Secondary Channels

Every schema names one **anchor channel** that defines `N_k(i)` and the graph
weights. In `corr_knn_v1`, the anchor is `rotmax_cc`, preserving the intended
rotational-max CC neighbourhood.

Secondary channels are first summarized over the anchor neighbourhood. For
example, an FRC distance channel can contribute the median FRC distance to the
same CC-selected neighbours, and a signal-statistics channel can report whether
those neighbours also agree in signal character. This has three advantages:

- all channel features refer to the same interpretable local family;
- the model can learn agreement or disagreement between relationship types;
- adding a matrix does not silently replace the graph used by Joe-derived
  features or neighbour-quality consensus.

A future schema may explicitly choose another anchor or define a second
channel-specific neighbourhood, but that is a new versioned schema. It must
not change the meaning of `corr_knn_v1`.

### Similarity and Distance Normalization

Reducers operate on a channel in its native semantics. A distance channel can
therefore expose local mean, median, quantile, margin, and cohesion distances
without pretending that it is a correlation.

If a distance channel is used for neighbour selection or weighting, its schema
must declare an affinity transform, for example:

```text
a_ij = exp(-d_ij / tau_channel)
```

The transform and scale are part of the channel parameters recorded in the
analysis and model artifact. There is no implicit `1 - distance` conversion.
After the raw per-class reducers are computed, the resulting scalar features
are robustly normalized per project in the same way as other model features.

### Reusable Reducers

The relational feature layer should provide channel-independent reducers for:

- top-`k` mean, median, quantiles, and nearest-versus-next-shell margin;
- pairwise cohesion within an anchor neighbourhood;
- value to a local neighbourhood representative;
- effective and population-weighted support;
- isolation or nearest valid relationship;
- similarity-weighted aggregation of another class-level value;
- cross-channel agreement on the anchor neighbourhood.

`corr_knn_v1` is then a declaration using these reducers, with the
Joe-derived centre/edge fraction implemented as a weighted aggregation of a
class-level indicator. Future schemas compose the same reducers with other
channels rather than adding one-off feature-extraction loops.

## Analyze-First Implementation Boundary

The first implementation should be deliberately broad in `analyze` and narrow
in `apply`. The boundary is between **observations that are expensive or
impossible to recover later** and **model choices that are cheap to revisit
from saved observations**.

The implementation has three layers:

```text
pairwise primitives       calculate once and persist
        |
candidate statistic bank derive broadly in analyze
        |
promoted model schema     select empirically; required by learn/apply
```

### 1. Pairwise Primitives: Implement and Persist Now

An analyze run should retain every useful quantity already produced by the
same class-average comparison passes. This avoids repeating the expensive
full-search alignment when a new reducer or feature combination is proposed.

The initial upper-triangular pair table should contain:

```text
class_i,class_j,valid,
rotmax_cc,frc05_distance,
best_rotation,best_shift_x,best_shift_y,best_mirror,
radial_power_distance,histogram_tvd,histogram_jsd,histogram_hellinger,
histogram_distance,signal_stats_distance
```

The first line comes from the same `match_imgs` result that currently feeds
`calc_cc_and_res_dmats`: CC, FRC, and the maximizing transform are already
available together. The power and histogram distances are already calculated
inside `calc_sigstats_dmats` and then collapsed into `dmat_sig`; expose the
components instead of discarding them. This remains close to the total work
performed by Joe's `cluster_cavgs` route and does not introduce another image
alignment objective.

Write channel metadata and parameter fingerprints beside the pair table. For
larger class sets, allow an equivalent compact binary matrix bundle, but keep
the text metadata and class-index map human-readable. The persisted data must
be sufficient to recompute neighbourhood membership and all reducers without
reading images again.

### 2. Candidate Statistic Bank: Compute Broadly in Analyze

The channel-independent reducers are cheap once the matrices exist. Analyze
mode should therefore emit a long-form candidate table rather than committing
all experimental statistics to the permanent model feature vector:

```text
class,channel,reducer,k,tau,raw_value,normalized_value,valid_count
```

For each available channel, calculate where meaningful:

- nearest value, top-`k` mean, median, lower/upper quartiles, and spread;
- nearest-versus-next-shell margin;
- local cohesion and value to the local representative;
- effective support, population-weighted support, and isolation;
- weighted mean, median, and dispersion of each existing per-class quality
  feature over the CC anchor neighbourhood;
- weighted fractions for named diagnostic predicates, initially including
  Joe's `bp_center_edge_var < 4.0` predicate;
- agreement between channel rankings and secondary-channel summaries over the
  CC anchor neighbourhood.

Use a small declared analysis grid for `k` and `tau`, recorded in metadata.
The default grid should cover small and moderate neighbourhoods without
creating hundreds of nearly duplicate columns; the persisted pair table is the
escape hatch for broader offline sweeps. Candidate statistics belong in the
long-form file and are not automatically accepted by the model learner.

Also emit one convenient wide per-class table for the current recommended
`corr_knn_v1` settings. This supports immediate plotting and comparison with
manual states while the long-form bank preserves alternatives.

### 3. Promoted Feature Schema: Keep Small and Versioned

`learn` consumes only a named, ordered feature schema promoted from the
candidate bank. It does not ingest every candidate column and allow an
uncontrolled feature search. Promotion records the selected channels,
reducers, `k`, `tau`, transforms, and predicates in a new schema version.

`apply` computes only the features named by that artifact. It does not write
the pair table or candidate bank unless explicitly requested through an
analysis/debug option. This keeps the production path predictable even though
analyze mode is deliberately rich.

### Deferred Beyond the First Analyze Implementation

Stop the initial implementation at statistics derivable from class-average
images, current class metadata, and the pairwise passes above. Defer:

- model-derived `neighbour_quality`, because it requires leakage-safe
  out-of-fold predictions;
- arbitrary pairwise feature interactions and automated feature-subset search;
- particle-level defocus/power trend fitting, which requires a new data pass;
- new alignment objectives or learned pairwise embeddings;
- any statistic that changes class state during analysis;
- all clustering, voting, and group-level decisions.

This boundary maximizes empirical information from the initial implementation
without coupling exploratory statistics to the trained artifact or production
apply path.

## Relational Feature Schema `corr_knn_v1`

The original cluster-derived ideas map naturally onto local, label-free graph
statistics:

| Original concept | Label-free replacement |
| --- | --- |
| cluster cohesion | cohesion within the local `k`-neighbourhood |
| medoid similarity | similarity to a local neighbourhood representative |
| assignment margin | density contrast between the nearest and next shells |
| cluster support | effective or particle-weighted local support |
| singleton evidence | low best-neighbour correlation and low effective support |

The initial schema contains the following raw features.

1. `corr_local_density`

   Mean of the top `k` rotational-max correlations:

   ```text
   corr_local_density_i = mean_{j in N_k(i)} s_ij
   ```

   This is the primary relational feature and the first feature to test on its
   own.

2. `corr_local_cohesion`

   Median pairwise correlation among the members of
   `{i} union N_k(i)`, excluding diagonal terms. Density asks whether `i` has
   close neighbours; cohesion asks whether those neighbours also resemble one
   another.

3. `corr_to_local_medoid`

   Within `{i} union N_k(i)`, choose the row with the highest mean correlation
   to the other rows and use `i`'s correlation to that representative. This is
   a local summary calculation only. It does not create a global medoid set or
   assign any other class to the representative. If `i` itself is the local
   representative, use its mean non-self correlation rather than a diagonal
   value of one.

4. `corr_neighbourhood_margin`

   Mean correlation to the nearest `k` neighbours minus mean correlation to
   the next `k` classes in rank order:

   ```text
   margin_i = mean(top-k s_ij) - mean(next-k s_ij)
   ```

   This is the label-free analogue of assignment confidence. It measures
   whether the local neighbourhood is distinct from the surrounding graph.
   Use all available rows in a short second shell and record when no second
   shell exists.

5. `corr_effective_support`

   Effective number of contributing neighbours:

   ```text
   n_eff_i = (sum_j w_ij)^2 / sum_j w_ij^2
   ```

   Use `log(1 + n_eff_i)` in the logistic feature vector. This distinguishes a
   broad family from a class supported by only one dominant match.

6. `corr_population_support`

   Similarity-weighted particle support:

   ```text
   log(1 + sum_{j in N_k(i)} w_ij * population_j)
   ```

   This is preferred to class count when the class populations differ greatly.
   Both support features are retained initially so held-out training can show
   which representation transfers better.

7. `corr_isolation`

   `1 - max_{j != i}(s_ij)`. This explicitly represents singleton-like or
   genuinely unusual classes without inventing a singleton cluster. It remains
   soft evidence: a strong individual feature vector may still support an
   isolated but valid class.

8. `corr_neighbour_bp_fail_fraction`

   Similarity-weighted fraction of neighbours for which
   `bp_center_edge_var < 4.0`:

   ```text
   sum_j w_ij * I(bp_center_edge_var_j < 4.0) / sum_j w_ij
   ```

   This is the continuous, label-free counterpart of Joe's 70% group test.
   The class itself is excluded. No 70% decision is made.

All relational features are robustly normalized per project using the same
median/MAD convention as the base features. Missing-neighbour cases receive
explicit deterministic handling and a recorded neighbour count; they never
inherit self-correlation.

`corr_knn_v1` consumes only the `rotmax_cc` channel. This keeps the initial
experiment close to the original proposal and provides the compatibility
baseline for future multi-channel schemas. Candidate additional matrices are
written as diagnostics before they are admitted into a trained schema.

## Stage 1 Logistic Model

The first model uses the existing base features plus `corr_knn_v1`:

```text
logit P_i = beta_0 + beta_x . x_i
          + gamma_1 corr_local_density_i
          + gamma_2 corr_local_cohesion_i
          + gamma_3 corr_to_local_medoid_i
          + gamma_4 corr_neighbourhood_margin_i
          + gamma_5 corr_effective_support_i
          + gamma_6 corr_population_support_i
          + gamma_7 corr_isolation_i
          + gamma_8 corr_neighbour_bp_fail_fraction_i
```

Fit a regularized logistic regression and learn its probability threshold and
calibration exactly as model evaluation requires. Begin feature ablations with:

1. base model only;
2. base plus `corr_local_density`;
3. base plus all image-only neighbourhood summaries;
4. base plus the Joe-derived neighbour failure fraction;
5. the complete `corr_knn_v1` schema.

An interaction between the class's own normalized centre/edge feature and
`corr_neighbour_bp_fail_fraction` is a candidate, because it is the smoothest
analogue of Joe's POC. It should enter the schema only if grouped validation
shows value beyond the two main effects. There is no automatic expansion to
all pairwise interactions.

## Stage 2 Neighbour-Quality Consensus

If label-free image-neighbourhood features improve held-out performance, add
one model-derived graph feature:

```text
neighbour_quality_i =
    sum_{j in N_k(i)} w_ij * logit(P_j^base)
    / sum_{j in N_k(i)} w_ij
```

The final model becomes:

```text
logit P_i^final = beta_0 + beta_x . x_i + gamma . g_i
                + delta * neighbour_quality_i
```

This remains soft evidence. It is neither majority voting nor a hard rescue.
Coherent fuzzy balls, ice, or other junk may support one another, so the final
coefficient must be learned and validated rather than assumed positive.

Training must generate `P_j^base` out of fold:

1. split by complete project or specimen;
2. fit the baseline model without the held-out group;
3. predict baseline logits for every row in the held-out group;
4. calculate neighbourhood consensus within that held-out project's graph;
5. fit and assess the final model using grouped or nested cross-validation.

Using in-sample fitted baseline probabilities would leak training fit into the
consensus feature and is prohibited. At apply time, the artifact's frozen base
model computes all baseline logits first; the graph feature and final logistic
probability are then computed in a second pass.

## Correlation Policy: Match Joe's POC

SIMPLE already contains the expensive part. `calc_cluster_cavgs_dmat` calls
`calc_cc_and_res_dmats`, which performs full in-plane matching and constructs a
symmetric correlation matrix before converting it to a distance matrix. The
refactor should expose that raw matrix through a reusable helper such as
`calc_cavg_rotmax_ccmat`.

The helper must reproduce the search policy used by Joe:

| Setting | `corr_knn_v1` policy | POC reference |
| --- | --- | --- |
| Objective | correlation (`cc`), `ctf=no` | `cluster_cavgs` commander defaults |
| Rotation | full in-plane rotational maximization | existing `match_imgs` path |
| Mirror | compare mirrored and unmirrored alternatives | existing matcher |
| High-pass | `100 Å` default | inherited `hp` default |
| Low-pass | explicit; coarse analysis uses `coarse_defaults%lpstop` | explicitly passed by Joe |
| Shift half-range | `10 px` default | `cluster_cavgs` commander default |
| Mask | caller's `mskdiam` | passed by Joe |
| Threads | runtime setting only | Joe currently uses `nthr=1` |

No clustering criterion or clustering method belongs in the relational schema.
The graph depends only on the rotational-max correlation matrix and the
neighbourhood parameters.

Expose the policy through named parameters:

```text
relational_features        = none | corr_knn_v1       (default: none)
relational_channels        = rotmax_cc                 (fixed by corr_knn_v1)
relational_anchor          = rotmax_cc                 (fixed by corr_knn_v1)
relational_corr_hp         = <Angstrom>               (default: 100.0)
relational_corr_lp         = <Angstrom>               (default: 15.0)
relational_corr_trs        = <pixels>                 (default: 10.0)
relational_knn             = <integer>                (default: 5)
relational_weight_tau      = <correlation units>      (schema default)
```

For `apply`, these values come from the learned model artifact. Command-line
overrides must match the artifact exactly or fail clearly. `analyze` is the
mode in which alternative low-pass limits, shift ranges, `k`, and `tau` are
deliberately explored and recorded.

Future schemas may expose multiple values in `relational_channels`. Channel
selection belongs to analysis/schema creation; `apply` always uses the exact
ordered channel list and provider parameters embedded in the trained artifact.

## Analyze, Learn, Evaluate, and Apply

The existing `model_cavgs_rejection` program remains the entry point. The
current `quality_context=sieve` hard-gate route remains unchanged when
`relational_features=none`. Selecting `corr_knn_v1` opts that context into
relational feature extraction and logistic scoring.

| Mode | Relational behaviour |
| --- | --- |
| `quality_mode=analyze` | Build and persist pairwise primitives, write the broad candidate statistic bank plus the recommended wide feature table, optionally score a supplied relational model, compare with manual states, and leave the project unchanged. |
| `quality_mode=learn` | Read a file table of schema-compatible analysis files, fit regularized logistic models, perform grouped cross-validation and feature ablations, and write a versioned artifact. |
| `quality_mode=evaluate` | Apply a frozen relational artifact to independent analysis files and report project-level and macro metrics without modifying projects. |
| `quality_mode=apply` | Recompute the artifact-defined graph, calculate probabilities, annotate class averages, map accepted classes to particles, and write the project. A trained relational artifact is required. |

This lifecycle also applies later to `chunk` and `pool`; the initial validation
target is the coarse sieve data for which Joe has the twelve-dataset test set.

## Artifact and Backward Compatibility

Do not increase `CAVG_QUALITY_NFEATS` in place or reinterpret current model
files. Existing models have fixed 14-element arrays and a fixed interaction
layout.

Introduce a versioned relational model artifact and dispatch by its header:

```text
model_version=10
model_family=relational_logistic
context=sieve
base_feature_schema=cavg_quality_base_v1
relational_feature_schema=corr_knn_v1
relation_channels=rotmax_cc
relation_anchor=rotmax_cc
relation_channel.rotmax_cc.semantics=similarity
relation_channel.rotmax_cc.provider=rotmax_cc_v1
relation_channel.rotmax_cc.parameters=hp,lp,trs,mskdiam,mirror,objfun
base_feature_names=...
relational_feature_names=...
neighbourhood_policy=knn,tau
base_coefficients=...
relational_coefficients=...
interactions=...
intercept=...
prob_threshold=...
regularization_lambda=...
calibration_temperature=...
```

The current model reader, writer, learner, and classifier continue to own
legacy artifacts through version 9. A relational header selects a separate
schema-aware model path. The trainer rejects file tables that mix feature
schemas, channel order, provider versions, channel parameters, low-pass limits,
shift ranges, masks, `k`, `tau`, or feature order.

The artifact for Stage 2 additionally embeds the frozen baseline model needed
to reproduce `neighbour_quality`. This keeps apply deterministic and prevents
the runtime from substituting a different baseline.

## Implementation Ownership

```text
simple_ui_params_common / simple_ui_cavgproc
    relational feature and correlation-policy parameters
                 |
simple_commanders_cavgs: model_cavgs_rejection
    mode validation and legacy/relational dispatch
                 |
simple_strategy2D_utils
    rotational-max CC/FRC relationship provider
                 |
simple_cavg_quality_relations
    current analyze provider result, reducers, and diagnostic output;
    evolve into channel contract plus schema declarations
                 |
relational model and learning code
    artifact I/O, training, scoring, grouped validation
```

Refactor `calc_cc_and_res_dmats` to consume the shared raw-correlation helper,
then verify that its resulting `dmat_cc` is unchanged. The relational extractor
uses the same helper directly. It must not invoke clustering code or a
commander, and must not alter class states while extracting features.

The CC helper is the first relationship provider. Keep channel construction
separate from reducers so the existing FRC matrix and future signal-statistics
or power-spectrum providers can populate the same validated contract.

## Analysis and Validation Outputs

Relational analysis files extend the current row-oriented output with schema
and policy metadata:

```text
# relational_feature_schema=corr_knn_v1
# relation_channels=rotmax_cc
# relation_anchor=rotmax_cc
# relation_channel.rotmax_cc.provider=rotmax_cc_v1
# relational_corr_hp=...
# relational_corr_lp=...
# relational_corr_trs=...
# relational_knn=...
# relational_weight_tau=...
class,manual_state,hard_reject,base_probability,relational_probability,
corr_local_density,corr_local_cohesion,corr_to_local_medoid,
corr_neighbourhood_margin,corr_effective_support,corr_population_support,
corr_isolation,corr_neighbour_bp_fail_fraction
```

Also write the actual neighbour indices, correlations, and normalized weights
for every class in a separate diagnostic table. This is essential for visual
inspection and for verifying that a false decision is supported by sensible
neighbours rather than a numerical summary alone.

When additional diagnostic channels are present, write their value and
validity for the same neighbour pairs in that table. This makes cross-channel
agreement inspectable before a new matrix becomes part of a trainable feature
schema.

The complete analyze artifact set is therefore:

```text
cavgs_quality_analysis.txt
cavgs_pairwise_relations.txt          upper-triangular pairwise primitives
cavgs_relational_candidates.txt       long-form reducer/channel/parameter bank
cavgs_relational_features.txt         recommended wide per-class schema
cavgs_relational_neighbours.txt       neighbour identities, values, weights
selected/rejected/ranked MRC stacks   existing visual diagnostics
```

Validation holds out complete projects or specimens and compares models on
identical hard-gate survivors. Report macro AUC, balanced accuracy, F1,
manual-good recall, manual-bad specificity, calibration, particle-weighted
effects, and worst-project behaviour. Inspect false-positive and false-negative
stacks, with particular attention to:

- the known dataset where Joe's POC rejects acceptable classes;
- overfit classes that survive Joe's hard group rule;
- coherent fuzzy-ball or ice neighbourhoods;
- genuinely good but isolated classes;
- sensitivity to low-pass limit, shift range, `k`, and `tau`;
- stability when classes are added to or removed from a project.

## Defocus–Power Diagnostic

The integrated CTF-corrected particle-power versus defocus idea can also be
made label-free: compare each class's fitted trend or residual summary with
those of its rotational-correlation neighbours. Keep this outside
`corr_knn_v1` initially. Add it only if it contributes held-out information
beyond the image-quality and graph features.

## Implementation Sequence

Completed in the analyze-only foundation:

1. Extract and test the shared symmetric rotational-max correlation matrix,
   reproducing Joe's POC search settings without running `cluster_cavgs`.
2. Introduce the validated pairwise-relation channel contract and wrap the CC
   and FRC outputs from the same match as the first provider.
3. Expose the power-spectrum, TVD, JSD, Hellinger, and merged signal-statistics
   distances already calculated by `calc_sigstats_dmats`.
4. Write the persisted pair table and neighbour diagnostic table in
   `quality_mode=analyze`.
5. Implement the reusable reducer bank, long-form candidate output, and the
   recommended wide `corr_knn_v1` table.

Next:

6. Analyze Joe's twelve datasets and promote a compact relational feature
   schema based on grouped empirical evidence.
7. Implement the versioned relational logistic learner, grouped evaluation,
   and feature-ablation report.
8. Implement `apply`, enforcing the artifact's channel and neighbourhood
   policies.
9. Validate the trained Stage 1 model on held-out projects. Add cross-fitted
   `neighbour_quality` only if the label-free model demonstrates held-out
   value.
