# Pairwise-Relationship Neighbourhood Features for Class-Average Logistic Rejection

Date: 2026-07-16

Status: Stage 1 implemented end to end. Analyze writes one canonical training
table, the pairwise learner reads it directly, model artifact version 10 stores
the relational schema and coefficient, and `chunk100mics` uses relational
scoring by default. Independent-dataset validation remains the next step.

## Implemented Stage 1

The implemented feature is deliberately smaller than the exploratory feature
bank:

```text
schema  = corr_knn_signal_v1
anchor  = rotational-max CC
k       = 5
feature = signal_stats_anchor_topk_mean
```

For every class surviving the unchanged hard gates, the code finds the five
other classes with the highest full-search rotational-max CC and averages the
existing composite signal-statistics distance to those neighbours. The result
is robust-normalized within the dataset and enters the current pairwise
logistic model as one additional scalar coefficient. There is no clustering,
cluster identifier, cluster-level decision, or relational interaction bank.

The model owns the search fingerprint: `ctf=no`, `objfun=cc`, `hp=100 Å`,
`lp=15 Å`, `trs=10 px`, and `k=5`. These values mirror Joe's POC defaults and
are recorded in both the canonical table and version-10 model artifact. They
are not user-facing tuning knobs in the production route.

The shared pairwise layer still exposes rotational-max CC/FRC plus the power,
TVD, JSD, Hellinger, histogram-composite, and final signal-statistics component
matrices. That extensible in-memory boundary is retained for future feature
experiments, but analyze no longer emits the exploratory pair, neighbour,
candidate, or wide-feature sidecars.

`quality_mode=analyze` now writes exactly one learner input per dataset:

```text
cavgs_quality_training.txt
```

It contains labels, hard-gate state, all fourteen raw and normalized base
features, the raw and normalized promoted relational feature, and the complete
schema/provider metadata needed by training. Existing visual MRC stacks remain
diagnostics, not learner inputs.

The one-time `scripts/pack_cavgs_quality_training.py` utility converts the old
five-file bundles to this canonical format and can remove the legacy text
artifacts only after the canonical file is written successfully.

### Relational Training Invariant

Relational evidence is now an invariant of newly generated training data, not
a property inherited accidentally from the model used during analyze:

- every `quality_mode=analyze` run computes and writes
  `corr_knn_signal_v1`, even when the scoring route is pool, sieve, linear, or
  a version-9 artifact;
- `quality_mode=learn` accepts only logistic training and requires every file
  to declare the same supported relational schema and policy;
- newly trained chunk and pool artifacts therefore write model version 10 and
  apply through the relational path;
- existing pool, sieve, linear, and version-9 artifacts remain valid for
  application and evaluation, but do not define the contents of new training
  tables.

## First Empirical Result: `chunk_training5`

The first analysis corpus contains twelve datasets and 1,200 class averages.
The unchanged hard gate removes 624 classes, leaving 576 trainable rows: 330
manual accepts and 246 manual rejects. All eligible pair tables are complete,
all eligible classes have the expected five neighbour rows, and the relational
tables join one-to-one by class. The `not` dataset has sixteen eligible accepts
and no eligible reject, so it can contribute training rows but cannot produce a
within-dataset AUC.

A project-held-out regularized logistic comparison was performed on the eleven
measurable datasets. The reference model was retrained from the existing
fourteen normalized base features for every fold; relational statistics were
calculated without labels and normalized within each dataset. Training projects
received equal total weight. The reference L2 penalty was 1.0, with sensitivity
checked from 0.1 through 100. The main result is:

| Feature set | Macro held-out AUC | Change from base |
| --- | ---: | ---: |
| fourteen base features | 0.836 | - |
| base + CC local density | 0.848 | +0.012 |
| base + FRC distance, CC-anchor top-five mean | 0.847 | +0.011 |
| base + signal-statistics distance, CC-anchor top-five mean | **0.860** | **+0.024** |
| base + all eight current CC-wide features | 0.821 | -0.015 |

The signal-statistics feature improved ten of eleven measurable held-out
datasets. A dataset-level bootstrap of the paired AUC change gave a 95% interval
of +0.012 to +0.035 at the reference regularization strength, and the direction
of improvement persisted across the tested regularization range. Adding
density or FRC on top of that feature did not improve macro AUC further.
Because the feature was selected on this corpus, the interval is descriptive;
it is not a substitute for evaluation on a new frozen holdout.

The exposed signal components justify retaining the existing composite rather
than promoting every component:

| Added CC-anchor top-five mean | Macro held-out AUC |
| --- | ---: |
| radial power distance | 0.837 |
| histogram TVD | 0.853 |
| histogram JSD | 0.852 |
| histogram Hellinger | 0.856 |
| all four component features | 0.849 |
| existing signal-statistics composite | **0.860** |

TVD, JSD, and Hellinger are highly redundant in this corpus (pairwise Spearman
correlations above 0.95 for their top-five median summaries), while radial
power is nearly independent but adds little by itself. The individual
components therefore remain analysis diagnostics rather than trainable schema
members.

The current frozen `chunk100mics` score already has macro AUC 0.971 on these
files. Recalibrating that score and adding relational statistics changes AUC by
at most about 0.002. This corpus is suitable for choosing a compact candidate
schema, but it is not sufficient evidence for final coefficients or a claim of
improvement over the frozen production model; the current score may already
encode closely related training information. Final acceptance requires a
frozen artifact and genuinely independent project-held-out evaluation.

### Promoted Schema

Promote one relational statistic for the first learner implementation:

```text
relational_feature_schema=corr_knn_signal_v1
relation_anchor=rotmax_cc
relational_knn=5
feature=signal_stats_distance|anchor_topk_mean
```

For class `i`, this is the unweighted mean signal-statistics distance to the
five other classes having the highest rotational-max CC with `i`. It is then
normalized across hard-gate survivors using the existing robust normalization.
Lower distance is evidence for acceptance; the fitted logistic coefficient,
not a hard rule, determines its contribution. The CC search fingerprint and
the complete signal-statistics provider fingerprint remain part of the schema.

This keeps rotational-max CC as the neighbourhood definition, uses the best
empirical secondary channel, performs no clustering, and adds only one degree
of freedom to the fourteen-feature model. FRC and the individual
power/histogram components remain available through the pairwise provider for
future datasets and ablations; they are not persisted by routine analyze runs.

### Full-Corpus Training and Promotion

After the schema was frozen, all twelve `chunk_training5` datasets were packed
to canonical files and used for the requested full-corpus fit, with no holdout
at this stage. The selected model has:

```text
model_family=pairwise_logistic
feature_policy=microchunk_plus_score_signal
prob_threshold=0.35
regularization_lambda=0.001
relational_coefficient=-0.3941003
macro_learn_score=0.50786
```

The negative coefficient is consistent with the feature semantics: larger
signal-statistics distance from the five CC-nearest neighbours is rejection
evidence. These coefficients are now the built-in `chunk100mics` preset.
Training and built-in evaluation reports agree exactly apart from the model
name. This is a training-set result, not independent validation; the other
assembled datasets should be used for that assessment.

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

The policy is model-owned rather than exposed as ordinary application
parameters. A model either declares `corr_knn_signal_v1` and its fixed
`hp/lp/trs/k` fingerprint or declares `none`. This prevents analyze, training,
and application from silently using different neighbourhoods. Future schemas
can add other named pairwise channels without changing the current artifact.

## Analyze, Learn, Evaluate, and Apply

The existing `model_cavgs_rejection` program remains the entry point. Relational
calculation is automatic when the selected model supports it.

| Mode | Relational behaviour |
| --- | --- |
| `quality_mode=analyze` | Compute the model-required neighbourhood statistic, compare with manual states, write `cavgs_quality_training.txt`, and leave the project unchanged. |
| `quality_mode=learn` | Read a file table of schema-compatible canonical files, fit the regularized pairwise logistic model including the relational scalar, and write a versioned artifact and report. |
| `quality_mode=evaluate` | Apply a frozen relational artifact to independent analysis files and report project-level and macro metrics without modifying projects. |
| normal application | Recompute the artifact-defined neighbourhood, calculate probabilities, and use the existing class/particle update path. |

The built-in `chunk100mics` preset declares the relational schema and therefore
uses this route by default. Existing pool, sieve, linear, and version-9
artifacts declare no relational schema and retain their existing base-only
application behaviour. A newly trained pool-context artifact carries the
schema and uses the relational route.

## Artifact and Backward Compatibility

Do not increase `CAVG_QUALITY_NFEATS` in place or reinterpret current model
files. Existing models have fixed 14-element arrays and a fixed interaction
layout.

Artifact version 10 extends the existing pairwise-logistic format without
changing the fourteen-feature arrays or interaction layout:

```text
model_version=10
model_family=pairwise_logistic
context=chunk
relational_feature_schema=corr_knn_signal_v1
relational_knn=5
relational_corr_hp=100.0
relational_corr_lp=15.0
relational_corr_trs=10.0
relational_coefficient=...
```

The current model reader, writer, learner, and classifier continue to own
legacy artifacts through version 9. A version-9 file has no relational schema
and is deliberately evaluated through the unchanged base-only path. The
trainer rejects file tables that mix relational schemas, neighbour counts, or
correlation-policy values.

Provider metadata distinguishes schema-fixed policy from dataset-resolved
geometry. For example, the twelve initial files share the signal mask policy
but have resolved mask radii from 56 to 136 pixels because particle diameter,
sampling, and box geometry differ. The learner must validate each resolved
value against the declared policy, not require the same pixel radius across
datasets. The artifact stores the policy; each canonical training file records
both its inputs and resolved value.

If Stage 2 neighbour quality is ever added, its artifact must additionally
embed the frozen baseline model needed to reproduce that feature.

## Implementation Ownership

```text
simple_commanders_cavgs: model_cavgs_rejection
    mode validation, canonical output, and model-driven dispatch
                 |
simple_strategy2D_utils
    rotational-max CC/FRC relationship provider
                 |
simple_cavg_quality_relations
    signal-component matrices, CC-neighbour selection, and promoted reducer
                 |
simple_cavg_quality_model / simple_cavg_quality_learn
    artifact I/O, training, scoring, evaluation, and compatibility fallback
```

Refactor `calc_cc_and_res_dmats` to consume the shared raw-correlation helper,
then verify that its resulting `dmat_cc` is unchanged. The relational extractor
uses the same helper directly. It must not invoke clustering code or a
commander, and must not alter class states while extracting features.

The CC helper remains the neighbourhood anchor. Channel construction stays
separate from reduction so future relationship matrices can populate the same
validated contract without reworking model I/O.

## Analysis and Validation Outputs

Analyze has one unambiguous tabular output:

```text
cavgs_quality_training.txt
```

This is the only per-dataset file required by `quality_mode=learn` and
`quality_mode=evaluate`. It contains:

- dataset identity, context, label semantics, base schema, relational schema,
  and normalization policy;
- the CC-anchor and signal-statistics provider fingerprints;
- one row per class with `class`, `manual_state`, and `hard_reject` fields;
- all fourteen base raw and normalized features in their declared order;
- the raw and normalized `signal_stats_anchor_topk_mean` columns.

Analyze writes this table directly from the in-memory base and relational
results. The learner receives a file table containing one canonical path per
dataset and validates the schema-fixed policy before fitting. Hard-rejected
rows remain in the file for diagnostics but do not enter the soft model fit.

The migration packer recognizes the old five-file layout, extracts the frozen
CC-anchor top-five signal-statistics statistic, writes the canonical table,
and can then remove the legacy text artifacts. The migrated `chunk_training5`
corpus now contains one canonical training table in each of its twelve dataset
directories.

Validation holds out complete projects or specimens and compares models on
identical hard-gate survivors. Report macro AUC, balanced accuracy, F1,
manual-good recall, manual-bad specificity, calibration, particle-weighted
effects, and worst-project behaviour. Inspect false-positive and false-negative
stacks, with particular attention to:

- the known dataset where Joe's POC rejects acceptable classes;
- overfit classes that survive Joe's hard group rule;
- coherent fuzzy-ball or ice neighbourhoods;
- genuinely good but isolated classes;
- sensitivity to low-pass limit, shift range, and `k` in a separately named
  experimental schema;
- stability when classes are added to or removed from a project.

## Defocus–Power Diagnostic

The integrated CTF-corrected particle-power versus defocus idea can also be
made label-free: compare each class's fitted trend or residual summary with
those of its rotational-correlation neighbours. Keep this outside
`corr_knn_v1` initially. Add it only if it contributes held-out information
beyond the image-quality and graph features.

## Implementation Sequence

Completed:

1. Extract and test the shared symmetric rotational-max correlation matrix,
   reproducing Joe's POC search settings without running `cluster_cavgs`.
2. Introduce the validated pairwise-relation channel contract and wrap the CC
   and FRC outputs from the same match as the first provider.
3. Expose the power-spectrum, TVD, JSD, Hellinger, and merged signal-statistics
   distances already calculated by `calc_sigstats_dmats`.
4. Analyze Joe's twelve datasets and promote a compact relational feature
   schema based on grouped empirical evidence: the
   CC-anchor top-five mean `signal_stats_distance` as `corr_knn_signal_v1`.
5. Write the self-contained `cavgs_quality_training.txt` table natively in
   analyze and provide a one-time packer for the current five-file corpus.
6. Remove the exploratory sidecar writers and obsolete relational parameter
   surface so routine analyze has one training artifact.
7. Implement version-10 artifact I/O, relational training/evaluation, automatic
   chunk application, and base-only fallback for unsupported models.
8. Pack and retrain the full `chunk_training5` corpus and promote the result as
   the built-in `chunk100mics` model.

Next:

9. Validate the trained Stage 1 model on independent datasets. Add cross-fitted
   `neighbour_quality` only if the label-free model demonstrates held-out
   value.
