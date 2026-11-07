# Multi-state solve3D: the state as part of the posterior draw

Proposal of 2026-11-03 (Hans's frame: hard state labels in solve3D, no
weights; `topk=3` kept with the pose weighting per state; the states sampled
probabilistically throughout the search from the same posterior the
single-state tables draw from), revised the same day after Hans's review, and
revised again on 2026-11-04 after the first autonomous run (`ms_posterior`,
dell2, stopped in its Phase 0 at Hans's request; what it measured is in
section 9). This revision changes how the work is done: the development is
driven by fast tests on small simulated data, a new test of a few minutes
(`multistate_workflows suite=labels`, which joins the nightly entry) is built
first, and the subset of Afan's benchmark made by the first run validates the
finished development against the truth alone, with no comparison arm of the
current code (Hans, 2026-11-04: "rather than assuming that whatever is in place
right now is of any value"). Sections 1-8 are the design, section 9 the first
run's measurements, sections 10-14 the plan of the second run (`ms_posterior2`,
dell2). Sections 1-8 were built in that run's Phases 1-2 (2026-11-05), on the
tree the batch alignment pipeline refactoring left (ruling 13), tested on small
data in Phase 3 and validated on the benchmark subset in Phase 4 (2026-11-06);
the Progress section holds the record, and the report beside this plan
(`multistate_posterior_state_draw_report.md`) the summary.

## 1. What is proposed, in one paragraph

A particle's state is one draw from its tempered posterior over state and
direction together, the same draw, primitive, temperature and annealing the
single-state tables and the 2D class draw use: one `draw_posterior_row` over
the particle's whole row (every active state's block) instead of the argmin
over states followed by a draw within the argmin's block. The label stays hard:
one state per particle per iteration, one half-map, the weighted insertion's
extras from that state's block alone with weights summing to 1 within it. The
dense fill screens and polishes shifts per state within the current per-particle
refinement budget, so every state's block carries polished entries at the same
cost as today. The neighbourhood table draws the same way over its pooled
candidates and takes its extras from the drawn state. Stages 1-2 (`shc`) stay
as they are. One tempering factor, `prob_beta`, flat at 0.3 and annealed to 1 in
the last stage's tail, governs states and directions alike, and a
full-population assignment at `prob_beta=1` before the final reconstruction
gives every active particle a label drawn at the physical factor against the
final maps. A state-draw readout shows whether these choices hold. No new
parameter of any kind.

## 2. Where the code stands

- **Dense table (`refine=prob`, solve3D's stages 3-5).** `ref_assign`
  (`simple_eul_prob_tab.f90`): the state is the particle's closest reference
  over all states (argmin, ties to the lowest index); the direction one draw
  from `exp(-beta T_i (L - L_min))` within that state's block. The weighted
  insertion's extras (`select_extras`) read the committed state from the
  assignment and come from its block, within a per-state cone (each state's
  FSC=0.5 on disk), renormalised over the committed entry and the extras:
  hard state membership, Hans's ruling of 2026-10-29 (`refine3D_policy.md`,
  "The weighted insertion").
- **The fill.** Both shift paths select globally. Under the window search
  `refine_best_refs` polishes the `projs_ns` best references *over all
  states* (`minnloc` over the whole row; `projs_ns` is the maximum of the
  per-state `calc_num2sample`). Under the polar shift map (the default)
  `process_particle_with_map` draws the coarse screen's pool from the whole
  row at the origin (`screen_pool_size` over `nrefs`) and polishes the
  `projs_ns` best of that screen. A state that trails gets few or no refined
  entries, so its block holds losses at the stored shift beside the leading
  state's polished losses (cause 4 of the regression plan of 2026-10-21, never
  tested). Changing only the polish selection would keep the starvation,
  because the screen already chose the pool.
- **The row is an approximation to the pose evidence.** Each entry of the
  dense row is the loss at *one* in-plane angle drawn from the posterior over
  the grid rotations (`gen_posterior_val`), not the loss marginalised over the
  angle. A marginal over such entries is a posterior over the constructed row,
  not the physical evidence: sampling an angle and weighing its likelihood
  again favours a sharp angular distribution over a flat one of equal
  integrated evidence (angular likelihoods (0.99, 0.01) and (0.5, 0.5) have
  equal evidence, and this procedure picks the first with probability about
  0.66, by direct enumeration; Hans's review). The in-plane marginal exists as
  a primitive (`gen_posterior_evid`, used by the 2D class draw under the shift
  map) and is not in the 3D row.
- **Neighbourhood table (`refine=prob_neigh`, stages 6-8 and
  `refine3D_states`).** `ref_assign_neigh`: `gather_state_row` chooses the
  state of the closest valid candidate over all states and gathers that
  state's candidates for the draw; `select_extras` calls `gather_state_row`
  *again* and so chooses the argmin's state on its own, and `which_pos` is the
  committed candidate's index in that gathering. Drawing the state elsewhere
  without changing the extras' gathering would take the extras from the wrong
  state, and gathering every state would admit extras from other states.
- **2D.** `ref_assign` and `ref_assign_sparse` in `simple_eul_prob_tab2D.f90`
  draw the class from the posterior over all active classes: a class in 2D is
  what a state is in 3D, and this is the template.
- **Stages 1-2.** `strategy3D_shc` walks a random order over
  `nspace x nstates` references and takes the first improvement: states and
  directions on an equal footing (ruling 5 of 2026-10-21).
- **Sampling and the final labels.** solve3D samples `nsample` (10,000 by
  default) per state, so `nsample x nstates` of the active particles per
  iteration, capped at 90 % of them (`update_frac` in
  `simple_commanders_solve3D.f90`). Only the sampled particles receive an
  iteration's draw; the others keep labels drawn at earlier iterations, against
  earlier maps, at that iteration's factor. Before the final reconstruction,
  `ensure_multistate_particle_assignments` only checks that every active
  particle was updated at least once and runs a greedy missing-update pass for
  the rest. So the annealing tail puts no more than its iterations' samples at
  the physical factor.
- **Convergence state.** The convergence object keeps the resolution history
  of the last two checks as scalars of the best populated state
  (`res05_hist(2)`, `res143_hist(2)`, `res_hist_iter(2)`), and
  `prev_iter_res` reads the best state's value back from `STATS_FILE`.
  `check_conv3D` makes a fresh statistics record each iteration and writes it,
  so a value computed at the assignment reaches the file only through the
  convergence object.
- **The regression evidence (2026-10-21/22).** On Cyril's simulated two-state
  set no build separated more than half its runs; failures either drained the
  minority state in stage 4 or 5 or froze the labels after one state locked on
  in stage 3. The regression plan's temperature on the state assignment was
  planned and never built.

## 3. The model

For a sampled particle `i` with noise temperature `T_i = pi wsqsum / (2 pftsz)`
and the row's values `L_{s,r}` (state `s`, direction `r`, each entry at its
drawn in-plane angle and its refined or stored shift, as the fill leaves them),
the tempered posterior over the pair is

    p(s, r | row_i)  proportional to  pi_s exp(-beta T_i (L_{s,r} - L_min)),

with a uniform state prior `pi_s` (no prior) and `beta = prob_beta`. Its state
marginal is `p(s | row_i) = sum_r p(s, r | row_i)`. One draw from the full row
is a draw of `s` from this marginal followed by a draw of `r` within `s`: the
implementation is the single-state draw over a longer row, and `nstates=1` is
the same code on the same random stream.

This is a posterior over the constructed row. It is exact for the row, and the
row approximates the pose evidence in two ways the readouts watch: each
entry's in-plane angle is a draw, not a marginal (section 2), and the entries of
a block mix polished and stored-shift losses (section 4.2). "The physical
posterior" at `beta = 1` below means this row's posterior at the noise
temperature; the label-recovery tests on true poses measure how far that is
from the truth.

**Why the marginal and not the best-per-state softmax** (ruling 6 of
2026-10-21). The two agree at `beta -> infinity` (both give the argmin's state)
and differ at the search's factor by an entropy term: a state whose map is
blurrier has a flatter block, so it holds more mass for the same best loss.
Under a uniform prior over directions that is the row's evidence, and it leans
against the drain. The in-plane bias of section 2 acts in the opposite
direction and is bounded by the same effective number of candidates; the
label-recovery tests show whether either matters, and the in-plane marginal is
the extension of section 6 if they do.

**Why one tempering factor.** The tempering compensates an overconfident
posterior (the temperature is the noise's, while the map at a stage is
band-limited, regularised and, early on, wrong). That applies to the state
entries as to the direction entries, and the argmin is the `beta = infinity`
limit of the state rule, the limit that locks on and drains. The tail anneals
both to `beta = 1` in the last stage (`prob_beta_end=1`, set by the controller),
and the full-population pass of section 4.5 then gives every active particle a
label drawn at `beta = 1` against the final maps.

**Symmetry breaking.** A balanced random start gives `nstates` nearly identical
mixture maps, a flat state posterior and labels re-drawn almost at random each
iteration. The separation grows from the bias the maps' small differences put
into the draws, which the next reconstruction amplifies: stochastic EM, slower
than the argmin's amplification, without its lock-on. The risk runs the other
way: at 0.3 the state posterior may be too flat to break the symmetry within a
stage's budget. Section 6 lists the fallback and section 7 the readout that
decides.

## 4. The changes

### 4.1 The dense assignment (`ref_assign`)

- Replace the argmin-over-states block and the within-block draw with one
  `draw_posterior_row` over the row's `nrefs` entries (the active states'
  blocks); the state is `ref_state(which)`, the direction `ref_proj(which)`.
  A row without a valid distance keeps reference 1, as now.
- The n_eff readout is offered the whole row; the state readout of section 7
  is computed here.
- The extras: unchanged in rule. `select_extras` reads the committed state
  from the assignment, so they come from the drawn state's block with that
  state's cone, and the weights renormalise within it. The
  refined-committed-entry rule stands.
- `nstates=1`: the same call on the same stream, bit for bit.

### 4.2 Per-state screen and polish at constant cost (ruling of 2026-11-04)

Every active state's block is built the same way, within the current
per-particle refinement budget, so that no block is starved and the cost stays
today's:

- the window path: `refine_best_refs` polishes the best
  `ceil(projs_ns / n_active)` entries of each active state's block, `projs_ns`
  as today (the maximum of the per-state `calc_num2sample`) and `n_active` the
  number of active states;
- the map path: the screen's pool is drawn per block (`minnloc` within the
  block) with today's pool size (`screen_pool_size` over the row) divided over
  the blocks, and the best `ceil(projs_ns / n_active)` of each block's screen
  are polished.

The per-particle refined capacity stays `projs_ns` (up to rounding: at most
`n_active - 1` entries more). With one state the split is the identity, so
`nstates=1` stays bit for bit. The part-file format is unchanged.

What this does not establish: each block holds fewer polished entries than the
leading state held before, and blocks mix polished and stored-shift losses
differently. The readout of the share of each block's tempered mass on
unrefined entries (section 7) bounds that; if it is not small, the candidate
rule is one polish count for every block. The first proposal gave every block
its own count; the first run measured that this would multiply the table time
by up to 2.5-2.8 at three states (section 9), hence the ruling.

### 4.3 The neighbourhood assignment and its extras (`ref_assign_neigh`)

Two gatherings, kept apart:

- for the draw, every valid candidate of every state (the pooled row); the
  draw runs once over it and the state is the drawn candidate's;
- for the extras, the candidates of the *committed* state only:
  `gather_state_row` takes the state as an argument instead of choosing it,
  and the committed candidate is identified by its store position (not its
  index in the pooled gathering) and remapped to its position in the committed
  state's row, so that `top_posterior_extras` excludes the committed entry and
  the cone is the committed state's.

`seed_fallback_if_empty` is unchanged. In `prob_neigh_mode=state` the pooled
neighbourhood is the same for every state (`build_pooled_neighborhood`) and in
`geom` one subspace serves every state, so the candidate sets are comparable.
In `prob_neigh_mode=shc|snhc` the candidates are hill-climb proposals in random
order over all states with an early exit, so the counts per state are unequal;
the draw over the store is then an acceptance step over the proposals and the
incumbent, and the state is part of that acceptance. No workflow runs those two
modes with `nstates > 1`; the unit tests cover them. The rule reaches
`refine3D_states`, which runs `prob_neigh` in state or geom mode.

### 4.4 Convergence and the readout's route

- The pose-overlap bypass for `prob`/`prob_neigh` stands.
- The annealing tail's stall test reads every populated state: per-state
  history (`res05_hist(2, nstates)`, `res143_hist(2, nstates)`, allocated at
  the first check; a state is its index, stable within a run), `prev_iter_res`
  read per state from `RESOLUTION_FSC05_STATEnn` / `RESOLUTION_STATEnn`, and the
  stall holds when no state gained: a state is "not gaining" when both its
  crossings are present in both records and neither improved; a state missing
  a value in either record is skipped; when every state is skipped the rule
  falls back to today's best-state test. The two-iteration margin of `topk > 1`
  applies per state.
- The state readout travels from the assignment to `STATS_FILE` through the
  convergence object: `prob_align_exec` / `prob_align_neigh_exec` return it, the
  refine3D strategy stores it on `conv` in both execution modes, and
  `check_conv3D` writes `STATE_NEFF_DRAW`, `STATE_NEFF_PHYS`,
  `STATE_MODE_FRAC` and the per-state `STATE_MASS_UNREFINEDnn` beside
  `STATE_OVERLAP`. These are statistics-file keys, not parameters.

### 4.5 The final labels: a full-population assignment at `beta = 1`

Before the final reconstruction of a multi-state solve3D, one `refine3D`
iteration over every active particle with `refine=prob`, `prob_beta=1` (no
tail), `update_frac=1.0`, `balance=none`, `fillin=no`, `trail_rec=no`,
`volrec=no`, `topk=1`, `maxits=1`, at the last stage's band and crop against
the last stage's maps, built as the greedy missing-update pass is today
(`run_multistate_missing_update`), which it replaces for `nstates > 1`: every
active particle is updated, so nothing is missing, and the coverage check that
follows must find nothing. Its draws are at the physical factor against the
final maps, so the final reconstruction uses labels of one provenance. Measured
cost (section 9): 1.2-1.3 late stage-5 iterations on the full benchmark, about
15 minutes on the subset.

### 4.6 What does not change

- The stage schedule: `shc` in 1-2, `prob` in 3-5, `prob_neigh` state mode from
  6; `MAXITS`; the 2,500-direction cap; the balanced stochastic sampling from
  stage 4; NU at stage 5; no stage omission or resolution exit for multi-state.
- `prob_beta` and `prob_beta_end` as the controller sets them (flat 0.3, the
  last stage's tail to 1), unless the fallback of section 6 is needed.
- Hard labels everywhere: one half-map per particle, `calc_3Drec`'s grouping by
  `(state, half)`, the trailing accounting, the final reconstruction from the
  labels of 4.5. `m_estimator=flex` stays refused on multi-state `refine3D`;
  the state weight set is `flex_pca`'s.
- Empty states: an emptied state leaves the table and does not return; no
  sampling or relabelling is driven by the maps (ruling of 2026-10-04).
- No parameter is added, internal or public (Hans, 2026-11-04: every field of
  `simple_parameters` becomes an accepted key of every program; ask before
  adding one).

## 5. Tests

The ladder the development runs on, fastest first. Section 12 says which rung
gates which step.

- **Unit tests** (seconds, the fast gate). `simple_eul_prob_tab_tester`
  (`unit_pftc_align2D3D`, sub-suite "3D probability table"):
  `check_posterior_assignment` loses "the state is the argmin's state under
  the draw" and gains: at a huge temperature a particle with a unique row
  minimum takes the argmin's state and reference; exact ties across states are
  drawn uniformly (the deterministic tie rule goes); at a finite temperature on
  a synthetic two-state row whose blocks hold known masses, the state
  frequencies over seeds within tolerance of the marginal; `nstates=1` draws
  identical to the base tree's at the same seed. `test_insert_extras_multistate`:
  every extra in the drawn state's block, with a case where the drawn state is
  not the argmin's. `test_ref_assign_neigh_posterior`: the pooled draw; a
  particle whose drawn state differs from the closest candidate's state, extras
  only from the drawn state and the committed candidate excluded (`state` and
  `geom`); the `shc`/`snhc` store with unequal candidate counts per state. The
  fill: after a two-state fill every active block carries a screened set and a
  polished entry per particle, and the total polished count is today's. The
  convergence tester: the per-state stall rule on synthetic records (one state
  gaining keeps the tail from starting; a state without an FSC is skipped; all
  skipped falls back); the readout keys in the file. New or changed OpenMP
  code: a tester at 16 threads, repeated, bit for bit against one thread.
- **The label-recovery route** (`multistate_workflows suite=labels`, minutes;
  built first, section 12): true poses, random labels, the state draw alone.
- **The multistate routes** (`states`, `cavgs`, `shc`, `consensus`, about
  5-20 minutes each at 8 threads; `hard_shc` and `shc_identical` by hand).
- **Single-state identity**: `simulated_workflow_1jxy` at one thread with the
  test environment's seed, base against new, through `test=identity` (the fills
  use dynamic scheduling and per-thread streams, so only one thread is
  bit-comparable).
- **Profiles**: `test=profile target=multistate_workflows suite=<route>
  nseeds=20` for `labels` (to set its floors) and `shc` (its floor, 6 of 20
  resolving every state, must hold; 13 of 20 was the count of 2026-10-19).
- **Validation on the subset S of Afan's benchmark** (hours): label recovery on
  true poses, then de novo solve3D, judged against the truth.

## 6. Risks, fallback and extensions

1. **Too hot for states.** The state n_eff stays near `nstates` while the
   labels stay at chance. The fallback is the tail to 1 at the end of every
   multi-state `prob` stage, not only the last (the controller sets
   `prob_beta_end=1` on stages 3-5; an existing key), which commits the labels
   at each stage boundary while keeping the exploration within the stage. A
   separate state temperature would need a parameter: if the tail does not
   suffice, the run stops and writes the question up for Hans.
2. **Cost**: the split of 4.2 keeps today's budget; the cheap variant (the
   global top-`projs_ns` with each block's mode guaranteed a polish) only if the
   readouts call for it.
3. **The row's in-plane bias** shows in the label-recovery tests as a state
   preference the truth does not support: the extension is the in-plane
   marginal per reference (`gen_posterior_evid`'s effective loss in the row for
   `nstates > 1`). It changes what the row holds, so it is Hans's decision.
4. **A state drains** when its map is garbage at any temperature: an emptied
   state's reseeding is outside this work; a population floor in the prior would
   be map-driven and is ruled out.
5. **Ambiguous particles' final labels** are draws at `beta = 1` against the
   final maps; a particle that fits two states alike lands in one of them at
   random. That is the honest hard label.

## 7. The readout

Per table assignment, beside `PROB_TAB_ASSIGN` and the n_eff line:
`>>> PROB_TAB_STATE_DRAW: n_eff over states <draw> / <beta 1>, <pct> % at the
state mode, mean max state posterior <p>`, the state n_eff being exp(entropy)
of the marginal over the same subsample of rows as the n_eff readout; and per
state `>>> PROB_TAB_STATE_MASS (state s): <frac> of the block's tempered mass on
unrefined entries`. Populations and `STATE OVERLAP` (`mi_state`) are already
printed by the convergence object. The values go to `STATS_FILE` (4.4). Read
beside per-iteration accuracy against the truth, they say whether 0.3 is too hot
(state n_eff near `nstates` while accuracy stays at chance), too cold (state
n_eff near 1 from the first iterations, the lock-on signature) or right, and
whether the unrefined mass is small enough to trust the marginal.

## 8. Decisions taken

1. The joint draw (the marginal over the constructed row), with the row's
   approximations named and measured (Hans, 2026-11-03).
2. One tempering factor for states and directions; the per-stage tail as the
   only fallback; no new parameter (Hans, 2026-11-04).
3. The per-state screen and polish at constant cost (ruling of 2026-11-04).
4. The same rule in the neighbourhood table, hence in `refine3D_states`, with
   the extras from the committed state and the `shc|snhc` modes covered.
5. The tail's stall test over every populated state.
6. Stages 1-2 and the stage schedule untouched.
7. The full-population assignment at `beta = 1` before the final
   reconstruction.
8. Development on fast tests first; validation on the subset against the truth,
   with no comparison arm of the current code (Hans, 2026-11-04).

## 9. What the first run measured (`ms_posterior`, dell2, 2026-11-03/04)

The first run's Phase 0, on the base tree fae5d8bf5 (current code). Its journal
is `/data2/multistate_bench/first_run_ms_posterior/journal_phase_0.md` on dell2.

- **The benchmark** (section 11) re-imported with true poses (PT) and zeroed
  poses (PZ), row identity checked against the label file.
- **Truth control**: maps reconstructed from the true poses and labels match
  their truth maps (above 0.9 to 12 A); with shuffled labels every map
  resembles the medium assembly (FSC about 0.95 against truth 2 at 20 A, 0.6-0.7
  against truths 1 and 3): a failed separation's maps.
- **Known-reference control** (true poses, truth maps as references, `refine3D
  refine=prob nstates=3 nspace=2500 lp=10`, two iterations): the argmin assigns
  every state correctly (accuracy 1.0000 at `prob_beta` 1 and 0.3). The state
  signal is there at 10 A.
- **The current code de novo on all 80,000 particles** (B_1, `solve3D
  nstates=3 nparts=2 nthr=8`): 14 h 23 min, failed: accuracy 0.453, ARI 0.078,
  no true state recovered. One found state locked onto the medium assembly from
  stage 4 (23,156 particles, 90 % true state 2, FSC=0.5 at 8.7 A; its final map
  docks onto truth 2 at correlation 0.96); the other two stayed 18 A mixtures
  of all three assemblies; the small assembly was never isolated; no state
  drained and the labels kept moving (state overlap 0.90-0.95). This is not
  the drain or freeze of Cyril's set: the labels move but never organise.
- **Where the time goes**: the probabilistic pre-step is 97-99 % of every prob
  stage; within the table fill the polar shift search (coarse screen and polish)
  is 79-88 % of the thread time and scoring the particle-reference pairs the
  rest. The screen computes 250-1,000 coarse shift maps per particle (the
  quantile pool triples when three near-identical states share every
  direction); the polish takes 54 references per particle in stages 3-4 and 138
  in stage 5 (the 10-degree threshold over 2,500 directions). Stage times 3 h
  (stage 3), 3 h (stage 4), 8.3 h (stage 5). This cost is the current code's
  and is outside this work.
- **The final full-population pass** (section 4.5) on B_1's finished project:
  3,595 s on all 80,000 particles, 1.2-1.3 late stage-5 iterations; about 15
  minutes on the subset. On a failed run it neither helps nor hurts (accuracy
  0.457, 10 % of labels changed).
- **The subset S** (section 11) made and checked; its 2D solution (P2S) was
  not finished when the run stopped.

# The run: `ms_posterior2` on dell2

## 10. Rulings

Fixed by the maintainer before the launch (Hans, 2026-11-03 and 2026-11-04,
and the standing rulings of 2026-10-04 to 2026-10-21 that hold):

1. State labels are hard: one state per particle per iteration, one half-map,
   no state weights anywhere in solve3D. The flex weight set is `flex_pca`'s
   and is neither read nor written.
2. `topk` stays 3; the weighted insertion's extras come from the committed
   state alone, with weights renormalised within it.
3. The state is part of the posterior draw: one draw over the row across
   states, at the one tempering factor `prob_beta` with the one annealing tail.
4. No parameter is added, internal or public. If the work seems to need one,
   the run stops and writes the question up.
5. `nstates=1` stays bit for bit identical at one thread; a change that alters
   it and cannot be made not to is a stop.
6. Stages 1-2 (`shc`), the stage schedule, budgets, the direction cap, NU at
   stage 5 and the speed changes stay.
7. Sampling and relabelling are never driven by the maps or poses; an emptied
   state's reseeding is outside this run.
8. State 0 is exclusion; a true state in state 0 is a failure. De novo routes
   are judged over several runs, never one. Bugs found are fixed in the run.
9. No floor, tolerance or timeout of an existing test is relaxed; new floors
   come only from a profile of the finished code by the floor rule (the
   observed count less three binomial standard deviations), never from a first
   run.
10. Development is driven by the fast tests (section 12); the subset validates
    the finished development against the truth; no comparison arm of the
    current code is run (Hans, 2026-11-04).
11. The `labels` route joins the default route list of `multistate_workflows`,
    so the nightly entry runs it (Hans, 2026-11-04).
12. The batch alignment pipeline refactoring
    (`doc/implementation_notes/planned/prob_align_batch_pipeline.md`, its own
    run on the dell, launched 2026-11-04) goes first: it edits the files of
    Phases 1-2 (the probability tables, their utilities and neighbourhood
    module, the pre-alignment strategy, the refine3D strategy, the
    pre-alignment commanders). This run therefore makes Phase 0 only now (the
    `labels` route and its variants, P2S, the tools, the builds); Phases 1-5
    start later on the tree that refactoring leaves, where the joint state draw
    becomes a selection rule of its selector (its section 3.3). Phase 0 edits
    only test, registry, policy and plan files, which the other run does not
    touch (Hans, 2026-11-04, after the plan was written; applied from Phase 0).
13. Phases 1-5 run on the tree the batch alignment pipeline refactoring left
    (`doc/implementation_notes/completed/prob_align_batch_pipeline.md` and its
    report; Hans, 2026-11-05, applied from Phase 1). There is no probabilistic
    pre-step any more: the matchers (`refine3D_exec`, `refine2D_exec`) fill the
    table and draw each particle's pose batch by batch, in shared memory and in
    every distributed part alike, and the distributed master merges the parts'
    readout records in part order. Read the design of sections 2-4 through it:
    - the dense joint draw of 4.1 lands in the per-row selection `assign_rows`
      of `simple_eul_prob_tab.f90` (where `ref_assign` stood; `select_extras`
      inside it, `assign_candidate` commits), the neighbourhood draw of 4.3 in
      `assign_rows_neigh` of `simple_eul_prob_tab_neigh.f90`; the batch frame
      around them (`open_assignment`, `assign_batch`, `close_assignment`, the
      kept batch rows) stays as it is;
    - every draw is keyed (`rnd_keyed`; the assignment's family `ASSIGN`): the
      joint draw uses the keyed variate the row selection already draws, so a
      draw depends on the seed, the pass, the iteration and the particle only;
    - the state readout of section 7 travels in the parts' readout records
      (`prob_readout` in `simple_eul_prob_tab_utils.f90`), summed by the master
      in part order (sums and counts, never means of means) as the `n_eff` and
      row-mode readouts are; `prob_align_exec` and its commanders are gone, so
      section 4.4's route to the statistics file runs from the matcher (shared
      memory) and the master's merge (distributed) to the refine3D strategy;
    - ruling 5 reads: `nstates=1` bit for bit identical at one thread and at
      eight in shared memory, and at the distributed cases' layout; the
      identity's base is this tree's HEAD (ea4c1c20a), exported again with
      `git archive`. A single-state comparison is made on run directories of
      the same build and layout (the shift search's simplex start still draws
      from the intrinsic stream);
    - known at this HEAD and not this run's to fix: `simulated_workflow_1jxy`
      stops at solve3D stage 3, `simulated_workflow_6vxx` fails its
      final-volume correlation floor;
    - the `labels` reports of Phase 0 are from the tree before the refactoring;
      Phase 1 reruns them once on this tree before any change.

## 11. The data

**The beta-galactosidase states of the multi-state fixture**
(`simple_test_state_fixture.f90`): box 96 at 2.2 A, mask 200 A; compositional
states (the tetramer, the tetramer less chain D, the A-B dimer; populations 50,
30, 20 %), conformational states (chain D rotated by 10 degrees) and identical
states; modes `search` (no poses, no states), `known` (the true poses, no
states) and `views`. Built and cached per build by the fixture itself
(`test=state_fixture` takes 7 s).

**Afan's Ribosembly benchmark** on dell2 at `/data2/multistate_bench`
(read-only): 80,000 simulated particles, 128 x 128 at 3.0 A, in one shuffled
stack `particles.mrcs`; true labels `ground_truth_labels.txt` (1-based, stack
order): state 1 the smallest assembly (8C9C, 24,000), state 2 the medium (8C96,
48,000), state 3 the largest (8C8X, 8,000). `reconstructed_volumes/` holds the
truth maps in the frame of the true poses. SNR 0.01. `mskdiam=300` (the maps
extend to 211, 224 and 267 A). The original project files point at a path that
does not exist on dell2 and are not used.

**The subset S** at `/data2/multistate_bench/subset20k` (read-only; made by the
first run, ruling of 2026-11-04): 20,000 particles with the benchmark's
proportions, 6,000 / 12,000 / 2,000, drawn at random without replacement within
each true state (`make_subset.py`, seed 20261104). `subset_indices.txt` holds the
drawn stack indices in ascending order, which is the subset projects' row order;
`subset_truth.txt` the true state per row; `subset_flags.txt` the 0/1 flags over
the full stack. `PZS/PZS.simple` (zeroed poses) and `PTS/PTS.simple` (true poses
and true labels) were made by `selection infile=<flags> prune=yes` on the
re-imported projects; pruning keeps the stack file and the surviving rows'
order, and every row keeps its physical image index (`indstk`) and its original
particle index (`pind`), which is what the iteration label files
(`states_iter<NNN>.txt`) print. Checked: every row equals its source row in
defocus, astigmatism, angles, shifts and half, 0 mismatches in 20,000.
solve3D samples 18,000 particles per iteration on S (update fraction 0.9, above
the full-sampling switch), so a de novo run takes roughly 9 hours on dell2.

**Tools** at `/data2/multistate_bench/tools`: `acc3.py` (accuracy under the
best matching of found to true states, recall per true state, ARI, confusion,
populations, verdict; `--rows subset_indices.txt` maps a `states_iter` file's
`pind` to S's rows), `run_metrics3.py`, `check_subset_rows.py`, `test_acc3.py`
(its tests, 12 of 12 passing), `dock_eval.sh` (both-hand docking onto the truth
maps; its run root is hard-coded to the first run: adapt a copy), and the first
run's scripts under `scripts_ms_posterior/` (adapt copies; never edit in place).

## 12. How the work is tested and judged

**The `labels` route** (built in Phase 0). A sub-suite of
`multistate_workflows` (`simple_commanders_test_highlevel_msworkflows.f90`)
that tests the state assignment alone: the compositional beta-galactosidase
fixture in `known` mode (true poses, 3,000 particles, SNR 1, populations 50 /
30 / 20 %), labels drawn at random with balanced populations by refine3D's own
fresh multi-state start (or, if refine3D needs starting maps, the same uniform
labelling followed by `reconstruct3D nstates=3` to make `vol1-3`), then
`refine3D refine=prob nstates=3` with `objfun=euclid ml_reg=yes` and a fixed
band, direction count and iteration budget chosen so the route takes at most
about five minutes at 8 threads, and `prob_beta_end=1` so the last iterations
anneal as solve3D's last stage does. It is judged by the suite's existing truth
metrics (accuracy, ARI, recall per true state, population error, discrimination
margin, the shuffled-label and duplicated-map controls) and reports per
iteration accuracy beside the state readout of section 7 once that exists. It
joins the default route list (`case('')`), so the CTest entry and the nightly run
it; its labels are reported, not gated, until Phase 3's profile sets floors by
the floor rule. Variants by hand: `labels_identical` (identical states: nothing
to recover; the populations and the maps' mutual correlations are reported) and
`hard_labels` (the conformational pair, 80 / 20 %, SNR 0.5). The suite's header
comment, the registry row's comment (the entry's time) and the test policy's
description of the entry are updated with it.

**Gates per step** (Phases 1-2): the fast gate and the touched unit suites; the
single-state identity; `labels` (all three variants), `states` and `cavgs`
once each at the end of each phase; `consensus` at the end of Phase 2 (it runs
`refine3D_states`, which the neighbourhood change reaches).

**Robustness** (Phase 3): `shc`, `hard_shc`, `shc_identical` and the full entry
once each; the profiles of `labels` and `shc` (section 5); `test=nightly`.

**Validation on S** (Phase 4), against the truth only:

- *Label recovery on true poses*: PTS with balanced random labels (drawn once,
  seed recorded; the starting maps by `reconstruct3D nstates=3` if refine3D
  does not make them), then `refine3D refine=prob nstates=3 pgrp=c1 mskdiam=300
  lp=10 nspace=2500 update_frac=0.9 balance=none maxits=8 objfun=euclid
  ml_reg=yes prob_beta_end=1 nthr=16 print_states=yes`, two runs. Read per
  iteration: accuracy, recall per state and the readout.
- *Known-reference control on S*: the same command on PTS with the true labels
  and the truth maps as `vol1-3`, two iterations, at `prob_beta` 1 and 0.3: the
  ceiling of the assignment at this band on S.
- *De novo*: three runs of `simple_exec prg=solve3D pgrp=c1 mskdiam=300
  nstates=3 nparts=2 nthr=8 print_states=yes projfile=<its copy of P2S>`
  (P2S: `solve2D ncls=30 mskdiam=300 nthr=16` on a copy of PZS, made in Phase
  0), no `SIMPLE_SEED`, one at a time on the machine.
- *Verdict per run*: **separated** when every true state is recovered (its
  matched found state holds at least 0.75 of its particles) and the accuracy is
  at least 0.90; **partial** when the accuracy is at least 0.75 or exactly two
  true states are recovered; **failed** otherwise (`acc3.py` defaults). If the
  known-reference control at `prob_beta=1` keeps less than 0.95, both accuracy
  thresholds drop by the shortfall. The thresholds are fixed in the journal
  before the first de novo run and never moved.
- *Recorded per run*: wall time per stage and the load; the labels after every
  stage and iteration; the log lines of the stages, the assignment
  (`PROB_TAB_ASSIGN`, `POSTERIOR DRAW N_EFF`, `PROB_TAB_STATE_DRAW`,
  `PROB_TAB_STATE_MASS`, `PROB_TAB_INSERT`), the per-state populations,
  overlaps and resolutions, the tail; `STATS_FILE` per iteration; the final maps
  docked onto the truth maps in both hands; the signature (where the separation
  forms, whether the small state drains, whether labels freeze).

## 13. Phases

### Phase 0: builds, the `labels` route, the subset's 2D solution

1. Build the tree in `repo/build` (`./compile_clean.sh`, Release, fast gate) and
   the base HEAD exported with `git archive` into `scratch/keep/src_base`, built
   in `scratch/keep/build_base` (for the single-state identity only).
2. The `labels` route and its two variants (section 12), joined to the default
   route list, with the suite's header, the registry comment and the test
   policy's description updated. Run `labels`, `labels_identical` and
   `hard_labels` once each on the unchanged code and record what they report:
   this is the route working, not a bar the development must beat.
3. P2S: `solve2D ncls=30 mskdiam=300 nthr=16` on a copy of PZS (base or tree:
   the 2D path is unchanged), alone on the machine, timed; record the particles
   it leaves inactive.
4. Copy and adapt the tools; check `test_acc3.py` passes.

Exit: both builds pass the fast gate; the route runs in the default list and by
hand, in at most about five minutes at 8 threads; its first reports are in
Progress; P2S exists.

### Phase 1: the dense joint draw, its readout, the tests, the identity

1. Section 4.1 in `ref_assign`; the state readout type in
   `simple_eul_prob_tab_utils.f90` and its route to `STATS_FILE` (section 4.4,
   the readout part); the dense table's unit tests of section 5.
2. The single-state identity, before anything else changes.
3. `labels` and its variants, iterated as needed (minutes each). If the
   readout says too hot, try the fallback of section 6 on `labels` here, not
   later. `states` and `cavgs` once at the end.

Exit: the fast gate and the unit suites pass; the identity holds; the `labels`
reports and the readout are in Progress with their reading (does the joint draw
recover the states from random labels on true poses, and what do the state
n_eff and the unrefined mass say).

### Phase 2: the fill's split, the neighbourhood table, the per-state convergence, the final pass

1. Section 4.2 in both fill paths, with the table's fill timers compared to
   Phase 1's on `labels` (same thread count, nothing else running).
2. Section 4.3 with the neighbourhood tests of section 5.
3. Section 4.4's per-state history and stall rule with its tester.
4. Section 4.5, the final pass, replacing the missing-update pass for
   `nstates > 1` in `simple_commanders_solve3D.f90`; the single-state and
   class-average routes untouched.
5. The identity again, the fast gate, the unit suites; `labels` (all
   variants), `states`, `cavgs` and `consensus` once each.

Exit: all pass; the identity holds; the fill's cost is within a tenth of
Phase 1's; the reports are in Progress.

### Phase 3: robustness on small data

1. `shc`, `hard_shc` and `shc_identical` once each, and the full
   `multistate_workflows` entry (`simple_test_exec test=multistate_workflows
   nthr=8`, the default route list with `labels`).
2. `test=profile target=multistate_workflows suite=labels nseeds=20` and its
   floors set in the suite by the floor rule (ruling 9), so the nightly gates
   it; `test=profile target=multistate_workflows suite=shc nseeds=20`, whose
   floor (6 of 20 resolving every state) must hold, its count reported beside
   the 13 of 20 of 2026-10-19 (the code has changed since: not paired).
3. `simple_test_exec test=nightly ncunits=16`.

Exit: the entry and the nightly pass with `labels` gated; the `shc` floor holds;
the profiles are in Progress.

### Phase 4: validation on the subset

The runs of section 12 in this order: the known-reference control on S, the two
label-recovery runs, then the three de novo runs, one at a time on the machine.
Each run's verdict and signature go into Progress as they finish, with a table
per run: accuracy and recall per true state after every stage, iterations per
stage, populations, the readout at the first, middle and last iteration of each
prob stage, the final maps against the truth.

Exit: the five runs and the control are made and judged; the reading (what the
joint draw does on the benchmark, and what is left) is in Progress.

### Phase 5: documents and report

1. The policy documents: `solve3D_policy.md` section 6 and the final-labels
   contract, `refine3D_policy.md`'s assignment sections (dense and
   neighbourhood), `refine3D_states_policy.md` where the pooled draw changes its
   description; skill text (`.claude/skills/simple-refine3d/references/bayesian-model.md`,
   `simple-refine3d` SKILL.md) proposed in the report, not edited.
2. The report beside the completed plan; the plan moved to `completed/`.

Exit: the documents are updated, the report is written, the plan is moved.

## 14. Common rules

- Nothing is committed; git is read-only; the base for the identity is a `git
  archive` export.
- Fix as you go: a defect found is fixed when found and recorded with its
  evidence.
- Thresholds are fixed before the runs they judge, never loosened; no run is
  repeated until it passes.
- New or changed OpenMP code gets a tester at 16 threads, repeated, compared
  bit for bit with one thread; no string temporaries inside parallel regions;
  work arrays sized by the team in use.
- No `goto`; THROW macros never continued with `//&`; no `.and.` short circuit
  on allocatable arrays; `use ..., only:`; no large type as a plain component
  (`doc/policies/compile_time_policy.md`); no new parameter (ruling 4).
- Unset `SIMPLE_SEED` for every benchmark run on S; the unit tests, the routes,
  the identity and the profiles use the seeds their test environment sets.
- Bulky outputs are deleted once their numbers are in the journal;
  `scratch/keep/runs/` keeps the labels, logs, statistics files and final maps
  of every validation run.
- Write prose for people: Progress entries and the report read without other
  documents, with SIMPLE terms explained and acronyms spelled out on first
  use; name things by what they do.

### File table

| file | why | phase |
| --- | --- | --- |
| `doc/implementation_notes/planned/multistate_posterior_state_draw.md` | this plan and its progress | each |
| `src/main/commanders/test/simple_commanders_test_highlevel_msworkflows.f90` | the `labels` route and its variants (0); their reports beside the readout (1, 2); their floors from the profile (3) | 0 1 2 3 |
| `src/main/commanders/test/simple_test_registry.f90` | the `multistate_workflows` row's comment (its time with the new route) | 0 3 |
| `doc/policies/test_environment_policy.md` | the entry's description and its routes | 0 3 5 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` | the joint draw and the state readout (1); the per-state screen and polish, the unrefined-mass readout (2) | 1 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` | the state readout type | 1 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_tester.f90` | the tables' tests | 1 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90` | the pooled draw; the extras gathered from the committed state | 2 |
| `src/main/commanders/test/simple_polar_matching_bench.f90` | the frozen-input fill harness records the fill's refined-reference budget, for the per-state polish check (`test=polar_matching_bench suite=fills`; the fill's thread and part identity at 16 threads) | 2 |
| `src/main/commanders/test/simple_commanders_test_highlevel_matchbench.f90` | the harness checks that every active state's block carries polished entries per particle within the budget | 2 |
| `src/main/strategies/search/simple_strategy3D_matcher.f90` | the matcher hands the dense table's readout record to its caller (shared memory; ruling 13: the pre-step and its strategy are gone) | 1 2 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | the readout stored on the convergence object in both execution modes (1); the distributed master's multi-state label start gated on the first iteration, as the shared-memory strategy's, a defect the final assignment exposed (2) | 1 2 |
| `src/main/simple_convergence.f90` | the readout keys (1); the per-state resolution history and the all-state stall test (2) | 1 2 |
| `src/main/simple_convergence_tester.f90` | tests of the per-state stall rule | 2 |
| `src/main/commanders/simple/simple_commanders_solve3D.f90` | the full-population assignment at prob_beta=1 before the final reconstruction | 2 |
| `src/main/solve/simple_solve3D_controller.f90` | the per-stage tail, only if the fallback of section 6 is needed | 1 2 |
| `doc/policies/3D/solve3D_policy.md` | section 6 and the final-labels contract | 5 |
| `doc/policies/3D/refine3D_policy.md` | the assignment sections of the dense and neighbourhood tables | 5 |
| `doc/policies/heterogeneity/refine3D_states_policy.md` | the pooled neighbourhood draw | 5 |
| `doc/implementation_notes/completed/multistate_posterior_state_draw_report.md` | the report | 5 |
| `doc/implementation_notes/completed/multistate_posterior_state_draw.md` | the plan, moved when the run completes | 5 |
| `doc/implementation_notes/completed/single_state_speed_week_closeout.md` | its link to this plan follows the plan's move to `completed/` | 5 |

Other files a fix needs are added to this table, with the reason, by the phase
that finds the defect.

## Progress

(one section per phase, written by the run)

### Phase 0 (2026-11-04, dell2): builds, the `labels` route, P2S, the tools

Under ruling 12 this is the only phase of the run for now; Phases 1-5 wait for
the batch alignment pipeline refactoring. The run directory is
`/home/hael/agent_runs/ms_posterior2` (below `R`).

**Files changed**: `simple_commanders_test_highlevel_msworkflows.f90` (the
route and its variants), `simple_test_registry.f90` (the comment of the
`multistate_workflows` row), `doc/policies/test_environment_policy.md` (the
entry's description), this plan (ruling 12, this section), and the new
report `doc/implementation_notes/completed/multistate_posterior_state_draw_report.md`
(this phase ends the run for now).

**Builds.** How SIMPLE builds on dell2: Ubuntu's gfortran 15.2 and CMake 4.2,
no modules, no conda; the login environment points `SIMPLE_PATH` and `PATH` at
`~/src/SIMPLE/build`, so every run of this phase sources an environment script
that points them at the build it means (`R/scratch/keep/scripts/env_{new,base}.sh`).
The tree built with `./compile_clean.sh` in 2 minutes and passed the fast gate
(15 of 15 entries, 6.2 s; `R/scratch/logs/build_tree.log`). The base HEAD
22a643642, exported with `git archive` into `R/scratch/keep/src_base` and built
in `R/scratch/keep/build_base` (Release, tests on), passed its fast gate too
(15 of 15; `R/scratch/logs/build_base.log`). After the route was added the
tree was rebuilt incrementally without compiler warnings in the changed file,
passed the fast gate again (15 of 15) and was installed.

**The `labels` route.** How it is built, and why:

- *The start.* `refine3D` draws the balanced random labels of a fresh
  multi-state start itself (`seed_multistate_startup_labels`, the uniform
  labelling of `gen_labelling`) in both execution modes, but only the
  distributed master (`nparts` given) also reconstructs the starting maps; the
  shared-memory strategy, which the other routes' `solve3D` runs use, stops
  without `vol1..N`. So the route draws the same uniform labelling itself (from
  the fixture's seed plus 2, so a profile's seed offsets give new labels),
  reconstructs the starting maps with `reconstruct3D nstates=3 objfun=cc`, and
  names them as the distributed master's fresh start names its own
  (`startvol_stateNN.mrc` with their halves and unfiltered half copies). With
  three state bins already in the project, `refine3D` keeps the drawn labels.
- *The refinement.* `refine3D refine=prob nstates=3 objfun=euclid ml_reg=yes
  update_frac=0.9 balance=none center=no bfac=0 prob_beta_end=1
  print_states=yes`, the keys solve3D's multi-state probabilistic stages set
  when there is no 2D solution, at a fixed band of 8 Å (`lp`; a given `lp`
  holds for every iteration outside the nonuniform filter mode), 1,000
  directions (solve3D's stages 3-4) and 10 iterations (`maxits`): five flat at
  `prob_beta`=0.3 and the annealing tail of five to 1 when the resolution has
  not stalled first. No new parameter: these are the test's constants
  (`LABELS_NSPACE`, `LABELS_LP`, `LABELS_MAXITS`, `LABELS_UPDATE_FRAC`).
- *The judgement.* The path check (one `prob` block over 1,000 directions,
  three states every iteration); a table `LABELS PER ITERATION` with the
  accuracy (under the best matching of found to true states), the adjusted Rand
  index, the smallest recall and the populations of the drawn start and of
  every iteration (from `states_iter<NNN>.txt`); the suite's `evaluate` with the
  final maps copied to `rec_final_stateNN.mrc` (labels, maps docked onto the
  truths, discrimination, poses, the shuffled-label and duplicated-map
  controls); and the populations and mutual correlations of the final maps.
  The labels start at random, so whether every true state is resolved is a
  `success` row, for now with a floor of 0; the labels and poses are reported,
  not gated, until Phase 3's profile sets floors by the floor rule.
- *Where it runs.* First in the default route list (`case('')`), so the CTest
  entry and the nightly run it; by hand as `suite=labels`, and the variants
  `suite=labels_identical` and `suite=hard_labels` through the suite's existing
  naming of identical and conformational variants.

**First reports on the unchanged code** (the final binary of the phase, eight
threads, alone on the machine; information about the route, not a bar).
Accuracy per iteration, the drawn start first, `|` where the annealing tail
starts:

| route | wall | accuracy per iteration | final ARI / smallest recall | populations (found) | maps |
| --- | --- | --- | --- | --- | --- |
| `labels` | 3 min 44 s | 0.34, 0.42, 0.61, 0.70, 0.75, 0.79, \| 0.81, 0.85, 0.92, 0.99, 0.998 | 0.993 / 0.992 | 0.20 / 0.30 / 0.50 | 0.97-0.99 to truth, margin 0.45, mutual correlations 0.70-0.85 |
| `labels_identical` | 3 min 18 s | 0.34, 0.34, 0.35, 0.37, \| 0.38, 0.39, 0.40, 0.41, 0.42 (null 99.9th percentile 0.415) | 0.002 / 0.035 | 0.04 / 0.44 / 0.52 | mutual correlations 0.93-0.98 |
| `hard_labels` | 2 min 50 s | 0.50, 0.51, 0.53, 0.56, 0.60, \| 0.63, 0.64, 0.66, 0.67, 0.68 (one state for all: 0.80) | 0.073 / 0.453 | 0.70 / 0.30 | mutual correlation 0.986 |

ARI is the adjusted Rand index (0 for labels unrelated to the truth, 1 for the
truth up to renaming); the smallest recall is the smallest fraction of a true
state's particles that its matched found state holds. A first run of each on
the binary before two small review edits (the starting maps' provenance file
now renamed with them, as the distributed master does, and comment text) gave
the same picture: 3 min 52 s and final accuracy 0.996 for
`labels`, 3 min 23 s and a state drained to 312 for `labels_identical`, 2 min
57 s and accuracy 0.66 for `hard_labels`; eight-thread runs are not
bit-reproducible.

What they say about the current assignment (the argmin over states): on the
compositional states it recovers the labels from random ones on true poses,
but not by the flat iterations alone: the middle state (the tetramer less chain
D, 900 particles) is starved to 220-320 particles in iterations 3-6 (smallest
recall 0.18) and is recovered only in the annealing tail. On identical states
it drains one state, to 116 of its 1,000 particles in the last iteration,
while the maps stay alike. The conformational pair is not separated. The
posterior over directions is sharp on this data (effective number of
candidates 1.00-1.05 at the search's factor 0.3), so the joint draw of Phase 1
will be close to the argmin here unless the state marginal is flatter than the
direction posterior. Evidence:
`R/scratch/keep/runs/phase0_{labels,labels_identical,hard_labels}_t{1,2}/`
(`route.log`, `metrics.tsv`).

**The default entry** with `labels` (`ctest -R '^multistate_workflows$'` in
`repo/build`, which runs `simple_test_exec test=multistate_workflows nthr=8`,
alone on the machine, on the binary before the two review edits): 2,517 s
(42 minutes; on 2026-11-02 the four routes took 54 minutes here at an older
commit). `labels` ran first and passed in 3 min 51 s (accuracy 0.996); then
`cavgs` (7.5 min, every state resolved), `shc` (16 min; it missed resolving
every state, which is a success row held only over a 20-seed profile, so not
a failure), `states` (3.5 min, every check passed) and `consensus` (11 min).
The entry **failed in the `consensus` route**: its `refine3D_states` march
emptied one state, merging the tetramer with the tetramer less chain D
(confusion rows [1500 891 0; 0 0 0; 0 9 600]), so the checks "every true
state resolved" and "no estimated state empty" failed (its accuracy 0.700 and
largest population error 0.300 sat exactly on their floor and ceiling).

This failure is not caused by the change, and is recorded rather than fixed:

- The `consensus` route's code path is untouched by the diff, and at one
  thread, where a run is deterministic under the test seed, it gives
  byte-identical `metrics.tsv` on the base build and on this tree (34 rows;
  both pass, accuracy 0.891, smallest recall 0.636; 54 minutes each; the logs
  differ only in paths, timestamps, the base build's first fixture build and
  the commit string, which the base export does not carry).
- At eight threads, runs are not reproducible, and the route is borderline on
  the current base: one more run on each build passed with the same partial
  merge (about 330-360 of true state 2's 900 particles in the tetramer's
  state; smallest recall 0.611 on the base, 0.598 here), against the full
  merge of the failing run. The route resolved every state in all 40 runs of
  its profile of 2026-10-19 and passed the nightly here on 2026-11-02 (at
  ef5fe84d5); 18 commits since (the polar shift search, the dynamic
  centring) changed the multi-state path. Making the route pass would need
  either a production fix outside this phase or relaxing its floors, which no
  ruling authorizes: it is left for the maintainer (the report says so).

Evidence: `R/scratch/keep/runs/phase0_ctest_default/` (the CTest log and
`metrics.tsv`), `R/scratch/keep/runs/phase0_consensus_{base,new}_{t1,t8}/`.

**P2S** (the subset's 2D solution): `solve2D ncls=30 mskdiam=300 nthr=16
mkdir=no` on a copy of PZS in `R/scratch/keep/projects/P2S/P2S.simple`, base
build, `SIMPLE_SEED` unset, alone on the machine: 8 min 33 s, peak memory
10.9 GB, six refine2D stages and 30 iterations, the final class resolution
(Fourier ring correlation 0.5) 22.6 Å for half the particles. It leaves no
particle inactive (all 20,000 in state 1 in both particle segments, all 30
classes in state 1). Its out segment names `cavgs_iter030.mrc` and `frcs.bin`
by absolute path inside its directory, so P2S is never modified: a de novo
run copies the directory. The earlier iterations' class averages were
deleted (29 MB left). Log: `R/scratch/keep/projects/P2S/solve2D.log`.

**The tools** are copied to `R/scratch/keep/tools`: `test_acc3.py` pointed at
the subset's files under `/data2/multistate_bench/subset20k` (the copy's
relative path to them did not exist) passes 12 of 12; `dock_eval.sh`'s run root
is this run's.

**Exit items.** Both builds pass the fast gate (above). The route runs in the
default list (it passed inside the CTest entry; the entry itself failed in the
`consensus` route, for the reason above, independent of this phase) and by
hand, in 3 min 44-52 s at eight threads (the variants 3 min 18-23 s and 2 min
50-57 s). Its first reports are above. P2S exists.

### Phase 1 (2026-11-05, dell2): the dense joint draw, its readout, the tests, the identity

On the tree the batch alignment pipeline refactoring left (ruling 13; base HEAD
ea4c1c20a, exported again into `R/scratch/keep/src_base` and built in
`R/scratch/keep/build_base`; both builds passed the fast gate, 15 of 15).

**Files changed**: `simple_eul_prob_tab.f90` (the joint draw, the state
readout, the kept record), `simple_eul_prob_tab_utils.f90` (the state readout
type, its place in the parts' readout record), `simple_eul_prob_tab_tester.f90`,
`simple_strategy3D_matcher.f90` (the record handed to the caller),
`simple_refine3D_strategy.f90` (the record stored on the convergence object in
both execution modes), `simple_convergence.f90` and its tester (the statistics
file keys), `simple_commanders_test_highlevel_msworkflows.f90` (the `labels`
route prints the readout beside its accuracy), this plan (ruling 13, the file
table: the rows of the retired pre-step strategy and commanders removed, the 3D
matcher added).

**What was built.**

- *The joint draw.* In the per-row selection `assign_rows` a particle's state
  and direction are one draw from the tempered posterior over the whole row
  (every active state's block), with the keyed variate the selection already
  drew (family `ASSIGN`, ordinal 0). The row's mode over all states stays the
  readouts' mode and the fallback of a row without a valid distance (its first
  reference). The extras of the weighted insertion come from the drawn state's
  block (`select_extras` reads the committed state). With one active state the
  row is the block, so the call, the variate and every output are the old ones.
- *The state readout* (`prob_state_readout`): over the rows the n_eff readout
  subsamples (the same rows whatever the number of parts), the effective number
  of states exp(entropy) of the row's state marginal at the draw's factor and at
  1, the largest state posterior and, per state, the share of the block's
  tempered mass on entries the fill did not refine; over every row, the draws
  in the row mode's state. It is carried in the parts' readout record (written,
  read and summed in part order as sums and counts) and printed after
  `PROB_TAB_ASSIGN` when more than one state is active:
  `>>> PROB_TAB_STATE_DRAW: n_eff over states <draw> at beta <b> / <phys> at beta 1.00, <pct>% at the state mode, mean max state posterior <p> (<n> rows)`
  and per state `>>> PROB_TAB_STATE_MASS (state s): <share> of the block's tempered mass on unrefined entries`.
- *To the statistics file.* The table keeps its last record; `refine3D_exec`
  returns it (shared memory) and the distributed master's merge returns the
  parts' sum; the refine3D strategy stores it on the convergence object, whose
  next check writes `STATE_NEFF_DRAW`, `STATE_NEFF_PHYS`, `STATE_MODE_FRAC` and
  `STATE_MASS_UNREFINEDnn` once (from as many rows as the log line needs). The
  neighbourhood table's readout is Phase 2's.

**Unit tests** (`unit_pftc_align2D3D` "3D probability table", 431 of 431;
`unit_cart_align3D` "pose statistics", 138 of 138; logs
`R/scratch/runs/unit/ut_prob2.log`, `ut_conv2.log`): at a huge temperature a
particle with a unique row minimum takes it; exact ties across states land on
the tied references only and each about half the time over 200 seeds (the
deterministic tie rule is gone); at temperature 2 some particles leave the
argmin's state; with one state the draw is the previous rule's draw within the
argmin's block with the same keyed variate, draw for draw, and with several it
is not; 2,000 identical two-state rows whose blocks hold known masses draw state
1 with the frequency of its marginal (0.429; the argmin rule would always take
state 1) and give the readout's exact values (the state n_eff, the largest
posterior, the mode fraction, the unrefined shares 0.269 and 0.667); the same
rows over two parts give records that add to the whole's, and a record
round-trips its file; in the extras test some particles draw a state other than
their argmin's and still get extras, all from the drawn state's block; the
convergence object writes the four keys once (not for a state without rows) and
not again without a new readout. No OpenMP region was added or changed.

**The single-state identity** (ruling 5 as ruling 13 reads it), base build
against this tree, same seed and layout:

- a small single-state `refine3D refine=prob` (1,500 simulated tetramer
  particles with true poses, 500 directions, 8 Å, three iterations;
  `R/scratch/keep/scripts/identity_ss_*.sh`) at one thread and at eight in
  shared memory and at two parts of one thread: every map (merged, halves,
  unfiltered, parts), FSC, sigma2 state, reprojection model and statistics file
  byte-identical, the particle records identical line for line, the logs'
  readout lines identical; the only differing files carry the run directory's
  path or process numbers, and differ between two runs of this tree at eight
  threads too (`R/scratch/keep/identity/phase1/small_refine3D/`). The check was
  repeated on the final code with the same result
  (`R/scratch/keep/identity/phase1/small_refine3D_final/`);
- `simulated_workflow_6vxx` at one thread (its solve3D runs 71 single-state
  probabilistic assignments): `test=identity`, 181 declared outputs the same, 0
  differing (`R/scratch/keep/identity/phase1/6vxx/`). `simulated_workflow_1jxy`
  was not used: at this HEAD it stops at solve3D stage 3, before any
  probabilistic stage, so it cannot exercise the change.

**The `labels` route.** First on the unchanged pipeline tree (the comparison arm
of ruling 13), then on the joint draw; eight threads, alone on the machine;
accuracy per iteration with the drawn start first and `|` where the annealing
tail starts:

| route | tree | wall | accuracy per iteration | final ARI / smallest recall | populations |
| --- | --- | --- | --- | --- | --- |
| `labels` | unchanged | 3 min 33 s | 0.34, 0.42, 0.60, 0.70, 0.77, 0.83, \| 0.94, 0.994, 0.999, 0.999, 0.999 | 0.998 / 0.998 | 0.20 / 0.30 / 0.50 |
| `labels` | joint draw | 3 min 41 s | 0.34, 0.42, 0.60, 0.70, 0.77, 0.83, \| 0.93, 0.990, 0.998, 0.999, 0.999 | 0.997 / 0.997 | 0.20 / 0.30 / 0.50 |
| `labels_identical` | unchanged | 3 min 48 s | 0.34 … 0.42 (no stall, ten iterations) | 0.003 / 0.032 | 0.03 / 0.46 / 0.50 |
| `labels_identical` | joint draw | 3 min 14 s | 0.34 … 0.42 (stall at 3, eight iterations) | 0.002 / 0.033 | 0.04 / 0.45 / 0.51 |
| `hard_labels` | unchanged | 3 min 04 s | 0.50, 0.51, 0.53, 0.56, 0.59, 0.62, \| 0.64 … 0.67 | 0.063 / 0.417 | 0.71 / 0.29 |
| `hard_labels` | joint draw | 2 min 34 s | 0.50, 0.51, 0.54, 0.57, 0.60, \| 0.64, 0.67, 0.68, 0.69 | 0.081 / 0.450 | 0.71 / 0.29 |

The state readout of the joint draw, per iteration (state n_eff at 0.3 / at 1,
draws at the state mode): `labels` 1.29 / 1.08, 91 %; 1.06 / 1.01, 98 %; 1.02,
99.4 %; from iteration 4 1.00, 99.8-100 %. `labels_identical` 1.55 / 1.16, 84 %;
1.13, 96 %; then 1.01-1.08, 98-99.7 %. `hard_labels` 1.38 / 1.13, 87 %; 1.12,
96 %; then 1.01-1.09, 97-99.9 %. The unrefined share of a block's tempered
mass stays at most 0.04 in the flat iterations and grows in the annealing tail
(`labels`' last iteration 0.64 / 0.40 / 0.25), as the refined-reference budget
shrinks with the angular moves and the fill's polish, chosen over the whole row,
concentrates on the drawn state's block. Two runs of each route on the final
code gave the same numbers to the last digit (the eight-thread route is
deterministic at a fixed seed). Evidence:
`R/scratch/keep/runs/p1base_*`, `p1try1_*`, `p1try2_*` (`route.log`,
`metrics.tsv`).

**Reading.** On true poses the joint draw recovers the compositional states
from random labels, as the argmin did, along the same path: the middle state
dips to about 300 particles in iterations 3-4 and comes back as the maps
separate. The state posterior at `prob_beta` 0.3 is almost as sharp as the
argmin on this data from the second iteration on (state n_eff near 1, 98-100 %
of the draws in the row mode's state), and at the first iteration (maps of
random labels) about 1.3-1.6: section 7's "too cold" reading, the lock-on
signature, not "too hot". The per-stage tail of section 6 is a remedy for too
hot and is not indicated; nothing in Phase 1 calls for a state temperature (a
parameter) because the states are recovered. The draw does not stop the drain
on identical states (a state falls to about 110 of 1,000) and does not separate
the conformational pair (accuracy 0.69, below the 0.80 of one state for all);
both are as before. The unrefined mass is small while the budget is wide and
large in the tail: the per-state screen and polish of Phase 2 is where that
goes.

**Known-reference run on the fixture** (true poses, the truth maps as the
references, the route's settings, two iterations;
`R/scratch/keep/scripts/kref_fixture.sh`, evidence
`R/scratch/keep/runs/p1_kref_fixture/`). Corrected in Phase 2: the import of
the particles resets every state to 1, so `refine3D`'s fresh multi-state start
drew random labels, and the run is not the known-label control it was meant
to be. Its first iteration assigned 2,700 particles against the truth maps,
0.847 of them to their true state; the second, against the maps of those
labels, 0.957 of all 3,000. Its log and statistics file show the readout:
iteration 2's `>>> PROB_TAB_STATE_DRAW: n_eff over states 1.00 at beta 0.30 /
1.00 at beta 1.00, 99.9% at the state mode, mean max state posterior 0.999
(1000 rows)` and per state 0.027 / 0.002 / 0.189 unrefined, written to
`simple_stats.txt` as `STATE_NEFF_DRAW=1.0044`, `STATE_NEFF_PHYS=1.0008`,
`STATE_MODE_FRAC=0.9985`, `STATE_MASS_UNREFINED01..03=0.0273, 0.0020, 0.1892`.

**`states` and `cavgs`** once each on the joint draw (eight threads):
`states` 3 min 32 s, every floor passed (adjusted Rand index 1.000, accuracy
1.000, margin 0.450; it runs the neighbourhood table, unchanged here); `cavgs` 7
min 08 s, every true state resolved, accuracy 0.935, adjusted Rand index 0.881,
every check and ceiling passed (`R/scratch/keep/runs/p1end_{states,cavgs}/`).

**Review.** An independent review of the diff found no defect; its two
consistency points were taken (the state readout keeps a row as the n_eff
readout does; the statistics file gets the readout only when the log prints
it) and its design note is Phase 2's (the joint draw compares the leading
state's polished entries with the other states' unpolished ones until the
per-state screen and polish exist).

**Exit items.** The fast gate (15 of 15) and the unit suites pass. The identity
holds at one thread and at eight in shared memory and at two parts, and on the
6vxx workflow at one thread. The `labels` reports and the readout, with their
reading, are above. `states` and `cavgs` pass.

### Phase 2 (2026-11-05, dell2): the fill's split, the neighbourhood table, the per-state convergence, the final pass

**Files changed**: `simple_eul_prob_tab.f90` (the per-state screen and polish),
`simple_eul_prob_tab_neigh.f90` (the pooled draw, the extras from the committed
state, its state readout), `simple_eul_prob_tab_utils.f90` (the readout's rows
with a state per entry), `simple_eul_prob_tab_tester.f90`,
`simple_convergence.f90` and its tester (the per-state stall rule),
`simple_commanders_solve3D.f90` (the final assignment),
`simple_refine3D_strategy.f90` (a defect the review found, below), the fill
harness `simple_polar_matching_bench.f90` and its driver
`simple_commanders_test_highlevel_matchbench.f90` (added to the file table for
this phase), this plan (the table, a correction of Phase 1's fixture run).

**What was built.**

- *The per-state screen and polish* (section 4.2). With more than one active
  state the refined-reference budget `projs_ns` is split over the states'
  blocks: `ceil(projs_ns / n_active)` per block (at most `n_active - 1`
  entries more in all). The window search polishes each block's best; the polar
  shift map keeps today's pool size, takes each block's share of it at the
  origin, tracks the screen's width cutoff per block and polishes each block's
  share of its own screen. With one state the old statements run unchanged.
- *The pooled neighbourhood draw* (section 4.3). The draw runs over the valid
  candidates of every state together; the extras gather the committed state's
  candidates only and find the committed candidate by its position in the
  candidate store, so the cone and the renormalisation are the committed
  state's. The neighbourhood table now gives the state readout too.
- *The per-state stall test* (section 4.4). The convergence object keeps every
  state's two last resolutions; the annealing tail starts when no populated
  state gained (a state is compared only when both its crossings, FSC 0.5 and
  0.143, are in this check and in the previous iteration's record; the
  weighted insertion's two-iteration margin applies per state); when no state
  can be compared, today's best-state test decides. With one state the old
  test runs unchanged; the resolution exit stays on the best state.
- *The final assignment* (section 4.5). Before the final reconstruction of a
  multi-state `solve3D`, one `refine3D` iteration over every active particle
  from the last stage's command line (its band, crop and maps) with
  `refine=prob prob_beta=1` (no tail), `update_frac=1`, `balance=none`,
  `fillin=no`, `trail_rec=no`, `volrec=no`, `topk=1`; the coverage check after
  it must find no active particle without an update. It replaces the greedy
  missing-update pass, deleted (it served only `nstates > 1`). The check is
  satisfied by the sample itself (every active particle is sampled at
  `update_frac=1` and every sampled particle is drawn), as the old check was.

**Defects found and fixed.** An independent review of the diff found two:

- the distributed master seeded and checked the multi-state labels at every
  `refine3D` start, unlike the shared-memory strategy, which does it only at a
  fresh start (the first iteration). A state drained in one stage would stop
  the next distributed stage, and now the final assignment, on "insufficient
  state population". The distributed path now has the same fresh-start gate
  (multi-state runs only);
- in the pooled neighbourhood row a particle whose valid candidates all carry
  a non-finite distance gave the draw no reference (index 0); the old code kept
  such a particle. It is kept again (the old condition exactly). The
  single-state `simulated_workflow_6vxx` met this case in its hill-climbing
  neighbourhood stage on the binary before the fix
  (`R/scratch/keep/identity/phase2/6vxx/new_prefix_failed/`): it would have broken
  single-state runs.

**Tests.** The fast gate (15 of 15) and the unit suites pass:
`unit_pftc_align2D3D` "3D probability table" 453 of 453 (twice) and
`unit_cart_align3D` "pose statistics" 142 of 142 (logs
`R/scratch/runs/unit/ut_prob_p2c_*.log`, `ut_conv_p2b.log`). New: the
neighbourhood draw leaves the closest candidate's state for some particles at
a finite temperature (the hill-climbing store with unequal candidate counts per
state); in the state and geometric modes every extra is a refined candidate of
the committed state and never the committed one, some particles commit a state
other than their closest candidate's and get extras, and the extras' parallel
selection gives the one-thread result bit for bit on the 16-thread team, twice;
a particle without a finite candidate distance is kept; the stall rule holds
the tail while one state still gains though the best state stalls, skips a
state without a resolution, and falls back to the best state when none can be
compared.

**The fill harness** (`simple_test_exec test=polar_matching_bench suite=fills
nthr=16`, twice; `R/scratch/keep/runs/harness_p2/`): every route (dense,
neighbourhood, 2D dense and sparse), one and two states, gives the one-thread
reference's rows, selections and extras bit for bit at 16 threads, in batches,
in three parts and with the earlier parts' outputs perturbed (48 of 48
configurations each time). Its new check of the dense polish: with two states
no block holds more than its share (6 of a budget of 11) and no particle more
than the budget. Deviation from the plan's wording ("every active block
carries a polished entry per particle"): the polar shift map records a polish
only where it beats the stored shift's distance, a rule kept because changing
it would change single-state results; the split guarantees each block its
share of screen and polish attempts, and 1 of the 64 particles ended without a
recorded polish in one state.

**The fill's cost**, on a replica of the `labels` route's setup that keeps the
bench files (the route deletes its directory when it passes;
`R/scratch/keep/scripts/fillbench_fixture.sh`: the fixture's particles, one
consensus map as the three starting maps and `refine3D`'s own random labels,
four iterations at eight threads, alone): Phase 1 70.4 and 74.8 s of table
fill, this phase 72.3 and 75.3 s (+1.7 % on the means, within the tenth); 57
polished references per particle against 55 (three blocks of 19), shift
searches +3.6 %, coarse maps +0.9 % (`R/scratch/keep/runs/fillbench_p{1,2}_r{1,2}/`).

**The `labels` route** on this phase's code (before the two late fixes, which
touch neither; eight threads, alone): `labels` 3 min 36 s, `labels_identical`
3 min 53 s, `hard_labels` 2 min 34 s (`R/scratch/keep/runs/p2_*`). `labels` and
`hard_labels` give Phase 1's accuracy and state readout per iteration to four
digits (`labels` 0.34 … 0.83 \| 0.93 … 0.999; `hard_labels` 0.69, unresolved);
the poses differ slightly (median rotation error 1.794° against 1.790°). The
unrefined share of each state's tempered mass, up to 0.64 / 0.40 / 0.25 in
Phase 1's annealing tail, now stays at or below 0.061: the split removes the
starvation of the blocks the draw does not take, but with the state posterior
as sharp as it is here those blocks never decide a draw, so the labels do not
move. `labels_identical` ran its full ten iterations this time and drained
state 1 to 5 particles (Phase 1: 114 at its eighth iteration).

**The single-state identity** on the final binary (base build against this
tree, same seed and layout): the small single-state `refine3D refine=prob` at
one thread and at eight in shared memory and at two parts of one thread gives
byte-identical maps, FSC, sigma2 state and statistics file, identical particle
records and identical log readouts, the differing files holding only the run
directory's path or process numbers as between two runs of one build
(`R/scratch/keep/identity/phase2/small_refine3D_final/`); and
`simulated_workflow_6vxx` at one thread, `test=identity`, 181 declared outputs
the same, 0 differing (`R/scratch/keep/identity/phase2/6vxx/`).

**The routes** on the final binary (eight threads; `R/scratch/keep/runs/p2end_*`):

- `states` 3 min 07 s: every floor passed (adjusted Rand index 1.000, accuracy
  1.000, margin 0.450, rotation error 1.25°); its neighbourhood stage now
  prints the state readout (n_eff 1.00, 100 % at the state mode at the end);
- `cavgs` 6 min 48 s: every true state resolved, accuracy 1.000, every check and
  ceiling passed;
- `consensus` 10 min 39 s: passed (accuracy 0.872, adjusted Rand index 0.650,
  the partial merge of about 380-520 of the middle state's particles seen in
  Phase 0's passing runs);
- `shc` 20 min 18 s, run beyond the plan's list because it is the one route
  that reaches the final assignment (particle `solve3D` with three states):
  every true state resolved, accuracy 1.000; the annealing tail started on the
  stall of every state in stage 5; the final assignment ran as iteration 67 over
  all 3,000 active particles at `prob_beta` 1 (2,500 directions x 3 states,
  state n_eff 1.00, the labels kept: state overlap 0.999-1.000), and the
  coverage check after it found none missing. (`shc` took 16 minutes inside
  Phase 0's entry on the tree before the batch alignment refactoring; the two
  are not paired.)

**Exit items.** All pass: the fast gate, the unit suites, the fill harness at
16 threads, the three `labels` routes, `states`, `cavgs`, `consensus` (and
`shc`). The identity holds at one and eight threads, at two parts and on the
6vxx workflow. The fill's cost is within a tenth of Phase 1's (+1.7 % on the
table fill). The reports are above.

### Phase 3 (2026-11-05, dell2): robustness on small data

**Files changed**: `simple_commanders_test_highlevel_msworkflows.f90` (the
`labels` route's floors), `doc/policies/test_environment_policy.md` (the gated
route, the profiles and the entry's time), `simple_test_registry.f90` (the
entry's time in its comment), this plan.

**The `labels` profile** (`simple_test_exec test=profile
target=multistate_workflows suite=labels nseeds=20 nthr=8 ncunits=16`, alone,
78 minutes; `R/scratch/keep/runs/p3_profile_labels/`): every true state
resolved in 13 of 20 runs (accuracy 0.998-1.000 in those; the shuffled-label
control in none). In the other 7 a random start lost a state: from the second
iteration one estimated state loses its particles to the other two, steadily
through the flat iterations at `prob_beta` 0.3 and sometimes into the
annealing tail, until two true states share one estimated state (accuracy
0.79-0.80, the smallest recall 0). In every run the state readout shows the
draw near-deterministic from the second iteration (about 1.0 effective
states).

**The drains diagnosed on `labels`** (scratch builds, deleted after;
`R/scratch/keep/runs/p3_diag_labels/`): at three of the draining seeds the code
before the draw (the closest-state rule) drains the same way; the joint draw at
`prob_beta` 0.1 (1.7-1.8 effective states at the first iteration) drains at all
three, and at 0.03 (2.6 effective states at first) at two of them. So the drain
is the reconstruction's feedback from a random start, not the draw, and no
tempering factor the rulings allow, or a separate one they exclude, removes it:
section 6's risk 4 ("a state drains when its map is garbage at any
temperature"), whose remedy, reseeding a drained state, is outside this work
(ruling 7). It is not a stop condition: the per-stage tail fallback treats a
draw too hot, the opposite of what the readout shows, and the joint draw is
not worse than the rule it replaced. The route is therefore gated as a
robustness profile.

**The `labels` floors** (by the floor rule, ruling 9, from that profile;
written with their evidence in the suite): every true state resolved in at
least 6 of 20 runs (13 - 3 x 2.13 = 6.6), adjusted Rand index at least 0.70
(the bound 0.739), accuracy at least 0.75 (0.792), largest population error at
most 0.30 (0.2015), median rotation error at most 2.0° (1.879) and median shift
error at most 0.10 px (0.0674); the smallest recall and the discrimination
margin reported (their worst runs, 0 and -2, admit no floor), the sign of the
maps' differences checked (above 0.50 wherever every map exists). On the
compositional states only; `hard_labels` reports with a success floor of 0,
`labels_identical` is the control. The gated route passes at a resolving and at
a draining seed (`R/scratch/keep/runs/p3_gated/`).

**The `shc` profile** (`simple_test_exec test=profile
target=multistate_workflows suite=shc nseeds=20 nthr=8 ncunits=16`, alone,
6.5 hours; `R/scratch/keep/runs/p3_profile_shc/`): every true state resolved in
15 of 20 runs, against 13 of 20 on 2026-10-19 (the code has changed since, so
the counts are not paired). That meets the floor of 6, and the shuffled-label
control resolved 0 of 20. Every execution passed and every check held in all 20
runs. In the five unresolved runs (seeds 1, 4, 8, 16 and 19) two true states
merged: accuracy 0.80, smallest recall 0, largest population error 0.20. The
resolved runs reached accuracy 0.965-1.000. Medians: accuracy 0.984, adjusted
Rand index 0.95, rotation error 1.51° (worst 1.90, ceiling 2.5), shift error
0.075 px (worst 0.088, ceiling 0.15), map correlations 0.975, 0.967 and 0.965.
The route's `prob` stages draw the state almost deterministically: 1.00-1.08
effective states, with 97.3-100 % of particles at the state mode over the 851
readout lines. `labels` behaves the same way.

**The routes once each** (eight threads, alone, on the final binary;
`R/scratch/keep/runs/p3_<route>/route.log`):

- `shc` (21 min 49 s): passed, every true state resolved. Accuracy 0.904,
  smallest recall 0.52 (288 of the 600 particles of the A-B dimer went to the
  state of the tetramer less chain D), rotation error 1.51°.
- `hard_shc` (12 min 47 s): passed; the route is reported, not gated. Its
  states were not resolved: all 3,000 particles went to one estimated state and
  the other stayed empty (accuracy 0.80, which equals the collapsed accuracy;
  adjusted Rand index 0); rotation error 1.09°. The code before the draw (the
  base build) collapses the same way at the same seed
  (`R/scratch/keep/runs/p3base_hard_shc/`), so the collapse predates this work.
- `shc_identical` (16 min 51 s): passed. Two estimated states, adjusted Rand
  index 0.000, accuracy 0.405 inside the null (median 0.407): the control
  behaves.

**The full entry** (`simple_test_exec test=multistate_workflows nthr=8`, all
five routes, alone; `R/scratch/keep/runs/p3_entry/route.log`): passed in 2,692 s.
With `labels` gated, it resolved every true state, with adjusted Rand index
0.997 (floor 0.70), accuracy 0.999 (0.75), population error 0.001 (ceiling
0.30), rotation error 1.79° (2.0) and shift error 0.058 px (0.10). `cavgs` and
`shc` resolved every state at accuracy 1.000. `states` reached 1.000 with
margin 0.450 (floor 0.40). `consensus` reached adjusted Rand index 0.650 (0.40)
and accuracy 0.872 (0.70), with every state resolved. The entry took 42 minutes
in Phase 0's run on an earlier tree; the registry comment and the policy now
give the new time.

**The nightly** (`simple_test_exec test=nightly ncunits=16`, alone, on the
final binary, 1 h 19 min; `R/scratch/keep/runs/p3_nightly/`, with the report
and every entry's log and metrics): every one of the 19 entries passed,
including `simulated_workflow_1jxy`, `simulated_workflow_6vxx` and
`multistate_workflows` (2,869 s beside the other entries). It gated `labels`
at the same values as the entry, since the seed is the same: adjusted Rand
index 0.997, accuracy 0.999, population error 0.001, rotation error 1.79°,
shift error 0.058 px, every true state resolved.

**Exit items.**

- The entry passes with `labels` gated: passed (2,692 s), and every gated row of
  `labels` is in the nightly's report with its floor.
- The nightly passes with `labels` gated: passed, 19 of 19 entries.
- The `shc` floor holds: 15 of 20 against the floor of 6 (13 of 20 on
  2026-10-19, not paired).
- The profiles are in Progress: the two profile paragraphs above.
- The single routes ran once each: `shc`, `hard_shc` and `shc_identical`.

`hard_shc` resolves nothing, but neither does the code before the draw at the
same seed; the route is reported only. The fast gate (15 of 15) passed after
the registry comment changed, and the build was installed before the nightly
started.

### Phase 4 (2026-11-05/06, dell2): validation on the subset

**Files changed**: this plan only (this section). The runs used the final
binary of Phase 3 (`repo/build`, installed 2026-11-05 19:27), `SIMPLE_SEED`
unset, one job at a time on the machine. The scripts are in
`R/scratch/keep/scripts/p4_*.sh` and `chain_p4_*.sh`; each run's labels, logs,
statistics files and final maps are in `R/scratch/keep/runs/p4_*`. The `refine3D` runs on S add `iterstats=yes` (a
record of every iteration's statistics) and `mkdir=no` to the plan's
commands; neither changes a result. Accuracy is
the fraction of particles in the found state matched to their true state, under
the best one-to-one matching of found to true states; the recall of a true
state is the fraction of its particles that its matched found state holds; the
adjusted Rand index (ARI) is 0 for labels unrelated to the truth and 1 for the
truth up to renaming. The true states are the smallest assembly (8C9C, 6,000
particles of S), the medium (8C96, 12,000) and the largest (8C8X, 2,000).

**The truth maps are in the frame of S's true poses**
(`R/scratch/keep/runs/p4_truthframe/`). Reconstructed from PTS, the subset
project with the true poses and labels, each state's map matches its own truth
map at Fourier shell correlation (FSC) 0.994-0.998 at 40 Å and 0.875-0.968 at
12 Å, crossing 0.5 only at Nyquist (6 Å). Against the other truths it reaches
0.29-0.54 at 40 Å. So the truth maps serve as references without docking.

**The known-reference control** (`refine3D refine=prob nstates=3 pgrp=c1
mskdiam=300 lp=10 nspace=2500 update_frac=0.9 balance=none maxits=2
objfun=euclid ml_reg=yes prob_beta_end=1 nthr=16 print_states=yes` on PTS with
the truth maps as `vol1-3`, at `prob_beta` 1 and 0.3;
`R/scratch/keep/runs/p4_kref_beta{1,0.3}/`): every particle assigned went to
its true state, accuracy 1.0000 at both factors. Each iteration assigned the
18,000 particles it sampled: the first against the truth maps, the second
against the maps of the first's labels. By the second, all 20,000 particles had
been assigned at least once. The state posterior is completely sharp at 10 Å:
1.00 effective states and 100 % of draws at the state mode even at 0.3, where
the posterior over the row spreads over 2.9-3.2 effective references. Each
run took 12-13 minutes, 4.4 GB: 8-9 minutes for the first iteration (18,000
particles x 7,500 references) and about 4 for the second. The ceiling at `prob_beta` 1 is 1.0000, so the verdict thresholds keep
their defaults. They were fixed in the journal before the first de novo run: a
run is **separated** when every true state is recovered (its matched found
state holds at least 0.75 of its particles) and the accuracy is at least 0.90;
**partial** when the accuracy is at least 0.75 or exactly two true states are
recovered; **failed** otherwise.

**Label recovery on true poses** (two runs from one start: PTS with a balanced
random labelling, drawn once by `random.Random(20261106)` and applied with
`selection`, accuracy 0.338; the starting maps by `reconstruct3D nstates=3`;
then the control's command with `maxits=8` and the starting maps as `vol1-3`;
`R/scratch/keep/runs/p4_labelrec_{1,2}/`). The annealing tail must finish by
the eighth iteration, so it started at the fourth: three flat iterations at
0.3, then 0.3, 0.41, 0.55, 0.74 and 1.0. Per iteration: a pair of values is
run 1, then run 2; a single value or a range covers both runs. The state n_eff
is exp(entropy) of the state marginal, at the factor and at 1. "Labels kept"
is the share of the iteration's sampled particles whose state did not change.

| iteration | factor | accuracy | recall smallest / medium / largest | found populations | state n_eff at the factor / at 1 | draws at the state mode | labels kept |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 0.30 | 0.350, 0.357 | 0.36/0.35/0.33, 0.37/0.35/0.33 | 6678/6703/6619, 6634/6782/6584 | 1.86/1.23, 1.86/1.24 | 76 % | 0.68 |
| 2 | 0.30 | 0.383, 0.414 | 0.41/0.38/0.30, 0.47/0.40/0.31 | 6584/6873/6543, 6474/7012/6514 | 1.81/1.18, 1.82/1.20 | 77 % | 0.73-0.74 |
| 3 | 0.30 | 0.478, 0.553 | 0.57/0.47/0.26, 0.72/0.51/0.31 | 6215/7392/6393, 6017/7678/6305 | 1.55/1.13, 1.47/1.12 | 85-86 % | 0.74-0.78 |
| 4 | 0.30 | 0.629, 0.667 | 0.80/0.62/0.20, 0.91/0.61/0.30 | 5347/8770/5883, 5514/8675/5811 | 1.30/1.08, 1.20/1.06 | 91-94 % | 0.76-0.81 |
| 5 | 0.41 | 0.691, 0.689 | 0.91/0.67/0.16, 0.96/0.63/0.27 | 4678/9606/5716, 5245/8986/5769 | 1.10/1.03, 1.07/1.02 | 97 % | 0.88-0.93 |
| 6 | 0.55 | 0.699, 0.683 | 0.94/0.68/0.09, 0.97/0.62/0.19 | 4329/10001/5670, 5061/9124/5815 | 1.04/1.02, 1.05/1.02 | 99 % | 0.94-0.96 |
| 7 | 0.74 | 0.691, 0.668 | 0.95/0.67/0.02, 0.97/0.61/0.08 | 4227/10062/5711, 4934/9227/5839 | 1.03/1.02 | 99 % | 0.97 |
| 8 | 1.00 | 0.685, 0.648 | 0.96/0.66/0.00, 0.98/0.59/0.01 | 4294/9946/5760, 5054/9086/5860 | 1.02/1.02 | 99 % | 0.97-0.98 |

Both runs **failed** by the thresholds (ARI 0.48 in both), and in the same way.
The smallest assembly separated first and almost whole (5,757 and 5,858 of its
6,000 particles in one found state). The medium assembly split over the two
other found states (4,069 / 7,931 and 4,926 / 7,074). The largest assembly went
almost whole (1,993 and 1,974 of 2,000) into the state holding the larger part
of the medium. The final maps say the same (direct FSC against the truths,
`fsc_truth/`): two maps of the medium assembly, one leaning towards the largest
(0.71-0.73 to its truth at 40 Å), and one of the smallest (0.96-0.97 at 40 Å,
0.43-0.45 at 12 Å). The unrefined share of each block's tempered mass stayed at
or below 0.047. Each run took 53 minutes (6-9 minutes per iteration) and 4.5 GB.

What the readout says: the draw is neither too hot nor locked from the start.
At the first iteration the state posterior is the flattest it gets (1.86
effective states of three). The accuracy climbs from the second iteration while
a fifth to a third of the labels move each iteration, and the posterior
sharpens to about one state as the maps separate. The partition forms at the
flat factor. By iterations 3-4, still at 0.3, 85-94 % of the draws are at the
state mode and most of the smallest assembly sits in one state. That
partition splits the dominant class (60 % of S) instead of isolating the rare
one (10 %). The largest assembly's recall never rises above the start's 0.33:
run 2 holds it at 0.30-0.31 to iteration 4, run 1 loses it from the first
iteration, and both end at 0.00-0.01. The annealing tail fixes the partition
that has formed.

**The de novo runs** (`simple_exec prg=solve3D pgrp=c1 mskdiam=300 nstates=3
nparts=2 nthr=8 print_states=yes` on each run's own copy of P2S, the subset's
2D solution; `SIMPLE_SEED` unset; `R/scratch/keep/runs/p4_denovo_{1,2,3}/`,
with `analysis/` holding the log and the label files). On S, solve3D plans
five stages for three independent states, with iteration budgets of 20, 20,
17, 15 and 12:
1. `shc` at 24 Å, 20 iterations;
2. `shc` at 22.6 Å, converged after 4 iterations;
3. `prob` at 22.6 Å over 1,000 directions per state, 17 iterations;
4. `prob` at 20.2 Å over 1,000 directions, 15 iterations;
5. `prob` at 6 Å over 2,500 directions with the nonuniform filter, 12
   iterations. The stage's budget starts the annealing tail at iteration 64,
   so the last five of the twelve are the tail.

The final assignment runs after stage 5 as iteration 69: all 20,000 particles at
`prob_beta` 1. Each run sampled 18,000 particles per iteration, and the master
peaked at 1.7 GB. The readout columns:
- the row n_eff: the effective number of references of the posterior over the
  whole row (directions x states: 3,000 references in stages 3-4, 7,500 in
  stage 5) at the draw's factor;
- the state n_eff at that factor and at 1;
- the share of draws at the state mode;
- the largest share of a state block's tempered mass on unrefined entries;
- the share of the iteration's sampled particles whose state did not change.

*Run 1* (5 h 44 min; **failed**: final accuracy 0.589, ARI 0.224):

| stage | Å | iterations | time | accuracy at the end | ARI | recall smallest / medium / largest | populations |
| --- | --- | --- | --- | --- | --- | --- | --- |
| start | | | | 0.338 | 0.000 | 0.34/0.34/0.33 | 6667/6667/6666 |
| 1 `shc` | 24 | 20 | 275 s | 0.451 | 0.050 | 0.42/0.50/0.21 | 5736/5553/8711 |
| 2 `shc` | 22.6 | 4 | 120 s | 0.446 | 0.049 | 0.43/0.49/0.22 | 5857/5681/8462 |
| 3 `prob` | 22.6 | 17 | 4,135 s | 0.403 | 0.025 | 0.41/0.42/0.29 | 6667/6064/7269 |
| 4 `prob` | 20.2 | 15 | 3,133 s | 0.505 | 0.121 | 0.54/0.51/0.40 | 6908/5907/7185 |
| 5 `prob` + final | 6 | 12 + 1 | 12,265 s | 0.589 | 0.224 | 0.60/0.60/0.51 | 6535/5687/7778 |

| iteration (stage) | accuracy | row n_eff | state n_eff | at the state mode | unrefined | labels kept |
| --- | --- | --- | --- | --- | --- | --- |
| 25 (3, first) | 0.403 | 64 | 2.39 / 1.49 | 59 % | 0.27 | 0.53 |
| 33 (3, middle) | 0.390 | 116 | 2.33 / 1.39 | 60 % | 0.53 | 0.53 |
| 41 (3, last) | 0.403 | 105 | 2.21 / 1.33 | 65 % | 0.59 | 0.58 |
| 42 (4, first) | 0.414 | 106 | 2.25 / 1.37 | 63 % | 0.53 | 0.56 |
| 49 (4, middle) | 0.446 | 73 | 2.01 / 1.28 | 70 % | 0.50 | 0.64 |
| 56 (4, last) | 0.505 | 45 | 1.60 / 1.17 | 80 % | 0.48 | 0.74 |
| 57 (5, first) | 0.538 | 61 | 1.63 / 1.21 | 79 % | 0.14 | 0.70 |
| 63 (5, middle) | 0.579 | 13 | 1.23 / 1.06 | 91 % | 0.04 | 0.89 |
| 68 (5, tail end) | 0.589 | 1.25 | 1.02 / 1.02 | 99 % | 0.00 | 0.99 |
| 69 (final, all) | 0.589 | 1.18 | 1.01 / 1.01 | 99 % | 0.01 | 0.99 |

Final confusion (true rows; found states 1, 2, 3): smallest 2,385 / 3,603 /
12, medium 3,131 / 1,707 / 7,162, largest 1,019 / 377 / 604. One found state is
92 % medium assembly but holds only 60 % of it. A second is smallest-rich (63 %
pure, holding 60 % of the smallest assembly), the third a mixture, and the
largest assembly is scattered.
The final maps docked onto the truths in both hands (`dock/dock_eval.txt`):
the medium-rich state's map docks onto the medium truth in the mirror hand at
correlation 0.969, FSC 0.991 / 0.968 / 0.911 at 40 / 20 / 12 Å, crossing 0.5 at
6.2 Å. The two other maps dock onto nothing (correlation at most 0.29).

*Run 2* (4 h 52 min; **failed**, above chance but far from any separation:
final accuracy 0.423, ARI 0.035; a balanced labelling unrelated to the truth
scores about 0.34 on S, as the random starts did):

| stage | Å | iterations | time | accuracy at the end | ARI | recall smallest / medium / largest | populations |
| --- | --- | --- | --- | --- | --- | --- | --- |
| start | | | | 0.337 | 0.000 | 0.34/0.34/0.34 | 6666/6667/6667 |
| 1 `shc` | 24 | 20 | 279 s | 0.387 | 0.008 | 0.38/0.37/0.51 | 6822/7103/6075 |
| 2 `shc` | 22.6 | 4 | 120 s | 0.394 | 0.011 | 0.40/0.37/0.52 | 6865/7082/6053 |
| 3 `prob` | 22.6 | 17 | 4,174 s | 0.370 | 0.004 | 0.40/0.35/0.38 | 6652/6741/6607 |
| 4 `prob` | 20.2 | 15 | 3,375 s | 0.394 | 0.016 | 0.45/0.37/0.39 | 6751/6840/6409 |
| 5 `prob` + final | 6 | 12 + 1 | 9,072 s | 0.423 | 0.035 | 0.51/0.39/0.33 | 6997/6481/6522 |

| iteration (stage) | accuracy | row n_eff | state n_eff | at the state mode | unrefined | labels kept |
| --- | --- | --- | --- | --- | --- | --- |
| 25 (3, first) | 0.386 | 65 | 2.49 / 1.52 | 57 % | 0.24 | 0.49 |
| 33 (3, middle) | 0.376 | 122 | 2.42 / 1.43 | 59 % | 0.54 | 0.51 |
| 41 (3, last) | 0.370 | 111 | 2.30 / 1.41 | 62 % | 0.58 | 0.55 |
| 42 (4, first) | 0.379 | 103 | 2.34 / 1.39 | 60 % | 0.52 | 0.52 |
| 49 (4, middle) | 0.388 | 87 | 2.24 / 1.37 | 65 % | 0.51 | 0.59 |
| 56 (4, last) | 0.394 | 65 | 2.10 / 1.29 | 69 % | 0.50 | 0.62 |
| 57 (5, first) | 0.392 | 94 | 2.22 / 1.40 | 64 % | 0.10 | 0.58 |
| 63 (5, middle) | 0.413 | 23 | 1.72 / 1.17 | 80 % | 0.04 | 0.76 |
| 68 (5, tail end) | 0.422 | 1.45 | 1.06 / 1.06 | 98 % | 0.00 | 0.97 |
| 69 (final, all) | 0.423 | 1.38 | 1.05 / 1.05 | 98 % | 0.00 | 0.97 |

Final confusion: smallest 1,420 / 1,501 / 3,079, medium 4,714 / 4,318 / 2,968,
largest 863 / 662 / 475: three mixtures. The `shc` stages left one true state
slightly enriched (the largest assembly's recall 0.52), and the probabilistic
stages never organised the labels.
No final map docks onto a truth: the best is 0.49 (FSC 0.57 at 40 Å).

*Run 3* (5 h 24 min; **failed**: final accuracy 0.582, ARI 0.187, one true state
recovered):

| stage | Å | iterations | time | accuracy at the end | ARI | recall smallest / medium / largest | populations |
| --- | --- | --- | --- | --- | --- | --- | --- |
| start | | | | 0.336 | 0.000 | 0.34/0.34/0.35 | 6666/6667/6667 |
| 1 `shc` | 24 | 20 | 275 s | 0.426 | 0.049 | 0.41/0.48/0.16 | 8302/5328/6370 |
| 2 `shc` | 22.6 | 4 | 120 s | 0.423 | 0.048 | 0.42/0.47/0.16 | 8220/5390/6390 |
| 3 `prob` | 22.6 | 17 | 4,096 s | 0.394 | 0.017 | 0.43/0.39/0.30 | 7197/6152/6651 |
| 4 `prob` | 20.2 | 15 | 3,167 s | 0.480 | 0.086 | 0.60/0.44/0.38 | 7215/6208/6577 |
| 5 `prob` + final | 6 | 12 + 1 | 11,038 s | 0.582 | 0.187 | 0.76/0.50/0.54 | 6476/7201/6323 |

| iteration (stage) | accuracy | row n_eff | state n_eff | at the state mode | unrefined | labels kept |
| --- | --- | --- | --- | --- | --- | --- |
| 25 (3, first) | 0.396 | 63 | 2.34 / 1.48 | 60 % | 0.26 | 0.53 |
| 33 (3, middle) | 0.387 | 109 | 2.31 / 1.37 | 60 % | 0.52 | 0.53 |
| 41 (3, last) | 0.394 | 100 | 2.19 / 1.34 | 65 % | 0.58 | 0.58 |
| 42 (4, first) | 0.394 | 98 | 2.24 / 1.36 | 63 % | 0.55 | 0.55 |
| 49 (4, middle) | 0.415 | 72 | 2.07 / 1.29 | 68 % | 0.51 | 0.63 |
| 56 (4, last) | 0.480 | 41 | 1.72 / 1.20 | 77 % | 0.47 | 0.71 |
| 57 (5, first) | 0.481 | 56 | 1.79 / 1.25 | 76 % | 0.13 | 0.68 |
| 63 (5, middle) | 0.579 | 10 | 1.15 / 1.03 | 94 % | 0.04 | 0.92 |
| 68 (5, tail end) | 0.582 | 1.28 | 1.01 / 1.01 | 100 % | 0.01 | 0.99 |
| 69 (final, all) | 0.582 | 1.23 | 1.01 / 1.01 | 100 % | 0.00 | 0.99 |

Final confusion: smallest 26 / 4,547 / 1,427, medium 6,013 / 2,175 / 3,812,
largest 437 / 479 / 1,084. There is a 93 % pure medium state holding half the
medium assembly, and a smallest-rich state (63 % pure) holding 76 % of the
smallest assembly. The third state is a mixture that holds 54 % of the largest.
Docked: the medium state's map onto the medium truth in the mirror hand at
correlation 0.972 (FSC 0.991 / 0.966 / 0.910, 0.5 at 6.2 Å). The mixture's map
docks onto the medium truth too, in the other hand, at 0.892 (FSC 0.93 / 0.85
/ 0.71, 0.5 at 9.6 Å); the states are refined independently, so their hands
may differ. The smallest-rich state's map docks onto nothing (at most 0.21):
it collects the smallest assembly's particles without being its map.

In every run the final assignment (iteration 69) assigned all 20,000 active
particles at `prob_beta` 1 in 8-12.5 minutes (722, 501 and 749 s; the plan
estimated 15). It kept 97-99 % of the labels, with 98-100 % of draws at the
state mode, and the coverage check after it found none missing (20,000 updated
of 20,000 active). No state drained in any run. Each de novo run took 4 h 52 min
to 5 h 44 min alone, against the plan's estimate of about 9 hours: stages 3 and
4 took about 70 and 55 minutes, stage 5 with the final assignment 2.5-3.4 hours.

**The reading: what the joint draw does on the benchmark.**

- *The state signal is there.* With the truth maps and true poses every
  particle goes to its true state at 10 Å, at the search's factor as at 1.
- *From random labels on true poses, the draw organises the labels, but not
  into the true partition.* In both runs the smallest assembly separates first
  and almost whole, already at the flat factor. The medium assembly, 60 % of S,
  is split over the two other states. The largest, 10 %, never gains a state
  and ends inside the medium-rich one. The tail fixes a partition that formed
  before it.
- *De novo, it does not separate the states.* The three runs ended at accuracy
  0.589, 0.423 and 0.582 (ARI 0.22, 0.04, 0.19), all failed. In runs 1 and 3:
  - one state converged to a clean medium map (92-93 % pure; correlation 0.97
    to the truth, FSC 0.91 at 12 Å), holding half to 60 % of the medium
    assembly's particles;
  - a second state gathered 60 % and 76 % of the smallest assembly at 63 %
    purity, but its map docks onto no truth;
  - the third stayed a mixture.

  In run 2 nothing separated. The largest assembly was never isolated.
- *What the readout shows.* The `shc` stages freeze the labels at accuracy
  0.39-0.45, with about 2 % changing per iteration by their end. The first
  probabilistic iteration keeps only 0.49-0.53 of them. Through stage 3 (22.6 Å,
  1,000 directions per state):
  - the posterior over the row spreads over 60-210 references;
  - the state posterior at 0.3 spreads over 2.2-2.7 of the three states;
  - 0.23-0.59 of every state block's tempered mass lies on entries the fill did
    not refine;
  - accuracy holds at 0.37-0.40 while 42-57 % of the labels change each
    iteration.

  That stage meets section 7's test for "too hot", the state n_eff near the
  number of states while the accuracy stays near chance; so does run 2's
  stage 4. Runs 1 and 3 leave it during stage 4, ending it at 1.6-1.7
  effective states and accuracy 0.48-0.51, still climbing. Two things keep this
  from being read as a draw too hot and nothing else. By section 7's own
  condition, the marginal cannot be trusted while that much of the mass sits on
  unrefined entries. And even at `prob_beta` 1 the state posterior spreads over
  1.3-1.7 states while the accuracy is at chance. That fits maps that cannot yet
  tell the states apart at 20-22 Å as well as it fits a draw that is too hot.

  The labels organise as the maps sharpen, at the end of stage 4 and in stage 5
  (6 Å). By stage 5's middle (iteration 63) the row posterior has narrowed to
  10-25 references, the state posterior to 1.2-1.7 states, and the unrefined
  share to 0.04 or less. Stage 5 adds 0.03-0.10 of accuracy. Its annealing
  tail (iterations 64-68) and the final assignment change it by less than
  0.01.
- For orientation only (not paired, ruling 10): the first run's de novo run of
  the code before this work, which gave each particle its closest state,
  failed on the full benchmark (section 9). Its failure looked like runs 1 and
  3: one medium-locked state and mixtures, without the smallest-rich state.

**What is left** (for the maintainer; none of it was tried here):

1. *The per-stage tail.* Section 6's fallback for a state draw too hot is the
   annealing tail at the end of every multi-state `prob` stage: the controller
   would set the existing key `prob_beta_end=1` on stages 3-5. Stage 3 meets
   the test that calls for it. These runs do not show that it would help:
   - in stages 3-4 the state posterior is broad at `prob_beta` 1 too;
   - the marginal rests on unrefined entries there;
   - the existing tail at the end of stage 5 moved accuracy by less than 0.01.

   It was not tried in this phase, for three reasons. Ruling 10 makes the subset
   validate the finished development, not drive it. The fast tests do not call
   for it: `labels` resolves its states in 13 of 20 seeds (the others lose a
   state, as the closest-state rule does), and its readout is never too hot.
   And a controller change alters the `shc` route, whose floor would need its
   6.5-hour profile and the nightly again before new de novo runs. A separate
   state factor would be a new parameter (ruling 4): a question for the
   maintainer.
2. *The polish budget.* In stages 3-4 a quarter to over half of each state
   block's tempered mass lies on stored-shift losses, with the per-particle
   refined-reference budget split over the three blocks. That share comes with
   the posterior's spread over 60-210 references, though within stage 3 it
   rises while the spread narrows. Section 4.2's candidate rule, one polish
   count for every block, would make each block's mass rest on refined
   entries, at up to 2.5-2.8 times the table time (section 4.2, from the first
   run's measurements): a cost decision.
3. *The rare state.* From random labels on true poses, the rare state does not
   get a state of its own: the partition that forms at the flat factor splits
   the dominant class instead. These runs do not say what would change that.
   They end in a partition the truth does not support, which is how section 6
   says the row's in-plane bias would show (risk 3). But they cannot tell that
   bias apart from the partition's own dynamics, so the in-plane marginal stays
   the open extension, the maintainer's decision.
4. *Reseeding* an emptied state (risk 4) is not what fails on S: no state
   drained.

**Exit items.**

- The known-reference control is made and judged: accuracy 1.0000 at
  `prob_beta` 1 and 0.3. The truth maps were first checked to be in the frame of
  the true poses. The verdict thresholds were fixed from the control, in the
  journal, before the first de novo run.
- The two label-recovery runs are made and judged: both failed, the same way.
- The three de novo runs are made and judged: all failed.
- Each run has its tables, its final maps against the truth and its signature,
  above. The tables give accuracy and recall after every stage or iteration,
  iterations per stage, populations, and the readout at the first, middle and
  last iteration of each probabilistic stage.
- The reading and what is left are above.

Evidence, all under `R/scratch/keep/runs/`:
- `p4_truthframe/`: the reconstruction log and the FSC logs against the
  truths;
- `p4_kref_beta{1,0.3}/` and `p4_labelrec_{1,2}/`: the refine3D log, the
  per-iteration labels and statistics, and the final maps (with `fsc_truth/`
  for the label-recovery runs);
- `p4_labels_start/`: the drawn labels and the start's logs;
- `p4_denovo_{1,2,3}/`: the solve3D log, labels, statistics snapshots, project,
  final maps and docking tables.

The journal is `R/state/journal_phase_4.md`.

### Phase 5 (2026-11-06, dell2): documents and report

**Files changed**:
- `doc/policies/3D/solve3D_policy.md`: section 3's tempering bullet, section 6's
  state assignment, state readout and validation, and section 9's final-labels
  contract (the final assignment in place of the greedy missing-update pass);
- `doc/policies/3D/refine3D_policy.md`: the dense assignment (the joint draw,
  the split refined-reference budget), the annealing tail's stall test over
  every populated state, the state readout, and the neighbourhood assignment
  (the pooled draw, the extras gathered from the committed state);
- `doc/policies/heterogeneity/refine3D_states_policy.md`: the pooled draw under
  both pose policies;
- `doc/policies/test_environment_policy.md` and
  `doc/implementation_notes/completed/single_state_speed_week_closeout.md` (added
  to the file table for this phase): their references to this plan follow its
  move to `completed/`. The test policy's `labels` paragraph also gives the
  hotter draws' drains as found (0.1 at all three seeds checked, 0.03 at two);
- this plan: its header (now built), this section, and the move to
  `completed/`. Two corrections: Phase 2's evidence path for the 6vxx failure
  now points where the log is kept, and Phase 4 names the run settings it added
  (`iterstats=yes`, `mkdir=no`);
- the report `doc/implementation_notes/completed/multistate_posterior_state_draw_report.md`,
  rewritten for the whole run. It replaces Phase 0's draft and proposes the skill
  text.

**Decisions.**
- The policies describe the rule as built, and give the rule before it in one
  sentence each. The validation results go into `solve3D_policy.md` in two
  sentences: the small-data profiles and the benchmark's failure with its
  readout. The detail stays here.
- The stream's multi-state `solve3D` stage runs `solve3D` and inherits the
  change unedited. No stream or GUI code reads the deleted pass or the new
  statistics keys, so neither needs an edit.
- No skill is edited (the run's rule); the report proposes the text.

**Checks**: no build (documents only); the fence check of the phase; `git diff
--check`; the links the edits added resolve (`../3D/refine3D_policy.md` from
`heterogeneity/`; `../completed/multistate_posterior_state_draw.md` from
`completed/`). An independent review checked the policy edits against the code
and the report against this record. It found no contradiction with the code. It
did find:
- ranges given for stages 3-4 that were stage 3's alone;
- an overstated "every block carries polished entries";
- the "too hot" reading stated without Phase 4's caveat (a broad posterior at
  factor 1 fits maps that cannot yet tell the states apart);
- the section headings still naming the retired `ref_assign` routines;
- several over-claims in the report.

All were corrected.

**Exit items.** The documents are updated (above). The report is written
beside the plan. The plan is moved to
`doc/implementation_notes/completed/multistate_posterior_state_draw.md`.
