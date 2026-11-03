# Multi-state solve3D: the state as part of the posterior draw

Proposal of 2026-11-03, revised the same day after Hans's review (the review's
points are worked into sections 2-5 and 8), and the plan of the autonomous run
that builds and measures it on dell2 (sections 10-14, arun run `ms_posterior`).
Nothing is built. Written after the single-state posterior work
(`pose_marginalization.md`, sections 16-21; the 2D class draw of
`solve2D_policy.md`) and the multi-state regression work of 2026-10-21/22
(`multistate_solve3D_regression.md` and its two reports,
`multistate_stage3_diagnosis_report.md` in `completed/`). The frame is Hans's,
2026-11-03: hard state labels in solve3D (no weights; weights are the flex
machinery's), `topk=3` kept with the pose weighting per state, and the states
sampled probabilistically throughout the search from the same posterior the
single-state tables draw from.

## 1. What is proposed, in one paragraph

A particle's state is one draw from its tempered posterior over state and
direction together, the same draw, primitive, temperature and annealing the
single-state tables and the 2D class draw use: one `draw_posterior_row` over
the particle's whole row (every active state's block) instead of the argmin
over states followed by a draw within the argmin's block. The label stays hard:
one state per particle per iteration, one half-map, the weighted insertion's
extras from that state's block alone with weights summing to 1 within it. The
dense fill screens and polishes shifts per state, so the state draw compares
losses of one kind with losses of the same kind. The neighbourhood table draws
the same way over its pooled candidates and takes its extras from the drawn
state. Stages 1-2 (`shc`) stay as they are. One tempering factor, `prob_beta`,
flat at 0.3 and annealed to 1 in the last stage's tail, governs states and
directions alike, and a full-population assignment at `prob_beta=1` before the
final reconstruction gives every active particle a label drawn at the
physical factor against the final maps. A state-draw readout is added so that
these choices can be checked on Afan's three-state benchmark. No new user
input.

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
  state's polished losses. This is cause 4 of the regression plan ("the
  leading state is favoured; a state that falls behind early cannot come
  back"), left open by every run so far. Changing only the polish selection
  would keep the starvation, because the screen already chose the pool.
- **The row is an approximation to the pose evidence.** Each entry of the
  dense row is the loss at *one* in-plane angle drawn from the posterior over
  the grid rotations (`gen_posterior_val`), not the loss marginalised over
  the angle. A marginal over such entries is a posterior over the
  constructed row, not the physical evidence: sampling an angle and weighing
  its likelihood again favours a sharp angular distribution over a flat one
  of equal integrated evidence (two states with angular likelihoods
  (0.99, 0.01) and (0.5, 0.5) have equal evidence, and this procedure picks
  the first with probability about 0.66, by direct enumeration; Hans's
  review). The in-plane marginal exists as a primitive
  (`gen_posterior_evid`, used by the 2D class draw under the shift map) and
  is not in the 3D row today.
- **Neighbourhood table (`refine=prob_neigh`, stages 6-8 and
  `refine3D_states`).** `ref_assign_neigh`: `gather_state_row` chooses the
  state of the closest valid candidate over all states and gathers that
  state's candidates for the draw; `select_extras` calls `gather_state_row`
  *again* and therefore chooses the argmin's state on its own, and
  `which_pos` is the committed candidate's index in that gathering. Drawing
  the state elsewhere without changing the extras' gathering would take the
  extras from the wrong state, and gathering every state would admit extras
  from other states.
- **2D.** `ref_assign` and `ref_assign_sparse` in `simple_eul_prob_tab2D.f90`
  draw the class from the posterior over all active classes: a class in 2D is
  what a state is in 3D, and this is the template.
- **Stages 1-2.** `strategy3D_shc` walks a random order over
  `nspace x nstates` references and takes the first improvement: states and
  directions on an equal footing, the stochastic start the rulings keep
  (ruling 5, 2026-10-21).
- **The final labels today.** solve3D's stages run fractional updates
  (`nsample` 10,000 of the active particles per iteration, balanced
  stochastic sampling from stage 4). Only the particles sampled in an
  iteration receive that iteration's draw; the others keep labels drawn at
  earlier iterations, against earlier maps, at the factor of that iteration.
  Before the final reconstruction, `ensure_multistate_particle_assignments`
  (`simple_commanders_solve3D.f90`) only checks that every active particle
  was updated at least once and runs a greedy missing-update pass for the
  rest (`refine=greedy`, `volrec=no`, `maxits=1`). So the annealing tail puts
  no more than its iterations' samples at the physical factor.
- **Convergence state.** The convergence object keeps the resolution history
  of the last two checks as scalars of the best populated state
  (`res05_hist(2)`, `res143_hist(2)`, `res_hist_iter(2)`), and
  `prev_iter_res` reads the best state's value back from `STATS_FILE`.
  `check_conv3D` makes a fresh statistics record each iteration and writes
  it (`STATE_OVERLAP`, the per-state resolutions, ...): a value computed at
  the assignment reaches the file only through the convergence object.
- **The regression evidence.** On Cyril's simulated two-state set no build
  separates more than half its runs (anchor 2 of 4; the stage-skipping fix
  1 of 5; `prob` in stage 5 1 of 4 plus 1 partial). The failures share two
  signatures: the minority state drains during stage 4 or 5, or the labels
  freeze after one state locks on in stage 3 (H5_3: 18 label changes over 14
  iterations). Where a run separates, it does so in stage 4 at 6.0 A, once the
  band reaches where the two maps differ, and stage 5 completes it. The
  regression plan's Phase 4, a temperature on the state assignment in stages
  3-4 (ruling 6, a softmax over the best score per state), was planned and
  never built; no state-temperature code exists in the tree.

## 3. The model

For a sampled particle `i` with noise temperature `T_i = pi wsqsum / (2 pftsz)`
and the row's values `L_{s,r}` (state `s`, direction `r`, each entry at its
drawn in-plane angle and its refined or stored shift, as the fill leaves them),
the tempered posterior over the pair is

    p(s, r | row_i)  proportional to  pi_s exp(-beta T_i (L_{s,r} - L_min)),

with a uniform state prior `pi_s` (no prior) and `beta = prob_beta`. Its state
marginal is `p(s | row_i) = sum_r p(s, r | row_i)`. One draw from the full row
is a draw of `s` from this marginal followed by a draw of `r` within `s`: the
two-step rule is the one-step rule, so the implementation is the single-state
draw over a longer row, and `nstates=1` is the same code on the same random
stream.

This is a posterior over the constructed row. It is exact for the row, and the
row approximates the pose evidence in two ways the readouts must watch: each
entry's in-plane angle is a draw, not a marginal (section 2), and the entries
of a block are a mixture of polished and stored-shift losses that depends on
the refinement budget (section 4.2). The claims below about "the physical
posterior" at `beta = 1` are claims about this row's posterior at the noise
temperature, and the readout of section 7 is how far that is from the truth
on data with known labels.

**Why the marginal and not the best-per-state softmax of ruling 6.** The two
agree at `beta -> infinity` (both give the argmin's state) and differ at the
search's factor by an entropy term: a state whose map is blurrier has a
flatter block, so it holds more mass for the same best loss. Under a uniform
prior over directions that is the row's evidence, and it leans against the
drain: the sharper state no longer wins every near-tie. At the n_eff the
single-state readouts show (1-5 candidates), the difference is small; the
marginal is the form with no rule of its own. The in-plane bias of section 2
acts in the opposite direction (it favours the sharper angular distribution)
and is bounded by the same n_eff; whether the two matter on real rows is what
the label-only experiment of Phase 1 measures, and the in-plane marginal is
the extension of section 6 if they do.

**Why one tempering factor.** The tempering compensates an overconfident
posterior: the temperature is the noise's, while the map at a stage is
band-limited, regularised and, early on, wrong, so the posterior at `beta = 1`
overstates what the data decide. That applies to the state entries as it does
to the direction entries, and the argmin is the `beta = infinity` limit of the
state rule, the limit that locks on and drains. The tail anneals both to
`beta = 1` together in the last stage (`prob_beta_end=1`, set by the
controller), and the full-population pass of section 4.5 then gives every
active particle a label drawn at `beta = 1` against the final maps; without
that pass the tail reaches only its iterations' samples (section 2).

**Symmetry breaking.** The balanced random start gives `nstates` nearly
identical mixture maps, a flat state posterior and labels re-drawn almost at
random each iteration. The separation grows from the bias the maps' small
difference puts into the draws, which the next reconstruction amplifies:
stochastic EM. Slower than the argmin's amplification, with no lock-on. On
Cyril's set nothing separates under either rule before the band reaches 6 A,
so what the draw changes is in stages 4-5: a lagging state keeps a fair share
of the particles that fit both states nearly alike, and the labels keep moving
while the posterior is flat. The risk runs the other way, that at 0.3 the
state posterior is too flat to break the symmetry within the stage budgets;
section 6 lists the fallbacks and section 7 the readout that decides.

## 4. The changes

### 4.1 The dense assignment (`ref_assign`)

- Replace the argmin-over-states block and the within-block draw with one
  `draw_posterior_row` over the row's `nrefs` entries (the active states'
  blocks); the state is `ref_state(which)`, the direction `ref_proj(which)`.
  A row without a valid distance keeps reference 1, as now.
- The n_eff readout is offered the whole row. The state readout of section 7
  is computed here.
- The extras: unchanged in rule. `select_extras` reads the committed state
  from the assignment, so they come from the drawn state's block with that
  state's cone, and the weights renormalise within it. The
  refined-committed-entry rule stands.
- `nstates=1`: the same call on the same stream, bit for bit.

### 4.2 Per-state screen and polish in the fill (the fairness contract)

The contract: every active state's block is built the same way, so the state
marginal compares like with like. Both shift paths select per state:

- the window path: `refine_best_refs` takes the `n_s` best of state `s`'s
  block, `n_s` the state's own `calc_num2sample` (already computed per state;
  today only their maximum is used);
- the map path: the screen's pool is drawn per block (`minnloc` within the
  block, `screen_pool_size` on the block's distances, at least `2 n_s`), and
  the `n_s` best of each block's screen are polished.

The per-particle refined capacity becomes the sum of the `n_s`. The part-file
format is unchanged (refined entries are a list with their reference).

What this does and does not establish. Every block then carries its own
polished mode and its own screened set, so the drained-state bias of section 2
is gone. It does not make the blocks' marginal masses comparable in every
respect: `n_s` differs between states (it follows each state's mean pose
change), so the blocks hold different mixtures of polished and stored-shift
losses. Two readouts bound that: the share of each block's tempered mass on
unrefined entries, per state (section 7), and, on the label-only fixture of
Phase 1, the sensitivity of the state marginals and of the accuracy to the
refinement budget (`prob_athres` doubled against the default). If the mass on
unrefined entries is not small, the candidate rule is one `n` for every block
(the maximum of the `n_s`, which is what the fill allocates today).

Cost: the refinement is 31-41 % of the table's thread time in stages 3-5 with
two states (Phase 7 profile), so this adds up to that much again at two
states and more at three. `n_s` shrinks as a state settles, which bounds it.
The bench record (`prob_fill_profile`) measures it. If a prob stage's table
time grows by more than half at three states, the cheaper variant is built
and compared: the global top-`projs_ns` with each block's mode guaranteed a
polish, which keeps the cost and loses only the screened set's fairness.

### 4.3 The neighbourhood assignment and its extras (`ref_assign_neigh`)

Two gatherings, kept apart:

- for the draw, every valid candidate of every state (the pooled row); the
  draw runs once over it and the state is the drawn candidate's;
- for the extras, the candidates of the *committed* state only:
  `gather_state_row` takes the state as an argument instead of choosing it,
  and the committed candidate's index is remapped from its position in the
  pooled row to its position in the committed state's row (`which_pos`
  carries the store position, not the gathered index), so that
  `top_posterior_extras` excludes the committed entry and the cone is the
  committed state's.

`seed_fallback_if_empty` is unchanged. In `prob_neigh_mode=state` the pooled
neighbourhood is the same for every state (`build_pooled_neighborhood` copies
the pooled subspaces to every state) and in `geom` the geometric
neighbourhood is one subspace for every state, so the candidate sets are
comparable across states. In `prob_neigh_mode=shc|snhc` the candidates are
the hill climb's proposals in a random order over all states with an early
exit, so the counts per state are unequal and random; the draw over the
pooled store is then an acceptance step over the proposals and the incumbent
(as the fill's own comment has it), and the state is part of that acceptance.
This changes the shared routine for those modes too; no workflow uses them
with `nstates > 1` (multi-state solve3D's stages 1-2 run the matcher's direct
`shc`, not the table; single-state solve3D's stages 1-2 have one state), so
the only affected configuration is a hand-run `refine3D refine=prob_neigh
prob_neigh_mode=shc|snhc nstates>1`, and the test of section 5 covers it.
The rule reaches `refine3D_states` (it runs `prob_neigh` in state or geom
mode from `flex_pca` labels): one rule for the tables.

### 4.4 Convergence and the readout's route

- The pose-overlap bypass for `prob`/`prob_neigh` stands.
- The annealing tail's stall test reads every populated state, not the best
  one: the history becomes per state (`res05_hist(2, nstates)`,
  `res143_hist(2, nstates)`, allocated at the first check; a state is its
  index, which is stable within a run), `prev_iter_res` reads the per-state
  keys `RESOLUTION_FSC05_STATEnn` / `RESOLUTION_STATEnn` of `STATS_FILE`
  per state, and the stall holds when no state gained: a state is "not
  gaining" when both its crossings are present in both records and neither
  improved; a state missing either value in either record is skipped (an
  emptied state, a state without an FSC yet); when every state is skipped the
  rule falls back to the scalar best-state test as it is today. The
  two-iteration margin of `topk > 1` applies per state the same way.
- The state readout travels from the assignment to the statistics file
  through the convergence object: `prob_align_exec` / `prob_align_neigh_exec`
  return it (a small type beside the fill profile), the refine3D strategy
  stores it on `conv` after the pre-alignment in both execution modes, and
  `check_conv3D` writes `STATE_NEFF_DRAW`, `STATE_NEFF_PHYS`,
  `STATE_MODE_FRAC` and the per-state `STATE_MASS_UNREFINEDnn` beside
  `STATE_OVERLAP`. Writing them at the assignment alone would lose them.

### 4.5 The final labels: a full-population assignment at `beta = 1`

Proposed (decision 7 of section 8): before the final reconstruction of a
multi-state solve3D, one `refine3D` iteration over every active particle with
`refine=prob`, `prob_beta=1` (no tail), `update_frac=1.0`, `balance=none`,
`fillin=no`, `trail_rec=no`, `volrec=no`, `topk=1`, `maxits=1`, at the last
stage's band and crop against the last stage's maps, exactly as the greedy
missing-update pass is built today (`run_multistate_missing_update`), which it
subsumes: every active particle is updated, so nothing is missing. Its draws
are at the physical factor against the final maps; the final reconstruction
then uses labels of one provenance. Cost: one dense table over all active
particles (eight times a stage's table at `nsample` 10,000 of 80,000), measured
in Phase 0 on the base tree with `prob_beta=1`. Without it the claim "final
labels at `beta = 1`" is false for every particle not sampled in the tail's
iterations (section 2), and the alternative contract is to say so in the
policy and ship the mixed-provenance labels.

### 4.6 What does not change

- The stage schedule: `shc` in 1-2, `prob` in 3-5, `prob_neigh` state mode
  from 6; `MAXITS`; the 2,500-direction cap; the balanced stochastic sampling
  from stage 4; NU at stage 5; no stage omission or resolution exit for
  multi-state (ruling of 2026-10-21).
- `prob_beta` and `prob_beta_end` as the controller sets them (flat 0.3, the
  last stage's tail to 1).
- Hard labels everywhere: one half-map per particle, `calc_3Drec`'s grouping
  by `(state, half)`, the trailing accounting, the final reconstruction from
  the labels of 4.5. `m_estimator=flex` stays refused on multi-state
  `refine3D`; the state weight set is `flex_pca`'s and nothing here reads or
  writes it.
- Empty states: an emptied state (below `MIN_POP`) leaves the table and does
  not return; ruling 7's balanced relabelling at a stage start is a separate
  item, not changed by this proposal. No sampling or relabelling is driven by
  the maps (ruling of 2026-10-04).
- No new user input. The A/B control is the base tree built from a `git
  archive` export the same day, on the same fixtures.

## 5. Tests

- `simple_eul_prob_tab_tester`, `check_posterior_assignment`: the assertion
  "the state is the argmin's state under the draw" goes. New: at a huge
  temperature a particle with a unique row minimum takes the argmin's state
  and reference (the old rule's limit); tied minima across states are drawn
  uniformly over the tied entries (the current deterministic tie rule is
  replaced: `draw_posterior_row` at a huge temperature puts equal weight on
  exact ties, so over seeds each tied state is drawn about equally); at a
  finite temperature on a synthetic two-state row whose blocks hold known
  masses, the state frequencies over seeds within tolerance of the marginal;
  `nstates=1` draws identical to the base tree's at the same seed (the stream
  contract). `test_insert_extras_multistate`: every extra in the drawn
  state's block, unchanged in rule, with a case where the drawn state is not
  the argmin's.
- `test_ref_assign_neigh_posterior`: the pooled draw; a particle whose drawn
  state differs from the closest candidate's state, with extras only from
  the drawn state and the committed candidate excluded (`state` and `geom`
  neighbourhoods), and the `shc`/`snhc` store with unequal candidate counts
  per state (no extras, the state the drawn candidate's).
- The fill: after a two-state fill every active block carries a screened set
  and a polished entry per particle; the per-state mass readout on synthetic
  rows against a direct computation.
- The convergence: the per-state stall rule on synthetic records (one state
  gaining keeps the tail from starting; a state without an FSC is skipped;
  all skipped falls back to the scalar rule); the readout keys in the file.
- Identity: `simulated_workflow_1jxy` at one thread with the test's fixed
  seed, base tree against the new tree, through `test=identity` (the fills
  use dynamic scheduling and per-thread streams, so only a one-thread run is
  bit-comparable; the policy's command, `doc/policies/test_environment_policy.md`
  section 1.5); the fast gate; the tables' unit entries.
- Behaviour: the label-only experiment and the de novo arms of sections
  12-13 against a freshly built base tree on the same fixtures (the historical
  four-run results on Cyril's set cannot isolate this change from the search
  changes since); the three-state profile `simple_test_exec test=profile
  target=multistate_workflows suite=shc nseeds=20` (13 of 20 the last count;
  `nseeds` passed to `test=multistate_workflows` runs no profile); the whole
  `multistate_workflows` entry (the class-average route `cavgs`, the `shc`
  route, and the `states` and `consensus` routes through `refine3D_states`
  and so the pooled neighbourhood draw); `test=nightly`. Timing per stage
  against the base tree's, the per-state refinement and the final pass being
  the costs this proposal adds.

## 6. Risks, fallbacks and extensions

1. **Too hot for states.** The state n_eff stays near `nstates` through
   stage 4 and the maps never diverge. Fallbacks, in order: the tail to 1 at
   the end of every multi-state `prob` stage, not only the last (the
   controller sets `prob_beta_end=1` on stages 3-5), which commits the labels
   at each stage boundary while keeping the exploration within the stage; a
   state factor of its own (an internal key beside `prob_beta`, not a UI
   input), swept on the benchmark; more iterations for multi-state stages
   4-5. None in the first build.
2. **Cost of the per-state screen and polish** beyond what the stage can
   afford: the cheaper variant of 4.2.
3. **The row's in-plane bias** (section 2) shows in the label-only
   experiment as a state preference the truth does not support: the
   extension is the in-plane marginal per reference (`gen_posterior_evid`'s
   effective loss in the row for `nstates > 1`), measured on the same
   fixture. Not in the first build: it changes what the row holds.
4. **A state still drains** when its map is garbage (every particle fits the
   other state better at any temperature): ruling 7's reseeding, outside this
   proposal; a population floor in the prior would be map-driven and is ruled
   out.
5. **Ambiguous particles' final labels** are draws from the row's posterior
   at `beta = 1` against the final maps (section 4.5); a particle that fits
   two states alike lands in one of them at random. That is the honest hard
   label, and the mode would be the argmin the tables no longer take.

## 7. The readout

Per table assignment, beside `PROB_TAB_ASSIGN` and the n_eff line:
`>>> PROB_TAB_STATE_DRAW: n_eff over states <draw> / <beta 1>, <pct> % at the
state mode, mean max state posterior <p>`, the state n_eff being
exp(entropy) of the marginal, over the same subsample of rows as the n_eff
readout; and per state `>>> PROB_TAB_STATE_MASS (state s): <frac> of the
block's tempered mass on unrefined entries`. The populations per state are
already printed by the convergence object (`STATE n JOINT DISTRIBUTION OVERLAP
... POPULATION`) and `STATE OVERLAP` (`mi_state`, the fraction keeping their
state) is the per-iteration flow. The values go to `STATS_FILE` through the
convergence object (section 4.4). Read beside the per-iteration accuracy on
the benchmark (`print_states=yes`), they say whether 0.3 is too hot (n_eff
near `nstates` through stage 4 while accuracy stays at chance), too cold
(n_eff near 1 from stage 3's first iterations, the lock-on signature) or
right, and whether the unrefined mass is small enough for the marginal to be
trusted.

## 8. Decisions to take before coding

1. The joint draw (the marginal over the constructed row) rather than ruling
   6's best-per-state softmax: section 3's argument, with the row's
   approximations named and measured.
2. One tempering factor for states and directions, with the state factor and
   the per-stage tail as fallbacks only.
3. The per-state screen and polish as the default, its cost measured and the
   cheap variant as the fallback; one `n` for every block if the unrefined
   mass is not small.
4. The same rule in the neighbourhood table, hence in `refine3D_states`, with
   the extras from the committed state and the `shc|snhc` modes covered.
5. The tail's stall test over every populated state, per-state history.
6. Stages 1-2 untouched, and no stage-schedule change in the first build.
7. The full-population assignment at `beta = 1` before the final
   reconstruction (section 4.5), recommended; the alternative is to document
   the mixed provenance of the final labels.

## 9. Order of work

The phases of section 13: the dense draw with its readout, tests and the
single-state identity, then the label-only experiment (Phase 1); the fairness
of the fill, the neighbourhood table, the per-state convergence and the final
pass (Phase 2); the de novo arms (Phase 3); the fallbacks if needed, the
profiles, the policy documents and the report (Phase 4).

# The run: `ms_posterior` on dell2

## 10. Rulings

Fixed by the maintainer before the launch (Hans, 2026-11-03, and the standing
rulings of 2026-10-04 to 2026-10-21 that hold):

1. State labels are hard: one state per particle per iteration, one half-map,
   no state weights anywhere in solve3D. The flex weight set is `flex_pca`'s
   and is neither read nor written.
2. `topk` stays 3; the weighted insertion's extras come from the committed
   state alone, with weights renormalised within it.
3. The state is part of the posterior draw: one draw over the row across
   states, at the one tempering factor `prob_beta` with the one annealing
   tail. No new user input; internal keys only where a fallback needs one.
4. `nstates=1` stays bit for bit identical at one thread; a change that
   alters it and cannot be made not to is a stop.
5. Stages 1-2 (`shc`) stay as they are. The stage schedule, budgets, the
   direction cap, NU at stage 5 and the speed changes stay.
6. Sampling and relabelling are never driven by the maps or poses; an
   emptied state's reseeding is outside this run.
7. State 0 is exclusion; a true state in state 0 is a failure. De novo
   routes are judged over many runs, never one. Bugs found are fixed in the
   run.
8. No floor, tolerance or timeout of an existing test changes; proposals only.
   Rulings taken during the run go to `spec/rulings.md` and are recorded in
   this plan where they belong.

## 11. What we know before the run

### The data

`/data2/multistate_bench` on dell2 (copied from the dell's
`~/multistate_bench`, Afan's set, 2026-11-03; read-only): the Ribosembly
three-state benchmark. 80,000 independently simulated particles, 128 x 128 at
3.0 A, in one shuffled stack `particles.mrcs` (5.2 GB); the true labels in
`ground_truth_labels.txt` (one 1-based label per line, in stack order, which
is the project row order): state 1 the smallest assembly (8C9C, 24,000, 30 %),
state 2 the medium (8C96, 48,000, 60 %), state 3 the largest (8C8X, 8,000,
10 %). `particles.star` holds the true poses, shifts and CTFs;
`bench.simple` is the imported project with the true poses and shifts in
`ptcl2D` and `ptcl3D` (all particles in state 1); `bench_zero.simple` the same
with the angles and shifts zeroed. Both refer to the stack at a beegfs path
that does not exist on dell2, so the run re-imports (section 12).
`reconstructed_volumes/` holds the three truth maps, reconstructed by SIMPLE
from the true poses and labels, in the frame of the true poses (reference FSC
at 12 A: 0.934, 0.985, 0.964). SNR 0.01 by global image variance; the
per-state SNR follows the assembly's signal. The README says the states are
recoverable by cryoSPARC's ab initio.

Checked on the dell, 2026-11-03: `print_project_field oritype=ptcl3D` on
either project prints exactly 80,000 record lines (one per particle, no other
lines), each with `dfx`, `dfy`, `angast` (microns, degrees), `e1-e3`, `x`,
`y`, `eo` and `state`, so the printout is an `oritab` for `import_particles`
(`ctf=yes smpd=3 kv=300 cs=2.7 fraca=0.07`, the units the program's defaults).
The truth maps' extents at 3 standard deviations above the mean are 211, 224
and 267 A (box 384 A), so `mskdiam=300`.

### What this set tests that Cyril's did not

Three states of different mass (compositional heterogeneity), one of them
small (10 %). Accuracy alone cannot see the small state: merging it into the
medium state scores 0.90. The verdict therefore requires every true state
recovered (section 12). The states differ at low resolution, so the
separation can form before the 6 A band Cyril's set needed; the draw's
symmetry breaking is then tested where it should be easiest, and the
small state's survival where the drain is most likely.

## 12. How runs are made and judged

**Machine.** dell2 (Ubuntu 26.04, 16 cores: 4 P, 8 E, 4 low-power E; 61 GB;
`/` and `/data2` about 935 GB each). One SIMPLE job at a time on the whole
machine when it is timed: distributed `nparts=2 nthr=8`, shared memory
`nthr=16`. Builds and analysis may run beside a job when the job is not
timed; the load is noted with every time. Timings are not comparable with the
dell's.

**Builds.** The run's tree in `repo/build` by `./compile_clean.sh` (Release,
tests on, about 2 min, fast gate 15 of 15). The base: `git archive <base HEAD>
| tar -x -C scratch/keep/src_base`, built in `scratch/keep/build_base` the
same way, the contemporaneous baseline of every comparison. Every change to
the tree is rebuilt incrementally (`make -C build -j16`, the fast gate,
`make install`). Builds no longer needed are deleted.

**Projects.** Made once in Phase 0 under `scratch/keep/projects/`:

- PT (true poses): `new_project`, then `import_particles stk=<the stack>
  oritab=<the printout of bench.simple's ptcl3D field> ctf=yes smpd=3
  kv=300 cs=2.7 fraca=0.07`; the labels then set from
  `ground_truth_labels.txt` by `selection oritype=ptcl3D infile=` (as the
  regression's truth project was), after checking that `import_particles`
  left every state at 1.
- PZ (zero poses): the same from `bench_zero.simple`'s printout, every
  particle in state 1.
- P2: `solve2D ncls=30 mskdiam=300 nthr=16` on a copy of PZ; every de novo
  run starts from its own copy of P2's directory.
- The project rows, the stack images and the label lines are the same index;
  Phase 0 checks it on the truth reconstruction.

**The de novo command** (the product under test; no `SIMPLE_SEED`, independent
repetitions):

```
simple_exec prg=solve3D pgrp=c1 mskdiam=300 nstates=3 nparts=2 nthr=8 print_states=yes projfile=<its copy of P2>
```

**The label-only fixture** (the state draw alone, Phase 1): PT with balanced
random labels and the maps reconstructed from them (the fresh multi-state
start of `refine3D`: all particles in state 1 and `nstates=3` draw uniform
labels; if the program does not reconstruct the starting maps itself, the
labels are drawn once with `refine3D`'s rule and `reconstruct3D nstates=3`
makes `vol1-3`), then `refine3D refine=prob nstates=3 pgrp=c1 mskdiam=300
lp=10 nspace=2500 update_frac=0.125 balance=none maxits=12 objfun=euclid
ml_reg=yes topk=3 nthr=16 print_states=yes`, shared memory, no `SIMPLE_SEED`,
three repetitions per build (`update_frac` 0.125 is the stages' `nsample` of
10,000 over 80,000 particles; `nsample` is a solve3D input that `refine3D`
does not consume; `prob_beta`, `print_states` and `update_frac` are registered
keys, which is what the command-line parser accepts). The poses are true, so the rows are sharp about
the truth and the state entries are the only open question: what the run
measures is the accuracy per iteration under the draw against the argmin, and
the readout beside it. The known-reference control is the same command on PT
with the true labels and the truth maps as `vol1-3` for two iterations, on
both builds, at `prob_beta=0.3` and at `prob_beta=1`: the accuracy each keeps
is the ceiling of the state assignment on this data at this band.

**Recorded per run:** the SIMPLE git-commit line (empty for an export: name
the build directory); wall time per stage and the machine's load; the labels
after every stage (`states_stage<N>.txt`) and iteration (`states_iter<NNN>.txt`,
`pind state updatecnt`); the final labels; the log lines `ENTERING REFINE3D
STAGE`, `STAGE <n> WITH`, `PROB_TAB_ASSIGN`, `POSTERIOR DRAW N_EFF`,
`PROB_TAB_STATE_DRAW`, `PROB_TAB_STATE_MASS`, `PROB_TAB_INSERT`, `STATE <s>
JOINT DISTRIBUTION OVERLAP ... POPULATION`, `STATE OVERLAP`, the per-state
resolution lines, `RESOLUTION NO LONGER GAINING`, `BUDGET END`, `ANNEALING
TAIL COMPLETE`, any empty-state or relabelling message; `STATS_FILE`; the
final state maps and their docking onto the truth maps (`dock_volpair`, both
hands); the `BENCH/` directory of a timed run (the table profile).

**Metrics:** accuracy under the best of the six matchings of found states to
true states (state 0 always wrong); the per-true-state recall under that
matching; ARI; the confusion matrix (true against found 0-3); the
populations against 24,000 / 48,000 / 8,000; the number of particles solve2D
left inactive (excluded from accuracy, reported); the label change rate per
iteration (from the iteration files); the state n_eff, mode fraction and
unrefined mass per iteration; each final map docked onto each truth map in
both hands (correlation and FSC against the truth); iterations per stage; the
signature (where the separation forms, whether the small state drains, whether
labels freeze).

**Verdict per run:** *separated* when every true state is recovered (its
matched found state holds at least 0.75 of its particles) and the accuracy is
at least 0.90; *partial* when the accuracy is at least 0.75 or exactly two
true states are recovered; *failed* otherwise. If the known-reference control
at the argmin (the base tree, `prob_beta=1`) keeps less than 0.95, both
accuracy thresholds drop by the shortfall. The thresholds are fixed in the
journal at the end of Phase 0 and never moved.

**Verdict per build:** equal counts of runs per build, interleaved, judged on
the count and the signature together; counts of six separate large differences
only, so the signatures carry the weight where the counts are close. A build is
*good* when most of its runs separate; the question of Phase 3 is whether the
new tree separates at least as often as the base and no longer shows the two
failure signatures.

**Tools.** SIMPLE programs do the SIMPLE work. The accuracy tool is one Python
script under `scratch/keep/tools/` (no numpy on dell2; `python3` is 3.14) that
reads a labels file and the truth and prints the metrics above; the regression
runs' `acc.py`, `run_metrics.py`, `pose_err.py`, `profile.py` and
`profile_table.py` are at `/data2/multistate_bench/tools_from_dell/` to adapt
(they assume two states and the dell's paths). It is tested in Phase 0 on the
true labels (1.0), on all-one-state (0.60), on the small state merged into the
medium (0.90 but state 3 not recovered) and on shuffled labels (about 0.46),
and used nowhere else than through its tests.

**Budget.** Phase 0 times one de novo run of the base tree and the journal
rescales the plan. The plan assumes 3-5 hours per de novo run (Cyril's set
took 140-165 minutes per run on the dell's 24 cores with two states and
30,000 particles; this set has 1.5 times the references, the same `nsample`
of 10,000, a smaller box, and dell2's 16 slower cores), about an hour per
label-only run, and two hours for the three-state profile. If the time runs
short, drop in this order: the budget-sensitivity run, the fourth and later
fallback runs, the de novo arms beyond four runs per build (never below four).

## 13. Phases

### Phase 0: builds, data, projects, tools, controls, the first base run

1. Build the tree and the base export (Release, fast gate). Record how SIMPLE
   builds on dell2 in the journal.
2. Check the data (sizes, the label counts 24,000 / 48,000 / 8,000, the
   printouts' line counts). Make PT, PZ and P2.
3. The accuracy tool, tested as section 12 says.
4. Truth control: `reconstruct3D nstates=3` on PT (true poses and labels);
   each map against its truth map by FSC in the same frame (no docking;
   expected above 0.9 to 12 A), and the three truth maps docked onto each
   other (the differences a separation must find). The same with the labels
   shuffled: the maps of a failed separation.
5. Known-reference control on the base tree at `prob_beta=1` (the argmin's
   ceiling) and at 0.3; fix the thresholds in the journal.
6. One de novo run of the base tree on P2 (arm B, run 1), timed, with its
   `BENCH/` kept: the per-stage times and the table's profile are the cost
   baseline. Then one full-population `refine3D refine=prob prob_beta=1
   update_frac=1.0 balance=none volrec=no maxits=1 topk=1` on its final
   project at its last stage's band (the cost of section 4.5), timed.
7. Rescale the plan in the journal.

Exit: both builds pass the fast gate; the projects exist with their row
identity checked; the tool passes its tests; the controls are recorded; the
thresholds are fixed; B_1 is judged and timed; the rescaled plan is written.

### Phase 1: the dense draw, its readout, the tests, the identity; the state draw alone

1. Section 4.1 in `ref_assign`; the state readout type in
   `simple_eul_prob_tab_utils.f90`; its route to `STATS_FILE` (section 4.4,
   the readout part only); the unit tests of section 5 for the dense table.
2. Identity: `simulated_workflow_1jxy` at one thread, base against new,
   `test=identity`, before anything else is changed.
3. The label-only experiment: three runs per build, interleaved; the
   known-reference control on the new tree at 0.3 and 1. Table in Progress:
   accuracy, recall per state, ARI and the readout per iteration, per run.
4. The budget-sensitivity run: one label-only run on the new tree with
   `prob_athres=20`, against the default's readout and accuracy.

Exit: the fast gate and the table's unit entry pass; the identity holds; the
label-only table is in Progress with its reading (does the draw classify at
least as well as the argmin on true poses, and what do n_eff and the unrefined
mass say).

### Phase 2: the fill's fairness, the neighbourhood table, the per-state convergence, the final pass

1. Section 4.2 in both fill paths, with the bench record read on the
   label-only fixture (table time against Phase 1's); the cheap variant if
   the rule of 4.2 fires.
2. Section 4.3 with the neighbourhood tests of section 5.
3. Section 4.4's per-state history and stall rule with its tester.
4. Section 4.5, the final pass, in the multi-state solve3D commander (it
   replaces the missing-update pass for `nstates > 1`); its cost from Phase 0.
5. Identity again (one thread, `simulated_workflow_1jxy`), the fast gate, the
   unit entries of the tables and the convergence.
6. The label-only experiment repeated on this tree (three runs), beside
   Phase 1's.

Exit: all tests pass, the identity holds, the costs are in the journal, the
repeated label-only table is in Progress.

### Phase 3: the de novo arms

Six runs of the new tree (arm N) and six of the base (arm B, B_1 of Phase 0
counting as the first), interleaved N, B, N, B, ... on copies of P2, judged as
section 12 says, with every readout of section 7 beside the accuracy per
iteration. The side-by-side table goes into Progress: per run the accuracy and
recall after every stage, iterations per stage, populations, the final maps
against the truth, the signature; per arm the counts.

Exit: the twelve runs are made and judged; the table and the reading (count,
signature, the readout's verdict on the factor) are in Progress.

### Phase 4: fallbacks if needed, profiles, documents, report

1. Only if Phase 3's reading says the draw is too hot or too cold (section
   7): the per-stage tail first (four runs), then a state factor (an internal
   key, four runs at one value chosen from the readout). Kept only if the
   count is not lower and the signatures improve; otherwise reported and left
   out of the tree.
2. `simple_test_exec test=profile target=multistate_workflows suite=shc
   nseeds=20 nthr=8 ncunits=16 keep_outputs=no`, alone on the machine,
   compared by count with the 13 of 20 of 2026-10-19 (the code has changed
   since, so runs cannot be paired). Fewer than 13 is reported with the
   discordant routes; 9 or fewer is a stop.
3. The whole `multistate_workflows` entry (`simple_test_exec
   test=multistate_workflows nthr=8`: the `cavgs`, `shc`, `states` and
   `consensus` routes) and `test=nightly ncunits=16` on the final tree.
4. Timings: the per-stage times of arm N against arm B; the table profile of
   one N run against B_1's.
5. The policy documents: `solve3D_policy.md` section 6, `refine3D_policy.md`'s
   assignment sections (dense and neighbourhood) and the final-labels
   contract, `refine3D_states_policy.md` where the pooled draw changes its
   description; skill text (`.claude/skills/simple-refine3d/references/bayesian-model.md`)
   proposed in the report, not edited.
6. The report beside the completed plan; the plan moved to `completed/`.

Exit: the profile and the entries pass as section 12 and this phase say; the
documents are updated; the report is written.

## 14. Common rules

- Nothing is committed; git is read-only; the base is a `git archive` export.
- Fix as you go: a defect found is fixed when found and recorded with its
  evidence; a change that moves results waits for the phase that owns it.
- Thresholds are fixed before the runs they judge, never loosened, no run
  repeated until it passes.
- New or changed OpenMP code gets a tester at 16 threads, repeated, compared
  bit for bit with one thread; no string temporaries inside parallel regions;
  work arrays sized by the team in use.
- No `goto`; THROW macros never continued with `//&`; no `.and.` short
  circuit on allocatable arrays; `use ..., only:`; no large type as a plain
  component (`doc/policies/compile_time_policy.md`).
- Unset `SIMPLE_SEED` for every benchmark run; the identity and the profile
  use the seeds their test environment sets.
- Bulky outputs of a run are deleted once its numbers are in the journal;
  `scratch/keep/runs/` keeps the labels, logs, statistics files and final
  maps of every run, `scratch/keep/projects/` the projects, `scratch/keep/`
  the base build.
- Write prose for people: Progress entries and the report read without other
  documents, with SIMPLE terms explained and acronyms spelled out on first
  use; name things by what they do.

### File table

| file | why | phase |
| --- | --- | --- |
| `doc/implementation_notes/planned/multistate_posterior_state_draw.md` | this plan and its progress | each |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` | the joint draw and the state readout (1); the per-state screen and polish, the unrefined-mass readout (2) | 1 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` | the state readout type (n_eff over states, mode fraction, unrefined mass) | 1 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_tester.f90` | the tables' tests | 1 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90` | the pooled draw; the extras gathered from the committed state | 2 |
| `src/main/strategies/parallelization/simple_prob_align_strategy.f90` | the readout handed back to the caller | 1 |
| `src/main/commanders/simple/simple_commanders_prob.f90` | the thin wrappers of the pre-alignment, if the routine's signature changes | 1 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | the readout stored on the convergence object in both execution modes | 1 |
| `src/main/simple_convergence.f90` | the readout keys (1); the per-state resolution history and the all-state stall test (2) | 1 2 |
| `src/main/simple_convergence_tester.f90` | tests of the per-state stall rule | 2 |
| `src/main/commanders/simple/simple_commanders_solve3D.f90` | the full-population assignment at prob_beta=1 before the final reconstruction | 2 |
| `src/main/solve/simple_solve3D_controller.f90` | the per-stage tail (the first fallback, only if Phase 3 calls for it) | 4 |
| `src/main/params/simple_parameters.f90` | a state tempering factor, internal key (the second fallback only) | 4 |
| `src/main/params/simple_parameters_parse.f90` | its registration (the second fallback only) | 4 |
| `src/main/params/simple_parameters_phases.f90` | its validation (the second fallback only) | 4 |
| `doc/policies/3D/solve3D_policy.md` | section 6 and the final-labels contract | 4 |
| `doc/policies/3D/refine3D_policy.md` | the assignment sections of the dense and neighbourhood tables | 4 |
| `doc/policies/heterogeneity/refine3D_states_policy.md` | the pooled neighbourhood draw | 4 |
| `doc/implementation_notes/completed/multistate_posterior_state_draw_report.md` | the report | 4 |
| `doc/implementation_notes/completed/multistate_posterior_state_draw.md` | the plan, moved when the run completes | 4 |

Other files a fix needs are added to this table, with the reason, by the phase
that finds the defect.

## Progress

(one section per phase, written by the run)
