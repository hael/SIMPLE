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
dell2). Nothing in sections 1-8 is built.

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
| `src/main/strategies/parallelization/simple_prob_align_strategy.f90` | the readout handed back to the caller | 1 |
| `src/main/commanders/simple/simple_commanders_prob.f90` | the thin wrappers of the pre-alignment, if the routine's signature changes | 1 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | the readout stored on the convergence object in both execution modes | 1 |
| `src/main/simple_convergence.f90` | the readout keys (1); the per-state resolution history and the all-state stall test (2) | 1 2 |
| `src/main/simple_convergence_tester.f90` | tests of the per-state stall rule | 2 |
| `src/main/commanders/simple/simple_commanders_solve3D.f90` | the full-population assignment at prob_beta=1 before the final reconstruction | 2 |
| `src/main/solve/simple_solve3D_controller.f90` | the per-stage tail, only if the fallback of section 6 is needed | 1 2 |
| `doc/policies/3D/solve3D_policy.md` | section 6 and the final-labels contract | 5 |
| `doc/policies/3D/refine3D_policy.md` | the assignment sections of the dense and neighbourhood tables | 5 |
| `doc/policies/heterogeneity/refine3D_states_policy.md` | the pooled neighbourhood draw | 5 |
| `doc/implementation_notes/completed/multistate_posterior_state_draw_report.md` | the report | 5 |
| `doc/implementation_notes/completed/multistate_posterior_state_draw.md` | the plan, moved when the run completes | 5 |

Other files a fix needs are added to this table, with the reason, by the phase
that finds the defect.

## Progress

(one section per phase, written by the run)
