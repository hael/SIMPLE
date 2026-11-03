# Multi-state solve3D: the state as part of the posterior draw

Proposal, 2026-11-03. Nothing here is built. Written after the single-state
posterior work (`pose_marginalization.md`, sections 16-21; the 2D class draw
of `solve2D_policy.md`) and the multi-state regression work of 2026-10-21/22
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
dense fill refines shifts per state, so the state draw compares refined losses
with refined losses. The neighbourhood table draws the same way over its
pooled candidates. Stages 1-2 (`shc`) stay as they are. One tempering factor,
`prob_beta`, flat at 0.3 and annealed to 1 in the last stage's tail, governs
states and directions alike; a state-draw readout is added so that this
choice can be checked on the two-state set. No new user input.

## 2. Where the code stands

- **Dense table (`refine=prob`, solve3D's stages 3-5).** `ref_assign`
  (`simple_eul_prob_tab.f90`): the state is the particle's closest reference
  over all states (argmin, ties to the lowest index); the direction one draw
  from `exp(-beta T_i (L - L_min))` within that state's block. The weighted
  insertion's extras (`select_extras`) come from the committed state's block,
  within a per-state cone (each state's FSC=0.5 on disk), renormalised over the
  committed entry and the extras: hard state membership, Hans's ruling of
  2026-10-29 (`refine3D_policy.md`, "The weighted insertion").
- **The fill.** The shift refinement (`refine_best_refs`, and the polish of the
  screen under the polar shift map, `process_particle_with_map`) takes the
  `projs_ns` best references *over all states* (`minnloc` over the whole row;
  `projs_ns` is the maximum of the per-state `calc_num2sample`). A state that
  trails gets few or no refined entries, so its block is scored at the stored
  shift and compared against the leading state's polished losses. This is
  cause 4 of the regression plan ("the leading state is favoured; a state that
  falls behind early cannot come back"), left open by every run so far.
- **Neighbourhood table (`refine=prob_neigh`, stages 6-8 and
  `refine3D_states`).** `ref_assign_neigh` (`gather_state_row`): the state of
  the closest valid candidate over all states, then the draw within that
  state's candidates.
- **2D.** `ref_assign` and `ref_assign_sparse` in `simple_eul_prob_tab2D.f90`
  draw the class from the posterior over all active classes: a class in 2D is
  what a state is in 3D, and this is the template.
- **Stages 1-2.** `strategy3D_shc` walks a random order over
  `nspace x nstates` references and takes the first improvement: states and
  directions on an equal footing, the stochastic start the rulings keep
  (ruling 5, 2026-10-21).
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
and the row's losses `L_{s,r}` (state `s`, direction `r`, each entry at its
drawn in-plane and its refined or stored shift, as the fill leaves them), the
tempered posterior over the pair is

    p(s, r | x_i)  proportional to  pi_s exp(-beta T_i (L_{s,r} - L_min)),

with a uniform state prior `pi_s` (no prior) and `beta = prob_beta`. Its state
marginal is `p(s | x_i) = sum_r p(s, r | x_i)`. One draw from the full row is
a draw of `s` from this marginal followed by a draw of `r` within `s`: the
two-step rule is the one-step rule, so the implementation is the single-state
draw over a longer row, and `nstates=1` is the same code on the same random
stream.

**Why the marginal and not the best-per-state softmax of ruling 6.** The two
agree at `beta -> infinity` (both give the argmin's state) and differ at the
search's factor by an entropy term: a state whose map is blurrier has a
flatter block, so it holds more posterior mass for the same best loss. That is
the evidence proper under a uniform prior over directions, and it leans
against the drain: the sharper state no longer wins every near-tie. At the
n_eff the single-state readouts show (1-5 candidates), the difference is
small; the marginal is the form with no rule of its own.

**Why one tempering factor.** The tempering compensates an overconfident
posterior: the temperature is the noise's, while the map at a stage is
band-limited, regularised and, early on, wrong, so the physical posterior
overstates what the data decide. That applies to the state losses as it does
to the direction losses, and the argmin is the `beta = infinity` limit of the
state rule, the limit that locks on and drains. The tail then anneals both to
the physical posterior together in the last stage (`prob_beta_end=1`, set by
the controller), so the final labels are draws from the posterior proper.

**Symmetry breaking.** The balanced random start gives `nstates` nearly
identical mixture maps, a flat state posterior and labels re-drawn almost at
random each iteration. The separation grows from the bias the maps' small
difference puts into the draws, which the next reconstruction amplifies:
stochastic EM. Slower than the argmin's amplification, with no lock-on. On
this set nothing separates under either rule before the band reaches 6 A, so
what the draw changes is in stages 4-5: a lagging state keeps a fair share of
the particles that fit both states nearly alike, and the labels keep moving
while the posterior is flat. The risk runs the other way, that at 0.3 the
state posterior is too flat to break the symmetry within the stage budgets;
section 6 lists the fallbacks and section 7 the readout that decides.

## 4. The changes

### 4.1 The dense assignment (`ref_assign`)

- Replace the argmin-over-states block and the within-block draw with one
  `draw_posterior_row` over the row's `nrefs` entries (the active states'
  blocks); the state is `ref_state(which)`, the direction `ref_proj(which)`.
  A row without a valid distance keeps reference 1, as now.
- The n_eff readout is offered the whole row. A state readout is added
  (section 7).
- The extras: unchanged. `select_extras` reads the committed state from the
  assignment, so they come from the drawn state's block with that state's
  cone, and the weights renormalise within it. The refined-committed-entry
  rule stands.
- `nstates=1`: the same call on the same stream, bit for bit.

### 4.2 Per-state refined set in the fill

`refine_best_refs` and the screen-and-polish of `process_particle_with_map`
select per active state: the `n_s` best of state `s`'s block, with `n_s` the
state's own `calc_num2sample` (already computed per state; today only their
maximum is used). The per-particle refined capacity becomes the sum over
states. The part-file format is unchanged (refined entries are a list with
their reference). Every state's block then carries its own polished mode, and
the state marginal is formed from comparable losses.

Cost: the refinement is 31-41 % of the table's thread time in stages 3-5 with
two states (Phase 7 profile), so this adds up to that much again at two
states and more at four. `n_s` shrinks as a state settles (it follows the
state's mean pose change), which bounds it. The bench record
(`prob_fill_profile`) measures it. A cheaper variant to fall back on: keep the
global top-`projs_ns` and guarantee only each state's block mode a refinement.

### 4.3 The neighbourhood assignment (`ref_assign_neigh`)

`gather_state_row` gathers all valid candidates over all states and the draw
runs once over them; the state is the drawn candidate's. `seed_fallback_if_empty`
and the extras (from the drawn state's candidates) are unchanged. In
`prob_neigh_mode=state` the neighbourhood is the same for every state, so the
candidate sets are comparable; `geom` likewise. This changes `refine3D_states`
too (it runs `prob_neigh` in state or geom mode from `flex_pca` labels): one
rule for the tables. The `multistate_workflows` test covers it.

### 4.4 Convergence

The pose-overlap bypass for `prob`/`prob_neigh` stands. The annealing tail's
stall test reads the best populated state's resolution; under the state draw
it should read every populated state (no state gaining over the margin's
iterations), so the tail does not start while a lagging state is still
climbing. `STATS_FILE` carries `RESOLUTION_FSC05_STATEnn` per state already.

### 4.5 What does not change

- The stage schedule: `shc` in 1-2, `prob` in 3-5, `prob_neigh` state mode
  from 6; `MAXITS`; the 2,500-direction cap; the balanced stochastic sampling
  from stage 4; NU at stage 5; no stage omission or resolution exit for
  multi-state (ruling of 2026-10-21).
- `prob_beta` and `prob_beta_end` as the controller sets them (flat 0.3, the
  last stage's tail to 1).
- Hard labels everywhere: one half-map per particle, `calc_3Drec`'s grouping
  by `(state, half)`, the trailing accounting, the final reconstruction from
  the last iteration's labels, the greedy missing-update pass. `m_estimator=flex`
  stays refused on multi-state `refine3D`; the state weight set is `flex_pca`'s
  and nothing here reads or writes it.
- Empty states: an emptied state (below `MIN_POP`) leaves the table and does
  not return; ruling 7's balanced relabelling at a stage start is a separate
  item, not changed by this proposal. No sampling or relabelling is driven by
  the maps (ruling of 2026-10-04).
- No new user input. The A/B control is the base tree built from a `git
  archive` export, as the regression runs did.

## 5. Tests

- `simple_eul_prob_tab_tester`: `check_posterior_assignment` loses "the state
  is the argmin's state under the draw" and gains: at a huge temperature the
  argmin's state (the old rule's limit); at a finite temperature on a
  synthetic two-state row whose blocks hold known masses, the state
  frequencies over seeds within tolerance of the marginal; ties as now;
  `nstates=1` draws identical to the current build's at the same seed (the
  stream contract). `test_ref_assign_neigh_posterior` likewise over pooled
  candidates. `test_insert_extras_multistate`: every extra in the drawn state's
  block, unchanged. A new check that every active state's block carries a
  refined entry per particle after the fill.
- Identity: single-state solve3D and refine3D bit for bit (`test=identity`);
  the fast gate; the tables' unit entry.
- Behaviour, on the two-state set (P2, the fixed command and the fixed
  threshold: separated at 0.90, partial 0.75-0.90): eight runs of the new
  build against the base tree's 1 separated and 1 partial of 4, judged on the
  count and on the signatures (does the minority still drain in stages 4-5?
  do labels still freeze after a stage-3 lock-on?), with `print_states=yes`
  for per-iteration accuracy. Then the three-state profile
  (`multistate_workflows suite=shc nseeds=20`, 13 of 20 the last count) and a
  `refine3D_states` run (the pooled neighbourhood draw). Timing per stage
  against Phase 7's (stage 5 is 66-91 minutes there; the per-state refinement
  is the only cost this proposal adds).

## 6. Risks and fallbacks

1. **Too hot for states.** The state n_eff stays near `nstates` through
   stage 4 and the maps never diverge. Fallbacks, in order: the tail to 1 at
   the end of every multi-state `prob` stage, not only the last (the
   controller sets `prob_beta_end=1` on stages 3-5), which commits the labels
   at each stage boundary while keeping the exploration within the stage; a
   state factor of its own (an internal key beside `prob_beta`, not a UI
   input), swept on the two-state set; more iterations for multi-state
   stages 4-5. None in the first build.
2. **Cost of the per-state refinement** beyond what the stage can afford: the
   cheaper variant of 4.2.
3. **A state still drains** when its map is garbage (every particle fits the
   other state better at any temperature): ruling 7's reseeding, outside this
   proposal; a population floor in the prior would be map-driven and is ruled
   out.
4. **Ambiguous particles' final labels** are draws from the physical posterior
   after the tail, as the directions are; a particle that fits two states
   alike lands in one of them at random. That is the honest hard label, and
   the mode would be the argmin the tables no longer take.

## 7. The readout

Per table assignment, beside `PROB_TAB_ASSIGN` and the n_eff line:
`>>> PROB_TAB_STATE_DRAW: n_eff over states <draw> / <beta 1>, <pct> % at the
state mode, mean max state posterior <p>`, the state n_eff being
exp(entropy) of the marginal, over the same subsample of rows as the n_eff
readout; the populations per state are already printed by the convergence
object (`STATE n JOINT DISTRIBUTION OVERLAP ... POPULATION`) and `STATE
OVERLAP` (`mi_state`, the fraction keeping their state) is the per-iteration
flow. `STATE_NEFF` goes to `STATS_FILE`. These lines, read beside the
per-iteration accuracy on the two-state set, say whether 0.3 is too hot
(n_eff near `nstates` through stage 4 while accuracy stays at 0.5), too cold
(n_eff near 1 from stage 3's first iterations, the lock-on signature) or
right.

## 8. Decisions to take before coding

1. The joint draw (the marginal) rather than ruling 6's best-per-state
   softmax: section 3's argument.
2. One tempering factor for states and directions, with the state factor and
   the per-stage tail as fallbacks only.
3. The per-state refined set as the default, with its cost measured and the
   cheap variant as the fallback.
4. The same rule in the neighbourhood table, hence in `refine3D_states`.
5. The tail's stall test over every populated state.
6. Stages 1-2 untouched, and no stage-schedule change in the first build.

## 9. Order of work

1. `ref_assign`: the joint draw, the state readout, the unit tests, the
   `nstates=1` identity.
2. The per-state refined set in both fill paths, with the bench record.
3. `ref_assign_neigh`.
4. The stall test over states.
5. The runs of section 5; then the policy documents (`solve3D_policy.md`
   section 6, `refine3D_policy.md`'s assignment sections, the refine3D
   skill's `bayesian-model.md`, and this note moved to `completed/` with the
   results).
