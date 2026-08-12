# Development `abinitio2D_sgd` commander plan

## Scope

Add an explicitly development-labelled `simple_exec prg=abinitio2D_sgd`
entry for the table-free streaming-SGD 2D-classification workflow. The
existing production `prg=abinitio2D` path must remain compatible and must not
expose SGD-only UI controls. The existing `simple_stream prg=abinitio2D_stream`
workflow is a separate acquisition/pooling workflow and is out of scope.

## Design

```text
simple_exec prg=abinitio2D
    standard UI and conventional/default behavior
    commander_abinitio2D

simple_exec prg=abinitio2D_sgd
    developer UI and explicit SGD controls
    commander_abinitio2d_sgd
        delegates staged work to commander_abinitio2D
```

The wrapper reuses the existing controller and strategies. It does not
duplicate stage setup, class averaging, checkpoint handling, or numerical
search code. Its public name intentionally follows the established
`abinitio2D` capitalization; its private Fortran module name remains
lowercase by local convention.

## Completed implementation

1. `simple_ui_cluster2D.f90` shares ordinary 2D inputs between commands, then
   appends the five developer-only controls solely to `abinitio2D_sgd`:
   `sgd_stage4_mode`, `sgd_diagnostic`, `sgd_eta_shift`,
   `sgd_update_frac`, and `sgd_shift_its`.
2. `simple_commanders_abinitio2d_sgd.f90` validates the Euclidean-only
   development route, rejects `inpl_cont=yes`, selects the internal stream
   route, and owns the checkpoint environment hooks.
3. `simple_exec_cluster2D.f90` dispatches `abinitio2D_sgd` through its own
   commander while preserving standard `abinitio2D` dispatch unchanged.
4. `simple_commanders_abinitio2D.f90` retains the shared staged lifecycle.
   The standard entry has no SGD-specific UI or environment handling. The
   wrapper passes its checkpoint request and a terminal command-line policy;
   that policy removes development stream fields before the conventional
   terminal greedy pass.
5. Focused UI visibility and dispatch tests verify that standard `abinitio2D`
   omits all five developer inputs, while `abinitio2D_sgd` registers the
   developer controls and rejects invalid objective, in-plane, and stage-mode
   combinations.

## Compatibility rules

- Do not rename or remove `abinitio2D`.
- Do not alter `simple_stream` or `abinitio2D_stream`.
- Keep generic parameter parsing available internally for delegated execution.
- Preserve stages 1--3 of the shared staged workflow.
- The wrapper terminal pass must restore conventional greedy behavior by
  removing its stream-only command-line fields before terminal execution.
- `sgd=yes` and `sgd_path=stream` are internal compatibility transport fields;
  users activate the development workflow only through `prg=abinitio2D_sgd`
  and, when desired, `sgd_stage4_mode=on|alternate`.

## Validation record

- 2026-08-04: Plan created after reviewing the UI refactor and streaming-SGD
  development work.
- 2026-08-05: Focused Oracle/Linux validation observed: UI visibility test
  passed 151 assertions, development-dispatch test passed 15 assertions, and
  the six-case base SGD regression suite passed.
- 2026-08-05: Independent fresh 100-class betagal baseline and stream jobs
  both reached `SIMPLE_ABINITIO2D NORMAL STOP`. The stream run used fresh
  60-percent mini-batches in its streamed stages. These independent runs are
  useful operational evidence, but they are not a matched scientific A/B
  comparison because their extracted particle counts differed.

## Remaining validation

1. Rebuild after the public-command capitalization change and rerun the UI and
   dispatch tests; verify `simple_exec prg=list` contains `abinitio2D_sgd`.
2. Run one full `prg=abinitio2D_sgd` workflow after that rebuild, including the
   terminal greedy pass.
3. For a scientific comparison, fork baseline and stream jobs from the same
   stage-3 checkpoint and compare runtime, resolution, support/occupancy, and
   zero-support classes.

For an initial full-workflow comparison, preprocessing through `5_extract`
should be performed once and copied into isolated baseline and SGD roots. This
holds input movies, picks, extracted particles, and the initial project state
fixed while allowing the two ab initio workflows to create their own
`6_abinitio2D` outputs.
