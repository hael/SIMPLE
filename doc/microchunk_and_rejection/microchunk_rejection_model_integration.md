# Microchunk and Model Rejection Relationship

## Scope

This note records the current relationship between the streaming microchunk rejection path and the shared `model_cavgs_rejection` backend.

## Current Code State

In this tree, `simple_microchunked2D` does not call the `src/main/cavg_quality` backend.

The live stream path is:

1. `collect_and_reject` detects `ABINITIO2D_FINISHED`.
2. `reject_cavgs` reads the chunk project and class averages.
3. `reject_cavgs` instantiates `cluster2D_rejector`.
4. `cluster2D_rejector` applies population, resolution, mask, and local-variance criteria.
5. `reject_cavgs` writes selected/rejected class-average stacks and propagates states into particles.

The shared class-average quality backend remains available through:

- `src/main/cavg_quality`
- `model_cavgs_rejection`

## Shared Concepts

The shared model backend mirrors several microchunk principles:

- population and resolution validity gates;
- foreground geometry and connected-component checks;
- local-variance/detail measurements;
- conservative state propagation from class averages to particles.

The two implementations are not identical:

- microchunk streaming currently uses hard, ordered scalar rejection criteria;
- `model_cavgs_rejection` uses hard validity gates followed by normalized scalar features, learned weights, clustering, and threshold controls.

## Ownership Boundary

No backend-selection flag exists in the current stream code.

The stream lifecycle is owned by `simple_microchunked2D`:

- sentinel handling;
- chunk failure/completion flags;
- selected/rejected stack writing;
- class-to-particle state propagation;
- match finalization and combined output.

The current class-average decision is owned by `cluster2D_rejector`.
