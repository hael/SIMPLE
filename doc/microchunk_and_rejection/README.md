# Microchunking and Class-Average Rejection

This folder collects the documentation for streaming microchunk class-average rejection and the `model_cavgs_rejection` backend.

## Documents

- `model_cavgs_rejection.md`: current command modes, feature bank, hard rejects, built-in models, classification, analysis output, learning, and promotion.
- `microchunk_rejection_model_integration.md`: current relationship between stream rejection and the shared model backend.
- `microchunk_rejection_policy.md`: current stream policy for class-average rejection, particle cleanup, lifecycle sentinels, and combined outputs.

## Implementation Pointers

- `src/main/cavg_quality`: shared class-average quality backend and `model_cavgs_rejection` implementation.
- `src/main/stream/simple_microchunked2D.f90`: microchunk workflow integration.
- `src/main/stream/simple_cluster2D_rejector.f90`: current microchunk rejection engine.
