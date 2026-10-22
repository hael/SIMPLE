# Stream reference picking (p04) at a glance

How the stream picks and extracts particles, for readers who know the technique but not the code.
p04 turns the picking references into templates once, picks every accepted micrograph by template
matching, extracts the particles, and hands them batch by batch to particle sieving (p05).
Contract and details: `reference_generation_policy.md` (p04's use of the references); this page
summarises the code.

```mermaid
flowchart TD
    REFS[/"<b>Picking references</b><br/>from the initial analysis (p03), or given by the user"/]
    MICS[/"<b>Micrographs from preprocessing (p01)</b><br/>sets of 5 movies; p01's rejections already applied"/]

    subgraph TPL ["Templates: made once, at the first accepted micrograph"]
        direction TB
        T1["<b>Reduce to 10 views</b><br/>more than 10 references are clustered;<br/>the 10 medoids are kept"]
        T2["<b>Prepare</b><br/>rescaled to the micrograph's pixel size, automasked, centred,<br/>low-passed to 15–30 Å (0.15 × the particle diameter)"]
        T3["<b>Particle size and boxes</b><br/>particle diameter from the masks;<br/>mask diameter = 1.2 × particle;<br/>extraction box ≈ 2 × particle (next FFT-friendly size)"]
        T4["<b>Rotate</b> each in 12 steps of 30°, plus mirrors<br/>(up to 240 templates)"]
        T1 --> T2 --> T3 --> T4
    end

    BATCH["<b>Batch</b><br/>each upstream set becomes one job (8 at once, 8 threads),<br/>accepted micrographs only"]

    subgraph JOB ["Per micrograph, in each job"]
        direction TB
        ROI["<b>Region of interest</b><br/>ice and carbon excluded; nothing picked<br/>if less than 1% of the area is usable"]
        COARSE["<b>Coarse search</b> at 4 Å per pixel, low-passed to 15 Å<br/>correlation with every template, best kept"]
        PEAKS{"<b>Peak selection</b><br/>non-maximum suppression (diameter / 3);<br/>automatic threshold (Otsu, refined)"}
        FINEP["<b>Refine positions</b> at 2 Å per pixel (±12 Å)"]
        EXT["<b>Extract</b> from the motion-corrected sum<br/>boxes fully inside the micrograph only;<br/>inverted, ramp removed, normalised on the background;<br/>per-particle defocus from the patch CTF"]
        ROI --> COARSE --> PEAKS --> FINEP --> EXT
        PEAKS -.->|no picks| NOP["micrograph deselected"]
    end

    MOL[/"<b>Mask diameter file</b><br/>what the sieve (p05) waits for and uses"/]
    OUT[("<b>Hand-off to sieving (p05)</b><br/>each finished set with particles,<br/>renamed into place as soon as its job is collected")]
    IDLE{"preprocessing idle or finished,<br/>nothing new, no job running?"}
    MARK[/"<b>Idle marker</b>: lets the sieve<br/>start its final ingestion"/]
    GUI[/"<b>GUI</b><br/>micrographs submitted / with particles, particles extracted, box;<br/>the templates; newest thumbnails with picks"/]

    REFS --> T1
    MICS --> BATCH
    T3 --> MOL
    T4 --> BATCH
    BATCH --> ROI
    EXT --> OUT
    OUT --> GUI
    BATCH --> IDLE -->|yes| MARK

    classDef step fill:#e5f5ee,stroke:#1b8a5f,color:#1c2230
    classDef once fill:#eaf0fc,stroke:#2c64c9,color:#1c2230
    classDef check fill:#fcf2e1,stroke:#b46e0c,color:#1c2230
    classDef fail fill:#fbeaea,stroke:#b93a34,color:#1c2230
    classDef ref fill:#d6e3fb,stroke:#2c64c9,stroke-width:2px,color:#1c2230
    classDef out fill:#f1ecfb,stroke:#6b46c1,color:#1c2230
    classDef input fill:#ffffff,stroke:#5a6374,color:#1c2230
    class T1,T2,T3,T4 once
    class BATCH,ROI,COARSE,FINEP,EXT step
    class PEAKS,IDLE check
    class NOP fail
    class OUT ref
    class MOL,MARK,GUI out
    class REFS,MICS input
```

## Worth knowing

- **At most 10 views become templates**, and they are made once, at the first micrograph's pixel
  size; nothing remakes them later. The initial analysis's 3D route publishes 50 reprojections,
  which are clustered down to 10.
- **The particle size comes from the references**, not from the user: a given particle diameter is
  ignored, and the mask diameter p04 measures is the one the sieve uses.
- **The extraction box is about twice the particle** and is not raised to a minimum; particles whose
  box crosses the micrograph edge are lost.
- **The picking low-pass is fixed at 15 Å**: the `lp_pick` setting currently has no effect.
- **A failed job is not retried**; its micrographs are only counted in the log.
- **p04 applies no micrograph thresholds of its own**, and a threshold the user changes in p01 later
  does not reach sets p04 has already taken.
