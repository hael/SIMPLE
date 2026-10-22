# Stream 2D pool (p06) at a glance

How the stream classifies all its particles in 2D, for readers who know the technique but not the
code. p06 keeps one growing pool of every particle the sieve accepts and classifies it one 2D
iteration at a time, shows the class averages in the GUI, writes snapshots of the classes the user
selects, and publishes classified particles for multistate 3D (p07). Contract and details:
`pool2D_policy.md` and `stream_3D_ingestion_policy.md` (sections 2, 3 and 9); this page
summarises the code.

```mermaid
flowchart TD
    IN[/"<b>From the sieve (p05)</b><br/>sets of ~10,000 particles with the sieve's 2D (classes,<br/>class averages, FRCs) and selection; an empty final set at the end"/]

    FIRST{"<b>First set for 3D due?</b><br/>100,000 accepted particles handed off,<br/>500 micrographs' worth, the final set,<br/>or pool iteration 10 done"}
    PUB1[("<b>First publication</b><br/>the sieve's own 2D of every set handed off,<br/>without waiting for a pool iteration")]

    IMP["<b>Import</b> (only between iterations)<br/>new particles join unclassified; the pool works<br/>downscaled to at most 2.67 Å per pixel (box ≥ 128)"]
    WAIT1{"enough particles?<br/>start: ~500 micrographs' worth;<br/>after a pause: 50 (≤ it. 20) or 500 micrographs' worth"}

    subgraph IT ["One 2D iteration (150 classes)"]
        direction TB
        SAMP["<b>Sample</b> up to 500,000 accepted particles<br/>(whole micrographs, new draw every iteration);<br/>never-classified particles get a random class"]
        REF2["<b>Classify</b><br/>low-pass ramp over iterations 1–20,<br/>then no fixed limit; shifts from iteration 6"]
        MSK["<b>After iteration 10</b>: switch to the sieve's mask diameter;<br/>the box grows when the resolution reaches Nyquist"]
        SAMP --> REF2 --> MSK
    end

    PUB[("<b>Publication for 3D</b><br/>after every iteration from 25: the classified particles,<br/>all classes with class averages and FRCs at full sampling;<br/>newest 2 kept")]
    FINAL{"<b>Sieve's final set in?</b><br/>run without pausing to iteration 25,<br/>flag that publication as the final one"}
    SNAP[/"<b>2D snapshots</b> on request<br/>particles of the selected classes, one of the last 5 iterations,<br/>with STAR files and class averages"/]
    GUI[/"<b>GUI</b><br/>class averages with resolution and population after every iteration;<br/>particles imported / accepted; mask, box, resolution"/]
    MASKGUI["<b>User changes the mask diameter</b><br/>applies from the next iteration"]

    IN --> FIRST
    FIRST -->|yes, once| PUB1
    IN --> IMP --> WAIT1
    WAIT1 -->|yes| SAMP
    WAIT1 -->|no: pause| IMP
    MSK --> PUB
    MSK -->|next iteration| IMP
    MSK --> GUI
    MSK --> SNAP
    FINAL -.-> PUB
    MASKGUI -.-> SAMP

    classDef step fill:#e5f5ee,stroke:#1b8a5f,color:#1c2230
    classDef check fill:#fcf2e1,stroke:#b46e0c,color:#1c2230
    classDef ref fill:#d6e3fb,stroke:#2c64c9,stroke-width:2px,color:#1c2230
    classDef out fill:#f1ecfb,stroke:#6b46c1,color:#1c2230
    classDef input fill:#ffffff,stroke:#5a6374,color:#1c2230
    class IMP,SAMP,REF2,MSK,MASKGUI step
    class FIRST,WAIT1,FINAL check
    class PUB1,PUB ref
    class SNAP,GUI out
    class IN input
```

## Worth knowing

- **The pool rejects nothing.** It classifies the particles the sieve accepted; class selection for
  3D happens in p07.
- **Above 500,000 particles each iteration uses a random subset** of about 500,000 (whole
  micrographs, a new draw every iteration); particles outside the sample keep their classes.
- **The first set for 3D comes from the sieve, not the pool**, and in practice appears about when
  the pool starts (both wait for roughly 500 micrographs' worth of particles). It is tried once; if
  it cannot be built, the pool publishes instead after iteration 10 or once 100,000 particles are
  classified.
- **Iterations 11 to 24 publish nothing**, so 3D sees a settled classification.
- **The pool pauses between imports**: it waits for more particles before the next iteration, and
  the sieve's final set lifts the pause until iteration 25, and past 25 until the final
  publication has been made.
- **Late particles never see the coarse low-pass limits**, which end at iteration 20.
- **A failed iteration is retried once**; a second failure stops the stage.
