# Stream preprocessing (p01) at a glance

How the stream preprocesses movies as they arrive, for readers who know the technique but not the
code. p01 watches the microscope's movie folder, motion-corrects and CTF-estimates every new movie
in batches, picks the first ~1,000 micrographs by segmentation, rejects micrographs that fail the
quality thresholds, and hands finished batches to the later stages. This page summarises the code;
the restart rules are in `restart_policy.md`.

```mermaid
flowchart TD
    IN[/"<b>Movies</b> as the microscope writes them<br/>(.mrc, .mrcs, .tif, .eer; flat folder or EPU layout);<br/>optional EPU metadata: beam shift, tilt group"/]
    WATCH["<b>1. Watch for new movies</b><br/>a movie is taken once untouched for 60 s;<br/>new grid-square folders are picked up while running"]
    GAIN{"<b>Gain reference</b><br/>given, flipped automatically,<br/>or generated from the first 1,000 movies"}
    BATCH["<b>2. Batch</b><br/>sets of 5 movies, one queued job each (16 at once, 4 threads);<br/>a leftover set goes after 3 min with no new movie"]

    subgraph JOB ["3. Per movie, in each job"]
        direction TB
        MC["<b>Motion correction</b><br/>global, then patch-based (~200 Å patches),<br/>falling back to global if the patch fit is poor;<br/>downscaled to 1.3 Å per pixel"]
        DW["<b>Dose weighting</b> when the total dose is given<br/>(CTF is fitted on the unweighted sum)"]
        CTF["<b>CTF estimation</b><br/>patch CTF, 30–5 Å, defocus 0.2–5 µm;<br/>records defocus, astigmatism, fit resolution, ice score"]
        PICK["<b>Segmentation picking</b><br/>first 1,000 micrographs only<br/>(used by the initial analysis, p03)"]
        MC --> DW --> CTF --> PICK
    end

    SEL{"<b>4. Micrograph selection</b><br/>rejected when worse than a threshold:<br/>CTF resolution 10 Å, astigmatism 10, ice score 1.0<br/>(the user can change them in the GUI)"}
    REJ["rejected micrographs stay in the project<br/>as deselected, left out of the STAR file"]
    OUT[("<b>5. Hand-off</b><br/>each finished set with at least one accepted micrograph,<br/>renamed into place as soon as its job is collected")]
    DOWN[/"<b>Read by</b><br/>optics groups (p02), initial analysis (p03),<br/>reference picking (p04)"/]
    IDLE{"no new movie for 15 min<br/>and no job running?"}
    MARK[/"<b>Idle marker</b>: lets the downstream stages<br/>finish their intake (removed when a movie arrives)"/]
    GUI[/"<b>GUI</b><br/>movies seen / processed / rejected, rate per hour;<br/>histograms and time plots of CTF resolution, ice, astigmatism, defocus;<br/>newest thumbnails with power spectra"/]

    IN --> WATCH --> GAIN --> BATCH --> MC
    PICK --> SEL
    SEL -->|accepted| OUT --> DOWN
    SEL -->|rejected| REJ
    OUT --> GUI
    WATCH --> IDLE -->|yes| MARK

    classDef step fill:#e5f5ee,stroke:#1b8a5f,color:#1c2230
    classDef check fill:#fcf2e1,stroke:#b46e0c,color:#1c2230
    classDef fail fill:#fbeaea,stroke:#b93a34,color:#1c2230
    classDef ref fill:#d6e3fb,stroke:#2c64c9,stroke-width:2px,color:#1c2230
    classDef out fill:#f1ecfb,stroke:#6b46c1,color:#1c2230
    classDef input fill:#ffffff,stroke:#5a6374,color:#1c2230
    class WATCH,BATCH,MC,DW,CTF,PICK step
    class GAIN,SEL,IDLE check
    class REJ fail
    class OUT ref
    class DOWN,MARK,GUI out
    class IN input
```

## Worth knowing

- **Nothing is processed until the gain is settled**: with a generated gain, that means waiting for
  the first 1,000 movies.
- **Rejection is one-way and applies from now on.** A threshold changed in the GUI applies to later
  batches; loosening it does not bring back micrographs already rejected, and batches already handed
  on are not revised. The thresholds survive a restart.
- **Astigmatism is relative** (|dfx − dfy| / mean defocus) while its default threshold is 10, so it
  very likely never rejects anything; worth checking.
- **A CTF fit worse than 30 Å is recorded at the threshold** and therefore accepted.
- **A failed job's movies count as rejected** and are not retried in the same run.
- **Two movies with the same file name in different folders**: the second is ignored.
- **An existing preprocessing folder can replace p01** entirely (the master's `dir_preprocess`).
