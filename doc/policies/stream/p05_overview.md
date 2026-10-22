# Stream particle sieving (p05) at a glance

How the stream removes junk particles continuously, for readers who know the technique but not the
code. p05 groups the particles from reference picking into chunks, classifies each chunk in 2D,
rejects bad class averages together with their particles, and does this twice (a coarse pass,
then a fine pass on merged chunks) before handing clean sets to the 2D pool (p06). Contract and
details: `doc/policies/sieving_and_rejection/ptcl_sieve_policy.md`; this page summarises the code.

```mermaid
flowchart TD
    IN[/"<b>From reference picking (p04)</b><br/>extracted particle sets, one record per micrograph;<br/>the mask diameter measured from the picking references"/]
    IMP["<b>Import</b> new sets<br/>(untouched for 2 s, up to 20 per pass)"]

    subgraph COARSE ["Coarse pass"]
        direction TB
        CC["<b>Coarse chunk</b><br/>whole micrographs until just over 5,000 particles"]
        C2D["<b>2D classification</b><br/>100 classes, 128 px box, to 15 Å,<br/>2,000 particles sampled per iteration"]
        CREJ{"<b>Reject junk classes</b><br/>hard gates: population, bad pixels, no object,<br/>object outside the mask;<br/>then the quality model (keep if p ≥ 0.25)<br/>and a size-compatibility filter"}
        CC --> C2D --> CREJ
    end

    subgraph FINE ["Fine pass"]
        direction TB
        FC["<b>Fine chunk</b><br/>coarse chunks merged until their accepted<br/>particles pass 10,000"]
        F2D["<b>2D classification</b> of the accepted particles<br/>100 classes, 128 px box, to 10 Å,<br/>10,000 sampled per iteration"]
        FREJ{"<b>Reject junk classes</b><br/>same gates and model as the coarse pass"}
        FC --> F2D --> FREJ
    end

    OUT[("<b>Hand-off to the 2D pool</b><br/>one set per fine chunk (~10,000 particles):<br/>all particles, rejected ones flagged (not deleted),<br/>class table, newest optics groups")]
    FIN{"<b>Final ingestion?</b><br/>reference picking idle or finished<br/>and no new sets since"}
    LAST["leftover particles go as a last chunk<br/>without the coarse 2D and rejection;<br/>then a last fine chunk below 10,000"]
    END[/"<b>Empty final set</b><br/>tells the 2D pool that nothing more comes"/]
    FAILC["a chunk that fails is retried once;<br/>a second failure drops its particles"]
    GUI[/"<b>GUI</b><br/>particles imported / accepted / rejected, failed chunks;<br/>the latest chunk's class averages with resolution,<br/>population and selection"/]

    IN --> IMP --> CC
    CREJ -->|kept classes' particles| FC
    CREJ -.->|rejected: state 0| FC
    FREJ --> OUT
    IMP --> FIN
    FIN -->|yes| LAST --> END
    C2D -.-> FAILC
    F2D -.-> FAILC
    OUT --> GUI

    classDef step fill:#e5f5ee,stroke:#1b8a5f,color:#1c2230
    classDef check fill:#fcf2e1,stroke:#b46e0c,color:#1c2230
    classDef fail fill:#fbeaea,stroke:#b93a34,color:#1c2230
    classDef ref fill:#d6e3fb,stroke:#2c64c9,stroke-width:2px,color:#1c2230
    classDef out fill:#f1ecfb,stroke:#6b46c1,color:#1c2230
    classDef input fill:#ffffff,stroke:#5a6374,color:#1c2230
    class IMP,CC,C2D,FC,F2D,LAST step
    class CREJ,FREJ,FIN check
    class FAILC fail
    class OUT ref
    class END,GUI out
    class IN input
```

Up to 4 chunk jobs run at once (16 threads each), fine chunks first; each 2D has a 29-minute
limit. The stage checks for new sets and finished chunks every 10 s.

## Worth knowing

- **Most particles are classified twice** (in chunks of ~5,000, then ≥ 10,000), but only the fine
  pass is handed on. The leftovers of the final ingestion skip the coarse pass by design.
- **The 2D is deliberately coarse:** a 128-pixel box, 15 Å (coarse) and 10 Å (fine).
- **Rejected particles are flagged, not deleted**, in the hand-offs; the 2D pool only uses the
  accepted ones.
- **The final set waits for reference picking to go idle**, which in turn waits until no new movie
  has arrived for 15 minutes.
- **"Accepted" and "rejected" in the GUI** count only finished fine chunks (and coarse chunks that
  came out empty).
- **A chunk that fails twice drops its particles**; only the GUI's failed-chunk count shows it.
- **The `sieve` quality model still needs a validation run** since it was refitted for the mask
  diameter.
