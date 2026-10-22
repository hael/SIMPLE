# Stream master (p00) at a glance

How the master runs a streaming session, for readers who know the technique but not the code. The
master does no image processing: it starts the seven processing stages as separate processes,
passes their progress to the GUI (NICE), and passes the user's actions back to them. The stages
hand data to each other through folders, not through the master. Contract and details:
`ipc_policy.md` and `restart_policy.md`; this page summarises the code.

## The pipeline the master starts

```mermaid
flowchart LR
    MOV[/"Movies<br/>(as the microscope writes them)"/]
    P01["<b>p01 Preprocessing</b><br/>motion correction, CTF,<br/>initial picking"]
    P02["<b>p02 Optics groups</b>"]
    P03["<b>p03 Initial analysis</b><br/>makes the picking references"]
    P04["<b>p04 Reference picking</b><br/>and extraction"]
    P05["<b>p05 Particle sieving</b><br/>2D in chunks, rejection"]
    P06["<b>p06 2D pool</b><br/>global 2D classification"]
    P07["<b>p07 Multistate 3D</b>"]
    MOV --> P01
    P01 --> P02
    P01 --> P03
    P01 --> P04
    P03 -->|picking references| P04
    P04 --> P05 --> P06 --> P07
    P02 -.->|optics map| P04
    P02 -.->|optics map| P05
    P02 -.->|optics map| P06

    classDef stage fill:#e5f5ee,stroke:#1b8a5f,color:#1c2230
    classDef opt fill:#eaf0fc,stroke:#2c64c9,color:#1c2230
    classDef input fill:#ffffff,stroke:#5a6374,color:#1c2230
    class P01,P03,P04,P05,P06,P07 stage
    class P02 opt
    class MOV input
```

All seven start at once; each waits for its upstream folder to fill. Preprocessing is skipped
when the user points at an existing preprocessing folder, and the initial analysis when the user
supplies picking references.

## What the master does

```mermaid
flowchart TD
    START["<b>Start</b> (from NICE)<br/>movies, gain, microscope and CTF settings, pixel size, dose,<br/>optional picking references or existing preprocessing"]
    RES["<b>Resources</b><br/>threads and parts per stage from one table<br/>(environment variables can override);<br/>a persistent worker server when the queue is a worker queue"]
    FORK["<b>Start the stages</b><br/>p01, p02, p04, p05, p06, p07, then p03,<br/>each in its own process with its own log"]
    FAIL0{"every stage<br/>running?"}
    ABORT["stop everything,<br/>kill after 60 s, exit with an error"]
    LISTEN["<b>Listener</b> (every 10 ms)<br/>collects each stage's progress messages"]

    subgraph LOOP ["Main loop: every 5 s"]
        direction TB
        HB["<b>Heartbeat to the GUI</b><br/>each stage's state (running, failed, finished, skipped)<br/>and its content (plots, images, volumes),<br/>only what changed since the last accepted heartbeat"]
        ANS["<b>GUI answer</b><br/>the user's actions"]
        ACT["<b>Act on them</b><br/>stop or restart one stage;<br/>CTF / astigmatism / ice thresholds → p01;<br/>reference selection → p03;<br/>2D mask diameter, 2D snapshots → p06;<br/>3D snapshots → p07"]
        HB --> ANS --> ACT --> HB
    end

    STOPQ{"stop requested?<br/>(GUI, SIGTERM or Ctrl-C)"}

    subgraph STOP ["Stop sequence"]
        direction TB
        S1["1. optics groups first (up to 60 s),<br/>so every stage gets the last map"]
        S2["2. preprocessing, initial analysis<br/>and sieving together"]
        S3["3. the large stages one at a time:<br/>3D, then 2D pool, then reference picking"]
        S4["a stage still running 600 s after being asked is killed<br/>and its queued jobs cancelled"]
        S5["4. last heartbeat (up to 3 tries), then exit"]
        S1 --> S2 --> S3 --> S5
        S3 -.-> S4
    end

    START --> RES --> FORK --> FAIL0
    FAIL0 -->|no| ABORT
    FAIL0 -->|yes| LISTEN --> HB
    ACT --> STOPQ
    STOPQ -->|no| HB
    STOPQ -->|yes| S1

    classDef step fill:#eaf0fc,stroke:#2c64c9,color:#1c2230
    classDef loop fill:#e5f5ee,stroke:#1b8a5f,color:#1c2230
    classDef check fill:#fcf2e1,stroke:#b46e0c,color:#1c2230
    classDef fail fill:#fbeaea,stroke:#b93a34,color:#1c2230
    classDef out fill:#f1ecfb,stroke:#6b46c1,color:#1c2230
    class START,RES,FORK,LISTEN step
    class HB,ANS,ACT loop
    class FAIL0,STOPQ check
    class ABORT,S4 fail
    class S1,S2,S3,S5 out
```

## Worth knowing

- **The master never ends on its own.** A session ends only when the user stops it (or the
  process gets SIGTERM / Ctrl-C). When every stage has finished, the master keeps sending
  heartbeats until then.
- **A stage that dies is not restarted automatically.** It shows as failed and the whole stream as
  "error"; the user can restart it from the GUI, and the stage restores what it can from its folder.
- **Only changes are sent.** The heartbeat carries only content that changed; if the GUI does not
  accept one, the next sends everything again.
- **Reference picking stops last**, after sieving, so it can keep picking while the rest stops.
- **A stage killed at the 600 s limit shows as failed**, which makes the session's final status
  "failed" even if everything else finished.
- **NICE marks a stream failed after a minute without a heartbeat**, which is why the last
  heartbeat is retried.
