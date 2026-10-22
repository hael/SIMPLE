# Stream optics groups (p02) at a glance

How the stream assigns optics groups, for readers who know the technique but not the code. p02
groups micrographs recorded with similar beam-image shift (and, if asked, the same EPU tilt
group), and publishes a map from micrograph to group that the downstream stages apply to their
micrographs and particles. This page summarises the code; the restart rules
are in `restart_policy.md`.

```mermaid
flowchart TD
    IN[/"<b>From preprocessing (p01)</b><br/>finished sets of 5 movies with at least one accepted micrograph;<br/>each micrograph carries its EPU beam shift (x, y)<br/>and the tilt group from its file name, when the EPU metadata were given"/]
    IMP["<b>1. Import new micrographs</b><br/>a set is taken once untouched for 2 s, up to 50 per pass;<br/>rejected micrographs are imported too"]
    NEW{"anything new?"}
    WAIT["wait 10 s"]
    SPLIT["<b>2. Split by tilt group</b> (optional, off by default)<br/>one subset per EPU tilt group"]
    LINK["<b>3. Group by beam shift</b><br/>single linkage: two micrographs are linked when their<br/>beam shifts differ by at most 0.05 (EPU units);<br/>each connected set is one optics group"]
    IDS["<b>4. Keep group ids stable</b><br/>a group takes the id most of its micrographs had before;<br/>new groups get new ids; ids are never reused"]
    MAP[("<b>5. Publish an optics map</b><br/>micrograph → group, and the optics table<br/>(pixel size, Cs, kV, amplitude contrast per group);<br/>written complete, numbered, newest 5 kept")]
    READ["<b>Applied downstream</b><br/>reference picking (p04): its micrographs and STAR files;<br/>sieving (p05): each set handed to the 2D pool;<br/>2D pool (p06): publications, snapshots, final project;<br/>a micrograph newer than the map goes to group 1"]
    GUI[/"<b>GUI</b><br/>micrographs imported and accepted, number of groups,<br/>each group's recent beam shifts"/]

    IN --> IMP --> NEW
    NEW -->|no| WAIT --> IMP
    NEW -->|yes| SPLIT --> LINK --> IDS --> MAP
    MAP --> READ
    MAP --> GUI
    MAP -.->|next pass, all micrographs regrouped| IMP

    classDef step fill:#e5f5ee,stroke:#1b8a5f,color:#1c2230
    classDef check fill:#fcf2e1,stroke:#b46e0c,color:#1c2230
    classDef out fill:#f1ecfb,stroke:#6b46c1,color:#1c2230
    classDef ref fill:#d6e3fb,stroke:#2c64c9,stroke-width:2px,color:#1c2230
    classDef input fill:#ffffff,stroke:#5a6374,color:#1c2230
    class IMP,SPLIT,LINK,IDS,WAIT step
    class NEW check
    class MAP ref
    class READ,GUI out
    class IN input
```

## Worth knowing

- **Optics groups are bookkeeping and export.** They end up as `rlnOpticsGroup` in the STAR files;
  no SIMPLE alignment or reconstruction uses them.
- **Without the EPU metadata every beam shift is 0**, so there is a single group.
- **Single linkage chains:** two shift clusters joined by a path of micrographs, each step within
  the threshold, become one group, and rejected micrographs can be the bridge. Within a session
  groups can only grow or merge, so merged groups leave gaps in the ids.
- **"Beam tilt" means the EPU tilt group** taken from the file name (`Data_…`), not a measured
  beam tilt; it is off by default.
- **All micrographs are regrouped on every pass** that imports something; the ids stay stable.
- **On a stop, p02 is stopped first** (up to 60 s) so every other stage gets the last map. After a
  restart it publishes nothing until every micrograph of its last map is back, so readers never
  apply a partial map.
