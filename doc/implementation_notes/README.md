# Implementation-note status

Every implementation note belongs to one lifecycle directory:

- [completed](completed/) — implemented features and finished investigations
- [planned](planned/) — proposed or partially implemented work
- [rejected](rejected/) — retired, superseded, or rejected designs

The directory is authoritative. Historical sections inside a note can retain
their original proposal language when that chronology is useful.

Two conventions:

- A note that mixed finished and open work is split into `<name>_completed.md`
  (under `completed/`, the full record) and `<name>_unfinished.md` (under
  `planned/`, only the open work, linking back to the record).
- `planned/next_release/` holds proposals deferred past the current release;
  the `rejected/` folder is top level only (no nested `planned/rejected/`).

For the cross-cutting audit summary and the most important open work, see the
[developer documentation index](../for_developers/index.md).
