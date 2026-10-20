# Refactoring-note status

Every refactoring note belongs to one lifecycle directory:

- [completed](completed/) — landed refactors and completed reviews
- [planned](planned/) — open refactors, validation, and test handovers
- [rejected](rejected/) — retired or rejected architectural directions

The directory is authoritative. A completed review may identify follow-up work;
that work belongs in a separate planned record rather than changing the review's
status.

Two conventions:

- A note that mixed finished and open work is split into `<name>_completed.md`
  (under `completed/`, the full record) and `<name>_unfinished.md` (under
  `planned/`, only the open work, linking back to the record).
- `planned/next_release/` holds proposals deferred past the current release, and
  `planned/need2basis/` holds proposals that need a decision or basis before work can start;
  the `rejected/` folder is top level only (no nested `planned/rejected/`).

For the cross-cutting audit summary and the most important open work, see the
[developer documentation index](../for_developers/index.md).
