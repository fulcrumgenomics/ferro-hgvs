# CodeRabbit review checklist — ferro-hgvs

Patterns CodeRabbit has flagged on this repository, mined from its review
history. Use it as a pre-push self-check: these are the misses that have
actually cost us review round-trips, not generic style advice.

Each entry names the **trigger** (what you changed), what to **verify**, and the
PRs where it was caught. A pattern earns a place here once CodeRabbit has raised
it on this repository **and** the trigger generalises — the entry has to describe
a class of future change, not a one-off bug. Every entry cites the PRs it was
caught on, so its evidence is checkable; a recurrence adds a citation to the
existing entry rather than a new entry.

This complements `.coderabbit.yaml` at the repo root, which tunes the review
itself. Nothing here is enforced by CI — it is a checklist, not a lint.

## Checklist

- **New path-bearing manifest field** (a `PathBuf`, `Option<PathBuf>` or
  `Vec<PathBuf>` added to a struct in `src/prepare/manifest.rs`):
  verify it is threaded into the path-enumerating helpers (`for_each_path`,
  `for_each_path_mut`, `deduplicate_paths`) and existence-checked in
  `src/check/mod.rs` like its sibling fields. A field that is parsed but never
  enumerated silently skips validation and path rewriting. Scalar and metadata
  fields (`prepared_at`, counts) are out of scope — the helpers walk paths, and
  threading a non-path through them is a different error.
  Flagged on `#527` (`canonical_overrides`) and `#530` (`protein_fastas`).
  → *registration-completeness*, Major.

- **New error/warning code** (`W####` / `E####`): verify `related_codes` in
  `src/error_handling/codes.rs` cross-links the sibling code **bidirectionally**,
  that the `docs/spec/error_code_audit.md` scope rows are updated (no row still
  says the new behaviour is "out of scope"), and that a strict/lenient/silent
  mode-matrix test exists. Audit-row drift is the usual miss — the code ships,
  the audit table still describes the old behaviour.
  Flagged on `#534`.
  → *doc-drift* + *test-gaps* + *sibling-pattern*.

- **New coordinate-system / axis token** (`g.` / `m.` / `o.` / `c.` / `n.` / `r.`):
  grep every doc comment and `# Errors` section that recites the axis set and
  confirm the new token is listed; add the validator `match` arm and a per-axis
  rejection/acceptance test. The axis set is spelled out in prose in many places,
  so adding a variant without sweeping them leaves the docs quietly wrong.
  Flagged on `#543` (missing `o.` axis in docs).
  → *doc-drift* + *registration-completeness* + *test-gaps*.

- **New normalize special-position path** (pter/qter/cen resolution in
  `src/normalize/mod.rs`): verify resolved endpoints are guarded against contig
  length and reversed order **before** the window-offset subtraction — the `u64`
  arithmetic underflows on reversed special ranges. That panics under
  `debug_assertions`, but `[profile.release]` sets no `overflow-checks`, so a
  release build wraps to a near-`u64::MAX` offset instead of crashing. Do not
  rely on the panic to surface it.
  Flagged on `#526`.
  → *boundary* / *silent-failure*, Major.

## Maintaining this file

Mine new entries from CodeRabbit's own review history on this repo. Before
adding one, confirm the paths and symbols it names still exist — a checklist
entry that points at a moved module is worse than no entry, because it sends the
reader somewhere irrelevant and erodes trust in the rest of the list.

Entries are written against this repo's real paths and symbols on purpose;
that specificity is what makes them useful. Re-verify them when the modules they
name are refactored.
