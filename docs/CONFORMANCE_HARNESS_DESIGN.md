# Conformance-harness redesign — annotated, self-correcting, hermetic divergence tracking

Parent umbrella: [#325](https://github.com/fulcrumgenomics/ferro-hgvs/issues/325) (biocommons-normalize fixture burn-down).

Related issues:
- [#471](https://github.com/fulcrumgenomics/ferro-hgvs/issues/471) — stronger diagnostic for FASTA-vs-cdot reconstruction mismatches (subsumed by §3).
- [#472](https://github.com/fulcrumgenomics/ferro-hgvs/issues/472), [#473](https://github.com/fulcrumgenomics/ferro-hgvs/issues/473) — residual bugs surfaced by the supplemental-FASTA wiring (PR #417); the kind of "known bug" this design makes self-correcting.
- [#324](https://github.com/fulcrumgenomics/ferro-hgvs/issues/324) / [#323](https://github.com/fulcrumgenomics/ferro-hgvs/issues/323) — corpus import (biocommons / mutalyzer).

Precedent: the mutalyzer harness (`tests/mutalyzer_normalize_tests.rs::AxisTally::record`) already carries a partial form of the annotation model proposed here. This design promotes that from an optional, mutalyzer-only convenience to the primary mechanism for both corpora.

## Problem

The biocommons (and mutalyzer) conformance suites pin every case where ferro's normalize output differs from the upstream reference implementation. Today those divergences live in a flat, hand-maintained text file — `tests/fixtures/biocommons-normalize/baseline-failures/normalized.txt` — that is **not read by any test and not checked by CI**. That single design choice is the root cause of a recurring class of maintenance failures.

### Symptoms observed (all real, all from the last burn-down cycle)

1. **Intent is conflated.** The file mixes two opposite meanings under one anonymous list:
   - *Accepted divergences* — ferro is intentionally, spec-correctly different from upstream (e.g. [#253](https://github.com/fulcrumgenomics/ferro-hgvs/issues/253) boundary-spanning `del` with `cross=false`; [#404](https://github.com/fulcrumgenomics/ferro-hgvs/issues/404) spanning-`dup` canonicalisation; the #293/Pattern-H case where biocommons throws a Python `AttributeError` and ferro emits the spec-canonical answer).
   - *Known bugs* — ferro is wrong, xfail until fixed (#472, #473).

   A reader cannot tell which row is which, why it is there, or who owns it.

2. **It drifts silently because nothing derives or enforces it.**
   - PR #410 fixed `NM_000051.3:c.1_2insCA` (a code fix) but never pruned its line; the row sat stale for a week.
   - PR #417 (supplemental-FASTA wiring) demoted **13 of 15** `NM_001166478.1` divergences by resolving GenBank-deposited bases — but the committed file still lists all 16, so it now overstates the divergence set by ~13 rows.
   - Three separate PRs (#405, #408, #413) hand-edited the file with inconsistent results.

3. **The authoritative test cannot run in CI.** `axis_normalized` requires a multi-GB reference manifest that exists only on a few machines, so CI skips it. The only artifact contributors see — the flat file — is the least trustworthy one, and there is no merge gate that would catch its drift.

4. **Correctness silently depends on un-pinned reference-data versions.** The "expected" answer for `NG_029146.1:g.6494delG` is `g.6496del`, which depends on `NG_029146.1`'s exact bases (a `GGG` tract at g.6494–6496). A locally-available manifest shipping `NG_029146.**2**` (different sequence, no tract) produces `g.6494del` and a *phantom* divergence, with no signal that the wrong accession version was loaded. Diagnosing one such phantom row cost roughly an hour of forensics. The same class bit `NM_212556.2` vs `.4` and `NM_001166478.1` deposited-vs-reconstructed.

### Root cause, one sentence

The committed divergence record is a **hand-maintained, CI-unchecked, undifferentiated flat list whose correctness depends on out-of-band reference data** — so it conflates intent, drifts silently, and cannot gate merges.

## Goals

- A reviewer can see, for every divergence, **what** ferro does, **why** it differs, and **who owns** closing it.
- A fixed bug **cannot silently linger** in the record — the harness forces its removal.
- A reference-data version mismatch is **detected at load**, never emitted as a phantom divergence.
- The merge-gating conformance test is **hermetic** — it runs in CI on every PR without an out-of-band manifest, and is reproducible by any contributor.
- The committed artifact is **derived and verified**, never hand-edited.

## Non-goals

- Changing ferro's normalization behavior. This is a test-infrastructure redesign only.
- Fixing the specific residual bugs (#472, #473) — they are the *clients* of the new model, not part of it.
- Replacing the full manifest-mode sweep. It remains, demoted to a scheduled (non-PR) tier (§5).

## Design

Five pillars. (1) and (2) are the high-value core and are independently shippable; (3)–(5) build on them.

### 1. Divergences become first-class annotated data, living with the cases

Each case in `cases.json` carries, **per axis** (direction × cross), an explicit disposition instead of a bare line in a side file:

- `match` *(default)* — ferro must equal the upstream `normalized` value. A mismatch fails the test.
- `accepted_divergence { reason, spec_citation, ferro_output }` — ferro is intentionally and **terminally** different (both forms are spec-allowed); ferro will not converge. The harness asserts ferro produces **exactly `ferro_output`** (the divergence is *pinned*, not merely "anything ≠ upstream"), and surfaces `reason` + `spec_citation` in failure output.
- `known_bug { tracking_issue, ferro_output }` — ferro is wrong, xfail until fixed. The harness asserts ferro still produces the recorded buggy `ferro_output`.
- `improvement { tracking_issue, section, ferro_output }` — ferro's output is valid HGVS but not the spec-*preferred* canonical form; `normalize()` should converge once `tracking_issue` lands. The middle disposition between `accepted_divergence` (terminal) and `known_bug` (wrong). Like `known_bug` it is an xfail and XPASS-guarded, but it requires a `section` citation substantiating why ferro's current form is non-preferred. (Added in #501 for mutalyzer; mirrored into the biocommons harness in #503.)
- `spec_citation { section }` *(mutalyzer harness; documentary)* — records that the upstream expected value in `cases.json` was corrected to ferro's spec-correct value, with the `section` that arbitrates it. Has no effect on tally bucketing (the corrected expectation makes the case a normal `match`).

The flat `baseline-failures/*.txt` file is deleted (or becomes a generated, checked view — §5).

### 2. XPASS detection — the self-correcting mechanism

For a `known_bug` or `improvement` axis, if ferro's output **starts matching** the upstream `normalized` value (or stops matching the recorded `ferro_output`), the harness **fails loudly**:

> `known_bug for NM_001166478.1:c.36_37insTC (#473) now matches upstream — the bug appears fixed. Remove the annotation and demote the row.`

This is the keystone. It converts burn-down from a manual, drift-prone chore into an invariant the test enforces. The exact failure mode that triggered this redesign — a fixed bug (`c.1_2insCA`) lingering stale in the file — becomes **structurally impossible**: a fixed bug cannot pass silently, because passing *is itself a test failure* until the annotation is removed. (This generalizes Rust's `#[should_panic]`/`xfail` XPASS semantics to per-row data.)

### 3. Reference-data version pinning + load-time spot-check

Each case (or the corpus header) records the accession **and version** its expected answers were validated against. When the harness loads a provider, it spot-checks that the provider actually serves that identity — a short reference window or sequence length — and reacts with a loud, specific message on mismatch:

> `manifest serves NG_029146.2 but case expects NG_029146.1 (seq differs at g.6494) — refusing to emit a divergence against the wrong reference version.`

Whether a mismatch is a hard FAIL or a SKIP is **not** left to the call site — it is fixed by the runtime context, so a gate can never silently downgrade itself to skip:

| Context | On version mismatch | Rationale |
| --- | --- | --- |
| **Per-PR hermetic gate (§4)** | **FAIL**, block the PR | The fixture provider commits the exact bases, so a mismatch is a genuine defect (stale fixture or wrong pin), never a missing-data artifact. |
| **Nightly full-data sweep (§5)** | **SKIP**, but record the incident and open/annotate a tracking issue | The full manifest is out-of-band; a mismatch usually means the manifest drifted, which must be triaged rather than silently passed *or* hard-failed against unreviewed data. |
| **Local dev (manifest mode)** | **SKIP** with a warning and reproduction guidance | A contributor's local manifest may legitimately differ; never block local iteration, but make the skip and its reason visible. |

The hermetic gate is the only FAIL path, and it is exactly the path where the bases are committed and reviewed — so the guard is strict precisely where strictness is sound. This is [#471](https://github.com/fulcrumgenomics/ferro-hgvs/issues/471) generalized from a soft warning into a context-typed hard guard, and it would have reduced the `NG_029146.2` forensics to a single line.

### 4. Hermetic, fixture-backed conformance (the deepest fix)

Today `axis_normalized` needs the full manifest, so it never runs in CI. Instead, for each case, capture the **exact reference window** ferro needs — the bases spanning the variant plus a normalization-shift margin, with the relevant CDS/exon metadata — into a small committed fixture, and run the conformance assertion against a fixture-backed provider (the `MockProvider` family already used by the per-issue `tests/issue_*.rs` tests).

Consequences:
- The merge gate runs in CI on every PR, fast, with zero out-of-band data.
- It is reproducible by any contributor on a fresh checkout.
- The reference-version-coupling problem (§3) *disappears for the gating test* because the exact bases are committed and reviewed.

The per-issue `MockProvider` tests (`issue_401.rs`, `issue_402.rs`, …) are already this pattern, hand-built one variant at a time. The work here is a one-time **extraction tool** that generalizes it across the whole corpus — analogous to the existing spec-fixture generator (`examples/generate_spec_fixture.rs`) with its `--check` mode.

### 5. Tier the suite; generate any human-readable view

- **Per-PR CI (hermetic):** fixture-backed conformance (§4) + the existing `MockProvider` regression pins. This is the gate.
- **Scheduled / nightly (full data):** the existing manifest-mode `axis_normalized` against complete reference data — catches anything the windows miss and regenerates the corpus + fixtures.
- **Human-readable summary, if still wanted:** *generated* from the annotations with a `--check` mode (mirroring `generate_spec_fixture`), so CI fails if the committed view is out of date. It can never drift because no human writes it.

## Schema sketch

Per-case, per-axis annotation in `cases.json` (illustrative — exact field names TBD in planning):

```json
{
  "input": "NM_001166478.1:c.36_37insTC",
  "reference": { "accession": "NM_001166478.1", "length": 2486, "spot_check": { "pos": 36, "bases": "..." } },
  "axes": {
    "3prime:cross": {
      "upstream": "NM_001166478.1:c.36_37dup",
      "disposition": "known_bug",
      "tracking_issue": 473,
      "ferro_output": "NM_001166478.1:c.35_36dup"
    }
  }
}
```

```json
{
  "input": "NM_000051.3:c.-2_-1insCA",
  "axes": {
    "3prime:no-cross": {
      "upstream": "NM_000051.3:c.-1_1insAC",
      "disposition": "accepted_divergence",
      "reason": "ferro applies edit-type priority (dup > ins) unconditionally; biocommons suppresses canon paths that cross the axis under cross=false. ferro's spanning-dup form is HGVS-canonical (biocommons emits it under cross=true).",
      "spec_citation": "HGVS DNA — priority of dup over ins",
      "ferro_output": "NM_000051.3:c.-1_1dup"
    }
  }
}
```

### Schema invariants

The loader **validates each annotation against these invariants before any comparison**, and surfaces a validation error (not a silent skip) on violation. Pinning them here prevents divergent implementations during migration:

- **Uniqueness.** `(input, axis)` is a primary key — at most one disposition per axis per case. Duplicate `(input, axis)` entries are a load error.
- **Per-disposition required / forbidden fields:**

  | `disposition` | Required fields | Forbidden fields |
  | --- | --- | --- |
  | `match` *(default; may be implicit — an axis with no annotation is `match`)* | none | `reason`, `spec_citation`, `ferro_output`, `tracking_issue` (a `match` axis has no divergence to describe) |
  | `accepted_divergence` | `reason`, `spec_citation`, `ferro_output` | `tracking_issue` (an accepted divergence is not a bug and has no owning issue) |
  | `known_bug` | `tracking_issue`, `ferro_output` | `reason`, `spec_citation` (the tracking issue is the source of truth for the why) |
  | `improvement` | `tracking_issue`, `section`, `ferro_output` | `reason` (the `section` + tracking issue substantiate it) |
  | `spec_citation` *(documentary)* | `section` | `ferro_output`, `tracking_issue` (no divergence is pinned — the corrected expectation makes it a `match`) |

- **Types.** `tracking_issue` is a positive integer; `ferro_output`, `reason`, `spec_citation`, `section`, and `upstream` are non-empty strings; `disposition` is one of the enum values above.
- **Reference identity.** Each case's `reference` block requires `accession` and `length`; `spot_check.pos` (if present) is a pure **1-based sequence offset** into the pinned reference window — an index into the reference bases, independent of the HGVS `input` text and its coordinate system (it is *not* an interpretation of HGVS expressions, intronic forms, or offset coordinates). It must satisfy `1 <= pos <= length`, and the window's bases must fit: `pos + len(bases) - 1 <= length`. The provider's underlying string is indexed `0`-based, so the loader reads from offset `pos - 1`. `spot_check.pos` and the HGVS coordinates in `input` are independent; they often coincide for simple coding cases — e.g. `pos: 36` aligns with `c.36_37insTC` — but `spot_check.pos` is never derived from or interpreted as an HGVS token.

`AxisTally::record` becomes annotation-aware: it looks up the disposition for `(case, axis)`, **validates the case entry against the invariants above**, compares ferro's actual output against the right expectation, and classifies the result as pass / fail / **xpass** (the §2 loud failure) — surfacing any validation error distinctly from a comparison failure.

## Migration / incremental rollout

No big-bang. Each step is independently valuable and independently reviewable:

1. **Annotation model + XPASS (pillars 1–2).** Add the schema, teach `AxisTally` to honor it with XPASS detection, seed the annotations from the current `baseline-failures` list + the documented dispositions (#253, #404, #472, #473, #400/#417). Delete or auto-generate the flat file. *Highest value, smallest blast radius — and it immediately reconciles the current stale state by construction, rather than via another hand-edit.*
2. **Version spot-check at load (pillar 3 / #471).**
3. **Window-extraction tool + hermetic provider (pillar 4); flip the gating test to it; demote manifest-mode to nightly (pillar 5).**

## Impact on the current burn-down state

- The stale `NM_000051.3:c.1_2insCA` row and PR #417's 13 undropped `NM_001166478.1` rows get reconciled **by the harness during step 1**, not by hand. That is the principled disposition of PR #477's one-line prune: fold it into the annotation migration rather than continue hand-patching the flat file.
- #325's definition of done is restated: the umbrella closes when the corpus is fully annotated (every divergence is `accepted_divergence` or a `known_bug` with an owning issue) and the gating test is hermetic — not when a flat file reaches a hand-counted row total.
- The optional "per-case `accepted_divergence` annotations" follow-up listed in #325 is **superseded and upgraded** by pillars 1–2 (it gains XPASS detection and applies to both corpora).

## Out of scope

- The mutalyzer corpus migration can follow the same model but is tracked separately; this design uses biocommons as the reference implementation and notes the shared `AxisTally` path.
- Hosting/caching the full manifest for the nightly tier (S3/Zenodo-pinned-by-hash) is a reasonable enhancement but not required for the per-PR hermetic gate, which is the correctness win.

## Open questions (resolve during planning)

- Exact annotation schema and whether axes nest under the case or live in a sibling file keyed by `(input, axis)`.
- Window margin for §4 extraction: fixed bp, or derived from each case's maximal shift distance? (Homopolymer/tandem cases shift farthest.)
- Whether `accepted_divergence` should additionally assert the upstream value is unchanged (to detect upstream corpus re-imports that silently alter the comparison).
- Do we keep a generated flat summary at all, or rely solely on `cases.json` diffs in review?
