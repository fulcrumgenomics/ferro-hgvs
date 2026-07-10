---
name: arbitrate-hgvs
description: Explain why ferro and another HGVS tool (Mutalyzer, VariantValidator, biocommons hgvs) disagree on a parse/normalize/projection result, and — only when ferro is actually wrong — help file a bug report. Use when the user is confused by a discrepancy between ferro's output and another tool's output for the same variant.
---

# Arbitrate an HGVS disagreement

The user has a variant where ferro's parse/normalize/projection output does
not match another tool's output (or they suspect it might not) and wants to
know which one — if either — is wrong. This skill drives `ferro arbitrate`
and `ferro bug-report` to get a grounded answer instead of a guess, and holds
the line on **not** blaming ferro just because a bug report would be
convenient, and **not** clearing ferro just because it is the home tool.

## Step 0: collect inputs

You need three things before running anything:

1. **The variant** — the HGVS description in question.
2. **A prepared reference directory** (`--reference <DIR>`). This is
   **required** — `ferro arbitrate` cannot normalize, project, or compare
   edited sequences without it. **If the user has no prepared reference**,
   say so plainly: without one, only parse-level reasoning is possible (does
   the string parse, is the syntax well-formed) — there is no way to
   adjudicate a normalization or projection disagreement. Point them at
   `ferro prepare` to build one, and stop here rather than pretending to
   arbitrate a normalization question you cannot actually check.
3. **The other tool's output.** Either:
   - let `ferro arbitrate` auto-fetch it from Mutalyzer (the default —
     confirm with the user that this is what they want, since it makes a
     network call), or
   - the user pastes it via `--other-output <STR>` (e.g. copied from
     VariantValidator or biocommons `hgvs`), optionally labeled with
     `--other-tool <NAME>` (default label: `provided`).

## Step 1: run `ferro arbitrate`

```
ferro arbitrate <VARIANT> --reference <DIR> \
  [--other-output <STR>] [--other-tool <NAME>] [--mutalyzer-url <URL>] \
  --format json
```

Always use `--format json` for this flow (the CLI's `--format text` verdict
line is for humans running it standalone; you want the structured fields).
Capture the JSON — you'll pipe it straight into `ferro bug-report` later if
needed, so keep it around (a file or a variable), not just its rendered
summary.

The JSON `Arbitration` bundle has the fields you branch on:

- `verdict`: `equivalent` | `different` | `basis_mismatch` | `other_unparseable` | `ferro_parse_error` | `inconclusive`
- `compliance`: `ferro` | `other` | `both_wrong` | `needs_interpretation` | `not_applicable`
- `category`: `ferro_correct` | `mutalyzer_correct` | `equivalent` | `both_incorrect` | `unknown`
- `ferro_output`, `other` (`{tool, status, output}`), `ferro_spdi`, `other_spdi`
- `spec_citations[]`: each `{ spec_version, file, heading, excerpt }` — the
  actual quoted spec text governing the decision

## Step 2: branch on `(verdict, compliance, category)`

**Read `compliance` even when `verdict == equivalent`.** This is the
single most important thing to get right in this flow: `equivalent` only
means the two outputs edit the reference identically — it says nothing
about whether the *notation* either tool used is spec-legal. A duplication
spelled correctly by one tool and spelled as a forbidden insertion by the
other reduces to the same edit (`equivalent`) but still has a real
compliance winner. Do not stop reading at `verdict` — always check
`compliance`/`category` too.

| verdict | compliance | category | meaning | action |
|---|---|---|---|---|
| `equivalent` | `not_applicable` | `equivalent` | Same variant; both spellings are independently spec-legal (e.g. two valid 3′-shift forms). | **No bug.** Explain using the quoted `spec_citations` excerpt: same edit, different valid spelling. |
| `equivalent` | `ferro` | `ferro_correct` | Same variant, but the spec requires one particular spelling (e.g. must be `dup`, not `ins`) — ferro used it, the other tool did not. | **No ferro bug.** This is the *other* tool's notation bug. Explain with the excerpt; do not offer `ferro bug-report`. |
| `equivalent` | `other` | `mutalyzer_correct` | Same variant, but ferro's spelling violates the spec's mandated form while the other tool's is correct. | **This IS a ferro bug** (a notation bug, not a wrong-variant bug). Offer `ferro bug-report`. |
| `different` | `ferro` | `ferro_correct` | Genuinely different edits; ferro's is the compliant one. | **No ferro bug.** Explain with the excerpt why ferro's edit is right and the other tool's is wrong. |
| `different` | `other` or `both_wrong` | `mutalyzer_correct` or `both_incorrect` | Genuinely different edits; ferro's is wrong (either the other tool is right, or neither is). | **Ferro is wrong.** Say so plainly, no hedging. Offer `ferro bug-report --from-arbitration - --category <category>`. |
| `different` | `needs_interpretation` | `unknown` | Genuinely different edits, but the deterministic oracle has no predicate to decide which is compliant. | **You must reason it out.** Read the quoted `spec_citations` excerpt yourself and apply it to this specific variant to reach a verdict. Only offer `ferro bug-report` if *your* reasoning concludes ferro is the wrong one — and say so with your reasoning shown, not as a bare assertion. |
| `basis_mismatch` | `not_applicable` | `unknown` | The two outputs are not on a shared basis — different accession, version, or coordinate axis (e.g. ferro on `NM_000059.4`, the other tool resolved `.3`, or one is `c.` and the other `g.`). | **Not a bug.** Explain the mismatch concretely (name the accessions/versions/axes involved) and help the user re-run both tools on the same basis before drawing any conclusion. Never file a bug report from a `basis_mismatch`. |
| `other_unparseable` | `not_applicable` | `unknown` | The other tool's output didn't parse as HGVS (or none was supplied). | State plainly which side failed. If the other tool's string looks like valid HGVS to you, ask the user to re-paste it via `--other-output`; this is not evidence ferro is wrong. |
| `ferro_parse_error` | `not_applicable` | `unknown` | Ferro itself failed to parse its own output — should not normally happen. | If the *original input* variant is valid HGVS, this is itself a ferro bug worth reporting (a parse failure on valid input is never correct behavior) — offer `ferro bug-report`. |
| `inconclusive` | `not_applicable` | `unknown` | The oracle/projection could not complete, so there is no verdict to draw — the `reason` field says why (e.g. the transcript is absent from the reference, or a bare `NM_:c.` intronic variant needs an explicit genomic parent). | **Not a bug.** Read out the `reason`, explain arbitration could not be completed, and stop. If the reason mentions needing an `NG/NC` parent, suggest re-running with the `NC_(NM_):c.…` form. Never file a bug report from `inconclusive`. |

## Step 3: offer to file, only when ferro is actually wrong

When the row above says to offer a bug report, pipe the arbitration JSON in
and let the CLI build the prefilled issue:

```
ferro bug-report --from-arbitration - --category <category>
```

(`<category>` is the `category` value from Step 1 — `ferro_correct`,
`mutalyzer_correct`, `equivalent`, `both_incorrect`, or your own reasoned
conclusion from the `needs_interpretation` case.) This shows the full issue
body in the terminal and requires confirmation before opening a browser; the
user reviews the prefilled title/body and submits it themselves. Never
submit on the user's behalf, and never file when the row says "no bug."

## Honesty guardrails

These are not optional flourishes — they are the point of this skill:

- **Start from the spec excerpt, never from tool identity.** Ferro being
  the home tool is not evidence it is right, and the other tool being
  well-known (Mutalyzer, VariantValidator, biocommons `hgvs`) is not
  evidence it is right either. Every "who is correct" statement must be
  traceable to a quoted `spec_citations` excerpt or to your own explicit
  application of one — never to "ferro usually gets this right" or "tool X
  is the reference implementation."
- **Trust the deterministic parts, reason carefully about the rest.** The
  `verdict` (same/different) and the dup-vs-ins `compliance`/`category`
  come from ferro's normalizer-independent oracle (it reduces both outputs
  to edited-reference-sequence comparisons and a local reference check —
  it never calls ferro's own shuffling/normalization code) — trust these
  as given. The 3′-shift compliance judgment and any `needs_interpretation`
  case are **not** independently verified; they require you to read the
  quoted excerpt and apply it yourself. Show your reasoning so the user can
  check it, rather than presenting a bare conclusion.
- **If the spec is silent or genuinely ambiguous, say `unknown` — do not
  manufacture a verdict.** A forced guess dressed up as a finding is worse
  than admitting the spec doesn't decide this case cleanly.
- **"Different-but-equivalent" is not "ferro is right."** `verdict:
  equivalent` only says the two outputs describe the same edit to the
  reference — it makes no claim about notation. Never report an
  `equivalent` verdict as "ferro is correct" without also checking
  `compliance`/`category` for a notation issue running the other way.
- **When ferro is wrong, say so plainly and route to bug-report without
  hedging.** Don't soften a `mutalyzer_correct`/`both_incorrect` finding,
  and don't let "ferro is usually right" bias the conclusion. The whole
  value of this tool is that it will tell the user when ferro is wrong.
