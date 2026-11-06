# Why do ferro and another tool disagree?

When ferro and another HGVS tool (Mutalyzer, VariantValidator, biocommons `hgvs`, …) give different answers for the same variant, `ferro arbitrate` explains *why*, grounded in the HGVS specification rather than in ferro's own opinion of itself. It names the spec-compliant tool only when one of its deterministic rules decides the case. When none does — a **Different** verdict whose compliance is `needs_interpretation` — it quotes the governing spec rule and leaves the judgment to you, and an **Inconclusive** result names no winner at all.

```bash
# Auto-fetch Mutalyzer's answer and compare:
ferro arbitrate "NM_000255.4:c.445_446insA" --reference ferro-reference/

# Or paste any other tool's output (works offline, for any tool):
ferro arbitrate "NM_000255.4:c.445_446insA" --reference ferro-reference/ \
  --other-output "NM_000255.4:c.445dup"
```

It reports one of: **Equivalent** (same variant, both spellings valid — or one spelling the spec mandates), **Different** (genuinely different edits, with the governing spec rule quoted), **BasisMismatch** (the two are on different transcript versions/axes — align them first), **OtherUnparseable** (the other tool's output could not be parsed as HGVS), **FerroParseError** (ferro's own output could not be parsed as HGVS — always a ferro bug), or **Inconclusive** (could not adjudicate, with the reason). The verdict comes from an oracle that reduces both answers to reference-anchored edits and shares no code with ferro's normalizer — so it can, and will, find ferro at fault.

If arbitration shows **ferro** is the one that's wrong (including a **FerroParseError**), file a bug in one step. ferro prints the issue title and body, then a URL to a prefilled GitHub issue form you review and submit. It always prints the URL, and offers to open it in a browser (after a `y/N` confirmation) only when `$BROWSER` is set, stdin is a terminal, and `--no-open` is not given. ferro needs no GitHub token to build the form, but GitHub may ask you to sign in before you can submit it:

```bash
# One command (only files when ferro is actually implicated):
ferro arbitrate "<variant>" --reference ferro-reference/ --other-output "<other>" --bug-report

# Or from a saved JSON result:
ferro arbitrate "<variant>" --reference ferro-reference/ --other-output "<other>" --format json \
  | ferro bug-report --from-arbitration -
```

> Requires a prepared reference (`ferro prepare`; see [Reference data](reference-data.md)) — the same one normalization uses. Claude Code users can also ask conversationally ("why do ferro and Mutalyzer disagree on …?") via the bundled `arbitrate-hgvs` skill, which drives these commands and applies the spec rule for you on the cases that need interpretation.

## Scripting `arbitrate`

The exit status says whether ferro delivered a verdict; the output says what the verdict was. A script reads `verdict` (and `compliance`) from `--format json`, and uses the status only to catch the runs that produced no verdict:

| status | meaning |
|---|---|
| 0 | a verdict was delivered: **Equivalent**, **Different**, **BasisMismatch** or **OtherUnparseable** |
| 1 | ferro could not run the arbitration: no usable reference, an input it cannot normalize, a failed fetch of the other tool's output |
| 2 | usage error |
| 3 | **Inconclusive**: frame resolution or the oracle itself failed; `reason` in the output says why |
| 4 | **FerroParseError**: ferro's own output did not parse, always a ferro bug |

The result is printed in full before a non-zero exit, so `ferro arbitrate … --format json > result.json` keeps the JSON for a 3 or a 4 as well. `--bug-report` does not change the status: a **FerroParseError** exits 4 after the report has been printed, because the status reports the arbitration, not whether a report was built.
