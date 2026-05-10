# Error code audit against the HGVS spec

This document maps every error and warning code defined in
[`src/error_handling/registry.rs`](../../src/error_handling/registry.rs) to the
section of the [HGVS nomenclature spec](https://hgvs-nomenclature.org/stable/)
it enforces, and identifies spec-mandated rules for which no code currently
exists.

Audited against the vendored spec at
[`assets/hgvs-nomenclature/`](../../assets/hgvs-nomenclature/) (version 21.0,
submodule pinned to `a32f970`).

The companion fixture [`tests/fixtures/error_code_audit.json`](../../tests/fixtures/error_code_audit.json)
mirrors this table machine-readably; the test
[`tests/error_code_audit.rs`](../../tests/error_code_audit.rs) asserts the
fixture and the registry stay in sync (every registry code must appear in the
fixture; every fixture code must exist in the registry) and that any W-row
classified `Enforced` corresponds to an `ErrorType::<Variant>` reference
emitted from a non-bookkeeping file under `src/`.

Tracking issue: [#81](https://github.com/fulcrumgenomics/ferro-hgvs/issues/81)
item **L1**. Follow-up gaps (codes the spec implies but ferro doesn't emit)
are tracked in [#115](https://github.com/fulcrumgenomics/ferro-hgvs/issues/115)
and the soft-validation issues
[#124](https://github.com/fulcrumgenomics/ferro-hgvs/issues/124),
[#125](https://github.com/fulcrumgenomics/ferro-hgvs/issues/125),
[#127](https://github.com/fulcrumgenomics/ferro-hgvs/issues/127),
[#128](https://github.com/fulcrumgenomics/ferro-hgvs/issues/128) under #81 L2.

## Status legend

- **Enforced** — the code reliably surfaces every input that violates the
  cited spec rule; no known false negatives at the boundary the code claims.
- **Partial** — the code surfaces a subset of inputs that violate the cited
  rule (e.g., one syntactic form but not all), or surfaces them only in some
  modes / coordinate systems. Concrete gap noted in the row.
- **Registered** — the code is present in `src/error_handling/registry.rs`
  and the `ErrorType` enum, but no preprocessor or parser site emits it. The
  spec rule is unenforced at runtime under this code (some Registered rows
  *are* enforced under a different code — e.g. `c.0A>G` is rejected as
  `E1003 InvalidPosition`, not as `W4002 PositionZero`; that overlap is
  noted in the row). Distinct from **Missing** because the code already
  exists in the registry; only the emission wiring is absent.
- **Missing** — the code is documented in the audit but does not exist in
  the registry yet. (Currently no rows use this status; reserved for future
  audits where the registry should grow.)
- **Infra** — the code covers an infrastructure failure (I/O, JSON, missing
  reference data) rather than a spec rule; included for completeness, no spec
  section applies.

## Error codes (E-prefix)

| Code | Description | Spec section | Status |
|------|-------------|--------------|--------|
| E1001 | InvalidAccession — accession does not match a recognized RefSeq/Ensembl/LRG prefix | [background/refseq.md — Reference sequence accession formats](https://hgvs-nomenclature.org/stable/background/refseq/) (lines 45–57: NC\_, NT\_, NW\_, NG\_, NM\_, NR\_, NP\_, LRG\_) | Enforced |
| E1002 | UnknownVariantType — coordinate-type prefix (g./c./n./r./p./m./o.) missing or invalid | [background/refseq.md — Sequence types and prefixes](https://hgvs-nomenclature.org/stable/background/refseq/) (lines 64, 81–141: c., g., m., n., o., p., r.); [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) (lines 20–28) | Enforced |
| E1003 | InvalidPosition — position is not a valid integer / offset position | [background/numbering.md — Position numbering](https://hgvs-nomenclature.org/stable/background/numbering/) (lines 16–28: c. numbering rules); [recommendations/DNA/](https://hgvs-nomenclature.org/stable/recommendations/DNA/substitution/) | Enforced |
| E1004 | InvalidEdit — edit type or format is invalid | [recommendations/DNA/substitution.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/substitution/), [deletion.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/deletion/), [insertion.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/insertion/), [duplication.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/duplication/), [inversion.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/inversion/), [delins.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/delins/) | Enforced |
| E1005 | UnexpectedEnd — input ended before description was complete | [recommendations/grammar.md](https://hgvs-nomenclature.org/stable/recommendations/grammar/) (full grammar rules) | Enforced |
| E1006 | UnexpectedChar — invalid character at this parser position | [recommendations/grammar.md](https://hgvs-nomenclature.org/stable/recommendations/grammar/); also catches retracted `c.IVS` notation per [background/numbering.md](https://hgvs-nomenclature.org/stable/background/numbering/) line 32 | Partial — generic `UnexpectedChar` is also raised for the retracted `c.IVS` notation, which the spec calls out specifically as "should not be used"; a dedicated code with an actionable hint is filed in [#115](https://github.com/fulcrumgenomics/ferro-hgvs/issues/115) item 2 |
| E1007 | InvalidBase — invalid nucleotide base (must be A/C/G/T/U or IUPAC code) | [background/standards.md — Nucleotide standards](https://hgvs-nomenclature.org/stable/background/standards/) (line 15: NC-IUB / IUBMB nucleotide codes) | Enforced |
| E1008 | InvalidAminoAcid — invalid amino acid (must be valid 1- or 3-letter IUPAC code) | [background/standards.md — Amino Acid Descriptions](https://hgvs-nomenclature.org/stable/background/standards/) (lines 210–214: IUPAC-IUB amino acid table) | Enforced |
| E2001 | ReferenceNotFound — accession not present in loaded reference data | n/a (infrastructure: depends on which reference set is loaded) | Infra |
| E2002 | SequenceNotFound — sequence data unavailable for accession with metadata | n/a (infrastructure) | Infra |
| E2003 | ChromosomeNotFound — chromosome / contig not in reference | n/a (infrastructure) | Infra |
| E3001 | PositionOutOfBounds — position beyond the length of the reference sequence | [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) line 15: "the reference sequence used must contain the residue(s) described to be changed" | Enforced |
| E3002 | ReferenceMismatch — stated reference base does not match the actual reference | [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) line 15 (residues must exist as described); [recommendations/DNA/substitution.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/substitution/) (substitution semantics) | Enforced |
| E3003 | InvalidRange — start > end or otherwise invalid coordinate range | [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) line 65: "`_` (underscore) is used to indicate a range" (range semantics imply start ≤ end) | Enforced |
| E3004 | ExonIntronBoundary — variant spans an exon/intron junction | [background/numbering.md — Coding DNA numbering](https://hgvs-nomenclature.org/stable/background/numbering/) (lines 21–28: intronic offsets are reckoned per-side); [recommendations/DNA/deletion.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/deletion/) (boundary cases) | Enforced |
| E3005 | UtrCdsBoundary — variant spans the UTR/CDS boundary | [background/numbering.md](https://hgvs-nomenclature.org/stable/background/numbering/) (UTR uses `-`/`*` markers, CDS uses bare integers; boundary requires special handling) | Enforced |
| E4001 | IntronicVariant — intronic variant cannot be normalized without genomic context | [background/numbering.md](https://hgvs-nomenclature.org/stable/background/numbering/) line 28: "a coding DNA reference sequence does not contain intron … sequences and can therefore not be used as a reference to describe variants in these regions"; [background/refseq.md](https://hgvs-nomenclature.org/stable/background/refseq/) (RNA/cDNA refs do not carry intronic sequence) | Enforced |
| E4002 | UnsupportedVariant — variant type not supported for this operation | n/a (implementation limit; surfaces incomplete coverage rather than a spec rule) | Infra |
| E5001 | ConversionFailed — coordinate conversion between reference systems failed | [background/refseq.md](https://hgvs-nomenclature.org/stable/background/refseq/) (cross-reference c./g./n./p. relations) | Enforced |
| E5002 | NoOverlappingTranscript — no transcript overlaps the genomic position | n/a (infrastructure: depends on transcript set loaded) | Infra |
| E9001 | IoError — file I/O failure | n/a (infrastructure) | Infra |
| E9002 | JsonError — JSON parse failure | n/a (infrastructure) | Infra |

## Warning codes (W-prefix)

| Code | Description | Spec section | Status |
|------|-------------|--------------|--------|
| W1001 | LowercaseAminoAcid — `val` instead of `Val` | [recommendations/protein/substitution.md](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/); [background/standards.md](https://hgvs-nomenclature.org/stable/background/standards/) (3-letter code preferred); [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) lines 41–43 | Registered — `correct_amino_acid_case` in `corrections.rs` returns this `ErrorType` but the preprocessor never invokes it; `NP_003997.1:p.val600glu` is auto-normalized in lenient mode with an empty corrections list. Tracked in [#124](https://github.com/fulcrumgenomics/ferro-hgvs/issues/124) (#81 L2 SVA-001) |
| W1002 | SingleLetterAminoAcid — `V600E` instead of `Val600Glu` | [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) lines 41–43: "three-letter amino acid code is preferred"; [recommendations/protein/substitution.md](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/) | Registered — `NP_003997.1:p.V600E` is auto-expanded to `p.Val600Glu` in lenient mode with no W1002 in the corrections list. Tracked in [#124](https://github.com/fulcrumgenomics/ferro-hgvs/issues/124) (#81 L2 SVA-002) |
| W1003 | LowercaseAccessionPrefix — `nm_000088.3` instead of `NM_000088.3` | [background/refseq.md](https://hgvs-nomenclature.org/stable/background/refseq/) lines 45–57 (accession prefixes shown in uppercase) | Enforced |
| W1004 | MixedCaseEditType — `Del` / `INS` instead of `del` / `ins` | [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) lines 87–104 (edit-type tokens shown in lowercase: `del`, `ins`, `dup`, `inv`, `delins`) | Registered — `correct_edit_type_case` in `corrections.rs` returns this `ErrorType` but no preprocessor phase calls it; `NM_000088.3:c.100Del` produces a generic parse error in lenient mode with no W1004 emission |
| W2001 | WrongDashCharacter — en-dash / em-dash instead of ASCII `-` | [recommendations/grammar.md](https://hgvs-nomenclature.org/stable/recommendations/grammar/) (ASCII grammar); [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) line 114 (minus sign semantics) | Enforced |
| W2002 | WrongQuoteCharacter — smart quotes instead of ASCII `"` | [recommendations/grammar.md](https://hgvs-nomenclature.org/stable/recommendations/grammar/) (ASCII grammar) | Partial — emitted by `correct_quote_characters` (preprocessor Phase 3), but standard HGVS does not use ASCII quotes inside any current production rule, so the catch-all parse failure shadows W2002 for most realistic inputs |
| W2003 | ExtraWhitespace — extraneous spaces in expression | [recommendations/grammar.md](https://hgvs-nomenclature.org/stable/recommendations/grammar/) (no whitespace in production rules) | Enforced |
| W2004 | InvalidUnicodeCharacter — generic non-ASCII character | [recommendations/grammar.md](https://hgvs-nomenclature.org/stable/recommendations/grammar/) (ASCII grammar) | Registered — real coverage of non-ASCII inputs is via W2001 (en/em-dash) and W2002 (smart quotes); the generic catch-all is registered but no preprocessor or parser site emits it on its own |
| W3001 | MissingVersion — accession without `.N` version | [background/refseq.md](https://hgvs-nomenclature.org/stable/background/refseq/) line 20: "variant descriptions lacking a version number are **not** valid" | Registered — `NM_000088:c.100A>G` (no `.N` version) parses successfully in lenient mode with no W3001 warning. Tracked in [#124](https://github.com/fulcrumgenomics/ferro-hgvs/issues/124) (#81 L2 SVA-017) |
| W3002 | ProteinSubstitutionArrow — `p.Val600>Glu` instead of `p.Val600Glu` | [recommendations/protein/substitution.md](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/) (substitution syntax; no `>` at protein level) | Enforced |
| W3003 | OldSubstitutionSyntax — `c.100_102>ATG` instead of `c.100_102delinsATG` | [recommendations/DNA/substitution.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/substitution/) line 16: "substitutions involving two or more consecutive nucleotides are described as deletion/insertion (delins) variants"; [recommendations/DNA/delins.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/delins/) | Registered — `c.100_102>ATG` (no ref bases) hits a generic parse error (no W3003); `c.79_80GC>TT` is silently rewritten to `c.79_80delinsTT` with no warning. The `ErrorType` variant exists but no parse-time hook emits it. Tracked in [#115](https://github.com/fulcrumgenomics/ferro-hgvs/issues/115) item 3 |
| W3004 | OldAlleleFormat — `[c.100A>G;c.200C>T]` instead of `c.[100A>G;200C>T]` | [recommendations/DNA/alleles.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/alleles/) (allele bracket placement) | Enforced |
| W3005 | TrailingAnnotation — `c.459A>G (p.Lys153=)` annotations stripped | n/a (defensive against ClinVar-style trailing protein annotations; tolerant input handling, not a spec rule) | Infra |
| W3006 | MissingCoordinatePrefix — `NC_000001.11:12345A>G` missing `g.` | [background/refseq.md](https://hgvs-nomenclature.org/stable/background/refseq/) line 64: "It is mandatory to indicate the type of reference sequence file using a prefix"; [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) lines 20–28 | Enforced |
| W4001 | SwappedPositions — `c.200_100del` corrected to `c.100_200del` | [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) line 65 (range semantics: start ≤ end implied by the underscore range syntax and reinforced throughout DNA/numbering examples) | Registered — `correct_swapped_positions` returns this `ErrorType` but the preprocessor's nine phases do not invoke it; `NM_000088.3:c.200_100del` produces a generic parse error in lenient mode with no W4001 emission |
| W4002 | PositionZero — `c.0A>G` rejected (no auto-correct) | [background/numbering.md — Coding DNA numbering](https://hgvs-nomenclature.org/stable/background/numbering/) line 16: "numbering starts with c.1 at the A of the ATG translation initiation codon"; lines 19–28 confirm there is no `c.0` (numbering goes `c.-1` → `c.1`) | Registered — rejection happens (`detect_position_zero` runs in preprocessor Phase 1) but the diagnostic uses `ErrorCode::InvalidPosition` (E1003), not W4002. The spec rule is enforced under E1003; W4002 has no distinct emission path |
| W5001 | RefSeqMismatch — stated reference base mismatches actual reference (warning sibling of `E3002`) | [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) line 15: "the reference sequence used must contain the residue(s) described to be changed" | Partial — emitted as `NormalizationWarning::RefSeqMismatch` from `src/normalize/mod.rs` when `validate_reference` detects a mismatch during normalization. Coverage is partial because emission requires loaded reference data; without reference data (the common CLI case) parsing succeeds with no warning. With reference data the sibling error `E3002` is raised in strict mode instead |
| W5002 | OverlapConflictingEdits — two or more cis-allele edits share identical reference bounds (e.g. `g.[100G>A;100A>C]`) | [recommendations/DNA/alleles.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/alleles/) (cis-allele descriptions); [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) line 15: "the reference sequence used must contain the residue(s) described to be changed" | Enforced — detected in `src/normalize/overlap.rs` by `detect_overlap_conflicts` and wired into `normalize_allele` after per-variant normalization. Variant output is preserved; the warning is advisory |

## Spec rules without an error code (gaps)

The audit identified spec rules that are either explicit prohibitions ("must
not", "not allowed") or retracted notations ("should not be used") for which
ferro-hgvs has no targeted error or warning code. The first four rows below
expand the [#115](https://github.com/fulcrumgenomics/ferro-hgvs/issues/115)
inventory; the remaining rows surfaced during the spot-check pass for this
audit and overlap with the upcoming #81 L2 work tracked in
[#124](https://github.com/fulcrumgenomics/ferro-hgvs/issues/124),
[#125](https://github.com/fulcrumgenomics/ferro-hgvs/issues/125),
[#127](https://github.com/fulcrumgenomics/ferro-hgvs/issues/127),
[#128](https://github.com/fulcrumgenomics/ferro-hgvs/issues/128). Per the
process for #81 item L1, **no new codes are added in this PR**; the table
captures the gap and the follow-up.

| Spec rule | Spec citation | Current ferro behavior | Follow-up |
|-----------|---------------|------------------------|-----------|
| Self-cancelling allele constructs ("removing part of a reference sequence replacing it with part of the same sequence are not allowed", e.g. `c.[762_768del;767_774dup]`) | [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) line 47 | Parses silently | [#115](https://github.com/fulcrumgenomics/ferro-hgvs/issues/115) item 1 |
| Retracted `c.IVS` intronic notation (e.g. `c.IVS2+2T>G`) | [background/numbering.md](https://hgvs-nomenclature.org/stable/background/numbering/) line 32: "should not be used" | Generic `E1006 UnexpectedChar` parse error with no actionable hint | [#115](https://github.com/fulcrumgenomics/ferro-hgvs/issues/115) item 2; also overlaps `deprecated::IVS_NOTATION` in `src/hgvs/version.rs` |
| Multi-base substitution with explicit ref bases (e.g. `c.79_80GC>TT`, `c.79GC>TT`) | [recommendations/DNA/substitution.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/substitution/) line 58: "this change can not be described as a substitution like c.79_80GC>TT" | Silently rewritten to `delins` (no warning) or generic parse error | [#115](https://github.com/fulcrumgenomics/ferro-hgvs/issues/115) item 3 (extends `W3003` coverage) |
| `con` (sequence conversion) edits | [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) (mentioned as a spec'd edit type) | Parses silently with no semantic validation; cannot distinguish "unrecognized syntax" from "recognized but unimplemented" | [#115](https://github.com/fulcrumgenomics/ferro-hgvs/issues/115) item 4 |
| Single-position range on `del` (e.g. `c.123_123del`) | [recommendations/checklist.md](https://hgvs-nomenclature.org/stable/recommendations/checklist/) line 27 | Silently collapsed to `c.123del` with no warning | [#127](https://github.com/fulcrumgenomics/ferro-hgvs/issues/127) SVA-008 |
| Single-position range on `dup` (e.g. `c.123_123dup`) | [recommendations/checklist.md](https://hgvs-nomenclature.org/stable/recommendations/checklist/) line 27 | Silently collapsed to `c.123dup` | [#127](https://github.com/fulcrumgenomics/ferro-hgvs/issues/127) SVA-009 |
| One-nucleotide inversion (e.g. `c.100_100inv`, `g.234inv`) | [recommendations/DNA/inversion.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/inversion/) line 15: "the description g.234inv is therefore not allowed" | Silently collapsed (single-position form) or accepted (`g.234inv` form) | [#127](https://github.com/fulcrumgenomics/ferro-hgvs/issues/127) SVA-014 |
| Size suffix on deletion (e.g. `c.123del6`, `g.123del3`) | [recommendations/checklist.md](https://hgvs-nomenclature.org/stable/recommendations/checklist/) line 34; [recommendations/DNA/deletion.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/deletion/) line 66: "not allowed" | Parses silently — the `6` is preserved as a size annotation | [#127](https://github.com/fulcrumgenomics/ferro-hgvs/issues/127) SVA-007 |
| Empty `delins` insert (e.g. `c.100_102delins`) | [recommendations/checklist.md](https://hgvs-nomenclature.org/stable/recommendations/checklist/) line 27 | Generic parse error | [#127](https://github.com/fulcrumgenomics/ferro-hgvs/issues/127) SVA-010 |
| Ambiguous insertion ("the format `c.52insT` is incomplete and not allowed") | [recommendations/checklist.md](https://hgvs-nomenclature.org/stable/recommendations/checklist/) line 26; [recommendations/DNA/insertion.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/insertion/) line 55 | Generic parse error; no targeted hint | New gap (file as #115 follow-up if not absorbed by L2) |
| Insertion size-only (e.g. `c.5439_5430ins6`) | [recommendations/checklist.md](https://hgvs-nomenclature.org/stable/recommendations/checklist/) line 27 | Generic parse error | New gap |
| Mixed-reference-types in allele (e.g. `c.[76A>C];g.[10091C>G]`) | [recommendations/DNA/alleles.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/alleles/) line 22: "should not be used" | Parse error (current allele grammar lacks coverage) | Already partially in #81 item H2 |
| `p.(Met1Val)` start-codon substitution | [recommendations/checklist.md](https://hgvs-nomenclature.org/stable/recommendations/checklist/) line 50: "the description p.(Met1Val) is not allowed" | Parses silently | New gap |
| Embedded whitespace in lenient mode (e.g. `c.100 A>G`) | [recommendations/grammar.md](https://hgvs-nomenclature.org/stable/recommendations/grammar/) (no whitespace in production rules) | `correct_whitespace` strips outer whitespace but not embedded; lenient parse hard-rejects rather than soft-warning W2003 | [#128](https://github.com/fulcrumgenomics/ferro-hgvs/issues/128) SVA-024 |

## Cross-reference: deprecated features in `src/hgvs/version.rs`

`src/hgvs/version.rs` already documents four deprecated HGVS notations under
`deprecated::*`. Each maps to a spec rule this audit also surfaces, but none
is wired to a warning code today. They are listed here so a future PR (likely
under #81 L2) can wire them to existing or new codes:

| Constant | Deprecated form | Spec citation | Audit row / follow-up |
|----------|-----------------|---------------|-----------------------|
| `STOP_AS_X` | `p.Arg97X` | [recommendations/checklist.md](https://hgvs-nomenclature.org/stable/recommendations/checklist/) line 48: "the X should not be used" | [#125](https://github.com/fulcrumgenomics/ferro-hgvs/issues/125) SVA-004 |
| `STOP_AS_STAR` | `p.Arg97*` | [recommendations/checklist.md](https://hgvs-nomenclature.org/stable/recommendations/checklist/) line 48: "'Ter' or '\*' should be used to indicate a translation stop codon; the X should not be used" | [#125](https://github.com/fulcrumgenomics/ferro-hgvs/issues/125) SVA-003 |
| `IVS_NOTATION` | `c.IVS2+2T>G` | [background/numbering.md](https://hgvs-nomenclature.org/stable/background/numbering/) line 32: "retracted and should not be used" | [#115](https://github.com/fulcrumgenomics/ferro-hgvs/issues/115) item 2 |
| `FS_X_NOTATION` | `p.Arg97fsX23` | [recommendations/protein/frameshift.md](https://hgvs-nomenclature.org/stable/recommendations/protein/frameshift/) | [#125](https://github.com/fulcrumgenomics/ferro-hgvs/issues/125) SVA-005 |

## Coordination with #81 L2 (soft-validation)

This audit (#81 L1) covers the **registry** of error and warning codes and
maps each to the spec section it enforces. The companion #81 L2 audit covers
the **soft-validation rules** the spec marks "should" or "should not"
(distinct from the "must" / "not allowed" rules tracked here). Several rows
in this audit are reclassified `Registered` — the code exists in the
registry but is not emitted at runtime. Wiring these up is the work of:

- [#124](https://github.com/fulcrumgenomics/ferro-hgvs/issues/124) — `W1001`,
  `W1002`, `W3001` (registered-but-unused soft-validation warnings; existing
  correctors in `corrections.rs` need wiring to the preprocessor).
- [#125](https://github.com/fulcrumgenomics/ferro-hgvs/issues/125) —
  deprecated stop-codon and frameshift forms (`p.Arg97X`, `p.Arg97*`,
  `p.Arg97fsX23`, `p.Arg97fs*23`).
- [#127](https://github.com/fulcrumgenomics/ferro-hgvs/issues/127) —
  non-canonical input forms (size suffixes, single-position ranges,
  1-nt inversions, redundant repeat labels).
- [#128](https://github.com/fulcrumgenomics/ferro-hgvs/issues/128) —
  embedded-whitespace lenient-mode soft-warn.

Boundary: `Enforced` rows in this audit are hard rejections (E-codes) or
warnings already emitted in lenient/silent mode. `Registered` and `Partial`
rows are exactly the surface where L2 work will land.

## Normalization-output canonicalization gaps

A separate, larger class of gaps — non-canonical inputs that ferro normalizes
without rewriting into the spec's preferred form (e.g. `ins` that should
collapse to `dup`, `delins` that should collapse to `inv`, `c.100A>A` that
should canonicalize to `=`, single-variant alleles emitted with bracket
wrapping) — is **out of scope** for this audit and is tracked under #81 items
A1, A2, A3, A4, A8, B1–B5, C1–C5, D1–D8, E1–E3, F1–F2, G1–G3, H1–H2, I1–I4,
J1–J3, K1–K2 and their downstream follow-up issues.

## How to keep this document in sync

1. Adding a code to `src/error_handling/registry.rs`: add a matching entry to
   `tests/fixtures/error_code_audit.json` (mirror status, spec citation, and
   any follow-up issues) and a row to this document.
2. Removing or renaming a code: update the fixture and document in the same
   commit.
3. Reclassifying a status: update both the fixture row's `status` and the
   markdown table row's status column. The drift-prevention test
   (`tests/error_code_audit.rs`) asserts that any W-row classified
   `Enforced` corresponds to an `ErrorType::<Variant>` substring in a
   non-bookkeeping file under `src/`; if you reclassify *to* `Enforced`
   without wiring an emission site, the test will fail. The inverse also
   holds: a `Registered` row whose `ErrorType` variant is now referenced in
   an emission-relevant file fails the symmetric drift check.
4. The test suite enforces (1)–(3) at test time
   (`cargo nextest run --features dev --test error_code_audit`).
