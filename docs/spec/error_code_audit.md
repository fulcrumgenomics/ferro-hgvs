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
  and the `ErrorType` enum, but no preprocessor or parser site emits it
  under this specific identity. The spec rule is unenforced at runtime
  under this code (some Registered rows *are* enforced under a different
  code — e.g. generic non-ASCII inputs are caught by the dash/quote
  correctors W2001 and W2002 rather than the catch-all W2004; that
  overlap is noted in the row). Distinct from **Missing** because the
  code already exists in the registry; only the emission wiring is
  absent.
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
| E3006 | SelfCancellingAllele — overlapping `del`+`dup` pair within a single cis allele (e.g. `c.[762_768del;767_774dup]`) | [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) line 47: "descriptions removing part of a reference sequence replacing it with part of the same sequence are not allowed" | Enforced |
| E4001 | IntronicVariant — intronic variant cannot be normalized without genomic context | [background/numbering.md](https://hgvs-nomenclature.org/stable/background/numbering/) line 28: "a coding DNA reference sequence does not contain intron … sequences and can therefore not be used as a reference to describe variants in these regions"; [background/refseq.md](https://hgvs-nomenclature.org/stable/background/refseq/) (RNA/cDNA refs do not carry intronic sequence) | Enforced |
| E4002 | UnsupportedVariant — variant type not supported for this operation | n/a (implementation limit; surfaces incomplete coverage rather than a spec rule) | Infra |
| E5001 | ConversionFailed — coordinate conversion between reference systems failed | [background/refseq.md](https://hgvs-nomenclature.org/stable/background/refseq/) (cross-reference c./g./n./p. relations) | Enforced |
| E5002 | NoOverlappingTranscript — no transcript overlaps the genomic position | n/a (infrastructure: depends on transcript set loaded) | Infra |
| E9001 | IoError — file I/O failure | n/a (infrastructure) | Infra |
| E9002 | JsonError — JSON parse failure | n/a (infrastructure) | Infra |

## Warning codes (W-prefix)

| Code | Description | Spec section | Status |
|------|-------------|--------------|--------|
| W1001 | LowercaseAminoAcid — `val` instead of `Val` | [recommendations/protein/substitution.md](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/); [background/standards.md](https://hgvs-nomenclature.org/stable/background/standards/) (3-letter code preferred); [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) lines 41–43 | Enforced |
| W1002 | SingleLetterAminoAcid — `V600E` instead of `Val600Glu` | [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) lines 41–43: "three-letter amino acid code is preferred"; [recommendations/protein/substitution.md](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/) | Enforced |
| W1003 | LowercaseAccessionPrefix — `nm_000088.3` instead of `NM_000088.3` | [background/refseq.md](https://hgvs-nomenclature.org/stable/background/refseq/) lines 45–57 (accession prefixes shown in uppercase) | Enforced |
| W1004 | MixedCaseEditType — `Del` / `INS` instead of `del` / `ins` | [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) lines 87–104 (edit-type tokens shown in lowercase: `del`, `ins`, `dup`, `inv`, `delins`) | Enforced — wired as InputPreprocessor Phase 6b1 via `correct_edit_type_case_full` (#268). Walks the input for case-insensitive edit-type tokens (`del`, `ins`, `dup`, `inv`, `delins`, `con`) anchored after a digit / `]` / `)` / `?` / `*` and lowercases them. Runs **before** Phase 6c (single-letter AA expansion) so `p.Arg8_Lys10DEL` is lowercased to `del` rather than mis-expanded to `AspGluLeu`. Lenient/silent mode emits W1004 and rewrites; strict rejects with `E1004 InvalidEdit` + hint |
| W2001 | WrongDashCharacter — en-dash / em-dash instead of ASCII `-` | [recommendations/grammar.md](https://hgvs-nomenclature.org/stable/recommendations/grammar/) (ASCII grammar); [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) line 114 (minus sign semantics) | Enforced |
| W2002 | WrongQuoteCharacter — smart quotes instead of ASCII `"` | [recommendations/grammar.md](https://hgvs-nomenclature.org/stable/recommendations/grammar/) (ASCII grammar) | Partial — emitted by `correct_quote_characters` (preprocessor Phase 3), but standard HGVS does not use ASCII quotes inside any current production rule, so the catch-all parse failure shadows W2002 for most realistic inputs |
| W2003 | ExtraWhitespace — extraneous spaces in expression | [recommendations/grammar.md](https://hgvs-nomenclature.org/stable/recommendations/grammar/) (no whitespace in production rules) | Enforced |
| W2004 | InvalidUnicodeCharacter — generic non-ASCII character | [recommendations/grammar.md](https://hgvs-nomenclature.org/stable/recommendations/grammar/) (ASCII grammar) | Registered — real coverage of non-ASCII inputs is via W2001 (en/em-dash) and W2002 (smart quotes); the generic catch-all is registered but no preprocessor or parser site emits it on its own |
| W3001 | MissingVersion — accession without `.N` version | [background/refseq.md](https://hgvs-nomenclature.org/stable/background/refseq/) line 20: "variant descriptions lacking a version number are **not** valid" | Enforced |
| W3002 | ProteinSubstitutionArrow — `p.Val600>Glu` instead of `p.Val600Glu` | [recommendations/protein/substitution.md](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/) (substitution syntax; no `>` at protein level) | Enforced |
| W3003 | OldSubstitutionSyntax — `c.100_102>ATG` instead of `c.100_102delinsATG` | [recommendations/DNA/substitution.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/substitution/) line 16: "substitutions involving two or more consecutive nucleotides are described as deletion/insertion (delins) variants"; [recommendations/DNA/delins.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/delins/) | Enforced — wired in preprocessor Phase 9 for the three deprecated multi-base substitution forms (`c.79_80GC>TT`, `c.79GC>TT`, `c.100_102>ATG`). Strict rejects with `E1004 InvalidEdit` + hint; lenient warns and rewrites to `delins`; silent rewrites without warning |
| W3004 | OldAlleleFormat — `[c.100A>G;c.200C>T]` instead of `c.[100A>G;200C>T]` | [recommendations/DNA/alleles.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/alleles/) (allele bracket placement) | Enforced |
| W3005 | TrailingAnnotation — `c.459A>G (p.Lys153=)` annotations stripped | n/a (defensive against ClinVar-style trailing protein annotations; tolerant input handling, not a spec rule) | Infra |
| W3006 | MissingCoordinatePrefix — `NC_000001.11:12345A>G` missing `g.` | [background/refseq.md](https://hgvs-nomenclature.org/stable/background/refseq/) line 64: "It is mandatory to indicate the type of reference sequence file using a prefix"; [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) lines 20–28 | Enforced |
| W3007 | DeprecatedStopCodonStar — `p.Arg97*` rewritten to `p.Arg97Ter` (`Ter` is the preferred stop-codon form) | [recommendations/checklist.md](https://hgvs-nomenclature.org/stable/recommendations/checklist/) (HGVS checklist: "Ter" is the preferred form for the translation stop codon; `*` is permitted but discouraged) | Enforced — `correct_deprecated_protein_forms` returns `ErrorType::DeprecatedStopCodonStar`; emitted by preprocessor Phase 6b in lenient/silent modes. Strict mode rejects the input. Closes [#125](https://github.com/fulcrumgenomics/ferro-hgvs/issues/125) |
| W3008 | DeprecatedStopCodonX — `p.Arg97X` rewritten to `p.Arg97Ter` (`X` is reserved for `Xaa` "any amino acid", not for stop) | [recommendations/checklist.md](https://hgvs-nomenclature.org/stable/recommendations/checklist/) ("the X should not be used" for the translation stop codon); [recommendations/protein/variant/substitution.md](https://hgvs-nomenclature.org/stable/recommendations/protein/variant/substitution/) (canonical `Ter`) | Enforced — `correct_deprecated_protein_forms` returns `ErrorType::DeprecatedStopCodonX`; emitted by preprocessor Phase 6b in lenient/silent modes. Strict mode rejects the input. Closes [#125](https://github.com/fulcrumgenomics/ferro-hgvs/issues/125) |
| W3009 | DeprecatedFrameshiftStar — `p.Arg97fs*23` rewritten to `p.Arg97fsTer23` | [recommendations/protein/variant/frameshift.md](https://hgvs-nomenclature.org/stable/recommendations/protein/variant/frameshift/) (canonical `fsTerN`; `fs*N` is permitted but discouraged) | Enforced — `correct_deprecated_protein_forms` returns `ErrorType::DeprecatedFrameshiftStar`; emitted by preprocessor Phase 6b in lenient/silent modes. Strict mode rejects the input. Closes [#125](https://github.com/fulcrumgenomics/ferro-hgvs/issues/125) |
| W3010 | DeprecatedFrameshiftX — `p.Arg97fsX23` rewritten to `p.Arg97fsTer23` | [recommendations/protein/variant/frameshift.md](https://hgvs-nomenclature.org/stable/recommendations/protein/variant/frameshift/) (canonical `fsTerN`; `X` is not used for HGVS frameshift termination) | Enforced — `correct_deprecated_protein_forms` returns `ErrorType::DeprecatedFrameshiftX`; emitted by preprocessor Phase 6b in lenient/silent modes. Strict mode rejects the input. Closes [#125](https://github.com/fulcrumgenomics/ferro-hgvs/issues/125) |
| W3011 | DelSizeSuffix — `g.123del6` instead of `g.123_128del` | [recommendations/DNA/deletion.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/deletion/): "a deletion of more than one residue should mention the first and last residue deleted, e.g. `NG_012232.1:g.123_128del` and not `NG_012232.1:g.123del6`" | Partial — lenient warns without rewriting (`warn_accept`) because synthesizing a safe end position depends on offset/intronic semantics not available at parse time; strict rejects. Wired as InputPreprocessor Phase 13 via `correct_del_size_suffix` (#127, #81 L2 SVA-007) |
| W3012 | EmptyDelinsInsert — `g.100_102delins` rewritten to `g.100_102del` | [recommendations/DNA/delins.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/delins/) (a `delins` whose inserted sequence is absent collapses to a plain `del`); #81 A3 canonicalization rule | Enforced — lenient/silent rewrite `delins` → `del` and warn once; strict rejects. Wired as InputPreprocessor Phase 11 via `correct_empty_delins` (#127, #81 L2 SVA-010) |
| W3013 | RedundantRepeatLabel — `r.100_102cug[4]` rewritten to `r.100_102[4]` | [recommendations/RNA/repeated.md](https://hgvs-nomenclature.org/stable/recommendations/RNA/repeated/): "the format `r.-125_-123cug[4]`, should not be used; it contains redundant information ('-125_-123' and 'cug')" | Enforced — lenient/silent strip the lowercase `a/c/g/u` label and warn once; strict rejects. Wired as InputPreprocessor Phase 12 via `correct_redundant_repeat_label`, scoped to the `r.` description so a later non-`r.` description in the same input is not stripped (#127, #81 L2 SVA-027) |
| W3014 | DeprecatedIvsNotation — retracted `c.IVS<n>+offset` / `n.IVS` / `r.IVS` intronic notation | [background/numbering.md](https://hgvs-nomenclature.org/stable/background/numbering/) line 32 (IVS notation retracted in favour of canonical `c.<exon-pos>+<offset>`) | Enforced — pre-parse reject; cannot be auto-rewritten without genomic intron metadata. Strict / lenient / silent all reject with `E1006 UnexpectedChar` + hint pointing at the canonical form. `ErrorOverride::Accept` lets a power-user opt out |
| W3015 | DeprecatedConSyntax — `<pos>_<pos>con<source>` rewritten to `<pos>_<pos>delins<source>` | [SVD-WG009](https://www.hgvs.org/varnomen/Conv-conv.html) (sequence conversion `con` retired in favour of `delins`) | Enforced — pre-parse rewrite in lenient/silent modes; strict rejects with `E1004 InvalidEdit` + hint. The legacy `parse_hgvs` parser path still accepts `NaEdit::Conversion` directly; only `parse_hgvs_with_config` / lenient / silent / CLI go through the preprocessor phase |
| W3016 | LengthMismatch — explicit reference sequence length does not match position range | [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) (range semantics: `len(reference) == end - start + 1` for `del` / `dup` / `inv` / `delins` with explicit ref) | Enforced — wired as InputPreprocessor Phase 13a via `detect_length_mismatch`. Lenient mode emits W3016 without rewriting (no safe auto-correction — either endpoint or the ref seq could be wrong); strict mode rejects with `E1004 InvalidEdit`. Detector handles ranges with simple integer endpoints; offset and `*N`/`-N` marker endpoints fall through (cannot compute length without provider). Covers `del<ref>`, `dup<ref>`, `inv<ref>`; `delins<ins>` is skipped (the inserted sequence is intentionally allowed to differ in length) |
| W3017 | AlleleFractionAnnotation — trailing `[level=NN%]` / `[heteroplasmy=NN%]` / `[mosaic=NN%]` / `(NN%)` rejected with a targeted SVA diagnostic | [recommendations/DNA/substitution.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/substitution/) (mosaic `/` and chimeric `//` markers are presence-only; the spec carries no allele fraction in the HGVS string) | Enforced — parser-level rejection at the trailing-content check in `src/hgvs/parser/variant.rs::parse_variant`. Closes the diagnostic gap left by [#133](https://github.com/fulcrumgenomics/ferro-hgvs/issues/133) work-item 4 (allele-fraction soft-warning code follow-up, [#278](https://github.com/fulcrumgenomics/ferro-hgvs/issues/278)). Cannot be auto-corrected: the fraction belongs in accompanying metadata (VCF `FORMAT/AF`, ClinVar heteroplasmy field, etc.), not in the HGVS string |
| W3018 | ClinVarProseMultiAllelic — ClinVar prose shorthand `m.<pos><ref>><alt>/<alt2>` (bare-base RHS) rejected with a targeted SVA diagnostic pointing at the three spec-supported alternatives (compound brackets, dual fully-qualified slash, compact mosaic) | [recommendations/DNA/substitution.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/substitution/) (mosaic `/` requires fully qualified RHS or the compact `<pos>=/<alt-edit>` form; the bare-base prose shape is not a spec form) | Enforced — parser-level rejection in `src/hgvs/parser/variant.rs::parse_phase_allele` when the bare-base prose RHS is detected after the spec-compact fallback fails. Closes the diagnostic gap left by [#133](https://github.com/fulcrumgenomics/ferro-hgvs/issues/133) work-item 4 (mosaic prose follow-up, [#278](https://github.com/fulcrumgenomics/ferro-hgvs/issues/278)) |
| W3019 | NonSpecMosaicForm — nested `/` + `//` at the same level (chimeric-of-mosaic / mosaic-of-chimeric) AND the bracketed `[a/b]` / `[a//b]` mosaic group; one diagnostic family for two detectors | [recommendations/DNA/substitution.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/substitution/) defines only compact mosaic `<pos>=/<edit>` and chimeric `<pos>=//<edit>`; HGVS v21 does **not** define nesting or `[a/b]` | Enforced (always-reject — no canonical alternative to auto-correct to). Two parser-level detectors in `src/hgvs/parser/variant.rs`: (A) `has_slash_inside_brackets` at the top of `parse_variant` catches `[a/b]` / `[a//b]`; (B) the existing chunk-has-inner-slash guard inside `parse_phase_allele` upgrades to the same family. Both call `non_spec_mosaic_form_error`, which inlines the spec alternatives (`[a;b]` cis; `acc:c.X/acc:c.Y` mosaic; `acc:c.X//acc:c.Y` chimeric; compact `=/<edit>` / `=//<edit>`) in `msg` so plain `to_string()` shows the migration path. Pinned by `tests/issue_281_mosaic_chimeric_nesting.rs` (closes #281, follow-up to #217 / #216) |
| W3020 | RnaThymineCanonicalized — `t`/`T` inside an `r.` (RNA) description rewritten to `u` | [recommendations/RNA/substitution.md](https://hgvs-nomenclature.org/stable/recommendations/RNA/substitution/) and adjacent RNA edit pages (RNA alphabet is `a/c/g/u`; `t` is non-canonical) | Enforced — pre-parse rewrite in lenient/silent modes via `correct_rna_thymine` (preprocessor Phase 16a); strict rejects with `E1006 UnexpectedChar` + hint. Closes the input-side loop on issue #232's deferred policy decision (closes #282); PR #293 / Display canonicalization handles the output side independently for raw `parse_hgvs` callers |
| W3021 | ProteinBracketedAaInsertion — `p.…ins[Ala;Pro]` rejected; canonical form is `p.…insAlaPro` | [recommendations/protein/variant/insertion.md](https://hgvs-nomenclature.org/stable/recommendations/protein/variant/insertion/) (protein insertions concatenate 3-letter codes with no separator; brackets are reserved for variant-level alleles) | Enforced — pre-parse reject; cannot be auto-rewritten because mixing 1-letter and 3-letter codes inside `[...]` is ambiguous. Strict / lenient / silent all reject with `E1004 InvalidEdit` + hint naming the canonical `insAlaPro` shape; the bare `parse_hgvs` path also rejects via a precheck inside `parse_variant`, so the diagnostic is identical whether the caller uses the preprocessor or not. Closes [#290](https://github.com/fulcrumgenomics/ferro-hgvs/issues/290) (follow-up to [#248](https://github.com/fulcrumgenomics/ferro-hgvs/pull/248)) |
| W4001 | SwappedPositions — `c.200_100del` corrected to `c.100_200del` | [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) line 65 (range semantics: start ≤ end implied by the underscore range syntax and reinforced throughout DNA/numbering examples) | Enforced — wired as InputPreprocessor Phase 15a via `correct_swapped_positions`. Lenient mode emits W4001 and rewrites the swapped pair; silent mode rewrites without emitting a warning. The endpoint tokenizer captures `(main_axis, offset)` sort keys per endpoint, so intronic offsets (`c.100+5_99+3del`), 3'UTR markers (`c.*5_*1del`), 5'UTR markers (`c.-3_-5del`), and cross-region swaps (`c.5_-3del`, `c.*1_100del`) are all detected and rewritten with offsets/markers preserved (#264 wired; #265 extended) |
| W4003 | SinglePositionRange — `c.123_123del` / `c.123_123dup` / `c.100_100inv` collapsed to the single-position form | [recommendations/DNA/deletion.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/deletion/): "position(s)\_deleted should contain two different positions, e.g. 123\_126 not 123\_123"; analogous wording in [recommendations/DNA/duplication.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/duplication/) and [recommendations/DNA/inversion.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/inversion/) | Enforced — lenient/silent collapse the redundant range and warn once; strict rejects. Wired as InputPreprocessor Phase 10 via `correct_single_position_range`, covering `del` / `dup` / `inv` (#127, #81 L2 SVA-008/009/014) |
| W4004 | PositionPastEnd — `c.<N>G>C` / `c.<N>dup` / `c.<a>_<b>del` where the position exceeds the resolved transcript's CDS-end or transcript-end (e.g. `c.946G>C` against an `NM_` whose CDS is 945 bases) | [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) (position numbering must reference an existing base in the transcript) | Enforced — emitted as `NormalizationWarning::PositionPastEnd` from `src/normalize/mod.rs::check_cds_pos_past_end`, called from `normalize_cds` before the `needs_normalization` short-circuit so substitutions are also checked. Strict mode converts the warning into `FerroError::InvalidCoordinates` in the outer `Normalizer::normalize` wrapper; lenient mode emits the warning and short-circuits to the canonical variant; silent mode skips the check entirely. Covers `c.<N>` against CDS length and `c.*<N>` against post-CDS transcript suffix. Intronic offsets, 5'UTR (`c.-N`), and `n.` variants are out of scope for the initial implementation (#336) |
| W5001 | RefSeqMismatch — stated reference base mismatches actual reference (warning sibling of `E3002`) | [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) line 15: "the reference sequence used must contain the residue(s) described to be changed" | Partial — emitted as `NormalizationWarning::RefSeqMismatch` from `src/normalize/mod.rs` when `validate_reference` detects a mismatch during normalization. Coverage is partial because emission requires loaded reference data; without reference data (the common CLI case) parsing succeeds with no warning. With reference data the sibling error `E3002` is raised in strict mode instead |
| W5002 | OverlapConflictingEdits — two or more cis-allele edits share identical reference bounds (e.g. `g.[100G>A;100A>C]`) | [recommendations/DNA/alleles.md](https://hgvs-nomenclature.org/stable/recommendations/DNA/alleles/) (cis-allele descriptions); [recommendations/general.md](https://hgvs-nomenclature.org/stable/recommendations/general/) line 15: "the reference sequence used must contain the residue(s) described to be changed" | Enforced — detected in `src/normalize/overlap.rs` by `detect_overlap_conflicts` and wired into `normalize_allele` after per-variant normalization. Variant output is preserved; the warning is advisory |
| W5003 | VariantExceedsReference — provider returned fewer bytes than the HGVS interval span; variant references positions outside the reference (e.g. `NG_032871.1:g.32476_53457delins...` when the reference is shorter than position 53457) | [background/refseq.md](https://hgvs-nomenclature.org/stable/background/refseq/#reference-sequences) §43: "the entirety of the variant sequence must be encompassed by the selected reference sequence" | Enforced — emitted as `NormalizationWarning::CanonicalSplitSkipped` from `apply_canonical_split` in `src/normalize/mod.rs` when `ref_bytes.len() != expected_span`. Strict mode promotes to `FerroError::VariantExceedsReference` via `Normalizer::normalize`; lenient warns + preserves; silent preserves without warning. Matches biocommons hgvs which raises `HGVSInvalidVariantError` (#355) |

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
