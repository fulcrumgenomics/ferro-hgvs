# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- *(loader)* unified `load_annotations` entry point with `LoaderConfig`/`LoaderReport`; auto-detects GFF3 vs GTF by extension and content (#191, #194, #195)
- *(loader)* exon-derivation ladder closes single-exon GFF3 without `exon` records (#183 fix) (#191)
- *(loader)* phase-aware CDS bounds with `start_codon` / `stop_codon` precedence; GFF3 input now extends `cds_end` to include the stop codon, matching ferro's downstream convention (fixes a latent off-by-one-codon bug in GFF3 protein conversion) (#194)
- *(loader)* UTR-merge ladder step: synthesizes exons from UTR + CDS when no explicit `exon` records (#194)
- *(loader)* `gene_symbol`, `mane_status`, `refseq_match`, `ensembl_match` extracted from attributes (#194)
- *(loader)* optional FASTA-aware validation: CDS length mod 3 and start codon (ATG/CTG/GTG/TTG) checks when `--fasta` is supplied
- *(error_handling)* 13 loader diagnostic codes (E-LOAD-*, W-LOAD-*) registered for `ferro explain` (#195)
- *(cli)* `ferro convert-gff --strict / --silent / --no-validate-fasta / --diagnostics-json` flags (#195)
- *(strand)* `Strand::Unknown` variant for GFF3 `.` and `?` strand values; transcripts with unknown strand are dropped at load with `E-LOAD-103` (#191)

### Changed

- *(strand)* `Strand`, `LoaderConfig`, and `LoaderReport` are marked `#[non_exhaustive]` for forward compatibility (#191)
- *(loader)* `load_gff3` and `load_gtf` have been removed; use `load_annotations` instead

## [0.5.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.4.1...v0.5.0) - 2026-05-12

### Added

- *(parser)* spec compact-mosaic form across sub/del/dup ([#133](https://github.com/fulcrumgenomics/ferro-hgvs/pull/133)) ([#153](https://github.com/fulcrumgenomics/ferro-hgvs/pull/153))
- *(normalize)* canonicalize con to delins per SVD-WG009 (#81 H1) ([#142](https://github.com/fulcrumgenomics/ferro-hgvs/pull/142))
- *(rna)* support r.spl, r.spl?, r.(spl), r.(spl?) splicing markers (#81 E2) ([#134](https://github.com/fulcrumgenomics/ferro-hgvs/pull/134))

### Fixed

- *(spdi)* emit inversion and repeat as SPDI delins ([#159](https://github.com/fulcrumgenomics/ferro-hgvs/pull/159))
- *(error-handling)* error codes for spec-mandated input rejections ([#115](https://github.com/fulcrumgenomics/ferro-hgvs/pull/115)) ([#152](https://github.com/fulcrumgenomics/ferro-hgvs/pull/152))
- *(spdi)* reference-aware HGVS→SPDI for del/dup/delins ([#158](https://github.com/fulcrumgenomics/ferro-hgvs/pull/158))
- *(parser)* preserve explicit deleted sequence in delins round-trip ([#154](https://github.com/fulcrumgenomics/ferro-hgvs/pull/154))
- *(error-handling)* soft-warn non-canonical input forms at parse time ([#127](https://github.com/fulcrumgenomics/ferro-hgvs/pull/127)) ([#151](https://github.com/fulcrumgenomics/ferro-hgvs/pull/151))
- *(ci)* drop redundant `#[allow(dead_code)]` on `mod common;` ([#178](https://github.com/fulcrumgenomics/ferro-hgvs/pull/178))
- *(spdi)* accept c./n./r./m. variants in HGVS->SPDI conversion ([#157](https://github.com/fulcrumgenomics/ferro-hgvs/pull/157))
- *(normalize)* apply 3' rule across phase-mismatched cyclic rotations for ins ([#132](https://github.com/fulcrumgenomics/ferro-hgvs/pull/132)) ([#155](https://github.com/fulcrumgenomics/ferro-hgvs/pull/155))
- *(error-handling)* surface deprecated stop-codon and frameshift forms as soft-warns ([#125](https://github.com/fulcrumgenomics/ferro-hgvs/pull/125)) ([#150](https://github.com/fulcrumgenomics/ferro-hgvs/pull/150))
- *(error-handling)* wire W1001/W1002/W3001 soft-validation warnings ([#124](https://github.com/fulcrumgenomics/ferro-hgvs/pull/124)) ([#149](https://github.com/fulcrumgenomics/ferro-hgvs/pull/149))
- *(spdi)* recover dup form on SPDI→HGVS for duplicated insertions ([#156](https://github.com/fulcrumgenomics/ferro-hgvs/pull/156))
- *(parser)* extend unknown-phase (;) support to g./n./m./o./p. coord systems ([#123](https://github.com/fulcrumgenomics/ferro-hgvs/pull/123)) ([#148](https://github.com/fulcrumgenomics/ferro-hgvs/pull/148))
- *(normalize)* detect cis-allele edits with coincident bounds (#81 A8) ([#147](https://github.com/fulcrumgenomics/ferro-hgvs/pull/147))
- *(parser)* soft-warn embedded whitespace and zero-width chars ([#128](https://github.com/fulcrumgenomics/ferro-hgvs/pull/128)) ([#145](https://github.com/fulcrumgenomics/ferro-hgvs/pull/145))
- *(normalize)* preserve r.*N UTR flag and translate via cds_end (closes #163) ([#164](https://github.com/fulcrumgenomics/ferro-hgvs/pull/164))
- *(parser)* correct LRG accession variant-type inference and compound-ref handling ([#122](https://github.com/fulcrumgenomics/ferro-hgvs/pull/122)) ([#141](https://github.com/fulcrumgenomics/ferro-hgvs/pull/141))
- *(normalize)* recognize revcomp inv sub-spans within delins ([#166](https://github.com/fulcrumgenomics/ferro-hgvs/pull/166))
- *(normalize)* apply 3'-rule to merged cis-allele deletions (closes #161) ([#162](https://github.com/fulcrumgenomics/ferro-hgvs/pull/162))
- *(normalize)* rewrite revcomp delins as inversion (#81 A2) ([#109](https://github.com/fulcrumgenomics/ferro-hgvs/pull/109))
- *(normalize)* enforce HGVS c. codon-frame exception for repeat notation (#81 B1) ([#110](https://github.com/fulcrumgenomics/ferro-hgvs/pull/110))
- *(normalize)* rewrite empty-insert delins as del (#81 A3) ([#113](https://github.com/fulcrumgenomics/ferro-hgvs/pull/113))
- *(normalize)* degenerate substitution (ref==alt) -> identity (#81 A4) ([#111](https://github.com/fulcrumgenomics/ferro-hgvs/pull/111))

### Other

- *(parser)* triage failure expectations (#174 phase 2) ([#176](https://github.com/fulcrumgenomics/ferro-hgvs/pull/176))
- *(parser)* per-input failure-expectations framework (#174 phase 1) ([#175](https://github.com/fulcrumgenomics/ferro-hgvs/pull/175))
- drop ferro_version from HGVS spec fixture ([#177](https://github.com/fulcrumgenomics/ferro-hgvs/pull/177))
- *(test)* cut Test job runtime ~398s → ~20s ([#173](https://github.com/fulcrumgenomics/ferro-hgvs/pull/173))
- *(test)* correct misleading comment on trans-allele expanded-form test ([#167](https://github.com/fulcrumgenomics/ferro-hgvs/pull/167))
- *(normalize)* tag r. positive bases as Region::Rna (closes #168) ([#169](https://github.com/fulcrumgenomics/ferro-hgvs/pull/169))
- *(allele)* pin trans-phase round-trip across all coord systems and merge-barrier (#81 C1) ([#146](https://github.com/fulcrumgenomics/ferro-hgvs/pull/146))
- *(compound)* pin cross-reference / cross-coord compound round-trip and merge-barrier (#81 H2) ([#143](https://github.com/fulcrumgenomics/ferro-hgvs/pull/143))
- *(error-handling)* audit error codes against HGVS spec sections (#81 L1) ([#137](https://github.com/fulcrumgenomics/ferro-hgvs/pull/137))
- *(parser)* pin gene-selector round-trip end-to-end with Display preservation (#81 I3) ([#135](https://github.com/fulcrumgenomics/ferro-hgvs/pull/135))
- *(mito)* audit heteroplasmy notation; tracking #133 (#81 F2) ([#139](https://github.com/fulcrumgenomics/ferro-hgvs/pull/139))
- *(normalize)* pin RNA path + edge cases for A9 substitution-after-trim ([#114](https://github.com/fulcrumgenomics/ferro-hgvs/pull/114))
- *(protein)* pin p.? unknown-effect round-trip across allele forms and edge cases (#81 D7) ([#136](https://github.com/fulcrumgenomics/ferro-hgvs/pull/136))
- *(protein)* pin p.0 no-product round-trip and adjacent guards (#81 D6) ([#130](https://github.com/fulcrumgenomics/ferro-hgvs/pull/130))
- *(protein)* pin silent `=` round-trip across allele forms and edge cases (#81 D5) ([#131](https://github.com/fulcrumgenomics/ferro-hgvs/pull/131))
- *(mito)* audit m. coord-system parse + wraparound behavior; tracking #129 (#81 F1) ([#138](https://github.com/fulcrumgenomics/ferro-hgvs/pull/138))
- pin normalize() against HGVS v21.0 spec fixture (closes #84) ([#105](https://github.com/fulcrumgenomics/ferro-hgvs/pull/105))

### Added

- *(normalize)* detect cis-allele edits with coincident reference bounds (e.g. `g.[100G>A;100A>C]`) and emit an advisory `OVERLAP_CONFLICTING_EDITS` (`W5002`) warning. Variant output is preserved unchanged. Addresses [#81](https://github.com/fulcrumgenomics/ferro-hgvs/issues/81) item A8.
- *(test)* `expected_warnings` field on the HGVS v21.0 spec-fixture row schema, pinning the warning set ferro emits per row.

### Fixed

- *(normalize)* apply 3' rule across phase-mismatched cyclic rotations for single-copy ins. When the inserted alt is a non-zero cyclic rotation of an adjacent reference repeat unit (e.g. `g.X_(X+1)insTG` against a `GT[3]` tract), shuffle's first-base check (`alt[0] == ref[ins_point]`) failed and the variant never moved. The new helper `insertion_to_duplication` mirrors `insertion_to_repeat`'s rotation iteration for the 1-copy case, so the variant now canonicalizes to a `dup` at the most-3' rotation-aligned position. Closes [#132](https://github.com/fulcrumgenomics/ferro-hgvs/issues/132).
- *(normalize)* recognize a reverse-complement sub-span within a `delins` (synthesized by cis-allele merge OR user-typed) and emit the spec-canonical `inv`, splitting the surrounding span into `[…; inv; …]`. The HGVS edit-priority rule places `inv` above `delins` in the priority order (`general.md:56`) and defines `delins` as the residual when no higher-priority form applies (`delins.md`). For example, `g.[1150T>G;1151C>A;1152C>G]` (over `TCC`) now normalizes to `g.[1150_1151inv;1152C>G]` instead of `g.1150_1152delinsGAG`; `g.[1092G>C;1093G>C]` (over `GG`) normalizes to `g.1092_1093inv` instead of `g.1092_1093delinsCC`. The same rule fires for a user-typed `g.1150_1152delinsGAG`, since the canonical form depends on `(ref, position, alt)` and not on input shape. Applies across all NA coord systems: `g.`, `m.`, `c.` (CDS-proper positions), `n.`, `r.` (T/U-equivalent comparison so `r.` alts with `U` align with transcript ref bytes that contain `T`). Sub-only decomposition (rewriting a non-inv multi-sub `delinsXY` to `[X>...; Y>...]`) is intentionally left out of scope and is a separate spec interpretation question. The codon-frame `c.` merge from issue [#79](https://github.com/fulcrumgenomics/ferro-hgvs/issues/79) is preserved automatically — the synthesized middle base in a codon-frame-merged delins makes a length-2 inv across the middle mathematically impossible. Closes issue [#160](https://github.com/fulcrumgenomics/ferro-hgvs/issues/160).
- *(normalize)* rewrite degenerate substitutions (ref == alt, e.g. `c.100A>A`) to identity (`=`) per HGVS v21 spec, which marks `c.X>X` as "not allowed" (`docs/recommendations/DNA/other.md`). The rule is purely syntactic on the edit's stated bases, so it fires in both the full-normalization path and the no-reference canonicalization path — `c.123C>C` rewrites to `c.123=` regardless of provider availability. ([#81](https://github.com/fulcrumgenomics/ferro-hgvs/issues/81) A4)
- *(normalize)* rewrite a `delins` whose inserted sequence is empty as a deletion, per the HGVS spec requirement that an empty insert is semantically a deletion and must be rendered as `del`. The rewritten deletion is then 3'-shifted under the standard del rule. Issue [#81](https://github.com/fulcrumgenomics/ferro-hgvs/issues/81) item A3.
- *(normalize)* rewrite a `delins` whose inserted sequence is the reverse complement of the deleted reference as an inversion, per the HGVS spec definition of `inv`. The complementary-outer-bases shortening rule applies to the result so that a `delins`-encoded inversion produces the same canonical output as a directly-encoded `inv`. Issue [#81](https://github.com/fulcrumgenomics/ferro-hgvs/issues/81) item A2.
- *(normalize)* canonicalize `delins` to the minimal HGVS form by trimming any shared prefix/suffix between the inserted sequence and the deleted reference, then reclassifying the residual edit as substitution / deletion / insertion / inversion / smaller `delins` per the sub > del > inv > dup > ins priority. For example, `c.1_4delinsAAGC` against ref `ATGC` collapses to `c.2T>A`; `c.1_4delinsAC` against ref `ATGC` collapses to `c.2_3del`; `c.1_3delinsATCG` against ref `ATG` collapses to `c.2_3insC`. Extension to issue [#81](https://github.com/fulcrumgenomics/ferro-hgvs/issues/81) item A2.
- *(normalize)* apply inversion shortening to direct `inv` inputs. `NaEdit::Inversion` was missing from `needs_normalization()`, so direct inversions bypassed `normalize_na_edit` entirely and the `shorten_inversion()` / identity-collapse logic was never exercised. Direct `inv` variants now also emit minimal notation after shortening (no stale explicit `sequence`/`length` from the input).
- *(normalize)* HGVS spec compliance ([#81](https://github.com/fulcrumgenomics/ferro-hgvs/issues/81) B1): repeat-notation rewrites in `c.` (coding DNA) context now enforce the spec's codon-frame exception — repeat notation requires `unit_len % 3 == 0`. Previously `insertion_to_repeat`, `deletion_to_repeat`, `duplication_to_repeat`, and `normalize_repeat` could emit `c.X_YA[N]`, `c.X_YAT[N]`, etc. for non-codon-aligned units, violating the spec (`docs/recommendations/DNA/repeated.md`). Variants in `c.` with non-codon-aligned units now retain the spec-prescribed alternative form: `dup` for 1 added copy, `ins<literal>` for ≥2 added copies, and plain `del` for contractions of ≥2 unit copies.

### Changed (public API)

- `insertion_to_repeat`, `duplication_to_repeat`, and `normalize_repeat` in `ferro_hgvs::normalize::rules` gain an `is_coding: bool` parameter to drive the codon-frame gate. (`deletion_to_repeat` is `pub(crate)` and gains the same parameter as an internal change.)
- `RepeatNormResult::Insertion { start, end, sequence }` and `DupToRepeatResult::GatedInsertion { start, end, sequence }` variants added so the rule layer can route gated rewrites to the spec-canonical literal `ins` form.
- *(normalize)* `NormalizationWarning` is now a sum type (`RefSeqMismatch` / `OverlapConflict`). Each warning code carries only the fields relevant to it. Read sites migrate from `.code` field access to `.code()` method and pattern-matching on the variant.

### Changed

- Internal: the four `delins_is_*` boolean helpers in `normalize::rules` are unified into one `canonicalize_delins()` function returning a `DelinsCanonical` enum, expressing HGVS edit-priority (sub > del > inv > dup > ins) in a single decision tree. The unreachable second delins arm in `normalize_na_edit` has been removed.

## [0.4.1](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.4.0...v0.4.1) - 2026-05-04

### Added

- *(python)* load reference data via Normalizer.from_manifest and extended from_json ([#86](https://github.com/fulcrumgenomics/ferro-hgvs/pull/86))
- *(normalize)* merge consecutive edits in cis alleles per HGVS spec ([#80](https://github.com/fulcrumgenomics/ferro-hgvs/pull/80))

### Fixed

- *(normalize)* codon-frame exception for c. SNV pairs separated by one nucleotide ([#79](https://github.com/fulcrumgenomics/ferro-hgvs/pull/79)) ([#104](https://github.com/fulcrumgenomics/ferro-hgvs/pull/104))
- *(normalize)* 5'UTR CDS↔tx off-by-one collapsed UTR del to c.? ([#97](https://github.com/fulcrumgenomics/ferro-hgvs/pull/97)) ([#102](https://github.com/fulcrumgenomics/ferro-hgvs/pull/102))
- *(normalize)* merge same-region UTR adjacency in cis alleles ([#89](https://github.com/fulcrumgenomics/ferro-hgvs/pull/89)) ([#103](https://github.com/fulcrumgenomics/ferro-hgvs/pull/103))
- *(normalize)* minus-strand intronic ref-base orientation ([#98](https://github.com/fulcrumgenomics/ferro-hgvs/pull/98)) ([#100](https://github.com/fulcrumgenomics/ferro-hgvs/pull/100))
- rewrite delins as identity when insert matches reference ([#78](https://github.com/fulcrumgenomics/ferro-hgvs/pull/78))
- rewrite single-base delins as substitution per HGVS priority ([#77](https://github.com/fulcrumgenomics/ferro-hgvs/pull/77))
- Emit single-variant alleles in bare spec form ([#76](https://github.com/fulcrumgenomics/ferro-hgvs/pull/76))

### Other

- unpin CI Rust toolchain and fix 1.95 clippy lints ([#108](https://github.com/fulcrumgenomics/ferro-hgvs/pull/108))
- dup 3'-shift coverage matrix (#81 A6) ([#107](https://github.com/fulcrumgenomics/ferro-hgvs/pull/107))
- del 3'-shift coverage matrix + tandem-repeat del canonical-form fix (#81 A5, B2) ([#106](https://github.com/fulcrumgenomics/ferro-hgvs/pull/106))
- tighten coverage_gap_tests assertions; restore intronic-ins coverage ([#94](https://github.com/fulcrumgenomics/ferro-hgvs/pull/94)) ([#99](https://github.com/fulcrumgenomics/ferro-hgvs/pull/99))
- *(coverage)* restore intronic insertion tests dropped in #93 ([#95](https://github.com/fulcrumgenomics/ferro-hgvs/pull/95))
- ins 3'-shift coverage matrix + tandem-repeat ins canonical-form fix (#81 A1, A7) ([#93](https://github.com/fulcrumgenomics/ferro-hgvs/pull/93))

### Added

- *(test)* HGVS v21.0 spec normalization regression fixture pinning ferro's current `normalize()` output for every variant string in the spec (823 rows). Each row carries a status (`preserved` / `diverges` / `correctly-rejected` / `false-acceptance` / `parse-error` / `needs-reference`) derived from `(parse_ok, normalize_ok, spec_expected, current == spec_expected)`. `spec_expected: null` is the sentinel for "spec rejects this," set automatically for inputs the spec marks via `<code class="invalid">…</code>` (35 occurrences in v21.0) and overridable by hand. Bare illustrative fragments like `c.1083A>C` get a default accession prepended per coord system so they parse — the prefixed form is recorded in a separate `input_prefixed` field. The fixture is generated by `examples/generate_spec_fixture.rs` from the vendored spec at `assets/hgvs-nomenclature/` (git submodule pinned to tag `21.0.0`). Hand overrides live in `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json` and accept `status`, `spec_expected` (string / null / absent), `input_prefixed`, `requires_reference`, and `todo`. CI runs the generator in `--check` mode to enforce byte-identical regeneration. Closes [#84](https://github.com/fulcrumgenomics/ferro-hgvs/issues/84); companion to [#83](https://github.com/fulcrumgenomics/ferro-hgvs/issues/83).
- 3'-shift coverage matrix for insertion variants
  (`tests/ins_shift_matrix.rs`): 84 rstest cases across all 7 nucleotide
  coord-system / strand combinations × 8 shuffle scenarios, with a
  shared `SyntheticBuilder` fixture helper at `tests/common/synthetic.rs`
  (reusable by future del / dup / repeat-notation matrices). Issue #81
  item A7.
- Tightened `tests/coverage_gap_tests.rs` from canary-shaped
  (`contains(...)`, `is_ok() || is_err()`) assertions to strict
  `assert_eq!` against the exact normalizer output. Builds on PR #95
  (intronic-insertion restoration) by extending the same tightening
  pass to the remaining gap-test buckets, and adds cross-references
  and procedural detail to the restored intronic-insertion comments.
  Several locked outputs are suspected-buggy and carry `FIXME(#NN)`
  comments pointing to follow-up tracking issues — most converge on
  #98 (minus-strand intronic reference-base orientation, identified
  as a likely common root cause behind A1 / A5 / A6 misfires on
  minus-strand intronic positions and the PR #78 identity-rewrite
  missing on minus strand), with #96 (wrong-strand repeat unit
  emission on minus-strand multi-base dup) and #97 (`c.?del` collapse
  on minus-strand 5'UTR del) as separate sub-symptoms. True
  boundary-spanning panic-canary tests are intentionally left as
  `is_ok() || is_err()` per the issue's scope. Issue #94.
- *(normalize)* Merge consecutive sub-variants in cis alleles into a single delins per HGVS spec. `g.[1000G>A;1001A>C]` now normalizes to `g.1000_1001delinsAC`; `g.[1000del;1001del]` to `g.1000_1001del`. Covers `g./c./n./r./m.` coordinate systems, sub/del/delins/ins edit combinations, chains, and same-boundary insertion pairs. Non-adjacent variants, intronic/UTR boundaries, uncertain edits, and non-`Literal` insertion payloads are barriers and pass through unchanged. The codon-frame exception (one-nt gap within a codon) is tracked separately in [#79](https://github.com/fulcrumgenomics/ferro-hgvs/issues/79). ([#80](https://github.com/fulcrumgenomics/ferro-hgvs/pull/80))

### Fixed

- *(normalize)* Cis-allele consecutive-edit merging now collapses
  adjacent sub-variants *within* a UTR / upstream / downstream region,
  not only within the CDS / transcript body. The PR #80 implementation
  rejected every `is_5utr() / is_3utr() / is_upstream() / is_downstream()`
  position outright, so inputs like `c.[-2A>G;-1C>T]` (both 5'UTR) and
  `c.[*1A>G;*2C>T]` (both 3'UTR) round-tripped unchanged even though
  HGVS allows ranges within those regions (`c.-2_-1`, `c.*1_*2`). The
  fix replaces the `Option<u64>` position keys with a `Region`-tagged
  `(Region, i64)` axis (covering CDS, 5'UTR, 3'UTR, transcript-body,
  upstream, downstream, and genomic / mitochondrial); merge eligibility
  becomes "same region + integer adjacency on the region's axis".
  `build_cds_merged` / `build_tx_merged` / `build_rna_merged` consume
  the region tag to reconstruct the right `CdsPos` / `TxPos` / `RnaPos`
  shape (negative base for 5'UTR / upstream, `utr3` / `downstream` flag
  for 3'UTR / downstream). Cross-region pairs (`c.[-1A>G;1A>T]` 5'UTR↔CDS,
  `c.[40C>T;*1A>G]` CDS↔3'UTR, …) still correctly do not merge.
  Issue #89.
- *(normalize)* CDS ↔ transcript coordinate mappings now respect the
  HGVS no-c.0 numbering rule for 5'UTR positions. The forward mapping
  (`Normalizer::cds_to_tx_pos`, `convert::coding::cds_to_transcript_pos`)
  previously computed `tx = cds_start + base - 1` for negative `base`,
  which double-counted the gap between c.-1 and c.1 and emitted a tx
  position one base 5' of the true location. The inverse mapping
  (`Normalizer::tx_to_cds_pos`) had the mirror bug: tx positions one
  base before `cds_start` mapped to `base = 0`, which `CdsPos::Display`
  renders as `c.?` (`CDS_BASE_UNKNOWN`). The most visible symptom was a
  5'UTR single-base deletion on a minus-strand transcript collapsing
  to `c.?del` instead of resolving to a real position (e.g. `c.-1del`
  after 3'-shifting within a UTR homopolymer). Forward and inverse are
  now `tx = cds_start + base` and `base = tx - cds_start` for negative
  / pre-cds_start positions, matching the spec's "c.-1 is one base 5'
  of c.1" rule and the existing exon-aware mapper at
  `convert::mapper::cds_to_tx` / `tx_to_cds`. Issue #97.
- *(normalize)* Cis-allele consecutive-edit merging now also collapses
  the codon-frame exception case: two `c.` exonic SNVs in the CDS
  proper, separated by exactly one nucleotide, that fall within the
  same codon merge into a single delins with the unchanged middle
  reference base preserved verbatim — per HGVS spec
  (`c.[145C>T;147C>G]` → `c.145_147delinsTGG`, where the middle base
  is the reference at `c.146`). Eligibility is narrow: same accession,
  both endpoints in `Region::Cds`, gap-of-one on the axis, both prev
  and next are single-base SUB anchors, and `(prev-1)/3 == (next-1)/3`
  (same codon). The unchanged middle base is fetched via the
  `ReferenceProvider` threaded into `merge_consecutive_edits`; if no
  transcript is registered or the position is out of range, the merge
  is silently declined and the variants pass through unchanged.
  Codon-frame–merged delins continue to participate in the
  strictly-consecutive walk, so a third SNV one base 3' of the pair
  (`c.[10A>G;12A>C;13A>T]`) still folds into the running delins.
  Cross-codon (`c.[3G>T;5A>C]`), gap-of-two (`c.[10A>G;13A>C]`), and
  non-CDS pairs (`g.`, UTR, `n.`, `r.`) all correctly do not merge.
  Issue #79.
- *(normalize)* Minus-strand intronic normalization now reads the
  reference window in transcript-view orientation. `normalize_intronic_cds`
  and `normalize_intronic_tx` previously passed the genomic-strand bytes
  fetched from `get_genomic_sequence` directly into `normalize_na_edit`
  alongside the variant's transcript-view edit alt; on minus-strand
  transcripts the two were mis-oriented, defeating every rule that
  compared the alt against the local reference window. The fix
  reverse-complements the genomic window and flips the relative
  positions / shuffle boundaries on minus strand before normalization,
  then maps the resulting positions back to the genomic frame for the
  CDS / tx coordinate conversion. As a single root-cause fix this
  resurfaces #81 A1 / A5 / A6 canonicalization, the PR #78 delins-as-
  identity rewrite, and the transcript-view repeat-unit letter on
  minus-strand intronic positions — all of which had been observed to
  misfire in #94's locking pass. Issue #98.
- Insertions that add ≥2 copies of a multi-base tandem repeat unit now
  emit repeat notation (`unit[N+k]`) instead of a duplication of the
  inserted sequence, per HGVS spec ("when more than one additional copies
  are inserted directly 3' of the original copy, the change is indicated
  using the format for Repeated sequences"). Single-unit additions remain
  `dup`. Issue #81 item A7.
- `MockProvider::get_sequence` now falls through to genomic contig
  lookup when the id is not a transcript, matching `FastaProvider`'s
  behavior. This unblocks 3'-shift normalization for genomic test
  fixtures that register only `add_genomic_sequence`. Issue #81 item A7.

## [0.4.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.3.0...v0.4.0) - 2026-04-30

### Added

- *(python)* add poethepoet task runner for Python dev workflow ([#53](https://github.com/fulcrumgenomics/ferro-hgvs/pull/53))
- Add variant accessor properties to Python bindings ([#49](https://github.com/fulcrumgenomics/ferro-hgvs/pull/49))

### Fixed

- accept gene selectors on non-RefSeq accessions ([#70](https://github.com/fulcrumgenomics/ferro-hgvs/pull/70))
- Use HGVS spec compact form for allele Display ([#48](https://github.com/fulcrumgenomics/ferro-hgvs/pull/48))

### Other

- *(readme)* vendor Fulcrum logo and use absolute URLs ([#71](https://github.com/fulcrumgenomics/ferro-hgvs/pull/71))
- bump pyo3 0.23 → 0.28 and add Python 3.14 wheels ([#57](https://github.com/fulcrumgenomics/ferro-hgvs/pull/57))
- *(python)* use dependency-groups, stricter mypy, and --locked in CI ([#52](https://github.com/fulcrumgenomics/ferro-hgvs/pull/52))
- publish Python wheels to PyPI via Trusted Publishing ([#58](https://github.com/fulcrumgenomics/ferro-hgvs/pull/58))
- *(python)* drop Python 3.8 and 3.9 support ([#55](https://github.com/fulcrumgenomics/ferro-hgvs/pull/55))
- *(python)* add uv lockfile for reproducible dev environment ([#50](https://github.com/fulcrumgenomics/ferro-hgvs/pull/50))
- build Python wheels and attach to GitHub Releases ([#39](https://github.com/fulcrumgenomics/ferro-hgvs/pull/39))
- *(prepare)* Modularize ReferenceManifest ([#44](https://github.com/fulcrumgenomics/ferro-hgvs/pull/44))
- switch reqwest from native-tls to rustls-tls ([#38](https://github.com/fulcrumgenomics/ferro-hgvs/pull/38))

### Changed

- Allele `Display` now emits HGVS spec-correct compact form (`ACC:c.[edit1;edit2]`) when sub-variants share an accession and coordinate type, instead of the expanded form (`[ACC:c.edit1;ACC:c.edit2]`). Mixed-accession alleles and alleles containing the per-variant unknown form (`c.?`, `r.?`, etc.) still emit the expanded form. Downstream consumers parsing the previous expanded output (including Python `str(variant)`) will see the new format. ([#48](https://github.com/fulcrumgenomics/ferro-hgvs/pull/48))

### Fixed

- Prevent panic when `NullAllele`/`UnknownAllele` are used as sub-variants in an allele ([#48](https://github.com/fulcrumgenomics/ferro-hgvs/pull/48))

## [0.3.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.2.0...v0.3.0) - 2026-03-30

### Added

- Bincode serialization for fast cdot loading ([#26](https://github.com/fulcrumgenomics/ferro-hgvs/pull/26))
- Compound reference syntax `NC_*(NM_*):c.…` support ([#21](https://github.com/fulcrumgenomics/ferro-hgvs/pull/21))
- GRCh37 cdot transcript metadata download in `ferro prepare` ([#24](https://github.com/fulcrumgenomics/ferro-hgvs/pull/24))

### Fixed

- Genomic-to-coding coordinate conversion ([#25](https://github.com/fulcrumgenomics/ferro-hgvs/pull/25))
- Coding-order positions after genomic-space normalization ([#20](https://github.com/fulcrumgenomics/ferro-hgvs/pull/20))
- Mutalyzer configuration across web service and benchmarks ([#19](https://github.com/fulcrumgenomics/ferro-hgvs/pull/19))
- Serde defaults for config structs and deploy workflow ([#27](https://github.com/fulcrumgenomics/ferro-hgvs/pull/27))
- Missing `in_memory` field in HgvsRsConfig construction ([#28](https://github.com/fulcrumgenomics/ferro-hgvs/pull/28))

### Changed

- Replaced vendored hgvs-rs with published crate v0.20.1 ([#22](https://github.com/fulcrumgenomics/ferro-hgvs/pull/22))

## [0.2.0](https://github.com/fulcrumgenomics/ferro-hgvs/compare/v0.1.0...v0.2.0) - 2026-03-23

### Fixed

- derive Python package version from Cargo.toml via maturin ([#14](https://github.com/fulcrumgenomics/ferro-hgvs/pull/14))
- allow valid HGVS characters in web service input validation ([#13](https://github.com/fulcrumgenomics/ferro-hgvs/pull/13))
- work-stealing for mutalyzer normalize shards ([#4](https://github.com/fulcrumgenomics/ferro-hgvs/pull/4))
- CIGAR-aware CDS-to-transcript coordinate mapping ([#7](https://github.com/fulcrumgenomics/ferro-hgvs/pull/7))
- prevent integer overflow and handle edge cases in normalization ([#6](https://github.com/fulcrumgenomics/ferro-hgvs/pull/6))
- prevent panics on malformed HGVS patterns during normalization ([#5](https://github.com/fulcrumgenomics/ferro-hgvs/pull/5))

### Other

- Add Bioconda and Zenodo badges ([#3](https://github.com/fulcrumgenomics/ferro-hgvs/pull/3))
- refresh test fixtures and data extraction scripts ([#8](https://github.com/fulcrumgenomics/ferro-hgvs/pull/8))
- release v0.1.0

## [0.1.0](https://github.com/fulcrumgenomics/ferro-hgvs/releases/tag/v0.1.0) - 2026-02-19

### Other

- Initial commit
