# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

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
