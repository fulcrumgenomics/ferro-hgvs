# Third-Party Notices

This file documents the sources and licenses of third-party data and code used in
ferro-hgvs.

## Test Data Sources

### Public Domain / CC0

The following test data sources are in the public domain and have no restrictions:

**NCBI ClinVar**
- Files: `tests/fixtures/bulk/clinvar_*.json`, `tests/fixtures/validation/clinvar.json`
- Source: https://www.ncbi.nlm.nih.gov/clinvar/
- License: Public Domain (NCBI places no restrictions on use or distribution)
- Citation: Landrum MJ, et al. ClinVar: public archive of interpretations of
  clinically relevant variants. Nucleic Acids Res. 2016.

**NCBI dbSNP**
- Files: `tests/fixtures/external/ncbi_dbsnp.txt`
- Source: https://www.ncbi.nlm.nih.gov/snp/
- License: Public Domain

**HGVS Nomenclature Specification**
- Files: `tests/fixtures/grammar/hgvs_spec_examples.json`, `tests/fixtures/grammar/hgvs_canonical_examples.json`
- Source: https://hgvs-nomenclature.org/
- License: CC0 1.0 Universal (Public Domain Dedication)

**MaveDB**
- Files: `tests/fixtures/external/mavedb_functional.json`
- Source: https://www.mavedb.org/
- License: CC0 (as of 2024)
- Citation: Esposito D, et al. MaveDB: an open-source platform to distribute and
  interpret data from multiplexed assays of variant effect. Genome Biol. 2019.

### MIT License

**mutalyzer/hgvs-parser**
- Files: `tests/fixtures/grammar/mutalyzer_github.json`, `tests/fixtures/normalization/mutalyzer.json`
- Source: https://github.com/mutalyzer/hgvs-parser
- License: MIT
- Note: Test patterns extracted from the MIT-licensed hgvs-parser repository,
  not from the main mutalyzer project.

### Apache 2.0 + CC-BY 4.0

**biocommons/hgvs**
- Files: `tests/fixtures/grammar/biocommons.json`, `tests/fixtures/external/biocommons_*.txt`
- Source: https://github.com/biocommons/hgvs
- License: Apache License 2.0 (code), Creative Commons CC-BY 4.0 (data)
- Attribution: biocommons contributors, https://github.com/biocommons/hgvs
- Citation: Hart RK, et al. A Python package for parsing, validating, mapping,
  and formatting sequence variants using HGVS nomenclature. Bioinformatics. 2015.

### Synthetic / Project-Generated Data

The following test data was generated synthetically by this project and is
released under the same MIT license as the project itself:

- `tests/fixtures/bulk/gnomad_*.json` - Synthetic patterns based on gnomAD variant
  type distributions (NOT actual gnomAD data)
- `tests/fixtures/transcripts/mock_transcripts.json` - Synthetic transcript data
- `tests/fixtures/normalization/normalization.json` - Project-generated test cases
- `tests/fixtures/edge_cases/*` - Project-generated edge case tests

## Dependencies

This project uses the following Rust crates as dependencies. See `Cargo.toml`
for the complete list and versions. All dependencies use permissive licenses
(MIT, Apache-2.0, or similar) compatible with this project's MIT license.

### Optional Comparison Dependencies

When built with the `hgvs-rs` feature, this project includes:

**hgvs-rs**
- Source: https://github.com/varfish-org/hgvs-rs
- License: Apache License 2.0
- Note: Used only for benchmarking comparison, not as a core dependency.

## External Tools (Not Included)

The benchmark feature can compare results against external tools. These tools are
NOT included in this repository and must be installed separately:

- **mutalyzer** (https://github.com/mutalyzer/mutalyzer) - AGPL v3.0
- **biocommons/hgvs** (https://github.com/biocommons/hgvs) - Apache 2.0
- **VariantValidator** (https://github.com/openvar/variantvalidator) - AGPL v3.0

These tools are called via HTTP API or subprocess and no code from them is
incorporated into ferro-hgvs.
