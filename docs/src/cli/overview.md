# CLI reference

The `ferro` command groups its work into subcommands. Run `ferro <command> --help` for the full,
authoritative options of any one — this page is a map, not a substitute.

## Core

| Command | What it does |
|---|---|
| `parse` | Parse and validate an HGVS description (no reference needed). |
| `normalize` | Rewrite descriptions into canonical form (needs a reference). See [Normalize variants](../guide/normalize-variants.md). |
| `project` | Re-express a variant on a chosen output axis (`g` / `c` / `n` / `p` / `r`). |
| `mutalyzer` | Show a Mutalyzer-style view: the normalized form and every axis at once. See [A Mutalyzer-style multi-axis view](../guide/mutalyzer-view.md). |
| `explain` | Explain an error or warning code (e.g. `ferro explain W3003`). |

## Reference data

| Command | What it does |
|---|---|
| `prepare` | Download and assemble reference data. See [Reference data](../guide/reference-data.md). |
| `check` | Verify a reference directory is ready (optionally pre-build the cache). |
| `convert-gff` | Convert a GFF3/GTF annotation to `transcripts.json`. |
| `build-transcript` | Build a single-exon `transcripts.json` from a FASTA + CDS coordinates. |

## VCF and interchange

| Command | What it does |
|---|---|
| `annotate-vcf` | Annotate a VCF file with HGVS notation. |
| `vcf-to-hgvs` / `hgvs-to-vcf` | Convert between VCF and HGVS. |
| `extract-hgvs` | Extract HGVS patterns from VEP-annotated VCF files. |
| `liftover` | Lift genomic coordinates between genome builds. |

## Generation and prediction

| Command | What it does |
|---|---|
| `generate` | Generate an HGVS description from components. |
| `describe` | Generate a description from reference and observed sequences. |
| `effect` | Predict the protein effect of a variant. |
| `backtranslate` | Backtranslate a protein variant to possible DNA variants. |

## Common options

Most commands accept `-i/--input <file>` (one description per line) and `-o/--output <file>`. The
output format is set with `-f/--format`, and the accepted values depend on the command: `text|json`
for `parse` and `project`, `text|json|tsv` for `normalize`, and `text|json|markdown` for `explain`.
Commands that normalize also accept `--reference <dir>` and `--error-mode <strict|lenient|silent>`.
