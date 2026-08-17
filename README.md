[![CI](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/ci.yml/badge.svg)](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/ci.yml)
[![Nightly reference-aware tests](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/nightly-mutalyzer.yml/badge.svg?branch=main)](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/nightly-mutalyzer.yml)
[![Codecov](https://codecov.io/gh/fulcrumgenomics/ferro-hgvs/branch/main/graph/badge.svg)](https://codecov.io/gh/fulcrumgenomics/ferro-hgvs)
[![Crates.io](https://img.shields.io/crates/v/ferro-hgvs.svg)](https://crates.io/crates/ferro-hgvs)
[![PyPI](https://img.shields.io/pypi/v/ferro-hgvs.svg)](https://pypi.org/project/ferro-hgvs/)
[![Python versions](https://img.shields.io/pypi/pyversions/ferro-hgvs.svg)](https://pypi.org/project/ferro-hgvs/)
[![Documentation](https://readthedocs.org/projects/ferro-hgvs/badge/?version=stable)](https://ferro-hgvs.readthedocs.io/en/stable/)
[![API docs](https://docs.rs/ferro-hgvs/badge.svg)](https://docs.rs/ferro-hgvs)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ferro-hgvs/README.html)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.18703103-blue.svg)](https://doi.org/10.5281/zenodo.18703103)

# ferro-hgvs

A high-performance HGVS variant nomenclature parser and normalizer written in Rust.

**WARNING: ALPHA SOFTWARE - USE AT YOUR OWN RISK**

This software is currently in **ALPHA**. While we have extensively tested it
across a wide variety of HGVS patterns, **no guarantees are made** regarding
correctness or stability.

<p>
<a href="https://fulcrumgenomics.com">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/fulcrumgenomics/ferro-hgvs/main/.github/logos/fulcrumgenomics-dark.svg">
    <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/fulcrumgenomics/ferro-hgvs/main/.github/logos/fulcrumgenomics-light.svg">
    <img alt="Fulcrum Genomics" src="https://raw.githubusercontent.com/fulcrumgenomics/ferro-hgvs/main/.github/logos/fulcrumgenomics-light.svg" height="100">
  </picture>
</a>
</p>

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-brightgreen.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-blue.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>

## Features

- **Full HGVS Parsing**: All coordinate systems (g/c/n/r/p/m/o) and edit types
- **Variant Normalization**: 3' shifting per HGVS specification
- **High Performance**: ~5M variants/sec single-threaded parsing (>12M/s parallel), zero-copy with nom
- **Type-Safe**: Leverages Rust's type system for correctness

## Installation

### Python

```bash
pip install ferro-hgvs
```

Pre-built wheels are available for Linux (x86_64, aarch64), macOS (x86_64, Apple Silicon), and Windows (x86_64) on Python 3.10+.

### Rust

Add to your `Cargo.toml`:

```toml
[dependencies]
ferro-hgvs = "0.1"
```

Or install the CLI:

```bash
cargo install ferro-hgvs
```

## Quick Start

### CLI

```bash
# Parse a variant
ferro parse "NM_000088.3:c.459A>G"

# Parse from file
ferro parse -i variants.txt -f json

# Prepare reference data (downloads RefSeq, genome, cdot — RefSeq-only by default)
ferro prepare --output-dir ferro-reference

# Verify reference data is ready
ferro check --reference ferro-reference

# (Optional) pre-build the on-disk cdot cache as a setup step, so the one-time
# cache build doesn't slow the start of a real (or timed/benchmarked) run.
ferro check --reference ferro-reference --build-cache

# Normalize with reference
ferro normalize "NM_000088.3:c.459del" --reference ferro-reference/
```

> **Read the warnings.** Normalization sometimes *repairs* a description in a way
> the normalized string does not record — separately reported cis members merged
> into one `delins` (`MEMBERS_COALESCED_FROM_REPORTED_FORM`), a `ins[100_110]`
> reference-range payload replaced by the bases it denotes
> (`INSERTED_SEQUENCE_EXPANDED`), a stated reference base that contradicted the
> reference and was accepted anyway (`REFSEQ_MISMATCH`). Those are reported as
> `warning[CODE]: message` on **stderr** (and in the `warnings` array under
> `--format json`, the `detail` column under `--format tsv`), so a pipeline
> reading only stdout will not see them. `--error-mode strict` is not a
> substitute: it rejects a specific ladder of conditions and reports the rest
> exactly as lenient does.

> **Throughput tip:** when normalizing many variants, feed them **sorted by transcript accession (or by genomic position)**. ferro caches each resolved transcript, so consecutive variants on the same transcript skip the (dominant) cost of re-reading and re-building it from the reference. Sorted input keeps the relevant transcripts resident in the cache and is markedly faster on large batches — see [Performance Comparison](#performance-comparison).

### Optional reference data

A bare `ferro prepare` builds a **RefSeq-only** reference (accessions `NM_`/`NR_`/`NP_`/`NG_`). Two opt-in flags provision additional data — pass them at prepare time; they are what a fully-provisioned ("blessed") reference is built with:

```bash
# Add Ensembl support (accessions ENST/ENSG/ENSP). Downloads the Ensembl cdot
# metadata and cDNA FASTAs (~1 GB+); off by default. Without it, an ENST/ENSG/ENSP
# input reports "Reference not found" and the message points back at this flag.
ferro prepare --output-dir ferro-reference --ensembl

# Derive version-independent NG_ placements and the NG_→transcript-version map
# (ng_hosted_transcripts) for a curated list of RefSeqGene accessions. Required to
# resolve legacy gene-symbol selectors (NG_(GENE):c.…) and bare-NG_ hosted lookups.
ferro prepare --output-dir ferro-reference \
  --derive-ng-placements path/to/ng_accessions.txt

# A fully-provisioned reference combines both in one run:
ferro prepare --output-dir ferro-reference --ensembl \
  --derive-ng-placements path/to/ng_accessions.txt
```

Both flags are incremental: re-running `ferro prepare` over an existing reference adds the requested data and preserves already-provisioned artifacts.

### Library

```rust
use ferro_hgvs::{parse_hgvs, HgvsVariant};

fn main() -> Result<(), ferro_hgvs::FerroError> {
    let variant = parse_hgvs("NM_000088.3:c.459A>G")?;

    match &variant {
        HgvsVariant::Cds(v) => println!("CDS variant: {}", v),
        HgvsVariant::Genome(v) => println!("Genomic variant: {}", v),
        _ => println!("Other: {}", variant),
    }

    Ok(())
}
```

### Python

```python
import ferro_hgvs

# Parse a variant
variant = ferro_hgvs.parse("NM_000088.3:c.459A>G")
print(variant.variant_type)  # "coding"
print(variant.reference)     # "NM_000088.3"
print(str(variant))          # "NM_000088.3:c.459A>G"

# Normalize with reference data
normalizer = ferro_hgvs.Normalizer(reference_json="ferro-reference/cdot.json")
normalized = normalizer.normalize("NM_000088.3:c.459del")

# `normalize` returns only the string, so it cannot tell you that normalization
# repaired something. `normalize_with_warnings` returns the same string plus the
# diagnostics — as a free function, or as a Normalizer method.
result = normalizer.normalize_with_warnings("NM_000088.3:c.459del")
print(str(result.result))                      # the same normalized string
print([(w.code, w.message) for w in result.warnings])
```

## Documentation

Full guides live in **[the documentation site](https://ferro-hgvs.readthedocs.io/)** (source under [`docs/src/`](docs/src/)):

- [Deriving a description from sequences](docs/src/guide/deriving-from-sequences.md) — turn a window of bases into one canonical description, no reference needed
- [Normalize variants](docs/src/guide/normalize-variants.md) · [Reference data](docs/src/guide/reference-data.md)
- [Error handling](docs/src/guide/error-handling.md) — strict / lenient / silent modes and warning codes
- [Comparing normalization rules (`FERRO_PARTITION`)](docs/src/guide/comparing-rules.md)
- [Benchmarking](docs/src/guide/benchmarking.md)
- [CLI reference](docs/src/cli/overview.md) · [Supported HGVS syntax](docs/src/reference/hgvs-syntax.md)
- [Interpreting the HGVS recommendations](docs/src/shadow-spec/index.md) — how ferro reads each clause


## Normalization rules

ferro's normalizer follows four rules about its output, and three about how it handles the gaps.

### The output contract

1. **Conformant.** Output follows the HGVS recommendations. *Absolute — never traded.*
   Scope: the `syntax.yaml` grammar, plus every prohibition, read from **prose force** rather than
   keyword casing — "not allowed", "not correct", "can not be used", "by definition", the
   `class="invalid"` markup, and the set `checklist.md` enumerates.

2. **Recommended form.** Where the spec prefers among conformant forms, ferro produces it.
   *Best effort.*
   Scope: "recommended", "preferred", lowercase "should", the 3' rule. A preference clause outranks
   maintainer judgment, but not rule 3: where it is not evaluable on what a normalizer holds, rule 3
   governs.

3. **Confluent.** Inputs denoting one variant produce one output. *Best effort.*
   Every rule — rule 2 included — is evaluated over the **resulting sequence**, never over the
   input's spelling. Reference context counts as sequence: transcript model, exon boundaries,
   reading frame, strand and topology are functions of the accession, not of the description.

4. **Deterministic.** Same input, same output. *Absolute.*
   Note that 4 does not imply 3 — a deterministic normalizer can be arbitrarily non-confluent.

### The procedure

5. **Where the spec is silent, ambiguous, or self-contradictory:** the issue is filed upstream
   **first** and cited, and only then does ferro ship a provisional choice.
    - **Self-contradictory** — two clauses that cannot both hold — is a **defect**. Every conforming
      tool must pick a side and none of them can be right, so filing is a bug report and is not
      optional.
    - **Silent** is merely incomplete. Ferro decides under rule 6 and violates nothing; filing is a
      feature request, worth making but never a reason to hold a release.

6. **Among multiple conformant forms:** the maintainers choose. **There are no user options for
   normalization form.** Error mode is an orthogonal axis and stays available. The 3'/5' shuffle
   direction was the one exception, and it is now removed from the public surface rather than
   excused: it is **not** orthogonal — it selects the frame every rule is evaluated in — so it was
   a user option for normalization form sitting inside this rule. ferro shifts 3', the only
   direction the HGVS recommendations describe, and no CLI flag, Python keyword or service config
   key selects otherwise. The 5' arm survives internally as a differential oracle over ferro's own
   test suite; an instrument is not a user option, and rules 2 and 3 are now claimed once rather
   than per direction.

7. **Disclosure.** Any change to these rules, and any different choice made under 5 or 6, is
   disclosed: in the changelog before v1, by a major version bump after. Output that *violates*
   rules 1-4 is a **bug**, not a disclosure.

### Why 2 and 3 are best effort, and 1 and 4 are not

Rules 1 and 4 are always achievable. Conformance is checkable against the spec text, and
determinism is a property of ferro's own code; nothing external can prevent us honouring them.

Rules 2 and 3 depend on the spec **determining an answer**, and sometimes it does not:

- **No preference exists.** The spec ranks substitution, deletion, inversion, duplication and
  insertion, but says nothing about competing `delins` forms.
- **Two preferences disagree.** `general.md` ranks duplication above insertion; `DNA/inversion.md`
  prefers an insertion for inverted copies. No single output satisfies both.
- **The same clause has two versions.** `general.md`'s current text and its forthcoming NOTE give
  opposite answers for variants separated by one nucleotide.
- **The variant's decomposition is not recoverable.** Recovering one means *choosing* an alignment,
  and the spec does not say which, so there is no derivable form to converge on. **Block length does
  not settle it**: an equal-length block can still carry a balanced del+ins pair, so its column
  correspondence need not be unique — `CAG -> AGA` is equal length with edit distance 2, not the
  position-wise 3. What decides whether a reference base is unchanged is whether **every minimal
  alignment** matches it, so the property to key on is edit distance against block length, never
  length alone. See `rulings[unchanged-is-read-over-every-minimal-alignment]`.
- **The preference keys on information a normalizer does not hold.** A "frequently occurring
  variant" (`RNA/delins.md:41`); a repeat "variable in the population" (`RNA/repeated.md:33`,
  `protein/repeated.md:22`); two variants "reported (or might occur) individually" (`DNA/delins.md:83`).
  The spec determines an answer; ferro cannot see what it keys on.

**"Best effort" is bounded by the spec's determinacy and by what a normalizer's inputs can decide —
not by ferro's implementation quality.** A failure of rule 2 or 3 caused by ferro's code is a
**bug** under rule 7; one caused by the spec not determining an answer triggers rule 5; one caused
by a clause keying on provenance or population data is a **declared deviation**, and a permanent
one — no upstream answer can put that information into a normalizer's hands.

Permanent is the only thing the last case adds, and it does **not** exempt the question from rule
5. The grain matters, because the two grains give opposite answers. Read **on its own**, such a
clause is not silent, ambiguous or self-contradictory — it determines an answer perfectly well, and
ferro is the one that cannot see what the answer keys on. Read **against the spec's own
re-derivation mandate** — `general.md:157-160`, which has a protein description derived by
comparing the variant and reference protein sequences and says knowledge of the underlying DNA
change "should not be used", with `general.md:13` extending the method to RNA — it is one half of a
pair that cannot both hold. That is rule 5's self-contradiction limb exactly, and the ruling ledger
says so in those words. What the carve-out states is that rule 5's escalation cannot *end*
the matter here: a choice made under rule 5 is provisional pending an upstream ruling, and no
ruling upstream makes provenance visible, so the deviation outlives the filing.

Where it is recorded is worth stating precisely, because the obvious pointer does not resolve.
`canonical-form-choice-when-both-legal` has **no `deviates_from` field**: it carries all four
clauses in its rationale as the recorded counter-evidence to re-derivation, opens that paragraph
with "The spec contradicts itself here", and adopts re-derivation over them deliberately. The one
of the four that is recorded *as* a deviation is `DNA/delins.md:83`, through the
`deviates_from: ["docs/recommendations/DNA/delins.md:79-84"]` on
`separation-is-a-property-of-the-spelling-not-of-the-variant` — the record whose ruling the rule 3
example below turns on.

### A worked example of reading force from prose

`DNA/duplication.md` says a variant that *can* be described as a duplication **must** be — but the
"must" is scoped by the preceding clause, which defines when a duplication *can* be used at all:
only when the additional copy is directly 3'-flanking the original. So the rule ranks the *label*
for one span; it does not require that a partition be chosen so as to produce a duplication.
Reading the force without the scope inverts the rule.

### What rule 3 excludes

"Never over the input's spelling" is narrower than it sounds. Most of a description is context.

| Carried by the description | Treatment |
|---|---|
| Accession, axis (`g.`/`c.`/`n.`), version | **Used** — it *is* the reference context |
| Which bases end up different | **Used** — this is the variant |
| Cis against trans (`[a;b]` against `[a];[b]`) | **Used** — different variants, not two spellings of one |
| Type label (`dup` against `ins`, `inv` against `delins`) | **Re-derived**, then ranked by `general.md` |
| How the edit set is cut into members | **Excluded** — a property of who wrote the string |
| Which copy in a run of identical residues a member names | **Excluded** — the 3' rule assigns this *"arbitrarily"* (`general.md:41`) |
| Repeat unit and phase, where several are equivalent | **Excluded** |

So three rows read **Excluded**, and they are three spellings of the one thing the input does not
get to decide: the **partition**, the run-position choice that feeds it, and — where several
unit-and-phase pairs describe one tract equally well — which pair a repeat member names.

`NC_000001.11:g.1001002_1001016` reads `ATGAGGGGCCACTGT`: a `GGGG` run at `1001006-1001009`, a `CC`
run at `1001010-1001011`, a lone `C` at `1001013`. Two spellings, one denoted sequence
(`ATGAGGGCATGT`), because `1001010` and `1001011` are both `C`:

```
g.[1001009del;1001010del;1001013del]     written gaps of 0 and 2
g.[1001009del;1001011del;1001013del]     written gaps of 1 and 1
```

`general.md:34`'s "two variants separated by one or more nucleotides should be described individually"
reads those gaps, so it answers twice for one variant. Rule 3 reads them off the partition ferro
derives instead: both give `g.[1001009_1001010del;1001013del]`. Rule 2 then keeps that over the
spanning `g.1001009_1001013delinsCA`, which merges across two unchanged nucleotides. Pinned in
[`tests/it/cis_confluence_adjudication.rs`](tests/it/cis_confluence_adjudication.rs).

### Known limitation

ferro cannot today guarantee that every input form is normalized according to these rules. They are
enforced intent, not a claim of current completeness.

Rule 7's disclosure mechanism — the `Representation-Change:` trailer and how it reaches the
changelog — is documented in
[CONTRIBUTING.md](CONTRIBUTING.md#declaring-a-representation-change).

### Where the individual decisions live

Rules 5 and 6 above say *how* a question is decided where the recommendations are silent,
ambiguous or self-contradictory. **What was decided, case by case**, is recorded separately, as
adjudication records — each naming the clauses in tension, which one governs, which is deviated
from, and why. Those records are published in full, with their clause quotes, in
[docs/NORMALIZATION_CONTRACT.md](docs/NORMALIZATION_CONTRACT.md).

That document is generated from the records and gated against them, so it cannot drift; it
deliberately does not restate the seven rules above, which are stated only here.

The inverse index — for each stage of the normalizer, *which* record or clause governs it, and
which decisions are governed by **nothing** — is
[docs/NORMALIZATION_STAGE_AUDIT.md](docs/NORMALIZATION_STAGE_AUDIT.md). Read that one if what
you want to know is whether a behaviour you are looking at was chosen or merely happened.

## Why ferro-hgvs?

ferro-hgvs provides the most comprehensive HGVS variant normalization across all pattern types, with performance orders of magnitude faster than alternatives.

### Normalization Capabilities Comparison

<!-- DO NOT EDIT — generated from docs/tool_support_matrix.json by `generate_tool_support_tables`. -->
<!-- BEGIN tool-support:normalization_capabilities -->
| Pattern Type | ferro | mutalyzer | biocommons | hgvs-rs |
|--------------|:-----:|:---------:|:----------:|:-------:|
| Genomic (g.) | ✓ | ✓ | ✓ | ✓ |
| Coding (c.) exonic | ✓ | ✓ | ✓ | ✓ |
| Coding (c.) intronic | ✓ | ✓** | ✗ | ✗ |
| Non-coding (n.) | ✓ | ✓ | ✓ | ✓ |
| RNA (r.) | ✓ | ✓ | ✓ | ✓ |
| Protein (p.) | ✓ | Net* | ✗ | ✗ |

\* mutalyzer protein normalization requires network access for NP_→NM_ lookups (cannot be cached locally).
\*\* mutalyzer intronic support is enabled by default via genomic-context rewriting; disable with --no-rewrite-intronic.
<!-- END tool-support:normalization_capabilities -->

### Performance Comparison

**All tools are benchmarked in ferro's offline configuration — best case for every tool.** Reference data is preloaded locally (a local UTA database and SeqRepo) and the network is disabled, so the figures below measure parse/normalize compute, not I/O. Out of the box, hgvs-rs, biocommons/hgvs, and mutalyzer resolve each variant against a remote UTA/SeqRepo or the Mutalyzer web API — a network round-trip per variant (~100–1000 ms), i.e. roughly **1–10 variants/sec, hundreds to thousands of times slower than shown here** (an order-of-magnitude estimate from per-call network latency, not separately benchmarked). That local, offline setup is exactly what ferro's `prepare` command builds; ferro needs no external service.

<!-- DO NOT EDIT — generated from data/benchmark/perf_results.json by `generate_perf_tables`. -->

_Median patterns/sec over 5 reps on an Apple M2 Max, local/offline. All tools draw from one stratified ClinVar population; per-tool sample sizes are calibrated so each tool is measured over a meaningful interval — fast cells (e.g. ferro/hgvs-rs parse) draw from millions of patterns, while slower cells (e.g. the per-tool normalize columns) draw from as few as tens to thousands. All tools exclude process/interpreter startup from the timed region — the mutalyzer/biocommons Python subprocesses are timed by their own internal startup-excluded timer, matching ferro/hgvs-rs. Only ferro parallelizes natively (rayon); the other tools are single-threaded libraries, so their normalize *@8 workers* figures come from the benchmark harness running 8 independent instances in parallel, while parsing is not sharded for them — hence the `single-threaded` label in their parse *@8 workers* column (mutalyzer normalize likewise shows no gain at 8 workers: per-call cache and IPC overhead dominate, so sharding does not help). Every tool runs fully offline against local reference data — a local UTA database and SeqRepo, with mutalyzer's network lookups disabled — the configuration ferro's `prepare` command enables; the figures therefore reflect compute throughput, not per-variant network latency. Reference-data load is excluded for all tools. ferro full-population peak: parse 20.0M/s, normalize 77.0k/s. See `docs/BENCHMARK_RUNBOOK.md` for the full method._

**Parse**

<!-- BEGIN perf:parse -->
| Tool | Throughput @ 1 worker | Throughput @ 8 workers | ferro speedup @ 8w |
|------|----------------------:|-----------------------:|-------------------:|
| ferro | 5.1M/s | 12.2M/s | — |
| mutalyzer | 352/s | single-threaded | 35,000× |
| biocommons | 3.9k/s | single-threaded | 3,100× |
| hgvs-rs | 3.6M/s | single-threaded | 3× |
<!-- END perf:parse -->

**Normalize**

<!-- BEGIN perf:normalize -->
| Tool | Throughput @ 1 worker | Throughput @ 8 workers | ferro speedup @ 8w |
|------|----------------------:|-----------------------:|-------------------:|
| ferro | 78.1k/s | 260.2k/s | — |
| mutalyzer | 4/s | 4/s | 73,000× |
| biocommons | 368/s | 818/s | 320× |
| hgvs-rs | 195/s | 1.3k/s | 200× |
<!-- END perf:normalize -->

**ferro thread scaling**

<!-- BEGIN perf:ferro_scaling -->
| Threads | 1 | 2 | 4 | 8 |
|---------|--:|--:|--:|--:|
| ferro parse | 5.1M/s | 9.4M/s | 16.0M/s | 12.0M/s |
<!-- END perf:ferro_scaling -->

**Input ordering matters for batch throughput.** Resolving a transcript (reading its full sequence from the reference and rebuilding its CDS/exon metadata) dominates per-variant cost. ferro memoizes resolved transcripts in a bounded in-memory cache, so repeated lookups of the same transcript are near-free. Providing variants **sorted by transcript accession — or by genomic position, which clusters variants onto the same transcripts** — maximizes the cache hit rate and can speed up large batches by an order of magnitude versus randomly-ordered input. Ordering matters most when the number of distinct transcripts in the run exceeds the cache capacity (very large or genome-wide inputs); below that, the working set stays resident regardless of order.

### Reference Data: What ferro Prepares

The `ferro prepare` command downloads and organizes all reference data needed for comprehensive normalization. This data is then shared with other tools (mutalyzer, biocommons, hgvs-rs) to enable their local operation.

| Data Type | Source | Size | Enables |
|-----------|--------|------|---------|
| **RefSeq transcripts** | NCBI | ~1GB | NM_/NR_/XM_ normalization |
| **cdot metadata** | MANE | ~200MB | Transcript-to-genome mappings |
| **GRCh38 + GRCh37 genomes** | NCBI | ~4GB | NC_ genomic normalization |
| **RefSeqGene** (sequences + genome alignments) | NCBI | ~600MB | NG_ gene-region normalization; projecting c./n. variants into an NG_ parent's own g. frame (via the RefSeqGene→genome alignment GFF3) |
| **LRG sequences + XML** | EBI | ~50MB | LRG_ stable-reference normalization; projecting c./n. variants into an LRG_ parent's own g. frame (via the LRG XML genomic mapping) |
| **Protein sequences** | Derived from CDS | ~200MB | NP_/XP_ protein normalization |
| **Legacy transcript versions** | NCBI | ~50MB | Historical ClinVar variants |

**Key insight**: Without ferro's reference preparation, other tools require network access for each variant lookup (adding 100-1000ms latency per variant). With ferro's cached reference data, all tools can operate fully offline with consistent, reproducible results.

#### Deriving version-independent NG_ placements (#728)

`ferro prepare --derive-ng-placements <accessions.txt>` derives genomic placements for the listed `NG_` versions (one exact accession per line, e.g. `NG_012337.3`; blank lines and `#` comments ignored), writing `derived_refseqgene_placements.json` into the reference directory and wiring the manifest's `derived_refseqgene_placements` field. This fills version gaps the archived RefSeqGene→genome GFF3 snapshots do not cover. It needs cdot + the genome in the same prepare run and uses NCBI EFetch per accession; accessions that cannot be validated are skipped with a warning. The field is preserved across subsequent `prepare` runs.

## Benchmark: Reference Data & Tool Comparison

The main `ferro` binary includes commands to prepare reference data (`ferro prepare`) and check its status (`ferro check`). The `ferro-benchmark` tool (build with `--features benchmark`) extends this for tool comparison benchmarks.

| Command | Description |
|---------|-------------|
| `prepare <tool>` | Prepare reference data for a tool |
| `check <tool>` | Verify tool configuration and dependencies |
| `parse <tool>` | Parse HGVS patterns with specified tool |
| `normalize <tool>` | Normalize HGVS patterns with specified tool |
| `compare results` | Compare parse/normalize results between tools |
| `extract` | Extract patterns from ClinVar, VCFs, or create samples |
| `setup` | Set up UTA database, SeqRepo, and other services |
| `generate` | Generate summary reports and configs |
| `collate` | Aggregate sharded results |

### Quick Start

```bash
# Prepare ferro reference (main binary - no special features needed)
ferro prepare --output-dir data/ferro

# Check reference data
ferro check --reference data/ferro

# Normalize with ferro
ferro normalize -i patterns.txt --reference data/ferro

# For tool comparison, build with benchmark support
cargo build --release --features benchmark

# Prepare other tools (uses ferro reference for transcript data)
ferro-benchmark prepare mutalyzer --ferro-reference data/ferro --output-dir data/mutalyzer
ferro-benchmark prepare biocommons --seqrepo-dir data/seqrepo --uta-dump uta_20210129b.pgd.gz --ferro-reference data/ferro

# Compare results between tools
ferro-benchmark normalize mutalyzer -i patterns.txt -o mutalyzer.json --mutalyzer-settings data/mutalyzer/mutalyzer_settings.conf
ferro-benchmark compare results normalize ferro.json mutalyzer.json -o comparison.json
```

**Supported tools**: ferro-hgvs, mutalyzer, biocommons/hgvs, hgvs-rs

> **Note**: The `pixi.toml` and `pixi.lock` files in this repository define a [pixi](https://pixi.sh) environment for the Python-based external tools (mutalyzer, biocommons/hgvs, seqrepo) used in benchmarking. Run `pixi shell` to activate it.

See [docs/BENCHMARK_GUIDE.md](docs/BENCHMARK_GUIDE.md) for detailed usage.

## Development

```bash
cargo build
cargo test                        # default features
cargo clippy -- -D warnings
```

The commands above use the default feature set, and CI keeps them compiling (see
the `build` job). They do not cover the whole suite — the feature-gated tests and
the integration tree need `dev`, which is what CI runs and what you want before
opening a PR:

```bash
cargo nextest run --features dev
cargo clippy --features dev --all-targets -- -D warnings
```

## License

Licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Disclaimer

This software is under active development.
While we make a best effort to test this software and to fix issues as they are reported, this software is provided as-is without any warranty (see the [license](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/LICENSE) for details).
Please submit an [issue](https://github.com/fulcrumgenomics/ferro-hgvs/issues), and better yet a [pull request](https://github.com/fulcrumgenomics/ferro-hgvs/pulls) as well, if you discover a bug or identify a missing feature.
Please contact [Fulcrum Genomics](https://www.fulcrumgenomics.com) if you are considering using this software or are interested in sponsoring its development.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.
