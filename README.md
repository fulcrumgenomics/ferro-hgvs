[![CI](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/ci.yml/badge.svg)](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/ci.yml)
[![Nightly reference-aware tests](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/nightly-mutalyzer.yml/badge.svg?branch=main)](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/nightly-mutalyzer.yml)
[![Codecov](https://codecov.io/gh/fulcrumgenomics/ferro-hgvs/branch/main/graph/badge.svg)](https://codecov.io/gh/fulcrumgenomics/ferro-hgvs)
[![Crates.io](https://img.shields.io/crates/v/ferro-hgvs.svg)](https://crates.io/crates/ferro-hgvs)
[![PyPI](https://img.shields.io/pypi/v/ferro-hgvs.svg)](https://pypi.org/project/ferro-hgvs/)
[![Python versions](https://img.shields.io/pypi/pyversions/ferro-hgvs.svg)](https://pypi.org/project/ferro-hgvs/)
[![Documentation](https://docs.rs/ferro-hgvs/badge.svg)](https://docs.rs/ferro-hgvs)
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
- **Variant Normalization**: 3'/5' shifting per HGVS specification
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
   maintainer judgment.

3. **Confluent.** Inputs denoting one variant produce one output. *Best effort.*
   Every rule is evaluated over the **resulting sequence**, never over the input's spelling.

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
   normalization form.** Options remain available for orthogonal axes — error mode, 3'/5' direction.

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
- **The variant's decomposition is not recoverable.** For an equal-length block the column
  correspondence is unique, so "two variants separated by N nucleotides" is a fact about the
  variant. For an unequal-length block there is no correspondence — recovering one means *choosing*
  an alignment, and the spec does not say which. There is no derivable form to converge on.

**"Best effort" is bounded by the spec's determinacy, not by ferro's implementation quality.** A
failure of rule 2 or 3 caused by ferro's code is a **bug** under rule 7. Only a failure caused by
the spec not determining an answer falls under best effort — and that triggers rule 5.

### A worked example of reading force from prose

`DNA/duplication.md` says a variant that *can* be described as a duplication **must** be — but the
"must" is scoped by the preceding clause, which defines when a duplication *can* be used at all:
only when the additional copy is directly 3'-flanking the original. So the rule ranks the *label*
for one span; it does not require that a partition be chosen so as to produce a duplication.
Reading the force without the scope inverts the rule.

### Known limitation

ferro cannot today guarantee that every input form is normalized according to these rules. They are
enforced intent, not a claim of current completeness.

Rule 7's disclosure mechanism — the `Representation-Change:` trailer and how it reaches the
changelog — is documented in
[CONTRIBUTING.md](CONTRIBUTING.md#declaring-a-representation-change).

## Supported HGVS Syntax

| Type | Prefix | Example |
|------|--------|---------|
| Genomic | `g.` | `NC_000001.11:g.12345A>G` |
| Coding DNA | `c.` | `NM_000088.3:c.459A>G` |
| Non-coding | `n.` | `NR_000001.1:n.100A>G` |
| RNA | `r.` | `NM_000088.3:r.459a>g` |
| Protein | `p.` | `NP_000079.2:p.Val600Glu` |
| Mitochondrial | `m.` | `NC_012920.1:m.3243A>G` |

### Edit Types

- Substitution: `A>G`, `Val600Glu`
- Deletion: `del`, `100_200del`
- Insertion: `100_101insATG`
- Deletion-Insertion: `100_102delinsATG`
- Duplication: `100_102dup`
- Inversion: `100_200inv`
- Repeat: `100CAG[20]`

## CLI Commands

The `ferro` CLI provides commands beyond parsing and normalization:

| Command | Description |
|---------|-------------|
| `prepare` | Download and prepare reference data for normalization |
| `check` | Verify reference data setup |
| `parse` | Parse and validate HGVS variants |
| `normalize` | Normalize HGVS variants (3'/5' shifting) |
| `explain` | Explain error/warning codes (e.g., `ferro explain W1001`) |
| `annotate-vcf` | Annotate VCF files with HGVS notation |
| `vcf-to-hgvs` | Convert VCF records to HGVS |
| `hgvs-to-vcf` | Convert HGVS to VCF format |
| `liftover` | Liftover coordinates between genome builds |
| `describe` | Generate HGVS from reference/observed sequences |
| `effect` | Predict protein effect from variant |
| `backtranslate` | Reverse translate protein to DNA variants |
| `convert-gff` | Convert GFF3/GTF to transcripts.json |
| `generate` | Generate HGVS descriptions from components |
| `extract-hgvs` | Extract HGVS from VEP-annotated VCFs |

## Error Handling

ferro-hgvs provides configurable error handling with three modes:

| Mode | Behavior |
|------|----------|
| `strict` | Reject non-conformant input (default) |
| `lenient` | Auto-correct with warnings |
| `silent` | Auto-correct silently |

```bash
# Use lenient mode to auto-correct common issues
ferro parse --error-mode lenient "p.val600glu"  # Corrects to p.Val600Glu

# Ignore specific warnings
ferro parse --ignore W1001,W2001 "p.val600glu"

# Get help on any error/warning code
ferro explain W1001
ferro explain --list
```

### Configuration File

Create `.ferro.toml` in your project directory:

```toml
[error-handling]
mode = "lenient"
ignore = ["W1001", "W2001"]  # Silently correct these
reject = ["W3003"]           # Always reject these
```

## Comparing normalization rules (`FERRO_PARTITION`)

> **Unstable.** `FERRO_PARTITION` is an evaluation switch, **not** a supported feature. It is not covered by semantic versioning, its values may change, and it is expected to be removed once the normalization rule is settled. Do not depend on it in production pipelines.

Normalization cuts each changed region of sequence into allele members. `FERRO_PARTITION` selects which rule does that cutting, so a candidate rule can be measured against the shipped one over a real corpus before anything changes for users.

| Value | Rule |
|-------|------|
| unset / empty / `live` | The shipped rule. This is what every normal invocation uses. |
| `shadow` | Cut only at alignment steps common to *every* minimal alignment. |
| `canonical` | The member-count-minimal minimal alignment. |
| `canonical-coalesced` | `canonical`, plus the `delins.md:44-47` merge: a split whose payload realigns as one block is re-spelled as a single `delins`. Applied after the downstream passes rather than at partition time, so it cuts identically to `canonical` and differs only in what survives. |

With the variable unset — or set to the empty string — output is byte-identical to a build with no switch at all.

### Running an A/B comparison

Run the same input twice and diff the normalized descriptions. In `tsv` format the columns are `line, input, normalized, changed, status, detail`, so column 3 is the normalized string:

```bash
# 1. the shipped behaviour
ferro normalize --input variants.txt --reference /path/to/reference \
  --format tsv --error-mode lenient -j 10 > shipped.tsv

# 2. the candidate rule
FERRO_PARTITION=canonical ferro normalize --input variants.txt --reference /path/to/reference \
  --format tsv --error-mode lenient -j 10 > candidate.tsv

# 3. what moved
diff <(cut -f3 shipped.tsv) <(cut -f3 candidate.tsv) | head

# 4. how many moved, counting only rows that succeeded on both sides
#    (columns 3/5 are `normalized`/`status`; a row that failed carries no
#     normalized string, so comparing column 3 alone would score two different
#     failures as identical)
paste <(cut -f3,5 shipped.tsv) <(cut -f3,5 candidate.tsv) \
  | awk -F'\t' '$2 == "ok" && $4 == "ok" && $1 != $3' | wc -l

# and, separately, rows whose status itself changed
paste <(cut -f5 shipped.tsv) <(cut -f5 candidate.tsv) | awk -F'\t' '$1 != $2' | wc -l
```

This works on a **stock release build** — the switch is not behind a build feature, so no special binary or wheel is needed.

### Two traps worth knowing

**A misspelled value is refused, loudly.** `FERRO_PARTITION=canonicl` aborts at the first variant ferro normalizes, naming the value you gave and the arms this build has:

```
FERRO_PARTITION="canonicl" is not a partitioner this build has. This build's arms
are: live, shadow, canonical, canonical-coalesced. Refusing rather than falling
back to `live`, because a bake-off served the shipped rule under a candidate's
name reports that the candidate changes nothing.
```

Naming *this build's* arms is the point: a value that exists on some other branch, or that used to exist, is reported as absent here rather than quietly answered as `live`.

Builds up to and including v0.13.1 fell back to `live` instead, and produced a clean, empty diff that read as "the candidate changes nothing". The fallback emitted a warning through the `log` facade, but the `ferro` CLI installs no logger, so that warning reached no stream and `RUST_LOG` could not surface it — there was no signal at all. If you are on an older build, treat every empty diff as unproven.

A positive control on an input **known** to differ is still worth running, because it catches the other way a comparison can be vacuous: a variable that never reached the process at all (a lost `export`, a `sudo` that scrubbed the environment, a container that did not forward it). That case is indistinguishable from `unset`, which is legitimately `live`, so no amount of validation inside ferro can catch it. On your own corpus a zero remains ambiguous between "the switch is not taking effect" and "this corpus has no affected variants".

These three inputs run against the built-in test data, so they need no `--reference` and no prepared reference directory. Between them they separate all four arms — every pair of arms disagrees on at least one row, so the control tells you *which* arm you got, not merely that something changed:

```bash
printf 'NM_001234.1:c.[5_6insAC;9del]\nNM_001234.1:c.[2del;5del]\nNM_001234.1:c.[2del;9del]\n' > control.txt

for arm in live shadow canonical canonical-coalesced; do
  echo "== $arm"
  if FERRO_PARTITION=$arm ferro normalize --input control.txt --format tsv \
       --error-mode lenient > "control.$arm.tsv"; then
    tail -n +2 "control.$arm.tsv" | cut -f3
  else
    echo "   FAILED (exit $?) -- ferro did not produce this arm's column"
  fi
done
```

The status check is not boilerplate. An unrecognised arm now aborts the process, and a pipeline reports the exit status of its *last* command — so `ferro … | tail | cut` reports `cut`'s success and the loop prints an empty column under the arm's heading. An empty column and an aborted run look identical, which is the same "a broken measurement reads as a result" failure this whole section is about. Redirecting first and testing the status makes the abort say so.

| input | `live` | `shadow` | `canonical` | `canonical-coalesced` |
|---|---|---|---|---|
| `c.[5_6insAC;9del]` | `c.6_9delinsACCAA` | `c.[5_6insAC;11del]` | `c.[5_6insAC;11del]` | `c.[5_6insAC;11del]` |
| `c.[2del;5del]` | `c.[2del;6del]` | `c.[2del;6del]` | `c.2_4delinsG` | `c.2_4delinsG` |
| `c.[2del;9del]` | `c.[2del;11del]` | `c.[2del;11del]` | `c.[2del;11del]` | `c.2_9delinsGCCCAA` |

Row 1 separates `live` from the other three, row 2 separates `{live, shadow}` from `{canonical, canonical-coalesced}`, and row 3 separates `canonical` from `canonical-coalesced`.

If the arm you selected does not produce its column above, the variable is not reaching ferro and any comparison you run is meaningless. Only once the control behaves is a zero on your own corpus informative.

> The control this section used to give did not work. It offered `c.[2del;9del]`, `c.[3del;9del]`, `c.[2del;9dup]` and claimed `canonical` answers `c.[2del;33del]`; on those three inputs **no arm differs from `live`**, so the documented check failed for a correct setup and would have been read as "the variable is not reaching ferro". The table above is verified against both a debug and a release build.

**From Python, it must be set before the first normalization.** The value is read once per process and cached, so:

```python
import os
os.environ["FERRO_PARTITION"] = "canonical"   # must precede the first normalize call
import ferro_hgvs
```

Setting it after any variant has been normalized silently does nothing, and it cannot be changed within a running process. Comparing two rules therefore means two separate processes (or two CLI runs, as above), not two calls in one script.

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
