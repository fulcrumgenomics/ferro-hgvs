# ferro-benchmark Comprehensive Guide

This guide walks through running a complete benchmark comparing ferro-hgvs against mutalyzer, biocommons/hgvs, and hgvs-rs.

## Table of Contents

1. [Overview](#overview)
2. [Quick Start](#quick-start)
3. [Prerequisites & Build](#prerequisites--build)
4. [Download ClinVar Data](#download-clinvar-data)
5. [Network Fetch Behavior](#network-fetch-behavior)
6. [Network vs Local Performance Comparison](#network-vs-local-performance-comparison)
7. [UTA Alignment Loading](#uta-alignment-loading)
8. [Recommended Workflow](#recommended-workflow)
9. [Prepare All Tools](#prepare-all-tools)
10. [Check All Tools](#check-all-tools)
11. [Parse Comparison: Small & Medium](#parse-comparison-small--medium)
12. [Normalize Comparison: Small & Medium](#normalize-comparison-small--medium)
13. [Parse Comparison: Full ClinVar](#parse-comparison-full-clinvar)
14. [Normalize Comparison: Full ClinVar](#normalize-comparison-full-clinvar)
15. [Troubleshooting](#troubleshooting)
16. [CLI Command Reference](#cli-command-reference)

---

## Overview

### Tools Compared

| Tool | Type | Language | Speed (local) | Speed (network) | ferro Speedup |
|------|------|----------|---------------|-----------------|---------------|
| **ferro-hgvs** | Parser + Normalizer | Rust | ~4M patterns/sec | N/A | N/A |
| **mutalyzer** | Parser + Normalizer | Python | ~20 patterns/sec | ~1 pattern/sec | **20x** |
| **biocommons/hgvs** | Parser + Normalizer | Python | ~20 patterns/sec | ~0.2 patterns/sec | **100x** |
| **hgvs-rs** | Parser + Normalizer | Rust | ~2 patterns/sec | ~0.2 patterns/sec | **10x** |

### Versions Used in This Guide

Results in this guide were generated with:

| Component | Version |
|-----------|---------|
| ferro-hgvs | 0.1.0 |
| mutalyzer-hgvs-parser | 0.3.9 |
| mutalyzer-algebra | 1.5.2 |
| biocommons hgvs | 1.5.6 |
| hgvs-rs | 0.19.1 |
| cdot | 0.2.32 (RefSeq GRCh38) |
| SeqRepo | 2021-01-29 |
| UTA | uta_20210129b |
| ClinVar | ~57.5M patterns (Dec 2024) |

> **Note**: Results may vary with different versions or ClinVar releases.

### Test Tiers

| Tier | Patterns | Tools | Approximate Runtime |
|------|----------|-------|---------------------|
| **Smoke** | 100 | All 4 | ~10 sec |
| **Small** | 1,000 | All 4 | ~1-2 min |
| **Medium** | 10,000 | All 4 | ~10-20 min |
| **Full ClinVar** | ~57.5M | ferro + mutalyzer + hgvs-rs | ~2-4 hours |

> **Note**: biocommons/hgvs doesn't scale to full ClinVar due to database query overhead (~27 days estimated).

### Normalization Capabilities by Tool

| Pattern Type | ferro | hgvs-rs | biocommons | mutalyzer |
|--------------|:-----:|:-------:|:----------:|:---------:|
| Genomic (g.) | ✓ | ✓ | ✓ | ✓ |
| Coding (c.) exonic | ✓ | ✓ | ✓ | ✓ |
| Coding (c.) intronic | ✓ | ✗ | ✗ | ✓* |
| Non-coding (n.) | ✓ | ✓ | ✓ | ✓ |
| RNA (r.) | ✓ | ✓ | ✓ | ✓ |
| Protein (p.) | ✓ | ✓ | ✗ | Net** |

\* mutalyzer intronic support enabled by default via genomic context rewriting.
Use `--no-rewrite-intronic` to disable.

\** mutalyzer protein normalization requires network access for NP_→NM_ lookups (see Design Limitations).

### Design Limitations

- **biocommons protein normalization**: Explicitly not supported by design. Use ferro or hgvs-rs instead.
- **hgvs-rs/biocommons intronic**: Neither tool normalizes intronic variants. Use ferro-hgvs or mutalyzer.
- **mutalyzer protein normalization**: Requires network access even with a local cache. Mutalyzer back-translates protein variants to DNA, which requires looking up NP_→NM_ transcript relationships from NCBI's API (`get_cds_to_mrna`). This lookup cannot be cached locally. Use ferro or hgvs-rs for offline protein normalization.

### Why ferro Has the Best Coverage

ferro-hgvs achieves comprehensive pattern support through its reference data preparation:

| Reference Data | Patterns Enabled | Without ferro |
|----------------|------------------|---------------|
| RefSeq transcripts (NM_, NR_, XM_) | All transcript-based variants | Network fetch per variant |
| GRCh38/GRCh37 genomes | NC_ genomic variants | Network fetch or failure |
| RefSeqGene (NG_) | Gene region variants | Network fetch or failure |
| LRG sequences | LRG_ stable references | Network fetch or failure |
| Protein sequences (NP_, XP_) | Protein normalization | Network fetch (mutalyzer) or failure (biocommons) |
| cdot transcript metadata | Intronic coordinate mapping | Not available (hgvs-rs/biocommons fail on intronic) |
| Legacy transcript versions | Historical ClinVar variants | Network fetch or failure |

When you run `ferro prepare`, all this reference data is downloaded and organized. Other tools
(mutalyzer, biocommons, hgvs-rs) can then use this cached data for local operation, dramatically
improving their performance from ~1 pattern/sec (network) to ~20 patterns/sec (local).

### Metrics Collected

At each comparison step, we collect:
- **Runtime**: Total time and patterns/second
- **Success rate**: Percentage of patterns successfully parsed/normalized
- **Cross-tool agreement**: How often tools produce the same result

---

## Quick Start

Minimal commands to run a small comparison (for experienced users).

> **⚠️ Quick Start Limitations**:
> - Only sets up ferro + mutalyzer (not biocommons/hgvs-rs)
> - For full ClinVar testing, see [Prepare All Tools](#prepare-all-tools)

```bash
# Clone (skip large LFS files for faster setup)
GIT_LFS_SKIP_SMUDGE=1 git clone <repo-url>
cd ferro-hgvs

# Build
cargo build --release --features benchmark,hgvs-rs

# 1. Prepare ferro reference data (main binary - no benchmark feature needed)
ferro prepare --output-dir data/ferro

# 2. Prepare mutalyzer cache (uses proteins from ferro)
ferro-benchmark prepare mutalyzer \
  --ferro-reference data/ferro \
  --output-dir data/mutalyzer \
  --shards 32

# 3. Create test patterns
echo "NM_000059.4:c.100del" > patterns.txt
echo "NM_001374258.1:c.1A>G" >> patterns.txt

# 4. Parse and normalize with ferro (for benchmark comparison, use ferro-benchmark)
ferro-benchmark parse ferro -i patterns.txt -o ferro_parse.json
ferro-benchmark normalize ferro -i patterns.txt -o ferro_norm.json --reference data/ferro

# 5. Parse and normalize with mutalyzer
ferro-benchmark parse mutalyzer -i patterns.txt -o mutalyzer_parse.json
ferro-benchmark normalize mutalyzer -i patterns.txt -o mutalyzer_norm.json \
  --mutalyzer-settings data/mutalyzer/mutalyzer_settings.conf -j 4

# 6. Compare results
ferro-benchmark compare results parse ferro_parse.json mutalyzer_parse.json -o parse_comparison.json
ferro-benchmark compare results normalize ferro_norm.json mutalyzer_norm.json -o norm_comparison.json
```

> **For biocommons/hgvs-rs**: These require UTA database setup. Download the UTA dump manually
> from https://dl.biocommons.org/uta/uta_20210129b.pgd.gz (requires human verification),
> then see [Prepare biocommons](#step-3-prepare-biocommons).

---

## Prerequisites & Build

### System Requirements

| Requirement | Minimum | Recommended |
|-------------|---------|-------------|
| **Disk Space** | 5GB | 50GB (full setup) |
| **RAM** | 4GB | 16GB |
| **CPU Cores** | 2 | 8+ |

### Software Dependencies

| Software | Version | Required For |
|----------|---------|--------------|
| **Rust** | 1.70+ | Building ferro-benchmark |
| **Python** | 3.10+ | mutalyzer, biocommons validators |
| **Docker** | Any | biocommons/hgvs-rs local UTA setup |
| **pixi** | Any | Python dependency management (recommended) |

### Installing Dependencies

```bash
# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ~/.cargo/env

# Install pixi (recommended for Python deps)
curl -fsSL https://pixi.sh/install.sh | bash

# Install Docker (macOS)
brew install --cask docker
# Start Docker Desktop from Applications
```

### NCBI API Key (Recommended)

The `prepare mutalyzer` command fetches sequences from NCBI. By default, NCBI limits unauthenticated requests to 3/second. With an API key, you get 10/second (3x faster).

```bash
# 1. Create an NCBI account at https://www.ncbi.nlm.nih.gov/account/
# 2. Go to Settings → API Key Management → Create API Key
# 3. Set the environment variable:
export NCBI_API_KEY="your-api-key-here"

# Add to shell profile for persistence:
echo 'export NCBI_API_KEY="your-api-key"' >> ~/.zshrc  # or ~/.bashrc
```

The `prepare mutalyzer --patterns` command will automatically use this key for faster fetching.

### Building ferro-benchmark

```bash
# Standard build (without hgvs-rs)
cargo build --release --features benchmark

# Full build with hgvs-rs support (recommended)
cargo build --release --features benchmark,hgvs-rs

# Verify build
./target/release/ferro-benchmark --help
```

### Using pixi Environment

```bash
pixi shell  # Activates Python environment with mutalyzer, biocommons
cargo build --release --features benchmark,hgvs-rs
```

---

## Download ClinVar Data

### Download and Extract

```bash
# Download ClinVar HGVS data
mkdir -p data/clinvar
curl -o data/clinvar/hgvs4variation.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/hgvs4variation.txt.gz

# Check size (~500MB compressed)
ls -lh data/clinvar/hgvs4variation.txt.gz

# Extract all patterns (~57.5M)
ferro-benchmark extract clinvar \
  -i data/clinvar/hgvs4variation.txt.gz \
  -o data/clinvar/clinvar_patterns.txt

# Check count
wc -l data/clinvar/clinvar_patterns.txt
```

### Create Sample Files

Sample files should be generated separately for parse and normalize tests:
- **Parse samples**: Include all pattern types (proteins can be parsed by all tools)
- **Normalize samples**: Exclude protein patterns (use `--exclude-protein`)

Protein patterns are excluded from normalize samples because protein normalization is handled
inconsistently across tools - ferro and hgvs-rs return proteins unchanged (no-op), biocommons
explicitly rejects them as "unsupported", and mutalyzer requires network access. Excluding
proteins ensures a fair comparison of actual normalization capabilities.

```bash
# Create sample directories
mkdir -p samples/parse samples/normalize

# ===== PARSE SAMPLES (include proteins) =====
# Smoke test (100 patterns)
ferro-benchmark extract sample -i data/clinvar/clinvar_patterns.txt \
  -o samples/parse/sample_100.txt --size 100 --seed 42

# Small (1K patterns)
ferro-benchmark extract sample -i data/clinvar/clinvar_patterns.txt \
  -o samples/parse/sample_1k.txt --size 1000 --seed 42

# Medium (10K patterns)
ferro-benchmark extract sample -i data/clinvar/clinvar_patterns.txt \
  -o samples/parse/sample_10k.txt --size 10000 --seed 42

# ===== NORMALIZE SAMPLES (exclude proteins) =====
# Smoke test (100 patterns, no proteins)
ferro-benchmark extract sample -i data/clinvar/clinvar_patterns.txt \
  -o samples/normalize/sample_100.txt --size 100 --seed 42 --exclude-protein

# Small (1K patterns, no proteins)
ferro-benchmark extract sample -i data/clinvar/clinvar_patterns.txt \
  -o samples/normalize/sample_1k.txt --size 1000 --seed 42 --exclude-protein

# Medium (10K patterns, no proteins)
ferro-benchmark extract sample -i data/clinvar/clinvar_patterns.txt \
  -o samples/normalize/sample_10k.txt --size 10000 --seed 42 --exclude-protein
```

> **Why exclude proteins from normalize samples?** Protein normalization is effectively a no-op
> (substitutions like `p.Glu483Lys` don't need 3' shifting). Different tools handle this
> inconsistently:
> - **ferro, hgvs-rs**: Accept and return unchanged → "success"
> - **biocommons**: Explicitly rejects as "unsupported" → "failure"
> - **mutalyzer**: Requires network for NP_→NM_ lookup
>
> Using `--exclude-protein` ensures normalize comparisons measure actual normalization capability.

> **Genomic patterns (NC_) and mutalyzer**: The mutalyzer cache is built for transcript-based
> normalization. Genomic patterns like `NC_000007.14:g.12345A>G` require either:
> - Full chromosome sequences in the cache (not included by default)
> - Network access (`--allow-network`) to fetch sequences from EBI
> - Rewriting to transcript context via `--rewrite-genomic` (if transcript mappings exist)
>
> For mutalyzer cache-only benchmarks, consider filtering NC_ patterns from normalize samples
> or using `--allow-network` and accepting slower performance.

---

## Network Fetch Behavior

Understanding what requires network access helps plan offline deployments and minimize external dependencies.

### Network Fetches by Command

| Command | Network Fetches | Notes |
|---------|-----------------|-------|
| `prepare ferro` | **Required** | Downloads RefSeq transcripts, cdot, genomes, RefSeqGene, LRG from NCBI/EBI |
| `prepare ferro --patterns` | **Required** | Additionally fetches legacy transcript versions detected in patterns |
| `prepare mutalyzer` | **Usually none** | Uses ferro reference; may fetch legacy versions not in patterns (see below) |
| `prepare biocommons` | **Required** | UTA Docker image + SeqRepo snapshot (~10GB each) |
| `prepare hgvs-rs` | **Usually none** | Shares UTA/SeqRepo with biocommons |

### How to Minimize Network Fetches

#### For Mutalyzer

The mutalyzer cache includes a **hardcoded list of 33 legacy transcript versions** (`LEGACY_TRANSCRIPT_VERSIONS` in `src/benchmark/cache.rs`). These are older RefSeq versions commonly found in ClinVar that have been superseded.

To ensure zero NCBI fetches during `prepare mutalyzer`:

```bash
# Option 1: Include all legacy versions in your patterns file
# Add these lines to your patterns file before running prepare ferro:
echo "NM_000051.3:c.1A>G" >> my_patterns.txt
echo "NM_000174.4:c.1A>G" >> my_patterns.txt
# ... (see full list in src/benchmark/cache.rs:942-977)

# Then prepare ferro with patterns and clinvar (derives proteins locally)
ferro-benchmark prepare ferro --output-dir data/ferro \
  --patterns my_patterns.txt \
  --clinvar data/clinvar/hgvs4variation.txt.gz

# Mutalyzer will use ferro's cached sequences and proteins
ferro-benchmark prepare mutalyzer --ferro-reference data/ferro --output-dir data/mutalyzer --shards 32
```

```bash
# Option 2: Accept that mutalyzer may fetch ~33 legacy versions
# This takes about 10 seconds with rate limiting
ferro-benchmark prepare mutalyzer --ferro-reference data/ferro --output-dir data/mutalyzer --shards 32
```

#### For Biocommons/hgvs-rs

These tools require one-time downloads that can be cached:

1. **UTA Database** (~10GB):
   - Download manually: https://dl.biocommons.org/uta/uta_20210129b.pgd.gz
   - Or use Docker: `docker pull biocommons/uta:uta_20210129b`

2. **SeqRepo Snapshot** (~10GB):
   - Downloads automatically during `prepare biocommons`
   - Cache at: `<seqrepo-dir>/2021-01-29/`

3. **Additional sequences from ferro** (automatic):
   ```bash
   # By default, all ferro transcripts (~270K) are loaded into SeqRepo
   ferro-benchmark prepare biocommons \
     --seqrepo-dir data/seqrepo \
     --ferro-reference data/ferro
   # Use --no-load-transcripts to skip transcript loading (only supplemental sequences)
   ```

### Legacy Transcript Versions

The `LEGACY_TRANSCRIPT_VERSIONS` constant is **manually maintained**. It was created from historical ClinVar test patterns and contains RefSeq versions that:
- Are referenced in ClinVar variant annotations
- Have been superseded by newer versions in the current cdot release
- Are not available in the standard RefSeq transcript downloads

**Current list (33 accessions):**
```
NM_000051.3, NM_000174.4, NM_000246.3, NM_000267.3, NM_000280.5,
NM_000426.3, NM_000449.3, NM_001127221.1, NM_001271223.2, NM_001458.4,
NM_002529.3, NM_002878.3, NM_004006.2, NM_004082.4, NM_004321.7,
NM_004364.4, NM_015125.4, NM_015243.2, NM_015247.2, NM_017890.4,
NM_018400.3, NM_018993.3, NM_019023.4, NM_021946.4, NM_030962.3,
NM_030984.5, NM_032790.3, NM_033517.2, NM_138387.3, NM_152383.4,
NM_172201.1, NM_177438.2, NM_203447.3
```

**To add new legacy versions:**
1. Edit `src/benchmark/cache.rs` and add to `LEGACY_TRANSCRIPT_VERSIONS`
2. Rebuild with `cargo build --release --features benchmark`
3. Re-run `prepare mutalyzer` to fetch the new versions

---

## Network vs Local Performance Comparison

ferro's reference preparation commands provide significant speedups by eliminating
network latency. Use the `--allow-network` flag to measure baseline performance
without local caching.

### Measuring the Speedup

#### Mutalyzer: Network vs Cache

```bash
# With ferro cache (fast, ~20 patterns/sec)
ferro-benchmark normalize mutalyzer -i sample_100.txt -o cached.json \
  --mutalyzer-settings data/mutalyzer/mutalyzer_settings.conf -j 1

# Without cache - network mode (slow, ~1 pattern/sec)
ferro-benchmark normalize mutalyzer -i sample_100.txt -o network.json \
  --mutalyzer-settings data/mutalyzer/mutalyzer_settings.conf --allow-network -j 1
```

#### Biocommons: Remote vs Local UTA

```bash
# With local UTA + SeqRepo (fast, ~20 patterns/sec)
ferro-benchmark normalize biocommons -i sample_100.txt -o local.json \
  --biocommons-settings data/ferro/biocommons_settings.conf

# With remote UTA (slow, ~0.2 patterns/sec)
# Simply omit --biocommons-settings to use remote uta.biocommons.org
ferro-benchmark normalize biocommons -i sample_100.txt -o remote.json
```

### Performance Impact Summary

| Tool | Network Mode | With ferro prep | Speedup |
|------|--------------|-----------------|---------|
| Mutalyzer | ~1 pattern/sec | ~20 patterns/sec | **20x** |
| biocommons/hgvs | ~0.2 patterns/sec | ~20 patterns/sec | **100x** |
| hgvs-rs | ~0.2 patterns/sec | ~2 patterns/sec | **10x** |

> **Note**: Network mode is rate-limited by NCBI (3-10 req/sec) and remote UTA server
> capacity. Actual speeds may vary based on network conditions and server load.

### When to Use Network Mode

- **Benchmarking**: Measure the actual speedup from ferro's preparation
- **Testing**: Verify that patterns work without local cache
- **Debugging**: Identify cache misses that cause network fallback

### Network Statistics

When using `--allow-network`, the output JSON includes network statistics:

```json
{
  "network_calls": 42,
  "network_stats": [
    {"url": "https://eutils.ncbi.nlm.nih.gov/...", "elapsed_seconds": 0.5, "success": true}
  ]
}
```

---

## UTA Alignment Loading

The UTA snapshot (2021-01-29) lacks transcript-to-genome alignments for newer transcripts.
ferro automatically loads missing alignments from cdot during `prepare hgvs-rs`.

### Why This Is Needed

hgvs-rs requires transcript alignments (exon structures) in UTA to normalize variants.
When a transcript is missing from UTA:

```
Error: replacing reference failed for NM_001407958.1
```

The cdot JSON file contains these alignments for newer transcripts. ferro bridges
the gap by loading them into UTA.

### How It Works

1. Loads cdot JSON to find all transcript alignments
2. Queries UTA for transcripts that already have alignments
3. Inserts missing transcripts and their exon structures
4. Creates exon alignments linking transcript to genomic coordinates

### Verifying Loaded Transcripts

```bash
# Check how many transcripts were loaded by ferro
docker exec ferro-uta psql -U anonymous -d uta -c \
  "SELECT COUNT(*) FROM uta_20210129b.transcript \
   WHERE origin_id = (SELECT origin_id FROM uta_20210129b.origin WHERE name = 'cdot-ferro');"

# Verify a specific transcript
docker exec ferro-uta psql -U anonymous -d uta -c \
  "SELECT tx_ac, alt_ac, alt_aln_method FROM uta_20210129b.exon_set \
   WHERE tx_ac = 'NM_001407958.1';"
```

### Disabling Alignment Loading

```bash
# Skip UTA alignment loading (not recommended)
ferro-benchmark prepare hgvs-rs \
  --seqrepo-dir data/seqrepo \
  --ferro-reference data/ferro \
  --no-load-alignments
```

> **Warning**: Without alignment loading, hgvs-rs will fail on transcripts
> added after the UTA snapshot date (January 2021).

---

## Recommended Workflow

For optimal performance with minimal network fetches, prepare ferro reference first, then use it for all other tool preparations.

### Optimization Summary

| Data Type | Source | Mutalyzer | Biocommons | hgvs-rs |
|-----------|--------|-----------|------------|---------|
| RefSeq transcripts | ferro | ✅ | ✅ | ✅ |
| Supplemental transcripts | ferro | ✅ | ✅ | ✅ |
| Legacy transcripts | ferro (if in patterns) | ✅ | ✅ | ✅ |
| Legacy GenBank | ferro | ✅ | ✅ | ✅ |
| Protein sequences | ferro (derived from CDS) | ✅ | N/A | N/A |
| LRG sequences | ferro | ✅ | ✅ | ✅ |
| LRG annotations | ferro XML | ✅ | N/A | N/A |
| Genomic (NC_) | ferro | Skipped | ✅ | ✅ |
| Transcript alignments | cdot | N/A | N/A | ✅ (UTA) |

> **Note**: Mutalyzer cannot normalize protein (p.) patterns - see [Design Limitations](#design-limitations). However, it uses protein sequences for reference validation.

### Expected Results

With a fully prepared ferro reference:
- **Zero NCBI nucleotide fetches** in mutalyzer prepare (if all legacy versions are in patterns)
- **Zero protein fetches** (derived from CDS translation)
- **Zero EBI fetches** for LRG data
- **Biocommons/hgvs-rs**: Zero NCBI fetches after initial SeqRepo/UTA setup

### What Gets Cached in Ferro

The `prepare ferro` command downloads and caches:
- **RefSeq transcripts**: Current transcript sequences from NCBI (~1GB)
- **cdot metadata**: Transcript-to-genome mappings (~200MB)
- **Genomes**: GRCh37 and GRCh38 chromosome sequences (~4GB)
- **RefSeqGene**: Gene-specific reference sequences (~600MB)
- **LRG sequences**: Locus Reference Genomic sequences (~50MB)
- **Legacy transcript versions** (with --patterns): Older RefSeq versions detected in patterns
- **Legacy GenBank sequences**: Historical non-RefSeq accessions (U31929.1, M75126.1, etc.)
- **Protein sequences** (with --clinvar): Derived from transcript CDS translation (~200MB)

---

## Prepare All Tools

### Quick: Prepare All at Once

The `prepare all` command sets up all four tools in sequence with sensible defaults:

```bash
# Prepare all tools (ferro → mutalyzer → biocommons → hgvs-rs)
# Requires Docker for biocommons/hgvs-rs UTA setup
ferro-benchmark prepare all \
  --output-dir data/benchmark \
  --seqrepo-dir data/seqrepo \
  --uta-dump path/to/uta_20210129b.pgd.gz \
  --patterns data/clinvar/clinvar_patterns.txt \
  --clinvar data/clinvar/hgvs4variation.txt.gz
```

This creates:
- `data/benchmark/ferro/` - Ferro reference data
- `data/benchmark/mutalyzer/` - Mutalyzer cache
- `data/benchmark/biocommons_settings.conf` - Biocommons settings
- `data/benchmark/hgvs-rs_settings.conf` - hgvs-rs settings
- `data/seqrepo/` - Shared SeqRepo for biocommons and hgvs-rs

After completion, verify with:

```bash
ferro-benchmark check all \
  --reference data/benchmark \
  --seqrepo-path data/seqrepo/2021-01-29
```

### Step 1: Prepare ferro

```bash
# Basic preparation (available in main ferro binary - no benchmark feature needed)
ferro prepare --output-dir data/ferro

# Full setup for ClinVar testing with supplemental accessions (~7GB total)
# Store patterns and clinvar for auto-detection by other tools
ferro-benchmark prepare ferro --output-dir data/ferro \
  --patterns data/clinvar/clinvar_patterns.txt \
  --clinvar data/clinvar/hgvs4variation.txt.gz
```

> **Note**: By default, `prepare` downloads all reference data:
> - GRCh38 + GRCh37 genomes (~4GB total)
> - RefSeqGene sequences (~600MB)
> - LRG sequences (~50MB)
>
> Use `--genome none` to skip genomes, `--no-refseqgene` to skip RefSeqGene, or `--no-lrg` to skip LRG.

What each option does:
- `--genome all` (default) - Both GRCh38 and GRCh37 genomes
- `--genome grch38` - GRCh38 only (~3GB)
- `--genome grch37` - GRCh37 only (~900MB)
- `--genome none` - Skip genome downloads
- `--no-refseqgene` - Skip RefSeqGene sequences
- `--no-lrg` - Skip LRG sequences
- `--patterns <file>` - Store patterns file for provenance and auto-detection
- `--clinvar <file>` - Store ClinVar reference and derive protein sequences from transcript CDS

> **Protein Derivation**: When `--clinvar` is provided, ferro automatically derives protein sequences
> from transcript CDS coordinates. For each NM_/XM_ transcript with CDS metadata in cdot,
> ferro translates the coding sequence to produce the corresponding NP_/XP_ protein sequence.
> This eliminates the need to fetch protein sequences from NCBI.

> **Auto-Detection**: When you use `--patterns` and `--clinvar`, these files are stored in the ferro reference:
> - `data/ferro/patterns/` - patterns file (copied)
> - `data/ferro/clinvar/` - ClinVar file (symlinked to avoid duplication)
>
> Subsequent `prepare mutalyzer`, `prepare biocommons`, and `prepare hgvs-rs` commands
> will automatically detect and use these files when you specify `--ferro-reference data/ferro`.

### Step 2: Prepare mutalyzer

> **Note**: The mutalyzer cache is built from ferro reference data. Prepare ferro first.

```bash
# Full setup for ClinVar testing (recommended)
# Patterns and ClinVar are auto-detected from ferro reference
ferro-benchmark prepare mutalyzer \
  --ferro-reference data/ferro \
  --output-dir data/mutalyzer \
  --shards 32
```

What each flag does:
- `--ferro-reference` - Ferro reference directory (patterns, clinvar, and proteins auto-detected)
- `--shards` - Number of logical shards for parallel workers (default: 32)
- `--no-lrg` - Skip LRG annotation enhancement (LRG is enhanced by default)
- `--no-proteins` - Skip copying proteins from ferro reference
- `--patterns` - Override auto-detected patterns with a specific file

> **Sharding for Memory Efficiency**: The `--shards` option enables memory-efficient parallel
> normalization. With 32 shards and 8 workers, each worker processes patterns from ~4 shards,
> loading only ~1/8 of accessions into memory. Mutalyzer loads sequences on-demand, so
> pattern routing (not cache partitioning) provides the memory savings.

> **Protein Sequences**: Mutalyzer uses protein sequences from the ferro reference. These are
> automatically derived from transcript CDS during `prepare ferro --clinvar`. If ferro has
> derived proteins, they are copied to the mutalyzer cache automatically.
>
> Use `--no-proteins` to skip this step (protein patterns will fail during normalization).

The prepare command automatically:
1. Populates the cache from ferro reference data
2. Auto-detects patterns from `{ferro-reference}/patterns/` directory
3. Enhances NC_ annotations with transcript mappings for intronic variants
4. Enhances LRG annotations (unless `--no-lrg`)
5. Pre-fetches accessions from patterns file
6. Copies protein sequences from ferro reference (unless `--no-proteins`)
7. Creates shard manifest for parallel normalization (if `--shards` > 1)
8. Creates `data/mutalyzer/mutalyzer_settings.conf`

### Step 3: Prepare biocommons

```bash
# Download UTA dump first (requires browser due to human verification):
# https://dl.biocommons.org/uta/uta_20210129b.pgd.gz

# Full setup (patterns auto-detected from ferro reference)
ferro-benchmark prepare biocommons --seqrepo-dir data/seqrepo \
  --uta-dump path/to/uta_20210129b.pgd.gz \
  --ferro-reference data/ferro

# If UTA is already running (e.g., from previous setup), omit --uta-dump
ferro-benchmark prepare biocommons --seqrepo-dir data/seqrepo \
  --ferro-reference data/ferro
```

What each option does:
- `--seqrepo-dir` - Where to store/find SeqRepo data (required)
- `--uta-dump` - Path to UTA dump file; auto-sets up UTA if not running
- `--ferro-reference` - Load sequences from ferro reference; auto-detects patterns from `patterns/` subdirectory
- `--patterns <file>` - Override auto-detected patterns with a specific file

This creates `data/ferro/biocommons_settings.conf` with UTA and SeqRepo paths.

> **Important**: The settings file contains **absolute paths**. If you move the reference data
> to a different location, you must either:
> 1. Edit the settings file to update `HGVS_SEQREPO_DIR` to the new location, OR
> 2. Use command-line overrides: `--seqrepo-path /new/path/to/seqrepo/2021-01-29`

> **Note**: The `--ferro-reference` option loads additional sequences into SeqRepo from the ferro reference data:
> - **NC_** (chromosomes): GRCh38 and GRCh37 genomes
> - **NG_** (RefSeqGene): Gene-specific reference sequences
> - **LRG_** (Locus Reference Genomic): Stable genomic references
>
> **LRG Support**:
> - **Genomic patterns** (`LRG_1:g.5000A>T`): Work automatically via SeqRepo aliases
> - **Transcript patterns** (`LRG_199t1:c.79del`): Use `--lrg-mapping` to translate to RefSeq:
>   ```bash
>   ferro-benchmark normalize biocommons -i patterns.txt -o results.json \
>     --biocommons-settings data/ferro/biocommons_settings.conf \
>     --lrg-mapping data/ferro/lrg/lrg_refseq_mapping.txt
>   ```
>   This translates `LRG_199t1` → `NM_004006.2` before normalization.
>
> Without these sequences, biocommons will attempt slow network fetches to NCBI/EBI for genomic patterns.

> **macOS Note**: The SeqRepo download uses `seqrepo pull` internally, which requires real `rsync` (not macOS's `openrsync`). If you encounter rsync errors, use pixi's rsync directly:
> ```bash
> pixi run seqrepo --rsync-exe $(pixi run which rsync) \
>   --root-directory data/seqrepo pull --instance-name 2021-01-29
> ```

### Step 4: Prepare hgvs-rs

hgvs-rs shares UTA and SeqRepo with biocommons:

```bash
# Full setup with cdot alignment loading (recommended)
# Patterns auto-detected from ferro reference, alignments loaded into UTA
ferro-benchmark prepare hgvs-rs --seqrepo-dir data/seqrepo \
  --ferro-reference data/ferro
```

> **Note**: hgvs-rs uses the same SeqRepo and UTA as biocommons. The `--ferro-reference` option:
> - Auto-detects patterns from `patterns/` subdirectory
> - Loads transcript alignments from cdot into UTA for transcripts missing from the snapshot (2021-01-29)

#### Transcript Alignment Loading

By default, `prepare hgvs-rs` loads transcript alignments from cdot into UTA for
transcripts missing from the UTA snapshot (2021-01-29). This enables hgvs-rs to
normalize variants for newer transcripts (e.g., NM_001407958.1).

```bash
# Default behavior: load cdot alignments into UTA
ferro-benchmark prepare hgvs-rs --seqrepo-dir data/seqrepo \
  --ferro-reference data/ferro

# Skip alignment loading (not recommended)
ferro-benchmark prepare hgvs-rs --seqrepo-dir data/seqrepo \
  --ferro-reference data/ferro --no-load-alignments
```

> **Warning**: Without alignment loading, hgvs-rs will fail on transcripts
> added after the UTA snapshot date (January 2021) with "replacing reference failed" errors.

---

## Check All Tools

After preparing all tools, verify everything is ready:

```bash
# Check all tools at once (recommended after 'prepare all')
ferro-benchmark check all \
  --reference data/benchmark \
  --mutalyzer-settings data/benchmark/mutalyzer/mutalyzer_settings.conf \
  --seqrepo-path data/seqrepo/2021-01-29

# Or check individual tools:
ferro-benchmark check ferro --reference data/ferro
ferro-benchmark check mutalyzer --mutalyzer-settings data/mutalyzer/mutalyzer_settings.conf
ferro-benchmark check biocommons --seqrepo-path data/seqrepo/2021-01-29
ferro-benchmark check hgvs-rs --seqrepo-path data/seqrepo/2021-01-29
```

### Expected Output

All tools should show "ready":

| Tool | Status | Notes |
|------|--------|-------|
| ferro | ✓ Ready | Transcripts, genome, refseqgene, lrg loaded |
| mutalyzer | ✓ Ready | Cache populated with ~290K entries |
| biocommons | ✓ Ready | UTA + SeqRepo connected |
| hgvs-rs | ✓ Ready | Provider initialization successful |

---

## Parse Comparison: Small & Medium

Run parsing comparisons with all 4 tools at increasing scale.

### Smoke Test (100 patterns)

```bash
SAMPLE=samples/parse/sample_100.txt
OUTDIR=results/parse/smoke
mkdir -p $OUTDIR

# Parse with all tools
ferro-benchmark parse ferro -i $SAMPLE -o $OUTDIR/ferro.json
ferro-benchmark parse mutalyzer -i $SAMPLE -o $OUTDIR/mutalyzer.json
ferro-benchmark parse biocommons -i $SAMPLE -o $OUTDIR/biocommons.json
ferro-benchmark parse hgvs-rs -i $SAMPLE -o $OUTDIR/hgvs_rs.json

# Compare all
ferro-benchmark compare results parse \
  $OUTDIR/ferro.json $OUTDIR/mutalyzer.json $OUTDIR/biocommons.json $OUTDIR/hgvs_rs.json \
  -o $OUTDIR/comparison.json
```

### Small (1K patterns)

```bash
SAMPLE=samples/parse/sample_1k.txt
OUTDIR=results/parse/small
mkdir -p $OUTDIR

ferro-benchmark parse ferro -i $SAMPLE -o $OUTDIR/ferro.json
ferro-benchmark parse mutalyzer -i $SAMPLE -o $OUTDIR/mutalyzer.json
ferro-benchmark parse biocommons -i $SAMPLE -o $OUTDIR/biocommons.json
ferro-benchmark parse hgvs-rs -i $SAMPLE -o $OUTDIR/hgvs_rs.json

ferro-benchmark compare results parse \
  $OUTDIR/ferro.json $OUTDIR/mutalyzer.json $OUTDIR/biocommons.json $OUTDIR/hgvs_rs.json \
  -o $OUTDIR/comparison.json
```

### Medium (10K patterns)

```bash
SAMPLE=samples/parse/sample_10k.txt
OUTDIR=results/parse/medium
mkdir -p $OUTDIR

ferro-benchmark parse ferro -i $SAMPLE -o $OUTDIR/ferro.json
ferro-benchmark parse mutalyzer -i $SAMPLE -o $OUTDIR/mutalyzer.json
ferro-benchmark parse biocommons -i $SAMPLE -o $OUTDIR/biocommons.json
ferro-benchmark parse hgvs-rs -i $SAMPLE -o $OUTDIR/hgvs_rs.json

ferro-benchmark compare results parse \
  $OUTDIR/ferro.json $OUTDIR/mutalyzer.json $OUTDIR/biocommons.json $OUTDIR/hgvs_rs.json \
  -o $OUTDIR/comparison.json
```

### Example Parse Results

Results from ClinVar patterns sampled with `--seed 42`:

#### Smoke Test (100 patterns)

| Tool | Success | Failed | Time | Patterns/sec |
|------|---------|--------|------|--------------|
| ferro | 100 (100%) | 0 | 0.0003s | ~390,000 |
| mutalyzer | 98 (98%) | 2 | 0.39s | ~250 |
| biocommons | 97 (97%) | 3 | 0.04s | ~2,700 |
| hgvs-rs | 97 (97%) | 3 | 0.0002s | ~490,000 |

#### Small (1K patterns)

| Tool | Success | Failed | Time | Patterns/sec |
|------|---------|--------|------|--------------|
| ferro | 1000 (100%) | 0 | 0.0008s | ~1,250,000 |
| mutalyzer | 990 (99%) | 10 | 3.4s | ~290 |
| biocommons | 991 (99.1%) | 9 | 0.80s | ~1,240 |
| hgvs-rs | 991 (99.1%) | 9 | 0.0008s | ~1,240,000 |

#### Medium (10K patterns)

| Tool | Success | Failed | Time | Patterns/sec |
|------|---------|--------|------|--------------|
| ferro | 10000 (100%) | 0 | 0.006s | ~1,630,000 |
| mutalyzer | 9936 (99.4%) | 64 | 33.6s | ~296 |
| biocommons | 9939 (99.4%) | 61 | 4.5s | ~2,200 |
| hgvs-rs | 9939 (99.4%) | 61 | 0.007s | ~1,370,000 |

### Why Ferro Has 100% Parse Success

Ferro parses patterns that other tools reject:

| Pattern Type | Example | Ferro | Others |
|--------------|---------|-------|--------|
| **Frameshift without terminator** | `p.Met847Glyfs` | ✓ Parses | ✗ Expects `*N` suffix |
| **Repeat notation** | `c.-963GCTACTGCTGCTGCT[3]` | ✓ Parses | ✗ Expects `=` |
| **Uncertain range deletions** | `g.(?_36993051)_(37050663_?)del` | ✓ Parses | ✗ Not supported |

### Failure Categories by Tool (10K sample)

| Tool | Failure Category | Count | Example |
|------|------------------|-------|---------|
| **mutalyzer** | Frameshift without terminator | 63 | `p.Pro1614Hisfs` |
| **mutalyzer** | Complex allele notation | 1 | `NM_006876.2:[c.1168A>G;c.1217C>T]` |
| **biocommons/hgvs-rs** | Repeat notation | 45 | `NC_000010.10:g.76788505AGG[2]` |
| **biocommons/hgvs-rs** | Uncertain range | 13 | `NC_000007.13:g.(?_26240172)_(26240197_?)del` |
| **biocommons/hgvs-rs** | Complex allele notation | 2 | `NM_033028.4:r.[118_261del;118_373del]` |

> **Note**: These are valid HGVS patterns found in ClinVar. Ferro's parser is more permissive to handle real-world variant annotations.

### Cross-Tool Output Agreement

When comparing the actual parse outputs across tools (patterns successful in all 4 tools):

| Sample | All Succeed | All Agree | Agreement Rate |
|--------|-------------|-----------|----------------|
| Smoke (100) | 96 | 93 | 96.9% |
| Small (1K) | 982 | 963 | 98.1% |
| Medium (10K) | 9,876 | 9,669 | 97.9% |

#### Output Difference Categories (10K sample)

| Category | Count | Example |
|----------|-------|---------|
| **Frameshift Ter suffix** | 140 | `p.Asn1070fs` (ferro/mutalyzer) vs `p.Asn1070fsTer` (biocommons/hgvs-rs) |
| **Deleted sequence stripped** | 50 | `c.4849_4851delAAG` (ferro/mutalyzer/hgvs-rs) vs `c.4849_4851del` (biocommons) |
| **Delins formatting** | 7 | `c.145_146delinsAA` (ferro/biocommons) vs `c.145_146delGCinsAA` (mutalyzer/hgvs-rs) |
| **Dup sequence stripped** | 7 | `c.3187dupT` (ferro/mutalyzer/hgvs-rs) vs `c.3187dup` (biocommons) |
| **Other** | 3 | Various minor formatting differences |

> **Key Insight**: Most differences are stylistic variations in HGVS representation:
> - **Ter suffix**: HGVS allows `p.Xfs` or `p.XfsTer` for frameshifts
> - **Explicit sequences**: Deleted/duplicated sequences can be included or omitted
> - **Star notation**: `*` and `Ter` are both valid for termination codons
>
> All outputs represent the same biological variant - they differ only in formatting preference.

---

## Normalize Comparison: Small & Medium

Run normalization comparisons with all 4 tools at increasing scale.

> **Parallel Workers**: Use `-j 4` for mutalyzer/biocommons to speed up normalization. Each Python worker uses ~12-14GB RAM, so 4 workers requires ~50GB RAM. Adjust based on available memory.

> **LRG Transcript Support**: ClinVar contains many LRG transcript patterns (e.g., `LRG_199t1:c.79del`).
> UTA doesn't have LRG transcript definitions, so biocommons/hgvs-rs will fail on these patterns.
> Use `--lrg-mapping data/ferro/lrg/lrg_refseq_mapping.txt` to translate LRG transcripts to RefSeq equivalents.

> **UTA Must Be Running**: biocommons and hgvs-rs require the UTA PostgreSQL database. If you've restarted
> your machine since running `prepare biocommons`, start UTA first:
> ```bash
> ferro-benchmark setup start-uta
> # Verify it's running:
> docker ps -f name=ferro-uta
> ```

### Smoke Test (100 patterns)

> **Note**: Normalize samples exclude protein patterns. See [Create Sample Files](#create-sample-files) for details.

```bash
SAMPLE=samples/normalize/sample_100.txt
OUTDIR=results/normalize/smoke
mkdir -p $OUTDIR

ferro-benchmark normalize ferro -i $SAMPLE -o $OUTDIR/ferro.json \
  --reference data/ferro
ferro-benchmark normalize mutalyzer -i $SAMPLE -o $OUTDIR/mutalyzer.json \
  --mutalyzer-settings data/mutalyzer/mutalyzer_settings.conf -j 4
ferro-benchmark normalize biocommons -i $SAMPLE -o $OUTDIR/biocommons.json \
  --biocommons-settings data/ferro/biocommons_settings.conf \
  --lrg-mapping data/ferro/lrg/lrg_refseq_mapping.txt
ferro-benchmark normalize hgvs-rs -i $SAMPLE -o $OUTDIR/hgvs_rs.json \
  --seqrepo-path data/seqrepo/2021-01-29 \
  --lrg-mapping data/ferro/lrg/lrg_refseq_mapping.txt

ferro-benchmark compare results normalize \
  $OUTDIR/ferro.json $OUTDIR/mutalyzer.json $OUTDIR/biocommons.json $OUTDIR/hgvs_rs.json \
  -o $OUTDIR/comparison.json
```

### Small (1K patterns)

```bash
SAMPLE=samples/normalize/sample_1k.txt
OUTDIR=results/normalize/small
mkdir -p $OUTDIR

ferro-benchmark normalize ferro -i $SAMPLE -o $OUTDIR/ferro.json \
  --reference data/ferro
ferro-benchmark normalize mutalyzer -i $SAMPLE -o $OUTDIR/mutalyzer.json \
  --mutalyzer-settings data/mutalyzer/mutalyzer_settings.conf -j 4
ferro-benchmark normalize biocommons -i $SAMPLE -o $OUTDIR/biocommons.json \
  --biocommons-settings data/ferro/biocommons_settings.conf \
  --lrg-mapping data/ferro/lrg/lrg_refseq_mapping.txt
ferro-benchmark normalize hgvs-rs -i $SAMPLE -o $OUTDIR/hgvs_rs.json \
  --seqrepo-path data/seqrepo/2021-01-29 \
  --lrg-mapping data/ferro/lrg/lrg_refseq_mapping.txt

ferro-benchmark compare results normalize \
  $OUTDIR/ferro.json $OUTDIR/mutalyzer.json $OUTDIR/biocommons.json $OUTDIR/hgvs_rs.json \
  -o $OUTDIR/comparison.json
```

### Medium (10K patterns)

```bash
SAMPLE=samples/normalize/sample_10k.txt
OUTDIR=results/normalize/medium
mkdir -p $OUTDIR

ferro-benchmark normalize ferro -i $SAMPLE -o $OUTDIR/ferro.json \
  --reference data/ferro
ferro-benchmark normalize mutalyzer -i $SAMPLE -o $OUTDIR/mutalyzer.json \
  --mutalyzer-settings data/mutalyzer/mutalyzer_settings.conf -j 4
ferro-benchmark normalize biocommons -i $SAMPLE -o $OUTDIR/biocommons.json \
  --biocommons-settings data/ferro/biocommons_settings.conf \
  --lrg-mapping data/ferro/lrg/lrg_refseq_mapping.txt --workers 4
ferro-benchmark normalize hgvs-rs -i $SAMPLE -o $OUTDIR/hgvs_rs.json \
  --seqrepo-path data/seqrepo/2021-01-29 \
  --lrg-mapping data/ferro/lrg/lrg_refseq_mapping.txt --workers 8

ferro-benchmark compare results normalize \
  $OUTDIR/ferro.json $OUTDIR/mutalyzer.json $OUTDIR/biocommons.json $OUTDIR/hgvs_rs.json \
  -o $OUTDIR/comparison.json
```

### Expected Normalize Metrics

| Tool | 100 patterns | 1K patterns | 10K patterns |
|------|--------------|-------------|--------------|
| ferro | ~0.01s | ~0.1s | ~1s |
| mutalyzer | ~5s | ~50s | ~8 min |
| biocommons | ~5s | ~50s | ~8 min |
| hgvs-rs | ~5s | ~50s | ~8 min |

---

## Parse Comparison: Full ClinVar

> **Important**: biocommons is excluded from full ClinVar comparisons due to speed limitations (~27 days estimated for 57.5M patterns).

```bash
INPUT=data/clinvar/clinvar_patterns.txt
OUTDIR=results/parse/full
mkdir -p $OUTDIR

# Parse with ferro (fastest)
echo "=== Parsing with ferro ==="
ferro-benchmark parse ferro -i $INPUT -o $OUTDIR/ferro.json

# Parse with mutalyzer
echo "=== Parsing with mutalyzer ==="
ferro-benchmark parse mutalyzer -i $INPUT -o $OUTDIR/mutalyzer.json

# Parse with hgvs-rs
echo "=== Parsing with hgvs-rs ==="
ferro-benchmark parse hgvs-rs -i $INPUT -o $OUTDIR/hgvs_rs.json

# Compare (3 tools)
ferro-benchmark compare results parse \
  $OUTDIR/ferro.json $OUTDIR/mutalyzer.json $OUTDIR/hgvs_rs.json \
  -o $OUTDIR/comparison.json
```

### Expected Full ClinVar Parse Metrics

| Tool | ~57.5M patterns | Patterns/sec |
|------|-----------------|--------------|
| ferro | ~12 sec | ~4,000,000 |
| mutalyzer | ~40 min | ~20,000 |
| hgvs-rs | ~4 min | ~200,000 |
| biocommons | ~27 days | ~20 (excluded) |

ferro parses 57,523,009 of 57,524,053 patterns successfully (99.998%).
The 1,044 parse failures are malformed syntax not conforming to any HGVS grammar
(e.g., multi-range uncertain positions, multi-repeat notation, non-standard suffixes).

---

## Normalize Comparison: Full ClinVar

> **Important**: biocommons is excluded from full ClinVar comparisons due to speed limitations.

```bash
INPUT=data/clinvar/clinvar_patterns.txt
OUTDIR=results/normalize/full
mkdir -p $OUTDIR

# Normalize with ferro (fastest)
echo "=== Normalizing with ferro ==="
ferro-benchmark normalize ferro -i $INPUT -o $OUTDIR/ferro.json \
  --reference data/ferro

# Normalize with mutalyzer (parallel workers)
echo "=== Normalizing with mutalyzer ==="
ferro-benchmark normalize mutalyzer -i $INPUT -o $OUTDIR/mutalyzer.json \
  --mutalyzer-settings data/mutalyzer/mutalyzer_settings.conf -j 24

# Normalize with hgvs-rs
echo "=== Normalizing with hgvs-rs ==="
ferro-benchmark normalize hgvs-rs -i $INPUT -o $OUTDIR/hgvs_rs.json \
  --seqrepo-path data/seqrepo/2021-01-29 --workers 24

# Compare (3 tools)
ferro-benchmark compare results normalize \
  $OUTDIR/ferro.json $OUTDIR/mutalyzer.json $OUTDIR/hgvs_rs.json \
  -o $OUTDIR/comparison.json
```

### Expected Full ClinVar Normalize Metrics

| Tool | Time | Patterns/sec |
|------|------|--------------|
| ferro | ~12 sec | ~4,000,000 |
| mutalyzer | ~40 min | ~20,000 |
| hgvs-rs | ~4 hours | ~200 |
| biocommons | ~27 days | ~20 (excluded) |

### ferro Full ClinVar Normalize Pass Rate

On 57,524,053 ClinVar patterns with a fully prepared ferro reference:

| Metric | Count | Rate |
|--------|-------|------|
| **Total patterns** | 57,524,053 | — |
| **Parse failures** | 1,044 | 0.002% |
| **Normalize success** | 57,522,070 | **99.997%** |
| **Normalize failures** | 1,983 | 0.003% |

The 1,983 normalize failures break down as:

| Category | Count | Description |
|----------|-------|-------------|
| Parse errors | 1,044 | Malformed HGVS syntax (inherently unparseable) |
| Intronic conversion | 909 | CDS-to-tx position doesn't match intron boundary (cdot coordinate discrepancy) |
| Boundary-spanning n. | 30 | n. variant spans exon-intron boundary (not yet supported) |

> **Note**: Excluding the 1,044 patterns that cannot be parsed by any tool,
> ferro normalizes 99.998% of all parseable patterns (57,522,070 / 57,523,009).

---

## Troubleshooting

### Clone is slow

LFS objects are large. Skip them for faster setup:

```bash
GIT_LFS_SKIP_SMUDGE=1 git clone <repo-url>
```

### Docker not found

```bash
docker --version
docker ps

# macOS: Start Docker Desktop from Applications
```

### SeqRepo rsync error on macOS

macOS ships with `openrsync`, but SeqRepo requires real `rsync`:

```bash
# Error: "Binary located at /usr/bin/rsync appears to be an `openrsync` instance"

# Solution: Use pixi's rsync
pixi run seqrepo --rsync-exe $(pixi run which rsync) \
  --root-directory data/seqrepo pull --instance-name 2021-01-29
```

### UTA connection failed

```bash
# Check container status
docker ps -f name=ferro-uta

# Start if stopped
ferro-benchmark setup start-uta
# or
docker start ferro-uta

# Verify connection
docker exec ferro-uta psql -U anonymous -d uta -c "SELECT 1"
```

### Mutalyzer cache miss errors

```bash
# Check cache coverage
ls data/mutalyzer/ | wc -l

# Re-run prepare with patterns to fetch missing accessions
ferro-benchmark prepare mutalyzer \
  --ferro-reference data/ferro \
  --output-dir data/mutalyzer \
  --shards 32 \
  --patterns your_patterns.txt
```

### SeqRepo not found

```bash
# Verify directory structure
ls data/seqrepo/2021-01-29/aliases.sqlite3
ls data/seqrepo/2021-01-29/sequences/
```

### Settings file has wrong SeqRepo path

```
Error: Unable to open SeqRepo directory /old/path/to/seqrepo/2021-01-29
```

The `biocommons_settings.conf` file contains absolute paths that may be outdated
if you moved the reference data.

**Solution 1**: Update the settings file:
```bash
# Check current path
cat data/ferro/biocommons_settings.conf

# Edit HGVS_SEQREPO_DIR to the correct path
```

**Solution 2**: Use command-line override:
```bash
ferro-benchmark normalize biocommons -i patterns.txt -o results.json \
  --biocommons-settings data/ferro/biocommons_settings.conf \
  --seqrepo-path /correct/path/to/seqrepo/2021-01-29
```

### Python validators not found

```bash
# Activate pixi environment
pixi shell

# Or install manually
pip install mutalyzer-hgvs-parser hgvs biocommons.seqrepo
```

### Build fails with hgvs-rs feature

```bash
# Ensure PostgreSQL client is available
# macOS
brew install postgresql

# Linux
apt install libpq-dev
```

---

## CLI Command Reference

### Core Commands

| Command | Description |
|---------|-------------|
| `prepare <tool>` | Download/prepare reference data for a tool |
| `check <tool>` | Verify tool configuration and dependencies |
| `parse <tool>` | Parse HGVS patterns with a tool |
| `normalize <tool>` | Normalize HGVS patterns with a tool |

Where `<tool>` is one of: `ferro`, `mutalyzer`, `biocommons`, `hgvs-rs`, or `all`

> **Note**: `prepare all` and `check all` set up and verify all four tools in sequence. `parse all` and `normalize all` are not supported - use individual tools instead.

### Compare Subcommands

| Command | Description |
|---------|-------------|
| `compare results parse` | Compare parsing results between tools |
| `compare results normalize` | Compare normalization results between tools |
| `compare arbitration` | View arbitration decisions |

### Extract Subcommands

| Command | Description |
|---------|-------------|
| `extract clinvar` | Extract patterns from ClinVar TSV |
| `extract json` | Extract patterns from JSON files |
| `extract proteins` | Extract protein accessions from ClinVar |
| `extract sample` | Create stratified random sample (use `--exclude-protein` for normalize) |
| `extract shard` | Split dataset for parallel processing |

### Setup Subcommands

| Command | Description |
|---------|-------------|
| `setup uta` | Set up UTA PostgreSQL database via Docker |
| `setup seqrepo` | Download SeqRepo data |
| `setup start-uta` | Start UTA Docker container |
| `setup stop-uta` | Stop UTA Docker container |

> **Note**: UTA and SeqRepo setup are also integrated into `prepare biocommons` and `prepare hgvs-rs` via the `--uta-dump` flag.

### Generate Subcommands

| Command | Description |
|---------|-------------|
| `generate summary` | Create aggregate statistics |
| `generate report` | Generate markdown report |
| `generate readme` | Generate README-ready tables |

### Collate Subcommands

| Command | Description |
|---------|-------------|
| `collate parsing` | Collate parsing results from shards |
| `collate normalization` | Collate normalization results from shards |
