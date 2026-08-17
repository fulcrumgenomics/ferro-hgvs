# Why ferro-hgvs?

ferro-hgvs provides the most comprehensive HGVS variant normalization across all pattern types, with performance orders of magnitude faster than alternatives.

## Normalization Capabilities Comparison

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

## Performance Comparison

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

## Reference Data: What ferro Prepares

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

### Deriving version-independent NG_ placements (#728)

`ferro prepare --derive-ng-placements <accessions.txt>` derives genomic placements for the listed `NG_` versions (one exact accession per line, e.g. `NG_012337.3`; blank lines and `#` comments ignored), writing `derived_refseqgene_placements.json` into the reference directory and wiring the manifest's `derived_refseqgene_placements` field. This fills version gaps the archived RefSeqGene→genome GFF3 snapshots do not cover. It needs cdot + the genome in the same prepare run and uses NCBI EFetch per accession; accessions that cannot be validated are skipped with a warning. The field is preserved across subsequent `prepare` runs.
