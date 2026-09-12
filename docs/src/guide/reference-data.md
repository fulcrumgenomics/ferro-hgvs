# Reference data

Parsing needs no reference data. Normalization does — it must read the transcript and genome
sequences a description refers to. `ferro prepare` downloads and assembles them into a reference
directory.

## Prepare a reference

```bash
ferro prepare --output-dir ferro-reference
```

This downloads RefSeq transcripts, genome FASTAs, and cdot metadata, and writes them under
`ferro-reference/`. A bare `prepare` builds a **RefSeq-only** reference — accessions `NM_` / `NR_` /
`NP_` / `NG_`.

## Verify it

```bash
ferro check --reference ferro-reference
```

Optionally pre-build the on-disk cdot cache as a setup step, so the one-time cache build does not
slow the start of a real (or timed) run:

```bash
ferro check --reference ferro-reference --build-cache
```

## Optional data

Two opt-in flags provision more than RefSeq. Pass them at prepare time; both are incremental, so
re-running `prepare` over an existing reference adds the requested data and preserves what is already
there.

**Ensembl support** (accessions `ENST` / `ENSG` / `ENSP`) — downloads Ensembl cdot metadata and cDNA
FASTAs (~1 GB+); off by default. Without it, an Ensembl input reports "Reference not found":

```bash
ferro prepare --output-dir ferro-reference --ensembl
```

**RefSeqGene placements** — derives version-independent `NG_` placements and the
`NG_`→transcript-version map, required to resolve legacy gene-symbol selectors (`NG_(GENE):c.…`) and
bare-`NG_` hosted lookups:

```bash
ferro prepare --output-dir ferro-reference \
  --derive-ng-placements path/to/ng_accessions.txt
```

A fully-provisioned reference combines both in one run:

```bash
ferro prepare --output-dir ferro-reference --ensembl \
  --derive-ng-placements path/to/ng_accessions.txt
```

## Genomic-parent placements

`ferro prepare` also downloads genomic-parent placement data used to project transcript-coordinate variants into a genomic parent's own frame (#480):
- **RefSeqGene→genome alignments** (`GCF_*_refseqgene_alignments.gff3`, from NCBI's archived `alignments/ARCHIVE/all/` — the feed stopped updating in 2024; latest GRCh38 is the RS-109/p13 snapshot, valid for p14 since primary `NC_` accessions are unchanged) → manifest `refseqgene_alignments`; parsed by `MultiFastaProvider` into per-`NG_` `GenomicPlacement`. The corresponding GRCh37 snapshot (`GCF_000001405.25_105.*`) → manifest `refseqgene_alignments_grch37`, merged with the GRCh38 file so an `NG_` resolves to its build-appropriate placement (#653/#713).
- **LRG XML** genomic mapping → per-`LRG_` `GenomicPlacement` (parsed on demand).

`ReferenceProvider::genomic_placement` exposes these; `VariantProjector::project_to_genomic` composes the `NM_`→`NC_` (cdot) step with the `NC_`→`NG_`/`LRG_` transform to re-express coordinates in the parent frame. With no placement it declines rather than emit chromosome coordinates under the parent accession. That transform is **affine only when the placement carries no alignment gaps** — an LRG record lists `<diff>` elements alongside its `<mapping_span>`, and since #1833 the indel ones are carried on the placement and honoured, so 46 of the prepared reference's 1294 LRG records map through a gap list rather than a pure offset (#1499).

## Using it

Point any normalizing command at the directory:

```bash
ferro normalize -i variants.txt --reference ferro-reference/
```

To guard against a reference that has drifted on disk, add `--strict-reference`, which hard-fails if
the reference's content no longer matches its recorded identity (the default warns and proceeds).

## Faster reads: the 2-bit sequence store

Normalization has to read reference bases constantly. Read straight from the FASTA text, every access
pays to strip newlines and upper-case the bytes — which profiling puts at roughly three-quarters of
`ferro normalize`'s runtime. `ferro prepare` avoids that by doing the work **once**: alongside the
usual files it writes a `sequence_store.pac` sidecar — the reference sequences decoded and packed two
bits per base (the same idea as bwa's `.pac`). At run time `normalize` and `project` memory-map that
file and read bases directly, so the per-access decode disappears. On a large batch this is roughly a
third faster, with **byte-for-byte identical output**.

The sidecar is written as part of the normal `prepare` step — there is no separate command:

```bash
# Writes ferro-reference/sequence_store.pac alongside the other files:
ferro prepare --output-dir ferro-reference

# Force a rebuild — for example after editing a reference FASTA in place:
ferro prepare --output-dir ferro-reference --force
```

There is nothing to turn on:

- **It is built automatically** by `ferro prepare`, and rebuilt when you re-prepare. `--force` rebuilds
  it even if one is already present; without `--force` a matching sidecar is kept.
- **It is loaded automatically** by `normalize` and `project` whenever a `sequence_store.pac` sits
  beside the reference. The one-time build cost is paid only by `prepare`, never by an ordinary run.
- **It fails safe.** If the sidecar is absent, was built for a different reference, or is unreadable,
  ferro silently falls back to reading the FASTA text. The result is the same either way — the store is
  a speed-up, never a source of different answers. If you edit a reference FASTA in place, re-run
  `ferro prepare` so the sidecar is rebuilt (a stale one is detected and ignored, not trusted).

To point at a store in a non-default location — or to force one to be built on first use — set
`FERRO_SEQUENCE_STORE` to its path:

```bash
FERRO_SEQUENCE_STORE=/path/to/sequence_store.pac \
  ferro normalize -i variants.txt --reference ferro-reference/
```
