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

## Using it

Point any normalizing command at the directory:

```bash
ferro normalize -i variants.txt --reference ferro-reference/
```

To guard against a reference that has drifted on disk, add `--strict-reference`, which hard-fails if
the reference's content no longer matches its recorded identity (the default warns and proceeds).
