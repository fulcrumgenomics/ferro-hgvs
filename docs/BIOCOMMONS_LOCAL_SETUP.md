# Local biocommons/hgvs Setup Guide

This guide explains how to set up local UTA (Universal Transcript Archive) and SeqRepo databases for high-performance biocommons/hgvs normalization.

## Overview

The biocommons/hgvs Python library requires access to:
1. **UTA (Universal Transcript Archive)**: PostgreSQL database with transcript data
2. **SeqRepo**: Local sequence repository for nucleotide/protein sequences

By default, biocommons/hgvs uses the remote UTA server at `uta.biocommons.org`, which is slow (~5 seconds per pattern). For bulk benchmarking, setting up local databases provides ~100x speedup.

## Resource Requirements

| Component | Size | Notes |
|-----------|------|-------|
| UTA Docker image | ~3GB compressed | ~15GB uncompressed PostgreSQL |
| SeqRepo data | ~10-13GB | Sequence repository |
| PostgreSQL RAM | 2GB+ | Recommended shared buffers |
| **Total disk** | **~25-30GB** | |
| **Setup time** | **1-2 hours** | Depends on network speed |

## Prerequisites

- Docker installed and running
- Python 3.10+ with `hgvs` and `seqrepo` packages
- `pixi` environment (recommended) or pip install

### Using pixi (recommended)

```bash
pixi shell
```

This activates the environment with all required Python dependencies.

### Manual pip install

```bash
pip install hgvs biocommons.seqrepo
```

## Quick Start (Automated Setup)

### 1. Download UTA Dump

First, download the UTA dump file (requires browser due to human verification):
- Download from: https://dl.biocommons.org/uta/uta_20210129b.pgd.gz

### 2. Prepare biocommons (with automatic UTA setup)

```bash
# Full setup with auto UTA setup (recommended)
# Patterns are auto-detected from ferro reference if prepared with --patterns
ferro-benchmark prepare biocommons --output-dir data/ferro \
  --seqrepo-dir data/seqrepo \
  --uta-dump path/to/uta_20210129b.pgd.gz \
  --ferro-reference data/ferro

# If UTA is already running (e.g., from a previous setup), omit --uta-dump
ferro-benchmark prepare biocommons --output-dir data/ferro \
  --seqrepo-dir data/seqrepo --ferro-reference data/ferro
```

This will:
- Set up UTA database in Docker (if `--uta-dump` provided and UTA not running)
- Download SeqRepo data (~10GB)
- Initialize the local sequence repository
- Load additional sequences from ferro reference (if `--ferro-reference` provided):
  - **NC_** (chromosomes): GRCh38 and GRCh37 genomes
  - **NG_** (RefSeqGene): Gene-specific reference sequences
  - **LRG_** (Locus Reference Genomic): Stable genomic references (genomic patterns only)
- Auto-detect patterns from `{ferro-reference}/patterns/` and fetch missing accessions from NCBI
- Generate settings file at `data/ferro/biocommons_settings.conf`

**Note**: You must specify `--seqrepo-dir`. Choose a location with sufficient disk space.

> **Tip**: The `--ferro-reference` option is recommended when you need to normalize genomic patterns (NC_, NG_, LRG_). Without these sequences, biocommons will attempt slow network fetches to NCBI/EBI.

> **LRG Support**:
> - **Genomic patterns** (`LRG_1:g.5000A>T`): Fully supported via SeqRepo aliases
> - **Transcript patterns** (`LRG_199t1:c.79del`): Supported via LRG-to-RefSeq translation. Use `--lrg-mapping`:
>   ```bash
>   ferro-benchmark normalize biocommons -i patterns.txt -o results.json \
>     --biocommons-settings data/ferro/biocommons_settings.conf \
>     --lrg-mapping data/ferro/lrg/lrg_refseq_mapping.txt
>   ```
>   This translates `LRG_199t1` → `NM_004006.2` before normalization. The translation is noted in the output JSON.

### 3. Verify Setup

```bash
ferro-benchmark check biocommons --seqrepo-path /path/to/seqrepo/2021-01-29
```

### 4. Normalize and Compare

```bash
# Normalize with biocommons
ferro-benchmark normalize biocommons -i patterns.txt -o biocommons_norm.json \
  --biocommons-settings data/ferro/biocommons_settings.conf

# Compare with ferro results
ferro-benchmark compare results normalize ferro_norm.json biocommons_norm.json -o comparison.json
```

## Manual Setup

### UTA Database (Docker)

1. Pull the UTA image:
```bash
docker pull biocommons/uta:uta_20210129b
```

2. Start the container:
```bash
docker run -d \
  --name ferro-uta \
  -p 5432:5432 \
  biocommons/uta:uta_20210129b
```

3. Verify connection:
```bash
psql -h localhost -U anonymous -d uta -c "SELECT 1"
# Password: anonymous
```

### SeqRepo

1. Install seqrepo:
```bash
pip install biocommons.seqrepo
```

2. Initialize and download:
```bash
seqrepo --root-directory /path/to/seqrepo pull -i 2021-01-29
```

3. Verify:
```bash
ls /path/to/seqrepo/2021-01-29/
# Should show: aliases.sqlite3  sequences/
```

## Container Management

### Start UTA (if stopped)
```bash
ferro-benchmark setup start-uta --container-name ferro-uta
# or
docker start ferro-uta
```

### Stop UTA
```bash
ferro-benchmark setup stop-uta --container-name ferro-uta
# or
docker stop ferro-uta
```

### Check status
```bash
docker ps -f name=ferro-uta
```

## Upgrading UTA for Newer Transcripts

The default UTA (uta_20210129b) may lack transcripts added after January 2021.

### Check Available Versions

Visit https://dl.biocommons.org/uta/ for available UTA dumps.

### Upgrade Steps

```bash
# 1. Stop and remove existing UTA
ferro-benchmark setup stop-uta
docker rm ferro-uta

# 2. Download newer UTA dump (check the website above for latest)
curl -o uta_20231004b.pgd.gz https://dl.biocommons.org/uta/uta_20231004b.pgd.gz

# 3. Setup with new dump
ferro-benchmark setup uta --uta-dump uta_20231004b.pgd.gz --image-tag uta_20231004b

# 4. Update your settings file to use new schema
# The UTA_DB_URL should point to the new schema:
# postgresql://anonymous:anonymous@localhost:5432/uta/uta_20231004b
```

> **Note**: Newer UTA versions include more recent transcript mappings from RefSeq. If you're getting "transcript not found" errors for newer transcripts (NM_ accessions from 2022+), upgrading UTA may help.

## Settings File Format

The settings file uses a simple key=value format:

```
UTA_DB_URL = postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b
HGVS_SEQREPO_DIR = /path/to/seqrepo/2021-01-29
```

## Environment Variables

You can also configure biocommons/hgvs using environment variables:

```bash
export UTA_DB_URL="postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b"
export HGVS_SEQREPO_DIR="/path/to/seqrepo/2021-01-29"
```

## Troubleshooting

### Docker not found
```
Error: Docker is not available
```
Install Docker and ensure the daemon is running:
```bash
docker --version
docker ps
```

### UTA connection failed
```
Error: Could not connect to UTA at localhost:5432
```
Check if the container is running:
```bash
docker ps -f name=ferro-uta
# If not running:
docker start ferro-uta
```

### SeqRepo not found
```
Error: SeqRepo directory does not exist or is invalid
```
Verify the directory and check for required files:
```bash
ls /path/to/seqrepo/2021-01-29/aliases.sqlite3
ls /path/to/seqrepo/2021-01-29/sequences/
```

### psycopg2 build errors
If you encounter build errors with psycopg2, ensure PostgreSQL development headers are installed:
- macOS: `brew install postgresql`
- Linux: `apt install libpq-dev` or `yum install postgresql-devel`

Or use pixi which includes PostgreSQL.

### SeqRepo download fails
SeqRepo downloads from NCBI. If downloads fail:
1. Check your internet connection
2. Try again later (NCBI may be overloaded)
3. Consider using a mirror if available

## Performance Comparison

| Setup | Speed | Use Case |
|-------|-------|----------|
| Remote UTA | ~5 sec/pattern | Quick testing (<100 patterns) |
| Local UTA only | ~0.5 sec/pattern | Medium testing |
| Local UTA + SeqRepo | ~0.05 sec/pattern | Bulk benchmarking (10K+ patterns) |

## See Also

- [biocommons/hgvs documentation](https://hgvs.readthedocs.io/)
- [UTA documentation](https://github.com/biocommons/uta)
- [SeqRepo documentation](https://github.com/biocommons/biocommons.seqrepo)
