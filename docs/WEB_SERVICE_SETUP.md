# Web Service Setup Guide

This guide explains how to set up the ferro-hgvs multi-tool web service for interactive HGVS variant normalization.

## Overview

The web service provides a REST API and web UI for normalizing HGVS variants using up to four tools:

| Tool | Implementation | Requirements |
|------|---------------|--------------|
| **ferro** | Native Rust | Reference data (~4-6GB) |
| **mutalyzer** | HTTP API or Python subprocess | Network access or Python package |
| **biocommons** | Python subprocess | UTA database + SeqRepo (~25GB) |
| **hgvs-rs** | Native Rust | UTA database + SeqRepo |

You can enable any combination of tools based on your needs and available resources.

## Quick Start (Ferro Only)

The simplest setup uses only the ferro tool:

```bash
# 1. Build with web-service feature
cargo build --release --features web-service

# 2. Prepare ferro reference data
cargo run --release -- prepare --output-dir data/ferro-reference

# 3. Generate config file
cargo run --release --features web-service --bin ferro-web -- config -o config/service.toml

# 4. Edit config to set your reference path
# Update tools.ferro.reference_dir in config/service.toml

# 5. Start the server
cargo run --release --features web-service --bin ferro-web -- serve -c config/service.toml --open
```

The `--open` flag automatically opens your browser to the web UI.

## Configuration

### Generate Sample Config

```bash
ferro-web config -o config/service.toml
```

### Configuration File Structure

```toml
[server]
host = "0.0.0.0"
port = 3000
request_timeout_seconds = 60
enable_cors = true

[tools.ferro]
enabled = true
reference_dir = "/path/to/ferro-reference"
shuffle_direction = "3prime"
error_mode = "lenient"

[tools.mutalyzer]
enabled = true
mode = "api"  # or "local"
api_url = "https://mutalyzer.nl"
timeout_seconds = 30
rate_limit_ms = 100

[tools.biocommons]
enabled = false
uta_url = "postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b"
seqrepo_path = "/path/to/seqrepo/2021-01-29"

[tools.hgvs_rs]
enabled = false
uta_url = "postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b"
uta_schema = "uta_20210129b"
seqrepo_path = "/path/to/seqrepo/2021-01-29"
```

## Tool Setup

### Ferro (Native Rust)

Ferro requires reference data for normalization:

```bash
# Download and prepare reference data (~4-6GB)
ferro prepare --output-dir /path/to/ferro-reference

# Verify setup
ferro check --reference /path/to/ferro-reference
```

Config:
```toml
[tools.ferro]
enabled = true
reference_dir = "/path/to/ferro-reference"
shuffle_direction = "3prime"  # or "5prime"
error_mode = "lenient"        # or "strict", "silent"
```

### Mutalyzer

Mutalyzer can run in two modes:

#### API Mode (Default)

Calls the public mutalyzer.nl API. Simple but requires network access and is rate-limited.

```toml
[tools.mutalyzer]
enabled = true
mode = "api"
api_url = "https://mutalyzer.nl"
timeout_seconds = 30
rate_limit_ms = 100  # Be respectful to public API

[tools.mutalyzer.connection_pool]
max_connections = 10

[tools.mutalyzer.circuit_breaker]
failure_threshold = 5
recovery_timeout_seconds = 60
```

#### Local Mode

Runs mutalyzer via Python subprocess, same as ferro-benchmark. Faster and works offline with cached data.

**Requirements:**
```bash
pip install mutalyzer
```

**Config:**
```toml
[tools.mutalyzer]
enabled = true
mode = "local"
api_url = "https://mutalyzer.nl"  # Ignored in local mode
timeout_seconds = 60
settings_file = "/path/to/mutalyzer_settings.txt"  # Optional
allow_network = true  # Set false for cache-only operation
```

**Settings file format (optional):**
```
MUTALYZER_CACHE_DIR=/path/to/mutalyzer/cache
MUTALYZER_FILE_CACHE_ADD=true
```

### Biocommons/HGVS

Biocommons requires UTA database and SeqRepo. See [BIOCOMMONS_LOCAL_SETUP.md](BIOCOMMONS_LOCAL_SETUP.md) for detailed setup.

**Quick setup:**
```bash
# 1. Start UTA database (Docker)
docker run -d --name ferro-uta -p 5432:5432 biocommons/uta:uta_20210129b

# 2. Download SeqRepo (~10GB)
pip install biocommons.seqrepo
seqrepo --root-directory /path/to/seqrepo pull -i 2021-01-29

# 3. Install hgvs package
pip install hgvs
```

**Config:**
```toml
[tools.biocommons]
enabled = true
uta_url = "postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b"
seqrepo_path = "/path/to/seqrepo/2021-01-29"
docker_container = "ferro-uta"  # Optional: auto-start container
```

### HGVS-RS

HGVS-RS is a native Rust implementation requiring the same databases as biocommons.

**Build with feature:**
```bash
cargo build --release --features "web-service hgvs-rs"
```

**Config:**
```toml
[tools.hgvs_rs]
enabled = true
uta_url = "postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b"
uta_schema = "uta_20210129b"
seqrepo_path = "/path/to/seqrepo/2021-01-29"
```

## Running the Server

### Basic Usage

```bash
# Start with config file
ferro-web serve -c config/service.toml

# Override host/port
ferro-web serve -c config/service.toml --host 127.0.0.1 --port 8080

# Auto-open browser
ferro-web serve -c config/service.toml --open

# Adjust log level
ferro-web serve -c config/service.toml --log-level debug
```

### Health Check

```bash
curl http://localhost:3000/health
```

Response:
```json
{
  "status": "healthy",
  "available_tools": ["ferro", "mutalyzer"],
  "unavailable_tools": [],
  "tools": [
    {"tool": "ferro", "available": true, "status": "Available"},
    {"tool": "mutalyzer", "available": true, "status": "Available"}
  ]
}
```

## API Endpoints

### Normalize Single Variant

```bash
curl -X POST http://localhost:3000/api/v1/normalize \
  -H "Content-Type: application/json" \
  -d '{"hgvs": "NM_000088.3:c.589G>T"}'
```

### Normalize with Specific Tools

```bash
curl -X POST http://localhost:3000/api/v1/normalize \
  -H "Content-Type: application/json" \
  -d '{
    "hgvs": "NM_000088.3:c.589G>T",
    "tools": ["ferro", "mutalyzer"]
  }'
```

### Batch Normalize

```bash
curl -X POST http://localhost:3000/api/v1/normalize/batch \
  -H "Content-Type: application/json" \
  -d '{
    "variants": [
      "NM_000088.3:c.589G>T",
      "NM_000088.3:c.590del"
    ],
    "tools": ["ferro"]
  }'
```

### Parse Variant

```bash
curl -X POST http://localhost:3000/api/v1/parse \
  -H "Content-Type: application/json" \
  -d '{"hgvs": "NM_000088.3:c.589G>T"}'
```

## Web UI

The web service includes an interactive UI at `http://localhost:3000/`.

Features:
- Enter HGVS variants for normalization
- Select which tools to use
- View results from all tools side-by-side
- See agreement/disagreement between tools
- View detailed error messages

## Resource Requirements

### Disk Space

| Component | Size | Notes |
|-----------|------|-------|
| **Ferro reference** | **11 GB** | Required for ferro tool |
| ├─ genome/ | 7.9 GB | GRCh37 + GRCh38 |
| ├─ transcripts/ | 1.3 GB | RefSeq transcripts |
| ├─ refseqgene/ | 823 MB | NG_ sequences |
| ├─ cdot/ | 546 MB | Transcript data |
| └─ lrg/ | 117 MB | LRG mappings |
| **UTA database** | **15 GB** | PostgreSQL for biocommons/hgvs-rs |
| **SeqRepo** | **10-13 GB** | Sequences for biocommons/hgvs-rs |
| **Mutalyzer cache** | **0-2 GB** | Optional, grows with usage |

### By Configuration

| Setup | Disk Space | Memory | Notes |
|-------|------------|--------|-------|
| Ferro only | ~11 GB | ~200 MB | Fastest, simplest |
| + Mutalyzer API | +0 | +0 | Requires network |
| + Mutalyzer local | +0-2 GB | +100 MB | Optional cache |
| + Biocommons/HGVS-RS | +25-28 GB | +2 GB | UTA + SeqRepo (shared) |
| **Full setup** | **~36-41 GB** | **~2.5 GB** | All four tools |

## Troubleshooting

### Ferro: "Reference directory does not exist"

Ensure you've run `ferro prepare` and the path in config matches:
```bash
ferro check --reference /path/to/ferro-reference
```

### Mutalyzer API: "HTTP 429 Too Many Requests"

Increase rate limiting in config:
```toml
rate_limit_ms = 200  # Increase delay between requests
```

Or switch to local mode.

### Mutalyzer Local: "mutalyzer not installed"

Install the Python package:
```bash
pip install mutalyzer
```

### Biocommons: "Could not connect to UTA"

Check if the Docker container is running:
```bash
docker ps -f name=ferro-uta
docker start ferro-uta  # If stopped
```

### Biocommons: "SeqRepo path does not exist"

Verify the path and instance name:
```bash
ls /path/to/seqrepo/2021-01-29/aliases.sqlite3
```

### HGVS-RS: "feature not enabled"

Rebuild with the feature:
```bash
cargo build --release --features "web-service hgvs-rs"
```

## Example Configurations

### Minimal (Ferro Only)

```toml
[server]
host = "0.0.0.0"
port = 3000

[tools.ferro]
enabled = true
reference_dir = "/data/ferro-reference"
```

### Ferro + Mutalyzer API

```toml
[server]
host = "0.0.0.0"
port = 3000

[tools.ferro]
enabled = true
reference_dir = "/data/ferro-reference"

[tools.mutalyzer]
enabled = true
mode = "api"
api_url = "https://mutalyzer.nl"
rate_limit_ms = 100
```

### Full Local Setup

```toml
[server]
host = "0.0.0.0"
port = 3000

[tools.ferro]
enabled = true
reference_dir = "/data/ferro-reference"

[tools.mutalyzer]
enabled = true
mode = "local"
allow_network = true

[tools.biocommons]
enabled = true
uta_url = "postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b"
seqrepo_path = "/data/seqrepo/2021-01-29"

[tools.hgvs_rs]
enabled = true
uta_url = "postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b"
uta_schema = "uta_20210129b"
seqrepo_path = "/data/seqrepo/2021-01-29"
```

## See Also

- [BIOCOMMONS_LOCAL_SETUP.md](BIOCOMMONS_LOCAL_SETUP.md) - Detailed UTA/SeqRepo setup
- [BENCHMARK_GUIDE.md](BENCHMARK_GUIDE.md) - Tool comparison benchmarking
- [mutalyzer documentation](https://mutalyzer.nl/documentation)
- [biocommons/hgvs documentation](https://hgvs.readthedocs.io/)
