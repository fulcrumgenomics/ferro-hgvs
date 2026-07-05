# Benchmark Runbook: Perf Matrix + Reference Stack Setup

This runbook guides a maintainer through the full setup and execution of the `benchmark matrix` command that produces `data/benchmark/perf_results.json`, the source of truth for the README performance tables.

Run this a few times a year when refreshing the published numbers, or when tooling versions change significantly.

**Audience:** A maintainer who already has the reference stack on a local scratch volume. The run is manual — it is never executed in CI.

---

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [One-Time / Per-Session Setup](#one-time--per-session-setup)
3. [Running the Matrix](#running-the-matrix)
4. [Refreshing the Published Tables](#refreshing-the-published-tables)
5. [Methodology Notes](#methodology-notes)
6. [Appendix: UTA Dump Recovery](#appendix-uta-dump-recovery)

---

## Prerequisites

### Software

| Requirement | Notes |
|-------------|-------|
| Docker | Running daemon required; used for the UTA PostgreSQL container |
| PostgreSQL client | `psql` + `pg_isready`; `setup uta` probes the UTA database over the published host port. Install e.g. `postgresql-client` (Debian/Ubuntu), `libpq` (Homebrew), or `postgresql` (conda/pixi) |
| pixi | Manages the Python tool environment (mutalyzer, biocommons/hgvs, hgvs-rs) |
| Rust toolchain | `cargo build --release` for the ferro-benchmark binary |
| System build libraries | The `benchmark` / `hgvs-rs` features link native C libraries, so a fresh machine needs a C toolchain plus dev headers. In particular the `benchmark` feature's `rusqlite` links **system SQLite** — without it the build fails at link time with `rust-lld: error: unable to find library -lsqlite3`. Install: **RHEL / Amazon Linux 2023** — `gcc gcc-c++ make cmake sqlite-devel zlib-devel bzip2-devel xz-devel openssl-devel`; **Debian / Ubuntu** — `build-essential cmake libsqlite3-dev zlib1g-dev libbz2-dev liblzma-dev libssl-dev`; **macOS** — Xcode Command Line Tools plus Homebrew `sqlite`. |

### Reference Stack

The reference stack lives on a machine-specific volume. All commands in this runbook refer to its root through the placeholder path:

```
/path/to/ferro-bench-data
```

Set it once and use throughout:

```bash
export D=/path/to/ferro-bench-data
```

Replace `/path/to/ferro-bench-data` with wherever the stack actually lives on your machine. The variable is not persisted by any config file — export it in every shell session before running benchmark commands. The `MUTALYZER_CACHE_DIR` / `HGVS_SEQREPO_DIR` settings written below use the same placeholder and must be edited to match.

The stack contains:

| Path | Contents |
|------|----------|
| `$D/ferro/` | ferro manifest, cdot, genomes, seqrepo-derived reference data |
| `$D/seqrepo/2021-01-29/` | SeqRepo instance (biocommons / hgvs-rs) |
| `$D/mutalyzer/` | Mutalyzer cache (~1.1 M files) |
| `$D/clinvar/clinvar_patterns.txt` | Full ClinVar population (~57.5 M lines) |

### Disk

The full stack requires ~50 GB of disk. Ensure adequate free space on the volume before running.

---

## One-Time / Per-Session Setup

Do all steps below at the start of any session where you will run the matrix. Steps 1 and 2 are one-time; steps 3–5 must be repeated after every machine reboot.

### Step 1: Build the binary

```bash
cargo build --release --features benchmark,hgvs-rs --bin ferro-benchmark
```

Verify:

```bash
./target/release/ferro-benchmark --help
```

### Step 2: Verify the Python environment

Confirm that all required Python packages are importable through the pixi environment:

```bash
pixi run python -c "import mutalyzer_hgvs_parser, hgvs, biocommons.seqrepo; print('ok')"
```

Expected output: `ok`. The pixi environment pins mutalyzer==3.1.1, mutalyzer-hgvs-parser, hgvs>=1.5.6, and biocommons-seqrepo. If the import fails, run `pixi install` to restore the environment.

### Step 3: Start the UTA Docker container

UTA is the PostgreSQL database used by biocommons/hgvs and hgvs-rs for normalization. Start it with:

```bash
pixi run ./target/release/ferro-benchmark setup uta \
  --container-name ferro-uta \
  --image-tag uta_20210129b
```

This pulls `biocommons/uta:uta_20210129b` if not already present, starts the container, and waits (polling the host port for the `uta_20210129b` schema) until the database is ready.

> **⚠️ The image no longer self-loads its data — you must supply a local dump.** The `biocommons/uta` image ships **empty**; on first boot its `load-uta.sh` init script downloads the ~2–3 GB dump from `https://dl.biocommons.org/uta/uta_20210129b.pgd.gz`. That URL is now behind a **human-verification wall** (`302 → /human-verify.html`), so the container's `curl` saves a 154-byte HTML page instead of the dump, `gzip` rejects it (`not in gzip format`), and the container comes up with an **empty schema** — `setup uta` then times out waiting for `uta_20210129b` (only the `public` schema exists). A plain `setup uta --image-tag …` therefore **cannot** provision UTA on a fresh host. You must obtain `uta_20210129b.pgd.gz` once by clicking through the verification page in a browser, then load it locally with **`setup uta --uta-dump /path/to/uta_20210129b.pgd.gz`** (see [Appendix: UTA Dump Recovery](#appendix-uta-dump-recovery)). Keep the dump on a stable local/scratch path; it is gitignored and not distributed in the repo.

**Readiness wait.** The wait is patient and configurable via `--uta-ready-timeout-secs` (default **300**); the biocommons/uta image can take a few minutes to initialize its schema on first start, so a heartbeat reports progress while it waits. If it does time out, the error is state-specific (container exited → check `docker logs`; postgres never came up; or schema never appeared) — re-run `setup uta` (it is **idempotent**) with a larger `--uta-ready-timeout-secs` if the disk is slow.

**Idempotency & re-runs.** Re-running `setup uta` is safe: a ready container is a no-op; a stopped one is started; one still initializing is waited on. If a container of the same name already exists with a **different published port or image** than requested, the command errors rather than silently using the wrong one — re-run with `--force` to recreate, or pass a matching `--uta-port` / `--uta-image-tag`.

**Port collisions.** The default `--uta-port` is 5432, commonly occupied by a local Postgres. If `docker run` reports the port is already allocated, pass `--uta-port <other>`.

**Offline / verification-gated dump.** `setup uta` pulls a Docker image; it never auto-downloads the dump, and (per the warning above) the image's *own* boot-time download is now broken by the biocommons human-verification wall. So on any fresh host you must pass a local dump with `--uta-dump <path>` (see the recovery appendix). The file must be a real gzip `.pgd.gz` — a saved HTML verification page is rejected with a clear error before any Docker work (this guard exists precisely because the wall serves HTML in place of the dump).

### Step 4: Discover the UTA host port and verify reachability

The container maps its internal port 5432 to a host port that can vary between sessions. Find the actual mapping:

```bash
docker ps --format '{{.Names}} {{.Ports}}' | grep ferro-uta
```

Example output:
```
ferro-uta 0.0.0.0:55432->5432/tcp
```

The host port in this example is `55432`. Substitute the actual value for `<PORT>` in all subsequent commands. (`5432` is also common — always check.)

Verify the database is reachable:

```bash
pixi run python -c "
import psycopg2
psycopg2.connect('postgresql://anonymous:anonymous@localhost:<PORT>/uta').close()
print('ok')
"
```

If this prints `ok`, UTA is up. If it raises a connection error, wait 10–15 seconds and retry — the database may still be initializing after a fresh container start.

**Schema-suffixed URL:** biocommons/hgvs and hgvs-rs both require the URL to include the schema name as a path suffix. Always use:

```
postgresql://anonymous:anonymous@localhost:<PORT>/uta/uta_20210129b
```

The unsuffixed URL (`…/uta`) connects to the database but does not set the schema search path and will cause lookup failures.

### Step 5: Write corrected settings files

The shipped settings files contain stale or relative paths. Write corrected copies to `/tmp/bench-settings/` — this directory is transient and safe to regenerate each session.

```bash
mkdir -p /tmp/bench-settings
```

**mutalyzer settings** (`/tmp/bench-settings/mutalyzer_settings.conf`):

```bash
cat > /tmp/bench-settings/mutalyzer_settings.conf << 'EOF'
MUTALYZER_CACHE_DIR=/path/to/ferro-bench-data/mutalyzer
MUTALYZER_FILE_CACHE_ADD=False
EOF
```

`MUTALYZER_FILE_CACHE_ADD=False` prevents the benchmark run from attempting to write new entries to the cache.

**biocommons settings** (`/tmp/bench-settings/biocommons_settings.conf`):

```bash
# Replace <PORT> with the actual host port discovered above (e.g. 55432)
cat > /tmp/bench-settings/biocommons_settings.conf << 'EOF'
UTA_DB_URL = postgresql://anonymous:anonymous@localhost:<PORT>/uta/uta_20210129b
HGVS_SEQREPO_DIR = /path/to/ferro-bench-data/seqrepo/2021-01-29
EOF
```

The shipped biocommons settings use a relative SeqRepo path that only resolves correctly when the working directory is `$D`. The corrected settings use the absolute path and the schema-suffixed UTA URL with the actual host port.

### Step 6: Check all four tools

Confirm each tool is ready before running the matrix:

```bash
D=/path/to/ferro-bench-data

# ferro
pixi run ./target/release/ferro-benchmark check ferro \
  --reference "$D/ferro"

# mutalyzer
pixi run ./target/release/ferro-benchmark check mutalyzer \
  --mutalyzer-settings /tmp/bench-settings/mutalyzer_settings.conf

# biocommons
pixi run ./target/release/ferro-benchmark check biocommons \
  --seqrepo-path "$D/seqrepo/2021-01-29" \
  --uta-db-url "postgresql://anonymous:anonymous@localhost:<PORT>/uta/uta_20210129b"

# hgvs-rs
pixi run ./target/release/ferro-benchmark check hgvs-rs \
  --seqrepo-path "$D/seqrepo/2021-01-29" \
  --uta-db-url "postgresql://anonymous:anonymous@localhost:<PORT>/uta/uta_20210129b"
```

All four should print a ready / ✓ status. Do not proceed to the matrix run until all checks pass.

---

## Running the Matrix

The `benchmark matrix` subcommand runs all four tools across parse and normalize operations, multiple worker counts, and multiple repetitions over a shared seeded stratified sample, then writes the result to a JSON file.

**The run is not part of CI** — it requires the full reference stack (UTA container, SeqRepo, mutalyzer cache) and takes tens of minutes, bounded by mutalyzer normalize.

### Smoke run (validation first)

Run a smoke pass before committing to the confident run. The smoke pass uses a small sample and a single repetition, and completes in a few minutes.

```bash
D=/path/to/ferro-bench-data

pixi run ./target/release/ferro-benchmark benchmark matrix \
  --population "$D/clinvar/clinvar_patterns.txt" \
  --output /tmp/perf_smoke.json \
  --sample-size 20 \
  --reps 1 \
  --workers 1,8 \
  --ferro-threads 1,2,4,8 \
  --operations parse,normalize \
  --reference "$D/ferro" \
  --mutalyzer-settings /tmp/bench-settings/mutalyzer_settings.conf \
  --biocommons-settings /tmp/bench-settings/biocommons_settings.conf \
  --seqrepo-path "$D/seqrepo/2021-01-29" \
  --uta-db-url "postgresql://anonymous:anonymous@localhost:<PORT>/uta/uta_20210129b" \
  --machine "Apple M2 Max (12 cores)"
```

Verify the output is well-formed and all four tools produced results (or recorded a `not_run` reason rather than crashing):

```bash
python3 -c "
import json, sys
d = json.load(open('/tmp/perf_smoke.json'))
print('schema_version:', d['schema_version'])
print('placeholder:', d['placeholder'])
for op, v in d['operations'].items():
    for w, wp in v['by_workers'].items():
        for tool, m in wp['tools'].items():
            status = 'not_run' if m.get('not_run') else f\"{m['median_pps']:.1f} pps\"
            print(f'  {op} W={w} {tool}: {status}')
"
```

### Confident run (for publication)

Once the smoke pass is clean, run the confident pass. This is the run whose output becomes the published `data/benchmark/perf_results.json`.

**Do not pass `--sample-size` for the confident run.** `--sample-size` pins one fixed N for every tool and *skips calibration* — it exists only for smoke runs. With it set, fast tools (ferro, hgvs-rs) are measured over a sample far too small to reflect their real throughput (e.g. ferro parse at N=60 reports a few hundred thousand p/s instead of millions), and the per-tool `sample_sizes` and `ferro_scaling` blocks are omitted, so `generate_perf_tables` cannot render the thread-scaling table.

Instead, let the matrix **calibrate per-(tool, op) sample sizes**: it estimates each tool's rate, then picks N so each tool runs for about `--target-seconds` per rep, clamped to `[--min-sample, --max-sample]`. This is what produces the "fast tools see millions of patterns, slow tools hundreds" property described in the README footnote. The defaults (`--target-seconds 3 --min-sample 50 --max-sample 2000000`) are what the published numbers use; they give a total run time of a few tens of minutes at `--reps 5`, bounded by mutalyzer normalize.

```bash
D=/path/to/ferro-bench-data

pixi run ./target/release/ferro-benchmark benchmark matrix \
  --population "$D/clinvar/clinvar_patterns.txt" \
  --output data/benchmark/perf_results.json \
  --target-seconds 3 \
  --min-sample 50 \
  --max-sample 2000000 \
  --reps 5 \
  --ferro-full-n 1000000 \
  --workers 1,8 \
  --ferro-threads 1,2,4,8 \
  --operations parse,normalize \
  --reference "$D/ferro" \
  --mutalyzer-settings /tmp/bench-settings/mutalyzer_settings.conf \
  --biocommons-settings /tmp/bench-settings/biocommons_settings.conf \
  --seqrepo-path "$D/seqrepo/2021-01-29" \
  --uta-db-url "postgresql://anonymous:anonymous@localhost:<PORT>/uta/uta_20210129b" \
  --machine "Apple M2 Max (12 cores)"
```

### Key flags explained

| Flag | Purpose |
|------|---------|
| `--target-seconds` | Calibration target: per-(tool, op) N is chosen so each tool runs ~this long per rep. Default 3. |
| `--min-sample` / `--max-sample` | Floor and ceiling on calibrated N. Defaults 50 / 2,000,000. The floor keeps very slow tools statistically meaningful; the ceiling bounds memory for very fast tools. |
| `--sample-size` | **Smoke runs only.** Pins one fixed N for all tools and skips calibration. Do not use for the published run (see warning above). |
| `--reps` | Repetitions per (tool, worker-count) cell. Median and min/max are computed across reps. |
| `--ferro-full-n` | Cap on the ferro full-population headline pass. `0` = whole population; `1000000` = first 1 M lines. The full-population pass measures ferro only — it does not run the slow tools. |
| `--workers` | Worker counts for the cross-tool comparison tables (e.g. `1,8` produces W=1 and W=8 columns). |
| `--ferro-threads` | Thread counts for the ferro-only thread-scaling table (e.g. `1,2,4,8`). |
| `--exclude-protein` | Omits protein variants from normalize samples (default: true). Protein normalization is inconsistent across tools; excluding them ensures a fair comparison. |
| `--machine` | Human-readable machine label written to the `provenance` block in the JSON. |

---

## Refreshing the Published Tables

After a confident run writes `data/benchmark/perf_results.json`, regenerate the README tables:

```bash
# Render the README tables from the new results
cargo run --features dev --example generate_perf_tables

# Confirm the on-disk file and the rendered output are in sync
cargo run --features dev --example generate_perf_tables -- --check
```

The `--check` invocation exits non-zero if the tables in README.md do not match what would be generated from the current `perf_results.json`. It must exit 0 before committing.

Inspect the diff to sanity-check the numbers:

```bash
git diff README.md
```

Things to verify before committing:

- No placeholder banner (the `placeholder` field must be `false`).
- ferro parse and normalize are substantially faster than the other tools.
- W=8 throughput is ≥ W=1 for every tool.
- ferro thread-scaling is roughly monotonic from 1→8 threads.
- biocommons and hgvs-rs show real numbers for normalize (not `—`); if they are `—`, confirm UTA was up during the run.

Commit the results and the updated README:

```bash
git add data/benchmark/perf_results.json README.md
git commit -m "perf(benchmark): publish measured performance tables"
```

---

## Methodology Notes

### Population and sampling

All tool comparisons use a single stratified ClinVar population (`$D/clinvar/clinvar_patterns.txt`, ~57.5 M lines). Stratification ensures representation across variant types (genomic, coding, intronic, non-coding, RNA, protein). For each repetition `r`, the sample is drawn with seed `base_seed + r`. The same seed-derived sample is used for every tool in that repetition, so per-tool differences reflect tool speed rather than input variation.

### ferro full-population headline

In addition to the per-rep sampled runs, ferro is run once across the full population (or `--ferro-full-n` lines if set). This produces the headline patterns-per-second figure that appears at the top of the README table. The full-population pass is not repeated across reps — it runs once per operation.

### Aggregation

For each (tool, operation, worker-count) cell, `--reps` independent measurements are collected. The reported value is the **median** patterns/second across reps; the min and max are also recorded and available in `perf_results.json`. Median is preferred over mean because it is robust to warm-up outliers on the first rep.

### Single provider per ferro worker

For ferro parallel runs, one reference provider is constructed per rayon worker thread before the timer starts. Provider construction is excluded from the timed region, matching the discipline used by hgvs-rs parallel runs.

### Timed region (startup-excluded for all tools)

All four tools exclude process/interpreter startup from the timed region. The mutalyzer and biocommons subprocesses are timed by their own internal (startup-excluded) `time.perf_counter()` timer for both parse and normalize, matching the discipline used by ferro and hgvs-rs — so Python interpreter startup and package import time are not folded into any published throughput figure. This is disclosed in the `provenance.notes` field of `perf_results.json` and in the README footnote.

For mutalyzer normalize at more than one worker under the pre-sharded cache, the recorded internal time is the slowest shard's critical path — a lower bound that still excludes the dominant per-process startup term. This closed issue #609 (the prior methodology folded Python startup into mutalyzer/biocommons normalize throughput).

---

## Appendix: UTA Dump Recovery

If the `ferro-uta` Docker container is fresh or its data volume is empty (the schema is unloaded), the database will not have the `uta_20210129b` schema and tool checks will fail with errors like `relation "uta_20210129b.transcript" does not exist`. As explained in [Step 3](#step-3-start-the-uta-docker-container), this is now the **default** state on a fresh host, because the image's boot-time auto-download is broken by the biocommons human-verification wall.

### Obtaining the dump

The dump is **not** distributed in this repo (it is gitignored). You must download `uta_20210129b.pgd.gz` (~281 MB) yourself:

```text
https://dl.biocommons.org/uta/uta_20210129b.pgd.gz
```

`curl -O`/`wget` on that URL **does not work** — it is behind a JavaScript human-verification wall that returns a 154-byte HTML redirect stub instead of the dump (which is why the container's own `load-uta.sh` fails silently). **Download it in a web browser**, which clears the challenge, and save the file to a stable local path (e.g. `/path/to/uta_20210129b.pgd.gz`). Verify it is a real gzip before using it:

```bash
file /path/to/uta_20210129b.pgd.gz   # => "gzip compressed data"; NOT "HTML document"
gzip -t /path/to/uta_20210129b.pgd.gz && echo OK
```

Then either pass it to `setup uta --uta-dump /path/to/uta_20210129b.pgd.gz` (preferred — it loads and readies the container in one step), or load it into an already-running container by hand with the steps below.

To restore it into a running `ferro-uta` container by hand:

**Step 1: Copy the dump into the container**

```bash
# Use the dump you downloaded in "Obtaining the dump" above.
docker cp /path/to/uta_20210129b.pgd.gz ferro-uta:/tmp/dump.gz
```

**Step 2: Load the schema and data (first psql pass)**

```bash
docker exec -i ferro-uta bash -c \
  "gunzip -c /tmp/dump.gz | psql -U anonymous -d uta"
```

This pass creates the `uta_20210129b` schema and loads all tables. It may take several minutes.

**Step 3: Set the search path and refresh materialized views (second psql pass)**

```bash
docker exec -i ferro-uta psql -U anonymous -d uta << 'EOF'
SET search_path TO uta_20210129b, public;
REFRESH MATERIALIZED VIEW uta_20210129b.tx_def_summary_mv;
REFRESH MATERIALIZED VIEW uta_20210129b.tx_exon_aln_v;
EOF
```

**Step 4: Verify**

```bash
pixi run python -c "
import psycopg2
conn = psycopg2.connect('postgresql://anonymous:anonymous@localhost:<PORT>/uta/uta_20210129b')
cur = conn.cursor()
cur.execute('SELECT COUNT(*) FROM uta_20210129b.transcript')
print('transcript count:', cur.fetchone()[0])
conn.close()
"
```

A count in the hundreds of thousands confirms the schema loaded correctly. Proceed to [Step 5 of the per-session setup](#step-5-write-corrected-settings-files) once this check passes.
