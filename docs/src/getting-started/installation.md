# Installation

ferro ships as a Python package, a Rust crate, and a standalone command-line tool. Install whichever
fits your workflow — they share the same engine.

## Python

```bash
pip install ferro-hgvs
```

Pre-built wheels are available for Linux (x86_64, aarch64), macOS (x86_64 and Apple Silicon), and
Windows (x86_64) on Python 3.10+.

## Command-line tool

```bash
cargo install ferro-hgvs
```

This builds the `ferro` binary from source and places it on your `PATH`. It needs a Rust toolchain
([rustup](https://rustup.rs/)).

## Rust library

Add the crate to your project:

```toml
[dependencies]
ferro-hgvs = "0.14"
```

## Verify the install

```bash
ferro parse "NM_000088.3:c.459A>G"
```

```python
import ferro_hgvs
print(ferro_hgvs.parse("NM_000088.3:c.459A>G"))
```

Parsing needs no reference data. Normalization against real transcripts does — see
[Reference data](../guide/reference-data.md).
