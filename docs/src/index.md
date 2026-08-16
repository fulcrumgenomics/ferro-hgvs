# ferro-hgvs

A high-performance HGVS variant nomenclature parser and normalizer, written in Rust with Python
bindings. It supports every HGVS coordinate system (`g` / `c` / `n` / `r` / `p` / `m` / `o`) and
edit type (substitution, deletion, insertion, duplication, inversion, repeat).

Use it as a Python package, a Rust crate, or a command-line tool. API docs are on
[docs.rs](https://docs.rs/ferro-hgvs); source is on
[GitHub](https://github.com/fulcrumgenomics/ferro-hgvs).

## Start here

- **[Installation](getting-started/installation.md)** — `pip install ferro-hgvs`,
  `cargo install ferro-hgvs`, or the Rust crate.
- **[Quickstart](getting-started/quickstart.md)** — parse a variant, then normalize one.

## Guides

- **[Normalize variants](guide/normalize-variants.md)** — batches, files, formats, error modes.
- **[Reference data](guide/reference-data.md)** — what `ferro prepare` downloads, and the optional
  Ensembl / RefSeqGene data.
- **[Benchmarking](guide/benchmarking.md)** — timing ferro and comparing it against other HGVS tools.
- **[CLI reference](cli/overview.md)** — every subcommand.

## Shadow Spec

The **[Shadow Spec](shadow-spec/index.md)** is a per-line reading of the HGVS recommendations that
records exactly how ferro interprets each clause — with runnable example spellings, the form ferro
normalizes each input to, and the reasoning behind each decision. It is under construction, one spec
page at a time.
