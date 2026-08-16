# ferro-hgvs

A high-performance HGVS variant nomenclature parser and normalizer, written in Rust with Python
bindings. It supports every HGVS coordinate system (`g` / `c` / `n` / `r` / `p` / `m` / `o`) and
edit type (substitution, deletion, insertion, duplication, inversion, repeat).

This site is the user-facing documentation. API docs are on [docs.rs](https://docs.rs/ferro-hgvs);
source is on [GitHub](https://github.com/fulcrumgenomics/ferro-hgvs).

## Shadow Spec

The **[Shadow Spec](shadow-spec/index.md)** is a per-line reading of the HGVS recommendations that
records exactly how ferro interprets each clause — with runnable example spellings, the form ferro
normalizes each input to, and links to the governing rulings. It is under construction, one spec
page at a time.
