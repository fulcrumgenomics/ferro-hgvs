# Normalize variants

Normalization rewrites a description into its canonical form — 3′-shifting, collapsing equivalent
spellings, and applying the recommendations' preferred representation. This guide covers batches,
files, output formats, and error modes.

Normalization needs reference sequences. Prepare a reference first (see
[Reference data](reference-data.md)) and pass `--reference`.

## A single variant

```bash
ferro normalize "NM_000088.3:c.459del" --reference ferro-reference/
```

Read from stdin instead:

```bash
echo "NC_000001.11:g.12345A>G" | ferro normalize --reference ferro-reference/
```

## A batch from a file

One description per line:

```bash
ferro normalize -i variants.txt --reference ferro-reference/ -o normalized.txt
```

**Sort the input by transcript accession (or genomic position) for large batches.** ferro caches
each resolved transcript, so consecutive variants on the same transcript skip the dominant cost of
re-reading it. Sorted input keeps the working set resident and is markedly faster.

Use several workers for large batches:

```bash
ferro normalize -i variants.txt --reference ferro-reference/ -j 8
```

## Output formats

`-f/--format` selects the output:

- `text` (default) — the normalized description, one per line.
- `json` — a structured record per input, including a `warnings` array.
- `tsv` — a table with header `line, input, normalized, changed, status, detail`, plus a summary
  line on stderr. This is the format for answering *"which of my variants changed?"*

```bash
ferro normalize -i variants.txt --reference ferro-reference/ -f tsv > normalized.tsv
```

## Error modes

`--error-mode` controls how strict ferro is about its input and output:

- `strict` (default) — validates that the input conforms to the recommendations and fails on a
  violation.
- `lenient` — accepts a wider range of inputs and repairs where it can; fails only when it cannot
  normalize.
- `silent` — lenient, but without the diagnostic messages.

Regardless of mode, normalization may *repair* a description in a way the output string does not
record. Those repairs are reported as `warning[CODE]: message` on **stderr** (and in the `warnings`
array under `-f json`, the `detail` column under `-f tsv`). A pipeline reading only stdout will not
see them, so capture stderr or use `-f json`/`-f tsv`. The `message` on each `warning[CODE]:` line
says what the code means; use `--ignore` / `--reject` to tune specific codes.
