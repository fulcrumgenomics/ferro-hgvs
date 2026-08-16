# Quickstart

This walks through parsing a variant, then normalizing one against reference data.

## Parse a variant

Parsing validates a description and returns its structure. It needs no reference data.

```bash
ferro parse "NM_000088.3:c.459A>G"
```

```python
import ferro_hgvs

variant = ferro_hgvs.parse("NM_000088.3:c.459A>G")
print(variant.variant_type)  # "coding"
print(variant.reference)     # "NM_000088.3"
print(str(variant))          # "NM_000088.3:c.459A>G"
```

## Normalize a variant

Normalization rewrites a description into its canonical form. It needs reference sequences, so first
prepare a reference (a one-time download — see [Reference data](../guide/reference-data.md)):

```bash
ferro prepare --output-dir ferro-reference
ferro check --reference ferro-reference
```

Then normalize:

```bash
ferro normalize "NM_000088.3:c.459del" --reference ferro-reference/
```

```python
import ferro_hgvs

normalizer = ferro_hgvs.Normalizer(reference_json="ferro-reference/cdot.json")
print(normalizer.normalize("NM_000088.3:c.459del"))
```

## Read the warnings

Normalization sometimes *repairs* a description in a way the normalized string does not itself
record — separately reported members merged into one `delins`
(`MEMBERS_COALESCED_FROM_REPORTED_FORM`), a reference-range insertion payload replaced by the bases
it denotes (`INSERTED_SEQUENCE_EXPANDED`), or a stated reference base that contradicted the reference
and was accepted anyway (`REFSEQ_MISMATCH`).

These are reported as `warning[CODE]: message` on **stderr** (and in the `warnings` array under
`--format json`, the `detail` column under `--format tsv`), so a pipeline reading only stdout will
not see them. Run `ferro explain <CODE>` for what any code means.

## Next steps

- [Normalize variants](../guide/normalize-variants.md) — batches, files, output formats, error modes.
- [Reference data](../guide/reference-data.md) — what `ferro prepare` downloads and the optional data.
- [CLI reference](../cli/overview.md) — every subcommand.
