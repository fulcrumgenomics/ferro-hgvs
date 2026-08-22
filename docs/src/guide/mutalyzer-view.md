# A Mutalyzer-style multi-axis view

`ferro mutalyzer` resolves one variant across every coordinate axis at once — the normalized form plus the genomic (`g`), coding (`c`), non-coding (`n`), RNA (`r`), and protein (`p`) descriptions — and prints them together, in the shape Mutalyzer's `/normalize` returns. Where `ferro project` asks you to pick one output axis, this shows all of them in a single call.

```bash
ferro mutalyzer --transcript NM_000059.4 --reference <dir> \
  'NC_000013.11:g.32316466T>G'
```

```text
Input:       NC_000013.11:g.32316466T>G
Normalized:  NC_000013.11:g.32316466T>G
Genomic:     NC_000013.11:g.32316466T>G
Coding:      NC_000013.11(NM_000059.4):c.6T>G
Noncoding:   NC_000013.11(NM_000059.4):n.205T>G
RNA:         NM_000059.4:r.(6u>g)
Protein:     NP_000050.3:p.(Pro2=)
Gene:        BRCA2
Warnings:    (none)
```

An axis that is not derivable from the input and reference is shown as `-` (for example, a bare-`NM_` coding input with no genome alignment has no `Genomic:` line). The `--transcript` accession is required for a bare genomic input and optional for a `c.`/`n.`/`r.` input that already names its transcript.

## Interface-compatible, not output-compatible

The *shape* is Mutalyzer-familiar; the *values* are ferro's own. Each axis is computed by ferro's spec-grounded [normalization](normalize-variants.md) and [projection](project-variants.md), so a description here can differ from what Mutalyzer would emit for the same input — that is expected, not a bug. Mutalyzer is not a spec oracle for ferro (its normalizer minimizes a weighted description length with constants dated 2014 and has no separation rule); see [Normalization rules](../reference/normalization-rules.md) for what ferro's output is contracted to be.

## JSON output

`--format json` (or `-f json`) emits one JSON object per input, with the axes as fields:

```bash
ferro mutalyzer --transcript NM_000059.4 --reference <dir> -f json \
  'NC_000013.11:g.32316466T>G'
```

```json
{"input":"NC_000013.11:g.32316466T>G","normalized":"NC_000013.11:g.32316466T>G","genomic":"NC_000013.11:g.32316466T>G","coding":"NC_000013.11(NM_000059.4):c.6T>G","noncoding":"NC_000013.11(NM_000059.4):n.205T>G","rna":"NM_000059.4:r.(6u>g)","protein":"NP_000050.3:p.(Pro2=)","gene_symbol":"BRCA2","transcript_id":"NM_000059.4","warnings":[]}
```

Pass `-i <file>` to read one variant per line and `-o <file>` to write the results.

## The same view from the APIs

The Rust and Python APIs expose the identical view through the variant projector.

```rust
use ferro_hgvs::project::VariantProjector;

let result = projector.mutalyzer("NC_000013.11:g.32316466T>G", "NM_000059.4")?;
println!("{result}"); // the aligned block above
let coding = result.coding.as_deref(); // Some("NC_000013.11(NM_000059.4):c.6T>G")
```

```python
import ferro_hgvs

projector = ferro_hgvs.VariantProjector.from_manifest("<dir>/manifest.json")
result = projector.mutalyzer("NC_000013.11:g.32316466T>G", transcript="NM_000059.4")
print(result.coding)   # NC_000013.11(NM_000059.4):c.6T>G
print(result.protein)  # NP_000050.3:p.(Pro2=)
print(result)          # the aligned block
```
