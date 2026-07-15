# Python guide: setting up a reference for normalization (without footguns)

This guide shows the **recommended** way to give `ferro_hgvs.Normalizer` a
reference in Python, and how to verify you got the capability you expect. The
short version: the easiest constructors can silently run in a reduced-capability
mode, so **always check `has_genomic_data()`** (or let the warning fail your
harness) before trusting a genome-aware result.

## Which path do I use?

| You have… | Use | Genome-aware rules run? |
|---|---|---|
| A downloaded RefSeq/Ensembl reference + genome | `ferro prepare` → `Normalizer.from_manifest("manifest.json")` | ✅ full |
| A **synthetic / local** construct (your own FASTA + annotation) and you need intronic / `g.` / exon-junction rules | `ferro_hgvs.convert_gff(...)` in-process (or `ferro convert-gff` / `ferro build-transcript`, all with `emit_genomic_sequences=True` / `--emit-genomic-sequences`) → `Normalizer(reference_json="transcripts.json")` | ✅ yes |
| Transcript-level normalization only (exonic SNV/indel shuffle; no introns, no `g.` axis) | a `transcripts.json` **without** `genomic_sequences` → `Normalizer(reference_json=...)` | ❌ by design (see the boundary below) |
| Just trying the API | `Normalizer()` (built-in toy data) | ❌ toy data |

"Genome-aware rules" = intronic positions (`c.10+5del`), the cross-exon/intron
3′-shift (#670), and the genomic (`g.`) axis / `project_to_genomic`. See
[the capability boundary](transcripts_json_schema.md#capability-boundary-transcriptsjson-vs-a-prepared-manifest).

## The one footgun to remember

`Normalizer(reference_json=...)` and `Normalizer()` are the most discoverable
constructors, and a transcript-only reference makes them run **reduced-capability**
— genome-dependent rules are skipped. This is not silent: a reduced-capability
build reports `has_genomic_data() == False` **and** emits a one-time
`UserWarning` at construction. There are two safe ways to handle it: either leave
the warning on as a safety net (if you don't otherwise check), or check capability
explicitly (recommended) — the two examples below show both.

```python
import ferro_hgvs

n = ferro_hgvs.Normalizer(reference_json="transcripts.json")
n.has_genomic_data()    # False → genome-dependent rules will NOT run
n.has_protein_data()    # protein sequences present?
n.reference_summary()   # {"provider_kind": "json"|"manifest"|"test_data", "has_genomic_data": ..., "has_protein_data": ...}
```

### Recommended: check capability explicitly and fail loud

Assert the capability you actually need, so a mis-provisioned reference fails at
load rather than producing quietly-degraded results across a whole batch. Because
you're checking `has_genomic_data()` yourself, the one-time reduced-capability
warning is redundant here and is silenced to keep logs clean — a `require_genome`
mismatch raises a clear error instead:

```python
import warnings
import ferro_hgvs


def load_normalizer(reference_json: str, *, require_genome: bool) -> ferro_hgvs.Normalizer:
    """Load a Normalizer, refusing a reference that lacks a capability you need."""
    with warnings.catch_warnings():
        # The reduced-capability UserWarning is a safety net for callers who do NOT
        # check; we check has_genomic_data() below, so silence the (redundant) warning.
        warnings.simplefilter("ignore", UserWarning)
        normalizer = ferro_hgvs.Normalizer(reference_json=reference_json)

    if require_genome and not normalizer.has_genomic_data():
        raise SystemExit(
            f"{reference_json} is not genome-capable: {normalizer.reference_summary()}. "
            "Re-emit it with `--emit-genomic-sequences`, or use Normalizer.from_manifest(...)."
        )
    return normalizer
```

This handles both cases correctly: `require_genome=True` on a transcript-only
reference raises, while `require_genome=False` (you *want* transcript-level only)
loads quietly.

The same `has_genomic_data()` / `reference_summary()` methods exist on
`VariantProjector`, `EquivalenceChecker`, `BatchProcessor`, and
`CoordinateMapper` — apply the same check. `VariantProjector.project_to_genomic`
in particular needs a genome.

## Example 1 — genome-aware normalization of a synthetic construct

Build a **genome-capable** `transcripts.json` from your own FASTA — the
`emit_genomic_sequences` option embeds the contig bytes so the genome-dependent
rules can run. You can do this entirely in-process (no `ferro` CLI required):

```python
import ferro_hgvs

# Many transcripts from an annotation + FASTA, without shelling out to the CLI:
ferro_hgvs.convert_gff(
    ferro_hgvs.ConvertGffConfig(
        gff="constructs.gff3",
        fasta="constructs.fa",
        output="constructs.json",
        emit_genomic_sequences=True,
    )
)
# Pass output=None instead to get the JSON back as report.transcripts_json.
```

or via the CLI:

```bash
# Single synthetic construct (one contig):
ferro build-transcript \
  --fasta construct.fa --cds-start 1 --cds-end 900 \
  --emit-genomic-sequences \
  -o construct.json

# Or many transcripts from an annotation + FASTA:
ferro convert-gff \
  --gff constructs.gff3 --fasta constructs.fa \
  --emit-genomic-sequences \
  -o constructs.json
```

```python
import ferro_hgvs

n = load_normalizer("construct.json", require_genome=True)   # asserts genome capability
assert n.reference_summary()["provider_kind"] == "json"

# Intronic / exon-junction / g.-axis normalization now runs. Use the transcript
# id from your reference (build-transcript defaults it to the FASTA contig name).
print(n.normalize("CONSTRUCT1:c.10+5del"))
```

If you forget `--emit-genomic-sequences`, `construct.json` carries the placement
but no genome; `load_normalizer(..., require_genome=True)` will refuse it (rather
than silently returning exon-confined results).

## Example 2 — a downloaded production reference

```bash
ferro prepare --output-dir ferro-reference     # downloads RefSeq + genome + cdot
ferro check   --reference ferro-reference       # sanity-check the prepared set
```

```python
import ferro_hgvs

n = ferro_hgvs.Normalizer.from_manifest("ferro-reference/manifest.json")
assert n.has_genomic_data()
print(n.normalize("NM_000088.3:c.589-1G>T"))    # intronic, genome-aware
```

## Example 3 — transcript-level only (no genome needed)

If your inputs are exonic SNVs/indels and you do **not** need introns or the `g.`
axis, a transcript-only reference is correct and lighter — just be explicit that
you expect reduced capability so it isn't mistaken for a genome-backed run:

```python
import ferro_hgvs

n = ferro_hgvs.Normalizer(reference_json="transcripts_only.json")
assert not n.has_genomic_data()                 # expected: transcript-level only
print(n.normalize("CONSTRUCT1:c.1152del"))      # exonic 3′-shuffle works
# n.normalize("CONSTRUCT1:c.10+5del")           # intronic → cannot normalize without a genome
```

## Verify from the shell

`ferro check` reports what a reference can do — point it at a prepared directory
or a standalone `transcripts.json`:

```bash
ferro check --reference construct.json
# === transcripts.json Check ===
#   Status: OK
#   Transcripts: 1
#   Genome-capable: yes
#   Protein data: no
```

## The capability boundary in one line

- **Any `transcripts.json`:** reference-allele checks, exonic 3′/5′ shuffle,
  dup/repeat canonicalization.
- **Needs `genomic_sequences` (genome-capable) or a prepared manifest:** intronic
  positions, cross-exon/intron 3′-shift, the `g.` axis / `project_to_genomic`.

Full details, the schema, and what the load-time validation does and does not
guarantee are in [`transcripts_json_schema.md`](transcripts_json_schema.md).
