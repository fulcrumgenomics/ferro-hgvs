# Normalize or rederive?

ferro exposes two entry points for canonicalizing an HGVS description, but they are different
operations with different guarantees. This guide says what each does, which axes each supports,
and when to reach for one over the other.

| | `normalize` | `rederive` |
|---|---|---|
| what it takes | an HGVS **description** | an HGVS **description** |
| what it does | shifts and re-spells the *given* description in place | throws the spelling away and re-derives from the **denoted bases** |
| axes | every axis — `g/c/n/r/p/m` | genomic only — `g.` (and `m.` on the two rCRS accessions) |
| confluence | best-effort (rule 3) — converges equivalent inputs | confluent **by construction**, over one window |
| needs a reference | yes | yes |

## `normalize` — canonicalize a description you were handed

Use `normalize` when you already have a description and want its canonical form: 3′/5′ shifting
and the recommendations' preferred representation, applied to the spelling you passed in. It works
on every axis, because a `c.`/`r.`/`p.` description carries information — reading frame, which
transcript, a protein consequence — that a genomic base window cannot recover.

```python
import ferro_hgvs

nz = ferro_hgvs.Normalizer(reference_json="…")
nz.normalize("NM_000088.3:c.100delA")   # -> "NM_000088.3:c.100del"
```

Converging equivalent inputs to one output is [confluence — rule 3 of the normalization
rules](../reference/normalization-rules.md), which `normalize` targets on a best-effort basis: a
handful of known residual defects remain, where two spellings of one variant still normalize to
two strings. Those are tracked bugs, not the intended contract.

## `rederive` — one canonical description per variant

Use `rederive` when the input's spelling should carry no weight — when you want *the* description
for a variant, regardless of how it happened to be written. `rederive` expresses the variant as
the bases it denotes (a reference/alternate window) and derives a description from those bases
alone, so two spellings of one variant reach one result:

```python
nz.rederive("NC_TEST.1:g.14_15insACGTACGT")   # -> "NC_TEST.1:g.7_14dup"
nz.rederive("NC_TEST.1:g.7_14dup")            # -> "NC_TEST.1:g.7_14dup"
```

Because it re-derives from a genomic sequence window, `rederive` is **genomic-only** and refuses a
transcript or protein axis. Project onto the axis you need afterwards if required.

`rederive` returns the alignment-derived form by default. Pass `recommended_form=True` to route
the result through `normalize` for ferro's recommended, reference-anchored form:

```python
nz.rederive("NC_TEST.1:g.11del", recommended_form=True)   # 3′-shifted per the recommendations
```

`rederive` is the single-call form of the `to_sequences` → `from_sequences` round trip described
in [Deriving a description from sequences](deriving-from-sequences.md); read that guide for the
mechanics and for the read-dependence cost of a window-local derivation.

## Which one?

- You have a description and want it in canonical form on its own axis → **`normalize`**.
- You want the same description for a variant no matter how two callers spelled it, on a genomic
  axis → **`rederive`**.
- You have raw bases (a reference/alternate window) rather than a description →
  [`from_sequences`](deriving-from-sequences.md).
