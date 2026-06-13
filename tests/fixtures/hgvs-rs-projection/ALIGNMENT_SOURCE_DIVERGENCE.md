# Alignment-source divergence in the hgvs-rs projection corpus

This document is the re-runnable provenance for the two **data-source**
divergence clusters in the hgvs-rs projection corpus
(`tests/fixtures/hgvs-rs-projection/cases.json`):

- `alignment-source-skew`
- `transcript-selection-vs-uta`

It records *why* a large fraction of the corpus diverges from ferro's output
without being a ferro bug, and it backs the **structural** quarantine in
`tests/hgvs_rs_projection_tests.rs` (the `classify_divergence` classifier plus
the `SOURCE_SKEW_TRANSCRIPTS` gate list). The quarantine is auditable here, not
spread across ~870 brittle per-case annotations.

## Root cause: two different transcript→genome alignments

ferro projects variants on **NCBI-canonical RefSeq-GFF alignments**
(`cdot-0.2.32.refseq`). The hgvs-rs corpus expectations were generated against
the **2021 biocommons UTA snapshot** (`uta_20210129`, which stores `splign`
alignments). Where the two alignment sources disagree, ferro returns a result
that is *correct for its alignment* but differs from the corpus oracle.

RefSeq-GFF is the more principled **default** for a canonical / MANE-aligned
engine: it is NCBI's own authoritative transcript-to-genome alignment, whereas
UTA is a 2021-frozen derived snapshot. So these divergences are "ferro tracks
NCBI-canonical; the corpus tracks a 2021 UTA snapshot," not ferro errors. We
*document and classify* the divergence rather than regenerating the corpus
against ferro's source (which would destroy the corpus's value as an
independent, published oracle).

## The gate: ungapped → 0% divergence, gapped/boundary → 39–100%

A gate over the corpus split the per-transcript divergence rate by whether the
two alignments agree:

- **Ungapped transcripts** (RefSeq-GFF and UTA splign produce the same exon
  structure with no indels): **0% divergence**. ferro and the corpus agree
  exactly. These never appear in the quarantine.
- **Gapped / boundary-skewed transcripts** (the alignments differ by a
  within-exon `D`/`I` gap or an exon-boundary off-by-one): **39–100%
  divergence**, concentrated entirely on the affected transcript.

This is the evidence that the divergence is *alignment-source*, not algorithmic:
it is perfectly predicted by whether the two sources' CIGARs differ, and it is
zero everywhere they agree.

## Category A — alignment-source skew (cluster `alignment-source-skew`)

ferro returns the **expected base transcript** but at a **different coordinate**
(the `c.`/`n.` position prefix differs). Each transcript below has a UTA
`uta_20210129` splign CIGAR (transcript vs genome) that contains an indel or a
boundary shift relative to ferro's RefSeq-GFF alignment, which is exactly the
coordinate offset observed in the diverging cases.

| Base accession | UTA `uta_20210129` splign CIGAR | Skew character |
|---|---|---|
| `NM_000682`      | `891=9D2375=`                          | 9 bp deletion gap → +9 coordinate shift downstream of the gap |
| `NM_000804`      | `150=2D37=`                            | 2 bp deletion gap → +2 coordinate shift downstream |
| `NM_001077527`   | `189=1X348=2D87=1X83=...`              | mismatch run + 2 bp deletion gap |
| `NM_001080519`   | `177=1D`                               | terminal 1 bp deletion gap |
| `NM_003777`      | `131=1I7=`                             | 1 bp insertion gap → −1 coordinate shift downstream |
| `NM_006158`      | `238=1I82=`                            | 1 bp insertion gap → −1 coordinate shift downstream |
| `NM_001277115`   | (ungapped per exon)                    | exon-**boundary** off-by-1 vs UTA: intronic offsets differ by 1 |

`NM_001277115` is the boundary case: each exon is internally ungapped, but the
exon boundaries are placed one base differently from UTA, so intronic
`+N`/`-N` offsets differ by 1. It produces a coordinate skew without a
within-exon indel.

These seven accessions are the `SOURCE_SKEW_TRANSCRIPTS` gate list in
`tests/hgvs_rs_projection_tests.rs`. A `CoordinateSkew` on one of them is
auto-quarantined as `divergence_accepted` (cluster `alignment-source-skew`). A
`CoordinateSkew` on **any other** transcript is deliberately **kept as a FAIL**
so an incomplete gate list surfaces for review instead of being silently
accepted.

## Category B — transcript selection/absence (cluster `transcript-selection-vs-uta`)

ferro does **not** return the expected base transcript at all (a
selection-coverage miss): the expected accession is absent from, or not
prioritized by, ferro's cdot set relative to the 2021 UTA snapshot. The clearest
example is the noncoding (`NR_`) rows, where ferro returns an `NM_` transcript
at the same locus rather than the corpus's expected `NR_`. The classifier routes
these structurally — no returned consequence carries the expected base
accession → `SelectionMiss` → `divergence_accepted`.

## Safety invariant (what is NOT quarantined)

A divergence where ferro returns the **expected base transcript** AND the
**same coordinate position** but a **different edit form** (e.g. `dup` vs `ins`,
`del` vs `delACCTT`, or the protein-prediction parenthesization `p.Arg268Trp`
vs `p.(Arg268Trp)`) is **not** data-source skew — it is a genuine
convention/algorithm delta and stays a **FAIL** (`FormOnly`). This preserves
"red = real": the residual FAILs are the genuine signal, dispositioned
individually in a later PR.

## Regenerating the CIGARs (needs live UTA)

The UTA splign CIGARs above come from the `tx_exon_aln_v` view of the 2021
biocommons UTA Postgres dump (`uta_20210129`). To regenerate them you need a
live UTA instance (the dump is large and access-gated), then per transcript:

```sql
-- uta_20210129 schema; one row per exon, ordered 5'→3'
SELECT tx_ac, alt_ac, alt_aln_method, ord, cigar
FROM   uta_20210129.tx_exon_aln_v
WHERE  tx_ac = 'NM_000682.5'
  AND  alt_aln_method = 'splign'
ORDER  BY ord;
```

Concatenate the per-exon `cigar` values (5'→3') to get the transcript-level
CIGAR shown above. ferro's RefSeq-GFF alignment for the same transcript comes
from the cdot bundle in the reference manifest (`cdot-0.2.32.refseq`); the
divergence is the difference between the two CIGARs.

## See also

- `tests/hgvs_rs_projection_tests.rs` — `classify_divergence`,
  `SOURCE_SKEW_TRANSCRIPTS`, and the per-axis structural routing.
- `tests/fixtures/hgvs-rs-projection/cases.json` — the `clusters` registry
  entries `alignment-source-skew` and `transcript-selection-vs-uta`.
- `tests/fixtures/CORPUS_LAYOUT.md` — corpus layout and consumer table.
