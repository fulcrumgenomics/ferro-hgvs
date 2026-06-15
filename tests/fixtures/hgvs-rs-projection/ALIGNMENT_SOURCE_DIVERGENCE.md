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
(the `c.`/`n.` position prefix differs), because ferro's RefSeq-GFF alignment
places the transcript on the genome differently from the UTA splign alignment.
This cluster has two sub-kinds: **gapped-CIGAR skew** (a within-exon indel or
boundary off-by-one, this section) and **genomic-anchor skew** (ungapped CIGAR,
agreeing CDS start, uniform whole-transcript offset — the next section).

### Gapped-CIGAR skew

Each transcript below has a UTA `uta_20210129` splign CIGAR (transcript vs
genome) that contains an indel or a boundary shift relative to ferro's
RefSeq-GFF alignment, which is exactly the coordinate offset observed in the
diverging cases.

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

These seven accessions, plus the three genomic-anchor accessions in the next
section, make up the `SOURCE_SKEW_TRANSCRIPTS` gate list in
`tests/hgvs_rs_projection_tests.rs`. A `CoordinateSkew` on one of them is
auto-quarantined as `divergence_accepted` (cluster `alignment-source-skew`). A
`CoordinateSkew` on **any other** transcript is deliberately **kept as a FAIL**
so an incomplete gate list surfaces for review instead of being silently
accepted.

### Sub-kind: genomic-anchor skew (ungapped CIGAR, agreeing CDS start)

A manifest run surfaces three low-count `CoordinateSkew` FAILs on transcripts
whose UTA splign CIGARs are **ungapped** — each a uniform **+2** `c.`-coordinate
shift:

| transcript | expected → ferro | shift |
|---|---|---|
| `NM_020451` (SELENON) | `c.943G>A` → `c.945G>A` | +2 |
| `NM_007199` (IRAK3) | `c.1A>G` → `c.3A>G` | +2 |
| `NM_002386` (MC1R) | `c.-11_19del` → `c.-9_21del` | +2 |

When PR2 first surfaced these it flagged them as *candidate* alignment skew
pending a UTA CIGAR check, and PR3 noted their CIGARs are pure `N=` (ungapped) —

| transcript | UTA `uta_20210129` splign CIGAR | gap? |
|---|---|---|
| `NM_002386` | `3099=` | none (single ungapped block) |
| `NM_007199` | all `N=` (every exon block ungapped) | none |
| `NM_020451` | `238=` / `118=` (per-exon, ungapped) | none |

The open question PR3 left was whether the +2 came from a **CDS-coordinate
disagreement** (where each source places the CDS start, the 0-point of `c.`) —
which could have been a ferro CDS-boundary bug — or from a genomic alignment
difference. **A CDS-frame probe has now answered it: it is alignment skew, of a
genomic-anchor sub-kind, not a CDS bug.** The probe compared ferro's CDS start
against the UTA `cds_start_i` for each transcript and found they **agree**:

| transcript | ferro CDS start (1-based) | UTA `cds_start_i` (0-based) | agree? | first codon |
|---|---|---|---|---|
| `NM_020451` | 56 | 55 | yes (1-vs-0-based) | `ATG` |
| `NM_007199` | 103 | 102 | yes (1-vs-0-based) | `ATG` |
| `NM_002386` | 1381 | 1380 | yes (1-vs-0-based) | `ATG` |

So the `c.` 0-point is **not** in dispute — ferro and UTA annotate the same CDS
start, and the start codon is a clean `ATG` in every case. The per-exon CIGARs
being ungapped only tells us there is no *within-exon* `D`/`I` indel; it does not
fix where the transcript as a whole is **anchored on the genome**. The uniform
`+2` on every `c.` position is exactly that: a 2 bp difference in the genomic
anchor (the first aligned genomic position / `alt_start_i`) between ferro's
RefSeq-GFF alignment and the UTA splign alignment, applied uniformly to the whole
transcript. This is the same root cause as the gapped-CIGAR rows above — the two
sources align the transcript to the genome differently — just expressed as a
whole-transcript offset rather than a downstream-of-gap one.

This corrects PR3's reading. The earlier "ungapped → 0% divergence" gate is about
*within-exon* CIGAR structure; it does not constrain the genomic anchor, so an
ungapped transcript can still carry a whole-transcript genomic-anchor skew. There
is no CDS-boundary bug here (the CDS starts agree and are valid `ATG`), so this is
**not** the second-bug lead PR3 hypothesized. (The genuine second bug surfaced
elsewhere — the CDS-frame / reference-consistency issue tracked as #625/#629 — and
is independent of these three.)

These three are therefore **added to `SOURCE_SKEW_TRANSCRIPTS`** as the
genomic-anchor sub-kind, so a `CoordinateSkew` on them is auto-quarantined as
`divergence_accepted` (cluster `alignment-source-skew`) like the gapped rows. With
this, the dashboard reaches **0 residual FAIL** — every divergence is now either a
classified data-source skew, a tracked improvement, a tracked bug, or an accepted
convention divergence.

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
