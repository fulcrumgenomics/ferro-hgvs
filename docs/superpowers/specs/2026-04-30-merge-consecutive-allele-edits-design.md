# Merge consecutive edits within alleles per HGVS spec

Issue: [#72](https://github.com/fulcrumgenomics/ferro-hgvs/issues/72)

## Problem

`normalize()` currently preserves the input form for an `AlleleVariant` whose
sub-variants describe changes to the same or nearby nucleotides. The HGVS
spec requires such changes to be expressed as a single combined operation
rather than as multiple bracketed edits. This document specifies a merging
pass that runs as part of `normalize_allele` to bring allele output into
spec compliance.

## Spec basis

From the HGVS recommendations on delins
(<https://hgvs-nomenclature.org/stable/recommendations/DNA/delins/>):

- **Consecutive nucleotide changes must collapse to a single delins.** Notes:
  *"changes involving two or more consecutive nucleotides are described as
  deletion/insertion (delins) variants."* The Discussion FAQ reinforces this
  with: di-nucleotide substitutions don't exist — the change `GC` to `TG` at
  positions 4-5 must be `g.4_5delinsTG`, not two adjacent SNVs.
- **Non-adjacent variants are described separately, with one exception.**
  Notes: *"two variants separated by one or more nucleotides should
  preferably be described individually and not as a 'delins'. Exception:
  two variants separated by one nucleotide, together affecting one amino
  acid, should be described as a 'delins'."* The Examples section makes the
  carve-out concrete: `c.[145C>T;147C>G]` is *not correct*; the required
  form is `c.145_147delinsTGG` because positions 145–147 fall in a single
  codon.
- **An insertion at the boundary of an adjacent change participates in the
  merge.** Discussion FAQ: `NM_007294.3:c.[2077G>A;2077_2078insTA]` is
  incorrect; the correct form is `NM_007294.3:c.2077delinsATA`. The
  substitution at 2077 and the insertion between 2077 and 2078 share the
  boundary at position 2077.

## Scope and phasing

To keep development and review tractable, the work is split into two phases,
each landing as a separate commit. Phase 1 covers strictly consecutive
merging and is sufficient to make several long-standing examples spec-compliant.
Phase 2 adds the codon-frame exception, which requires reference-provider
access for the unchanged middle nucleotide.

### Phase 1: strictly consecutive merging

Two sub-variants merge iff their anchors (defined below) are exactly
adjacent — `A.end + 1 == B.start` — with no intervening unchanged
nucleotide. Applies to Genome (`g.`), Cds (`c.`), Tx (`n.`), Rna (`r.`),
and Mt (`m.`) variants.

### Phase 2: codon-frame exception for `c.` variants

Two CDS sub-variants whose anchors are separated by exactly one nucleotide
merge when their positions fall within a single codon. The unchanged middle
nucleotide is fetched from the reference provider. If no reference provider
is available or the lookup fails, the variants pass through unchanged (with
a warning surfaced via the existing `NormalizationWarning` channel) — Phase
1 behavior, not an error. Phase 2 applies only to `c.` variants in the
exonic CDS region (positions ≥ 1, not in 5'UTR or 3'UTR, not intronic).

Whether the same exception applies to `r.` is left for a follow-up; the
spec example uses `c.` and `r.` codon framing aligns with `c.` only when
the RNA region is coding. This is captured as a deferred test, not a
deferred phase of the same issue.

## Common scope (both phases)

### Phase

- **Cis only.** Trans, Mosaic, Chimeric, and Unknown phases never merge —
  they encode distinct chromosomes or distinct cell populations rather than
  a single contiguous coding sequence.

### Variant types

- Sub-variants whose edit is `NaEdit` participate: `Genome`, `Cds`, `Tx`,
  `Rna`, `Mt`.
- Other variant types (`Protein`, `Circular`, `RnaFusion`, `NullAllele`,
  `UnknownAllele`) pass through unchanged.

### Edit types

| Edit | Participates? | Notes |
| --- | --- | --- |
| `Substitution`, `SubstitutionNoRef` | yes | Single-base anchor |
| `Deletion` | yes | Contributes ref-range, empty alt |
| `Delins` | yes | Ref-range + alt sequence; only when `InsertedSequence` is `Literal` |
| `Insertion` | yes | Boundary-only edit; only when `InsertedSequence` is `Literal` |
| `Duplication` | no | Spec describes duplications as a distinct event; this exclusion is a deliberate scope choice, not a spec mandate |
| `Inversion` | no | Spec describes inversions as a distinct event; the merge pass also does *not* try to detect "this delins is actually an inversion" — that would be a separate canonicalization |
| `Repeat`, `MultiRepeat`, `Identity`, `Conversion`, `Unknown`, `Methylation`, `CopyNumber`, `Splice`, `NoProduct`, `PositionOnly` | no | Don't fit the consecutive-nucleotide model |

**Non-`Literal` `InsertedSequence` is a hard exclusion.**
`InsertedSequence` admits non-literal forms (`Count`, `Range`, `Repeat`,
`SequenceRepeat`, `Complex`, `Named`, `Reference`, `PositionRange`,
`PositionRangeInv`, `Uncertain`). These cannot be safely concatenated into
a single delins. If any participant `Insertion` or `Delins` carries a
non-`Literal` payload, the merge for that pair is refused; both inputs
pass through unchanged.

### Edit certainty

- Only `Mu::Certain` edits merge. `Mu::Uncertain` (parens-wrapped) and
  `Mu::Unknown` sub-variants pass through.

### Position constraints

A position participates in merging only if it is:

- not unknown (`?`),
- not a special marker (`pter`/`qter`/`cen`),
- not intronic (per the existing `is_intronic()` predicate on `CdsPos`,
  `TxPos`, etc.),
- not crossing a coordinate-system boundary. For `c.` this means 5'UTR
  (`-N`), CDS (`1..=N`), and 3'UTR (`*N`) are mutually exclusive merge
  regions — there is no valid HGVS range syntax that spans these (e.g.,
  `c.-1_1` does not exist), so adjacency cannot be expressed even when
  positions are physically consecutive. The same logic applies to `n.`
  (upstream `-N` / transcript / downstream `*N`) and `r.` (UTR boundaries).

### Accession and variant type

Two sub-variants are merge-eligible only if they share both:

- the same accession string (`Accession::full()` equality), and
- the same `HgvsVariant` discriminant (`Cds` with `Cds`, `Genome` with
  `Genome`, etc.).

Gene symbol is not separately checked; in practice it is determined by the
accession. The merged result inherits the accession, gene symbol, and
discriminant of the first input.

## Algorithm

### Anchor model

Each merge-eligible sub-variant has an *anchor*: an inclusive position
range it occupies plus the alt content it contributes.

| Edit | Anchor start | Anchor end | Alt |
| --- | --- | --- | --- |
| `sub` at `p` (`X>Y`) | `p` | `p` | `Y` |
| `del` at `[p, q]` | `p` | `q` | `""` |
| `delins` at `[p, q]` (alt `S`) | `p` | `q` | `S` |
| `ins` between `p` and `p+1` (alt `S`) | `p+1` | `p` | `S` |

For an insertion, `start = end + 1` — the range is empty but locates the
boundary precisely. Two anchors `A` and `B` (in position order) are
*Phase-1 adjacent* iff `A.end + 1 == B.start`. This produces the right
answer for every combination:

- sub/del/delins at `[…, p]` then sub/del/delins at `[p+1, …]` → adjacent
- sub/del/delins at `[…, p]` then `ins` between `p` and `p+1` → adjacent
- `ins` between `p` and `p+1` then sub/del/delins at `[p+1, …]` → adjacent
- two `ins` at the same boundary `p|p+1` → adjacent
- `ins` between `p` and `p+1` then `ins` between `p+1` and `p+2` → not
  adjacent (one intervening unchanged nucleotide; spec rule keeps separate)

Two anchors are *Phase-2 adjacent* iff they are CDS exonic SNVs (not
ranges, not insertions, both `Substitution` or `SubstitutionNoRef`, both
position ≥ 1 with no offset and no UTR), with `A.end + 2 == B.start`, and
`(A.end - 1) / 3 == (B.start - 1) / 3` (same codon, 0-indexed). Phase 2
extends only the SNV-pair case; broader forms (`sub`+`ins`+`sub` straddling
a codon, etc.) are not in scope for issue #72.

### Walk

In `normalize_allele`, after the existing `resolve_overlaps` call:

1. Return the input unchanged unless `phase == Cis` and `len() >= 2`.
2. Walk left-to-right with a single output cursor, preserving input order
   (input order should already match position order; no resort, so
   ordering is stable for tie-cases like two insertions at the same
   boundary).
3. For each next variant, attempt to merge with the head of the output
   under Phase-1 rules. If that fails and Phase 2 is in scope, attempt the
   codon-frame merge. Otherwise push.
4. Sub-variants that aren't merge-eligible (different accession, different
   type, non-NaEdit, uncertain edit, disqualifying position, non-Literal
   insertion payload, etc.) act as merge barriers — they're emitted in
   place and reset the cursor.

### Merge result

For a merged group spanning anchor `[start, end]` with concatenated alt
`alt_seq`:

- If `start > end` (all inputs were insertions at one boundary) → emit
  `Insertion` between `start - 1` and `start` with content `alt_seq`.
- Else if `alt_seq.is_empty()` (only deletions contributed) → emit
  `Deletion` at `[start, end]`. Drop any per-input ref-sequence; emit the
  no-sequence form (`g.<start>_<end>del`), which is the spec-recommended
  shape and avoids the ambiguity of partial ref-sequence info.
- Else → emit `Delins` at `[start, end]` with
  `InsertedSequence::Literal(alt_seq)`.

Phase-2 codon-frame merges always produce the `Delins` branch. The middle
nucleotide is fetched from the reference provider and embedded in
`alt_seq` between the contributions of the two SNVs.

### Two insertions at the same boundary

`g.[100_101insT;100_101insA]` → `g.100_101insTA`. Order is the input order
of the allele bracket. Two insertions at the *same* boundary merge into a
single `Insertion`; the `Delins` branch is not used because there is no
covered range.

## Placement

The merge pass lives in a new module `src/normalize/merge.rs`, called from
`normalize_allele` in `src/normalize/mod.rs` after the existing
`resolve_overlaps` call.

Note that `resolve_overlaps` is currently a no-op for the variants
themselves — it walks for overlap detection but does not modify the input
(see the comment block in `mod.rs` near the end of `resolve_overlaps`).
The merge pass therefore sees the post-individual-normalization positions
directly. As a defensive measure, the merge pass refuses to merge any pair
where `A.end >= B.start` (true overlap), so that an unresolved overlap
silently degrades to "leave both inputs as-is" rather than producing a
malformed merged variant.

The new module exposes:

```rust
pub(super) fn merge_consecutive_edits(
    variants: &[HgvsVariant],
    phase: AllelePhase,
    provider: Option<&dyn ReferenceProvider>, // None disables Phase 2
) -> Vec<HgvsVariant>;
```

Private helpers handle anchor extraction, adjacency checks, and merged-edit
construction. The position-typed pieces (`GenomePos`, `CdsPos`, `TxPos`,
`RnaPos`) are bridged via small per-type helpers; the core walk is shared.

## Tests

All assertions use `parse → normalize → display` round-trip and check the
rendered HGVS string, which is the user-visible contract.

### Phase 1 tests

#### Issue examples

- `g.[1000G>A;1001A>C]` → `g.1000_1001delinsAC`
- `g.[1000del;1001del]` → `g.1000_1001del`

#### Coordinate systems

One sub+sub case per merge-eligible coordinate system (g/c/n/r/m). RNA
case verifies lowercase nucleotide preservation in the merged alt
sequence.

#### Spec FAQ example (sub+ins)

- `NM_007294.3:c.[2077G>A;2077_2078insTA]` → `NM_007294.3:c.2077delinsATA`
- Mirror: `c.[2077_2078insTA;2078G>A]` (ins-then-sub at the shared
  boundary `2077|2078`).

#### Mixed mergeable types

- `sub+del`, `del+sub` → delins
- `sub+delins`, `delins+sub` → delins
- `del+delins`, `delins+del`, `delins+delins`
- `del+ins` and `ins+del` at a shared boundary

#### Chains

- Three consecutive substitutions → single delins
- `sub`, `ins`, `sub` chained at a single shared boundary triple → single
  delins
- Long chain (≥ 4) of mixed sub/del/ins → single delins

#### Two insertions at the same boundary

- `g.[100_101insT;100_101insA]` → `g.100_101insTA` (preserve input order)

#### Deletion ref-sequence handling

- `g.[100delA;101delC]` → `g.100_101del` (sequence dropped)
- `g.[100del;101del]` → `g.100_101del`
- `g.[100del3;101del2]` → `g.100_101del` (length info dropped; emit
  no-sequence form)

#### Round-trips that must remain unchanged

- One-nucleotide gap between variants in non-CDS context (`g.[100G>A;102C>T]`)
- Same position twice (zero-gap but overlap, not adjacency)
- Different accessions in the bracket
- Different variant-type discriminants in the bracket
- Trans, Mosaic, Chimeric, Unknown phase
- Sub at intronic position (`c.79+1A>G` adjacent to `c.79+2T>G`)
- Sub crossing 5'UTR/CDS boundary (`c.[-1A>G;1G>T]`)
- `Mu::Uncertain` edit
- Duplication adjacent to substitution
- Inversion adjacent to substitution
- Two insertions at *different* boundaries (`g.[100_101insT;101_102insA]`)
- `Delins` with non-`Literal` `InsertedSequence` adjacent to `sub`
- Reverse-input-order pair that would not be position-order-adjacent
  (`g.[1001A>C;1000G>A]` — input order preserved, no merge in Phase 1
  since the walk checks position-order adjacency; this pins the no-resort
  decision)

### Phase 2 tests

- Spec example: `NM_…:c.[145C>T;147C>G]` (with a reference provider that
  returns `G` at position 146) → `c.145_147delinsTGG`
- Same-codon SNV pair, position 1+3 of a codon, exhaustive over the codon
  positions (1+3, 2+? — Phase 2 only handles the 1+3 / 1-intervening case)
- Cross-codon SNV pair (`c.[3C>T;5G>A]` straddles codon 1/codon 2) → no
  merge
- One-nt-gap pair in `g.` context → no merge (Phase 2 is `c.` only)
- One-nt-gap pair with no reference provider → no merge, warning
  surfaced
- One-nt-gap pair with reference-provider failure → no merge, warning
  surfaced
- Phase-1 merge still wins for strictly consecutive `c.` SNVs (Phase 2
  rule is only consulted when Phase 1 doesn't apply)

### Deferred (follow-up issue, not in this PR)

- Codon-frame exception for `r.` coding regions
- Codon-frame exception for chains involving more than two SNVs (`sub`,
  `ins`, `sub` straddling a codon)
- Codon-frame exception for SNV+`del` pairs separated by one unchanged
  nucleotide

These get a single `#[ignore]`-d test each so the gap is visible in the
test suite rather than only in this design doc.
