<!-- GENERATED — do not edit; regenerate via `cargo run --features dev --example generate_conformance_summary` -->

# hgvs-rs-projection conformance summary

Generated from `cases.json` — do not edit by hand; regenerate with the example above. Every row below is a tracked disposition. The live divergence set against full reference data is emitted only by the manifest run and is never committed (it is non-hermetic).

## Root-cause clusters

### ferro RefSeq-GFF alignment vs 2021 UTA splign snapshot

Spec: `recommendations/general/refseq.md`

ferro projects on NCBI-canonical RefSeq-GFF alignments (cdot-0.2.32.refseq); the corpus expectations come from the 2021 biocommons UTA snapshot (uta_20210129, splign). Where the two alignments disagree (within-exon D/I gaps or exon-boundary off-by-1) ferro returns the expected base transcript at a different coordinate. Auto-quarantined structurally for the gate-listed transcripts only; see tests/fixtures/hgvs-rs-projection/ALIGNMENT_SOURCE_DIVERGENCE.md.

_No seeded member rows yet (manifest-gated seeding, #325)._

### expected transcript absent/not-selected in ferro cdot-0.2.32.refseq vs UTA

Spec: `recommendations/general/refseq.md`

ferro does not return the expected base transcript at all (selection-coverage miss): the expected accession is absent from, or not prioritized by, ferro cdot-0.2.32.refseq relative to the 2021 UTA snapshot. Auto-quarantined structurally when no returned consequence carries the expected base accession.

_No seeded member rows yet (manifest-gated seeding, #325)._

### ferro emits the spec-preferred parenthesized predicted protein form (p.(X)) vs hgvs-rs bare p.X

Spec: `recommendations/protein/substitution.md`

On a (c., p.) pair where the coding component matches version-insensitively, ferro renders the DNA-predicted protein consequence parenthesized (p.(Arg268Trp)) per the HGVS predicted-consequence rule (a consequence predicted from a DNA change is uncertain and rendered p.(...); docs/recommendations/protein/substitution.md). hgvs-rs carries the bare legacy form (p.Arg268Trp). ferro is spec-correct, so these are routed structurally to `improvement` in tests/hgvs_rs_projection_tests.rs (is_predicted_parenthesization_pairs) rather than FAIL. #615-related: the same predicted-consequence machinery now renders the C-terminal stop-loss extension correctly after the #615 fix (PR #621).

_No seeded member rows yet (manifest-gated seeding, #325)._

### ferro emits the spec-preferred three-letter extension glyph (extTer) vs hgvs-rs legacy ext*

Spec: `background/standards.md`

For a C-terminal stop-loss extension, ferro canonicalizes to the three-letter form p.(...extTer<n>) (#224); the corpus carries the legacy single-character p.(...ext*<n>). Three-letter codes (including Ter) are the preferred canonical form, with `*` an accepted alternative (background/standards.md). Both are valid HGVS, so these rows are dispositioned per-case as `improvement` (#224), XPASS-guarded.

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `NM_000051.3:c.9170_9171delGA` | protein_description | spec_citation | — | — |
| `NM_000249.3:c.2266_2269dupTGTT` | protein_description | spec_citation | — | — |
| `NM_000249.3:c.2269dupT` | protein_description | spec_citation | — | — |
| `NM_000425.3:c.3772dupT` | protein_description | spec_citation | — | — |
| `NM_000488.3:c.1391dupA` | protein_description | spec_citation | — | — |
| `NM_152263.2:c.855delA` | protein_description | spec_citation | — | — |

### whole-CDS-deletion protein form: corpus p.0? vs ferro p.(Met1?)

Spec: `recommendations/protein/deletion.md`

For a deletion spanning the entire coding sequence (including the initiation codon), the corpus emits p.0? (no protein produced) while ferro emits p.(Met1?) (uncertain start-loss). Both are spec-allowed predicted forms for a variant removing the start codon; ferro reports the start-loss form it derives from the normalized variant. Terminal convention difference (accepted_divergence, policy ferro-policy-whole-cds-del-met1), not a bug; p.0? is arguably preferable as a future refinement.

| input | axis | disposition | ferro output | tracking |
|---|---|---|---|---|
| `NM_000249.3:c.-7_*46del` | protein_description | accepted_divergence | — | — |

## Disposition tallies

| axis | accepted_divergence | known_bug | improvement | reference_unavailable | spec_citation |
|---|---:|---:|---:|---:|---:|
| protein_description | 1 | 0 | 0 | 0 | 6 |
