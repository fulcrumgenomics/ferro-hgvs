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

## Disposition tallies

| axis | accepted_divergence | known_bug | improvement | spec_citation |
|---|---:|---:|---:|---:|
