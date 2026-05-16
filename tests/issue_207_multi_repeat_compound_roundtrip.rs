//! Audit (item B3 of tracking issue #81 / closes #207): multi-repeat
//! compound rendering `[A];[B]` round-trip.
//!
//! HGVS DNA `repeated.md` example (spec line cited in the tests):
//!
//! ```
//! NC_000014.8:g.[101179660_101179695TG[14]];[101179660_101179695TG[18]]
//! ```
//!
//! A trans-phase allele where each allele carries the same repeat locus
//! at a different copy count. `[X];[Y]` is the spec's canonical form for
//! "X on one allele, Y on the other"; the two alleles are siblings in
//! trans phase, not a cis bracket.
//!
//! Goal: pin parse + normalize round-trip across every coord system that
//! admits the repeat shape (`g.`, `c.`, `n.`, `r.`, `m.`), and cover the
//! single-allele baseline plus the cis-phase form `[X; Y]` as a foil so
//! a regression flips one without the others.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

fn round_trip(input: &str) -> String {
    let normalizer = Normalizer::new(MockProvider::new());
    let variant =
        parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for {:?}: {:?}", input, e));
    let normalized = normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize failed for {:?}: {:?}", input, e));
    format!("{}", normalized)
}

// =============================================================================
// Spec textbook example — genomic trans-phase compound repeat.
// =============================================================================

#[test]
fn genomic_trans_compound_repeat_spec_example() {
    // `assets/hgvs-nomenclature/docs/recommendations/DNA/repeated.md`:
    //   NC_000014.8:g.[101179660_101179695TG[14]];[101179660_101179695TG[18]]
    // "a repeated TG di-nucleotide sequence, starting at position
    //  g.101179660 on human chromosome 14, is present with 14 TG copies on
    //  one allele and 18 TG copies on the other allele."
    let s = "NC_000014.8:g.[101179660_101179695TG[14]];[101179660_101179695TG[18]]";
    assert_eq!(round_trip(s), s);
}

#[test]
fn genomic_single_allele_repeat_baseline() {
    // The single-allele baseline (no bracket wrapper) must round-trip
    // unchanged — if this fails the compound case can't possibly work.
    let s = "NC_000014.8:g.101179660_101179695TG[14]";
    assert_eq!(round_trip(s), s);
}

#[test]
fn genomic_cis_compound_repeat_at_two_loci() {
    // Cis form: a single allele bracket carrying two repeat edits at
    // distinct loci. Different from B3's trans case but worth pinning so
    // a parser regression that broke the trans bracket separator also
    // shows up against the cis form.
    let s = "NC_000001.11:g.[100_105TG[3];200_205AC[4]]";
    assert_eq!(round_trip(s), s);
}

// =============================================================================
// Per-coord-system trans-phase compound repeat round-trips.
// =============================================================================

#[test]
fn cds_trans_compound_repeat() {
    let s = "NM_000088.3:c.[10_15TG[3]];[10_15TG[5]]";
    assert_eq!(round_trip(s), s);
}

#[test]
fn noncoding_trans_compound_repeat() {
    let s = "NR_046018.2:n.[10_15TG[3]];[10_15TG[5]]";
    assert_eq!(round_trip(s), s);
}

#[test]
fn rna_trans_compound_repeat_lowercase() {
    // RNA spec requires lowercase bases; the compound trans form must
    // preserve case end-to-end.
    let s = "NM_000088.3:r.[10_15ug[3]];[10_15ug[5]]";
    assert_eq!(round_trip(s), s);
}

#[test]
fn mt_trans_compound_repeat() {
    let s = "NC_012920.1:m.[100_105TG[3]];[100_105TG[5]]";
    assert_eq!(round_trip(s), s);
}

// =============================================================================
// Range counts and uncertainty inside compound repeats.
// =============================================================================

#[test]
fn genomic_trans_compound_repeat_range_count() {
    // `TG[14_18]` carries an uncertain `(min_max)` count. Each allele in
    // the trans pair may itself express uncertainty independently.
    let s = "NC_000014.8:g.[101179660_101179695TG[14_18]];[101179660_101179695TG[20]]";
    assert_eq!(round_trip(s), s);
}

#[test]
fn genomic_trans_compound_repeat_unknown_count() {
    // `TG[?]` is the spec's "present, count unknown" shape; the compound
    // round-trip must preserve it on either allele.
    let s = "NC_000014.8:g.[101179660_101179695TG[?]];[101179660_101179695TG[18]]";
    assert_eq!(round_trip(s), s);
}

// =============================================================================
// Single-unit repeats inside the compound trans form.
// =============================================================================

#[test]
fn genomic_trans_compound_single_base_repeat() {
    // Single-base repeat (`A[N]`) inside the trans compound. The repeat
    // alphabet is intentionally narrow — `A` vs `T` vs `C` vs `G` — so
    // any single-letter parser regression shows up here.
    let s = "NC_000001.11:g.[100A[6]];[100A[8]]";
    assert_eq!(round_trip(s), s);
}

// =============================================================================
// Structural variants — distinct loci / distinct units. These must
// round-trip without the trans bracket quietly normalising the two
// halves toward each other.
// =============================================================================

#[test]
fn genomic_trans_compound_repeat_distinct_loci() {
    // Two alleles at different loci, each its own repeat. Spec admits
    // this as a valid trans description; ferro must round-trip both
    // halves without merging or reordering.
    let s = "NC_000001.11:g.[100_105TG[3]];[200_205AC[5]]";
    assert_eq!(round_trip(s), s);
}

#[test]
fn genomic_trans_compound_repeat_distinct_units() {
    // Two alleles at the same locus but different repeat units (TG vs
    // GT). Round-trip must preserve the units verbatim.
    let s = "NC_000001.11:g.[100_105TG[3]];[100_105GT[5]]";
    assert_eq!(round_trip(s), s);
}

// =============================================================================
// Coverage for the remaining parsers touched by the fix.
// =============================================================================

#[test]
fn circular_trans_compound_repeat() {
    // The `o.` (circular genome) trans-allele parser
    // (`parse_circular_trans_allele_shorthand`) is part of the
    // 8-site fix; lock in a regression test against the same
    // `[A];[B]` shape with nested repeat counts.
    let s = "NC_001422.1:o.[100_105TG[3]];[100_105TG[5]]";
    assert_eq!(round_trip(s), s);
}

#[test]
fn top_level_trans_compound_repeat_fully_qualified() {
    // The top-level `parse_trans_allele` (driven via
    // `detect_allele_type` for ClinVar-style fully-qualified shapes
    // where every member carries its own `ACC:coord.edit` prefix) is
    // the eighth call site touched by the fix and deserves its own
    // direct regression test. Each member shares an accession so
    // Display normalizes to the compact-prefix `ACC:g.[X];[Y]` form;
    // re-parsing that canonical output must round-trip stably.
    let fully_qualified =
        "[NC_000014.8:g.101179660_101179695TG[14]];[NC_000014.8:g.101179660_101179695TG[18]]";
    let canonical = "NC_000014.8:g.[101179660_101179695TG[14]];[101179660_101179695TG[18]]";
    assert_eq!(round_trip(fully_qualified), canonical);
    assert_eq!(round_trip(canonical), canonical);
}

#[test]
fn delins_reference_with_nested_repeat_count_round_trip() {
    // Companion fix in `parser/edit.rs::parse_bracketed_inserted_sequence`:
    // a `delins[…]` whose bracketed payload is a reference range that
    // *itself* contains an `[N]` repeat-count bracket previously had
    // its close-bracket lookup truncated. The reference payload
    // (e.g. `NC_000022.11:g.100_200TG[3]`) is now scanned with the
    // depth-aware variant of the helper used by the trans-allele
    // parsers above.
    let s = "NC_000001.11:g.100_105delins[NC_000022.11:g.100_200TG[3]]";
    assert_eq!(round_trip(s), s);
}

// =============================================================================
// Trans-allele branches whose sibling is a null/unknown marker.
// =============================================================================

#[test]
fn trans_compound_repeat_alongside_null_allele() {
    // `[0]` is the spec's "no allele" marker. The trans helper has a
    // dedicated branch for it (`content == "0"`); confirm it still
    // works when the *other* allele contains a nested repeat
    // count's `[N]`. The Display path falls back to the
    // fully-qualified form `[ACC:variant];[0]` (compact-prefix is
    // only used when both alleles are real variants sharing an
    // accession).
    let input = "NC_000001.11:g.[100_105TG[3]];[0]";
    let canonical = "[NC_000001.11:g.100_105TG[3]];[0]";
    assert_eq!(round_trip(input), canonical);
    assert_eq!(round_trip(canonical), canonical);
}

#[test]
fn trans_compound_repeat_alongside_unknown_allele() {
    // `[?]` is the spec's "unknown allele" marker. Same Display shape
    // as the null-allele test above.
    let input = "NC_000001.11:g.[100_105TG[3]];[?]";
    let canonical = "[NC_000001.11:g.100_105TG[3]];[?]";
    assert_eq!(round_trip(input), canonical);
    assert_eq!(round_trip(canonical), canonical);
}

// =============================================================================
// Empty-member rejection — `[];[A]` is malformed.
// =============================================================================

#[test]
fn empty_first_member_is_rejected() {
    // The trans helper accepts a content of `"0"` or `"?"` (the null
    // and unknown markers) and otherwise routes to the coord-system
    // variant parser; an empty content fails fast inside that parser.
    // Locks in the rejection so a future "fix" that silently accepts
    // it shows up as a regression.
    let normalizer = Normalizer::new(MockProvider::new());
    let result = parse_hgvs("NC_000001.11:g.[];[100_105TG[3]]");
    assert!(
        result.is_err(),
        "expected parse error for empty first member, got: {:?}",
        result.map(|v| normalizer.normalize(&v).map(|n| n.to_string()))
    );
}
