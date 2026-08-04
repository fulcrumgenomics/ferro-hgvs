//! #1362 — an identity edit's SPDI triple must span the bases it claims.
//!
//! `hgvs_to_spdi` used to convert an identity whose bases the input did not
//! spell (`g.100=`) into a zero-width triple that claimed nothing. That is not
//! merely lossy: a zero-width triple at an interbase is indistinguishable from
//! an *insertion* at that interbase, so an applier walking an allele's members
//! cannot tell `g.[261_262dup;263=]` apart from two insertions competing for
//! one position and has to decline the whole description.
//!
//! The unit tests in `src/spdi/convert.rs` cover the conversion arm directly.
//! This file covers it through the public API at the two places an
//! inclusive-vs-exclusive coordinate error would surface — the **last base of
//! the contig** and one past it — because the fetch converts HGVS's 1-based
//! inclusive end into the provider's own end convention, and an off-by-one
//! there is invisible on an interior span.

use ferro_hgvs::spdi::convert::{hgvs_to_spdi, hgvs_to_spdi_simple, spdi_to_hgvs};
use ferro_hgvs::{parse_hgvs, ConversionError, MockProvider, SpdiVariant};

/// A 28-base `ACGT…` contig: 1-based position `p` holds `"ACGT"[(p - 1) % 4]`,
/// so the final base (28) is `T` and 27..28 is `GT`.
const CONTIG: &str = "ACGTACGTACGTACGTACGTACGTACGT";
const CONTIG_LEN: u64 = 28;

fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NC_000001.11", CONTIG.to_string());
    provider
}

/// The identity that ends exactly on the contig's last base must fetch that
/// base, not one past it and not one short of it.
#[test]
fn a_terminal_identity_claims_the_last_base() {
    assert_eq!(CONTIG.len() as u64, CONTIG_LEN, "contig length assumption");
    let spdi = hgvs_to_spdi(&parse_hgvs("NC_000001.11:g.28=").unwrap(), &provider()).unwrap();
    assert_eq!(spdi.position, 27, "0-based interbase of 1-based base 28");
    assert_eq!(spdi.deletion, "T", "base 28 of an ACGT tandem");
    assert_eq!(spdi.insertion, "T", "an identity changes nothing");
}

/// A multi-base identity ending on the last base pins the inclusive end
/// against the provider's exclusive one: an off-by-one yields `"G"` (end
/// treated as exclusive) or an out-of-range error (end over-extended).
#[test]
fn a_terminal_identity_span_claims_both_of_its_bases() {
    let spdi = hgvs_to_spdi(&parse_hgvs("NC_000001.11:g.27_28=").unwrap(), &provider()).unwrap();
    assert_eq!(spdi.position, 26);
    assert_eq!(
        spdi.deletion, "GT",
        "HGVS 27_28 is inclusive of both ends, so the span is two bases"
    );
    assert_eq!(spdi.insertion, "GT");
}

/// One past the end must be reported, not clamped to the last base.
#[test]
fn an_identity_past_the_contig_end_is_reported() {
    let result = hgvs_to_spdi(&parse_hgvs("NC_000001.11:g.29=").unwrap(), &provider());
    assert!(
        matches!(result, Err(ConversionError::MissingReferenceData { .. })),
        "expected MissingReferenceData for position 29 on a 28-base contig, got {result:?}"
    );
}

/// The same span through the provider-less entry point reports rather than
/// inventing an empty span. This is the breaking half of #1362.
#[test]
fn the_provider_less_path_reports_instead_of_claiming_nothing() {
    let result = hgvs_to_spdi_simple(&parse_hgvs("NC_000001.11:g.28=").unwrap());
    assert!(
        matches!(result, Err(ConversionError::MissingReferenceData { .. })),
        "expected MissingReferenceData without a provider, got {result:?}"
    );
}

/// A whole-entity `=` names no interval, so it has no triple at all — rather
/// than one at position 0, which read back as `g.1=` and narrowed a statement
/// about the whole sequence to one about its first base.
#[test]
fn a_whole_entity_identity_has_no_triple() {
    let result = hgvs_to_spdi(&parse_hgvs("NC_000001.11:g.=").unwrap(), &provider());
    assert!(
        matches!(result, Err(ConversionError::UnsupportedEditType { .. })),
        "expected UnsupportedEditType for a whole-entity identity, got {result:?}"
    );
}

/// A zero-width triple names no bases, so it cannot come back as an identity
/// over one. `spdi_to_hgvs` used to read `NC_000001.11:99::` as `g.100=`, an
/// identity asserting a base the triple never claimed — the same
/// invent-a-base error as #1362's forward direction, in the reverse direction.
///
/// The triple is built directly rather than parsed so the case survives a
/// parser that later rejects the spelling: the public constructors and
/// `SpdiVariant`'s own fields admit it regardless, and this arm is what must
/// refuse it.
#[test]
fn a_zero_width_triple_is_not_an_identity() {
    let spdi = SpdiVariant::new("NC_000001.11", 99, "", "");
    assert!(
        spdi.is_identity(),
        "deletion == insertion, so it takes the identity arm"
    );
    let result = spdi_to_hgvs(&spdi);
    assert!(
        matches!(result, Err(ConversionError::UnsupportedEditType { .. })),
        "a triple naming zero bases must be refused, not read as `g.100=`, got {result:?}"
    );
}

/// Round-tripping a multi-base identity keeps its span. Before #1362 the
/// triple was zero-width, so the reverse direction never saw a multi-base
/// identity from this input at all.
#[test]
fn a_multi_base_identity_keeps_its_span_through_spdi() {
    let spdi = hgvs_to_spdi(&parse_hgvs("NC_000001.11:g.27_28=").unwrap(), &provider()).unwrap();
    let back = spdi_to_hgvs(&spdi).expect("round trip");
    assert_eq!(
        back.to_string(),
        "NC_000001.11:g.27_28GT=",
        "a two-base identity must come back as a two-base interval, not a point"
    );
}
