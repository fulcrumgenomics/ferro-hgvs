//! A cancelled identity member must not survive beside a sibling that claims
//! its bases (#1416).
//!
//! `g.[11_12inv;11_12insAA]` normalized to `g.[11_12=;10_11A[4]]`, whose two
//! members both claim base 11 — so the allele denotes **no sequence at all**,
//! and the independent applier declines it in either member order. The whole
//! variant is the repeat: `11_12inv` over `AT` is its own reverse complement and
//! cancels to `=`, while `11_12insAA` grows the `A` tract at `10_11` from two
//! copies to four. `g.10_11A[4]` alone denotes exactly the input's bases.
//!
//! This is #1235's criterion 2 — "no normalized cis allele contains overlapping
//! or out-of-order members" — on a shape its confluence proptest cannot
//! generate, since that model never puts an inversion beside a co-located
//! insertion.
//!
//! Pinned with the apply oracle rather than by re-normalizing, deliberately: the
//! bad output also fires `FERRO_ASSERT_IDEMPOTENT` (#1414), so a
//! re-normalization test would report the ordering symptom instead of the
//! defect. What matters is that the output denotes something.

use crate::common::cis_apply_oracle::apply;
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};

/// 18 bases; `10` and `11` are `A`, `12` is `T`, so `11_12` is `AT` — its own
/// reverse complement — and the `A` tract at `10_11` has two copies.
const REFERENCE: &str = "CGCGCGCGCAATCGCGCG";

fn normalize(input: &str) -> String {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("TEMPLATE", REFERENCE.to_string());
    Normalizer::new(provider)
        .normalize(&parse_hgvs(input).expect("parse input"))
        .expect("normalize")
        .to_string()
}

/// The output must denote a sequence.
///
/// The primary assertion: a description the applier declines is malformed in the
/// strongest sense available, whatever else is true of it.
#[test]
fn the_normalized_allele_denotes_a_sequence() {
    let out = normalize("TEMPLATE:g.[11_12inv;11_12insAA]");
    // The whole accession-qualified string: `apply` re-parses it with
    // `parse_hgvs`, which needs the accession. Passing the bare `g.…` makes it
    // return `None` for *every* input, which looks exactly like the defect.
    assert!(
        apply(REFERENCE, &out).is_some(),
        "normalized to a description the applier declines — its members overlap, \
         so it denotes no sequence: {out}"
    );
}

/// …and it must denote the *right* one.
///
/// Stated separately because `is_some()` alone would accept any appliable
/// answer, including one that quietly changed the variant's meaning while
/// fixing the overlap.
#[test]
fn the_normalized_allele_denotes_the_inputs_bases() {
    let out = normalize("TEMPLATE:g.[11_12inv;11_12insAA]");
    // The whole accession-qualified string: `apply` re-parses it with
    // `parse_hgvs`, which needs the accession. Passing the bare `g.…` makes it
    // return `None` for *every* input, which looks exactly like the defect.
    // `11_12insAA` grows the two-copy `A` tract at `10_11` to four copies; the
    // inversion cancels. Spliced by hand so the expectation does not depend on
    // the normalizer that produced the description.
    assert_eq!(
        apply(REFERENCE, &out).as_deref(),
        Some("CGCGCGCGCAAAATCGCGCG"),
        "normalized description denotes different bases than the input: {out}"
    );
}

/// A reversed `<high>_<low>` span must not be read as an ordinary interval.
///
/// SVD-WG006 admits the reversed range for a circular deletion or duplication,
/// and `cis_axis_parts` passes the endpoints through raw, so a member with
/// `start > end` reaches the coverage test above. Widening containment to
/// overlap is what made that matter: `s.start <= span.end && s.end >= span.start`
/// is satisfied by a reversed sibling against a forward identity, where the
/// containment form refused it — so without the forward-span guard this rule
/// would start dropping identity members on the circular axis on the strength of
/// an interval comparison that means nothing there.
///
/// The assertion is deliberately weak on the output *form* and strong on the
/// invariant: whatever `m.` normalization does with this allele, it must not
/// silently lose a member to a wraparound comparison.
#[test]
fn a_wraparound_sibling_does_not_drop_an_identity_member() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NC_012920.1", "CGCGCGCGCAATCGCGCG".to_string());
    let input = "NC_012920.1:m.[16569_1del;5_6=]";
    let Ok(variant) = parse_hgvs(input) else {
        return; // the parser is entitled to refuse this spelling; nothing to pin
    };
    let output = Normalizer::new(provider)
        .normalize(&variant)
        .expect("lenient normalization must not reject")
        .to_string();
    assert!(
        output.contains("5_6=") || output.contains("5="),
        "`{input}` -> `{output}`: the identity member was dropped by a \
         wraparound span comparison"
    );
}
