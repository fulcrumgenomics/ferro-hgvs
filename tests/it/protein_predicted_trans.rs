//! Trans protein alleles whose members are predicted-change wrappers (#544).
//!
//! `p.[(Ser68Arg)];[(Ser73Arg)]` — each trans member is a predicted protein
//! change `(...)`. Previously the protein trans bracket-member parser only
//! handled bare members and whole-entity markers, not the predicted wrapper.

use ferro_hgvs::hgvs::AllelePhase;
use ferro_hgvs::{parse_hgvs, HgvsVariant};

#[test]
fn protein_predicted_trans_round_trips() {
    for s in [
        "NP_003997.1:p.[(Ser68Arg)];[(Ser73Arg)]",
        "NP_003997.1:p.[(Ser73Arg)];[(Asn103del)]",
    ] {
        let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
        let HgvsVariant::Allele(a) = &v else {
            panic!("expected Allele for `{s}`, got {v:?}");
        };
        assert_eq!(a.phase, AllelePhase::Trans, "phase for `{s}`");
        assert_eq!(a.variants.len(), 2, "member count for `{s}`");
        assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
    }
}

/// A predicted member paired with the unknown-allele marker `[?]` must at
/// least parse to a trans allele (its compact round-trip is tracked with the
/// other trans `[?]` rendering work).
#[test]
fn protein_predicted_trans_with_unknown_member_parses() {
    let v = parse_hgvs("NP_003997.1:p.[(Ser68Arg)];[?]").expect("must parse");
    let HgvsVariant::Allele(a) = &v else {
        panic!("expected Allele, got {v:?}");
    };
    assert_eq!(a.phase, AllelePhase::Trans);
    assert_eq!(a.variants.len(), 2);
    assert!(matches!(a.variants[1], HgvsVariant::UnknownAllele));
}
