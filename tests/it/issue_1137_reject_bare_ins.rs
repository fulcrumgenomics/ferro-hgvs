//! Issue #1137: the parser must reject an insertion with no inserted sequence.
//!
//! An insertion is defined by what it inserts (`DNA/insertion.md:22`): every
//! form the spec offers supplies a payload — literal bases (`insAGC`), a
//! same-reference range (`ins858_895`), or another reference. There is no
//! production for a bare `ins` in `syntax.yaml`'s `dna.ins`. Before this fix
//! the parser accepted `…ins` and round-tripped it unchanged, minting a
//! `NaEdit::Insertion { sequence: Empty }` — a semantically empty edit that
//! asserts no change and forces downstream code to defend against it (see the
//! internally-constructed empty-insertion guard from #1136).
//!
//! A bare `ins` has no repair (the bases are simply absent), so it is rejected
//! in every error mode — standard, lenient, and strict — at parse time. These
//! tests need no prepared reference, so they gate every PR.

use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::{parse_hgvs_lenient, parse_hgvs_with_config};
use ferro_hgvs::parse_hgvs;

/// Bare-`ins` inputs across coordinate systems, each a two-position anchor with
/// the payload omitted.
const BARE_INS: &[&str] = &[
    "NC_000001.11:g.100_101ins",
    "NM_004006.3:c.179_180ins",
    "NM_000088.3:n.100_101ins",
    "NM_004006.2:r.100_101ins",
    "NC_012920.1:m.100_101ins",
];

/// The default (standard) parse path must reject a bare `ins`.
#[test]
fn standard_parse_rejects_bare_ins() {
    for input in BARE_INS {
        assert!(
            parse_hgvs(input).is_err(),
            "expected a bare `ins` to be rejected in standard mode: {input:?}"
        );
    }
}

/// Lenient mode has no repair for a missing inserted sequence, so it rejects
/// rather than round-tripping the empty edit.
#[test]
fn lenient_parse_rejects_bare_ins() {
    for input in BARE_INS {
        assert!(
            parse_hgvs_lenient(input).is_err(),
            "expected a bare `ins` to be rejected in lenient mode: {input:?}"
        );
    }
}

/// Strict mode rejects a bare `ins`.
#[test]
fn strict_parse_rejects_bare_ins() {
    for input in BARE_INS {
        assert!(
            parse_hgvs_with_config(input, ErrorConfig::strict()).is_err(),
            "expected a bare `ins` to be rejected in strict mode: {input:?}"
        );
    }
}

/// Regression guard: insertions that DO carry a payload still parse. The fix
/// keys on the resulting empty sequence, so every legitimate inserted-sequence
/// form must be unaffected.
#[test]
fn insertions_with_a_payload_still_parse() {
    for input in [
        "NC_000001.11:g.100_101insA",
        "NC_000001.11:g.100_101insATG",
        "NC_000001.11:g.100_101ins10",
        "NC_000001.11:g.100_101ins[A[10];T]",
        "NM_000088.3:c.100_101insA",
        "NC_000001.11:g.849_850ins858_895",
    ] {
        assert!(
            parse_hgvs(input).is_ok(),
            "a payloaded insertion must still parse: {input:?}"
        );
    }
}

/// Regression guard: a `delins` with no inserted sequence is a distinct case
/// (it collapses to a plain `del`, handled elsewhere as W3012) and must not be
/// swept up by the bare-`ins` rejection. A plain `del` is likewise unaffected.
#[test]
fn del_and_bare_delins_are_unaffected() {
    // A bare `delins` still parses — the bare-`ins` rejection lives in
    // `parse_insertion`, which `delins` never reaches. The standard path parses
    // it as-is; lenient rewrites it to the equivalent plain `del` (W3012).
    let standard = parse_hgvs("NC_000001.11:g.100_101delins")
        .expect("a bare `delins` must still parse in standard mode");
    assert_eq!(standard.to_string(), "NC_000001.11:g.100_101delins");

    let lenient = parse_hgvs_lenient("NC_000001.11:g.100_101delins")
        .expect("a bare `delins` must still parse in lenient mode");
    assert_eq!(lenient.result.to_string(), "NC_000001.11:g.100_101del");

    // A plain deletion is valid and stays valid.
    assert!(parse_hgvs("NC_000001.11:g.100_101del").is_ok());
}
