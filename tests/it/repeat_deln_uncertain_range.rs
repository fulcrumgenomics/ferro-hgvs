//! N-padded deletion over an uncertain range: `(start_end)delN[k]` (#545).
//!
//! A fragment containing a repeat was amplified and its size determined to be
//! `k` nucleotides smaller than the reference, without sequencing — so the
//! deleted size is `N[k]` (k unknown bases) over the uncertain amplified range
//! `(start_end)`. Spec: `recommendations/uncertain.md` (lines 38-40),
//! `recommendations/DNA/repeated.md` (line 51). The matching `insN[k]` form
//! already parses; `delN[k]` was the gap.

use ferro_hgvs::parse_hgvs;

#[test]
fn deln_over_uncertain_range_round_trips() {
    for s in [
        "NM_000333.3:c.(4_246)delN[15]",
        "NC_000003.12:g.(63912602_63912844)delN[15]",
        "NM_000333.3:c.(4_246)delN[(150_180)]",
    ] {
        let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
        assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
    }
}
