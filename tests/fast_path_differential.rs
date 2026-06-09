//! Differential identity test: `parse_hgvs_fast` must agree with `parse_hgvs`.
//!
//! `parse_hgvs` routes through the full nom parser (`parse_variant`).
//! `parse_hgvs_fast` first tries a specialized substitution parser
//! (`fast_path::try_fast_path`) and falls back to `parse_variant` for
//! everything it does not recognize.
//!
//! Making the fast path the *default* (`C1`) is only sound if the two entry
//! points are observationally identical for **every** input — i.e. for any
//! string they either both reject, or both accept and yield the same
//! [`HgvsVariant`]. A single divergence (one accepts where the other rejects,
//! or they produce different ASTs) would be a behavior change, not a perf win.
//!
//! This test pins that invariant two ways:
//!   1. a deterministic table covering each documented fast-path pattern and
//!      its fallback counterparts (substitutions on every accession flavor,
//!      plus intronic / UTR / non-coding / RNA / del / ins / dup / delins /
//!      inv that must route to the generic parser), and
//!   2. a `proptest` fuzzer that assembles structurally-varied HGVS strings
//!      and asserts agreement on each.
//!
//! It is intentionally a pure correctness gate: it does not flip the default,
//! it proves whether the default *could* be flipped safely.

use ferro_hgvs::{parse_hgvs, parse_hgvs_fast};
use proptest::prelude::*;

/// Assert that the standard and fast entry points agree on `input`.
///
/// Agreement means: both `Err`, or both `Ok` with equal variants. Returns a
/// human-readable description of the divergence on mismatch (so both the
/// table test and the proptest fuzzer can surface the exact offending input)
/// or `None` when they agree.
fn divergence(input: &str) -> Option<String> {
    let standard = parse_hgvs(input);
    let fast = parse_hgvs_fast(input);
    match (standard, fast) {
        (Ok(a), Ok(b)) if a == b => None,
        (Err(_), Err(_)) => None,
        (Ok(a), Ok(b)) => Some(format!(
            "both parsed but differ:\n  input:    {input:?}\n  standard: {a:?}\n  fast:     {b:?}"
        )),
        (Ok(a), Err(e)) => Some(format!(
            "standard accepted, fast rejected:\n  input:    {input:?}\n  standard: {a:?}\n  fast err: {e:?}"
        )),
        (Err(e), Ok(b)) => Some(format!(
            "fast accepted, standard rejected:\n  input:    {input:?}\n  fast:     {b:?}\n  std err:  {e:?}"
        )),
    }
}

/// Representative inputs spanning every documented fast-path pattern and the
/// fallback patterns it must defer on. If the fast path ever diverges from the
/// generic parser on one of these, flipping the default is unsafe.
const DIFFERENTIAL_CASES: &[&str] = &[
    // --- fast-path substitutions (every accession flavor) ---
    "NC_000001.11:g.12345A>G",
    "NC_000023.11:g.1A>T",
    "NM_000088.3:c.459A>G",
    "NM_000088.3:c.1G>C",
    "ENST00000357033.8:c.100A>G",
    "ENSG00000139618.15:g.500C>T",
    "LRG_1:g.5000A>G",
    "LRG_1t1:c.100A>G",
    "GRCh38(chr1):g.12345A>G",
    "GRCh37(chrX):g.999G>A",
    // --- non-coding / RNA substitutions ---
    "NR_003051.3:n.100A>G",
    "NM_000088.3:r.100a>g",
    // --- intronic / UTR (must fall back) ---
    "NM_000088.3:c.100+5A>G",
    "NM_000088.3:c.100-5A>G",
    "NM_000088.3:c.-50A>G",
    "NM_000088.3:c.*100A>G",
    "NM_000088.3:c.*1G>A",
    // --- non-substitution edits (must fall back) ---
    "NM_000088.3:c.459del",
    "NM_000088.3:c.459_460del",
    "NM_000088.3:c.459dup",
    "NM_000088.3:c.459_460insACGT",
    "NM_000088.3:c.459_460delinsAC",
    "NM_000088.3:c.459_460inv",
    "NC_000001.11:g.12345_12350del",
    // --- whitespace / trimming edges ---
    "  NC_000001.11:g.12345A>G  ",
    "\tNM_000088.3:c.459A>G\n",
    // --- malformed (both must reject) ---
    "",
    "not an hgvs string",
    "NM_000088.3:c.A>G",
    "NC_000001.11:g.12345A>",
    "NM_000088.3:c.459X>G",
    "NM_000088.3:p.Arg100Gly",
];

#[test]
fn fast_path_matches_standard_on_table() {
    let mut failures = Vec::new();
    for &input in DIFFERENTIAL_CASES {
        if let Some(msg) = divergence(input) {
            failures.push(msg);
        }
    }
    assert!(
        failures.is_empty(),
        "fast path diverged from standard parser on {} case(s):\n\n{}",
        failures.len(),
        failures.join("\n\n")
    );
}

/// A nucleotide base for substitution generation.
fn base() -> impl Strategy<Value = &'static str> {
    prop_oneof![Just("A"), Just("C"), Just("G"), Just("T")]
}

/// An accession prefix spanning the flavors the fast path special-cases plus
/// ones it does not, so both `Success` and `Fallback` branches are exercised.
fn accession() -> impl Strategy<Value = &'static str> {
    prop_oneof![
        Just("NC_000001.11"),
        Just("NM_000088.3"),
        Just("NR_003051.3"),
        Just("ENST00000357033.8"),
        Just("ENSG00000139618.15"),
        Just("LRG_1"),
        Just("LRG_1t1"),
        Just("GRCh38(chr1)"),
        Just("GRCh37(chrX)"),
    ]
}

/// Coordinate-system letter.
fn coord_system() -> impl Strategy<Value = &'static str> {
    prop_oneof![Just("g"), Just("c"), Just("n"), Just("r")]
}

prop_compose! {
    /// A bare substitution string, e.g. `NM_000088.3:c.459A>G` — the fast
    /// path's core domain.
    fn substitution()(
        acc in accession(),
        sys in coord_system(),
        pos in 1u64..200_000,
        from in base(),
        to in base(),
    ) -> String {
        format!("{acc}:{sys}.{pos}{from}>{to}")
    }
}

prop_compose! {
    /// An intronic / UTR-offset substitution, e.g. `c.100+5A>G` or
    /// `c.-50A>G` — structurally a substitution but outside the fast path,
    /// so it must route to the generic parser identically.
    fn offset_substitution()(
        acc in accession(),
        pos in 1u64..50_000,
        offset in prop_oneof![Just(""), Just("+5"), Just("-5"), Just("*"), Just("-")],
        from in base(),
        to in base(),
    ) -> String {
        format!("{acc}:c.{offset}{pos}{from}>{to}")
    }
}

prop_compose! {
    /// A non-substitution edit (del/ins/dup/delins/inv) that the fast path
    /// never claims — both parsers must agree (typically both accept).
    fn non_sub_edit()(
        acc in accession(),
        pos in 1u64..50_000,
        span in 1u64..20,
        edit in prop_oneof![
            Just("del"),
            Just("dup"),
            Just("inv"),
            Just("insACGT"),
            Just("delinsAC"),
        ],
    ) -> String {
        format!("{acc}:c.{pos}_{end}{edit}", end = pos + span)
    }
}

prop_compose! {
    /// A truncated / boundary-length input: `{acc}:{sys}.{tail}` where `tail`
    /// is a short string from the HGVS punctuation alphabet. This deliberately
    /// stresses the fast path's byte-index arithmetic with bodies too short to
    /// hold a position + edit (it is how the `c.A>G` slice-underflow panic was
    /// found). The two entry points must still agree — and neither may panic.
    fn truncated_boundary()(
        acc in accession(),
        sys in coord_system(),
        tail in proptest::collection::vec(
            prop_oneof![
                Just('A'), Just('C'), Just('G'), Just('T'),
                Just('0'), Just('1'), Just('9'),
                Just('>'), Just('+'), Just('-'), Just('*'), Just('_'), Just('.'),
            ],
            0..6,
        ),
    ) -> String {
        let tail: String = tail.into_iter().collect();
        format!("{acc}:{sys}.{tail}")
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(5000))]

    #[test]
    fn fast_path_matches_standard_substitution(input in substitution()) {
        prop_assert!(divergence(&input).is_none(), "{}", divergence(&input).unwrap());
    }

    #[test]
    fn fast_path_matches_standard_offset(input in offset_substitution()) {
        prop_assert!(divergence(&input).is_none(), "{}", divergence(&input).unwrap());
    }

    #[test]
    fn fast_path_matches_standard_non_sub(input in non_sub_edit()) {
        prop_assert!(divergence(&input).is_none(), "{}", divergence(&input).unwrap());
    }

    #[test]
    fn fast_path_matches_standard_truncated(input in truncated_boundary()) {
        prop_assert!(divergence(&input).is_none(), "{}", divergence(&input).unwrap());
    }
}
