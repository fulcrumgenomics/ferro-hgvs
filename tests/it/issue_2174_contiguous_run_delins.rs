//! #2174 — a **contiguous** changed run (no interior reference base retained
//! *in place*) must be described as the single spanning `delins` that
//! `DNA/delins.md:16` recommends, not fragmented into a `del`+`ins` / `del`+`dup`
//! shift-anchoring.
//!
//! The minimal-edit alignment of a balanced run (e.g. `GACT -> ACTG`) is
//! `del`+3·match+`ins` at cost 2, cheaper than the 4-substitution position-wise
//! alignment — so before `#2174` ferro emitted `g.[11del;14_15insG]`. The three
//! interior "matched" bases are matched only under a *shifted* alignment; none
//! is a reference base retained at its own coordinate. With no in-place boundary
//! the run is one member, typed `delins`.
//!
//! All blocks are equal-length and every position differs, so applying either
//! form to the reference yields the same resulting sequence (a pure
//! representation change).

use ferro_hgvs::{from_sequences, FromSequencesOptions};

fn derive(refe: &str, alt: &str) -> String {
    from_sequences("TEMPLATE", 1, refe, alt, &FromSequencesOptions::default())
        .expect("from_sequences must derive this block")
        .to_string()
}

/// `g.11_14 GACT -> ACTG`: was `g.[11del;14_15insG]`; `coalesce_solid_run`
/// collapses the contiguous run into one `delins`.
#[test]
fn a_contiguous_equal_length_run_collapses_to_one_delins() {
    let got = derive("ACGTTCAGGTGACTTTAGCTAGCTAG", "ACGTTCAGGTACTGTTAGCTAGCTAG");
    assert_eq!(got, "TEMPLATE:g.11_14delinsACTG");
}

/// `g.11_13 GAC -> ACT`: was `g.[11del;15dup]` — a single-base `dup` two bases
/// 3' of the change, the shift anchor `#2174` keeps merged. `coalesce_solid_run`
/// collapses the contiguous run into one `delins`.
#[test]
fn a_del_dup_shift_anchor_collapses_to_one_delins() {
    let got = derive("ACGTTCAGGTGACTTAGCTAGCTAG", "ACGTTCAGGTACTTTAGCTAGCTAG");
    assert_eq!(got, "TEMPLATE:g.11_13delinsACT");
}
