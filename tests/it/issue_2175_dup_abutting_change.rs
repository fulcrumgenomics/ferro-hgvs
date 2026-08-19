//! #2175 — a tandem duplication that abuts an adjacent change must still be
//! surfaced as a `dup` (`DNA/duplication.md:18` is a MUST), with the adjacent
//! change kept as its own member — not spread into anonymous insertions or
//! flattened into one `delins`.
//!
//! The isolated expansion `CACA -> CACACA` correctly derives `g.26_27dup`. Once
//! the immediately-3' base also changes (`g.28A>C`), the minimal-edit derivation
//! spreads the tandem copy across two single-base insertions
//! (`g.[27_28insC;28_29insC]`) and the duplication vanishes. The dup must be
//! peeled out and the substitution kept separate.

use ferro_hgvs::{from_sequences, FromSequencesOptions, ShuffleDirection};

fn derive(refe: &str, alt: &str) -> String {
    from_sequences("TEMPLATE", 1, refe, alt, &FromSequencesOptions::default())
        .expect("from_sequences must derive this block")
        .to_string()
}

fn derive_5p(refe: &str, alt: &str) -> String {
    from_sequences(
        "TEMPLATE",
        1,
        refe,
        alt,
        &FromSequencesOptions::default().with_direction(ShuffleDirection::FivePrime),
    )
    .expect("from_sequences must derive this block")
    .to_string()
}

// A `CA` tandem array at g.24_27 (`CACA`), then `AAT`.
const REF: &str = "TAGTAAACCATTTTACGGAGGATCACAAATTCCTCCTTAT";

/// The isolated one-unit expansion is a clean `dup` — the control that must not
/// regress.
#[test]
fn an_isolated_tandem_expansion_is_a_dup() {
    let alt = "TAGTAAACCATTTTACGGAGGATCACACAAATTCCTCCTTAT";
    assert_eq!(derive(REF, alt), "TEMPLATE:g.26_27dup");
}

/// The same expansion plus the immediately-3' `g.28A>C`: the duplication must
/// survive, with the substitution as its own member.
#[test]
fn a_tandem_expansion_abutting_a_substitution_keeps_the_dup() {
    let alt = "TAGTAAACCATTTTACGGAGGATCACACACATTCCTCCTTAT";
    assert_eq!(derive(REF, alt), "TEMPLATE:g.[26_27dup;28A>C]");
}

/// The same case under 5' placement. `peel_tandem_dup_beside_change` takes no
/// direction — its descending scan always selects the 3'-most duplication — so
/// the peeled shape is pinned here to prove the 5' shuffle leaves it where the
/// 3' pass does, and to give the `#2175` peel a regression under the non-default
/// direction that the rest of the coverage never exercises.
#[test]
fn a_tandem_expansion_abutting_a_substitution_keeps_the_dup_under_five_prime() {
    let alt = "TAGTAAACCATTTTACGGAGGATCACACACATTCCTCCTTAT";
    assert_eq!(derive_5p(REF, alt), "TEMPLATE:g.[26_27dup;28A>C]");
}
