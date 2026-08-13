//! One definition of the unknown-offset sentinel pair, guarded against
//! re-divergence.
//!
//! `+?` / `-?` are stored in band, as two extreme `i64` values. Those two values
//! used to be declared **twice**: once by the parser that produces them
//! (`parser::position::OFFSET_UNKNOWN_{POSITIVE,NEGATIVE}`) and once, as bare
//! `i64::MAX` / `i64::MIN` literals, by `location::GENOME_OFFSET_UNKNOWN_*`.
//! Nothing tied the two declarations together, so changing either would leave
//! the other behind: the parser would store one value and every `Display` — all
//! of which key off the parser's pair — would fail to recognise it and print a
//! raw 19-digit integer where a `?` belongs. That is the leak #1762 fixed for
//! `TxPos`/`RnaPos`, re-armed by a different mechanism.
//!
//! The collapse makes the divergence unrepresentable: `GENOME_OFFSET_UNKNOWN_*`
//! is now a re-export of the parser's pair, so there is one definition and one
//! literal `i64::MAX` / `i64::MIN` in the crate. The names are kept because they
//! are public API of a published crate.
//!
//! **Why these assertions are not tautological.** Comparing the two constants to
//! each other would be — they are the same item today, so `assert_eq!` holds by
//! definition and would keep holding if the collapse were undone with two
//! matching literals. So nothing here compares a constant to a constant. Every
//! assertion compares a **named constant against a value the parser actually
//! produced** from the text `+?` / `-?`, or against what the renderer actually
//! prints. The parser is driven by the canonical constant and knows nothing of
//! the alias; so if the alias is ever re-declared independently and drifts, the
//! parsed value stops matching the name this test asserts it against, and the
//! test fails. That is true of a re-declaration that is *initially* equal, too —
//! it only has to be checked once it drifts, which is exactly when it matters.

use ferro_hgvs::hgvs::location::{
    CdsPos, GenomePos, RnaPos, TxPos, GENOME_OFFSET_UNKNOWN_NEGATIVE,
    GENOME_OFFSET_UNKNOWN_POSITIVE,
};
use ferro_hgvs::hgvs::parser::position::{
    parse_cds_pos, parse_genome_pos, parse_rna_pos, parse_tx_pos, OFFSET_UNKNOWN_NEGATIVE,
    OFFSET_UNKNOWN_POSITIVE,
};
use ferro_hgvs::parse_hgvs;

/// The genome axis is the one that carries the alias, so it is the axis where a
/// re-divergence would first show. `parse_genome_pos` resolves `+?` through the
/// shared `parse_offset`, which is keyed off the *canonical* constant — so this
/// asserts the alias against a value produced without reference to it.
#[test]
fn the_genome_axis_alias_names_the_value_the_parser_produces() {
    let (rest, pos) = parse_genome_pos("100+?").expect("g. position with an unknown offset");
    assert!(rest.is_empty(), "the whole position must be consumed");
    assert_eq!(
        pos.offset,
        Some(GENOME_OFFSET_UNKNOWN_POSITIVE),
        "the parser stored an offset that GENOME_OFFSET_UNKNOWN_POSITIVE does \
         not name — the alias has drifted from the definition it re-exports"
    );

    let (rest, pos) = parse_genome_pos("100-?").expect("g. position with an unknown offset");
    assert!(rest.is_empty(), "the whole position must be consumed");
    assert_eq!(
        pos.offset,
        Some(GENOME_OFFSET_UNKNOWN_NEGATIVE),
        "the parser stored an offset that GENOME_OFFSET_UNKNOWN_NEGATIVE does \
         not name — the alias has drifted from the definition it re-exports"
    );
}

/// The same question on the three transcript-relative axes, under the canonical
/// name. Every axis stores the one pair, which is the fact that makes a
/// genome-specific spelling of it wrong.
#[test]
fn every_axis_stores_the_one_sentinel_pair_the_parser_produces() {
    let (rest, cds) = parse_cds_pos("100+?").expect("c. position with an unknown offset");
    assert!(rest.is_empty(), "c.100+? must be consumed whole");
    assert_eq!(cds.offset, Some(OFFSET_UNKNOWN_POSITIVE), "c.100+?");
    let (rest, cds) = parse_cds_pos("100-?").expect("c. position with an unknown offset");
    assert!(rest.is_empty(), "c.100-? must be consumed whole");
    assert_eq!(cds.offset, Some(OFFSET_UNKNOWN_NEGATIVE), "c.100-?");

    let (rest, tx) = parse_tx_pos("100+?").expect("n. position with an unknown offset");
    assert!(rest.is_empty(), "n.100+? must be consumed whole");
    assert_eq!(tx.offset, Some(OFFSET_UNKNOWN_POSITIVE), "n.100+?");
    let (rest, tx) = parse_tx_pos("100-?").expect("n. position with an unknown offset");
    assert!(rest.is_empty(), "n.100-? must be consumed whole");
    assert_eq!(tx.offset, Some(OFFSET_UNKNOWN_NEGATIVE), "n.100-?");

    let (rest, rna) = parse_rna_pos("100+?").expect("r. position with an unknown offset");
    assert!(rest.is_empty(), "r.100+? must be consumed whole");
    assert_eq!(rna.offset, Some(OFFSET_UNKNOWN_POSITIVE), "r.100+?");
    let (rest, rna) = parse_rna_pos("100-?").expect("r. position with an unknown offset");
    assert!(rest.is_empty(), "r.100-? must be consumed whole");
    assert_eq!(rna.offset, Some(OFFSET_UNKNOWN_NEGATIVE), "r.100-?");
}

/// The other half of the seam: a position *built from the named constant* must
/// render as `?`. Together with the two tests above this closes the loop —
/// parser -> constant -> renderer — without ever comparing a constant to a
/// constant. A drifted alias fails here by printing its 19-digit value, which is
/// the exact symptom the sentinel pair exists to prevent.
#[test]
fn a_position_built_from_the_named_sentinel_renders_as_a_question_mark() {
    assert_eq!(
        GenomePos::with_offset(100, GENOME_OFFSET_UNKNOWN_POSITIVE).to_string(),
        "100+?"
    );
    assert_eq!(
        GenomePos::with_offset(100, GENOME_OFFSET_UNKNOWN_NEGATIVE).to_string(),
        "100-?"
    );
    assert_eq!(
        CdsPos::with_offset(100, OFFSET_UNKNOWN_POSITIVE).to_string(),
        "100+?"
    );
    assert_eq!(
        TxPos::with_offset(100, OFFSET_UNKNOWN_NEGATIVE).to_string(),
        "100-?"
    );
    assert_eq!(
        RnaPos::with_offset(100, OFFSET_UNKNOWN_POSITIVE).to_string(),
        "100+?"
    );
}

/// End to end through the public entry point, on the genome axis — the axis
/// #1762's round-trip test did not cover, and the one whose sentinel name is now
/// an alias. A drifted alias would not fail this on its own; it is here so the
/// alias's axis has a whole-description guard rather than only a position-level
/// one.
#[test]
fn a_parsed_genomic_unknown_offset_round_trips() {
    for input in ["NC_000001.11:g.100+?del", "NC_000001.11:g.100-?del"] {
        let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
        assert_eq!(
            parsed.to_string(),
            input,
            "an unknown offset must render back as `?`, not as the i64 sentinel"
        );
    }
}

/// A control, so the assertions above cannot pass by the renderer printing `?`
/// for everything: an ordinary offset one step inside the sentinel still renders
/// as a number. This is what pins the sentinel band at exactly two values.
#[test]
fn an_ordinary_offset_next_to_the_sentinel_still_renders_numerically() {
    let near_positive = OFFSET_UNKNOWN_POSITIVE - 1;
    let near_negative = OFFSET_UNKNOWN_NEGATIVE + 1;

    assert_eq!(
        GenomePos::with_offset(100, near_positive).to_string(),
        format!("100+{near_positive}"),
        "only the sentinel itself may render as `+?`"
    );
    assert_eq!(
        CdsPos::with_offset(100, near_negative).to_string(),
        format!("100{near_negative}"),
        "only the sentinel itself may render as `-?`"
    );
}
