//! A normalized cis allele renders its members in ascending coordinate order.
//!
//! `sort_cis_members_by_genomic_order` (#1098) exists to make the rendered order
//! canonical and input-order-independent, and it keys on
//! `cis_member_order_key`. That key's only positional discriminator was the
//! member's **start**; when two members shared a start it fell through to a
//! lexical compare of the rendered descriptor, which is arbitrary with respect
//! to where the members actually sit:
//!
//! ```text
//! reference   ("ACGT" x 64) + "CAGTATGCAGGCAA" + ("ACGT" x 64)
//! in   g.[258_259insA;259_260insAG]
//! out  g.[258_259dup;258dup]          <- descending
//! ```
//!
//! Both members start at 258, so every positional field tied and
//! `"258_259dup" < "258dup"` decided it (`_` is 0x5F, `d` is 0x64). But a
//! duplication places its copy at the junction *after* its last base, so those
//! two sit at interbase 259 and 258 respectively — ordered by their **end**,
//! the one field the key omitted.
//!
//! This is the ordering-key twin of #1276/#1277, which found the same blind spot
//! in the overlap detector: a duplication occupies a junction, not just a span.
//!
//! The fix adds the member's end to the key between the offset and the
//! descriptor. Members whose spans do not overlap already agree on start-order,
//! end-order and junction-order, so this can only reorder members that share a
//! start — which is precisely the case the descriptor tie-break was deciding
//! arbitrarily. The descriptor stays as the final tie-break, so the key remains
//! total and the result input-order-independent.
//!
//! The sequence is correct on both sides here; this is an ordering violation
//! only, so every assertion below pairs the order check with the apply oracle.

use crate::common::cis_apply_oracle::{apply, assert_members_ascending, normalize};

/// `("ACGT" x 64) + core + ("ACGT" x 64)`, so `core[i]` is at 1-based position
/// `257 + i`. Padding keeps every shift well away from a contig boundary.
fn padded(core: &str) -> String {
    format!("{}{}{}", "ACGT".repeat(64), core, "ACGT".repeat(64))
}

/// Assert that `input` normalizes to members in ascending interbase order and
/// that normalization preserved the denoted sequence.
fn assert_ordered_and_preserving(core: &str, input: &str) -> String {
    let seq = padded(core);
    let output = normalize(&seq, input);
    assert_members_ascending(&seq, &output);
    // Both sides are unwrapped rather than compared as `Option`s: `apply`
    // declines an overlapping description, and `assert_eq!(None, None)` would
    // report preservation it never checked.
    let from_input = apply(&seq, input).expect("input applies");
    let from_output = apply(&seq, &output).unwrap_or_else(|| {
        panic!("`{input}` normalized to `{output}`, which has no well-defined resulting sequence")
    });
    assert_eq!(
        from_output, from_input,
        "`{input}` -> `{output}` changed the sequence"
    );
    output
}

#[test]
fn two_duplications_sharing_a_start_render_in_junction_order() {
    // The #1261 reproduction. `258dup` sits at interbase 258 and `258_259dup`
    // at 259, so the narrower one must come first.
    let output =
        assert_ordered_and_preserving("CAGTATGCAGGCAA", "TEMPLATE:g.[258_259insA;259_260insAG]");
    assert_eq!(output, "TEMPLATE:g.[258dup;258_259dup]");
}

#[test]
fn the_rendered_order_does_not_depend_on_the_authored_order() {
    // The whole point of #1098's sort. Both spellings must reach one string.
    let forward =
        assert_ordered_and_preserving("CAGTATGCAGGCAA", "TEMPLATE:g.[258_259insA;259_260insAG]");
    let reverse =
        assert_ordered_and_preserving("CAGTATGCAGGCAA", "TEMPLATE:g.[259_260insAG;258_259insA]");
    assert_eq!(forward, reverse);
}

#[test]
fn members_with_distinct_starts_are_unaffected() {
    // Non-overlapping spans already agree on start- and end-order, so adding
    // the end to the key must not move them. Authored both ways round.
    for input in ["TEMPLATE:g.[258A>C;266G>T]", "TEMPLATE:g.[266G>T;258A>C]"] {
        let output = assert_ordered_and_preserving("CAGTATGCAGGCAA", input);
        assert_eq!(output, "TEMPLATE:g.[258A>C;266G>T]");
    }
}

#[test]
fn a_duplication_beside_a_later_substitution_keeps_its_place() {
    // A dup's junction is its end, which for a disjoint sibling is still 5' of
    // that sibling's start — so the dup stays first.
    let output = assert_ordered_and_preserving("CAGTATGCAGGCAA", "TEMPLATE:g.[266G>T;258_259dup]");
    assert!(
        output.starts_with("TEMPLATE:g.[258"),
        "the duplication must render first, got `{output}`"
    );
}

#[test]
fn protein_members_sharing_a_start_render_in_end_order() {
    // The protein axis has the same shared-start shape: `p.Ala3dup` and
    // `p.Ala3_Ala4dup` both start at residue 3. Ordering by end puts the
    // narrower first, matching the nucleotide axes. No reference needed —
    // ordering is positional, so a bare parse/normalize round trip shows it.
    use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};
    let normalizer = Normalizer::new(MockProvider::new());
    for input in [
        "NP_000001.1:p.[Ala3_Ala4dup;Ala3dup]",
        "NP_000001.1:p.[Ala3dup;Ala3_Ala4dup]",
    ] {
        let variant = parse_hgvs(input).expect("parse");
        let output = normalizer
            .normalize(&variant)
            .expect("normalize")
            .to_string();
        assert_eq!(
            output, "NP_000001.1:p.[Ala3dup;Ala3_Ala4dup]",
            "authored as `{input}`"
        );
    }
}
