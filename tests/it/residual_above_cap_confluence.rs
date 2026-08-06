//! A block larger than `MAX_SPLIT_BLOCK` is still partitionable when its two
//! sides are the same length.
//!
//! `partition_block` short-circuits to the whole block on length alone:
//!
//! ```text
//! if reference.len() > MAX_SPLIT_BLOCK || result.len() > MAX_SPLIT_BLOCK { whole() }
//! ```
//!
//! The cap is there to bound `best_alignment`'s gap-placement sweep, which is
//! O(n²) in the block length — at 1024 that is already ~1e6 column comparisons.
//! But that sweep only runs when the two sides differ in length. For an
//! **equal-length** block `best_alignment` returns immediately:
//!
//! ```text
//! if reference.len() == result.len() {
//!     return Some((0..reference.len()).map(|i| (Some(i), Some(i))).collect());
//! }
//! ```
//!
//! There is no gap to place, so there is nothing to search and nothing for the
//! cap to protect — the partition is position-wise and linear. Capping it anyway
//! cost confluence, because the un-partitioned whole block is then refused by the
//! weight bound whenever the input was spelled as its individual changes:
//!
//! ```text
//! a 1100 nt near-palindrome at 257_1356, differing from its own reverse
//! complement at exactly two mirror pairs (257/1356 and 267/1346)
//!
//! g.257_1356inv                                 stable
//! g.[257C>A;267A>C;1346G>T;1356T>G]             stable, and denotes the same bases
//! ```
//!
//! One variant, two stable normal forms. The cap boundary was exact — 1024
//! confluent, 1026 not — which is the signature of a length short-circuit rather
//! than anything about the sequence.
//!
//! The four substitutions are the canonical form here, and not merely the
//! lighter one: an inversion whose interior is unchanged is the #1230 shape from
//! #1235's own table ("inv over-recognized (unchanged interior) -> separate
//! subs"). A 1100 nt inversion that alters four bases is not an inversion.
//!
//! ## What still bounds the cost
//!
//! Lifting the cap here does not open an unbounded path, because
//! `MAX_CANONICAL_WINDOW` (4096) already refuses a wider group before any block
//! is partitioned — it bounds the reference window the canonicalizer will fetch
//! at all. Measured on whole-span inversions: 2,000 nt and 3,000 nt partition
//! (1.6 ms and 0.3 ms), 4,000 nt and beyond are refused at the window and stay
//! `inv`. So this change moves the equal-length confluence ceiling from 1024 up
//! to that window bound, and leaves the guard that actually caps cost in place.

use crate::common::synthetic::assert_padded_preserving;

/// Length of the near-palindromic span, chosen to sit above `MAX_SPLIT_BLOCK`
/// (1024) so the short-circuit is what the test exercises.
const SPAN: usize = 1100;

/// Offsets into the span whose bases are perturbed away from palindromy. Each
/// one breaks its own column *and* its mirror, so two perturbations yield four
/// changed positions.
const PERTURBED: [usize; 2] = [0, 10];

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

/// A span that is its own reverse complement except at [`PERTURBED`] and their
/// mirrors, built deterministically rather than committed as a fixture.
fn near_palindrome() -> String {
    const UNIT: &[u8] = b"ACGTTGCAAGCT";
    let mut bases = vec![0u8; SPAN];
    // Fill the 5' half, then mirror it as the reverse complement so the span
    // inverts to itself.
    for i in 0..SPAN / 2 {
        let base = UNIT[i % UNIT.len()];
        bases[i] = base;
        bases[SPAN - 1 - i] = complement(base);
    }
    // Break palindromy at two positions. `inv` then differs from the reference
    // at those two columns and at their two mirrors: four changed bases.
    for &offset in &PERTURBED {
        bases[offset] = match bases[offset] {
            b'A' => b'C',
            b'C' => b'A',
            b'G' => b'T',
            _ => b'G',
        };
    }
    String::from_utf8(bases).expect("ASCII bases")
}

/// The four positions at which the span differs from its own reverse
/// complement, as 1-based coordinates in the padded reference (the core starts
/// at 257).
fn changed_positions() -> Vec<usize> {
    let mut positions: Vec<usize> = PERTURBED
        .iter()
        .flat_map(|&offset| [257 + offset, 257 + SPAN - 1 - offset])
        .collect();
    positions.sort_unstable();
    positions
}

#[test]
fn the_inversion_spelling_reduces_to_its_four_substitutions() {
    let core = near_palindrome();
    let output = assert_padded_preserving(&core, "NC_TEST.1:g.257_1356inv");
    assert!(
        !output.contains("inv"),
        "a 1100 nt inversion that alters four bases must not stay an inversion: {output}"
    );
    // Not pinned as a literal string — the bases come from `near_palindrome`, so
    // pinning one would restate the generator rather than check it. Pinned
    // against the member list that generator implies instead, which checks the
    // substituted bases and not only the positions. Containment on the bare
    // positions was the weaker earlier form: `257` also matches inside `1257`,
    // and nothing there read the payload at all.
    assert_eq!(
        output,
        substitution_spelling(&core),
        "the inversion must reduce to exactly the four substitutions it performs"
    );
}

/// The same variant spelled as its individual substitutions, built from the
/// generated span so the two spellings cannot drift apart.
///
/// `inv` replaces the span with its reverse complement, so the base the
/// inversion puts at offset `i` is `complement(span[SPAN - 1 - i])`.
fn substitution_spelling(core: &str) -> String {
    let bases = core.as_bytes();
    let members: Vec<String> = changed_positions()
        .into_iter()
        .map(|position| {
            let offset = position - 257;
            let from = bases[offset] as char;
            let to = complement(bases[SPAN - 1 - offset]) as char;
            format!("{position}{from}>{to}")
        })
        .collect();
    format!("NC_TEST.1:g.[{}]", members.join(";"))
}

#[test]
fn both_spellings_reach_the_same_normal_form() {
    let core = near_palindrome();
    let from_inversion = assert_padded_preserving(&core, "NC_TEST.1:g.257_1356inv");
    let from_substitutions = assert_padded_preserving(&core, &substitution_spelling(&core));
    assert_eq!(
        from_inversion, from_substitutions,
        "the two spellings of one variant must converge on one string"
    );
}
