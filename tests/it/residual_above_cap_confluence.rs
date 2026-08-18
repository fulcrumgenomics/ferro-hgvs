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
//! cost confluence, because the un-partitioned whole block was then handled
//! differently from the same variant spelled as its individual changes:
//!
//! ```text
//! a 1100 nt near-palindrome at 257_1356, differing from its own reverse
//! complement at exactly two mirror pairs (257/1356 and 267/1346)
//!
//! g.257_1356inv                                 stable
//! g.[257C>A;267A>C;1346G>T;1356T>G]             denotes the same bases
//! ```
//!
//! One variant, and both spellings must reach one normal form. The cap boundary
//! was exact — 1024 confluent, 1026 not — which is the signature of a length
//! short-circuit rather than anything about the sequence.
//!
//! ## The canonical form is the inversion, not the four substitutions
//!
//! `g.257_1356inv` replaces its whole span with its exact reverse complement, so
//! `rulings[whole-span-reverse-complement-types-as-inv]` (`DNA/inversion.md:5`,
//! 2026-08-13, closing #1703) types it as one `inv` — and it does so **uniformly,
//! with no length or density bound**: a 1100 nt span whose reverse complement
//! coincides with the reference at all but four columns is one `inv` exactly as a
//! 4 nt span is. The ruling considered and **refused** a density floor by name
//! ("supplies no length or density bound … The obvious bound was already refuted
//! by measurement … A later bound would be a new ruling"), and its own withdrawn-
//! costings paragraph names this exact `long_block_inversion`/`near_palindromic_core`
//! shape (four one-column runs at length 1100) as a case it reaches. This
//! overturns the earlier #1230 reading ("inv over-recognized (unchanged interior)
//! -> separate subs") that a prior version of this test asserted.
//!
//! So what this file still pins above the split cap is the cap-lifting itself and
//! confluence: an equal-length block wider than `MAX_SPLIT_BLOCK` is handled (not
//! refused or broken), and both spellings of the one variant converge on the same
//! `inv`.
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
fn the_whole_span_reverse_complement_types_as_one_inv() {
    let core = near_palindrome();
    let output = assert_padded_preserving(&core, "NC_TEST.1:g.257_1356inv");
    // Per `rulings[whole-span-reverse-complement-types-as-inv]` a span replaced by
    // its exact reverse complement is one `inv`, uniformly and with no length or
    // density bound — so this 1100 nt near-palindrome, altering only four bases,
    // stays a single inversion rather than reducing to those four substitutions.
    // Pinned as the exact string: the input already IS the canonical form, so the
    // normalizer must return it unchanged.
    assert_eq!(
        output, "NC_TEST.1:g.257_1356inv",
        "a whole-span reverse complement is one `inv`, not its interior substitutions"
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
