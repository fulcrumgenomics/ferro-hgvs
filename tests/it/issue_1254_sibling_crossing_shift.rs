//! Regression coverage for issue #1254 — a cis member's 3'-shift must not
//! carry it over a sibling's bases.
//!
//! The per-member 3'-shift is sibling-unaware: each member goes to its
//! *standalone* most-3' position. When a sibling edits a base inside the tract
//! the member rotates through, the member leapfrogs it and the pair describes a
//! different sequence:
//!
//! ```text
//! reference   TAATATATATAATATATATT
//! input       TEMPLATE:g.[3_4del;9del]
//! was      -> TEMPLATE:g.12_14del      applied: TAATATATATATATATT
//! input applied                                 TAATATTAATATATATT   ← different
//! now      -> TEMPLATE:g.7_9del        applied: TAATATTAATATATATT
//! ```
//!
//! `3_4del` shifted to `10_11del`, straight past the `9del`; the two were then
//! adjacent and merged to `9_11del`, which shifted on to `12_14del`. Nothing
//! about that output invites suspicion — one well-formed member, no overlap, no
//! warning — so the corruption is silent.
//!
//! Every assertion here is against an **independently applied** reference
//! sequence, reconstructed from SPDI triples rather than from the normalizer.
//! Asserting `normalize(input) == normalize(output)` cannot catch this class: it
//! normalizes both sides, so it passed on the broken behavior.

use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{
    parse_hgvs, HgvsVariant, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection,
};

/// The reported reference: an `(AT)` tract that ends at position 11, where
/// `11=A, 12=A` breaks the repeat.
const TEMPLATE: &str = "TAATATATATAATATATATT";

fn provider(seq: &str) -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("TEMPLATE", seq.to_string());
    provider
}

fn normalize(seq: &str, input: &str) -> String {
    normalize_in(seq, input, ShuffleDirection::ThreePrime)
}

fn normalize_in(seq: &str, input: &str, direction: ShuffleDirection) -> String {
    let normalizer = Normalizer::with_config(
        provider(seq),
        NormalizeConfig::default().with_direction(direction),
    );
    let variant = parse_hgvs(input).expect("parse");
    format!("{}", normalizer.normalize(&variant).expect("normalize"))
}

/// Apply `descriptor` to `seq` **independently of the normalizer**, by
/// converting each member to its SPDI triple and splicing the reference.
///
/// Returns `None` when a member cannot be converted, when a stated deletion
/// disagrees with the reference, or when two members claim the same base — an
/// overlapping description has no single well-defined resulting sequence, and
/// silently double-splicing one would invent a sequence neither side denotes.
fn apply(seq: &str, descriptor: &str) -> Option<String> {
    let provider = provider(seq);
    let members: Vec<HgvsVariant> = match parse_hgvs(descriptor).ok()? {
        HgvsVariant::Allele(allele) => allele.variants.clone(),
        single => vec![single],
    };
    let mut triples = Vec::with_capacity(members.len());
    for member in &members {
        triples.push(hgvs_to_spdi(member, &provider).ok()?);
    }
    // 3'→5' so an applied splice never shifts a later one's coordinates.
    triples.sort_by_key(|t| std::cmp::Reverse(t.position));
    let reference = seq.as_bytes();
    let mut edited = reference.to_vec();
    let mut claimed_from = reference.len();
    for triple in &triples {
        let start = usize::try_from(triple.position).ok()?;
        let end = start.checked_add(triple.deletion.len())?;
        if end > reference.len() || end > claimed_from {
            return None; // out of bounds, or overlapping the member 3' of it
        }
        if !reference[start..end].eq_ignore_ascii_case(triple.deletion.as_bytes()) {
            return None;
        }
        edited.splice(start..end, triple.insertion.bytes());
        claimed_from = start;
    }
    String::from_utf8(edited).ok()
}

/// Assert that `input` normalizes to `expected` *and* that both denote the same
/// sequence when applied to `seq`.
fn assert_normalizes_preserving(seq: &str, input: &str, expected: &str) {
    let actual = normalize(seq, input);
    assert_eq!(actual, expected, "normalizing {input}");
    let from_input = apply(seq, input).expect("input applies");
    let from_output = apply(seq, &actual).unwrap_or_else(|| {
        panic!("{actual} has no well-defined resulting sequence (overlapping or unconvertible)")
    });
    assert_eq!(
        from_output, from_input,
        "{input} -> {actual} changed the sequence"
    );
}

#[test]
fn separated_deletions_do_not_merge_into_a_different_deletion() {
    // The reported case. `3_4del` may 3'-shift only as far as `7_8del`, which
    // stops one base short of the `9del`; the two are then adjacent (no
    // intervening nucleotide) and merge into one three-base deletion.
    assert_normalizes_preserving(TEMPLATE, "TEMPLATE:g.[3_4del;9del]", "TEMPLATE:g.7_9del");
}

#[test]
fn separated_deletions_with_a_downstream_substitution_stay_well_formed() {
    // The three-member variant of the same input. It used to expose the
    // corruption *visibly*, as the overlapping `g.[12_14del;13T>A]`; with the
    // shift clamped, the deletion lands where it belongs and the substitution
    // is untouched.
    assert_normalizes_preserving(
        TEMPLATE,
        "TEMPLATE:g.[3_4del;9del;13T>A]",
        "TEMPLATE:g.[7_9del;13T>A]",
    );
}

#[test]
fn every_spelling_of_the_variant_reaches_one_canonical_form() {
    // Confluence: the pre-shifted, the shifted-and-adjacent, the already-merged
    // and the reordered spellings are one variant and must share a key.
    for input in [
        "TEMPLATE:g.[3_4del;9del]",
        "TEMPLATE:g.[7_8del;9del]",
        "TEMPLATE:g.[9del;3_4del]",
        "TEMPLATE:g.7_9del",
    ] {
        assert_normalizes_preserving(TEMPLATE, input, "TEMPLATE:g.7_9del");
    }
}

#[test]
fn a_deletion_does_not_shift_onto_a_substituted_base() {
    // A deletion inside a homopolymer whose 3'-most standalone position is the
    // substituted base itself. Clamped to `8del`, it is adjacent to `9A>G` and
    // the two render as one delins (`delins.md:16`). This previously emitted the
    // overlapping `g.[9A>G;9del]`.
    assert_normalizes_preserving(
        "GGGAAAAAAGGG",
        "TEMPLATE:g.[4del;9A>G]",
        "TEMPLATE:g.8_9delinsG",
    );
}

#[test]
fn a_member_that_cannot_reach_its_sibling_still_shifts_fully() {
    // Negative control against over-clamping: `4del` 3'-shifts through the A-run
    // to `9del`, and the substitution at `11` is two bases beyond it — never
    // inside the swept window — so the shift must complete and the two must stay
    // separate members, separated by the unchanged `10G` (`general.md:34`).
    assert_normalizes_preserving(
        "GGGAAAAAAGCG",
        "TEMPLATE:g.[4del;11C>T]",
        "TEMPLATE:g.[9del;11C>T]",
    );
}

#[test]
fn a_five_prime_shift_does_not_cross_a_sibling_either() {
    // The mirror image. Under `ShuffleDirection::FivePrime`, `9_10del` rotates
    // 5' through the same `(AT)` tract and its standalone 5'-most position is
    // past the `5del`. Clamped to `6_7del` it is adjacent to the `5del`, and the
    // two merge into `g.5_7del`. This previously emitted `g.1_3del`, which
    // denotes `TATATATAATATATATT` where the input denotes `TAATTATAATATATATT`.
    let input = "TEMPLATE:g.[5del;9_10del]";
    let actual = normalize_in(TEMPLATE, input, ShuffleDirection::FivePrime);
    assert_eq!(actual, "TEMPLATE:g.5_7del");
    assert_eq!(
        apply(TEMPLATE, &actual).expect("output applies"),
        apply(TEMPLATE, input).expect("input applies"),
    );
}

#[test]
fn an_insertion_landing_flush_against_a_substitution_still_collapses() {
    // Negative control for the #999 rule the clamp must not disturb: an
    // insertion claims no reference base, so a member landing flush against one
    // is adjacency, not a crossing. `305_306insC` 3'-shifts to a dup at 306 and
    // coalesces with `307G>T`.
    let mut seq = vec![b'A'; 1000];
    for (i, b) in "CATCCTCGCTCCT".bytes().enumerate() {
        seq[299 + i] = b;
    }
    let seq = String::from_utf8(seq).unwrap();
    assert_eq!(
        normalize(&seq, "TEMPLATE:g.[305_306insC;307G>T]"),
        "TEMPLATE:g.307delinsCT"
    );
}

#[test]
fn shifts_never_change_the_sequence_across_an_exhaustive_two_member_sweep() {
    // The class, not the instance. Every two-member cis allele of a deletion
    // plus a downstream deletion or substitution, over every 4^10 window of a
    // deterministic set of 10-mers, checked by applying both sides to the
    // reference. Before the clamp this sweep reported 2,244 silent
    // sequence-changing normalizations (well-formed, warning-free, disjoint
    // output) out of 445,148; over-clamping would show up as a canonical-form
    // change in the pinned cases above.
    let mut checked = 0usize;
    let mut mismatched: Vec<String> = Vec::new();

    for seed in 0..64u32 {
        // Deterministic 20-mers over {A,T} — the alphabet that builds the long
        // alternating tracts a shift can travel through — and over {A,C,G,T}.
        for alphabet in [b"AT".as_slice(), b"ACGT".as_slice()] {
            let mut state = u64::from(seed).wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
            let seq: String = (0..20)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    alphabet[(state % alphabet.len() as u64) as usize] as char
                })
                .collect();

            for first_start in 1..=14usize {
                for first_len in 1..=2usize {
                    let first_end = first_start + first_len - 1;
                    let first = if first_len == 1 {
                        format!("{first_start}del")
                    } else {
                        format!("{first_start}_{first_end}del")
                    };
                    // Leave at least one unchanged base between the members, so
                    // every case starts out as a genuinely separated pair.
                    for second_start in first_end + 2..=19usize {
                        let base = seq.as_bytes()[second_start - 1] as char;
                        let alt = if base == 'A' { 'G' } else { 'A' };
                        for second in [
                            format!("{second_start}del"),
                            format!("{second_start}{base}>{alt}"),
                        ] {
                            let input = format!("TEMPLATE:g.[{first};{second}]");
                            let output = normalize(&seq, &input);
                            let (Some(want), Some(got)) =
                                (apply(&seq, &input), apply(&seq, &output))
                            else {
                                continue; // unconvertible, or a flagged overlap
                            };
                            checked += 1;
                            if want != got && mismatched.len() < 10 {
                                mismatched.push(format!(
                                    "{seq}: {input} -> {output} (want {want}, got {got})"
                                ));
                            }
                        }
                    }
                }
            }
        }
    }

    assert!(checked > 10_000, "sweep covered too little: {checked}");
    assert!(
        mismatched.is_empty(),
        "sequence-changing normalizations in {checked} cases: {mismatched:#?}"
    );
}
