//! Shared oracle and corpus for cis-allele normalization tests.
//!
//! Every sibling-crossing test asks the same question — did normalizing this
//! allele change the sequence it denotes? — and answering it requires applying
//! a description to the reference *without* going through the normalizer, or
//! the check is circular.
//!
//! [`apply_with`] is the load-bearing piece: its `claimed_from` walk is what
//! makes "declines an overlapping description" true, and every
//! sequence-preservation assertion in `cis_junction_crossing_shift.rs`,
//! `issue_1234_sibling_clamped_shift.rs`, `issue_1254_sibling_crossing_shift.rs`,
//! `issue_1261_cis_member_order.rs`, `issue_1281_reducing_member_shift.rs` and
//! `repeat_span_sibling_overlap.rs` rests on it. It lives here rather than
//! being copied into each of them so a change cannot silently weaken one file's
//! oracle while the others keep their own version — it is generic over the
//! provider so the `JsonProvider`-backed file shares it too.
//!
//! [`assert_normalizes_preserving`] and [`assert_normalizes_preserving_in`] are
//! the assertion half of that, one per shuffle direction, for the same reason.
//!
//! [`sweep_sequences`] is here for the same reason: the three exhaustive sweeps
//! assert case-count floors, and one of them an exact residual count, over the
//! sequences it generates. Three private copies of that generator would be
//! three ways for one sweep's corpus to drift out from under its pinned
//! numbers unnoticed.

use ferro_hgvs::reference::ReferenceProvider;
use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{
    parse_hgvs, HgvsVariant, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection,
};

/// The deterministic 20-mers every cis-allele sweep enumerates over.
///
/// `2 * seeds` sequences: for each seed, one over `{A,T}` — the alphabet that
/// builds the long alternating tracts a shift can travel through — and one over
/// `{A,C,G,T}`. The generator is a plain xorshift64 seeded by the golden-ratio
/// odd multiple of the seed, so the corpus is reproducible without a dependency
/// and identical across the sweeps that share it.
///
/// The corpus is load-bearing: each sweep's case total is asserted against a
/// floor, and `FIVE_PRIME_DUP_DEL_SEQUENCE_CHANGES` in
/// `cis_junction_crossing_shift.rs` is an exact count over *these* sequences.
/// Changing the generator or the seed count re-rolls the corpus and moves those
/// numbers. The order does not: it decides only which few failing cases a
/// failure message samples.
pub fn sweep_sequences(seeds: u32) -> Vec<String> {
    let mut sequences = Vec::with_capacity(2 * seeds as usize);
    for seed in 0..seeds {
        for alphabet in [b"AT".as_slice(), b"ACGT".as_slice()] {
            let mut state = u64::from(seed).wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
            sequences.push(
                (0..20)
                    .map(|_| {
                        state ^= state << 13;
                        state ^= state >> 7;
                        state ^= state << 17;
                        alphabet[(state % alphabet.len() as u64) as usize] as char
                    })
                    .collect(),
            );
        }
    }
    sequences
}

/// A single-contig genomic provider under the accession `TEMPLATE`.
pub fn provider(seq: &str) -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("TEMPLATE", seq.to_string());
    provider
}

/// Normalize `input` against `seq` in the default 3' direction.
pub fn normalize(seq: &str, input: &str) -> String {
    normalize_in(seq, input, ShuffleDirection::ThreePrime)
}

/// Normalize `input` against `seq` in an explicit direction.
pub fn normalize_in(seq: &str, input: &str, direction: ShuffleDirection) -> String {
    let normalizer = Normalizer::with_config(
        provider(seq),
        NormalizeConfig::default().with_direction(direction),
    );
    let variant = parse_hgvs(input).expect("parse");
    format!("{}", normalizer.normalize(&variant).expect("normalize"))
}

/// Apply `descriptor` to `seq` **independently of the normalizer**, against the
/// single-contig `TEMPLATE` provider this module builds.
///
/// See [`apply_with`] for the contract; this is the convenience wrapper for the
/// common case where the reference is just a string.
pub fn apply(seq: &str, descriptor: &str) -> Option<String> {
    apply_with(&provider(seq), seq, descriptor)
}

/// Apply `descriptor` to `reference` **independently of the normalizer**, by
/// converting each member to its SPDI triple and splicing the reference.
///
/// Takes the provider and reference explicitly so every caller shares this one
/// walk: `claimed_from` below is what makes "declines an overlapping
/// description" true, and a second copy of it is a second thing to drift. (The
/// generic is presentational today — `MockProvider` is a type alias of
/// `JsonProvider`, so all callers pass one concrete type — but the parameters
/// are what let the padded-synthetic and tempfile-JSON tests reuse it.)
///
/// Returns `None` when a member cannot be converted, when a stated deletion
/// disagrees with the reference, or when two members claim the same base — an
/// overlapping description has no single well-defined resulting sequence, and
/// silently double-splicing one would invent a sequence neither side denotes.
pub fn apply_with<P: ReferenceProvider + ?Sized>(
    provider: &P,
    reference: &str,
    descriptor: &str,
) -> Option<String> {
    let members: Vec<HgvsVariant> = match parse_hgvs(descriptor).ok()? {
        HgvsVariant::Allele(allele) => allele.variants.clone(),
        single => vec![single],
    };
    let mut triples = Vec::with_capacity(members.len());
    for member in &members {
        triples.push(hgvs_to_spdi(member, provider).ok()?);
    }
    // 3'→5' so an applied splice never shifts a later one's coordinates.
    //
    // The tie-break matters. `sort_by_key` is stable, so members sharing a
    // position would otherwise keep input order, and a zero-length member (a
    // pure insertion) applied first would set `claimed_from` to that position
    // and make a deletion starting there look like an overlap. Longer deletions
    // first means a span-claiming member is always applied before the
    // zero-length one it is flush against, which is not an overlap.
    triples.sort_by_key(|t| {
        (
            std::cmp::Reverse(t.position),
            std::cmp::Reverse(t.deletion.len()),
        )
    });
    let reference = reference.as_bytes();
    let mut edited = reference.to_vec();
    let mut claimed_from = reference.len();
    let mut insertion_at: Option<usize> = None;
    for triple in &triples {
        let start = usize::try_from(triple.position).ok()?;
        let end = start.checked_add(triple.deletion.len())?;
        if end > reference.len() || end > claimed_from {
            return None; // out of bounds, or overlapping the member 3' of it
        }
        // Two pure insertions at one interbase have no defined order between
        // them, so the sequence they jointly denote is ambiguous — decline it
        // rather than let the arbitrary sort order pick a winner. The
        // `claimed_from` walk above cannot see this: a zero-length member never
        // advances the claim, so both would pass it.
        if triple.deletion.is_empty() && insertion_at == Some(start) {
            return None;
        }
        if !reference[start..end].eq_ignore_ascii_case(triple.deletion.as_bytes()) {
            return None;
        }
        edited.splice(start..end, triple.insertion.bytes());
        if triple.deletion.is_empty() {
            insertion_at = Some(start);
        }
        claimed_from = start;
    }
    String::from_utf8(edited).ok()
}

/// The half-open interbase spans of `descriptor`'s members, in the order they
/// are rendered.
///
/// Spans come from `hgvs_to_spdi` rather than the printed HGVS endpoints
/// because the two disagree for the edits that consume no base: `a_b ins` names
/// two positions while occupying the junction at `a`, and `a_b dup` places its
/// copy at the junction *after* `b`. Interbase coordinates are the frame the
/// ordering and disjointness contracts are actually stated in — reading the
/// printed endpoints instead is what let #1261 render as ascending.
pub fn member_interbase_spans(seq: &str, descriptor: &str) -> Vec<(u64, u64)> {
    let template = provider(seq);
    let members: Vec<HgvsVariant> = match parse_hgvs(descriptor).expect("parse") {
        HgvsVariant::Allele(allele) => allele.variants.clone(),
        single => vec![single],
    };
    members
        .iter()
        .map(|member| {
            let triple = hgvs_to_spdi(member, &template)
                .unwrap_or_else(|e| panic!("`{descriptor}`: member `{member}` has no SPDI: {e}"));
            (
                triple.position,
                triple.position + triple.deletion.len() as u64,
            )
        })
        .collect()
}

/// Assert that `descriptor`'s members are rendered in ascending interbase
/// order (#1098's contract, and #1235's second acceptance criterion).
pub fn assert_members_ascending(seq: &str, descriptor: &str) {
    let spans = member_interbase_spans(seq, descriptor);
    assert!(
        spans.windows(2).all(|pair| pair[0] <= pair[1]),
        "`{descriptor}` renders its members out of order: interbase spans {spans:?}"
    );
}

/// Assert that `input` normalizes to `expected` in the default 3' direction
/// *and* that both denote the same sequence when applied to `seq`.
pub fn assert_normalizes_preserving(seq: &str, input: &str, expected: &str) {
    assert_normalizes_preserving_in(seq, input, expected, ShuffleDirection::ThreePrime);
}

/// Assert that `input` normalizes to `expected` in `direction` *and* that both
/// denote the same sequence when applied to `seq`.
///
/// The 5' direction gets the same single home as the 3' one for the reason
/// [`apply_with`] does: a hand-rolled copy per file is a way for one file's
/// assertion to weaken while the others keep theirs.
pub fn assert_normalizes_preserving_in(
    seq: &str,
    input: &str,
    expected: &str,
    direction: ShuffleDirection,
) {
    let actual = normalize_in(seq, input, direction);
    assert_eq!(actual, expected, "normalizing {input} in {direction:?}");
    let from_input = apply(seq, input).expect("input applies");
    let from_output = apply(seq, &actual).unwrap_or_else(|| {
        panic!("{actual} has no well-defined resulting sequence (overlapping or unconvertible)")
    });
    assert_eq!(
        from_output, from_input,
        "{input} -> {actual} changed the sequence"
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Pins [`sweep_sequences`]'s shape and its first sequence, so a change to
    /// the generator fails here — naming the generator — rather than surfacing
    /// as a moved residual count in a 76-second sweep that names normalization.
    #[test]
    fn sweep_sequences_is_a_stable_corpus() {
        let sequences = sweep_sequences(4);
        assert_eq!(sequences.len(), 8, "two alphabets per seed");
        assert!(
            sequences.iter().all(|s| s.len() == 20),
            "every sequence is a 20-mer: {sequences:?}"
        );
        // Alphabets alternate within a seed: even index {A,T}, odd {A,C,G,T}.
        assert!(
            sequences
                .iter()
                .step_by(2)
                .all(|s| s.bytes().all(|b| b == b'A' || b == b'T')),
            "even indices are the {{A,T}} corpus: {sequences:?}"
        );
        assert!(
            sequences
                .iter()
                .all(|s| s.bytes().all(|b| b"ACGT".contains(&b))),
            "every sequence is over {{A,C,G,T}}: {sequences:?}"
        );
        // Seed 0's pair, hardcoded. Checking the alphabets by predicate alone
        // would not distinguish the two halves — a generator that drew every
        // sequence from {A,T} satisfies both `all` assertions above. These two
        // values are what actually pin the xorshift and the alphabet split.
        // The first is the nine-`T` tract the pinned cases use as `TRACT`, which
        // is how the sweeps reach repeat notation at all.
        assert_eq!(sequences[0], "TTTTTTTTTAATATATTTTA");
        assert_eq!(sequences[1], "CCCCCCCCTGACGTATCCTA");
    }

    /// The sweeps re-derive their corpus on every run, so the generator must be
    /// a pure function of its seed count.
    #[test]
    fn sweep_sequences_is_deterministic() {
        assert_eq!(sweep_sequences(8), sweep_sequences(8));
        // A larger seed count extends the corpus rather than re-rolling it,
        // which is what lets one sweep use 48 seeds and another 64.
        assert_eq!(sweep_sequences(8)[..], sweep_sequences(16)[..16]);
    }
}
