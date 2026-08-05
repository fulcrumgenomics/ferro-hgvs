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

/// Seeds a sweep should draw, given the count it wants when run in full (#1295).
///
/// The exhaustive sweeps are 92% of a local `cargo nextest run`: measured on this
/// commit's parent, `cis_junction_crossing_shift` alone is 79.6s of an 86.6s
/// suite, and `repeat_span_sibling_overlap` adds 33.1s. That is a poor edit-test
/// loop for work that touches none of it, and it is the only cost — CI shards the
/// `Test` job, so it does not gate CI wall-clock.
///
/// So the default is a small prefix of the corpus and CI asks for all of it:
///
/// | `FERRO_SWEEP_SEEDS` | seeds drawn |
/// |---|---|
/// | unset | [`DEFAULT_SWEEP_SEEDS`], or `full` if that is smaller |
/// | `full` | `full` |
/// | a number | that number |
///
/// **Why a prefix, and why this does not weaken the pinned counts.**
/// [`sweep_sequences`] is prefix-stable — `sweep_sequences(8)` is exactly the
/// first sixteen entries of `sweep_sequences(16)`, asserted by
/// `sweep_sequences_is_prefix_stable` below. So a reduced run enumerates a strict
/// *subset* of the full run's cases, which is what makes the assertions survive
/// unchanged rather than becoming profile-dependent:
///
/// - `overlapping.is_empty()` / `changed.is_empty()` — a property that holds over
///   a set holds over any subset, so a reduced run cannot pass one the full run
///   would fail. It can only find *fewer* failures, never different ones.
/// - `FIVE_PRIME_DUP_DEL_SEQUENCE_CHANGES` — pinned at **zero**, and zero is
///   subset-stable for the same reason. Had it been a non-zero count, this knob
///   would have had to make it profile-dependent, which #1295 flags as its own
///   hazard; it is only safe because the residual was driven to zero first.
///
/// Two assertions are **not** covered by that subset argument, and it is worth
/// being exact about which:
///
/// - each sweep's case-count floor is genuinely seed-dependent, which is why
///   those are expressed per-seed rather than as absolutes;
/// - the `skipped * 10 < checked` share is an **aggregate ratio**, and a ratio
///   is not subset-stable. Adding seeds changes which cases are skipped, so a
///   prefix run clearing the bound does not prove the full run will: later seeds
///   could be skip-heavier and push the aggregate over. It is not profile-*noisy*
///   — both runs assert the same bound — but a green local prefix is not a
///   guarantee here the way it is for the three above. CI's `full` run is what
///   actually settles it.
///
/// The escape hatch is the numeric form: `FERRO_SWEEP_SEEDS=12` when bisecting a
/// failure that a reduced run reproduces but a full one takes too long to iterate
/// on.
pub fn sweep_seeds(full: u32) -> u32 {
    let setting = match std::env::var("FERRO_SWEEP_SEEDS") {
        Ok(value) => Some(value),
        Err(std::env::VarError::NotPresent) => None,
        // Set, but not valid UTF-8. Folding this into "unset" would silently run
        // the 4-seed prefix for someone who asked for the full corpus and typoed
        // it — the same "ran less than you think" failure the panic below exists
        // to prevent, so it gets the same treatment rather than a quiet default.
        Err(std::env::VarError::NotUnicode(raw)) => panic!(
            "FERRO_SWEEP_SEEDS is set to a non-UTF-8 value ({raw:?}); it must be \
             `full` or a seed count. Unset it for the default \
             {DEFAULT_SWEEP_SEEDS}-seed prefix."
        ),
    };
    sweep_seeds_from(setting.as_deref(), full)
}

/// The selection rule of [`sweep_seeds`], as a pure function of the setting.
///
/// Split out so the rule can be tested **without touching the environment**.
/// The tests used to `set_var`/`remove_var` around their assertions, guarded by
/// a SAFETY note arguing that nextest gives each test its own process. That
/// argument is sound under nextest and *false* under `cargo test`, which
/// `CLAUDE.md` lists as a supported alternative and which runs tests as threads
/// in one process: there, setting `FERRO_SWEEP_SEEDS=12` mid-run races every
/// concurrently-executing sweep reading the same variable, silently shrinking
/// its corpus. A test that can quietly hollow out the sweeps it shares a process
/// with is not worth the coverage, and none is lost by parsing a string instead.
fn sweep_seeds_from(setting: Option<&str>, full: u32) -> u32 {
    match setting {
        Some(value) if value.eq_ignore_ascii_case("full") => full,
        Some(value) => value.trim().parse().unwrap_or_else(|_| {
            panic!(
                "FERRO_SWEEP_SEEDS must be `full` or a seed count, got {value:?}. \
                 Unset it for the default {DEFAULT_SWEEP_SEEDS}-seed prefix."
            )
        }),
        None => DEFAULT_SWEEP_SEEDS.min(full),
    }
}

/// Seeds drawn when `FERRO_SWEEP_SEEDS` is unset.
///
/// Four keeps every *shape* the sweeps enumerate — the loop bounds over member
/// spellings, positions and directions are untouched — and reduces only sequence
/// diversity. That is the right axis to cut: each sweep's own notes record that
/// every blocking defect found so far lived in a shape the generator could not
/// emit rather than in a sequence it did not draw.
pub const DEFAULT_SWEEP_SEEDS: u32 = 4;

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
    apply_reason(provider, reference, descriptor).ok()
}

/// Why [`apply_with`] declined a descriptor.
///
/// `apply_with` collapses all of these to `None`, which is the right shape for a
/// sweep asking "does this denote a sequence at all". It is the wrong shape for
/// reporting: #1268 records a three-member measurement whose 61,765 "unapplicable"
/// cases conflated genuine overlap with conversion failure and with parse
/// rejection, and notes that no number should be quoted until the causes are
/// separated. These are those causes.
///
/// Ordered by severity as a sweep should read them: [`Self::Overlapping`] and
/// [`Self::CoincidentInsertions`] mean the description denotes no single
/// sequence, which for a *normalizer output* is a defect. The rest are limits of
/// the harness or of the input.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApplyFailure {
    /// `parse_hgvs` rejected it. Suite-wide, `FERRO_ASSERT_REPARSE` already
    /// covers this for normalizer output (#1259), so a sweep seeing it here is
    /// usually looking at a generated *input*.
    Unparseable,
    /// A member has no SPDI triple — `hgvs_to_spdi` declined it.
    Unconvertible,
    /// A member's span runs past the end of the reference.
    OutOfBounds,
    /// A member claims a base another member already claimed.
    Overlapping,
    /// Two pure insertions at one interbase, whose relative order — and so the
    /// sequence they jointly denote — is undefined.
    CoincidentInsertions,
    /// A member's stated deleted bases disagree with the reference.
    StatedBasesMismatch,
}

/// [`apply_with`], but reporting *why* it declined.
///
/// The two share this one walk rather than the walk being duplicated, for the
/// reason `apply_with`'s own note gives about `claimed_from`: a second copy is a
/// second thing to drift.
pub fn apply_reason<P: ReferenceProvider + ?Sized>(
    provider: &P,
    reference: &str,
    descriptor: &str,
) -> Result<String, ApplyFailure> {
    let members: Vec<HgvsVariant> = match parse_hgvs(descriptor) {
        Ok(HgvsVariant::Allele(allele)) => allele.variants.clone(),
        Ok(single) => vec![single],
        Err(_) => return Err(ApplyFailure::Unparseable),
    };
    let mut triples = Vec::with_capacity(members.len());
    for member in &members {
        triples.push(hgvs_to_spdi(member, provider).map_err(|_| ApplyFailure::Unconvertible)?);
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
        let start = usize::try_from(triple.position).map_err(|_| ApplyFailure::OutOfBounds)?;
        let end = start
            .checked_add(triple.deletion.len())
            .ok_or(ApplyFailure::OutOfBounds)?;
        // Split, where `apply_with` used to test both in one condition and
        // report them alike. They are different severities: running off the
        // reference is a coordinate that does not exist, while claiming a base
        // the member 3' of it already claimed is the overlap #1268 asks to have
        // counted on its own.
        if end > reference.len() {
            return Err(ApplyFailure::OutOfBounds);
        }
        if end > claimed_from {
            return Err(ApplyFailure::Overlapping);
        }
        // Two pure insertions at one interbase have no defined order between
        // them, so the sequence they jointly denote is ambiguous — decline it
        // rather than let the arbitrary sort order pick a winner. The
        // `claimed_from` walk above cannot see this: a zero-length member never
        // advances the claim, so both would pass it.
        if triple.deletion.is_empty() && insertion_at == Some(start) {
            return Err(ApplyFailure::CoincidentInsertions);
        }
        if !reference[start..end].eq_ignore_ascii_case(triple.deletion.as_bytes()) {
            return Err(ApplyFailure::StatedBasesMismatch);
        }
        edited.splice(start..end, triple.insertion.bytes());
        if triple.deletion.is_empty() {
            insertion_at = Some(start);
        }
        claimed_from = start;
    }
    String::from_utf8(edited).map_err(|_| ApplyFailure::Unconvertible)
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
    }

    /// A larger seed count **extends** the corpus rather than re-rolling it.
    ///
    /// This is what lets one sweep use 48 seeds and another 64, and since #1295
    /// it is load-bearing for [`sweep_seeds`]: a reduced run must enumerate a
    /// strict subset of the full run's cases, or the sweeps' `is_empty()` and
    /// pinned-zero assertions could pass at the default seed count while failing
    /// at the full one. Split out of the determinism test and named, because a
    /// property the seed knob's soundness rests on should fail under its own name.
    #[test]
    fn sweep_sequences_is_prefix_stable() {
        assert_eq!(sweep_sequences(8)[..], sweep_sequences(16)[..16]);
        // The default prefix specifically, since that is the one every local run
        // takes.
        let full = sweep_sequences(64);
        let default = sweep_sequences(DEFAULT_SWEEP_SEEDS);
        assert_eq!(default[..], full[..default.len()]);
    }

    /// The selection rule, exercised through [`sweep_seeds_from`] so nothing
    /// here mutates a process-wide variable the sweeps are concurrently reading.
    #[test]
    fn sweep_seeds_defaults_to_the_prefix_and_honours_full() {
        assert_eq!(sweep_seeds_from(None, 64), DEFAULT_SWEEP_SEEDS);
        // A sweep asking for fewer than the default gets its own count, never more.
        assert_eq!(sweep_seeds_from(None, 2), 2);

        assert_eq!(sweep_seeds_from(Some("full"), 64), 64);
        assert_eq!(
            sweep_seeds_from(Some("FULL"), 48),
            48,
            "the sentinel is case-insensitive"
        );

        assert_eq!(
            sweep_seeds_from(Some("12"), 64),
            12,
            "an explicit count overrides both"
        );
        assert_eq!(
            sweep_seeds_from(Some("  12  "), 64),
            12,
            "surrounding whitespace is tolerated, as a shell can easily add it"
        );
    }

    /// A value that is neither `full` nor a count must be loud, not defaulted —
    /// silently running 4 seeds for someone who asked for 48 is the failure this
    /// whole knob has to avoid.
    #[test]
    #[should_panic(expected = "FERRO_SWEEP_SEEDS must be `full` or a seed count")]
    fn an_unparseable_sweep_seed_setting_panics() {
        let _ = sweep_seeds_from(Some("most"), 64);
    }
}
