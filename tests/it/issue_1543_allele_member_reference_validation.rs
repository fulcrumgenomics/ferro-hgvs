//! #1543 — a cis-allele member's stated reference bases must be validated
//! whether or not a sibling is close enough to merge with it.
//!
//! # The question
//!
//! Strict mode rejects a substitution that misstates its reference base. It did
//! **not** reject the same substitution when a sibling sat close enough for the
//! two to merge. Measured on `main` at `5616cdb9` against the prepared
//! reference, where `c.1 = A`, `c.3 = G` and `c.500 = G`:
//!
//! ```text
//! NM_000143.3:c.1G>T            ->  normalize_error  Reference mismatch at 64-64
//! NM_000143.3:c.[1G>T;500G>C]   ->  normalize_error  Reference mismatch at 64-64
//! NM_000143.3:c.[1G>T;3T>A]     ->  NM_000143.3:c.1_3delinsTTA   status=ok
//! ```
//!
//! The discriminator was **merge distance and nothing else**: the identical
//! false `G` at `c.1` was caught with a sibling 499 nt away and missed with one
//! 2 nt away. Worse than an ordinary validation gap, because the accepted output
//! is built from the *real* bases — so ferro silently rewrote a description whose
//! reference assertions were false into a well-formed, correct-looking one, with
//! no warning, in the mode whose whole purpose is to reject such input.
//!
//! # The ruling: the verdict is a property of the member, not of its neighbours
//!
//! `DNA/substitution.md:26` reads a stated reference base as a claim about the
//! reference — "a substitution of the `C` nucleotide at `g.33038255` by an `A`"
//! — so a member whose stated base contradicts the reference is not a
//! description of that reference at all, however its siblings are spelled.
//!
//! The operative authority here is not a spec conflict, though; it is ferro's
//! own already-adjudicated contract, recorded verbatim at
//! `merge::canonicalize_from_sequence` since #1052/#1097:
//!
//! > A member that misstates its reference base is a reference mismatch, not a
//! > canonicalization problem: strict mode must still reject it and lenient mode
//! > must still warn.
//!
//! This file pins that contract against the one input property that had been
//! deciding it: the distance to the nearest sibling.
//!
//! # The mechanism, which is *not* the one the issue names
//!
//! The issue attributes the defect to `stated_reference_bases_match` being
//! vacuous at `merge.rs`'s canonicalization guard. That guard is real and its
//! doc comment does predict this class of failure, but it is not what let these
//! inputs through — measured directly, and reported on the PR.
//!
//! `Normalizer::normalize_allele` runs `merge::collapse_overlapping_cis_edits`
//! and `merge::merge_consecutive_edits` over the **raw** members *before* any
//! member reaches `normalize_na_edit`, which is the only place `RefSeqMismatch`
//! (W5001) is raised. A merge keeps each member's *alt* bases and nothing else —
//! `merge::Anchor` has no channel for an asserted reference — so a member
//! consumed by a merge had its assertions discarded before anything looked at
//! them. The merged anchor was then rebuilt from the real bases, which is why
//! the output looked correct. Everything downstream, the canonicalization guard
//! included, was seeing a member that no longer stated anything.
//!
//! Two consequences the issue's framing understates, both pinned below:
//!
//! * it is not only the codon-frame merge (`general.md:35`, gap of one). Plain
//!   **strict adjacency** loses the assertions just the same, so
//!   `c.[1G>T;2A>C]` was accepted too;
//! * it is not only substitutions, and not only `c.`. Every channel that can
//!   carry an asserted reference — `delACG`, `delACinsGT`, `dupCA`, `invCCCCC` —
//!   and every axis that merges — `g.`, `m.`, `c.`, `n.`, `r.` — was affected.
//!
//! # The fix, and why it is not at the merge site
//!
//! This is the sixth defect in one family: #1052 (substitutions), #1068 (the
//! `m.` axis), #1092 (the parser discarded the stated bases), #1097 (multi-base
//! range substitutions), #1352 (`NaEdit::Identity` omitted while the doc comment
//! claimed the channel list was complete). Each of the five was fixed one shape
//! at a time.
//!
//! So the question is asked on the **authored input**, in `normalize_core`,
//! before any pass has had the chance to strip anything
//! (`merge::authored_member_reference_mismatches`). That is ordering-immune by
//! construction — no rearrangement of the merge pipeline can make it vacuous —
//! and it reuses `normalize::validate::validate_reference`, the same function
//! the per-member pipeline uses, so the channel list cannot drift between the
//! two sites the way #1352's did.
//!
//! It adds warnings and changes nothing else. A correctly-spelled allele merges
//! exactly as it did before, which is what
//! [`the_correctly_spelled_sibling_still_merges`] exists to hold.

use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::normalize::{NormalizationWarning, NormalizeConfig};
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Test sequence, shared by every fixture below.
///
/// ```text
///   axis:  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 ...  30
///   base:  A T G C A A A A A  C  C  C  C  C  G  G ...   C
/// ```
///
/// The first three bases reproduce the issue's transcript exactly — `c.1 = A`,
/// `c.3 = G` — so `1G>T` and `3T>A` misstate the reference in the same way the
/// reported `NM_000143.3` rows do, and `c.1_3delinsTTA` is the same accepted
/// output. CDS = `c.1..=c.60`, so `c.N`, `n.N` and `r.N` all name transcript
/// position `N`; codon 1 is `c.1..3`, which is what makes `c.1` and `c.3`
/// eligible for the `general.md:35` codon-frame merge.
const CORE: &str = "ATGCAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGT";

/// The transcript accession every `c.`/`n.`/`r.` fixture uses.
const TX: &str = "NM_TEST.1";

/// A hermetic provider carrying [`CORE`] as a coding transcript, a plain contig
/// (`CTG`, for `g.`) and a mitochondrial contig (`MITO`, for `m.`).
///
/// Synthetic on purpose: a test that needed `FERRO_MANIFEST` would be skipped on
/// every PR and so would gate nothing.
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let length = CORE.len() as u64;
    provider.add_transcript(Transcript::new(
        TX.to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        CORE.to_string(),
        Some(1),
        Some(60),
        vec![Exon::new(1, 1, length)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    provider.add_genomic_sequence("CTG", CORE.to_string());
    provider.add_genomic_sequence("MITO", CORE.to_string());
    provider
}

/// Whether strict mode rejects `input`.
fn strict_rejects(input: &str) -> bool {
    Normalizer::with_config(
        provider(),
        NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
    )
    .normalize(&parse_hgvs(input).expect("fixture must parse"))
    .is_err()
}

/// `input`'s lenient normalization, as a string.
fn lenient(input: &str) -> String {
    Normalizer::new(provider())
        .normalize(&parse_hgvs(input).expect("fixture must parse"))
        .expect("lenient mode must not reject")
        .to_string()
}

/// Every `(stated_ref, actual_ref, position)` reference mismatch lenient mode
/// reports for `input`.
fn reference_mismatches(input: &str) -> Vec<(String, String, String)> {
    Normalizer::new(provider())
        .normalize_with_diagnostics(&parse_hgvs(input).expect("fixture must parse"))
        .expect("lenient mode must not reject")
        .warnings
        .iter()
        .filter_map(|w| match w {
            NormalizationWarning::RefSeqMismatch {
                stated_ref,
                actual_ref,
                position,
                ..
            } => Some((stated_ref.clone(), actual_ref.clone(), position.clone())),
            _ => None,
        })
        .collect()
}

// ---------------------------------------------------------------------------
// The reported defect, and the two controls that keep a fix from being credited
// to breaking something else.
// ---------------------------------------------------------------------------

/// The merged case: `c.[1G>T;3T>A]` is rejected in strict mode.
///
/// `c.1` is `A` and `c.3` is `G`, so both members misstate. They are one
/// nucleotide apart inside codon 1, which is what `general.md:35` merges — and
/// what used to make the two false assertions disappear into
/// `c.1_3delinsTTA` with `status=ok`.
#[test]
fn the_merged_case_is_rejected_in_strict_mode() {
    assert!(
        strict_rejects(&format!("{TX}:c.[1G>T;3T>A]")),
        "both members misstate the reference (c.1 is A, c.3 is G); strict mode \
         must reject them whether or not they are close enough to merge (#1543)"
    );
}

/// The unmerged control: `c.[1G>T;30G>C]` was *already* rejected, and still is.
///
/// The stand-in for the issue's `c.[1G>T;500G>C]`. It passed before the fix, so
/// a fix cannot be credited to it — its job is to fail if the change ever
/// stops rejecting the case that always worked.
#[test]
fn the_unmerged_control_is_still_rejected() {
    assert!(
        strict_rejects(&format!("{TX}:c.[1G>T;30G>C]")),
        "a member 29 nt from its sibling reaches the per-member validator and \
         was already rejected; that must not regress (#1543)"
    );
}

/// The correctly-spelled sibling still merges, and still merges to the same
/// string.
///
/// This is the regression the fix is most likely to cause: over-rejecting, or
/// refusing the merge, for input whose reference assertions are *true*.
/// `c.1 = A` and `c.3 = G`, so `c.[1A>T;3G>A]` states them correctly and must
/// reach `c.1_3delinsTTA` — byte for byte the same output the defective build
/// produced for the *mis*-spelled pair, which is exactly why the two are
/// indistinguishable downstream and why the mis-spelled one has to be rejected.
#[test]
fn the_correctly_spelled_sibling_still_merges() {
    let input = format!("{TX}:c.[1A>T;3G>A]");
    assert!(
        !strict_rejects(&input),
        "`{input}` states c.1 = A and c.3 = G, which is what the reference \
         holds; strict mode must accept it (#1543)"
    );
    assert_eq!(
        lenient(&input),
        format!("{TX}:c.1_3delinsTTA"),
        "a correctly-spelled allele must merge exactly as it did before the \
         #1543 guard was added"
    );
    assert!(
        reference_mismatches(&input).is_empty(),
        "`{input}` misstates nothing, so it must carry no RefSeqMismatch (#1543)"
    );
}

/// Lenient mode warns and continues rather than rejecting.
///
/// Both halves matter. The warning is the one the per-member pipeline already
/// raises for a lone misstating substitution (W5001 `RefSeqMismatch`, reused
/// rather than reinvented), reported once per false member at the member's own
/// position. And normalization still *completes*, producing the same string it
/// produced before — the fix adds a diagnostic, it does not move output.
#[test]
fn lenient_mode_warns_and_continues() {
    let input = format!("{TX}:c.[1G>T;3T>A]");
    assert_eq!(
        reference_mismatches(&input),
        vec![
            ("G".to_string(), "A".to_string(), "1-1".to_string()),
            ("T".to_string(), "G".to_string(), "3-3".to_string()),
        ],
        "lenient mode must report both false assertions, each at its own \
         position, exactly as the lone-substitution path does (#1543)"
    );
    assert_eq!(
        lenient(&input),
        format!("{TX}:c.1_3delinsTTA"),
        "lenient mode warns and continues; the normalized form is unchanged"
    );
}

/// A member that reaches the per-member validator is not reported twice.
///
/// The guard runs on the authored members, and the per-member pipeline reports
/// the same finding for any member a merge did not consume. `c.[1G>T;30G>C]`
/// merges nothing, so both members reach both — and must still yield exactly two
/// warnings.
#[test]
fn an_unmerged_member_is_not_reported_twice() {
    assert_eq!(
        reference_mismatches(&format!("{TX}:c.[1G>T;30G>C]")),
        vec![
            ("G".to_string(), "A".to_string(), "1-1".to_string()),
            ("G".to_string(), "C".to_string(), "30-30".to_string()),
        ],
        "one mismatch per false member; the authored-member guard must not \
         duplicate what the per-member pipeline already reported (#1543)"
    );
}

// ---------------------------------------------------------------------------
// The invariant: distance decides nothing.
// ---------------------------------------------------------------------------

/// Sibling positions spanning the merge threshold, with the base each one
/// truthfully states.
///
/// `c.2 = T` (strictly adjacent, merges), `c.3 = G` (gap of one inside codon 1,
/// merges via `general.md:35`), `c.4 = C` (gap of two, `general.md:34` splits),
/// `c.5 = A`, `c.10 = C`, `c.30 = C` (far).
const SIBLINGS: [(u32, char); 6] = [(2, 'T'), (3, 'G'), (4, 'C'), (5, 'A'), (10, 'C'), (30, 'C')];

/// The verdict on a false `c.1G>T` does not depend on how far away its sibling
/// sits.
///
/// This is the actual invariant — the defect's whole discriminator was merge
/// distance — and it is what keeps the guard alive through a refactor of the
/// merge threshold: a change to `general.md:34`/`:35` eligibility moves which
/// rows merge, and this test does not care which those are.
///
/// The sweep is paired with a correctly-spelled control at every distance, for
/// two reasons. It proves the rejections come from the false `c.1`, not from
/// something about the sibling; and the controls' outputs are what show the
/// sweep genuinely *spans* the threshold — a sweep entirely on one side of it
/// would assert the invariant without ever testing it, which is how a distance
/// sweep quietly stops covering anything.
#[test]
fn the_verdict_does_not_depend_on_merge_distance() {
    let mut merged = 0usize;
    let mut split = 0usize;
    for (position, reference_base) in SIBLINGS {
        // `c.1` is `A`, so `1G>T` is false at every distance.
        let false_at_one = format!("{TX}:c.[1G>T;{position}{reference_base}>N]");
        assert!(
            strict_rejects(&false_at_one),
            "`{false_at_one}` misstates c.1 as G when the reference holds A; \
             a sibling {} nt away must not change that verdict (#1543)",
            position - 1
        );

        // The same allele with `c.1` stated truthfully. Accepted at every
        // distance, which is what makes the rejections above attributable to
        // the false base rather than to the pair.
        let true_at_one = format!("{TX}:c.[1A>T;{position}{reference_base}>N]");
        assert!(
            !strict_rejects(&true_at_one),
            "`{true_at_one}` states every reference base correctly and must be \
             accepted at every distance (#1543)"
        );

        if lenient(&true_at_one).contains('[') {
            split += 1;
        } else {
            merged += 1;
        }
    }
    assert!(
        merged > 0 && split > 0,
        "the sweep must straddle the merge threshold, or it asserts the \
         distance-independence invariant without exercising it — got \
         {merged} merged and {split} split of {} distances (#1543)",
        SIBLINGS.len()
    );
}

/// Strict adjacency loses the assertions too, not only the codon-frame merge.
///
/// `c.[1G>T;2A>C]` is the shape the issue's reproducer does not cover: `c.1 = A`
/// and `c.2 = T`, both misstated, and the two members are *consecutive*, so they
/// merge under `delins.md:16` rather than under the `general.md:35` codon
/// exception. It was accepted as `c.1_2delinsTC` before the fix.
#[test]
fn strictly_adjacent_members_are_rejected_too() {
    assert!(
        strict_rejects(&format!("{TX}:c.[1G>T;2A>C]")),
        "consecutive members merge under delins.md:16 and lost their reference \
         assertions the same way the codon-frame pair did (#1543)"
    );
    assert!(
        !strict_rejects(&format!("{TX}:c.[1A>T;2T>C]")),
        "the same pair spelled correctly must still be accepted (#1543)"
    );
}

// ---------------------------------------------------------------------------
// Coverage: every channel that can assert a reference, and every axis that
// merges. #1352's lesson was that a channel list which claims to be complete
// and is not fails silently, so each one is pinned rather than assumed.
// ---------------------------------------------------------------------------

/// Every `NaEdit` channel that can carry an asserted reference is validated.
///
/// Each row pairs a member misstating its reference with a sibling close enough
/// to merge. The `del`/`delins`/`dup`/`inv` channels were all accepted before
/// the fix; a substitution is included as the shape the issue reported.
#[test]
fn every_stated_reference_channel_is_validated() {
    // `c.1_2` is `AT`, `c.4_5` is `CA`, `c.10_14` is `CCCCC`.
    for input in [
        format!("{TX}:c.[1G>T;3T>A]"),           // substitution
        format!("{TX}:c.[1_2delAA;3T>A]"),       // deletion, stated bases
        format!("{TX}:c.[1_2delAAinsGG;3T>A]"),  // delins, stated deleted bases
        format!("{TX}:c.[4_5dupAA;6A>G]"),       // duplication, stated bases
        format!("{TX}:c.[10_14invGGGGG;15G>T]"), // inversion, stated bases
    ] {
        assert!(
            strict_rejects(&input),
            "`{input}` asserts reference bases the reference contradicts; every \
             channel that can carry them must be validated, not just \
             substitutions (#1352, #1543)"
        );
    }

    // The same five shapes spelled truthfully, which must all still be accepted.
    for input in [
        format!("{TX}:c.[1A>T;3G>A]"),
        format!("{TX}:c.[1_2delAT;3G>A]"),
        format!("{TX}:c.[1_2delATinsGG;3G>A]"),
        format!("{TX}:c.[4_5dupCA;6A>G]"),
        format!("{TX}:c.[10_14invCCCCC;15G>T]"),
    ] {
        assert!(
            !strict_rejects(&input),
            "`{input}` states its reference bases correctly and must be \
             accepted (#1543)"
        );
    }
}

/// Every axis whose members merge is covered — including `m.` (#1068) and `r.`.
///
/// #1068 is in this family precisely because a fix applied on one axis and not
/// another is the same defect twice. The guard runs on the shared cis-member
/// path, so it reaches all five; each is pinned so a future axis-specific
/// narrowing fails here.
#[test]
fn every_merging_axis_is_covered() {
    for (input, corrected) in [
        (format!("{TX}:c.[1G>T;2A>C]"), format!("{TX}:c.[1A>T;2T>C]")),
        (format!("{TX}:n.[1G>T;2A>C]"), format!("{TX}:n.[1A>T;2T>C]")),
        (format!("{TX}:r.[1g>u;3u>a]"), format!("{TX}:r.[1a>u;3g>a]")),
        (
            "CTG:g.[1G>T;2A>C]".to_string(),
            "CTG:g.[1A>T;2T>C]".to_string(),
        ),
        (
            "MITO:m.[1G>T;2A>C]".to_string(),
            "MITO:m.[1A>T;2T>C]".to_string(),
        ),
    ] {
        assert!(
            strict_rejects(&input),
            "`{input}` misstates its reference; the axis it is written on must \
             not decide whether that is caught (#1068, #1543)"
        );
        assert!(
            !strict_rejects(&corrected),
            "`{corrected}` states its reference correctly and must be accepted \
             on every axis (#1543)"
        );
    }
}

/// The `r.` axis is judged in its own alphabet, so a truthful `u` over a `T` is
/// not a mismatch.
///
/// `r.` spells its bases as `a/c/g/u` while the transcript is served as DNA, so
/// a validator comparing the two byte for byte would reject `r.2u>a` over a `T`
/// — over-rejecting *correct* input, which is worse than the hole being closed.
/// `normalize_rna` rewrites `u` to `T` before it validates (#736); the
/// authored-member guard makes the same rewrite, and this pins it.
#[test]
fn the_rna_alphabet_is_not_read_as_a_mismatch() {
    // `c.2` is `T`, i.e. `r.2` is `u`.
    for input in [format!("{TX}:r.2u>a"), format!("{TX}:r.[2u>a;4c>g]")] {
        assert!(
            !strict_rejects(&input),
            "`{input}` states the reference truthfully in the RNA alphabet; \
             folding u/T is what keeps the #1543 guard from over-rejecting"
        );
        assert!(
            reference_mismatches(&input).is_empty(),
            "`{input}` must carry no RefSeqMismatch (#736, #1543)"
        );
    }
}
