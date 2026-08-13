//! The coding codon exception, exercised where the block and the window do
//! **not** share an origin.
//!
//! # The coordinate basis this guards
//!
//! `canonicalize_from_sequence` fetches a reference *window* `[w_lo, w_hi]`
//! padded by `CANONICAL_PAD` around the members, applies the edits, and then
//! trims to the minimal changed *block* with `trim_common_flanks`. That trim
//! yields an offset `lo`, and from it two different coordinate bases coexist in
//! the same function:
//!
//! * **block-relative** — offsets into `ref_bytes[lo..hi_ref]`, whose axis
//!   coordinate is `w_lo + lo + offset`;
//! * **window-relative** — offsets into `ref_bytes`, whose axis coordinate is
//!   `w_lo + offset`.
//!
//! Pieces start block-relative and are shifted into the window (`ref_start +=
//! lo`) before `apply_coding_codon_exception` and
//! `coalesce_coding_frame_separation` run, so those two passes must convert with
//! `w_lo` and nothing else. Handing either of them `w_lo + lo` — the *block*
//! basis — reads the codon arithmetic `lo` positions too far along. That is a
//! real defect that has been written: it was found in a partition-preserving
//! arm whose `split_codon_incompatible_triplets` call passed block coordinates
//! for pieces that had already been offset into the window, an off-by-`lo` on
//! both the bases read and the codon phase.
//!
//! # Why every pre-existing codon fixture is blind to it
//!
//! `lo` is `min(CANONICAL_PAD, c_lo - 1)`: the window starts at `c_lo -
//! CANONICAL_PAD`, clamped to 1. Every codon-exception fixture in this suite
//! sits within a hundred nucleotides of `c.1`, so `w_lo` clamps to `1` and `lo`
//! collapses to `c_lo - 1`.
//!
//! That regime is not merely a weak test, it is **structurally** blind. The
//! codon exception's positive case needs a codon-aligned first position, i.e.
//! `c_lo % 3 == 1`; and in the clamped regime `lo == c_lo - 1`, which is then
//! always a multiple of 3. A phase error of `lo` is a phase error of zero, so
//! the block basis and the window basis agree by construction and a test there
//! cannot tell a correct implementation from the broken one — whatever
//! assertion it makes.
//!
//! Escaping it requires a variant far enough into the transcript that `w_lo` is
//! set by the pad rather than by the clamp. Then `lo == CANONICAL_PAD == 128`,
//! and `128 % 3 == 2`, so the two bases disagree by two codon phases and the
//! codon verdict genuinely differs between them.
//!
//! **This guard discriminates only while `CANONICAL_PAD % 3 != 0`.** It is 128
//! today. Were the pad ever set to a multiple of 3, every assertion here would
//! still pass under the broken basis, exactly as the `c.1` fixtures do — so
//! change the pad and this file needs re-deriving, not re-running. Within that
//! constraint the guard is insensitive to the pad's exact value: the merging
//! case sits at `c_lo % 3 == 1`, and shifting it by either non-zero residue
//! lands it across a codon boundary, so *any* `lo % 3 != 0` is caught.
//!
//! # Which assertions here discriminate, and which are controls
//!
//! Measured, not asserted: by reintroducing the defect at both window-relative
//! call sites (`w_lo` -> `w_lo + lo`) and re-running, and separately by making
//! `apply_coding_codon_exception` a no-op to find what reaches it at all. The
//! split is written down because a file in which every assertion *reads* like a
//! guard, while most of them cannot fail, is the exact failure this fixture
//! exists to end.
//!
//! Two tests **discriminate**, one for each window-relative pass, and they fail
//! in opposite directions:
//!
//! | test | pass | wrong-basis output |
//! |---|---|---|
//! | `a_length_changing_straddling_pair_deep_in_the_transcript_stays_split` | `coalesce_coding_frame_separation` | `c.449_451delinsTT` — a merge across a codon boundary |
//! | `a_literal_delins_holding_an_in_codon_triplet_stays_whole` | `apply_coding_codon_exception` | `c.[301G>A;303_304delinsCT]` — a split inside one codon |
//!
//! "Forbids" would overstate the authority: `general.md:34` and `:35` are
//! lowercase prose graded by a **modal**, and the ruling record
//! `separation-rule-force-modal-or-negation` settles that the clause states a
//! preference rather than a prohibition. Each wrong-basis output above is a
//! deviation this file pins as a tripwire, not a rule the spec bans.
//!
//! The rest are **controls**, and were measured not to move under that same
//! sabotage. Two reasons, worth keeping apart:
//!
//! * The cis-pair tests (`an_in_codon_pair_*`, `a_codon_straddling_pair_*`,
//!   `a_length_changing_in_codon_pair_*`) are merged by
//!   `merge_consecutive_edits`, which works in axis coordinates from end to end
//!   and never converts a window offset — so no coordinate-basis error *can*
//!   move them. Verified separately by making `apply_coding_codon_exception` a
//!   no-op, which leaves all of them green while reddening three tests in
//!   `issue_165_delins_sub_only_decompose`.
//! * `a_literal_delins_straddling_a_codon_boundary_is_split` reaches the right
//!   pass and its codon verdict does invert under the wrong basis, but the
//!   output does not: something downstream declines the merge anyway. It is kept
//!   as the negative half of the literal-`delins` pair, not offered as evidence
//!   about the basis.
//!
//! They hold the fixture's amino-acid claims and pin that the deep positions
//! behave like the shallow ones. Do not read them as basis guards.
//!
//! # Spec
//!
//! > "two variants separated by one or more nucleotides should be described
//! > individually and **not** as a 'delins'."
//! >   — `assets/hgvs-nomenclature/docs/recommendations/general.md:34`
//! >
//! > "**exception**: two variants separated by one nucleotide, together
//! > affecting one amino acid, should be described as a 'delins'."
//! >   — `general.md:35`

use ferro_hgvs::backtranslate::codon::{Codon, CodonTable};
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Transcript length. Long enough that a site in the middle is more than
/// `CANONICAL_PAD` (128) from either end, so the canonical window is set by the
/// pad on both sides and `lo` is the pad rather than a clamp.
const TRANSCRIPT_LEN: usize = 600;

/// The canonicalizer's window pad, hand-copied from `merge.rs`'s private
/// `CANONICAL_PAD` (`merge.rs:1822`).
///
/// It is a **mirror, not a reference**: the crate's constant is module-private,
/// so nothing here can read it and nothing here can verify this copy against
/// it. The value is load-bearing rather than decorative — it is read by the
/// site-depth and codon-phase assertions in
/// [`the_deep_sites_are_far_enough_in_to_escape_the_window_clamp`] as well as
/// by the assertion messages — which is exactly why the divergence matters. See
/// that test's doc for what a change to the crate's pad would and would not do
/// here.
const CANONICAL_PAD: i64 = 128;

/// A 600 nt single-exon coding transcript with `cds_start = 1`, so `c.p` is
/// transcript position `p` and codon `k` spans `c.3k-2 ..= c.3k`.
///
/// The body is `ATG` followed by a `GCT` tandem, which has two properties the
/// fixture needs and which a hand-typed sequence would not reliably give:
///
/// * **no homopolymer anywhere**, so no member can 3'-shift and two sites can
///   differ in codon phase and in nothing else;
/// * **a uniform codon**, `GCT` (Ala) for every codon from 2 on, so the
///   reference amino acid is the same at each site and the residue counts below
///   are attributable to the variant rather than to the local sequence.
///
/// Generated rather than committed, per the repository's test-data rule.
fn deep_coding_transcript() -> MockProvider {
    let sequence = format!("ATG{}", "GCT".repeat((TRANSCRIPT_LEN - 3) / 3));
    assert_eq!(
        sequence.len(),
        TRANSCRIPT_LEN,
        "the codon arithmetic in this module assumes a {TRANSCRIPT_LEN} nt transcript",
    );
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NM_DEEP.1".to_string(),
        Some("DEEP".to_string()),
        Strand::Plus,
        sequence,
        Some(1),
        Some(len),
        vec![Exon::new(1, 1, len)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    let mut provider = MockProvider::new();
    provider.add_transcript(transcript);
    provider
}

/// The transcript sequence the provider above serves, as bytes, 0-based.
fn deep_coding_sequence() -> Vec<u8> {
    format!("ATG{}", "GCT".repeat((TRANSCRIPT_LEN - 3) / 3)).into_bytes()
}

/// The reference base at `c.position` on the fixture (`cds_start = 1`).
fn base_at(position: usize) -> char {
    deep_coding_sequence()[position - 1] as char
}

fn normalize(input: &str) -> String {
    let normalizer = Normalizer::new(deep_coding_transcript());
    let variant = parse_hgvs(input).expect("parse failed");
    normalizer
        .normalize(&variant)
        .expect("normalize failed")
        .to_string()
}

/// The 1-based codon number holding CDS position `p` — the numbering
/// `merge::same_codon` is defined on, restated here so the fixture's claims are
/// checked against an independent derivation rather than against the crate's.
fn codon_of(position: i64) -> i64 {
    (position - 1) / 3 + 1
}

/// How many residues the full-CDS translation of `mutant` differs from that of
/// `reference` in.
///
/// The spec clause the exception rests on is about amino acids, so the fixture
/// states its precondition in amino acids: an in-codon pair changes **one**
/// residue and a boundary-straddling pair changes **two**. Both sequences must
/// be the same length — this is only meaningful for equal-length (substitution)
/// pairs, since a length-changing pair reframes everything downstream.
fn residues_changed(reference: &[u8], mutant: &[u8]) -> usize {
    assert_eq!(
        reference.len(),
        mutant.len(),
        "residue counting is only defined for an equal-length pair",
    );
    let table = CodonTable::standard();
    let translate = |cds: &[u8]| -> Vec<Option<_>> {
        cds.chunks_exact(3)
            .map(|triplet| {
                let codon = Codon::parse_bytes(triplet).expect("fixture CDS must be ACGT");
                table.amino_acid_for(&codon)
            })
            .collect()
    };
    let (before, after) = (translate(reference), translate(mutant));
    before
        .iter()
        .zip(after.iter())
        .filter(|(a, b)| a != b)
        .count()
}

/// Apply substitutions (1-based CDS position, alternate base) to the fixture.
fn with_substitutions(substitutions: &[(usize, u8)]) -> Vec<u8> {
    let mut sequence = deep_coding_sequence();
    for &(position, alternate) in substitutions {
        assert_ne!(
            sequence[position - 1],
            alternate,
            "c.{position} must actually change, else the fixture asserts nothing",
        );
        sequence[position - 1] = alternate;
    }
    sequence
}

/// The geometry claim this whole file rests on, asserted rather than assumed.
///
/// A site is discriminating only if the window's origin and the block's origin
/// differ by an amount that is **not** a multiple of 3. In the pad-limited
/// regime that difference is `CANONICAL_PAD`.
///
/// **What these assertions can and cannot catch.** They pin that the four sites
/// are deep enough for the pad-limited regime, and that the pad *this file
/// assumes* is not a multiple of 3. Both read the local mirror, which is the
/// only value a test outside the crate can see, so a change to `merge.rs`'s own
/// `CANONICAL_PAD` is **invisible here**: the mirror would still say 128 and
/// every assertion in this file would still pass.
///
/// That gap cannot be closed from a test. A pad that became a multiple of 3
/// would leave every *output* in this file unchanged — that is precisely the
/// blindness the module doc describes — so there is no observable for an
/// assertion to key on. The guard against a pad change is the module doc's
/// instruction to re-derive this file, not anything asserted below. Do not read
/// the `assert_ne!` as one; per this file's own standard, an assertion that
/// cannot fail must not be presented as one that can.
#[test]
fn the_deep_sites_are_far_enough_in_to_escape_the_window_clamp() {
    for site in [301, 353, 397, 449] {
        assert!(
            site as i64 - CANONICAL_PAD >= 1,
            "c.{site} must be more than {CANONICAL_PAD} nt from c.1, or the window \
             clamps to 1 and `lo` collapses to `c_lo - 1`",
        );
        assert!(
            site as i64 + 2 + CANONICAL_PAD <= TRANSCRIPT_LEN as i64,
            "c.{site} must be more than {CANONICAL_PAD} nt from the transcript end, \
             so the window is set by the pad on both sides",
        );
    }
    assert_ne!(
        CANONICAL_PAD % 3,
        0,
        "the block/window origins differ by `lo` == {CANONICAL_PAD}; if that is a \
         multiple of 3 the two coordinate bases agree in codon phase and every \
         assertion in this file passes under the broken basis too",
    );
}

/// **The exception applies, at a non-zero window offset.** `c.301` and `c.303`
/// are one unchanged nucleotide apart and both lie in codon 101, so they
/// together affect one amino acid and `general.md:35` asks for a `delins`.
///
/// This is the discriminating assertion of the file. `c.301 % 3 == 1`, so the
/// pair is codon-aligned in the window basis the pass must use; add `lo` (128,
/// i.e. two codon phases) and `c.429_431` straddles codons 143 and 144, the
/// exception is refused, and the output falls back to the two individually
/// described substitutions. Same reference triplet, same payloads, same axis —
/// only the coordinate basis differs.
#[test]
fn an_in_codon_pair_deep_in_the_transcript_merges() {
    assert_eq!(
        (base_at(301), base_at(302), base_at(303)),
        ('G', 'C', 'T'),
        "the fixture's codon 101 must be the uniform `GCT`",
    );
    assert_eq!(
        codon_of(301),
        codon_of(303),
        "c.301 and c.303 must share a codon, else this tests the wrong clause",
    );

    // The amino-acid half of `general.md:35`, stated in amino acids: the pair
    // changes exactly one residue, Ala101 -> Thr101.
    assert_eq!(
        residues_changed(
            &deep_coding_sequence(),
            &with_substitutions(&[(301, b'A'), (303, b'C')]),
        ),
        1,
        "`together affecting one amino acid` is the exception's second conjunct",
    );

    assert_eq!(
        normalize("NM_DEEP.1:c.[301G>A;303T>C]"),
        "NM_DEEP.1:c.301_303delinsACC",
        "codon 101 holds both changes, so `general.md:35` asks for a delins — a \
         split here means the codon arithmetic ran in the block basis, {} codon \
         phases off",
        CANONICAL_PAD % 3,
    );
}

/// **The exception does not apply, at the same depth.** `c.353` and `c.355` are
/// one unchanged nucleotide apart — the distance conjunct of `general.md:35` is
/// met, identically to the merging case — but they lie in codons 118 and 119,
/// so the amino-acid conjunct is not, and `general.md:34` governs.
///
/// The mirror of the test above: `c.353 % 3 == 2`, so adding `lo` (two codon
/// phases) lands `c.481` codon-aligned and merges a pair the spec requires be
/// described individually. The two tests fail in opposite directions under the
/// same wrong basis.
///
/// The alternate bases are chosen so that **both** codons change residue.
/// `c.345`/`c.347` was the other straddling phase available on this fixture and
/// was dropped: `c.345` is a third codon position, where every substitution in
/// the uniform `GCT` background is silent (`GCN` is Ala throughout), so the pair
/// changes one residue and cannot state the two-amino-acid precondition.
#[test]
fn a_codon_straddling_pair_deep_in_the_transcript_stays_split() {
    assert_ne!(
        codon_of(353),
        codon_of(355),
        "c.353 and c.355 must NOT share a codon, else this tests the wrong clause",
    );
    assert_eq!(
        (base_at(353), base_at(355)),
        ('C', 'G'),
        "the fixture's bases at the two changed positions",
    );

    // Codon 118 `GCT` (Ala) -> `GAT` (Asp); codon 119 `GCT` (Ala) -> `TCT`
    // (Ser). Two amino acids, which is precisely what puts the pair outside
    // `general.md:35`.
    assert_eq!(
        residues_changed(
            &deep_coding_sequence(),
            &with_substitutions(&[(353, b'A'), (355, b'T')]),
        ),
        2,
        "a pair spanning two codons affects two amino acids, which is what \
         puts it outside `general.md:35`",
    );

    assert_eq!(
        normalize("NM_DEEP.1:c.[353C>A;355G>T]"),
        "NM_DEEP.1:c.[353C>A;355G>T]",
        "two amino acids are affected, so `general.md:34` governs — a delins \
         here means the codon arithmetic ran in the block basis",
    );
}

/// The length-changing route to the same clause, which reaches
/// `coalesce_coding_frame_separation` rather than `apply_coding_codon_exception`
/// — a second pass with a second `w_lo` argument, and so a second place the two
/// bases can be confused.
///
/// `c.397del` and `c.399T>A` are separated by the unchanged `c.398`, and
/// `c.397_399` is codon 133 exactly. The clause names no edit type, so the pair
/// is a `delins`; add `lo` and the span reads as `c.525_527`, which straddles
/// codons 175 and 176 and is refused.
#[test]
fn a_length_changing_in_codon_pair_deep_in_the_transcript_merges() {
    assert_eq!(
        codon_of(397),
        codon_of(399),
        "c.397 and c.399 must share a codon",
    );
    assert_eq!(
        normalize("NM_DEEP.1:c.[397del;399T>A]"),
        "NM_DEEP.1:c.397_399delinsCA",
        "the deleted `G`, the carried-through `C` at c.398 and the substituted \
         `A` are one codon's worth of change",
    );
}

/// The other side of the length-changing pair. `c.449` is in codon 150 and
/// `c.451` in codon 151, so `general.md:34` governs; under the block basis
/// `c.577_579` is codon-aligned and the two would merge.
#[test]
fn a_length_changing_straddling_pair_deep_in_the_transcript_stays_split() {
    assert_ne!(
        codon_of(449),
        codon_of(451),
        "c.449 and c.451 must NOT share a codon",
    );
    assert_eq!(
        normalize("NM_DEEP.1:c.[449del;451G>T]"),
        "NM_DEEP.1:c.[449del;451G>T]",
        "two amino acids are affected, so the two variants are described \
         individually",
    );
}

/// The same clause reached through the **other** window-relative pass.
///
/// A cis pair of members is re-merged by `merge_consecutive_edits`, which works
/// in axis coordinates throughout and so cannot express this defect — that is
/// what serves the two substitution pairs above, and it is why they are controls
/// rather than discriminators (see the module doc). A **literal `delins`** is
/// different: it arrives as one member, is decomposed against the reference, and
/// is put back together by `apply_coding_codon_exception`, which converts window
/// offsets with `w_lo` and is exactly where the basis can be wrong.
///
/// `c.301_304delinsACCT` decomposes to `[Sub@301; Identity@302; Sub@303;
/// Sub@304]`. The `(301, 303)` pair is codon 101, so `general.md:35` merges it,
/// and `c.304` joins because it touches `c.303` (`delins.md:16` — two
/// consecutive changed nucleotides are one `delins`). The whole thing is
/// therefore a fixed point. Read the pair's codon in the block basis instead and
/// it is `c.429`/`c.431`, codons 143 and 144, the exception is refused, and the
/// member is emitted split.
///
/// The shallow twin is `issue_165_delins_sub_only_decompose::
/// cds_embedded_codon_frame_triplet_from_literal_delins`, at `c.10_13` where
/// `lo` is 9 and the two bases agree.
#[test]
fn a_literal_delins_holding_an_in_codon_triplet_stays_whole() {
    assert_eq!(
        (base_at(301), base_at(302), base_at(303), base_at(304)),
        ('G', 'C', 'T', 'G'),
    );
    assert_eq!(codon_of(301), codon_of(303), "the merged pair is one codon");

    assert_eq!(
        normalize("NM_DEEP.1:c.301_304delinsACCT"),
        "NM_DEEP.1:c.301_304delinsACCT",
        "`general.md:35` merges c.301 with c.303 and `delins.md:16` pulls c.304 \
         in — a split here means the codon arithmetic ran in the block basis",
    );
}

/// The mirror of the test above, at the codon phase where the *wrong* basis
/// says "one amino acid".
///
/// `c.353_356delinsATTA` decomposes to `[Sub@353; Identity@354; Sub@355;
/// Sub@356]`. `c.353` and `c.355` are in codons 118 and 119, so `general.md:35`
/// does not reach them and `general.md:34` keeps `c.353` individual; `c.355` and
/// `c.356` are consecutive and are one `delins` by `delins.md:16`.
///
/// **A control, not a basis guard.** In the block basis the pair reads as
/// `c.481`/`c.483`, codon 161 both, so `apply_coding_codon_exception` *would*
/// merge it — but the output was measured unchanged under that sabotage, so
/// something later declines the merge regardless and this assertion cannot see
/// the defect. Its sibling above can, in the other direction. Recorded rather
/// than deleted because the expectation is correct and the pair reads better
/// together; recorded rather than left implied because an assertion that cannot
/// fail must not be presented as one that can.
#[test]
fn a_literal_delins_straddling_a_codon_boundary_is_split() {
    assert_eq!(
        (base_at(353), base_at(354), base_at(355), base_at(356)),
        ('C', 'T', 'G', 'C'),
    );
    assert_ne!(
        codon_of(353),
        codon_of(355),
        "the pair the exception would merge must span two codons",
    );

    assert_eq!(
        normalize("NM_DEEP.1:c.353_356delinsATTA"),
        "NM_DEEP.1:c.[353C>A;355_356delinsTA]",
        "two amino acids are affected, so `general.md:34` keeps c.353 \
         individual — one whole member here means the codon arithmetic ran in \
         the block basis",
    );
}

/// The axis control for the assertions above.
///
/// Without it, the splits could be produced by any axis-blind decline — a window
/// clamp, the weight bound, the round-trip guard — rather than by the codon
/// test. `NM_DEEP.1` has `cds_start = 1`, so `n.p` and `c.p` address the same
/// base and each pair below differs from its `c.` sibling in nothing but the
/// axis letter. On `n.` there is no amino acid to affect, so `general.md:35`
/// does not apply at either phase and *neither* position merges — including
/// `c.301`, the one position where the `c.` axis does merge.
#[test]
fn without_a_reading_frame_no_deep_position_merges() {
    for input in [
        "NM_DEEP.1:n.[301G>A;303T>C]",
        "NM_DEEP.1:n.[353C>A;355G>T]",
        "NM_DEEP.1:n.[397del;399T>A]",
        "NM_DEEP.1:n.[449del;451G>T]",
        "NM_DEEP.1:n.[353C>A;355_356delinsTA]",
    ] {
        assert_eq!(
            normalize(input),
            input,
            "`general.md:34`'s plain rule governs an axis with no amino acid to affect",
        );
    }

    // The literal `delins` is the sharpest of the five, because it is the one
    // input whose `c.` and `n.` answers differ in *shape* rather than merely in
    // whether a merge happened: with no codon to hold them together, `c.301` is
    // described individually and only the two consecutive positions stay one
    // member.
    assert_eq!(
        normalize("NM_DEEP.1:n.301_304delinsACCT"),
        "NM_DEEP.1:n.[301G>A;303_304delinsCT]",
        "with no reading frame the c.301/c.303 pair has no amino acid to \
         share, so `general.md:34` splits what the `c.` axis keeps whole",
    );
}
