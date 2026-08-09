//! #1453 — a cis allele on the `r.` axis of a **non-coding** transcript emitted
//! the same member twice.
//!
//! ```text
//! NR_TEST.1:r.[9dup;10dup;11dup]  ->  NR_TEST.1:r.[9dup;11dup;11dup]
//! NR_TEST.1:n.[9dup;10dup;11dup]  ->  NR_TEST.1:n.9_10insTAA
//! ```
//!
//! Two members claiming one interbase point denote no sequence at all, so this
//! is a well-formedness defect rather than a representation preference — and it
//! is invisible to two of the three seam oracles. `FERRO_ASSERT_REPARSE` accepts
//! `r.[9dup;11dup;11dup]` because it is syntactically well-formed, and
//! `FERRO_ASSERT_IN_BOUNDS` accepts it because every coordinate is in range.
//!
//! ## The mechanism, which is narrower than "the `r.` axis"
//!
//! The issue framed this as `r.` versus `n.`. It is not: `r.` on a **coding**
//! transcript collapses the same three members correctly
//! (`NM_TEST.1:r.[9dup;10dup;11dup]` -> `NM_TEST.1:r.9_10insuaa`). The
//! discriminator is whether the transcript has a CDS.
//!
//! `coalesce_members_at_one_junction` (#1286) is the pass that merges two
//! members which shifted onto one junction, and it does reach this input — it
//! groups the two `11dup`s and reads both payloads. It then calls
//! `respell_at_gap`, whose `r.` arm resolves the gap *through the transcript*
//! via `ordered_cds_bounds` (#1284), because on a coding transcript the `r.`
//! axis is CDS-relative and `+1` across the CDS start would name the
//! non-existent `r.0`. `ordered_cds_bounds` requires **both** `cds_start` and
//! `cds_end`, so on a non-coding record it returns `None`, `respell_at_gap`
//! reverts the member, and the pass abandons the group — leaving exactly the
//! repeated member it exists to remove.
//!
//! The fallback is the one `axis_frame` and `region_sequence_delta` already
//! make for `Region::Rna`: with no CDS the `r.` axis numbers the transcript
//! directly, exactly as `n.` does, so the gap is `[junction, junction + 1]` on
//! the axis itself.
//!
//! ## Two exits, not one
//!
//! `respell_at_gap` resolves through the CDS in two places, and both had to be
//! fixed. The second is the **terminal** path: when the landing junction is one
//! past the last base, the repair is re-spelled as the equivalent boundary
//! `delins` (#1327), and `respell_at_sequence_end` resolves the last base
//! through `cds_axis_end` — `ordered_cds_bounds` again. It calls it *twice*,
//! once to place the write and once to verify it landed, so a fallback applied
//! to only one of the two leaves the guard rejecting a correct write. The first
//! cut of this fix did exactly that, and the census below still showed one
//! surviving row until `rna_axis_end` was shared between them.
//!
//! ## Census
//!
//! Over 3,200 rows — 6 deterministic homopolymer-enriched 63 nt cores plus the
//! issue's, × two- and three-member cis alleles of nearby junction-occupying
//! edits, × both shuffle directions, on the non-coding `r.` axis:
//!
//! | | on `main` | with the fix |
//! |---|---|---|
//! | outputs repeating a member | 849 (2 families) | **0** |
//! | outputs disagreeing with the `n.` axis of the same record | 1,046 (3 families) | **0** |
//!
//! ## What these tests pin
//!
//! The `n.`/`r.` pair is asserted as **agreement**, not as two independent
//! expected strings. A test that pinned only the `r.` string would go green if
//! a future change moved the `n.` axis onto the same defect, which is the
//! opposite of what this issue is about.

use ferro_hgvs::conformance::spec_corpus::{denotation_of, Denotation};
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// The issue's 63 nt core. Positions 1-9 are `T`, 10 and 11 are `A`, 12 is `T`,
/// so `9dup`, `10dup` and `11dup` each land inside a homopolymer run and shift
/// independently — `10dup` and `11dup` both to the junction 3' of position 11.
const CORE: &str = "TTTTTTTTTAATATATTTTAATATAATTAAAAAAATAATTTTTATAAATATATTATTTTAAAAA";

/// A core whose 3'-most homopolymer runs to the **last** base (positions 58-63
/// are `G` on this 63 nt sequence), so a member shifting into that run lands its
/// copy at the junction one past the end.
///
/// That is the terminal path: `respell_at_gap` hands off to
/// `respell_at_sequence_end`, which resolves the last base through
/// `cds_axis_end` rather than through `cds_relative_gap`. Drawn from the census
/// rather than authored — it is the one row that survived the first cut of the
/// fix.
const TERMINAL_CORE: &str = "CGGGTGGGGTTTTCCCCCCGAAAAAAAAAAAAAAAGAGGGGGTTTTTAAAATTTAAAGGGGGG";

/// One transcript carrying [`CORE`], coding or not.
///
/// The coding arm sets the CDS to the whole transcript (`1..=63`) so that `c.N`,
/// `n.N` and `r.N` all name transcript position `N`. That is what makes the
/// coding/non-coding pair below a controlled comparison: the two records differ
/// only in whether a CDS is annotated, never in which bases a position names.
fn provider(accession: &str, coding: bool) -> MockProvider {
    provider_for(accession, coding, CORE)
}

fn provider_for(accession: &str, coding: bool, core: &str) -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence = core.to_string();
    let length = sequence.len() as u64;
    let (cds_start, cds_end) = if coding {
        (Some(1), Some(length))
    } else {
        (None, None)
    };
    provider.add_transcript(Transcript::new(
        accession.to_string(),
        Some("TESTGENE".to_string()),
        Strand::Plus,
        sequence,
        cds_start,
        cds_end,
        vec![Exon::new(1, 1, length)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

fn normalize(provider: &MockProvider, descriptor: &str) -> String {
    let variant = parse_hgvs(descriptor).expect("fixture must parse");
    Normalizer::new(provider.clone())
        .normalize(&variant)
        .expect("fixture must normalize")
        .to_string()
}

/// Normalize `input`, assert it reaches `expected`, and assert the two
/// properties every test in this file *claimed* to have measured: that the
/// output denotes the same bases as the input, and that it is a fixed point.
///
/// This file has no apply oracle of its own, so the denotation goes through
/// `spec_corpus::denotation_of` — i.e. `hgvs_to_spdi`, a path the normalizer
/// does not consult. Five doc comments here stated `Denotation::Sequence`
/// equality as a measurement while every test asserted an exact string and
/// nothing else; a pinned string cannot separate a re-spelling from a
/// corruption, which for this file is the whole question — a repeated member
/// denotes **no** sequence, and that is the defect #1453 is about.
///
/// `Denotation::Sequence` is asserted explicitly rather than left implied,
/// because two `Inexpressible` or two `NoSequence` values compare equal and
/// would make the check pass while measuring nothing.
fn assert_reaches_preserving(
    provider: &MockProvider,
    served: &str,
    input: &str,
    expected: &str,
    message: &str,
) {
    let output = normalize(provider, input);
    assert_eq!(output, expected, "{message}");

    let from_input = denotation_of(provider, served, input);
    assert!(
        matches!(from_input, Denotation::Sequence(_)),
        "`{input}` denotes no sequence ({from_input:?}), so the comparison below \
         would measure nothing"
    );
    assert_eq!(
        denotation_of(provider, served, &output),
        from_input,
        "`{input}` -> `{output}` changed the sequence"
    );

    let again = normalize(provider, &output);
    assert_eq!(
        again, output,
        "`{input}` -> `{output}` is not a fixed point: a second pass reaches `{again}`"
    );
}

/// The `n.` rendering of a description, re-spelled in the `r.` alphabet: `r.`
/// lower-cases its sequences and writes `u` for `T`.
///
/// Only the axis letter and the inserted bases differ between the two axes on a
/// non-coding transcript, which is the whole point of the parity assertion —
/// the coordinates must be identical because the two axes number the same
/// transcript.
fn as_rna_spelling(tx_rendering: &str) -> String {
    let (accession, description) = tx_rendering
        .split_once(":n.")
        .expect("an `n.` rendering to convert");
    let description: String = description
        .chars()
        .map(|c| match c {
            'T' => 'u',
            other if other.is_ascii_uppercase() => other.to_ascii_lowercase(),
            other => other,
        })
        .collect();
    format!("{accession}:r.{description}")
}

/// The reproducer: three duplications inside two homopolymer runs must not
/// produce a member twice, on the `r.` axis of a non-coding transcript.
///
/// Pinned as an exact string rather than as "no repeated member": a predicate
/// over the member set is satisfied by a refusal or by a dropped member just as
/// well as by the correct collapse.
///
/// RE-BLESSED under `partition-is-the-unit-of-normalization` (DECIDED,
/// 2026-08-08, `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`).
///
/// ```text
/// was  NR_TEST.1:r.9_10insuaa            one member, re-derived
/// now  NR_TEST.1:r.[9dup;10_11a[4]]      two members
/// ```
///
/// **What the ruling changed here is only the first member.** `9dup` sits in
/// the `u` run at 1-9 and writes at the junction `9|10`; the other two write at
/// `11|12`, two unchanged bases away, which `general.md:34` says to describe
/// individually. The old single insertion fused all three by re-deriving the
/// partition from the resulting sequence.
///
/// **What the ruling did NOT change is the collapse this file exists for.**
/// `10dup` and `11dup` both shift onto the junction 3' of `r.11` — separation
/// zero, below the axis floor — so the ruling's MERGE move applies and they
/// become the one member `10_11a[4]`. The #1453 defect was emitting them as
/// `r.[11dup;11dup]`, two members claiming one interbase point; that is still
/// gone, and the assertion is still an exact string rather than a predicate.
///
/// **Sequence unchanged, and this file has no apply oracle, so it was measured
/// rather than assumed.** `spec_corpus::denotation_of` (via `hgvs_to_spdi`, not
/// the normalizer) gives `Denotation::Sequence` for both input and output and
/// the two are equal — `u`×10 then `a`×4 then `u`… — where a repeated member
/// would have given `Denotation::NoSequence`. The output is a fixed point.
#[test]
fn a_noncoding_rna_allele_collapses_instead_of_repeating_a_member() {
    let provider = provider("NR_TEST.1", false);
    assert_reaches_preserving(
        &provider,
        CORE,
        "NR_TEST.1:r.[9dup;10dup;11dup]",
        "NR_TEST.1:r.[9dup;10_11a[4]]",
        "the two members shifted onto the junction 3' of `r.11` must be \
         coalesced into one; emitting `r.[9dup;11dup;11dup]` claims one \
         interbase point twice and so denotes no sequence (#1453)",
    );
}

/// The same three members on the `n.` axis of the *same* transcript, asserted as
/// agreement with the `r.` axis rather than as an independent expectation.
#[test]
fn the_noncoding_rna_and_tx_axes_agree_on_one_transcript() {
    let provider = provider("NR_TEST.1", false);
    let tx = normalize(&provider, "NR_TEST.1:n.[9dup;10dup;11dup]");
    let rna = normalize(&provider, "NR_TEST.1:r.[9dup;10dup;11dup]");
    assert_eq!(
        rna,
        as_rna_spelling(&tx),
        "`n.` and `r.` number the same non-coding transcript, so the same three \
         members must reach the same coordinates on both — only the rendering \
         alphabet may differ (`n.` gave `{tx}`)"
    );
}

/// The discriminating case, and the one that says the fix is scoped to the
/// missing non-coding fallback rather than to the `r.` axis as a whole: a
/// **coding** transcript's `r.` axis already collapsed this correctly, and must
/// keep doing so.
///
/// Its `r.` gap is still resolved through the CDS (`respell_at_gap`'s
/// `cds_relative_gap`), which is what names a junction on the CDS start
/// `r.-1_1` instead of the non-existent `r.0` (#1284).
///
/// **This fixture cannot detect an over-broad fix on its own.** Its CDS spans
/// the whole transcript, so `cds_start = 1`, the CDS delta is 0, and the
/// CDS-relative and transcript-relative resolutions coincide — measured:
/// replacing the predicate with a bare `span.region == Region::Rna`, applying
/// the non-coding fallback to *every* `r.` record, leaves this test and the
/// terminal one below passing. That is the price of the controlled comparison
/// (`c.N`, `n.N` and `r.N` naming one base) and it is worth paying here, but
/// the scoping claim has to be pinned somewhere the two conventions differ:
/// [`a_coding_transcripts_terminal_rna_repair_resolves_through_the_cds`].
///
/// RE-BLESSED under `partition-is-the-unit-of-normalization` (DECIDED,
/// 2026-08-08, `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`).
///
/// ```text
/// r.  was NM_TEST.1:r.9_10insuaa    now NM_TEST.1:r.[9dup;10_11dup]
/// c.  was NM_TEST.1:c.9_10insTAA    now NM_TEST.1:c.[9dup;10_11dup]
/// ```
///
/// Same movement as the non-coding reproducer above and for the same reason:
/// `9dup` is two unchanged bases from the pair, so `general.md:34` describes it
/// individually, while `10dup` and `11dup` land on one junction and merge. The
/// coding record spells that merged member `10_11dup` rather than the
/// non-coding record's `10_11a[4]`, which is the per-axis repeat-vs-dup
/// preference and not a difference in the partition — both denote the same two
/// added `a`s.
///
/// **The discriminating claim this row carries is untouched.** It is that the
/// coding `r.` axis is decided by the CDS-relative gap resolution rather than
/// by #1453's non-coding fallback, and that is a statement about *which*
/// resolution runs, not about how many members come out. The `r.`/`c.` pair
/// still agrees, which is what says the axes number the same transcript.
///
/// Sequence unchanged, measured with `spec_corpus::denotation_of` (via
/// `hgvs_to_spdi`, not the normalizer) on both axes; both outputs are fixed
/// points.
#[test]
fn a_coding_transcripts_rna_axis_is_unchanged() {
    let provider = provider("NM_TEST.1", true);
    assert_reaches_preserving(
        &provider,
        CORE,
        "NM_TEST.1:r.[9dup;10dup;11dup]",
        "NM_TEST.1:r.[9dup;10_11dup]",
        "the coding `r.` axis must still coalesce the two members that land on \
         one junction — it is the CDS-relative gap resolution that #1453 leaves \
         in place",
    );
    assert_reaches_preserving(
        &provider,
        CORE,
        "NM_TEST.1:c.[9dup;10dup;11dup]",
        "NM_TEST.1:c.[9dup;10_11dup]",
        "and so must the `c.` axis of the same record",
    );
}

/// The terminal path (`respell_at_sequence_end`), where the landing junction is
/// one past the last base and the repair becomes the boundary `delins` of
/// #1327.
///
/// This is a *second* exit from `respell_at_gap`, not the one the reproducer
/// above takes, and it declines through `cds_axis_end` rather than
/// `cds_relative_gap`. It is pinned separately because a fix applied to only
/// the junction path leaves it broken — which is how the first cut of this fix
/// measured, with 1 of the census's 849 rows surviving.
#[test]
fn a_noncoding_rna_repair_at_the_sequence_end_also_takes() {
    let provider = provider_for("NR_TEST.1", false, TERMINAL_CORE);
    let tx = normalize(&provider, "NR_TEST.1:n.[57dup;58dup;59dup]");
    let rna = normalize(&provider, "NR_TEST.1:r.[57dup;58dup;59dup]");
    assert_eq!(
        rna, "NR_TEST.1:r.[57dup;58_63g[8]]",
        "the two members that shifted into the terminal `G` run must be \
         coalesced and re-spelled at the boundary; `r.[57dup;63dup;63dup]` \
         claims one interbase point twice (#1453)"
    );
    assert_eq!(
        rna,
        as_rna_spelling(&tx),
        "and must still agree with the `n.` axis of the same record (`n.` gave \
         `{tx}`)"
    );
}

/// The discriminating case for the terminal path, mirroring
/// [`a_coding_transcripts_rna_axis_is_unchanged`]: on a **coding** record the
/// last base is still resolved through the CDS, which is what lets
/// `respell_at_sequence_end` name it `r.*N` when it falls outside the CDS
/// (#1284). Measured on `origin/main` before the fix and unchanged by it.
///
/// RE-BLESSED under `partition-is-the-unit-of-normalization` (DECIDED,
/// 2026-08-08, `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`).
///
/// ```text
/// r.  was NM_TERM.1:r.57_58insagg   now NM_TERM.1:r.[57dup;62_63dup]
/// c.  was NM_TERM.1:c.57_58insAGG   now NM_TERM.1:c.[57dup;62_63dup]
/// ```
///
/// `57dup` and the `58dup`/`59dup` pair are separated by unchanged bases, so
/// `general.md:34` keeps them apart; the pair itself lands on one junction in
/// the terminal `G` run and merges to `62_63dup`. The old single insertion fused
/// all three by re-deriving the partition from the resulting sequence.
///
/// **The claim this row carries — that the coding record still resolves its last
/// base through the CDS — is untouched, and is what the discriminating test
/// below actually measures.** Sequence unchanged on both axes, measured with
/// `spec_corpus::denotation_of` (via `hgvs_to_spdi`); both outputs are fixed
/// points.
#[test]
fn a_coding_transcripts_terminal_rna_repair_is_unchanged() {
    let provider = provider_for("NM_TERM.1", true, TERMINAL_CORE);
    assert_reaches_preserving(
        &provider,
        TERMINAL_CORE,
        "NM_TERM.1:r.[57dup;58dup;59dup]",
        "NM_TERM.1:r.[57dup;62_63dup]",
        "the coding `r.` axis must still resolve its terminal repair through the CDS",
    );
    assert_reaches_preserving(
        &provider,
        TERMINAL_CORE,
        "NM_TERM.1:c.[57dup;58dup;59dup]",
        "NM_TERM.1:c.[57dup;62_63dup]",
        "and so must the `c.` axis of the same record",
    );
}

/// One transcript carrying `core`, coding, with an untranslated 5' region — so
/// `cds_start != 1` and the CDS-relative axis is *numerically distinct* from
/// the transcript-relative one.
///
/// [`provider`]'s coding arm deliberately sets the CDS to the whole transcript
/// so that `c.N`, `n.N` and `r.N` name the same base, which is what makes its
/// coding/non-coding pair a controlled comparison. The cost is that its CDS
/// delta is **0**, so on those records the two conventions this fix chooses
/// between produce identical output — see
/// [`a_coding_transcripts_terminal_rna_repair_resolves_through_the_cds`].
fn provider_with_five_prime_utr(accession: &str, core: &str, cds_start: u64) -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence = core.to_string();
    let length = sequence.len() as u64;
    provider.add_transcript(Transcript::new(
        accession.to_string(),
        Some("TESTGENE".to_string()),
        Strand::Plus,
        sequence,
        Some(cds_start),
        Some(length),
        vec![Exon::new(1, 1, length)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

/// A **coding** record whose CDS does not start at transcript position 1 still
/// resolves its terminal `r.` repair through the CDS.
///
/// This is the guard that actually discriminates. The two
/// `..._is_unchanged` tests above use a CDS spanning the whole transcript, so
/// their CDS delta is 0 and the CDS-relative and transcript-relative
/// resolutions coincide — replacing the predicate with a bare
/// `span.region == Region::Rna`, i.e. applying the non-coding fallback to
/// *every* `r.` record, leaves all of them passing. On this fixture it does
/// not: the fallback yields `r.[53dup;*4delinsggg]`, which both fails to
/// coalesce the two members and names the last base under the wrong
/// convention.
///
/// The terminal path is the discriminating one because that is where the two
/// conventions disagree about *which position is the last base*
/// (`cds_axis_end` vs the transcript length). The junction path agrees on this
/// shape, which is why it is not also pinned here.
///
/// RE-BLESSED under `partition-is-the-unit-of-normalization` (DECIDED,
/// 2026-08-08, `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`):
/// was `NM_UTRT.1:r.53_54insagg`, now `NM_UTRT.1:r.[53dup;58_59dup]` — `53dup`
/// is separated from the pair by unchanged bases (`general.md:34`), and the pair
/// lands on one junction and merges.
///
/// **This row's discriminating content survives the move intact, and that is
/// the thing to check rather than the member count.** The wrong-convention
/// answer it exists to exclude is one that names the last base under the
/// transcript axis instead of the CDS axis — `*4delinsggg` in the note below.
/// The new answer names no `*` position at all, so the CDS resolution is still
/// what ran, and the non-coding fallback firing here would still be visible.
///
/// Sequence unchanged, measured with `spec_corpus::denotation_of` (via
/// `hgvs_to_spdi`, not the normalizer); the output is a fixed point.
#[test]
fn a_coding_transcripts_terminal_rna_repair_resolves_through_the_cds() {
    let provider = provider_with_five_prime_utr("NM_UTRT.1", TERMINAL_CORE, 5);
    assert_reaches_preserving(
        &provider,
        TERMINAL_CORE,
        "NM_UTRT.1:r.[53dup;58_59dup]",
        "NM_UTRT.1:r.[53dup;58_59dup]",
        "a coding record with a 5' UTR must still resolve its terminal repair \
         through the CDS; the non-coding fallback firing here names the last \
         base under the transcript axis instead (`*4delinsggg`, #1453)",
    );
    assert_reaches_preserving(
        &provider,
        TERMINAL_CORE,
        "NM_UTRT.1:r.[53dup;54dup;55dup]",
        "NM_UTRT.1:r.[53dup;58_59dup]",
        "and the authored spelling must reach it",
    );
}

/// One transcript carrying [`CORE`] annotated with a **stop codon but no start
/// codon** — `cds_end` present, `cds_start` absent.
fn provider_five_prime_incomplete(accession: &str) -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence = CORE.to_string();
    let length = sequence.len() as u64;
    provider.add_transcript(Transcript::new(
        accession.to_string(),
        Some("TESTGENE".to_string()),
        Strand::Plus,
        sequence,
        None,
        Some(length),
        vec![Exon::new(1, 1, length)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

/// The fallback keys on `cds_start` alone, so a record with `cds_end` but no
/// `cds_start` is transcript-relative too.
///
/// This is the slice where a predicate requiring *both* bounds absent still
/// reproduced the issue verbatim, emitting `r.[9dup;11dup;11dup]`. It is
/// reachable from real data rather than only by direct construction:
/// `CdotTranscript`'s CDS guard is an **or**
/// (`start_codon.is_some() || stop_codon.is_some()`) and then assigns
/// `(self.start_codon, self.stop_codon)` straight through, so a 5'-incomplete
/// transcript (`cds_start_NF`, #972) arrives with exactly this shape.
///
/// Pinned as agreement with the `n.` axis for the same reason the non-coding
/// pair above is: `region_sequence_delta` and `axis_frame` already treat this
/// record's `r.` axis as transcript-relative, so a repair that resolved it
/// through the CDS would read its bases under one convention and name them
/// under another.
#[test]
fn a_five_prime_incomplete_records_rna_axis_is_transcript_relative() {
    let provider = provider_five_prime_incomplete("NR_NOSTART.1");
    let tx = normalize(&provider, "NR_NOSTART.1:n.[9dup;10dup;11dup]");
    let rna = normalize(&provider, "NR_NOSTART.1:r.[9dup;10dup;11dup]");
    // RE-BLESSED under `partition-is-the-unit-of-normalization` (DECIDED,
    // 2026-08-08,
    // `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`): was
    // `NR_NOSTART.1:r.9_10insuaa`, now `NR_NOSTART.1:r.[9dup;10_11a[4]]`,
    // byte-identical to the non-coding reproducer's answer on the same core —
    // which is the point of this slice. `9dup` stays its own member
    // (`general.md:34`); `10dup` and `11dup` land on one junction and merge.
    //
    // The claim this row carries is unchanged: the write still *takes*, so the
    // repeated member is still gone. Sequence unchanged, measured with
    // `spec_corpus::denotation_of` (via `hgvs_to_spdi`, not the normalizer) —
    // `Denotation::Sequence`, equal to the input's, where a repeated member
    // gives `Denotation::NoSequence`. Output is a fixed point.
    assert_eq!(
        rna, "NR_NOSTART.1:r.[9dup;10_11a[4]]",
        "a record with `cds_end` but no `cds_start` has no CDS to resolve \
         through, so its `r.` axis numbers the transcript — resolving it \
         through the CDS instead reverts the write and leaves \
         `r.[9dup;11dup;11dup]`, one interbase point claimed twice (#1453)"
    );
    assert_eq!(
        rna,
        as_rna_spelling(&tx),
        "and it must agree with the `n.` axis of the same record (`n.` gave \
         `{tx}`)"
    );
    // The two properties the comment above states, asserted rather than
    // described. Run on the `n.` axis: `as_rna_spelling` is a re-spelling of
    // that same answer, so checking it here would measure the renderer twice
    // and the normalization once.
    assert_reaches_preserving(
        &provider,
        CORE,
        "NR_NOSTART.1:n.[9dup;10dup;11dup]",
        &tx,
        "the `n.` axis answer is the one the `r.` spelling above is derived from",
    );
}
