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
#[test]
fn a_noncoding_rna_allele_collapses_instead_of_repeating_a_member() {
    let provider = provider("NR_TEST.1", false);
    assert_eq!(
        normalize(&provider, "NR_TEST.1:r.[9dup;10dup;11dup]"),
        "NR_TEST.1:r.9_10insuaa",
        "two members shifted onto the junction 3' of `r.11` must be coalesced \
         into one insertion; emitting `r.[9dup;11dup;11dup]` claims one \
         interbase point twice and so denotes no sequence (#1453)"
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
#[test]
fn a_coding_transcripts_rna_axis_is_unchanged() {
    let provider = provider("NM_TEST.1", true);
    assert_eq!(
        normalize(&provider, "NM_TEST.1:r.[9dup;10dup;11dup]"),
        "NM_TEST.1:r.9_10insuaa",
        "the coding `r.` axis already coalesced these members and must be \
         unaffected — it is the CDS-relative gap resolution that #1453 leaves \
         in place"
    );
    assert_eq!(
        normalize(&provider, "NM_TEST.1:c.[9dup;10dup;11dup]"),
        "NM_TEST.1:c.9_10insTAA",
        "and so must the `c.` axis of the same record"
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
#[test]
fn a_coding_transcripts_terminal_rna_repair_is_unchanged() {
    let provider = provider_for("NM_TERM.1", true, TERMINAL_CORE);
    assert_eq!(
        normalize(&provider, "NM_TERM.1:r.[57dup;58dup;59dup]"),
        "NM_TERM.1:r.57_58insagg"
    );
    assert_eq!(
        normalize(&provider, "NM_TERM.1:c.[57dup;58dup;59dup]"),
        "NM_TERM.1:c.57_58insAGG"
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
#[test]
fn a_coding_transcripts_terminal_rna_repair_resolves_through_the_cds() {
    let provider = provider_with_five_prime_utr("NM_UTRT.1", TERMINAL_CORE, 5);
    assert_eq!(
        normalize(&provider, "NM_UTRT.1:r.[53dup;54dup;55dup]"),
        "NM_UTRT.1:r.53_54insagg",
        "a coding record with a 5' UTR must still resolve its terminal repair \
         through the CDS; the non-coding fallback firing here gives \
         `r.[53dup;*4delinsggg]` (#1453)"
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
    assert_eq!(
        rna, "NR_NOSTART.1:r.9_10insuaa",
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
}
