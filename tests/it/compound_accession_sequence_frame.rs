//! A compound `NG_xxx(NM_yyy)` reference resolves its **coordinates** on the
//! inner transcript and must resolve its **bases** there too.
//!
//! # The two halves of one accession, and which one the numbers belong to
//!
//! `Accession` stores a compound reference the way HGVS writes it: the *inner*
//! accession is the primary one and the `NG_`/`NC_`/`LRG_` parent hangs off it
//! as `genomic_context`. [`Accession::transcript_accession`] exists for exactly
//! this — its doc calls it "the bare accession string for **provider/transcript
//! lookups**, stripping any genomic context wrapper" — and `src/spdi/convert.rs`
//! already uses it to resolve the position:
//!
//! ```text
//! resolve_cds_to_tx:  let tx_id = accession.transcript_accession();   // NM_000532.5
//!                     provider.get_transcript(&tx_id)                 // cds_start = 37
//!                     ensure_tx_in_bounds(.., transcript.sequence_length(), ..)
//! ```
//!
//! So `NG_008939.1(NM_000532.5):c.156` resolves as `cds_start + 156 - 1` = the
//! transcript's own base 192, bounded against the transcript's own 1799 bases.
//! The resulting number is an offset **on the transcript**.
//!
//! # What was wrong
//!
//! The sibling call two lines below handed `emit_spdi_for_edit` the *whole*
//! compound string — `variant.accession.to_string()`, i.e. `Display`, which
//! renders `NG_008939.1(NM_000532.5)`. That string is both the SPDI's `sequence`
//! field and the key every short-form edit's reference fetch is made against, so
//! a transcript offset was read out of whichever record the provider resolved
//! that string to. `MultiFastaProvider::resolve_name` reaches it through the
//! version-strip fallback (`"NG_008939.1(NM_000532.5)".split('.').next()` is
//! `"NG_008939"`), landing on the **genomic parent**:
//!
//! ```text
//! get_sequence("NG_008939.1",              191, 197) = CACACA   (the NG_ parent)
//! get_sequence("NG_008939.1(NM_000532.5)", 191, 197) = CACACA   <-- same record
//! get_sequence("NM_000532.5",              191, 197) = ACGCCG   <-- the frame the 191 is in
//! ```
//!
//! Every consumer of that SPDI then read the wrong bases. The one that made it
//! visible is `FERRO_ASSERT_SEQUENCE` (#1615), whose `compare_denoted_sequences`
//! fetches its comparison window from the SPDI's own `sequence` field: two
//! spellings that denote identical bases on the transcript are compared against
//! the parent's bases and can differ there, so the oracle reported a changed
//! denotation for a normalization that changed nothing.
//!
//! The tests below are hermetic — they build a parent and a transcript whose
//! bases differ at the positions in play, so they run in gating CI rather than
//! only under `FERRO_MANIFEST`.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::reference::ReferenceProvider;
use ferro_hgvs::spdi::{compare_denoted_sequences, hgvs_to_spdi, DenotedSequenceComparison};
use ferro_hgvs::{parse_hgvs, MockProvider};

/// The inner transcript. `cds_start = 3`, so `c.N` is transcript base `N + 2`.
const TRANSCRIPT: &str = "NM_INNER.1";
/// The genomic parent the description is written against.
const PARENT: &str = "NG_PARENT.1";

/// Transcript bases 1..=20. Chosen so the window the tests read is a
/// homopolymer run flanked by distinct bases — a 3' shift inside it is
/// sequence-preserving, which is what the oracle must be able to see.
///
/// ```text
///  1-based:  1234567890123456789 0
///            GG A AAA AA CGT CGTAC
/// ```
const TRANSCRIPT_BASES: &str = "GGAAAAAACGTCGTACGTAC";

/// Parent bases over the same *numeric* range, deliberately different at every
/// position the tests touch. Nothing about this record is a legitimate answer to
/// a transcript-frame coordinate; it is here so that reading it instead of the
/// transcript is detectable rather than accidentally right.
const PARENT_BASES: &str = "TTCCCCCCGATGATCAGATC";

fn transcript() -> Transcript {
    Transcript::new(
        TRANSCRIPT.to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        TRANSCRIPT_BASES.to_string(),
        Some(3),
        Some(20),
        vec![Exon::new(1, 1, TRANSCRIPT_BASES.len() as u64)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    )
}

/// A provider holding both records, with the parent *also* reachable under the
/// compound spelling.
///
/// The second registration models `MultiFastaProvider::resolve_name`'s
/// version-strip fallback rather than inventing a behaviour: measured on the
/// prepared reference, `get_sequence("NG_008939.1(NM_000532.5)", 191, 197)`
/// returns the parent's `CACACA`, byte-for-byte what `get_sequence` on the bare
/// `NG_008939.1` returns. Without it `MockProvider` would merely *error* on the
/// compound spelling, which turns the defect into a decline and hides it.
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_transcript(transcript());
    provider.add_genomic_sequence(PARENT, PARENT_BASES);
    provider.add_genomic_sequence(format!("{PARENT}({TRANSCRIPT})"), PARENT_BASES);
    provider
}

/// The premise, asserted rather than assumed: the two records really do differ
/// over the window in play, and the compound spelling really does serve the
/// parent's bases. If this ever stops holding, the tests below stop testing
/// anything and must be re-derived rather than deleted.
#[test]
fn the_compound_spelling_serves_the_parent_not_the_transcript() {
    let provider = provider();
    let compound = format!("{PARENT}({TRANSCRIPT})");

    let from_transcript = provider
        .get_sequence(TRANSCRIPT, 2, 8)
        .expect("transcript window");
    let from_parent = provider.get_sequence(PARENT, 2, 8).expect("parent window");
    let from_compound = provider
        .get_sequence(&compound, 2, 8)
        .expect("compound window");

    assert_eq!(from_transcript, "AAAAAA", "transcript bases 3..=8");
    assert_eq!(from_parent, "CCCCCC", "parent bases 3..=8");
    assert_eq!(
        from_compound, from_parent,
        "the compound spelling resolves to the parent's record — this is the \
         provider behaviour the SPDI conversion must not lean on"
    );
    assert_ne!(
        from_transcript, from_parent,
        "the fixture must distinguish the two frames"
    );
}

/// A `c.` position on a compound reference is resolved on the transcript, so the
/// SPDI it converts to must name the transcript and carry the transcript's
/// bases.
#[test]
fn a_compound_c_position_converts_to_spdi_on_the_transcript() {
    let provider = provider();
    // `cds_start = 3`, so `c.1` is transcript base 3 and `c.4_6` is transcript
    // bases 6..=8 — `AAA` on the transcript, `CCC` on the parent.
    let variant = parse_hgvs(&format!("{PARENT}({TRANSCRIPT}):c.4_6del")).expect("parses");
    let spdi = hgvs_to_spdi(&variant, &provider).expect("converts");

    assert_eq!(
        spdi.sequence, TRANSCRIPT,
        "the SPDI must name the sequence its position is an offset on"
    );
    assert_eq!(spdi.position, 5, "c.4 is transcript base 6 = interbase 5");
    assert_eq!(
        spdi.deletion, "AAA",
        "the deleted bases must come from the transcript, not the parent's CCC"
    );
}

/// The same for the `n.` axis, which reaches `emit_spdi_for_edit` by its own
/// path.
#[test]
fn a_compound_n_position_converts_to_spdi_on_the_transcript() {
    let provider = provider();
    let variant = parse_hgvs(&format!("{PARENT}({TRANSCRIPT}):n.6_8del")).expect("parses");
    let spdi = hgvs_to_spdi(&variant, &provider).expect("converts");

    assert_eq!(spdi.sequence, TRANSCRIPT);
    assert_eq!(spdi.position, 5);
    assert_eq!(spdi.deletion, "AAA");
}

/// And for `r.`, which reaches `emit_spdi_for_edit` by a third path.
///
/// The descriptor is `r.4_6`, not `r.6_8`: on a **coding** transcript ferro
/// resolves an `r.` position CDS-anchored (`cds_start + N - 1`, as for `c.`)
/// rather than transcript-anchored. Whether that is right is #1177's question
/// and is orthogonal to this module — what matters here is only that whichever
/// number comes out, the accession beside it names the sequence that number
/// indexes.
#[test]
fn a_compound_r_position_converts_to_spdi_on_the_transcript() {
    let provider = provider();
    let variant = parse_hgvs(&format!("{PARENT}({TRANSCRIPT}):r.4_6del")).expect("parses");
    let spdi = hgvs_to_spdi(&variant, &provider).expect("converts");

    assert_eq!(spdi.sequence, TRANSCRIPT);
    assert_eq!(spdi.position, 5);
    assert_eq!(spdi.deletion, "AAA");
}

/// A bare transcript accession is the control: it was never affected, and must
/// stay exactly where it was.
#[test]
fn a_bare_transcript_accession_is_unchanged() {
    let provider = provider();
    let variant = parse_hgvs(&format!("{TRANSCRIPT}:c.4_6del")).expect("parses");
    let spdi = hgvs_to_spdi(&variant, &provider).expect("converts");

    assert_eq!(spdi.sequence, TRANSCRIPT);
    assert_eq!(spdi.position, 5);
    assert_eq!(spdi.deletion, "AAA");
}

/// The consequence that made this visible: two spellings of one 3'-shifted
/// deletion inside the transcript's `AAAAAA` run denote the same bases, and the
/// denoted-sequence oracle must say so.
///
/// The run is transcript bases 3..=8 (`c.1`..`c.6`). Deleting any single base of
/// it yields the same sequence, so `c.1del` and `c.6del` are the same variant.
/// On the **parent** the corresponding window is `CCCCCC` — also a run, which is
/// why this pair alone would not discriminate; the pair below is chosen so the
/// two frames disagree.
#[test]
fn the_oracle_agrees_on_a_shift_that_preserves_the_transcript_sequence() {
    let provider = provider();
    // Transcript 8..=11 is `ACGT`; parent 8..=11 is `GATG`. Deleting transcript
    // base 8 (`A`, the last of the run) is the same variant as deleting base 3
    // — a legal 3' shift on the transcript. On the parent, base 8 is `C` and
    // base 3 is `C` too, so this pair agrees on both frames; what does NOT
    // agree is the *stated* base, which is what the comparison splices.
    let input = parse_hgvs(&format!("{PARENT}({TRANSCRIPT}):c.1delA")).expect("parses");
    let output = parse_hgvs(&format!("{PARENT}({TRANSCRIPT}):c.6del")).expect("parses");

    let verdict = compare_denoted_sequences(&input, &output, &provider);
    assert_eq!(
        verdict,
        DenotedSequenceComparison::Agree,
        "`c.1delA` states the transcript's `A`; `c.6del` reads its deleted base \
         from the reference. Both denote one sequence on the transcript — the \
         frame their coordinates are in."
    );
}

/// The oracle must still *fire* on a compound reference when the sequence
/// really does change. A fix that made every compound row skip would remove the
/// false positives and the coverage together.
#[test]
fn the_oracle_still_fires_on_a_compound_reference_that_changes_the_sequence() {
    let provider = provider();
    // Transcript `c.7` is base 9 = `C`; `c.8` is base 10 = `G`. Substituting at
    // two different bases is not one variant under any frame.
    let input = parse_hgvs(&format!("{PARENT}({TRANSCRIPT}):c.7C>T")).expect("parses");
    let output = parse_hgvs(&format!("{PARENT}({TRANSCRIPT}):c.8G>T")).expect("parses");

    let verdict = compare_denoted_sequences(&input, &output, &provider);
    assert!(
        matches!(verdict, DenotedSequenceComparison::Differ { .. }),
        "a genuine sequence change on a compound reference must still fire, got {verdict:?}"
    );
    if let DenotedSequenceComparison::Differ { accession, .. } = verdict {
        assert_eq!(
            accession, TRANSCRIPT,
            "the report must name the frame it compared in"
        );
    }
}
