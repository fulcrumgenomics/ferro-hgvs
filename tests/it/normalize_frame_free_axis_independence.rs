//! An axis with no reading frame normalizes without consulting transcript
//! annotation.
//!
//! # What this asserts
//!
//! Normalize one `g.` (and one `m.`) description twice — once against a provider
//! that can resolve the covering transcript *and* its CDS, once against a
//! provider that cannot — and the two outputs are identical. A genomic string
//! must be derivable from its own reference; if the answer moves when annotation
//! appears, the normalizer is importing information the reference does not
//! carry, and the same input would normalize differently for two callers holding
//! different reference bundles.
//!
//! Two seams decide it, and each is `false` for the genomic axes by
//! construction:
//!
//! - `Normalizer::apply_canonical_split`'s `codon_frame_aware` — `true` for
//!   `HgvsVariant::Cds`, asked of the provider for `HgvsVariant::Rna` (#1241: an
//!   `r.` on a non-coding transcript has no frame), `false` for everything else.
//! - `merge::axis_frame`'s `reading_frame` — `false` for `CisKind::Genome` and
//!   `CisKind::Mt` with no provider round-trip at all.
//!
//! Nothing enforces either, which is why this is worth a test: a rule that
//! reached for a covering transcript on the genomic axis would compile.
//!
//! **The two seams cross-check each other, and that bounds what this module can
//! catch.** Measured by mutating each in turn to consult the annotated provider:
//! with only the first leaking, `sequence_first_pass` re-derives the merged
//! delins back into the split and the output is unchanged; with only the second
//! leaking, `apply_canonical_split` has already split and the codon exception
//! does not re-merge. Under either single-seam leak all five tests here still
//! pass. Only when **both** leak does the answer move.
//!
//! **Which tests then fire depends on whether the leak consults the provider,
//! and that is the sharper bound.** Both figures below are measured, not
//! reasoned:
//!
//! - A **provider-mediated** leak on both genomic seams — one that reaches the
//!   annotated transcript's CDS — turns
//!   `a_genomic_description_normalizes_the_same_…` red, together with
//!   `the_genomic_case_is_one_a_leaked_reading_frame_would_have_changed`. This
//!   is the shape the module is written against and the A/B design catches it.
//! - A leak that simply hard-codes a reading frame for the genomic axis, without
//!   consulting the provider at all, changes **both** arms of the A/B
//!   identically, so every equality assertion here still passes. Only the exact
//!   string pinned by `the_genomic_case_…` catches that shape.
//!
//! So the value pin is not redundant with the A/B — it is the only assertion
//! covering a non-provider-mediated leak, and the A/B is the only one covering a
//! leak whose *effect* depends on the reference bundle. The `m.` assertion
//! answers to a leak on the `Mt` arm and to nothing else (measured: mutating
//! `CisKind::Mt` and `HgvsVariant::Mt` fires it alone).
//!
//! Note also what no mutation here can reach: the `ReferenceProvider` trait
//! exposes no genomic-region-to-transcript lookup, so a leak in `src/normalize/`
//! cannot find the covering transcript of a `g.` description by position
//! through the public trait at all. The provider-mediated mutation above had to
//! name the transcript accession outright. That is a statement about today's
//! trait surface, and it is exactly the thing that would change if such a lookup
//! were ever added — which is when this module starts earning its keep.
//!
//! The seams are named because they are where a leak would come from, not because
//! this module reaches into them. Every assertion below is on public behaviour:
//! parse, normalize twice, compare the strings.
//!
//! # The case is chosen so that a leak would be visible
//!
//! The fixture places `CGC` at genomic positions that are exactly one codon of a
//! transcript the annotated provider knows. Replacing them with `TGG` changes
//! the first and third base and leaves the middle one, so:
//!
//! ```text
//!   frame-free  g.265_267delinsTGG  ->  g.[265C>T;267C>G]     split
//!   framed      c.[10C>T;12C>G]     ->  c.10_12delinsTGG      merged
//! ```
//!
//! `general.md:34` splits two changes separated by an unchanged base; `:35`'s
//! one-amino-acid exception is what licenses the merge, and it needs a reading
//! frame. A genomic normalizer that leaked the frame would therefore emit the
//! merged form, which is a *different string* — the assertion is not vacuous.
//! The shape is `LRG_199:g.494841_494843delinsTGG`, whose bases are `c.145`-`c.147`
//! of `LRG_199t1`, rebuilt hermetically at codon 4 of a synthetic transcript.
//!
//! # The negative control is load-bearing
//!
//! An "identical" result proves nothing unless the harness can detect
//! annotation-dependence at all — two providers that differ in a way the
//! normalizer never sees would produce identical outputs for every input,
//! including ones that *should* differ. So the same provider pair is run against
//! the `c.` axis, where the dependence is legitimate and expected, and the
//! outputs must differ. Both sides return `Ok` there: the difference is the rule
//! that fired, not a provider error, which is the strongest form the control can
//! take.
//!
//! # What this does NOT claim
//!
//! **Nothing about the `n.` axis.** Whether `n.` on a *coding* transcript should
//! carry that transcript's reading frame is an open question — normalization and
//! projection disagree about it, and the disagreement is recorded but unruled.
//! The omission is deliberate: asserting either side here would pre-empt an
//! operator ruling by freezing it in a test. (`axis_frame` currently answers
//! `reading_frame: false` for `CisKind::Tx`; that is the status quo, not a
//! ruling, and this module does not pin it.)
//!
//! It also claims nothing about *which* output is correct. The property is
//! stability of the answer across reference bundles, not the answer itself; the
//! one exact string it does pin is the genomic split, and only to establish that
//! a leak would have been visible.

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer, ReferenceProvider};

/// Transcript sequence, 10 codons. Bases 10-12 are `CGC` — codon 4 — and are the
/// three the test edits. Their neighbours (`A` at 9, `T` at 13) stop the edit
/// shuffling out of the codon, so the case stays where the fixture puts it.
const TRANSCRIPT_SEQUENCE: &str = "ATGGCCTTACGCTAGGATCACTTGGACTAA";

/// Bases of `ACGT` flanking the transcript on each side, so the normalizer's
/// 100bp window stays in bounds. Same idea as `common::synthetic::padded`,
/// spelled out here because this module needs the genomic accession and the
/// transcript placement to be its own.
///
/// **255, not 256, and the value is load-bearing.** `merge::same_codon` asks
/// `(pos - 1) / 3` of the *description's own* coordinates, so a leaked reading
/// frame merges only when the two edited endpoints are codon-aligned in the frame
/// the leak happens to use. A flank of 255 makes the edited bases codon-aligned
/// on **both** frames at once — `c.10`/`c.12` on the transcript and `g.265`/
/// `g.267` on the chromosome — so the case detects a leak that projects into the
/// transcript *and* one that merely turns the coding flag on and reads genomic
/// positions. With 256 the second is invisible: no codon of the transcript can
/// land on a genomically codon-aligned position, and a mutation enabling the
/// coding path for `HgvsVariant::Genome` passes this module unchanged (measured).
const FLANK_LEN: u64 = 255;

/// 1-based genomic position of transcript base 1.
const TX_GENOMIC_START: u64 = FLANK_LEN + 1;

/// 1-based genomic position of `c.10` (== transcript base 10, since the CDS
/// starts at transcript base 1).
const CODON_GENOMIC_START: u64 = FLANK_LEN + 10;

const CHROMOSOME: &str = "NC_TEST.1";
const MITOCHONDRION: &str = "NC_012920.1";
const TRANSCRIPT: &str = "NM_TEST.1";

fn genomic_sequence() -> String {
    let mut flank = "ACGT".repeat(FLANK_LEN as usize / 4 + 1);
    flank.truncate(FLANK_LEN as usize);
    format!("{flank}{TRANSCRIPT_SEQUENCE}{flank}")
}

/// A provider holding the genomic sequence and nothing else: it cannot resolve
/// any transcript, so it cannot know where a CDS is.
fn without_annotation(contig: &str) -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence(contig, genomic_sequence());
    provider
}

/// The same genomic sequence, plus a transcript placed on `contig` that covers
/// the edited bases and declares a CDS putting them on a codon boundary.
///
/// The *only* difference from [`without_annotation`] is the transcript. The
/// genomic accession, its sequence and its coordinates are identical, so any
/// difference in a `g.` result is attributable to the annotation and to nothing
/// else.
fn with_annotation(contig: &str) -> MockProvider {
    let mut provider = without_annotation(contig);
    let tx_len = TRANSCRIPT_SEQUENCE.len() as u64;
    let tx_genomic_end = FLANK_LEN + tx_len;
    provider.add_transcript(Transcript::new(
        TRANSCRIPT.to_string(),
        Some("SYNTH".to_string()),
        Strand::Plus,
        TRANSCRIPT_SEQUENCE.to_string(),
        Some(1),
        Some(tx_len),
        vec![Exon::with_genomic(
            1,
            1,
            tx_len,
            TX_GENOMIC_START,
            tx_genomic_end,
        )],
        Some(contig.to_string()),
        Some(TX_GENOMIC_START),
        Some(tx_genomic_end),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

fn normalize(provider: MockProvider, input: &str) -> Result<String, String> {
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("`{input}` must parse: {e}"));
    Normalizer::new(provider)
        .normalize(&variant)
        .map(|normalized| normalized.to_string())
        .map_err(|e| e.to_string())
}

/// The fixture's own claims, checked rather than asserted in prose.
///
/// If the transcript stopped covering the edited bases, or the CDS stopped
/// putting them on a codon boundary, every "identical output" below would still
/// pass while measuring nothing — the annotated provider would simply have no
/// frame to leak.
#[test]
fn the_fixture_puts_the_edited_bases_on_a_codon_of_the_annotated_transcript() {
    let provider = with_annotation(CHROMOSOME);
    let transcript = provider
        .get_transcript(TRANSCRIPT)
        .expect("the annotated provider resolves the transcript");
    let cds_start = transcript.cds_start.expect("…and it declares a CDS");
    assert_eq!(cds_start, 1, "c.<n> maps 1:1 onto transcript base <n>");

    let exon = transcript.exons.first().expect("one exon");
    let (exon_g_start, exon_g_end) = (
        exon.genomic_start.expect("exon carries genomic coords"),
        exon.genomic_end.expect("exon carries genomic coords"),
    );
    assert!(
        exon_g_start <= CODON_GENOMIC_START && CODON_GENOMIC_START + 2 <= exon_g_end,
        "the edited bases must lie inside the annotated exon"
    );

    // c.10 is the first base of codon 4: (10 - 1) % 3 == 0.
    let cds_position = CODON_GENOMIC_START - TX_GENOMIC_START + 1;
    assert_eq!(cds_position, 10);
    assert_eq!(
        (cds_position - cds_start) % 3,
        0,
        "the three edited bases must be a whole codon, or there is no frame to leak"
    );

    // The edited endpoints must ALSO be codon-aligned in the genomic frame, or a
    // leak that turns the coding flag on without projecting is invisible here —
    // `merge::same_codon` reads `(pos - 1) / 3` of whatever coordinates the
    // description carries. See [`FLANK_LEN`].
    assert_eq!(
        (CODON_GENOMIC_START - 1) / 3,
        (CODON_GENOMIC_START + 1) / 3,
        "the genomic endpoints must share a codon under `same_codon`'s arithmetic too"
    );

    // And the bases really are `CGC`, so replacing them with `TGG` leaves the
    // middle one unchanged — the shape whose split-vs-merge answer depends on
    // whether a reading frame is in scope.
    let genome = genomic_sequence();
    let start = (CODON_GENOMIC_START - 1) as usize;
    assert_eq!(&genome[start..start + 3], "CGC");

    // The provider without annotation must genuinely lack it, or the A/B is
    // between two identical inputs.
    assert!(
        without_annotation(CHROMOSOME)
            .get_transcript(TRANSCRIPT)
            .is_err(),
        "the unannotated provider must not resolve the transcript"
    );
}

/// The invariant, on `g.`.
#[test]
fn a_genomic_description_normalizes_the_same_with_and_without_transcript_annotation() {
    for input in [
        &format!(
            "{CHROMOSOME}:g.{CODON_GENOMIC_START}_{}delinsTGG",
            CODON_GENOMIC_START + 2
        ),
        // The split spelling of the same variant, so the invariant is checked in
        // both directions rather than only on the input that happens to move.
        &format!(
            "{CHROMOSOME}:g.[{CODON_GENOMIC_START}C>T;{}C>G]",
            CODON_GENOMIC_START + 2
        ),
        // A `g.` carrying a gene selector — the genomic shape that most plausibly
        // hands the normalizer a route to the annotation.
        &format!(
            "{CHROMOSOME}(SYNTH):g.{CODON_GENOMIC_START}_{}delinsTGG",
            CODON_GENOMIC_START + 2
        ),
    ] {
        let bare = normalize(without_annotation(CHROMOSOME), input);
        let annotated = normalize(with_annotation(CHROMOSOME), input);
        assert_eq!(
            bare, annotated,
            "`{input}` normalized differently once the covering transcript and its CDS \
             became resolvable — a genomic description must be derivable from its own reference"
        );
    }
}

/// The same, on `m.`. The mitochondrial axis is the other frame-free DNA axis
/// `axis_frame` answers without a provider, and nothing else in this module
/// reaches it.
#[test]
fn a_mitochondrial_description_normalizes_the_same_with_and_without_transcript_annotation() {
    let input = format!(
        "{MITOCHONDRION}:m.{CODON_GENOMIC_START}_{}delinsTGG",
        CODON_GENOMIC_START + 2
    );
    assert_eq!(
        normalize(without_annotation(MITOCHONDRION), &input),
        normalize(with_annotation(MITOCHONDRION), &input),
        "`{input}` normalized differently once transcript annotation became resolvable"
    );
}

/// The genomic answer is the *split* form, which is what makes the equality
/// above worth asserting: a leaked reading frame would have merged it, so the
/// two providers would have disagreed.
///
/// It is also the only assertion in this module that survives a leak which
/// never consults the provider — such a leak moves both arms of the A/B
/// identically and every equality above stays green. See the module doc; do not
/// delete this as a redundant value pin.
#[test]
fn the_genomic_case_is_one_a_leaked_reading_frame_would_have_changed() {
    let input = format!(
        "{CHROMOSOME}:g.{CODON_GENOMIC_START}_{}delinsTGG",
        CODON_GENOMIC_START + 2
    );
    let expected = format!(
        "{CHROMOSOME}:g.[{CODON_GENOMIC_START}C>T;{}C>G]",
        CODON_GENOMIC_START + 2
    );
    assert_eq!(
        normalize(with_annotation(CHROMOSOME), &input).expect("normalizes"),
        expected,
        "the frame-free split is the whole point of the case; if this becomes the merged \
         delins, the equality assertions above stop distinguishing anything"
    );
}

/// **Negative control.** The same provider pair, on the axis where dependence on
/// the CDS is legitimate, must produce *different* answers — otherwise "the
/// outputs are identical" above is indistinguishable from a harness that cannot
/// see the provider difference at all.
///
/// Both sides return `Ok`, and both return a normalized description: with the
/// CDS, `general.md:35`'s one-amino-acid exception merges the two substitutions
/// into one `delins`; without it, there is no frame and `general.md:34` keeps
/// them apart. So what the control demonstrates is that a difference in
/// resolvable annotation reaches the *output* — not merely that one provider
/// errors where the other does not, which would say nothing about whether the
/// normalizer consults annotation at all.
#[test]
fn the_control_shows_the_harness_detects_annotation_dependence_on_the_coding_axis() {
    let input = "NM_TEST.1:c.[10C>T;12C>G]";
    let bare = normalize(without_annotation(CHROMOSOME), input);
    let annotated = normalize(with_annotation(CHROMOSOME), input);

    assert_ne!(
        bare, annotated,
        "the control did not fire: the coding axis answered identically with and without \
         the CDS, so this harness cannot detect annotation-dependence and the genomic \
         equality assertions prove nothing"
    );
    assert_eq!(
        bare.as_deref(),
        Ok("NM_TEST.1:c.[10C>T;12C>G]"),
        "with no reading frame in scope the two substitutions stay separate"
    );
    assert_eq!(
        annotated.as_deref(),
        Ok("NM_TEST.1:c.10_12delinsTGG"),
        "with the CDS resolvable the codon exception merges them"
    );
}
