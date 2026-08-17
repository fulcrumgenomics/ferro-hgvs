//! `to_sequences` and `from_sequences` are each other's inverse, up to the
//! derivation the second one performs.
//!
//! The pair is also its own oracle, which is what makes the corpus tests beside
//! this file possible without a single new fixture: any HGVS corpus becomes a
//! `from_sequences` corpus by passing it through `to_sequences` first.

use crate::common::synthetic::{padded, SyntheticBuilder, PAD_OFFSET};
use ferro_hgvs::{from_sequences, FromSequencesOptions, Normalizer};

/// A synthetic genomic normalizer over `core`, plus the contig it serves.
fn genomic(core: &str) -> (Normalizer<ferro_hgvs::JsonProvider>, String) {
    let provider = SyntheticBuilder::genomic(core).build();
    (Normalizer::new(provider), padded(core))
}

/// `to_sequences` reports **1-based** positions, the convention
/// `from_sequences` consumes.
///
/// Asserted against a hand-computed coordinate rather than against the other
/// half of the pair: an off-by-one that is consistent between them would be
/// invisible to a round-trip test and would shift every derived description by
/// one base.
#[test]
fn to_sequences_reports_one_based_positions() {
    let core = "ACGTACGTAAGCATGCTAGCATCG";
    let (normalizer, contig) = genomic(core);
    // `PAD_OFFSET + 1` is the first 1-based position of the core.
    let first_core_base = PAD_OFFSET + 1;
    let base = contig.as_bytes()[PAD_OFFSET as usize] as char;
    let substitution = format!(
        "NC_TEST.1:g.{first_core_base}{base}>{}",
        if base == 'A' { 'T' } else { 'A' }
    );

    let variant = ferro_hgvs::parse_hgvs(&substitution).expect("parses");
    let pair = normalizer.to_sequences(&variant, 0).expect("to_sequences");

    assert_eq!(pair.accession, "NC_TEST.1");
    assert_eq!(pair.position, first_core_base);
    assert_eq!(pair.reference, base.to_string());
}

/// The pair composes: taking a description to its sequences and deriving from
/// them yields a description denoting the same bases over the same window.
///
/// The comparison is on the **denoted bases**, not on the strings. Deriving is
/// allowed to change the string — that is what it is for — and asserting string
/// equality here would pin today's normalizer output rather than the invariant.
#[test]
fn a_description_survives_a_round_trip_through_its_sequences() {
    let core = "ACGTACGTAAGCATGCTAGCATCG";
    let (normalizer, _contig) = genomic(core);
    let input = format!("NC_TEST.1:g.[{}del;{}del]", PAD_OFFSET + 5, PAD_OFFSET + 8);

    let variant = ferro_hgvs::parse_hgvs(&input).expect("parses");
    let pair = normalizer
        .to_sequences(&variant, 128)
        .expect("to_sequences");

    let derived = from_sequences(
        &pair.accession,
        pair.position,
        &pair.reference,
        &pair.alternate,
        &FromSequencesOptions::default(),
    )
    .expect("from_sequences");

    let back = normalizer
        .to_sequences(&derived, 128)
        .expect("to_sequences again");
    assert_eq!(
        (back.position, &back.reference, &back.alternate),
        (pair.position, &pair.reference, &pair.alternate),
        "the derived description must denote the same bases over the same window"
    );
}

/// Deriving twice changes nothing.
///
/// Idempotence is not implied by determinism: a deterministic function can
/// still move its own output on a second pass, and that is precisely what the
/// per-member pipeline does for the shapes this surface exists for.
#[test]
fn deriving_twice_reaches_the_same_description() {
    let core = "ACGTACGTAAGCATGCTAGCATCG";
    let (normalizer, _contig) = genomic(core);
    let input = format!("NC_TEST.1:g.{}_{}delinsGC", PAD_OFFSET + 5, PAD_OFFSET + 8);

    let variant = ferro_hgvs::parse_hgvs(&input).expect("parses");
    let pair = normalizer
        .to_sequences(&variant, 128)
        .expect("to_sequences");
    let first = from_sequences(
        &pair.accession,
        pair.position,
        &pair.reference,
        &pair.alternate,
        &FromSequencesOptions::default(),
    )
    .expect("derives");

    let after = normalizer.to_sequences(&first, 128).expect("re-applies");
    let second = from_sequences(
        &after.accession,
        after.position,
        &after.reference,
        &after.alternate,
        &FromSequencesOptions::default(),
    )
    .expect("derives again");

    assert_eq!(first.to_string(), second.to_string());
}

/// The method refuses an interval running past the end of the sequence — a
/// validity failure the free function cannot see, because seeing it needs the
/// reference and that would make the provider a hidden input.
#[test]
fn the_method_refuses_an_interval_past_the_sequence_end() {
    let core = "ACGTACGTAAGCATGCTAGCATCG";
    let (normalizer, contig) = genomic(core);
    let past_end = contig.len() as u64 + 10;

    let err = normalizer
        .from_sequences(
            "NC_TEST.1",
            past_end,
            "ACGT",
            "A",
            &FromSequencesOptions::default(),
            false,
        )
        .expect_err("refuses");
    assert!(
        err.to_string().contains("exceeds available reference"),
        "{err}"
    );
}

/// Post-normalizing a derived description is a no-op **on this material**.
///
/// Scoped in the name and meant literally: an internal sweep of many synthetic
/// shapes found `normalize` moves a meaningful share of derived descriptions, so
/// this is a property of the rows below and never of the surface. See
/// `Normalizer::from_sequences`.
///
/// Pinned as a **measured property of today's normalizer**, not asserted as a
/// law: the derivation and the normalizer are separate rules, and the value of
/// this test is that it goes red if they stop agreeing, which is a fact worth
/// knowing rather than a contract worth defending.
#[test]
fn post_normalizing_a_derived_description_is_a_no_op_here() {
    let core = "ACGTACGTAAGCATGCTAGCATCG";
    let (normalizer, contig) = genomic(core);
    let start = PAD_OFFSET + 5;
    let window = &contig[(start - 1) as usize..(start + 3) as usize];
    let observed: String = window
        .chars()
        .enumerate()
        .filter(|(i, _)| *i != 0 && *i != 3)
        .map(|(_, c)| c)
        .collect();

    let bare = normalizer
        .from_sequences(
            "NC_TEST.1",
            start,
            window,
            &observed,
            &FromSequencesOptions::default(),
            false,
        )
        .expect("derives");
    let normalized = normalizer
        .from_sequences(
            "NC_TEST.1",
            start,
            window,
            &observed,
            &FromSequencesOptions::default(),
            true,
        )
        .expect("derives and normalizes");

    assert_eq!(bare.to_string(), normalized.to_string());
}

/// A provider that serves the **front** of any `get_sequence` request but at
/// most `limit` bases, delegating everything else.
///
/// It is not a hypothetical: `get_sequence`'s contract does not promise the full
/// range, and a real provider near a sequence start, at a FASTA record boundary,
/// or with a partial index short-reads exactly this way. `fetch_window`
/// (`spdi/apply.rs`) already guards against it by requiring
/// `bases.len() == end - start`.
struct ShortReadProvider<P> {
    inner: P,
    limit: usize,
}

impl<P: ferro_hgvs::ReferenceProvider> ferro_hgvs::ReferenceProvider for ShortReadProvider<P> {
    fn get_transcript(
        &self,
        id: &str,
    ) -> Result<std::sync::Arc<ferro_hgvs::reference::Transcript>, ferro_hgvs::FerroError> {
        self.inner.get_transcript(id)
    }

    fn get_sequence(
        &self,
        id: &str,
        start: u64,
        end: u64,
    ) -> Result<String, ferro_hgvs::FerroError> {
        let full = self.inner.get_sequence(id, start, end)?;
        Ok(full.chars().take(self.limit).collect())
    }

    fn get_genomic_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<String, ferro_hgvs::FerroError> {
        self.inner.get_genomic_sequence(contig, start, end)
    }

    fn has_genomic_data(&self) -> bool {
        self.inner.has_genomic_data()
    }

    fn get_sequence_length(&self, id: &str) -> Result<u64, ferro_hgvs::FerroError> {
        self.inner.get_sequence_length(id)
    }
}

/// **A short 5' read must not move the window's reported start.**
///
/// `prepend_five_prime_flank` asks for `[start - wanted, start)` and then
/// prepends whatever comes back. A provider that serves fewer bases has served
/// the **front** of that range, but labelling them as the *last* `served` bases
/// — which is what `start - served` does — shifts every coordinate the
/// subsequent derivation names by `wanted - served`, silently: the pair is still
/// well-formed and the description still parses, so nothing downstream can tell.
///
/// The fix is `fetch_window`'s own discipline — serve it exactly or not at all —
/// so a short read costs the 5' flank and nothing else. Asserted against the
/// full-service window rather than against a hand-computed number, since the
/// property is "the two agree about where the bases are", not a particular pad.
#[test]
fn a_short_five_prime_read_does_not_shift_the_window() {
    let core = "ACGTACGTAAGCATGCTAGCATCG";
    let input = format!("NC_TEST.1:g.{}del", PAD_OFFSET + 8);
    let variant = ferro_hgvs::parse_hgvs(&input).expect("parses");

    let honest = Normalizer::new(SyntheticBuilder::genomic(core).build());
    let full = honest.to_sequences(&variant, 32).expect("to_sequences");

    // Serves 8 of the 32 flank bases asked for.
    let short = Normalizer::new(ShortReadProvider {
        inner: SyntheticBuilder::genomic(core).build(),
        limit: 8,
    });
    let clipped = short.to_sequences(&variant, 32).expect("to_sequences");

    // The clipped window is allowed to be narrower — the flank is best effort —
    // but wherever it starts, that position must name its own first base.
    let full_first = full.position;
    let clipped_first = clipped.position;
    assert!(
        clipped_first >= full_first,
        "a short read cannot widen the window: {clipped_first} < {full_first}"
    );
    let offset = (clipped_first - full_first) as usize;
    assert_eq!(
        clipped.reference,
        full.reference[offset..],
        "the clipped window's position does not name its own first base — a short 5' read \
         shifted the coordinate frame"
    );
    assert_eq!(clipped.alternate, full.alternate[offset..]);
    // And the derivation over each names the same locus.
    let options = FromSequencesOptions::default();
    assert_eq!(
        from_sequences(
            &clipped.accession,
            clipped.position,
            &clipped.reference,
            &clipped.alternate,
            &options
        )
        .expect("derives")
        .to_string(),
        from_sequences(
            &full.accession,
            full.position,
            &full.reference,
            &full.alternate,
            &options
        )
        .expect("derives")
        .to_string(),
    );
}

/// **The round trip survives a soft-masked reference.**
///
/// `fetch_window` does not case-fold and `prepend_five_prime_flank` upper-cased
/// only the flank, so a masked region produced a mixed-case pair no caller
/// wrote; `apply_edits_to_window` then copied the reference verbatim while
/// splicing an upper-case payload, and `verify_round_trip`'s byte comparison
/// refused a correct derivation with "expected atgt, got aTgt".
///
/// A substitution is the shape that broke — `del`/`dup` carry no payload and
/// passed throughout, which is why a masked reference looked fine on the corpora.
#[test]
fn a_soft_masked_reference_round_trips() {
    let core = "acgtacgtaagcatgctagcatcg";
    let (normalizer, _contig) = genomic(core);
    let position = PAD_OFFSET + 5;
    let input = format!("NC_TEST.1:g.{position}A>T");
    let variant = ferro_hgvs::parse_hgvs(&input).expect("parses");

    let pair = normalizer
        .to_sequences(&variant, 128)
        .expect("to_sequences");
    assert_eq!(
        pair.reference,
        pair.reference.to_ascii_uppercase(),
        "to_sequences must hand back one alphabet, not a mixed-case window"
    );

    let derived = from_sequences(
        &pair.accession,
        pair.position,
        &pair.reference,
        &pair.alternate,
        &FromSequencesOptions::default(),
    )
    .expect("a soft-masked window must derive");
    assert_eq!(derived.to_string(), input);

    // And the caller's own lower-case window derives the same description, which
    // is the half `to_sequences` folding does not cover.
    let lower = from_sequences(
        &pair.accession,
        pair.position,
        &pair.reference.to_ascii_lowercase(),
        &pair.alternate.to_ascii_lowercase(),
        &FromSequencesOptions::default(),
    )
    .expect("a lower-case window must derive");
    assert_eq!(lower.to_string(), input);
}
