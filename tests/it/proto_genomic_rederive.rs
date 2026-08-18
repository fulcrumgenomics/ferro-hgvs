//! Genomic rederive prototype (issue #2155).
//!
//! These tests exercise the prototype whose goal is to make genomic (`g.`/`m.`)
//! `normalize` re-derive from the resulting sequence through the confluent
//! sequence-first core, so every reference-aware entry point converges on one
//! description. The tests assert **convergence + idempotence**, never *which*
//! form the core picks — the form is what Phase 4 inspects and what the operator
//! may re-adjudicate (`// TODO(readjudicate #2155)`).

use ferro_hgvs::{parse_hgvs, FromSequencesOptions, MockProvider, Normalizer};

/// A genomic provider backed by a single contig's bases.
fn provider(contig: &str, seq: &str) -> MockProvider {
    let mut p = MockProvider::new();
    p.add_genomic_sequence(contig, seq);
    p
}

/// Reverse-complement of a DNA string.
fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            other => other,
        })
        .collect()
}

/// Normalize a parsed HGVS string and render it.
fn norm(p: &MockProvider, s: &str) -> String {
    Normalizer::new(p.clone())
        .normalize(&parse_hgvs(s).unwrap())
        .unwrap()
        .to_string()
}

#[test]
fn genomic_rederive_is_idempotent_on_the_2155_main_case() {
    // TEMPLATE contig from issue #2155. g. is 1-based inclusive, so g.10_17 is
    // the 0-based slice [9, 17) = "CTTAGTTA" (pos 10 = C). The #2155 main case
    // replaces that window with "AAACAAAC" (CTTAGTTA -> AAACAAAC).
    let refseq = "AGAACCCCCCTTAGTTAAGAACAAAAGCAACAATCTTCGTGGTCCTGG";
    assert_eq!(
        &refseq[9..17],
        "CTTAGTTA",
        "fixture sanity: g.10_17 == CTTAGTTA"
    );
    // Splice the observed full sequence so nothing is hand-transcribed.
    let altseq = format!("{}{}{}", &refseq[..9], "AAACAAAC", &refseq[17..]);
    let contig = "NC_TEST.1";
    let p = provider(contig, refseq);
    let opts = FromSequencesOptions::default();

    // (a) Two spellings via `normalize` converge, and the result is a fixed point.
    let split = norm(
        &p,
        &format!("{contig}:g.[10_12delinsAA;14_16delinsCAA;17_18insC]"),
    );
    let span = norm(&p, &format!("{contig}:g.10_17delinsAAACAAAC"));
    assert_eq!(
        split, span,
        "two spellings of one variant must converge under normalize"
    );
    assert_eq!(
        norm(&p, &split),
        split,
        "the canonical form must be a fixed point"
    );

    // (b) The OTHER reference-aware entry points must land on the SAME form.
    // Confluence-by-construction is a claim ACROSS entry points, not within one.
    let seq_norm = Normalizer::new(p.clone())
        .sequence_normalize(
            &parse_hgvs(&format!("{contig}:g.10_17delinsAAACAAAC")).unwrap(),
            &opts,
            true,
        )
        .unwrap()
        .to_string();
    assert_eq!(
        seq_norm, span,
        "sequence_normalize must converge to the same form as normalize"
    );

    // from_sequences — the (ref, obs) constructor the issue reproduces through.
    // Use the issue's reproduction: FULL sequences at pos 1 (the path that has
    // flanks to shift within — a windowed sub-call tests strictly less).
    let from_seq = Normalizer::new(p.clone())
        .from_sequences(contig, 1, refseq, &altseq, &opts, true)
        .unwrap()
        .to_string();
    assert_eq!(
        from_seq, span,
        "from_sequences must converge to the same form as normalize"
    );

    // FINDING (Phase 1): on current main all four converge to the min-edit SPLIT
    //   NC_TEST.1:g.[10_12delinsAA;14_16delinsCAA;17_18insC]
    // (NOT the spanning delins) — already confluent-by-construction via normalize's
    // post-hoc core. The form is a Phase-4 / re-adjudication question.
}

/// Phase 2: a whole-span reverse-complement block must type as ONE `inv`, in the
/// RAW `from_sequences` derivation (normalize=false) — the path the free Python
/// `ferro_hgvs.from_sequences` uses. Today this SHREDS the block into
/// per-column subs/indels because inv typing is per-member (post-cut) and cannot
/// coalesce an already-split block; the in-cut predicate types the whole block
/// before the partition so it is never split.
///
/// (`normalize` already produces `inv` for this via `coalesce_whole_block_inversion`,
/// so a `normalize`-level assertion would not fail — the raw path is the reproducer.)
#[test]
fn a_whole_span_inversion_is_not_shredded_in_the_raw_derivation() {
    // 16-mer block with interior self-coincidence, so the min-edit cut splits it.
    let block = "TTCAAGTTTTCAAGTT";
    let refseq = format!("AAAAAAAAAA{block}AAAAAAAAAA");
    let altseq = format!("AAAAAAAAAA{}AAAAAAAAAA", revcomp(block));
    let contig = "NC_TEST.1";
    let p = provider(contig, &refseq);
    let opts = FromSequencesOptions::default();
    // block occupies g.11..=g.26 (1-based inclusive).
    let end = 10 + block.len();
    let raw = Normalizer::new(p.clone())
        .from_sequences(contig, 1, &refseq, &altseq, &opts, false)
        .unwrap()
        .to_string();
    assert_eq!(
        raw,
        format!("{contig}:g.11_{end}inv"),
        "whole-span revcomp must be one inv in the raw derivation, not shredded"
    );
}

/// Phase 3 safety: the in-cut inversion predicate must NOT mis-type a
/// duplication as `inv`. A dup is a pure insertion (empty reference block), and
/// the predicate requires `block_ref.len() >= 2`, so it structurally cannot
/// fire — verified here on both `normalize` and the raw `from_sequences` path,
/// including the adversarial case of a duplicated unit that is its own
/// reverse-complement palindrome.
#[test]
fn a_duplication_survives_and_is_not_mistyped_as_inversion() {
    let contig = "NC_TEST.1";
    let opts = FromSequencesOptions::default();

    // Isolated tandem dup of ACGT at g.11_14.
    {
        let refseq = "GGGGGGGGGGACGTCCCCCCCCCC";
        let altseq = "GGGGGGGGGGACGTACGTCCCCCCCCCC";
        let p = provider(contig, refseq);
        let n = Normalizer::new(p.clone());
        assert_eq!(
            norm(&p, &format!("{contig}:g.11_14dup")),
            format!("{contig}:g.11_14dup")
        );
        let raw = n
            .from_sequences(contig, 1, refseq, altseq, &opts, false)
            .unwrap()
            .to_string();
        assert_eq!(
            raw,
            format!("{contig}:g.11_14dup"),
            "raw derivation of a dup is a dup"
        );
    }

    // Adversarial: the duplicated unit GAATTC is its own reverse complement, so a
    // naive whole-block revcomp test could mis-fire. It must still be a dup.
    {
        let unit = "GAATTC";
        assert_eq!(
            revcomp(unit),
            unit,
            "fixture: GAATTC is a revcomp palindrome"
        );
        let refseq = format!("GGGGGGGGGG{unit}CCCCCCCCCC");
        let altseq = format!("GGGGGGGGGG{unit}{unit}CCCCCCCCCC");
        let p = provider(contig, &refseq);
        let n = Normalizer::new(p.clone());
        assert_eq!(
            norm(&p, &format!("{contig}:g.11_16dup")),
            format!("{contig}:g.11_16dup")
        );
        let raw = n
            .from_sequences(contig, 1, &refseq, &altseq, &opts, false)
            .unwrap()
            .to_string();
        assert_eq!(
            raw,
            format!("{contig}:g.11_16dup"),
            "a palindrome-unit dup must stay a dup, not become an inv"
        );
    }
}
