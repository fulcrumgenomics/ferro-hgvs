//! Regression coverage for issue #1004.
//!
//! Two separate insertions at the **same reference junction** in a cis allele
//! (`g.[4_5insT;4_5insA]`) are an order-ambiguous *overlap conflict*, not a
//! normalizable form. Every reference agrees:
//!   - the HGVS spec expresses "insert both" as a single ordered compound
//!     payload `ins[T;A]` (general.md:79 / DNA/insertion.md), not two members;
//!   - mutalyzer rejects it (`EOVERLAP`);
//!   - VariantValidator rejects it (`AlleleSyntaxError: ... ranges overlap`);
//!   - ferro's strict mode rejects it (W5002 / OverlapConflictingEdits, #486).
//!
//! So in strict mode the input is rejected (covered by #486). In the non-strict
//! modes ferro must **warn and preserve** the allele as authored rather than
//! fabricating a merged/collapsed `delins` — the old behavior invented an
//! order-dependent form no other tool produces and, worse, was **not
//! idempotent** (`normalize(normalize(x)) != normalize(x)`): the fabricated
//! insertion 3'-shifted flush against a neighbour on the second pass and
//! collapsed. Preserving the overlap is both spec-consistent and idempotent.
//!
//! A pair of insertions at *different* junctions flanking an edit is a distinct
//! case — two changes less than two nucleotides apart (general.md:34-39) — and
//! still collapses to a single `delins`.

use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// Build a MockProvider whose genomic contig `NC_000001.11` carries
/// `substantive` starting at 1-based `start`, padded to 1000 bases of filler
/// `A` so the default normalize window always lands inside.
fn provider_with(start: usize, substantive: &str) -> MockProvider {
    let mut seq = vec![b'A'; 1000];
    for (i, b) in substantive.bytes().enumerate() {
        seq[start - 1 + i] = b;
    }
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NC_000001.11", String::from_utf8(seq).unwrap());
    provider
}

/// Non-strict (default) 3'-shifting normalizer over `substantive` at `start`.
fn normalize(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::with_config(
        provider,
        NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
    );
    let variant = parse_hgvs(input).expect("parse");
    format!("{}", normalizer.normalize(&variant).expect("normalize"))
}

/// Non-strict normalize that also returns the emitted warning codes, so the
/// *warn* half of "warn-and-preserve" can be asserted alongside the output.
fn normalize_with_warnings(provider: MockProvider, input: &str) -> (String, Vec<String>) {
    let normalizer = Normalizer::with_config(
        provider,
        NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
    );
    let variant = parse_hgvs(input).expect("parse");
    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("normalize");
    let codes = result
        .warnings
        .iter()
        .map(|w| w.code().to_string())
        .collect();
    (format!("{}", result.result), codes)
}

fn ref_mono() -> MockProvider {
    provider_with(300, "CATCCTCGCTCCT")
}
fn ref_di() -> MockProvider {
    provider_with(600, "TGACTTCAGTCACCTGACTGACTG")
}

// --- idempotency: normalize(normalize(x)) == normalize(x) -----------------

#[test]
fn samegap_insertion_pair_with_sub_is_idempotent() {
    let input = "NC_000001.11:g.[305_306insC;305_306insA;307G>T]";
    let pass1 = normalize(ref_mono(), input);
    let pass2 = normalize(ref_mono(), &pass1);
    assert_eq!(pass1, pass2, "normalize must be idempotent; pass1 != pass2");
}

#[test]
fn samegap_dinucleotide_insertion_pair_with_sub_is_idempotent() {
    let input = "NC_000001.11:g.[611_612insCG;611_612insGC;613C>A]";
    let pass1 = normalize(ref_di(), input);
    let pass2 = normalize(ref_di(), &pass1);
    assert_eq!(pass1, pass2, "normalize must be idempotent; pass1 != pass2");
}

#[test]
fn samegap_insertion_pair_with_deletion_is_idempotent() {
    let input = "NC_000001.11:g.[305_306insC;305_306insA;306del]";
    let pass1 = normalize(ref_mono(), input);
    let pass2 = normalize(ref_mono(), &pass1);
    assert_eq!(pass1, pass2, "normalize must be idempotent; pass1 != pass2");
}

// --- warn-and-preserve: the co-located overlap is left as authored, never
// fabricated into a merged insertion or a `delins` ---------------------------

#[test]
fn samegap_insertion_pair_with_sub_is_preserved_as_authored() {
    let input = "NC_000001.11:g.[305_306insC;305_306insA;307G>T]";
    let (out, warnings) = normalize_with_warnings(ref_mono(), input);
    assert_eq!(
        out, input,
        "a co-located insertion overlap must be preserved as authored, not canonicalized"
    );
    assert!(
        !out.contains("delins"),
        "must not collapse to a delins; got: {out}"
    );
    assert!(
        warnings.iter().any(|c| c == "OVERLAP_CONFLICTING_EDITS"),
        "warn-and-preserve: expected an OVERLAP_CONFLICTING_EDITS warning, got {warnings:?}"
    );
}

#[test]
fn samegap_insertion_pair_with_deletion_is_preserved_as_authored() {
    let input = "NC_000001.11:g.[305_306insC;305_306insA;306del]";
    let (out, warnings) = normalize_with_warnings(ref_mono(), input);
    assert_eq!(
        out, input,
        "a co-located insertion overlap must be preserved as authored"
    );
    assert!(
        !out.contains("delins"),
        "must not collapse to a delins; got: {out}"
    );
    assert!(
        warnings.iter().any(|c| c == "OVERLAP_CONFLICTING_EDITS"),
        "warn-and-preserve: expected an OVERLAP_CONFLICTING_EDITS warning, got {warnings:?}"
    );
}

// --- strict mode still rejects the overlap (spec / mutalyzer / VV consensus,
// #486). This locks the contract against my non-strict change. --------------

#[test]
fn strict_mode_still_rejects_samegap_insertion_overlap() {
    let normalizer = Normalizer::with_config(ref_mono(), NormalizeConfig::strict());
    let v = parse_hgvs("NC_000001.11:g.[305_306insC;305_306insA;307G>T]").expect("parse");
    let err = normalizer
        .normalize(&v)
        .expect_err("strict mode must reject a co-located insertion overlap");
    let msg = format!("{err}");
    assert!(
        msg.contains("W5002")
            || msg.contains("OverlapConflictingEdits")
            || msg.contains("coincident"),
        "expected an overlap-conflict rejection; got: {msg}"
    );
}

// --- regression guard: a pair at DIFFERENT junctions flanking an edit is a
// distinct-position case and still collapses to a single `delins` (spec
// general.md:34-39), and is idempotent. -------------------------------------

#[test]
fn different_junction_flanking_insertions_still_collapse_to_delins() {
    // insA at gap 304, insC at gap 305, deletion at 305: two changes less than
    // two nucleotides apart -> a single delins over the affected span.
    let input = "NC_000001.11:g.[304_305insA;305_306insC;305del]";
    let out = normalize(ref_mono(), input);
    assert!(
        out.contains("delins"),
        "distinct-junction insertions flanking an edit must collapse to a delins; got: {out}"
    );
    let again = normalize(ref_mono(), &out);
    assert_eq!(out, again, "the collapsed delins must itself be idempotent");
}
