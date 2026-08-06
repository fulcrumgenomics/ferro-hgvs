//! W5005 — provenance for members the input reported separately.
//!
//! `DNA/delins.md:79-84` gives the spec's reason for describing nearby variants
//! individually: *"the two variants may have been reported (or might occur)
//! individually"*, and names the concrete risk — "an overlap with the
//! description of the combined variant might be missed in the annotation step
//! (database queries)". That is **provenance**: it lives in the input's
//! spelling, and neither the reference nor the observed bases carry any trace
//! of it.
//!
//! The normalized string cannot answer it. Making the output depend on which
//! spelling arrived is the definition of non-confluence, which is what #1235
//! exists to remove — and the cost is measured rather than assumed: a single
//! input-relative comparand in the normalizer's weight bound cost 427 converged
//! classes per shuffle direction over the 11,272-class corpus.
//!
//! So the spelling-dependence goes where it is harmless. The canonical string
//! stays a function of the denoted sequence; this diagnostic carries how the
//! variant was reported. Two spellings of one variant still normalize to one
//! string — they simply carry different diagnostics, and that is what lets a
//! downstream query recover the individually-reported form.
//!
//! The assertions below pin both halves of that contract: the warning appears
//! when members are coalesced, and **the normalized string is identical whether
//! or not it does**.

use ferro_hgvs::{parse_hgvs, JsonProvider, Normalizer};
use std::io::Write;

/// A genome-capable provider over one contig of cyclic `ACGT`, with `payload`
/// written 1-based at `pos1`.
fn provider(contig: &str, len: usize, pos1: usize, payload: &str) -> JsonProvider {
    let mut bases: Vec<u8> = "ACGT".bytes().cycle().take(len).collect();
    for (i, b) in payload.bytes().enumerate() {
        bases[pos1 - 1 + i] = b;
    }
    let doc = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [],
        "genomic_sequences": { contig: String::from_utf8(bases).unwrap() },
    });
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(doc.to_string().as_bytes()).unwrap();
    JsonProvider::from_json(f.path()).unwrap()
}

/// `(normalized string, W5005 codes emitted)`.
fn normalize(p: &JsonProvider, input: &str) -> (String, Vec<String>) {
    let variant = parse_hgvs(input).unwrap();
    let result = Normalizer::new(p.clone())
        .normalize_with_diagnostics(&variant)
        .unwrap();
    let codes = result
        .warnings
        .iter()
        .filter(|w| w.code() == "MEMBERS_COALESCED_FROM_REPORTED_FORM")
        .map(|w| w.to_string())
        .collect();
    (result.result.to_string(), codes)
}

/// How many `;`-separated members a rendered allele carries.
fn members(rendered: &str) -> usize {
    match (rendered.find('['), rendered.rfind(']')) {
        (Some(open), Some(close)) if close > open => rendered[open + 1..close].split(';').count(),
        _ => 1,
    }
}

/// Two adjacent substitutions the normalizer merges into one `delins`: the
/// input reported two members, the output describes one.
#[test]
fn coalescing_two_reported_members_emits_the_diagnostic() {
    let p = provider("c", 400, 200, "GCA");
    let (rendered, codes) = normalize(&p, "c:g.[200G>T;201C>G]");

    assert_eq!(
        members(&rendered),
        1,
        "precondition: this input must coalesce ({rendered})"
    );
    assert_eq!(codes.len(), 1, "expected exactly one W5005 ({codes:?})");
    let message = &codes[0];
    assert!(
        message.contains("input described 2 cis members"),
        "message must state the reported count: {message}"
    );
    assert!(
        message.contains("describes 1"),
        "message must state the normalized count: {message}"
    );
}

/// The load-bearing property: the diagnostic is *informational*. The normalized
/// string is what it was; nothing about emitting a warning moves it.
///
/// Asserted by normalizing the same input through both public exits — the one
/// that collects diagnostics and the one that does not — and requiring them to
/// agree byte for byte.
#[test]
fn the_diagnostic_does_not_change_the_normalized_string() {
    let p = provider("c", 400, 200, "GCA");
    for input in [
        "c:g.[200G>T;201C>G]",
        "c:g.200_201delinsTG",
        "c:g.[100T>C;300T>G]",
        "c:g.200G>T",
    ] {
        let variant = parse_hgvs(input).unwrap();
        let plain = Normalizer::new(p.clone())
            .normalize(&variant)
            .unwrap()
            .to_string();
        let (with_diagnostics, _) = normalize(&p, input);
        assert_eq!(
            plain, with_diagnostics,
            "the diagnostics exit moved the string for {input}"
        );
    }
}

/// Two spellings of one variant: the strings converge (confluence), and the
/// diagnostics differ (provenance). That is the whole design in one assertion —
/// the spelling-dependence lives in the metadata, never in the string.
#[test]
fn equivalent_spellings_converge_while_their_provenance_differs() {
    let p = provider("c", 400, 200, "GCA");
    let (split_form, split_codes) = normalize(&p, "c:g.[200G>T;201C>G]");
    let (span_form, span_codes) = normalize(&p, "c:g.200_201delinsTG");

    assert_eq!(
        split_form, span_form,
        "the two spellings must normalize to one string"
    );
    assert_eq!(split_codes.len(), 1, "the split spelling was coalesced");
    assert!(
        span_codes.is_empty(),
        "the span spelling reported one member and was not coalesced ({span_codes:?})"
    );
}

/// A member count that does not fall emits nothing. Guards against the warning
/// firing on every multi-member allele, which would make it useless noise.
#[test]
fn an_allele_whose_members_survive_emits_nothing() {
    let p = provider("c", 400, 200, "GCA");
    let (rendered, codes) = normalize(&p, "c:g.[100T>C;300T>G]");
    assert_eq!(
        members(&rendered),
        2,
        "precondition: these members are far apart and must survive ({rendered})"
    );
    assert!(codes.is_empty(), "unexpected W5005: {codes:?}");
}

/// A single-member input can never coalesce, so it can never warn — whatever
/// the normalizer does to its spelling.
#[test]
fn a_single_member_input_never_warns() {
    let p = provider("c", 400, 200, "GCA");
    for input in ["c:g.200G>T", "c:g.200_202delinsTGC", "c:g.200_201del"] {
        let (_, codes) = normalize(&p, input);
        assert!(codes.is_empty(), "unexpected W5005 for {input}: {codes:?}");
    }
}

/// A `trans` group is not a coalescing candidate: its members sit on different
/// molecules, so a change in their number would mean something else entirely.
#[test]
fn a_trans_group_is_not_reported_as_coalesced() {
    let p = provider("c", 400, 200, "GCA");
    let (_, codes) = normalize(&p, "c:g.[200G>T];[201C>G]");
    assert!(
        codes.is_empty(),
        "a trans group must not report coalescing: {codes:?}"
    );
}

/// The W-code is stable and resolves both ways through the public registry, so
/// a caller can select on `"W5005"` and get this warning back.
#[test]
fn the_w_code_round_trips_through_the_registry() {
    use ferro_hgvs::error_handling::ErrorType;

    let error_type = ErrorType::MembersCoalescedFromReportedForm;
    assert_eq!(error_type.code(), "W5005");
    assert_eq!(
        ErrorType::from_code("W5005"),
        Some(error_type),
        "W5005 must resolve back to its ErrorType"
    );
    assert!(
        !error_type.is_correctable(),
        "there is nothing to correct: the coalesced form IS the canonical one"
    );
}

/// The `p.` axis coalesces too, and must carry the same provenance.
///
/// `accession_and_axis` originally returned `None` for `Protein`, so the
/// warning was dropped silently on exactly the axis where a database query is
/// most likely to be looking for the individually reported form. Measured
/// before the fix: `p.[Ala100Gly;Ala101Gly]` normalized to
/// `p.Ala100_Ala101delinsGlyGly` — two reported members, one normalized — with
/// `warnings=[]`.
#[test]
fn a_protein_cis_allele_carries_the_provenance_warning() {
    use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};
    let provider = MockProvider::with_test_data();
    let variant = parse_hgvs("NP_000079.2:p.[Ala100Gly;Ala101Gly]").expect("fixture parses");
    let result = Normalizer::new(provider)
        .normalize_with_diagnostics(&variant)
        .expect("normalization succeeds");
    assert_eq!(
        result.result.to_string(),
        "NP_000079.2:p.Ala100_Ala101delinsGlyGly",
        "two adjacent protein substitutions coalesce into one delins"
    );
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code() == "MEMBERS_COALESCED_FROM_REPORTED_FORM"),
        "the coalesce must be reported on the `p.` axis as it is on the \
         nucleotide axes; got {:?}",
        result.warnings.iter().map(|w| w.code()).collect::<Vec<_>>()
    );
}
