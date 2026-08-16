//! Insertions whose payload is given as a **range of the reference** rather
//! than as literal bases — `ins201_212`, `delins201_212`, `ins201_212inv`,
//! `ins[201_212]`, `ins[OTHER:g.a_b]`.
//!
//! # What the spec says, and what it conspicuously does not
//!
//! The form is permitted, in one sentence, as one of two alternatives —
//! `DNA/insertion.md:22`: the inserted sequence "can be given as the nucleotides
//! inserted (e.g., `insAGC`) or, for larger insert sequences, by referring to
//! the sequence in the reference sequence (e.g., `c.849_850ins858_895`)".
//!
//! "larger" is defined nowhere in the corpus, and `open-issues.md:27` cuts
//! against reading a threshold into it: "the formats for description are
//! general, **irrespective of the length of the variant**". The committee then
//! states the gap outright at `open-issues.md:77-78`: "Since the current
//! recommendations do not specify when to use which of these formats, one
//! variant can be described using different formats. This is undesired; HGVS
//! recommendations should be extended by specifying when to use which format."
//!
//! So: **neither form is preferred, mandated, or forbidden** for a
//! fully-sequenced insert. Two clauses that look like they bear on it do not —
//! `checklist.md:33` (`ins6`) and `DNA/insertion.md:119` (`ins24`) forbid a bare
//! *count*, which does not define the inserted sequence at all, whereas a range
//! does. Do not cite either against a copy-range.
//!
//! # What this module therefore pins
//!
//! That ferro **reads** every published spelling correctly and resolves it to
//! the right bases — which is the part that is unambiguously required, because a
//! copy-range denotes a specific sequence and getting it wrong is a wrong
//! answer, not a style choice.
//!
//! That ferro **writes** the literal form, always, and says so via
//! `INSERTED_SEQUENCE_EXPANDED`. That is a deliberate divergence (tracked as
//! `ferro-policy-654-positional-insert-literal-expansion` in the Mutalyzer
//! fixture dispositions; Mutalyzer preserves the copy-range on `c.` and
//! inconsistently on `g.`), and it follows the ledger's
//! `canonical-form-choice-when-both-legal`: derive from the resulting sequence,
//! do not preserve the input's spelling. It is pinned as *policy*, not as
//! compliance, because the spec ranks neither form.
//!
//! The one place the choice is **not** free is an inverted copy — see
//! `insertion_adjacency_defects::inverted_copy_range_is_the_recommended_form`.

use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::{parse_hgvs, HgvsVariant, MockProvider, NormalizeConfig, Normalizer};

use crate::common::hg38_window::{HG38_WINDOW, LOCAL_CONTIG};

/// A second contig, so the cross-reference spelling has somewhere to point.
///
/// Its sequence is the window reversed — not reverse-complemented, just
/// reversed — so that a payload fetched from it can never be confused with the
/// same coordinates fetched from the primary contig, and a resolver that
/// silently ignored the accession would produce visibly wrong bases rather than
/// coincidentally right ones.
const OTHER_CONTIG: &str = "NC_TESTALT.1";

fn other_sequence() -> String {
    HG38_WINDOW.chars().rev().collect()
}

fn two_contig_provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence(LOCAL_CONTIG, HG38_WINDOW.to_string());
    provider.add_genomic_sequence(OTHER_CONTIG, other_sequence());
    provider
}

/// Normalize and return `(rendered body, warning codes)`.
fn normalize_with_codes(body: &str) -> (String, Vec<String>) {
    let input = format!("{LOCAL_CONTIG}:g.{body}");
    let variant: HgvsVariant = parse_hgvs(&input).expect("parse");
    let normalizer = Normalizer::with_config(
        two_contig_provider(),
        NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
    );
    let out = normalizer
        .normalize_with_diagnostics(&variant)
        .unwrap_or_else(|e| panic!("{body}: {e}"));
    let rendered = out
        .result
        .to_string()
        .strip_prefix(&format!("{LOCAL_CONTIG}:g."))
        .expect("accession and axis survive rendering")
        .to_string();
    let codes = out
        .warnings
        .iter()
        .map(|w| w.code().to_string())
        .collect::<Vec<_>>();
    (rendered, codes)
}

#[track_caller]
fn assert_expands_to(body: &str, expected: &str) {
    let (got, codes) = normalize_with_codes(body);
    assert_eq!(got, expected, "expansion drifted for {body}");
    assert!(
        codes.iter().any(|c| c == "INSERTED_SEQUENCE_EXPANDED"),
        "{body} expanded silently; codes were {codes:?}"
    );
}

/// The bases the copy-range cases below name, stated once.
///
/// Written out rather than sliced so the expected payload in each case is
/// checkable against the assembly by eye, and so a fixture edit that shifted the
/// window would fail here rather than diffusely.
const COPY_201_212: &str = "TTTCTGAGCCAG";

#[test]
fn the_copy_source_is_the_bases_the_cases_assume() {
    assert_eq!(&HG38_WINDOW[200..212], COPY_201_212);
}

// ---------------------------------------------------------------------------
// Reading: every published spelling resolves to the same bases.
// ---------------------------------------------------------------------------

/// The bare form, `ins<start>_<end>`.
#[test]
fn bare_copy_range_resolves_to_the_named_bases() {
    assert_expands_to("302_303ins201_212", &format!("302_303ins{COPY_201_212}"));
}

/// The bracketed form. `general.md:78` gives `[ ]` for inserted-sequence lists,
/// and a single-element list is the same content — so it must resolve
/// identically, not merely similarly.
#[test]
fn bracketed_copy_range_resolves_identically_to_the_bare_form() {
    let bare = normalize_with_codes("302_303ins201_212").0;
    let bracketed = normalize_with_codes("302_303ins[201_212]").0;
    assert_eq!(bare, bracketed);
}

/// `delins` with a copy-range payload — the SVD-WG009 shape.
///
/// SVD-WG009 (accepted) retired `con` with the instruction to "simply replace
/// 'con' by 'delins'", and every example it publishes keeps a copy-range
/// payload (`SVD-WG009.md:22`, `DNA/delins.md:54`). So this spelling is the
/// current recommendation for the whole conversion class and has to read.
#[test]
fn delins_with_a_copy_range_payload_resolves() {
    let (got, codes) = normalize_with_codes("302_313delins201_212");
    assert!(codes.iter().any(|c| c == "INSERTED_SEQUENCE_EXPANDED"));
    assert_eq!(got, "[303del;305_306del;310_311insA;313_314insAG]");
}

/// The inverted copy-range resolves to the reverse complement.
///
/// A resolver that returned the forward bases would still produce a
/// length-correct, plausible-looking answer, which is why this asserts the
/// actual complement rather than just "differs from the forward form".
#[test]
fn inverted_copy_range_resolves_to_the_reverse_complement() {
    let (got, _) = normalize_with_codes("302_303ins201_212inv");
    let rc: String = COPY_201_212
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            other => other,
        })
        .collect();
    let payload = got.split("ins").nth(1).expect("an ins payload");
    assert_eq!(payload.len(), rc.len(), "inverted payload changed length");
    assert_ne!(payload, COPY_201_212, "inv resolved to the forward bases");
    // The rendered insertion may sit at a 3'-shifted junction, so the payload is
    // a rotation of the reverse complement rather than equal to it.
    assert!(
        format!("{rc}{rc}").contains(payload),
        "{payload} is not a rotation of the reverse complement {rc}"
    );
}

/// A copy-range far from the edit it feeds.
///
/// Distance matters because the resolver and the normalizer's own 100-base
/// window are different mechanisms; a payload fetch that reused the shift window
/// would fail here and pass everywhere else.
#[test]
fn a_distant_copy_range_resolves() {
    let (got, _) = normalize_with_codes("302_303ins500_511");
    assert_eq!(got, "302_303insGAAAGCCGCGGG");
    assert_eq!(&HG38_WINDOW[499..511], "GAAAGCCGCGGG");
}

/// A cross-reference payload resolves against the **named** accession.
#[test]
fn cross_reference_copy_range_resolves_against_the_named_accession() {
    let (got, codes) = normalize_with_codes(&format!("302_303ins[{OTHER_CONTIG}:g.201_212]"));
    assert!(codes.iter().any(|c| c == "INSERTED_SEQUENCE_EXPANDED"));
    let payload = got.split("ins").nth(1).expect("an ins payload");
    let expected = &other_sequence()[200..212];
    assert!(
        format!("{expected}{expected}").contains(payload),
        "cross-reference payload {payload} does not come from {OTHER_CONTIG} ({expected})"
    );
    assert_ne!(
        payload, COPY_201_212,
        "cross-reference resolved against the wrong contig"
    );
}

// ---------------------------------------------------------------------------
// Writing: always literal, always disclosed. POLICY — the spec ranks neither.
// ---------------------------------------------------------------------------

/// Output never carries a copy-range, whatever the input carried.
///
/// Pinned as policy. If ferro ever gains a "preserve the authored payload" mode
/// this test is the one to change, and changing it should be a decision rather
/// than a side effect.
#[test]
fn policy_output_is_always_literal_never_a_copy_range() {
    for body in [
        "302_303ins201_212",
        "302_303ins[201_212]",
        "302_313delins201_212",
        "302_303ins500_511",
    ] {
        let (got, _) = normalize_with_codes(body);
        assert!(
            !got.contains("ins201_") && !got.contains("ins500_"),
            "{body} kept a copy-range in its output: {got}"
        );
    }
}

/// The expansion is disclosed rather than silent.
///
/// The warning is the only signal that the emitted string is a different
/// *description* from the one authored, so a consumer diffing input against
/// output has something to key on.
#[test]
fn policy_expansion_is_always_disclosed() {
    for body in ["302_303ins201_212", "302_313delins201_212"] {
        let (_, codes) = normalize_with_codes(body);
        assert!(
            codes.iter().any(|c| c == "INSERTED_SEQUENCE_EXPANDED"),
            "{body} expanded without INSERTED_SEQUENCE_EXPANDED"
        );
    }
}

/// A copy-range and the literal it denotes converge.
///
/// This is the confluence property for gap 1: two legal spellings of one variant
/// must reach one normalized form, or downstream equality is spelling-sensitive.
#[test]
fn a_copy_range_and_its_literal_converge() {
    let from_range = normalize_with_codes("302_303ins201_212").0;
    let from_literal = normalize_with_codes(&format!("302_303ins{COPY_201_212}")).0;
    assert_eq!(from_range, from_literal);
}

/// A reversed copy-range is refused rather than silently reordered.
#[test]
fn a_reversed_copy_range_is_refused() {
    let input = format!("{LOCAL_CONTIG}:g.302_303ins212_201");
    let variant: HgvsVariant = parse_hgvs(&input).expect("parse");
    let normalizer = Normalizer::with_config(
        two_contig_provider(),
        NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
    );
    let err = normalizer
        .normalize(&variant)
        .expect_err("a reversed copy-range names no span");
    assert!(
        format!("{err}").contains("start <= end") || format!("{err}").contains("range"),
        "unexpected refusal for a reversed copy-range: {err}"
    );
}
