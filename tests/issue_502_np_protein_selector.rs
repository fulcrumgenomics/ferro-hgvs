//! Issue #502: emit the `NP_` protein-accession selector for `p.` variants on
//! a genomic reference.
//!
//! For `NG_/LRG_/NC_(NM_):p.` inputs the parenthetical selector is a coding
//! transcript, but a `p.` coordinate is a protein-level statement — so the
//! spec-preferred selector is the paired protein accession (`NP_`), per HGVS
//! `background/refseq.md` L38-42. `normalize` resolves the transcript's
//! `protein_id` and rewrites `NM_`→`NP_` (keeping the genomic-context wrapper),
//! while leaving `c./n./r.` selectors and bare/already-`NP_` inputs untouched.

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{ManeStatus, Strand};
use ferro_hgvs::reference::Transcript;
use ferro_hgvs::{parse_hgvs, Normalizer};

/// A `MockProvider` carrying the `NM_003002.4` transcript (SDHD) used by the
/// #502 corpus row, optionally with its paired protein accession populated.
/// The sequence is an all-`A` filler long enough to validate the `c.` coordinate
/// used below; protein bases are intentionally absent so `has_protein_data()` is
/// false and protein-reference validation is skipped (the selector rewrite is
/// independent of it).
fn provider_with(protein_id: Option<&str>) -> MockProvider {
    let mut provider = MockProvider::new();
    let tx = Transcript::new(
        "NM_003002.4".to_string(),
        Some("SDHD".to_string()),
        Strand::Plus,
        Some("A".repeat(300)),
        Some(1),
        Some(300),
        Vec::new(),
        None,
        None,
        None,
        Default::default(),
        ManeStatus::Select,
        None,
        None,
    )
    .with_protein_id(protein_id.map(str::to_string));
    provider.add_transcript(tx);
    provider
}

fn normalize_str(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
    let normalized = normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize {input}: {e}"));
    format!("{normalized}")
}

#[test]
fn ng_nm_protein_input_emits_np_selector() {
    let out = provider_with(Some("NP_002993.1"));
    assert_eq!(
        normalize_str(out, "NG_012337.3(NM_003002.4):p.(Asp92Tyr)"),
        "NG_012337.3(NP_002993.1):p.(Asp92Tyr)"
    );
}

#[test]
fn np_selector_rewrite_is_idempotent() {
    // Normalizing the rewritten output again must be a no-op: the selector is
    // already `NP_`, so there is nothing to resolve.
    let once = normalize_str(
        provider_with(Some("NP_002993.1")),
        "NG_012337.3(NM_003002.4):p.(Asp92Tyr)",
    );
    let twice = normalize_str(provider_with(Some("NP_002993.1")), &once);
    assert_eq!(once, twice);
    assert_eq!(twice, "NG_012337.3(NP_002993.1):p.(Asp92Tyr)");
}

#[test]
fn already_np_selector_is_preserved() {
    // An input that already carries the `NP_` selector is left unchanged (the
    // prefix is not a coding transcript, so no resolution is attempted).
    assert_eq!(
        normalize_str(
            provider_with(Some("NP_002993.1")),
            "NG_012337.3(NP_002993.1):p.(Asp92Tyr)"
        ),
        "NG_012337.3(NP_002993.1):p.(Asp92Tyr)"
    );
}

#[test]
fn bare_nm_protein_without_genomic_context_is_preserved() {
    // No genomic-context wrapper → out of scope; the input selector is
    // preserved verbatim (#121).
    assert_eq!(
        normalize_str(
            provider_with(Some("NP_002993.1")),
            "NM_003002.4:p.(Asp92Tyr)"
        ),
        "NM_003002.4:p.(Asp92Tyr)"
    );
}

#[test]
fn missing_protein_id_preserves_nm_selector() {
    // The reference resolves the transcript but it carries no `protein_id`:
    // fall back to preserving the input `NM_` selector rather than erroring.
    assert_eq!(
        normalize_str(provider_with(None), "NG_012337.3(NM_003002.4):p.(Asp92Tyr)"),
        "NG_012337.3(NM_003002.4):p.(Asp92Tyr)"
    );
}

#[test]
fn unknown_transcript_preserves_nm_selector() {
    // The provider has no such transcript (empty provider): preserve the input
    // selector; never error on the missing-data path.
    let out = normalize_str(MockProvider::new(), "NG_012337.3(NM_003002.4):p.(Asp92Tyr)");
    assert_eq!(out, "NG_012337.3(NM_003002.4):p.(Asp92Tyr)");
}

#[test]
fn non_protein_protein_id_preserves_nm_selector() {
    // The provider carries a `protein_id` that parses but is not a protein
    // accession (e.g. a transcript accession from malformed reference data).
    // Emitting it as a `p.` selector would violate the NP_/XP_ contract, so the
    // input `NM_` selector is preserved instead.
    assert_eq!(
        normalize_str(
            provider_with(Some("NM_999999.9")),
            "NG_012337.3(NM_003002.4):p.(Asp92Tyr)"
        ),
        "NG_012337.3(NM_003002.4):p.(Asp92Tyr)"
    );
}

#[test]
fn coding_selector_on_genomic_reference_is_not_rewritten() {
    // The rewrite is `p.`-only: a `c.` coordinate on the same `NG_(NM_)`
    // reference must keep its transcript selector (must NOT become `NP_`).
    let out = normalize_str(
        provider_with(Some("NP_002993.1")),
        "NG_012337.3(NM_003002.4):c.92A>T",
    );
    assert!(
        out.contains("(NM_003002.4)"),
        "c. selector must stay NM_: {out}"
    );
    assert!(
        !out.contains("NP_"),
        "c. selector must not become NP_: {out}"
    );
}

#[test]
fn predicted_xm_transcript_selector_emits_xp() {
    // The rewrite covers predicted models too: an `XM_` coding-transcript
    // selector resolves to its `XP_` protein accession.
    let mut provider = MockProvider::new();
    let tx = Transcript::new(
        "XM_011512873.2".to_string(),
        Some("SDHD".to_string()),
        Strand::Plus,
        Some("A".repeat(300)),
        Some(1),
        Some(300),
        Vec::new(),
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    )
    .with_protein_id(Some("XP_011511175.1".to_string()));
    provider.add_transcript(tx);
    assert_eq!(
        normalize_str(provider, "NG_012337.3(XM_011512873.2):p.(Asp92Tyr)"),
        "NG_012337.3(XP_011511175.1):p.(Asp92Tyr)"
    );
}
