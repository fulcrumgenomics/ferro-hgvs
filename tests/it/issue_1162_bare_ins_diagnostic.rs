//! Issue #1162: a bare `ins` must be *diagnosed*, not merely rejected.
//!
//! #1137 made an insertion with no inserted sequence a parse error, but the
//! rejection lives in `parse_insertion` at the nom level, so the user-facing
//! message was the generic fall-through — byte-for-byte what unrecognized
//! garbage produces, leaking raw nom `Debug` output:
//!
//! ```text
//! Parse error at position 13: Failed to parse variant: Error(Error { input: "ins", code: Tag })
//! ```
//!
//! The class now has its own registry code (`W3027`), so the message names the
//! problem and `ferro explain W3027` documents it. The diagnosis is surfaced on
//! both entry points: the preprocessor phase (config-driven paths: CLI,
//! service) and a pre-parse check inside `parse_variant` (the bare `parse_hgvs`
//! path, which never reaches the preprocessor).
//!
//! These tests need no prepared reference, so they gate every PR.

use ferro_hgvs::error::{ErrorCode, FerroError};
use ferro_hgvs::error_handling::{get_code_info, ErrorConfig, InputPreprocessor, ModeAction};
use ferro_hgvs::hgvs::parser::{parse_hgvs_lenient, parse_hgvs_with_config};
use ferro_hgvs::parse_hgvs;

/// Bare-`ins` inputs across coordinate systems, each a two-position anchor with
/// the payload omitted. Mirrors `issue_1137_reject_bare_ins::BARE_INS` plus the
/// protein axis and the uncertain-range anchor from the bulk fixture.
const BARE_INS: &[&str] = &[
    "NC_000001.11:g.100_101ins",
    "NM_004006.3:c.179_180ins",
    "NM_000088.3:n.100_101ins",
    "NM_004006.2:r.100_101ins",
    "NC_012920.1:m.100_101ins",
    "NP_000079.2:p.Arg97_Trp98ins",
    // The shape recorded in the bulk failure-expectations fixture under
    // `[insertion_without_inserted_sequence]` — the anchor ends in `)`.
    "NM_001134407.3:c.(414+1_415-1)_(1007+1_1008-1)ins",
];

// ---------------------------------------------------------------------------
// The registry code itself.
// ---------------------------------------------------------------------------

/// `ferro explain W3027` must resolve, cite the spec, and declare that the
/// class is rejected in every mode (there is no derivable repair — the bases
/// are simply absent).
#[test]
fn w3027_is_registered_and_always_rejects() {
    let info = get_code_info("W3027").expect("W3027 must be registered for `ferro explain`");
    assert_eq!(info.name, "InsertionWithoutInsertedSequence");

    let behavior = info
        .mode_behavior
        .expect("a W-code must declare a mode behavior");
    for (mode, resolved) in [
        ("strict", behavior.strict),
        ("lenient", behavior.lenient),
        ("silent", behavior.silent),
    ] {
        assert_eq!(
            resolved,
            ModeAction::Reject,
            "W3027 must reject in {mode} mode (no repair is derivable); got {resolved:?}"
        );
    }

    let prose = format!("{} {}", info.summary, info.explanation);
    assert!(
        prose.contains("insertion.md"),
        "W3027 must cite the insertion recommendation; got: {prose}"
    );
    assert!(
        prose.contains("syntax.yaml"),
        "W3027 must cite the absence of a payload-less production in syntax.yaml; got: {prose}"
    );
    assert!(
        info.bad_examples.iter().any(|e| e.ends_with("ins")),
        "W3027 needs a bare-`ins` bad example; got: {:?}",
        info.bad_examples
    );
    assert!(
        info.good_examples.iter().any(|e| e.contains("insA")),
        "W3027 needs a payloaded good example; got: {:?}",
        info.good_examples
    );
}

// ---------------------------------------------------------------------------
// The bare `parse_hgvs` path (never reaches the preprocessor).
// ---------------------------------------------------------------------------

/// Every bare-`ins` shape must carry the dedicated W3027 code on the structured
/// `Diagnostic.code` channel, so downstream tooling can discriminate this
/// rejection from a generic parse error.
#[test]
fn bare_ins_carries_the_w3027_code() {
    for input in BARE_INS {
        let err = parse_hgvs(input).expect_err("a bare `ins` must be rejected");
        assert_eq!(
            err.code(),
            Some(ErrorCode::InsertionWithoutInsertedSequence),
            "{input:?} must reject with the dedicated W3027 code; got: {err:?}"
        );
        assert_eq!(err.code().map(|c| c.as_str()), Some("W3027".to_string()));
    }
}

/// The message must name the problem instead of dumping a nom `ErrorKind`.
#[test]
fn bare_ins_message_names_the_class() {
    for input in BARE_INS {
        let msg = parse_hgvs(input)
            .expect_err("a bare `ins` must be rejected")
            .to_string();
        assert!(
            msg.to_ascii_lowercase().contains("inserted sequence"),
            "{input:?} must self-identify as a missing inserted sequence; got: {msg}"
        );
        assert!(
            !msg.contains("ErrorKind") && !msg.contains("Error(Error {"),
            "{input:?} must not leak raw nom Debug output; got: {msg}"
        );
    }
}

/// The diagnostic must point at the spec-offered payload forms so the user can
/// repair the description.
#[test]
fn bare_ins_message_points_at_the_payload_forms() {
    let msg = parse_hgvs("NC_000001.11:g.100_101ins")
        .expect_err("a bare `ins` must be rejected")
        .to_string();
    assert!(
        msg.contains("insertion.md"),
        "the diagnostic should cite the spec section; got: {msg}"
    );
    assert!(
        msg.contains("insA") || msg.contains("insAGC"),
        "the diagnostic should show a literal-bases repair; got: {msg}"
    );
    assert!(
        msg.contains("858_895") || msg.to_ascii_lowercase().contains("range"),
        "the diagnostic should mention the same-reference range form; got: {msg}"
    );
}

/// The span must point at the offending `ins` token, not at position 0.
#[test]
fn bare_ins_span_covers_the_ins_token() {
    let input = "NC_000001.11:g.100_101ins";
    let err = parse_hgvs(input).expect_err("a bare `ins` must be rejected");
    let FerroError::Parse {
        diagnostic: Some(diagnostic),
        ..
    } = &err
    else {
        panic!("the W3027 rejection must be a diagnostic-carrying parse error; got: {err:?}");
    };
    let span = diagnostic
        .span
        .as_ref()
        .expect("the W3027 diagnostic must carry a span");
    assert_eq!(
        &input[span.start..span.end],
        "ins",
        "the span must cover the bare `ins` token; got {span:?}"
    );
}

/// Lenient and strict modes reject with the same diagnosis — there is no
/// derivable repair, so the mode makes no difference to the outcome.
#[test]
fn lenient_and_strict_reject_with_the_same_diagnosis() {
    for input in BARE_INS {
        for (mode, err) in [
            (
                "lenient",
                parse_hgvs_lenient(input).err().map(|e| e.to_string()),
            ),
            (
                "strict",
                parse_hgvs_with_config(input, ErrorConfig::strict())
                    .err()
                    .map(|e| e.to_string()),
            ),
        ] {
            let msg = err.unwrap_or_else(|| panic!("{input:?} must be rejected in {mode} mode"));
            assert!(
                msg.to_ascii_lowercase().contains("inserted sequence"),
                "{input:?} in {mode} mode must name the class; got: {msg}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// The preprocessor path (CLI / service).
// ---------------------------------------------------------------------------

/// The preprocessor detects the class up-front in every mode, so config-driven
/// callers surface the diagnosis before the parser is reached.
#[test]
fn preprocessor_rejects_bare_ins_in_every_mode() {
    for input in BARE_INS {
        for (mode, preprocessor) in [
            ("strict", InputPreprocessor::strict()),
            ("lenient", InputPreprocessor::lenient()),
            ("silent", InputPreprocessor::silent()),
        ] {
            let result = preprocessor.preprocess(input);
            assert!(
                !result.success,
                "{input:?} must be rejected by the {mode} preprocessor"
            );
            let err = result
                .error
                .unwrap_or_else(|| panic!("{input:?}: a failed result must carry an error"));
            assert_eq!(
                err.code(),
                Some(ErrorCode::InsertionWithoutInsertedSequence),
                "{input:?} in {mode} mode must carry W3027; got: {err:?}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Regression guards: the neighbouring classes keep their own diagnoses.
// ---------------------------------------------------------------------------

/// Payloaded insertions are untouched — the detector keys on the *absence* of
/// an inserted sequence, so every legitimate payload form must still parse.
#[test]
fn payloaded_insertions_still_parse() {
    for input in [
        "NC_000001.11:g.100_101insA",
        "NC_000001.11:g.100_101insATG",
        "NC_000001.11:g.100_101ins10",
        "NC_000001.11:g.100_101ins[A[10];T]",
        "NM_000088.3:c.100_101insA",
        "NC_000001.11:g.849_850ins858_895",
        "NM_000088.3:c.100_101ins[NM_004006.2:c.1_10]",
        "NM_000088.3:c.100_101insN[10]",
        // A cis allele whose members all carry payloads.
        "NM_000088.3:c.[100_101insA;200_201insT]",
    ] {
        assert!(
            parse_hgvs(input).is_ok(),
            "a payloaded insertion must still parse: {input:?}"
        );
    }
}

/// A bare `delins` is a distinct class (it collapses to a plain `del`, W3012)
/// and must not be swept up by the W3027 detector. Likewise `dupins`, which
/// keeps its own spec-cited rejection.
#[test]
fn delins_and_dupins_keep_their_own_diagnoses() {
    // `delins` with no insert: still parses (standard) / rewrites to `del`
    // (lenient) — exactly as before W3027 existed.
    let standard = parse_hgvs("NC_000001.11:g.100_101delins")
        .expect("a bare `delins` must still parse in standard mode");
    assert_eq!(standard.to_string(), "NC_000001.11:g.100_101delins");
    let lenient = parse_hgvs_lenient("NC_000001.11:g.100_101delins")
        .expect("a bare `delins` must still parse in lenient mode");
    assert_eq!(lenient.result.to_string(), "NC_000001.11:g.100_101del");

    // `dupins` keeps its `duplication.md:92` diagnosis.
    let msg = parse_hgvs("NC_000001.11:g.100_101dupinsA")
        .expect_err("`dupins` is rejected")
        .to_string();
    assert!(
        msg.contains("dupins"),
        "`dupins` must keep its own diagnosis; got: {msg}"
    );
}

/// Unrecognized garbage still gets the generic fall-through — W3027 must not
/// widen into a catch-all for "the edit did not parse".
#[test]
fn unrecognized_edits_are_not_claimed_by_w3027() {
    for input in [
        "NC_000001.11:g.100_101zzz",
        "NC_000001.11:g.100_101in",
        "NM_000088.3:c.100_101sni",
    ] {
        let err = parse_hgvs(input).expect_err("garbage must be rejected");
        assert_ne!(
            err.code(),
            Some(ErrorCode::InsertionWithoutInsertedSequence),
            "{input:?} is not a bare-`ins` case; got: {err:?}"
        );
    }
}

/// A gene symbol spelled `INS` (insulin) sits next to the very token the
/// detector scans for. Case matters: the edit keyword is lowercase, the symbol
/// is not, and the symbol is never preceded by a position.
#[test]
fn ins_gene_symbol_is_not_mistaken_for_a_bare_edit() {
    for input in [
        "NM_000207.3(INS):c.100A>G",
        "NC_000011.10(INS):g.2159830A>G",
        "NM_000207.3(INS):c.100_101insA",
    ] {
        assert!(
            parse_hgvs(input).is_ok(),
            "the INS gene symbol must not trip the bare-`ins` detector: {input:?}"
        );
    }
}
