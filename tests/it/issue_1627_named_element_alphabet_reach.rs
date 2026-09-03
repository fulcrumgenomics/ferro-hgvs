//! Issue #1627 / #1630 — how far the named-element fallback reaches, measured.
//!
//! Sources: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/1627>,
//! <https://github.com/fulcrumgenomics/ferro-hgvs/issues/1630>.
//!
//! # This file pins the GRAMMAR, which #1627 deliberately did not change
//!
//! What ferro should *do* about `standards.md:39` is settled — the decided
//! `rulings[alignment-only-symbol-in-a-description]` says neither `X` nor `-`
//! may appear in a description and ferro must refuse both. **Where** it refuses
//! is the decided `rulings[absolute-prohibition-enforcement-stage]`: strict
//! fails at PARSE, lenient does not validate input conformance and fails only
//! when it cannot NORMALIZE, silent is lenient without messages.
//!
//! **That two-stage schedule describes `X` only.** `-` never reaches the AST,
//! so no conformance rule is consulted for it at either stage: the grammar
//! refuses it in every mode, for grammar reasons, which is measured below in
//! `the_dash_is_refused_by_the_grammar_rather_than_by_a_rule`. Reading the
//! schedule as covering both symbols is the misattribution this note exists to
//! prevent — it would make `-` look mode-dependent when it is not.
//!
//! **So a grammar-level refusal would be the wrong fix — it would refuse in
//! lenient too — and #1627 does not make one.** The alphabet of
//! `parse_inserted_sequence` is untouched; the refusal is a post-parse rule on
//! the AST (`src/hgvs/alignment_symbols.rs`, applied in
//! `parse_hgvs_with_config` and in `Normalizer::normalize_core_checked`).
//!
//! Every assertion below therefore still describes ferro today, and that is the
//! point of keeping the file: it reads through the bare `parse_hgvs`, which
//! applies no `ErrorConfig` (#1632), so it observes the grammar alone. If a
//! later change moves the refusal into the grammar, this file goes red and says
//! why that is wrong. The three-mode, two-stage behaviour is asserted
//! separately, in `corpus_prohibited_inputs::`
//! `an_alignment_only_symbol_is_refused_in_every_mode_for_both_x_and_dash`.
//!
//! What it measured, and what that measurement was for: the surface the fix had
//! to cover, because the record understated it.
//!
//! # The correction
//!
//! `corpus_alphabet_and_intronic_validity.rs`'s header analysis says `X`
//! survives via the named-element arm because "a lone uppercase letter with
//! nothing alphanumeric following matches it", and adds "this is not
//! `X`-specific — any lone uppercase letter that is not otherwise consumed hits
//! the same arm". Both sentences are true and both are about **one** of the two
//! arms that produce `InsertedSequence::Named`.
//!
//! `parse_inserted_sequence` (`src/hgvs/parser/edit.rs`) dispatches on the
//! *first* byte. A non-IUPAC first byte (`X…`) takes the trailing
//! `c if c.is_ascii_uppercase()` arm, which is the one the record describes. But
//! an **IUPAC** first byte takes the `is_iupac_base` arm, and that arm walks the
//! whole run setting `has_non_iupac` on any uppercase byte that is not an IUPAC
//! base — so a sequence that is literal apart from one stray symbol is
//! reclassified wholesale as a mobile-element name.
//!
//! `NC_TEST.1:g.10delinsACGTX` is therefore `Named("ACGTX")`, not
//! `Literal(ACGT)` plus a rejected `X`, and nothing in the corpus or in the
//! committed analysis covers that shape. It mattered for the fix in two ways:
//! the prohibited-row counts are drawn from lone-`X` inputs only, and a check
//! keyed on "the insert is exactly one non-IUPAC letter" would close the narrow
//! case and leave this one open. `alignment_only_symbol` accordingly tests for
//! the *presence* of a daggered symbol anywhere in the name.
//!
//! **The corpus figure is a STRUCTURAL zero, not a clean bill of health.** All
//! 24 `standards.md:39` rows spell the insert `delinsX`; re-measured at
//! `78c43230`, the embedded count over the corpus is **0**, because the
//! generator never varies the insert's literal run. Quoting that zero as
//! "embedded shapes do not occur" would be exactly the mistake `CONTRIBUTING.md`
//! warns about, so the embedded cases are asserted here by construction.

use ferro_hgvs::hgvs::edit::{InsertedSequence, NaEdit};
use ferro_hgvs::hgvs::variant::HgvsVariant;
use ferro_hgvs::parse_hgvs;

/// The `InsertedSequence` a single genomic `delins` carries, or `None` when the
/// input does not parse.
fn inserted(input: &str) -> Option<InsertedSequence> {
    let variant = parse_hgvs(input).ok()?;
    let HgvsVariant::Genome(g) = variant else {
        panic!("{input}: expected a genomic variant");
    };
    match g.loc_edit.edit.inner() {
        Some(NaEdit::Delins { sequence, .. }) => Some(sequence.clone()),
        other => panic!("{input}: expected a delins, got {other:?}"),
    }
}

/// True when the insert was classified as a named mobile element.
fn is_named(input: &str) -> bool {
    matches!(inserted(input), Some(InsertedSequence::Named(_)))
}

/// The shape the record already describes: a lone non-IUPAC uppercase letter is
/// taken as a mobile-element name. Pinned here as the baseline the correction
/// below is measured against, and to show it is not `X`-specific.
#[test]
fn a_lone_non_iupac_uppercase_letter_is_taken_as_a_named_element() {
    for input in [
        "NC_TEST.1:g.10delinsX",
        "NC_TEST.1:g.10delinsZ",
        "NC_TEST.1:g.10delinsQ",
        "NC_TEST.1:g.10delinsE",
    ] {
        assert!(is_named(input), "{input} should currently parse as Named");
    }
}

/// **The correction.** An otherwise-literal IUPAC run carrying one stray
/// non-IUPAC uppercase byte is reclassified *in its entirety* as a named
/// element, whichever end the stray sits at. This path is the `is_iupac_base`
/// arm, not the `c.is_ascii_uppercase()` arm the committed analysis describes,
/// and no corpus row exercises it.
#[test]
fn a_literal_run_with_one_stray_symbol_is_reclassified_whole() {
    for input in [
        "NC_TEST.1:g.10delinsACGTX",  // stray at the end
        "NC_TEST.1:g.10delinsXACGT",  // stray at the start
        "NC_TEST.1:g.10delinsACXGT",  // stray in the middle
        "NC_TEST.1:g.10delinsACGTZQ", // two strays
    ] {
        assert!(
            is_named(input),
            "{input} is currently Named, not Literal — this is the unrecorded half"
        );
    }
}

/// The genuine mobile-element spellings the arm exists for. They must keep
/// parsing as `Named` through whatever #1630 does, so they are pinned here as
/// the constraint on any narrowing of the arm — the fix cannot simply be "an
/// insert containing a non-IUPAC uppercase letter is refused", because `AluYb8`
/// and `LINE1` both contain several.
#[test]
fn real_mobile_element_names_are_named_and_must_stay_so() {
    for input in [
        "NC_TEST.1:g.10delinsAluYb8",
        "NC_TEST.1:g.10delinsLINE1",
        "NC_TEST.1:g.10delinsL1",
        "NC_TEST.1:g.10delinsAlu",
    ] {
        assert!(is_named(input), "{input} is a legitimate named element");
    }
}

/// A pure IUPAC run is `Literal`, so the reclassification above really is
/// triggered by the stray byte and not by the arm always winning.
#[test]
fn a_pure_iupac_run_is_still_literal() {
    assert!(matches!(
        inserted("NC_TEST.1:g.10delinsACGT"),
        Some(InsertedSequence::Literal(_))
    ));
}

/// `standards.md:39`'s other daggered symbol, for contrast — and the reason the
/// two halves of one footnote needed different work.
///
/// `-` matches no arm at all, so it is left unconsumed and the variant-level
/// trailing-character check refuses it. That is a *grammar* refusal, not a
/// conformance rule, which is why it is refused in every mode and why it needed
/// no rule of its own. It is measured in every position — leading, trailing and
/// interior, `ins` and `delins` — because "refused" for the lone symbol would
/// not have entailed "refused when embedded", which is precisely the asymmetry
/// `X` turned out to have.
#[test]
fn the_dash_is_refused_by_the_grammar_rather_than_by_a_rule() {
    for input in [
        "NC_TEST.1:g.10delins-",
        "NC_TEST.1:g.10delins-ACGT",
        "NC_TEST.1:g.10delinsACGT-",
        "NC_TEST.1:g.10delinsAC-GT",
        "NC_TEST.1:g.10_11ins-",
    ] {
        let err = parse_hgvs(input).unwrap_err_or_else_msg(input);
        assert!(
            err.contains("trailing") || err.contains("Parse error"),
            "{input}: the refusal is a grammar one, not an alphabet finding; got: {err}"
        );
    }
}

/// Helper: the error message, or a panic naming the input that was accepted.
trait UnwrapErrMsg {
    fn unwrap_err_or_else_msg(self, input: &str) -> String;
}

impl<T: std::fmt::Debug> UnwrapErrMsg for Result<T, ferro_hgvs::error::FerroError> {
    fn unwrap_err_or_else_msg(self, input: &str) -> String {
        match self {
            Ok(v) => panic!("{input} must be refused by the grammar, but parsed to {v:?}"),
            Err(e) => e.to_string(),
        }
    }
}
