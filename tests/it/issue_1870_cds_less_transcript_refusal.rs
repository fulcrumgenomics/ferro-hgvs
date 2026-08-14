//! Issue #1870 — a `c.` description against a transcript whose CDS the
//! reference cannot resolve is REFUSED, not normalized silently.
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/1870>.
//! Ruling: `rulings[c-description-against-an-unresolvable-cds-is-refused]`
//! in `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`.
//!
//! # What this replaces
//!
//! `normalize` used to return the input verbatim with `status: ok`,
//! `changed: false` and an EMPTY warning vector, byte-identically in strict,
//! lenient and silent. On such an accession `normalize` was the **identity
//! function**, not a canonicalising one — measured on `NM_000546.3` (TP53), all
//! five of `c.528del`…`c.532del` (one deletion from one `CCCCC` run, one
//! variant) came back as five distinct strings, every one `changed=false` /
//! `ok`. A caller deduplicating by normalized string got five variants where
//! there is one. That is the shape `the_five_spellings_of_one_deletion_no_longer_diverge`
//! pins, on a hermetic fixture with the same geometry.
//!
//! # The three things the ruling settles, and which row pins each
//!
//! | ruling | guard |
//! |---|---|
//! | refuse, rather than answer with the input | [`the_refusal_is_unconditional_across_every_mode`] |
//! | at NORMALIZE, not at parse — the input is conformant, the *reference* is not | [`the_input_still_parses_because_the_spelling_is_conformant`] |
//! | the `n.` axis on the same record is untouched | [`the_n_axis_still_shifts_on_the_same_record`] |
//!
//! # Why NORMALIZE and not PARSE
//!
//! `absolute-prohibition-enforcement-stage` rules that **input conformance** is
//! answered at parse, before the input is accepted. This is not input
//! conformance: `NM_NOCDS.1:c.10del` is a conformant description — correct
//! prefix, correct axis for a coding accession, a position and an edit the
//! grammar admits — and the identical string against a reference that carries
//! the CDS normalizes correctly. What fails is the *reference*. A fact that
//! changes with the provider cannot be settled at parse, and `parse_hgvs` holds
//! no provider to see it with, which is decisive.
//!
//! The contrast with the existing prefix rule is the point:
//! `src/hgvs/parser/variant.rs` already refuses `c.` on an `NR_`/`XR_`
//! accession "which has no CDS" — at parse, in every mode, correctly, because
//! the accession's own *prefix* carries the fact and no reference is needed to
//! read it. This ruling makes the refusal follow the fact rather than the
//! naming convention; it does not move the prefix rule, which is pinned here by
//! [`the_prefix_rule_is_untouched_and_still_refuses_at_parse`].
//!
//! # Why unconditional across strict / lenient / silent
//!
//! Not a new mode policy — an application of the existing one.
//! `absolute-prohibition-enforcement-stage` states lenient's contract as "input
//! conformance is NOT validated. It fails only if it cannot NORMALIZE, which
//! can happen", and silent's as "lenient with no error messages, but an exit
//! code". Ferro genuinely cannot normalize here: with nothing locating the
//! `ATG`, there is no origin to number from, so the 3'-shift has no window and
//! the W4004 bounds gate has no bound. Lenient's own stated behaviour is
//! therefore to fail.
//!
//! # The error is `project`'s, deliberately
//!
//! `ferro project 'NM_000546.3:c.528del' --axis c` already exited 1 with
//! `Coordinate conversion error: transcript NM_000546.3 has no CDS start` while
//! `ferro normalize` on the identical input returned `status: ok` and exit 0.
//! Reusing `FerroError::ConversionError` and that message text is what stops
//! the two surfaces drifting apart again, and is asserted rather than assumed
//! by [`the_refusal_reuses_projections_own_error_and_words`].

use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, FerroError, MockProvider, NormalizeConfig, Normalizer};

/// The accession whose provider record carries **no** CDS bounds.
///
/// An `NM_` on purpose: a `c.` description against an `NR_`/`XR_` is refused at
/// parse by the prefix rule, so an `NM_` is the only way to reach the
/// reference-derived path this file is about.
const NO_CDS: &str = "NM_NOCDS.1";

/// The control: the identical sequence, served with a resolvable CDS.
const WITH_CDS: &str = "NM_HASCDS.1";

/// A bare **genomic parent**, served with no CDS — as every `NG_` is, by
/// construction. Out of the ruling's scope; see
/// [`a_bare_genomic_parent_is_out_of_scope_and_still_answers`].
const GENOMIC_PARENT: &str = "NG_000001.1";

/// 5'UTR (10 nt) + a CDS whose bases mirror the TP53 geometry the issue
/// measured: a `CCCCC` run five bases into the coding sequence, so a `del`
/// anywhere in it denotes one variant with one canonical spelling.
///
/// `ATGAA CCCCC ...` — `c.1` is the `A` of the `ATG`, so the run is `c.6`–`c.10`
/// and the canonical (3'-most) spelling of deleting one `C` is `c.10del`.
const SEQUENCE: &str = "GGTTAAGGTTATGAACCCCCGGTTAAGGTTAA";

/// CDS start in **transcript** coordinates, 1-based: the `A` of the `ATG`.
const CDS_START: u64 = 11;
/// CDS end in transcript coordinates, 1-based.
const CDS_END: u64 = 22;

/// Build a provider serving both records: `NO_CDS` with
/// `cds_start`/`cds_end` absent, `WITH_CDS` with both resolved. Same bases,
/// same exon table, same strand — the *only* difference between the two is the
/// CDS annotation, which is what makes the pair a controlled comparison rather
/// than two unrelated fixtures.
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let len = SEQUENCE.len() as u64;
    for (id, cds) in [
        (NO_CDS, None),
        (WITH_CDS, Some((CDS_START, CDS_END))),
        (GENOMIC_PARENT, None),
    ] {
        provider.add_transcript(Transcript::new(
            id.to_string(),
            Some("TEST".to_string()),
            Strand::Plus,
            SEQUENCE.to_string(),
            cds.map(|(s, _)| s),
            cds.map(|(_, e)| e),
            vec![Exon::new(1, 1, len)],
            None,
            None,
            None,
            Default::default(),
            ManeStatus::None,
            None,
            None,
        ));
    }
    provider
}

fn config_for(mode: ErrorMode) -> NormalizeConfig {
    match mode {
        ErrorMode::Strict => NormalizeConfig::strict(),
        ErrorMode::Lenient => NormalizeConfig::lenient(),
        ErrorMode::Silent => NormalizeConfig::silent(),
    }
}

fn normalize_in(mode: ErrorMode, input: &str) -> Result<String, FerroError> {
    let normalizer = Normalizer::with_config(provider(), config_for(mode));
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} must parse: {e}"));
    normalizer.normalize(&variant).map(|v| v.to_string())
}

/// THE RULING. A `c.` description against a record whose CDS does not resolve
/// is refused — in strict, in lenient and in silent alike.
///
/// The three-mode loop is the assertion, not decoration: the previous behaviour
/// was byte-identical across all three (`warnings: []` in every one), so a guard
/// that checked only strict would pass on a mode-gated implementation that left
/// two thirds of the finding in place.
#[test]
fn the_refusal_is_unconditional_across_every_mode() {
    for mode in [ErrorMode::Strict, ErrorMode::Lenient, ErrorMode::Silent] {
        let outcome = normalize_in(mode, &format!("{NO_CDS}:c.6del"));
        let err = match outcome {
            Err(e) => e,
            Ok(out) => {
                panic!("{mode:?} must refuse a c. description on a CDS-less record, got: {out}")
            }
        };
        let msg = err.to_string();
        assert!(
            msg.contains(NO_CDS) && msg.contains("has no CDS start"),
            "{mode:?}: the refusal must name the accession and the missing CDS start, got: {msg}"
        );
    }
}

/// The refusal reuses the error `project` already returns for this exact fact,
/// so `normalize` and `project` cannot disagree about one input again.
#[test]
fn the_refusal_reuses_projections_own_error_and_words() {
    let err = normalize_in(ErrorMode::Strict, &format!("{NO_CDS}:c.6del"))
        .expect_err("a CDS-less record must refuse");
    assert!(
        matches!(err, FerroError::ConversionError { .. }),
        "expected FerroError::ConversionError (E5001), the variant `project` already returns \
         for `has no CDS start`; got: {err:?}"
    );
    assert_eq!(
        err.code(),
        Some(ferro_hgvs::error::ErrorCode::ConversionFailed),
        "the code must match `project`'s"
    );
}

/// THE SECOND PUBLIC EXIT. `Normalizer` has exactly two public normalizing
/// methods — `normalize()` and `normalize_with_diagnostics()` — and every other
/// row in this file goes through the first.
///
/// That is not a stylistic gap. The defect this file replaces was *reported* in
/// the second one's vocabulary: `status: ok`, `changed: false` and an empty
/// warning vector are fields of [`ferro_hgvs::normalize::NormalizeResult`] and of
/// the CLI rendering built on it, not of anything `normalize()` returns. So
/// pinning only `normalize()` would leave the surface the finding was actually
/// observed on untested.
///
/// Both exits do reach the check today — it sits in `normalize_cds`, dispatched
/// from `normalize_core`, which `normalize_with_diagnostics` routes through via
/// `normalize_core_checked` (#1382). But that is an implementation fact, and it
/// has already been false once: before #1382 this exit called
/// `normalize_core_canonical` directly and applied none of the rejection ladder.
/// A guard is what stops it becoming false again silently.
#[test]
fn the_refusal_reaches_the_diagnostics_exit_too() {
    for mode in [ErrorMode::Strict, ErrorMode::Lenient, ErrorMode::Silent] {
        let normalizer = Normalizer::with_config(provider(), config_for(mode));
        let input = format!("{NO_CDS}:c.6del");
        let variant = parse_hgvs(&input).unwrap_or_else(|e| panic!("{input} must parse: {e}"));
        let err = match normalizer.normalize_with_diagnostics(&variant) {
            Err(e) => e,
            Ok(r) => panic!(
                "{mode:?}: normalize_with_diagnostics must refuse a c. description on a \
                 CDS-less record, got {} with warnings {:?} and infos {:?} — which is \
                 exactly the `status: ok`, empty-warnings shape #1870 reported",
                r.result, r.warnings, r.infos
            ),
        };
        assert!(
            matches!(err, FerroError::ConversionError { .. }),
            "{mode:?}: the diagnostics exit must refuse with the same ConversionError as \
             `normalize()`, got: {err:?}"
        );
        assert!(
            err.to_string().contains(NO_CDS) && err.to_string().contains("has no CDS start"),
            "{mode:?}: the refusal must name the accession and the missing CDS start, got: {err}"
        );
    }
}

/// A MULTI-MEMBER CIS ALLELE on the CDS-less record is refused as a whole.
///
/// This is the shape the change's only real-corpus measurement is made of:
/// `multi_member_cis_axis`'s census moves `declined: 0 -> 6`, and all six rows
/// are `c.[…;…]` alleles. That census is reference-aware — it resolves a prepared
/// reference from `FERRO_MANIFEST` — and PR CI does not provision one, so it
/// skips green there and the allele path would otherwise ship with no guard that
/// PR CI actually runs. This row is that guard, hermetically.
///
/// `normalize_allele` recurses through `normalize_core` per member, so refusing
/// one member refuses the allele; asserting it here pins that composition rather
/// than leaving it inferred from the single-member rows above.
#[test]
fn a_multi_member_allele_on_the_same_record_is_refused_whole() {
    for mode in [ErrorMode::Strict, ErrorMode::Lenient, ErrorMode::Silent] {
        let err = match normalize_in(mode, &format!("{NO_CDS}:c.[6del;10del]")) {
            Err(e) => e,
            Ok(out) => panic!(
                "{mode:?}: a cis allele whose members are all c. descriptions on a CDS-less \
                 record must be refused, got: {out}"
            ),
        };
        assert!(
            err.to_string().contains("has no CDS start"),
            "{mode:?}: unexpected refusal reason for the allele: {err}"
        );
    }

    // The control, so the row above cannot pass by the allele grammar being
    // rejected outright: the identical allele on the record whose CDS resolves
    // still normalizes.
    let out = normalize_in(ErrorMode::Strict, &format!("{WITH_CDS}:c.[6del;10del]"))
        .expect("the same allele on the CDS-resolving control must normalize");
    assert!(
        out.starts_with(WITH_CDS),
        "the control allele must normalize against its own accession, got: {out}"
    );
}

/// The control, and the thing that makes every row above non-vacuous: the
/// SAME bases with a resolvable CDS still normalize, and still shift.
///
/// Without this row the whole file would pass on an implementation that refused
/// every `c.` description outright.
#[test]
fn a_record_whose_cds_resolves_is_untouched_and_still_shifts() {
    for mode in [ErrorMode::Strict, ErrorMode::Lenient, ErrorMode::Silent] {
        let out = normalize_in(mode, &format!("{WITH_CDS}:c.6del"))
            .unwrap_or_else(|e| panic!("{mode:?}: the control must normalize, got: {e:?}"));
        assert_eq!(
            out,
            format!("{WITH_CDS}:c.10del"),
            "{mode:?}: the control must 3'-shift across the CCCCC run"
        );
    }
}

/// The consumer-facing fact, stated as a property rather than as one row: on a
/// CDS-less record `normalize` used to be the IDENTITY, so every spelling of one
/// deletion survived as its own "normalized" string and a caller deduplicating
/// by normalized string got five variants where there is one.
///
/// Both halves are asserted. On the control the five spellings collapse to one;
/// on the CDS-less record none of them is answered at all. What must never
/// again be true is the third state — five distinct answers.
#[test]
fn the_five_spellings_of_one_deletion_no_longer_diverge() {
    let spellings = ["c.6del", "c.7del", "c.8del", "c.9del", "c.10del"];

    let collapsed: std::collections::BTreeSet<String> = spellings
        .iter()
        .map(|s| {
            normalize_in(ErrorMode::Strict, &format!("{WITH_CDS}:{s}"))
                .unwrap_or_else(|e| panic!("{s} on the control: {e:?}"))
        })
        .collect();
    assert_eq!(
        collapsed,
        std::collections::BTreeSet::from([format!("{WITH_CDS}:c.10del")]),
        "five spellings of one deletion must collapse to one normalized string when the CDS \
         resolves; got {collapsed:?}"
    );

    for s in spellings {
        let err = match normalize_in(ErrorMode::Strict, &format!("{NO_CDS}:{s}")) {
            Err(e) => e,
            Ok(out) => panic!("{s} on the CDS-less record must be refused, got: {out}"),
        };
        assert!(
            err.to_string().contains("has no CDS start"),
            "{s}: unexpected refusal reason: {err}"
        );
    }
}

/// The `n.` axis on the SAME record is untouched, and still shifts.
///
/// This is a requirement of the ruling rather than a side effect:
/// `background/numbering.md:52` numbers the non-coding axis "from the first to
/// the last nucleotide of the reference sequence", which has no CDS term in it,
/// so a missing CDS annotation cannot reach it. The bases are present and
/// correct on these records — only the annotation is missing — which is exactly
/// why the two axes must part company here.
///
/// `n.16`–`n.20` is the same `CCCCC` run as `c.6`–`c.10` (`CDS_START` is 11), so
/// this is the identical edit to [`the_refusal_is_unconditional_across_every_mode`]'s,
/// spelled on the axis that survives.
#[test]
fn the_n_axis_still_shifts_on_the_same_record() {
    for mode in [ErrorMode::Strict, ErrorMode::Lenient, ErrorMode::Silent] {
        let out = normalize_in(mode, &format!("{NO_CDS}:n.16del"))
            .unwrap_or_else(|e| panic!("{mode:?}: the n. axis must keep working, got: {e:?}"));
        assert_eq!(
            out,
            format!("{NO_CDS}:n.20del"),
            "{mode:?}: the n. axis must still 3'-shift across the CCCCC run"
        );
    }
}

/// The stage. The input's SPELLING is conformant, so it parses — in every mode.
/// The refusal happens later, when a provider exists to be asked.
///
/// This is what separates the ruling from `absolute-prohibition-enforcement-stage`,
/// whose subject is a spelling the spec prohibits and which is therefore
/// answered at parse. Nothing about `NM_NOCDS.1:c.6del` is wrong as a
/// string.
#[test]
fn the_input_still_parses_because_the_spelling_is_conformant() {
    let input = format!("{NO_CDS}:c.6del");
    let parsed = parse_hgvs(&input).expect("a conformant c. spelling must still parse");
    assert_eq!(
        parsed.to_string(),
        input,
        "parse must round-trip the description untouched — the refusal is not the parser's"
    );
}

/// The scope boundary. A bare genomic parent (`NG_`/`NC_`/`LRG_` with nothing in
/// its transcript-selector slot) is **not** a transcript, so the ruling does not
/// reach it and it still answers.
///
/// This is not defensive trimming — it was measured. An unscoped check looks up
/// the `NG_` record itself, which carries no CDS *by construction* (a RefSeqGene
/// genomic record is not a transcript and has no CDS annotation of its own), and
/// so refuses a `c.` description whose implied transcript the legacy-selector
/// machinery (#500 / #793 / #923) is what resolves. Unscoped, that newly failed
/// `NG_012337.1(10683):c.274G>T` in the mutalyzer comparator — a row whose
/// accepted value is ferro's own pass-through — and declined six real-corpus cis
/// rows for the same reason.
#[test]
fn a_bare_genomic_parent_is_out_of_scope_and_still_answers() {
    for mode in [ErrorMode::Strict, ErrorMode::Lenient, ErrorMode::Silent] {
        let out = normalize_in(mode, &format!("{GENOMIC_PARENT}:c.6del")).unwrap_or_else(|e| {
            panic!("{mode:?}: a bare genomic parent is outside the ruling's scope, got: {e:?}")
        });
        assert_eq!(
            out,
            format!("{GENOMIC_PARENT}:c.6del"),
            "{mode:?}: the pre-existing pass-through on a genomic parent must be unchanged"
        );
    }
}

/// The prefix rule is untouched: a `c.` description on a non-coding accession is
/// still refused at PARSE, in every mode, with its own message.
///
/// Recorded here because this ruling is easy to misread as replacing it. It does
/// not: an `NR_`/`XR_` carries the fact in its own accession, so no reference is
/// needed and parse is the right stage. What the ruling adds is the case where
/// the fact is the *reference's* and the prefix says nothing.
#[test]
fn the_prefix_rule_is_untouched_and_still_refuses_at_parse() {
    let err = parse_hgvs("NR_000005.1:c.100del")
        .expect_err("`c.` on a non-coding accession is refused at parse");
    let msg = err.to_string();
    assert!(
        msg.contains("has no CDS"),
        "the parse-stage prefix refusal must still name the missing CDS, got: {msg}"
    );
}
