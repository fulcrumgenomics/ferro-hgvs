//! Issue #1789 — an insertion sized by a number states no inserted sequence.
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/1789>.
//!
//! # The clause
//!
//! `docs/recommendations/checklist.md` item 3 ("Insertions") has two bullets.
//! The second, `:32`–`:33`, quoted from spec checkout `6f85311`:
//!
//! ```text
//! 32: - **do you provide the inserted sequence?**<br>
//! 33:   Describing a variant as <code class="invalid">c.5439_5430ins6</code>
//!       is not allowed, the inserted sequence (for `ins6`, e.g., `TGCCAT`)
//!       should be specified.
//! ```
//!
//! The force of `:33` is not a keyword question. `recommendations/style.md:9`
//! binds the RFC 2119 reading to the ten UPPERCASE key words, and `:33` uses
//! none of them — nor does almost any other clause in this corpus, so keyword
//! strength cannot rank them. What grounds `:33` is the checklist's own frame
//! and its worked example: `checklist.md:5` introduces the list as descriptions
//! that "**do not correctly follow HGVS nomenclature**", and `:33` marks
//! `c.5439_5430ins6` `class="invalid"` — the spec's verdict on a concrete
//! description rather than prose about one. What ferro does about it is then
//! fixed by the decided `rulings[absolute-prohibition-enforcement-stage]`,
//! which names `:33` in its own clause list (see "The stage" below).
//!
//! Its sibling `:49` ("Descriptions like `g.123del3` are not allowed, correct
//! is `g.123_125del`") is the same construction — a numbered checklist bullet
//! whose example carries `class="invalid"` — and is the deletion case ferro has
//! enforced since #1079. So before this change one such bullet was enforced for
//! `del` and not for `ins`.
//!
//! **The axis clauses say the same thing constructively**, which matters because
//! `checklist.md` is not molecule-scoped and a `DNA/` clause cannot scope `r.`.
//! Both enumerate what an inserted sequence may be, and neither list contains a
//! bare count:
//!
//! - `DNA/insertion.md:22` — "the **\"inserted_sequence\"** can be given as the
//!   nucleotides inserted (e.g., `insAGC`) or, for larger insert sequences, by
//!   referring to the sequence in the reference sequence (e.g.,
//!   `c.849_850ins858_895`) or another reference".
//! - `RNA/insertion.md:20` — the same sentence for `r.`, with `insagc` and
//!   `r.849_850ins858_895`.
//!
//! And both name the spelling that *is* conformant when the number of inserted
//! nucleotides is known but their identity is not: `DNA/insertion.md:77`
//! (`g.32717298_32717299insN[100]`, "the insertion of 100 nucleotides (not
//! specified)") and `RNA/insertion.md:41` (`r.1149_1150insn[100]`). So the rule
//! costs no expressive power — `ins6` has an exact conformant translation,
//! `insN[6]`, that ferro already parses and re-emits.
//!
//! `DNA/insertion.md:119` states the same verdict about a real reported
//! description: "`c.23ins24` is not correct since the position of the insertion
//! is not described properly and because \"ins24\" does not define the sequence
//! inserted."
//!
//! # Why it stayed invisible
//!
//! The spec's own example at `:33` carries an **inverted range** —
//! `c.5439_5430ins6`, 5439 > 5430 — so ferro refuses it for the *range*
//! (`Interval positions are swapped`) and the conformance corpus scores the row
//! `correctly-rejected`. The `ins6` acceptance was masked by an unrelated
//! rejection. Upstream corrected the range in `14fcbf9`, which is in the
//! unmerged `6f85311..565b973` submodule range.
//!
//! **Every input below therefore uses a VALID range of this file's own**, so the
//! guard holds at both spec pins and cannot be masked again by the example's
//! defect. [`the_guard_does_not_depend_on_the_spec_example`] asserts that
//! property directly.
//!
//! # The stage, which is not this file's to choose
//!
//! The decided `rulings[absolute-prohibition-enforcement-stage]` sets it, and
//! names `checklist.md:33` in its own clause list:
//!
//! | mode | parse | normalize |
//! |---|---|---|
//! | strict | **fails**, `W3029` | (unreached) |
//! | lenient | accepts, warns `W3029` | **fails** |
//! | silent | accepts, quiet | **fails** |
//!
//! The normalize rung is **not** mode-gated, on that record's own reasoning:
//! the mode gate governs whether the INPUT is judged and "does not, and cannot,
//! govern whether the output conforms". Lenient fails there on the ground the
//! record gives it — it cannot normalize — because `ins6` denotes no sequence
//! (`corpus_prohibited_inputs::an_insertion_stating_a_length_denotes_no_sequence`
//! measures exactly that). That record also predicted this outcome in as many
//! words: the `ins6` rows "genuinely denote nothing … But lenient does NOT fail
//! on them. The normalizer passes the offending member through untouched and
//! emits it back — in all three modes, byte-identically, with an EMPTY warning
//! vector. Normalization is not impossible here, it is VACUOUS."
//!
//! # Scope, and the two neighbours deliberately left alone
//!
//! **`InsertedSequence::Range` (`ins(10_20)`) is NOT matched.** In ferro's AST
//! that variant is ambiguous between a count range and an uncertain *position*
//! range: `LRG_308:g.?_?ins(23632682_23625413)_(23625324_23619334)` is a
//! position range the spec sanctions (`DNA/complex.md`), and its single-bracket
//! sibling parses to the same `Range`. Refusing it would risk refusing a
//! conformant description, so the widening needs its own measurement rather than
//! being folded in here. [`a_count_range_insert_is_deliberately_untouched`] pins
//! that boundary.
//!
//! **`delins<number>` is NOT matched.** `checklist.md:33` sits under item 3,
//! "Insertions", and `DNA/delins.md` states no equivalent prohibition — its
//! `class="invalid"` notes are about restating the deleted sequence (`:29`,
//! `:34`), not about the insert's form. Extending there would be an
//! adjudication no clause supports.
//!
//! **`inv<number>` is NOT matched, and is not a defect.** See
//! [`an_inversion_sized_by_a_number_is_a_different_shape`] for the measurement.

use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::edit::{InsertedSequence, NaEdit};
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::hgvs::variant::HgvsVariant;
use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

// ---------------------------------------------------------------------------
// Harness
// ---------------------------------------------------------------------------

/// A single-exon transcript with `cds_start = 1`, so `c.N == n.N == r.N` and the
/// axis is the only thing that differs between the probes below.
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_U.1".to_string(),
        Some("UG".to_string()),
        Strand::Plus,
        "AATTTGCCAATTTGCCAATTTGCC".to_string(),
        Some(1u64),
        Some(24u64),
        vec![Exon::new(1, 1, 24)],
        None,
        None,
        None,
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

/// Normalize `input` in `config`'s mode, from parse through normalize, as a
/// caller would. `Err` carries whichever stage refused.
fn round_trip(input: &str, config: ErrorConfig) -> Result<String, String> {
    let parsed = parse_hgvs_with_config(input, config).map_err(|e| e.to_string())?;
    Normalizer::new(provider())
        .normalize(&parsed.result)
        .map(|v| v.to_string())
        .map_err(|e| e.to_string())
}

/// Insertions whose payload is a bare count, one per shape that reaches
/// [`InsertedSequence::Count`] and one per nucleic-acid axis that can carry it.
///
/// Both surface spellings are here because they are different parser arms —
/// `ins6` comes from `parse_simple_count` and `ins(6)` from the parenthesized
/// arm — and a check keyed on the rendered text would have caught one and missed
/// the other. Both land on the same AST node, which is the point of keying the
/// rule on the AST.
const SIZE_COUNT_INSERTS: &[&str] = &[
    "NM_U.1:c.10_11ins6",
    "NM_U.1:c.10_11ins(6)",
    "NM_U.1:n.10_11ins6",
    "NM_U.1:r.10_11ins6",
    "NC_TEST.1:g.10_11ins6",
    "NC_TEST.1:m.10_11ins6",
];

/// The conformant translations of the same descriptions, which must stay
/// accepted — otherwise the rule above is refusing "an insertion of six
/// unspecified nucleotides" rather than the spelling that fails to state one.
const CONFORMANT_CONTROLS: &[&str] = &[
    "NM_U.1:c.10_11insTGCCAT",
    "NM_U.1:c.10_11insN[6]",
    "NC_TEST.1:g.10_11insN[6]",
    "NM_U.1:c.10_11ins1_6",
    // The `r.` row is not decoration: it is the exact spelling the RNA
    // diagnostic *offers*, so without it the suggestion could name a
    // description ferro refuses and nothing here would notice.
    "NM_U.1:r.10_11insn[6]",
];

// ---------------------------------------------------------------------------
// The refusal
// ---------------------------------------------------------------------------

/// Strict refuses at parse; lenient and silent accept at parse and fail at
/// normalize. The schedule is `rulings[absolute-prohibition-enforcement-stage]`'s
/// and is asserted in the same shape
/// `issue_1715_rna_alignment_symbol_reach::`
/// `a_lowercase_masked_nucleotide_on_the_rna_axis_is_refused_at_the_ruled_stage`
/// uses for `W3028`, its nearest sibling.
#[test]
fn a_size_count_insertion_is_refused_at_the_ruled_stage() {
    for input in SIZE_COUNT_INSERTS {
        // STRICT — refused at PARSE. Strict validates input conformance, not
        // merely parseability.
        let err = parse_hgvs_with_config(input, ErrorConfig::strict())
            .map(|p| panic!("{input}: strict must refuse at parse, got {:?}", p.result))
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("W3029"),
            "{input}: strict must refuse with W3029; got: {err}"
        );
        assert!(
            err.contains("checklist.md:33"),
            "{input}: the refusal must cite the clause; got: {err}"
        );

        // LENIENT and SILENT — accepted at parse, because neither validates
        // input conformance.
        let lenient = parse_hgvs_with_config(input, ErrorConfig::lenient())
            .expect("lenient accepts at parse");
        assert!(
            lenient
                .warnings
                .iter()
                .any(|w| w.error_type.code() == "W3029"),
            "{input}: lenient must say why it will fail later"
        );
        let silent =
            parse_hgvs_with_config(input, ErrorConfig::silent()).expect("silent accepts at parse");
        assert!(
            !silent
                .warnings
                .iter()
                .any(|w| w.error_type.code() == "W3029"),
            "{input}: silent is lenient without messages"
        );

        // All three modes refuse to NORMALIZE. Output conformance is rule 1 of
        // the README ruleset and has no mode escape, so this rung is not
        // mode-gated: lenient fails on the ruling's own ground, because a count
        // names no bases to normalize.
        for (label, config) in [
            ("strict", ErrorConfig::strict()),
            ("lenient", ErrorConfig::lenient()),
            ("silent", ErrorConfig::silent()),
        ] {
            let err =
                round_trip(input, config).expect_err("every mode must refuse a size-count insert");
            assert!(
                err.contains("W3029"),
                "{input}: the {label} refusal must name W3029; got: {err}"
            );
        }
    }
}

/// The rule refuses the **spelling**, not the variant class. Every conformant
/// way of saying the same thing must still round-trip, or the change has
/// deleted expressive power rather than a prohibited form.
#[test]
fn the_conformant_spellings_are_untouched() {
    for input in CONFORMANT_CONTROLS {
        assert!(
            parse_hgvs_with_config(input, ErrorConfig::strict()).is_ok(),
            "{input}: strict must still accept the conformant spelling"
        );
        assert!(
            round_trip(input, ErrorConfig::lenient()).is_ok(),
            "{input}: the conformant spelling must still normalize"
        );
    }
}

/// A member of a composite is reached on its own, so the spelling cannot hide
/// behind a bracketed description. All three shapes, because "the cis case
/// works" does not entail the other two: a trans pair nests one `Allele` inside
/// another, and an uncertain group sets `uncertain` on it.
#[test]
fn a_member_of_a_composite_is_refused_too() {
    for input in [
        "NM_U.1:c.[10_11ins6;20del]",
        "NM_U.1:c.[10_11ins6];[20del]",
        "NM_U.1:c.(10_11ins6)",
    ] {
        let err = parse_hgvs_with_config(input, ErrorConfig::strict())
            .map(|_| panic!("{input}: strict must refuse"))
            .unwrap_err()
            .to_string();
        assert!(err.contains("W3029"), "{input}: got {err}");
        let err = round_trip(input, ErrorConfig::lenient())
            .expect_err("lenient cannot normalize a size-count insert");
        assert!(err.contains("W3029"), "{input}: got {err}");
    }
}

/// **The diagnostic is axis-aware in both halves, and both are asserted here.**
///
/// `general.md:50` puts an `r.` description's nucleotides in lower case, so an
/// RNA author must be offered `insn[6]` and not `insN[6]`. And a `DNA/` clause
/// cannot scope `r.` — the reason this rule cites two axis clauses rather than
/// `checklist.md:33` alone — so the citations swap with the spelling:
/// `RNA/insertion.md:20`/`:41` rather than `DNA/insertion.md:22`/`:77`, with the
/// worked range example spelled on the author's own axis.
///
/// Without this test the axis branch is a no-op guard: every other assertion in
/// this file matches on `W3029` alone, so swapping
/// `conformant_spelling_for_rna` and `conformant_spelling` at both call sites
/// leaves the suite green. Written to fail on that swap specifically — each
/// side asserts the *absence* of the other axis's spelling as well as the
/// presence of its own.
#[test]
fn the_diagnostic_cites_the_axis_the_author_is_writing() {
    // (input, expected spelling, expected clauses, the other axis's spelling)
    let cases: [(&str, &str, [&str; 3], &str); 2] = [
        (
            "NM_U.1:r.10_11ins6",
            "insn[6]",
            ["RNA/insertion.md:20", "RNA/insertion.md:41", "r.849_850ins"],
            "insN[6]",
        ),
        (
            "NC_TEST.1:g.10_11ins6",
            "insN[6]",
            ["DNA/insertion.md:22", "DNA/insertion.md:77", "c.849_850ins"],
            "insn[6]",
        ),
    ];

    for (input, spelling, clauses, other_axis_spelling) in cases {
        // STRICT, at parse — the message carrying the citations.
        let err = parse_hgvs_with_config(input, ErrorConfig::strict())
            .map(|_| panic!("{input}: strict must refuse"))
            .unwrap_err()
            .to_string();
        assert!(
            err.contains(spelling),
            "{input}: the strict refusal must offer `{spelling}`; got: {err}"
        );
        assert!(
            !err.contains(other_axis_spelling),
            "{input}: the strict refusal must not offer the other axis's \
             `{other_axis_spelling}`; got: {err}"
        );
        for clause in clauses {
            assert!(
                err.contains(clause),
                "{input}: the strict refusal must cite `{clause}`; got: {err}"
            );
        }

        // LENIENT, at parse — the W3029 warning says what to write instead.
        let warning = parse_hgvs_with_config(input, ErrorConfig::lenient())
            .expect("lenient accepts at parse")
            .warnings
            .into_iter()
            .find(|w| w.error_type.code() == "W3029")
            .expect("lenient warns")
            .message;
        assert!(
            warning.contains(spelling) && !warning.contains(other_axis_spelling),
            "{input}: the lenient warning must offer `{spelling}`; got: {warning}"
        );

        // NORMALIZE — the rung that is not mode-gated carries it too.
        let err = round_trip(input, ErrorConfig::lenient())
            .expect_err("normalize refuses a size-count insert");
        assert!(
            err.contains(spelling) && !err.contains(other_axis_spelling),
            "{input}: the normalize refusal must offer `{spelling}`; got: {err}"
        );
    }
}

// ---------------------------------------------------------------------------
// The scope lines
// ---------------------------------------------------------------------------

/// **The guard is independent of the spec's own example**, which is the whole
/// reason the defect survived: `c.5439_5430ins6` at `:33` has 5439 > 5430, so
/// ferro refuses it for the *range* and the corpus scores the row
/// `correctly-rejected` while the `ins6` acceptance goes unmeasured.
///
/// Asserted rather than asserted-about: the example is shown to be refused for a
/// reason that is **not** this rule, and the file's own valid-range inputs are
/// shown to be refused for a reason that **is**. So a submodule bump to
/// `565b973` (which corrects the range upstream, `14fcbf9`) cannot silently
/// change what this file proves in either direction.
#[test]
fn the_guard_does_not_depend_on_the_spec_example() {
    let example = parse_hgvs_with_config("NM_004006.2:c.5439_5430ins6", ErrorConfig::strict())
        .map(|_| panic!("the spec's example must be refused"))
        .unwrap_err()
        .to_string();
    assert!(
        !example.contains("W3029"),
        "the spec example is refused for its inverted range, not for `ins6` — if that \
         changes, this file's inputs are the ones carrying the guard; got: {example}"
    );

    for input in SIZE_COUNT_INSERTS {
        assert!(
            !input.contains("5439"),
            "{input}: no input here may be the spec's example"
        );
    }
}

/// **`ins(10_20)` is deliberately untouched.** `InsertedSequence::Range` is
/// ferro's node for both "a count range" and "an uncertain position range", and
/// the second is spec-sanctioned. Refusing the node would refuse conformant
/// descriptions, so the widening is a separate question with its own
/// measurement.
///
/// Pinned so the exclusion is a decision rather than an oversight: if a later
/// change splits the node, this test fails and the boundary gets re-argued.
#[test]
fn a_count_range_insert_is_deliberately_untouched() {
    let input = "NM_U.1:c.10_11ins(6_9)";
    let parsed = parse_hgvs(input).expect("the grammar admits a parenthesized range");
    let HgvsVariant::Cds(cds) = &parsed else {
        panic!("{input}: expected a c. variant");
    };
    assert!(
        matches!(
            cds.loc_edit.edit.inner(),
            Some(NaEdit::Insertion {
                sequence: InsertedSequence::Range(_, _)
            })
        ),
        "{input}: expected the ambiguous Range node, got {:?}",
        cds.loc_edit.edit.inner()
    );
    assert!(
        parse_hgvs_with_config(input, ErrorConfig::strict()).is_ok(),
        "{input}: Range is out of scope, so strict must still accept it"
    );
}

/// **`inv<number>` is not the same defect, measured three ways.**
///
/// The issue filed alongside this one observed that `g.123_125inv3` is also
/// accepted, and asked whether it is a parser laxity. It is neither that nor a
/// live conformance defect:
///
/// 1. **It is not a trailing token being absorbed.** `NaEdit::Inversion` has a
///    first-class `length: Option<u64>` field, documented "Explicit length of
///    inverted region (e.g., `inv3`)", and `parse_inversion` fills it
///    deliberately. So there is nothing accidental to fix.
/// 2. **No clause reaches it.** `checklist.md:49` is about naming "the first and
///    last residue involved in a deletion", and `g.123_125inv3` names both
///    endpoints already — the number is redundant, not a substitute for the
///    range. That is the class ferro already treats as soft
///    (`variant.rs::validate_no_point_size_suffix` says so of `c.100_102del3`),
///    and `DNA/inversion.md` states no prohibition on a size suffix at all.
/// 3. **The output already conforms.** Normalization drops the count — the rung
///    `checklist.md:33` needed and did not have. That is asserted below rather
///    than argued, because it is the whole difference between the two shapes.
///
/// The single-position spelling `g.123inv3`, which *would* be `:49`'s shape, is
/// already refused in every mode by the "an inversion covers more than one
/// nucleotide" check.
#[test]
fn an_inversion_sized_by_a_number_is_a_different_shape() {
    // The count is dropped on output, in every mode — so no prohibited-looking
    // token is re-emitted and rule 1 is not engaged. Asserted on the *token*
    // rather than on a fixed string, because what an inversion normalizes to
    // depends on the bases under it (over this transcript's sequence
    // `c.10_12inv` reduces further, to a substitution) and that reduction is
    // not what is being measured here.
    for config in [
        ErrorConfig::strict(),
        ErrorConfig::lenient(),
        ErrorConfig::silent(),
    ] {
        let out =
            round_trip("NM_U.1:c.10_12inv3", config.clone()).expect("an inversion normalizes");
        assert!(
            !out.contains("inv3"),
            "normalization must drop the redundant size from an inversion; got {out}"
        );
        // The control: the same description without the size normalizes to the
        // same thing, so the size is genuinely inert rather than the input being
        // refused on some other ground.
        assert_eq!(
            round_trip("NM_U.1:c.10_12inv", config).as_deref(),
            Ok(out.as_str()),
            "`inv3` and `inv` over the same range must agree"
        );
    }

    // The single-position spelling — the one `:49`'s reasoning would reach — is
    // already refused, in every mode, by the grammar rather than by a rule.
    for config in [
        ErrorConfig::strict(),
        ErrorConfig::lenient(),
        ErrorConfig::silent(),
    ] {
        assert!(
            parse_hgvs_with_config("NM_U.1:c.10inv3", config).is_err(),
            "a single-position inversion is refused in every mode"
        );
    }
}
