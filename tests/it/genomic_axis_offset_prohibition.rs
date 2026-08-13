//! `background/numbering.md:6`/`:8`/`:11` — no `+`/`-` offset on a `g.`, `o.`
//! or `m.` position, enforced at the stage the ruling names.
//!
//! # The question, and why it was not already answered
//!
//! `docs/recommendations/checklist.md:16` names three symbols in one sentence —
//! a genomic reference "can not have nucleotides with additions like a `+`, `-`,
//! or `*`" — and ferro enforced exactly one of them. Measured on `9fb126ba`,
//! through every entry point:
//!
//! ```text
//!                        parse_hgvs   parse strict/lenient/silent   normalize strict/lenient/silent
//! g.100+2A>G             accepted     accepted / accepted / accepted   re-emitted verbatim, all three
//! g.100-2A>G             accepted     accepted / accepted / accepted   re-emitted verbatim, all three
//! g.*100A>G              refused      refused  / refused  / refused    (does not parse)
//! g.-100A>G              refused      refused  / refused  / refused    (does not parse)
//! m.100+2A>G             accepted     accepted / accepted / accepted   re-emitted verbatim, all three
//! o.100+2del             accepted     accepted / accepted / accepted   re-emitted verbatim, all three
//! ```
//!
//! **The `*` column is a grammar accident, not enforcement.** `*` and a leading
//! `-` match no production of the genomic position parser, so they are refused
//! for the same reason `g.@100` would be — which is why the clause read as
//! "one-third enforced" rather than "not enforced". [`the_star_and_leading_hyphen_forms_are_grammar_refusals`]
//! pins that half so the distinction cannot be lost again.
//!
//! # Authority
//!
//! The clause is stated **per axis, three times**, in three consecutive bullets
//! of `background/numbering.md` — `:6` for `g.`, `:8` for `o.`, `:11` for `m.`
//! — each ending "Nucleotide numbers based on a … reference sequence **do not
//! include** `+`, `-`, `*`, or other prefixes." So covering `m.` and `o.` is not
//! a widening of a `g.`-only clause; refusing on `g.` alone would be the
//! under-reading. `checklist.md:16` says the same of `g.` in the checklist's
//! register, and `checklist.md:45` supplies the shape that actually occurs — a
//! range written with `-` where the separator is `_`, which it marks
//! `Not correct`.
//!
//! # The stage
//!
//! The decided `rulings[absolute-prohibition-enforcement-stage]` (2026-08-10)
//! governs: **strict** validates input conformance and so fails at PARSE;
//! **lenient** does not validate input conformance and fails only when it cannot
//! NORMALIZE; **silent** is lenient without messages. The output half is not a
//! mode question — rule 1 of the README ruleset ("Output follows the HGVS
//! recommendations. Absolute — never traded.") is about OUTPUT and has no mode
//! escape — so normalize refuses in every mode, and that is what
//! [`no_mode_may_emit_the_prohibited_spelling`] asserts.
//!
//! `tests/it/corpus_prohibited_inputs.rs::the_decided_target_is_a_mode_gated_refusal`
//! is the shared acceptance criterion for all four clauses in that record's
//! census; this file closes its two `g.`-offset rows. It stays `#[ignore]`d
//! because its `checklist.md:33` `ins6` row is still open (#1627).
//!
//! # What is NOT claimed here
//!
//! Ferro never **manufactures** a genomic offset. `GenomePos::with_offset` has
//! exactly one caller in the crate — `parse_genome_pos` — so no projection,
//! conversion or normalization path can mint one; every prohibited output was a
//! propagated input. That is a lower severity than a self-emitted violation, and
//! it is stated rather than glossed. What made it a rule-1 defect anyway is that
//! `normalize` **re-emitted** the propagated spelling as a fixed point, in all
//! three modes, with an empty warning vector.

use ferro_hgvs::error_handling::{ErrorConfig, ErrorMode};
use ferro_hgvs::hgvs::parser::{parse_hgvs, parse_hgvs_with_config};
use ferro_hgvs::normalize::NormalizeConfig;
use ferro_hgvs::{MockProvider, Normalizer};

/// Every genomic-family accession this file names, with enough sequence that a
/// normalization attempt reaches the conformance refusal rather than dying on a
/// short contig.
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    for accession in [
        "NC_000023.10",
        "NC_012920.1",
        "NC_001416.1",
        "NC_000001.11",
        "NC_TEST.1",
    ] {
        provider.add_genomic_sequence(accession, "ACGT".repeat(100_000));
    }
    provider
}

fn normalize_in(mode: ErrorMode, input: &str) -> Result<String, String> {
    let variant = parse_hgvs(input).map_err(|e| e.to_string())?;
    let config = match mode {
        ErrorMode::Strict => NormalizeConfig::strict(),
        ErrorMode::Lenient => NormalizeConfig::lenient(),
        ErrorMode::Silent => NormalizeConfig::silent(),
    };
    Normalizer::with_config(provider(), config)
        .normalize(&variant)
        .map(|v| v.to_string())
        .map_err(|e| e.to_string())
}

/// Every shape that states an offset on a genomic-family position, with the
/// clause governing its axis. One list, so no test can cover a subset by
/// accident.
const OFFENDING: &[(&str, &str)] = &[
    // `checklist.md:16`'s `+`, on each of the three axes.
    ("NC_000023.10:g.100+2A>G", "background/numbering.md:6"),
    ("NC_012920.1:m.100+2A>G", "background/numbering.md:11"),
    ("NC_001416.1:o.100+2del", "background/numbering.md:8"),
    // `checklist.md:45`'s hyphen range, which the grammar reads as an offset —
    // on each of the three axes, for the same reason the `+` rows above are:
    // `:6`, `:8` and `:11` state the prohibition separately, so a matrix that
    // carried the hyphen on `g.` alone would be covering a subset by accident,
    // which the comment above says this list exists to prevent.
    ("NC_000023.10:g.100-2A>G", "background/numbering.md:6"),
    ("NC_012920.1:m.100-2A>G", "background/numbering.md:11"),
    ("NC_001416.1:o.100-2del", "background/numbering.md:8"),
    // The real-world instance of that shape: two rows of `clinvar_hgvs_unique`
    // spell a repeat span with a hyphen. They are the only genuine genomic
    // offsets in 9,949,738 committed corpus rows.
    (
        "NC_000023.10:g.146993569-146993571[6]",
        "background/numbering.md:6",
    ),
    // The unknown-offset sentinels render `+`/`-` exactly as a numeric offset
    // does; `:6` prohibits the symbol, not a magnitude.
    ("NC_000023.10:g.100+?del", "background/numbering.md:6"),
    ("NC_000023.10:g.100-?del", "background/numbering.md:6"),
    // Inside a cis allele — a per-description check that missed this would
    // report the same verdict for the standalone member and not for the allele.
    ("NC_TEST.1:g.[265del;266+2del]", "background/numbering.md:6"),
    // A complex `(a_b)` boundary. `UncertainBoundary::inner()` returns `None`
    // for a `Range`, so a walker built on it alone skips this shape — which is
    // why the detector enumerates both of a range's endpoints.
    (
        "NC_000001.11:g.(100+1_101-1)_(200_300)del",
        "background/numbering.md:6",
    ),
];

/// **Question.** Where does strict refuse a genomic offset?
///
/// **At parse**, per `rulings[absolute-prohibition-enforcement-stage]`: strict
/// validates whether the INPUT conforms, which is answered before the input is
/// accepted rather than part-way through normalizing it.
///
/// Keyed on the message, never on a bare `is_err()` — a bare check would count a
/// short contig or a reference mismatch as a successful refusal.
#[test]
fn strict_refuses_a_genomic_offset_at_parse() {
    for (input, clause) in OFFENDING {
        let err = parse_hgvs_with_config(input, ErrorConfig::strict())
            .expect_err(&format!("STRICT must refuse `{input}` at parse"))
            .to_string();
        assert!(
            err.contains(clause),
            "`{input}` must be refused citing its own axis's clause `{clause}`; got: {err}"
        );
        assert!(
            err.contains("W4009"),
            "`{input}` must carry the W4009 tag; got: {err}"
        );
    }
}

/// **Question.** Do lenient and silent refuse at parse too?
///
/// **No — and that is the ruling, not an oversight.** Those modes do not
/// validate input conformance, so the description parses; the refusal falls to
/// normalize. The mode gate exists so a caller round-tripping a real-world
/// corpus is not blocked at the door.
#[test]
fn lenient_and_silent_accept_a_genomic_offset_at_parse() {
    for (input, _) in OFFENDING {
        for (mode, config) in [
            ("lenient", ErrorConfig::lenient()),
            ("silent", ErrorConfig::silent()),
        ] {
            let parsed = parse_hgvs_with_config(input, config)
                .unwrap_or_else(|e| panic!("{mode} must accept `{input}` at parse; got {e}"));
            assert_eq!(
                parsed.result.to_string(),
                *input,
                "{mode} must accept `{input}` unchanged — there is no repair arm"
            );
        }
    }
}

/// **Question.** May a permissive mode hand the description back?
///
/// **No.** Rule 1 of the README ruleset is about OUTPUT and carries no mode
/// escape, so `normalize` refuses in every mode. There is nothing to be lenient
/// *toward*: a genomic accession has no exon table, so the offset is measured
/// from nothing and the position names no nucleotide.
///
/// This is the assertion that would have failed before the fix. Measured on
/// `9fb126ba`, all three modes returned the input byte-identically with an empty
/// warning vector — normalization was not impossible, it was **vacuous**, which
/// is worse, because the output looks normalized.
#[test]
fn no_mode_may_emit_the_prohibited_spelling() {
    for (input, clause) in OFFENDING {
        for mode in [ErrorMode::Strict, ErrorMode::Lenient, ErrorMode::Silent] {
            let err = normalize_in(mode, input).expect_err(&format!(
                "{mode:?} must refuse to normalize `{input}` rather than re-emit it"
            ));
            assert!(
                err.contains(clause) && err.contains("W4009"),
                "{mode:?}'s refusal of `{input}` must cite `{clause}` and W4009; got: {err}"
            );
        }
    }
}

/// **Question.** Is `g.*100A>G` / `g.-100A>G` evidence that the clause was
/// already partly enforced?
///
/// **No, and reading it that way is what hid the defect.** Neither matches any
/// production of the genomic position grammar, so both are refused the way any
/// unparseable text is — in every mode, with a nom error rather than a clause
/// citation. Pinned here so the difference between a grammar accident and a
/// conformance rule stays visible.
#[test]
fn the_star_and_leading_hyphen_forms_are_grammar_refusals() {
    for input in ["NC_000023.10:g.*100A>G", "NC_000023.10:g.-100A>G"] {
        assert!(
            parse_hgvs(input).is_err(),
            "`{input}` has no production and must stay refused at the bare entry"
        );
        for (mode, config) in [
            ("strict", ErrorConfig::strict()),
            ("lenient", ErrorConfig::lenient()),
            ("silent", ErrorConfig::silent()),
        ] {
            let err = parse_hgvs_with_config(input, config)
                .expect_err(&format!("`{input}` must be refused in {mode}"))
                .to_string();
            assert!(
                !err.contains("W4009"),
                "`{input}` is a grammar refusal, not this rule's: {err}"
            );
        }
    }
}

/// **Question.** Does the bare `parse_hgvs` entry refuse?
///
/// **No.** It applies no `ErrorConfig` at all (#1632), so it runs neither the
/// repairs lenient performs nor the refusals that make strict strict. It is the
/// same answer W3026 and W3028 give at that entry, and it is why the normalize
/// half has to be unconditional: a caller who never touches a config still must
/// not be handed a prohibited spelling back.
#[test]
fn the_bare_parse_entry_is_neither_mode_and_still_accepts() {
    for (input, _) in OFFENDING {
        let parsed = parse_hgvs(input)
            .unwrap_or_else(|e| panic!("the bare entry must still accept `{input}`; got {e}"));
        assert_eq!(parsed.to_string(), *input);
    }
}

/// **Question.** Does the check reach anything it should not?
///
/// Three controls, each for a way a text-keyed gate would misfire: a gene
/// symbol carries `-`, an accession carries digits and dots, and the transcript
/// axes are where an offset is *legitimate* (`numbering.md:53` grants intronic
/// offsets by name). The check is AST-keyed on `GenomePos::offset`, so none of
/// these is touched.
///
/// The obvious fourth control — a repeat count spelled as a range,
/// `NG_008047.1:g.17267CAG[(54-68)]` — is **not** here because ferro's grammar
/// does not accept that form at all ("Unexpected trailing characters:
/// `[(54-68)]`"), so it would pin the grammar rather than this rule. Six such
/// rows exist in `clinvar_hgvs_unique`; they are already refused, for an
/// unrelated reason, and this change does not move them.
#[test]
fn legitimate_neighbours_are_untouched() {
    for input in [
        // Plain genomic positions.
        "NC_000023.10:g.100A>G",
        "NC_000023.10:g.100_200del",
        "NC_012920.1:m.3243A>G",
        // A gene symbol containing a hyphen, beside a genomic position — the
        // shape a text-keyed gate would refuse.
        "NC_012920.1(MT-ND1):m.3460G>A",
        // The transcript axes. `numbering.md:53` grants these explicitly.
        "NM_000088.3:c.100+5A>G",
        "NM_000088.3:c.101-3del",
        "NR_003051.3:n.100+5A>G",
    ] {
        for (mode, config) in [
            ("strict", ErrorConfig::strict()),
            ("lenient", ErrorConfig::lenient()),
            ("silent", ErrorConfig::silent()),
        ] {
            let parsed = parse_hgvs_with_config(input, config)
                .unwrap_or_else(|e| panic!("`{input}` must still parse in {mode}; got {e}"));
            assert_eq!(parsed.result.to_string(), *input, "in {mode}");
        }
    }
}
