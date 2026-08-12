//! Issue #1748 — `n.*N` and `n.-N` name zones the non-coding axis does not have.
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/1748>.
//!
//! # The clause, and why `:53` is what makes it bite
//!
//! `background/numbering.md:50`–`:54` enumerates the whole `n.` axis:
//!
//! ```text
//! 52: - nucleotide numbering is `n.1`, `n.2`, `n.3`, ..., etc., from the first
//!      to the last nucleotide of the reference sequence.
//! 53: - nucleotides in introns are numbered as for coding DNA reference
//!      sequences (see above), although proceeded by `n.` (not `c.`).
//! 54: - it is **not** allowed to describe variants in nucleotides beyond the
//!      boundaries of a transcript reference sequence, using that transcript
//!      reference sequence.
//! ```
//!
//! `:52` is exhaustive — first nucleotide to last — and names no `*` zone and no
//! negative one. On its own that could be read as terseness. `:53` is what rules
//! that reading out: the spec knew how to add a zone to this axis, and added
//! exactly one (introns). `:54` then forbids naming anything past the boundary
//! outright, and `numbering.md:45` records that the sibling proposal to mark
//! non-transcribed nucleotides "[has] been made but [was] rejected".
//!
//! # The one reading that would save these spellings, and why it does not
//!
//! `:53` incorporates the coding numbering by reference — "numbered as for
//! coding DNA reference sequences (see above)" — and the coding numbering is
//! where `-N` and `*N` are defined. So the obvious counter-argument is that
//! `:53` pulls those two zones onto `n.` along with the intronic ones, which is
//! very likely the premise `issue_255_noncoding_markers.rs` was written on.
//!
//! **It does not, and the structure of the section it points at is what settles
//! it.** The coding numbering (`:18`–`:44`) is three named bullets:
//!
//! | bullet | line | defines |
//! |---|---|---|
//! | **protein coding region** | `:20` | `c.1` … the stop codon |
//! | **untranslated region (UTR)** | `:28` | `c.-1`, `c.-2`, … and `c.*1`, `c.*2`, … |
//! | **introns** | `:33` | `c.87+1`, `c.88-1`, … |
//!
//! `:53`'s subject is *"nucleotides **in introns**"*, so its "see above" is
//! scoped by its own sentence to the **introns** bullet. The `-`/`*` zones are
//! in the **UTR** bullet, which `:53` never reaches. There is no clause tension
//! here to adjudicate — one clause is narrower than it first reads.
//!
//! The coding section's own closing bullet is the clincher, because it is
//! `:54`'s sentence verbatim: `:43`–`:44` **transcript flanking**. So `c.*1` and
//! "nothing beyond the boundary" coexist on the coding axis *because the UTR is
//! inside the transcript*. `n.` has no UTR bullet, so nothing puts `-`/`*`
//! inside its boundary and `:54` stands unopposed there.
//!
//! # ONE READING, TWO STAGES — and the split is measured, not argued
//!
//! This is the part to read first, because it is the only thing here that is not
//! derived from the clause:
//!
//! | marker | stage | code | corpus rows |
//! |---|---|---|---|
//! | `n.*N` | refused at parse in **every mode**, bare entry included | `E1003` | **0** of 103,762 |
//! | `n.-N` | refused at **strict** parse only | `W4008` | **5** of 103,762 |
//!
//! Counted over ferro's four committed corpora — `clinvar_hgvs_500k`,
//! `clinvar_hgvs_unique`, `cmrg_genes_exhaustive` and
//! `paraphase_genes_exhaustive`. The `n.*N` zero is a **real** zero rather than a
//! structural one, and the `n.-N` five is the proof: the same scan over the same
//! rows finds the sibling shape five times, so the counter can see this class
//! when it is there.
//!
//! Those five are `NR_003051.3:n.-57T>C`, `NR_003051.3:n.-30_-7dup` and
//! `LRG_163t1:n.-5delins17` — all RMRP, whose upstream promoter variants are the
//! clinically conventional case for this spelling — plus
//! `NR_029595.1:n.-4771G>T` (MIR208A) and `NR_033294.1:n.-6G>A` (SNORD118).
//! They are descriptions NCBI publishes today.
//!
//! So: `n.*N` is refused because the axis has no such zone **and nothing in the
//! wild uses it**; `n.-N` is deviation-tolerated in the permissive modes
//! precisely **because real clinical data uses it**, on the identical clause
//! reading. Whoever revisits this needs to see that the split was measured.
//!
//! # The `n.*N` arm DEPARTS from a decided ruling. That is disclosed, not hidden.
//!
//! `rulings[absolute-prohibition-enforcement-stage]` (decided 2026-08-10) rules
//! enforcement **mode-dependent, uniformly**. An unconditional refusal does not
//! fit that schedule, and two tempting reconciliations are weaker than they look:
//!
//! - *"`n.*5` names no position at all, so it is a grammar matter and outside
//!   that record."* Arguable, **not conclusive**: `checklist.md:16` says a
//!   genomic reference "can not have nucleotides with additions like a `+`, `-`,
//!   or `*`", so `g.*10` is exactly parallel — and that record's census places
//!   `g.*10` under mode-gating and says lenient should newly *accept* it. A
//!   scoping argument that also carves out `g.*10` proves too much.
//! - *"The record's question names only `checklist.md`."* True of the question
//!   sentence; its rationale generalises to "every absolute prohibition", and
//!   `:54` is one.
//!
//! What does hold is narrower and empirical. The record rejected unconditional
//! refusal for **one stated reason** — it "would newly refuse inputs ferro
//! accepts today, with no escape for a caller round-tripping a real-world
//! corpus". For `n.*N` that objection is measurably empty. So the unconditional
//! arm is a **maintainer's decision** under rule 6 of the `README.md` ruleset,
//! disclosed under rule 7, and **revisitable on user demand**: a real `n.*N`
//! corpus is grounds to move it onto the `n.-N` schedule, not to defend it.
//!
//! An amendment to `rulings[absolute-prohibition-enforcement-stage]` recording
//! the carve-out is drafted and **awaiting operator sign-off**. It is not written
//! into the ledger by this change, because that record is a signed operator
//! ruling. `the_unconditional_arm_is_a_disclosed_departure` below is the
//! committed record in the meantime, so the next reader meets a disclosed
//! decision rather than a contradiction.
//!
//! # What this is NOT about — the coding axis
//!
//! `c.*N` and `c.-N` are legal HGVS and nothing here touches them. The zones
//! exist on `c.` because they are anchored to the **CDS**, not to the sequence:
//! `c.-1` is the base before the ATG and `c.*1` the base after the stop, both
//! still inside the transcript. `n.` has no CDS to anchor them to. Ferro also
//! deliberately entertains a past-the-end `c.*` in the #797 poly-A carve-out,
//! and that is untouched; `the_coding_axis_is_not_touched` below is the guard.
//!
//! # And NOT about the `r.` axis, on the spec's own words
//!
//! `numbering.md:58` — "nucleotide numbering for an RNA reference sequence
//! follows that of the associated **coding or non-coding** DNA reference
//! sequence" — makes the `r.` zone set a property of the underlying
//! **reference**, not of the prefix. `:61` gives a coding RNA reference the
//! coding numbering (so `r.*5` is legal there) and `:60` gives a non-coding one
//! the `n.` numbering (so it is not). The parser holds no provider and cannot
//! tell which it has, so refusing `r.*5` at parse would refuse a conformant
//! description. `n.` is decidable because `general.md:26` makes the prefix a
//! claim about the reference **type**.

use ferro_hgvs::error_handling::{ErrorConfig, ErrorMode};
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::parse_hgvs;

/// The transcript #255's audit used. 1786 nt, so `n.1786` is its last base.
const NR: &str = "NR_037639.1";

#[track_caller]
fn strict_refusal(input: &str) -> String {
    match parse_hgvs_with_config(input, ErrorConfig::strict()) {
        Ok(parsed) => panic!("strict accepted {input:?}, emitting {}", parsed.result),
        Err(e) => e.to_string(),
    }
}

#[track_caller]
fn strict_accepts(input: &str) {
    parse_hgvs_with_config(input, ErrorConfig::strict())
        .unwrap_or_else(|e| panic!("strict refused {input:?}: {e}"));
}

/// Refused on **every** door: the bare entry and all three modes.
#[track_caller]
fn refused_everywhere(input: &str) -> String {
    let bare = parse_hgvs(input)
        .map(|v| v.to_string())
        .expect_err(&format!("bare parse_hgvs accepted {input:?}"))
        .to_string();
    for (name, config) in [
        ("strict", ErrorConfig::strict()),
        ("lenient", ErrorConfig::lenient()),
        ("silent", ErrorConfig::silent()),
    ] {
        if let Ok(parsed) = parse_hgvs_with_config(input, config) {
            panic!("{name} accepted {input:?}, emitting {}", parsed.result);
        }
    }
    bare
}

// =============================================================================
// `n.*N` — refused at parse, in EVERY mode
// =============================================================================

/// `numbering.md:52` puts no `*` zone on the `n.` axis, so `n.*5` names no
/// position. Refused everywhere, as `E1003` — the same code `n.0` carries.
#[test]
fn a_downstream_marker_is_refused_in_every_mode() {
    for input in [
        format!("{NR}:n.*1A>G"),
        format!("{NR}:n.*5A>G"),
        format!("{NR}:n.*1_*5del"),
        format!("{NR}:n.*5+10A>G"),
    ] {
        let msg = refused_everywhere(&input);
        assert!(
            msg.contains("E1003"),
            "refusal of {input:?} should name E1003, got: {msg}"
        );
        assert!(
            msg.contains("numbering.md:52"),
            "refusal of {input:?} should cite the clause, got: {msg}"
        );
    }
}

/// `n.*N` and `n.0` are the two positions this axis can spell and does not have,
/// so they must carry one code between them rather than two. `W4008` is now
/// solely the `n.-N` diagnosis and must not appear on either.
#[test]
fn one_input_gets_one_diagnosis() {
    let star = refused_everywhere(&format!("{NR}:n.*5A>G"));
    assert!(
        !star.contains("W4008"),
        "`n.*5` must not carry the `n.-N` code as well: {star}"
    );

    // `n.0` is the pre-existing member of the class; it stays where it was.
    let zero = parse_hgvs(&format!("{NR}:n.0A>G"));
    assert!(zero.is_err(), "`n.0` must still be refused");
}

/// A range with the `*` at either end is still outside, and a cis allele must
/// not launder it past the check.
#[test]
fn a_downstream_marker_cannot_hide_behind_a_composite() {
    refused_everywhere(&format!("{NR}:n.40_*3del"));
    refused_everywhere(&format!("{NR}:n.*3_40del"));
    refused_everywhere(&format!("{NR}:n.[*5A>G;100del]"));
    refused_everywhere(&format!("{NR}:n.[100del;*5A>G]"));
    refused_everywhere(&format!("{NR}:n.[*5A>G];[100del]"));
    refused_everywhere(&format!("{NR}:n.100_*1insAAA"));
}

/// The refusal reaches the selector form too, where the transcript is only the
/// numbering selector — `:52` binds on the prefix, not on the reference.
#[test]
fn the_selector_form_is_refused_too() {
    refused_everywhere("NG_007485.1(NR_003529.3):n.*40000del");
}

// =============================================================================
// `n.-N` — refused at STRICT parse only. The five ClinVar rows depend on this.
// =============================================================================

/// Same clause, the other end — but the stage is different, and the five rows
/// below are the whole reason. `n.-1` is a nucleotide before the first, which
/// `:52`'s enumeration does not reach.
#[test]
fn an_upstream_marker_is_refused_at_strict_parse() {
    for input in [
        format!("{NR}:n.-1A>G"),
        format!("{NR}:n.-5A>G"),
        format!("{NR}:n.-100_-50del"),
        format!("{NR}:n.-3_5del"),
    ] {
        let msg = strict_refusal(&input);
        assert!(
            msg.contains("W4008"),
            "refusal of {input:?} should name W4008, got: {msg}"
        );
        assert!(
            msg.contains("numbering.md:52"),
            "refusal of {input:?} should cite the clause, got: {msg}"
        );
    }
}

/// A cis allele must not launder the prohibited member past the strict check.
#[test]
fn an_upstream_allele_member_cannot_hide_the_marker() {
    strict_refusal(&format!("{NR}:n.[100del;-5A>G]"));
    strict_refusal(&format!("{NR}:n.[-5A>G;100del]"));
}

/// **The required check.** These are the five real `n.-N` rows in ferro's
/// committed corpora, and they are the reason `n.-N` is not refused
/// unconditionally. All four corpus consumers, and every `ferro` subcommand,
/// parse through the bare `parse_hgvs` entry — so this is the door that matters.
///
/// If this test goes red, the conformance figures moved and five descriptions
/// NCBI publishes today stopped parsing.
#[test]
fn the_five_real_clinvar_rows_still_parse_on_the_bare_entry() {
    for input in [
        "NR_003051.3:n.-57T>C",
        "NR_003051.3:n.-30_-7dup",
        "NR_029595.1:n.-4771G>T",
        "LRG_163t1:n.-5delins17",
        "NR_033294.1:n.-6G>A",
    ] {
        let parsed =
            parse_hgvs(input).unwrap_or_else(|e| panic!("real ClinVar row {input:?} failed: {e}"));
        assert_eq!(
            parsed.to_string(),
            input,
            "{input:?} must round-trip as authored"
        );
    }
}

/// …and on the two permissive modes, which is what a caller validating a real
/// corpus actually uses. Lenient says so with `W4008`; silent is quiet.
#[test]
fn the_five_real_clinvar_rows_still_parse_in_lenient_and_silent() {
    for input in [
        "NR_003051.3:n.-57T>C",
        "NR_003051.3:n.-30_-7dup",
        "NR_029595.1:n.-4771G>T",
        "LRG_163t1:n.-5delins17",
        "NR_033294.1:n.-6G>A",
    ] {
        let lenient = parse_hgvs_with_config(input, ErrorConfig::lenient())
            .unwrap_or_else(|e| panic!("lenient refused {input:?}: {e}"));
        assert_eq!(lenient.result.to_string(), input);
        assert!(
            lenient
                .warnings
                .iter()
                .any(|w| w.error_type.code() == "W4008"),
            "lenient should warn W4008 for {input:?}, got {:?}",
            lenient.warnings
        );

        let silent = parse_hgvs_with_config(input, ErrorConfig::silent())
            .unwrap_or_else(|e| panic!("silent refused {input:?}: {e}"));
        assert_eq!(silent.result.to_string(), input);
        assert!(
            !silent
                .warnings
                .iter()
                .any(|w| w.error_type.code() == "W4008"),
            "silent should emit no message for {input:?}, got {:?}",
            silent.warnings
        );
    }
    assert_eq!(ErrorConfig::silent().mode, ErrorMode::Silent);
}

// =============================================================================
// The message defect: `:54`'s clause is conditional, and the selector form
// does not meet its condition
// =============================================================================

/// `numbering.md:54` prohibits describing a nucleotide beyond a transcript's
/// boundaries **"using that transcript reference sequence"** — the clause is
/// conditioned on its own subject.
///
/// For `NG_007485.1(NR_003529.3):n.-40000del` the reference sequence is the
/// genomic `NG_007485.1`, which *does* contain that nucleotide; the transcript
/// is only the numbering selector. So `:54` is not the operative clause there
/// and the message must not assert it. The refusal still stands, on `:52` alone
/// — that axis has no zone to number `-40000` from.
#[test]
fn the_boundary_clause_is_not_asserted_for_the_selector_form() {
    let bare = strict_refusal("NR_003529.3:n.-40000del");
    assert!(
        bare.contains(":54"),
        "on a bare transcript `:54` IS operative and should be cited: {bare}"
    );

    let selector = strict_refusal("NG_007485.1(NR_003529.3):n.-40000del");
    assert!(
        !selector.contains(":54"),
        "the genomic reference contains the nucleotide, so `:54`'s \"using that \
         transcript reference sequence\" is not the operative clause: {selector}"
    );
    assert!(
        selector.contains("numbering.md:52"),
        "the refusal still stands on `:52`: {selector}"
    );
}

// =============================================================================
// What must keep parsing — the axis the spec DOES describe
// =============================================================================

/// `:52`'s own enumeration: `n.1` to the last nucleotide. `NR_024540.1` is
/// 1786 nt, so `n.1786` is a real last base rather than an invented one.
#[test]
fn in_transcript_positions_are_untouched() {
    strict_accepts("NR_024540.1:n.1A>G");
    strict_accepts("NR_024540.1:n.1786A>G");
    strict_accepts(&format!("{NR}:n.100_200del"));
}

/// `:53` grants introns on this axis **explicitly** — "numbered as for coding
/// DNA reference sequences … although proceeded by `n.`". An offset is a zone
/// the spec added; refusing it would be the opposite error to the one this
/// change fixes. Note the negative offset in particular: the `-3` in `n.6-3` is
/// a distance from position 6, not a zone, and the predicate must not read it
/// as one.
#[test]
fn intronic_offsets_stay_legal_on_the_noncoding_axis() {
    for input in [
        "NG_012337.1(NR_037639.1):n.5+3A>G",
        "NG_012337.1(NR_037639.1):n.6-3A>G",
        "NG_012337.1(NR_037639.1):n.100+10del",
        "NG_012337.1(NR_037639.1):n.100-5_100-1del",
    ] {
        strict_accepts(input);
        parse_hgvs(input).unwrap_or_else(|e| panic!("bare parse refused {input:?}: {e}"));
    }
}

/// The coding axis keeps both zones. `NM_000492.4` (CFTR) is 6070 nt with CDS
/// `[70, 4513)` 0-based, so `c.*1557` is its last 3'UTR base — a real
/// coordinate, not a synthetic one.
#[test]
fn the_coding_axis_is_not_touched() {
    for input in [
        "NM_000492.4:c.*1A>G",
        "NM_000492.4:c.-1A>G",
        "NM_000492.4:c.*1557A>G",
        "NM_000492.4:c.*1_*10del",
        "NM_000492.4:c.-20_-1del",
    ] {
        strict_accepts(input);
        parse_hgvs(input).unwrap_or_else(|e| panic!("bare parse refused {input:?}: {e}"));
    }
}

/// The `r.` axis is left alone, and this test is the record of that being a
/// decision. See the module docs: `numbering.md:58` makes `r.`'s zone set a
/// property of the underlying reference, which the parser cannot resolve.
///
/// If a future change gives the parser reference context and enforces this on
/// `r.`, this test is the one to revisit — deliberately, rather than by
/// discovering it went red.
#[test]
fn the_rna_axis_is_left_alone() {
    strict_accepts("NM_003002.4:r.*5a>g");
    strict_accepts("NM_003002.4:r.-5a>g");
    parse_hgvs("NM_003002.4:r.*5del").expect("`r.*5` is conformant on a coding reference");
}

// =============================================================================
// The adjudication record
// =============================================================================

/// **Adjudication record — maintainer decision, disclosed, revisitable.**
///
/// Question: at which stage is an `n.`-axis out-of-zone marker refused?
///
/// Ruling: **the two markers are refused at different stages, on measured cost
/// rather than on a difference in the clause reading.**
///
/// - `n.*N` — refused at parse in **every** mode, as `E1003`. Authority for the
///   *invalidity* is `background/numbering.md:52`; authority for the *stage* is
///   the maintainer, under rule 6 of the `README.md` ruleset and disclosed under
///   rule 7. This **departs** from the decided
///   `rulings[absolute-prohibition-enforcement-stage]`, which makes enforcement
///   mode-dependent uniformly. The departure is justified empirically and not by
///   scoping: that record rejected unconditional refusal because it "would newly
///   refuse inputs ferro accepts today, with no escape for a caller
///   round-tripping a real-world corpus", and for this shape there is no such
///   caller — **0 of 103,762** committed `n.`-axis corpus rows. An amendment to
///   that record is drafted and awaiting operator sign-off.
/// - `n.-N` — refused at **strict** parse only, as `W4008`, exactly as that
///   ruling schedules. **5 of 103,762** rows use it and NCBI publishes all five.
///
/// The scoping argument that would have avoided the departure — "`n.*5` names no
/// position, so it is a grammar matter outside that record's reach" — is
/// recorded here as **considered and not relied on**: `checklist.md:16` makes
/// `g.*10` exactly parallel, and that record's own census mode-gates it. An
/// argument that also carves out `g.*10` proves too much.
///
/// This test asserts the *shape* of the decision, so it fails if either arm is
/// quietly moved to the other's schedule.
#[test]
fn the_unconditional_arm_is_a_disclosed_departure() {
    // `n.*N`: refused on every door, including the two the ruling says must
    // accept a non-conformant input.
    for mode in [ErrorConfig::lenient(), ErrorConfig::silent()] {
        assert!(
            parse_hgvs_with_config(&format!("{NR}:n.*5A>G"), mode).is_err(),
            "the `n.*N` arm is deliberately unconditional — see this test's docs"
        );
    }

    // `n.-N`: on the ruling's schedule, unchanged.
    assert!(parse_hgvs_with_config(&format!("{NR}:n.-5A>G"), ErrorConfig::strict()).is_err());
    assert!(parse_hgvs_with_config(&format!("{NR}:n.-5A>G"), ErrorConfig::lenient()).is_ok());
    assert!(parse_hgvs_with_config(&format!("{NR}:n.-5A>G"), ErrorConfig::silent()).is_ok());
}

/// The bare entry applies no `ErrorConfig` (#1632), so it sits on none of the
/// three mode arms — which is why the `n.-N` half is invisible to it and the
/// `n.*N` half is not. Asserted rather than reasoned about, because this is the
/// blast-radius claim the whole change rests on.
#[test]
fn the_bare_parse_entry_sees_only_the_unconditional_arm() {
    assert!(
        parse_hgvs(&format!("{NR}:n.*5A>G")).is_err(),
        "the unconditional arm reaches the bare entry"
    );
    for input in [
        format!("{NR}:n.-5A>G"),
        format!("{NR}:n.-3_5del"),
        format!("{NR}:n.-1_1insAAA"),
    ] {
        let parsed = parse_hgvs(&input).unwrap_or_else(|e| panic!("parse {input:?}: {e}"));
        assert_eq!(parsed.to_string(), input);
    }
}
