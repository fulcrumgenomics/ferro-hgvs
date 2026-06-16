//! Parser-behavior tests for HGVS input hygiene (F1) and explicit
//! spec-rejected forms (F2).
//!
//! Source: phase4_final_triage.md families F1 (input hygiene — whitespace,
//! keyword case) and F2 (spec-explicit rejections — `dupins`, `dupSEQ`,
//! `dupN`, single-position `ins`, `^` separator, pre-2016 `del3` size
//! suffix).
//!
//! Each test pins ferro's CURRENT behavior at the `parse_hgvs` API level.
//! For forms the spec rejects but ferro accepts, the test documents the
//! divergence so any future tightening is intentional. For forms ferro
//! already rejects, the test locks the rejection so a regression to
//! accept-and-warn is caught.
//!
//! Policy decisions in `phase4_final_triage.md`:
//! - D1 — `dupT`/`dupSEQ`/`dupN` (soft) are currently ACCEPTED by ferro
//!   at the API level; spec says reject. Tracked as divergence; lenient
//!   mode may later accept these soft forms with a W code. The `dupins`
//!   (absolute) form is now REJECTED outright per `DNA/duplication.md:92`
//!   (#445) — its rejection is pinned in `tests/reject_dupins.rs`, so it
//!   is intentionally not duplicated here.
//! - D2 — whitespace inside variant body rejected per `general.md:96`;
//!   ferro currently accepts spaces inside allele brackets (stripped on
//!   Display).

use ferro_hgvs::parse_hgvs;

fn assert_rejected(input: &str, why: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_err(),
        "ferro should reject {input:?} ({why}); got: {:?}",
        result.map(|v| v.to_string())
    );
}

fn assert_accepted_diverges_from_spec(input: &str, displays_as: &str, spec_says: &str) {
    let parsed = parse_hgvs(input)
        .unwrap_or_else(|e| panic!("ferro currently accepts {input:?} but failed: {e}"));
    assert_eq!(
        parsed.to_string(),
        displays_as,
        "Pinning current divergent behavior. Spec says: {spec_says}"
    );
}

// =====================================================================
// F1 — Input hygiene (whitespace + keyword case)
// =====================================================================

// Whitespace inside the variant body — ferro rejects these (good; matches
// general.md:96).

#[test]
fn rejects_whitespace_around_substitution_operator() {
    // dna-substitution:S3
    assert_rejected("NM_004006.2:c. 76 A > G", "spaces inside variant body");
}

#[test]
fn rejects_whitespace_inside_deletion() {
    // dna-deletion:S2
    assert_rejected("NM_004006.2:c.123 _ 125 del", "spaces inside del");
}

#[test]
fn rejects_whitespace_inside_duplication() {
    // dna-duplication:S16
    assert_rejected("NM_004006.2:c.20 _ 21 dup", "spaces inside dup");
}

#[test]
fn rejects_whitespace_inside_insertion() {
    // dna-insertion:S2
    assert_rejected("NM_004006.2:c.123_124 ins ATG", "spaces around ins");
}

#[test]
fn rejects_whitespace_inside_rna_delins() {
    // rna-delins:S11
    assert_rejected("NM_004006.2:r.775 _ 777 delins ga", "spaces in RNA delins");
}

#[test]
fn rejects_whitespace_inside_rna_insertion() {
    // rna-insertion:S10
    assert_rejected("NM_004006.2:r.756_757 ins uacu", "spaces around RNA ins");
}

// Whitespace inside allele bracket forms — ferro INTENTIONALLY ACCEPTS
// spaces after `;`, silently strips on Display. Diverges from
// `general.md:96` but is documented behavior (see lenient-mode comment in
// `parse_variant`). Pin the current shape.

#[test]
fn accepts_whitespace_inside_allele_bracket_dna_strips_on_display() {
    // dna-alleles:S7
    assert_accepted_diverges_from_spec(
        "NM_004006.2:c.[76A>C; 80T>G]",
        "NM_004006.2:c.[76A>C;80T>G]",
        "general.md:96 forbids; ferro's lenient parser strips on display",
    );
}

#[test]
fn accepts_whitespace_inside_allele_bracket_rna_strips_on_display() {
    // rna-alleles:S20
    assert_accepted_diverges_from_spec(
        "NM_004006.2:r.[76a>c; 80u>g]",
        "NM_004006.2:r.[76a>c;80u>g]",
        "general.md:96 forbids; ferro's lenient parser strips on display",
    );
}

#[test]
fn accepts_whitespace_inside_allele_bracket_protein_strips_on_display() {
    // protein-alleles:S20
    // The protein form must be a predicted-wrap allele (p.[(a);(b)]) per
    // protein/alleles.md. Pin behavior here.
    // ferro currently rejects this protein allele form (the bracket
    // member is not a predicted-wrap), which is consistent with
    // general.md:96 forbidding whitespace inside the variant body. Pin
    // the rejection so a regression to accept-and-strip is caught.
    assert_rejected(
        "NP_003997.1:p.[Trp24Cys; Leu25Pro]",
        "protein allele form with embedded whitespace should remain rejected",
    );
}

// Leading / trailing whitespace in the identifier — ferro INTENTIONALLY
// accepts and trims (lenient mode). Pin the divergence.

#[test]
fn accepts_leading_whitespace_in_identifier_trims_on_display() {
    // dna-substitution:S26
    assert_accepted_diverges_from_spec(
        "  NM_004006.2:c.76A>G",
        "NM_004006.2:c.76A>G",
        "general.md:96 forbids; ferro's lenient parser trims outer whitespace",
    );
}

#[test]
fn accepts_trailing_whitespace_in_identifier_trims_on_display() {
    // dna-insertion:S9
    assert_accepted_diverges_from_spec(
        "NM_004006.2:c.76A>G  ",
        "NM_004006.2:c.76A>G",
        "general.md:96 forbids; ferro's lenient parser trims outer whitespace",
    );
}

// Keyword case sensitivity — ferro rejects uppercase + title-case keywords.
// Per Sonnet's nit (v1 + v2 reviews): one cross-construct policy, not
// seven separate corners — but we keep one probe per affected construct so
// a future change is grep-able.

#[test]
fn rejects_uppercase_del_keyword() {
    // dna-deletion:S1
    assert_rejected("NM_004006.2:c.123_125DEL", "uppercase del");
    assert_rejected("NM_004006.2:c.123_125Del", "title-case Del");
}

#[test]
fn rejects_uppercase_dup_keyword() {
    // dna-duplication:S15, rna-duplication:S10/S14
    assert_rejected("NM_004006.2:c.20_21DUP", "uppercase dup");
    assert_rejected("NM_004006.2:c.20_21Dup", "title-case Dup");
}

#[test]
fn rejects_uppercase_rna_dup_keyword() {
    // rna-duplication:S10
    assert_rejected("NM_004006.2:r.20_21DUP", "uppercase dup on RNA");
}

// =====================================================================
// F2 — Spec-explicit rejections
// =====================================================================
//
// Forms the spec marks `<code class="invalid">` or explicitly says "No" /
// "is not used" / "is not allowed". ferro CURRENTLY ACCEPTS most of these
// at the parse_hgvs API level — pin the divergence.

#[test]
fn accepts_dup_size_suffix_pre_2016_diverges_from_spec() {
    // checklist:C-5 — checklist.md:49 explicitly lists `g.123del3` as
    // invalid (pre-2016 size-suffix form; must be `g.123_125del`).
    // The probe uses `del3` not `dup3` per the checklist citation.
    assert_accepted_diverges_from_spec(
        "NC_000023.11:g.123del3",
        "NC_000023.11:g.123del3",
        "checklist.md:49 — pre-2016 size-suffix form, must be g.123_125del",
    );
}

#[test]
fn accepts_dup_with_explicit_sequence_diverges_from_spec() {
    // dna-duplication:C-C18 — "the recommendation is not to" (soft).
    assert_accepted_diverges_from_spec(
        "NM_004006.2:c.20_21dupGT",
        "NM_004006.2:c.20_21dupGT",
        "DNA/duplication.md — the recommendation is not to describe the duplicated sequence",
    );
}

#[test]
fn accepts_dup_with_explicit_count_diverges_from_spec() {
    // dna-duplication:C-C19 — explicit count after `dup` is non-canonical.
    assert_accepted_diverges_from_spec(
        "NM_004006.2:c.20_21dup2",
        "NM_004006.2:c.20_21dup2",
        "DNA/duplication.md — dup<N> count suffix is non-canonical",
    );
}

// Single-position insertion (`g.123insG`) is invalid per
// insertion.md:95-101 and is now REJECTED at parse time (closes #446);
// the former accept-divergence probe is obsolete. Its rejection is
// pinned in `tests/reject_single_position_insertion.rs`, so it is
// intentionally not duplicated here.

#[test]
fn rejects_caret_separator_in_insertion() {
    // dna-insertion:C-C10 — insertion.md:103-106 Q&A "Can I use the `^`
    // character?" Answer: No, for insertions. ferro already rejects.
    assert_rejected("NC_000023.11:g.123^124insG", "^ separator in ins");
    assert_rejected("NC_000023.11:g.123^124G", "^ as anchor separator");
}
