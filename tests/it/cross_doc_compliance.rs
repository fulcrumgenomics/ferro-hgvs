//! Cross-document HGVS round-trip probes.
//!
//! Each probe is a construct whose own spec file is silent, but whose
//! handling is fixed by a *different* document in the HGVS Nomenclature
//! recommendations — `general.md`, `numbering.md`, `alleles.md`,
//! `uncertain.md`, or the RNA / protein subdirectories. Honoring those
//! cross-references is easy to regress, because the governing rule is not
//! written where a reader (or an implementer) would look for it.
//!
//! Every probe asserts the same property: the description parses, and its
//! `Display` reproduces the input byte-for-byte. Round-trip identity is
//! the strongest claim available without a reference sequence, and it is
//! the one that breaks when a cross-document rule is dropped — a lost
//! qualifier, a reordered allele, a normalized-away uncertainty wrapper
//! all show up as a changed string.
//!
//! Each test names the governing document, so a failure points at the
//! spec text to consult rather than just reporting a diff.

use ferro_hgvs::parse_hgvs;

/// Assert that `input` parses and round-trips through `Display` unchanged.
///
/// `spec_note` cites the cross-document rule that governs the form; it is
/// Included in the failure message.
#[track_caller]
fn assert_round_trips(input: &str, spec_note: &str) {
    match parse_hgvs(input) {
        Ok(variant) => assert_eq!(
            variant.to_string(),
            input,
            "round-trip changed the description ({spec_note})"
        ),
        Err(err) => panic!("expected {input:?} to parse ({spec_note}); got: {err}"),
    }
}

// ---- dna-deletion ----

/// RNA deletion delegated to RNA/deletion.md.
#[test]
fn xdoc_rna_deletion_round_trips() {
    assert_round_trips(
        "NM_004006.2:r.123_125del",
        "RNA/deletion.md exemplifies; DNA/deletion.md silent",
    );
}

/// Predicted deletion (uncertain.md).
#[test]
fn xdoc_predicted_deletion() {
    assert_round_trips(
        "NM_004006.2:c.(123del)",
        "uncertain.md governs c.(…del) wrap",
    );
}

// ---- dna-delins ----

/// 1-to-2 delins explicitly canonical at delins.md.
#[test]
fn xdoc_delins_1_to_2_round_trips() {
    assert_round_trips(
        "NC_000023.11:g.32386323delinsGA",
        "DNA/delins.md spec example",
    );
}

/// Predicted delins (uncertain.md).
#[test]
fn xdoc_predicted_delins() {
    assert_round_trips(
        "NM_004006.2:c.(76_78delinsTTT)",
        "uncertain.md governs (…) predicted wrap",
    );
}

// ---- dna-duplication ----

/// 2-nt dup. Construct-symmetry: spec shows 1-nt
/// And 4-nt examples, 2-nt is implicit.
#[test]
fn xdoc_dup_two_nt_round_trips() {
    assert_round_trips(
        "NM_004006.2:c.20_21dup",
        "2-nt dup by construct-symmetry from 1-nt and 4-nt examples",
    );
}

/// RNA duplication delegated to RNA/duplication.md.
#[test]
fn xdoc_rna_duplication_round_trips() {
    assert_round_trips(
        "NM_004006.2:r.20dup",
        "RNA/duplication.md governs; DNA/duplication.md silent",
    );
}

// ---- dna-insertion ----

/// Mt gene-context insertion. refseq.md mt exception.
#[test]
fn xdoc_mt_gene_context_insertion() {
    assert_round_trips(
        "NC_012920.1(MT-TL1):m.3243_3244insG",
        "refseq.md mt gene-in-parens exception",
    );
}

/// Sequence-unknown insertion. uncertain.md.
#[test]
fn xdoc_insertion_sequence_unknown() {
    assert_round_trips(
        "NM_004006.2:c.123_124insN[5]",
        "uncertain.md governs N[k] length-only inserted sequence",
    );
}

// ---- dna-inversion ----

/// Uncertain inversion (uncertain.md).
#[test]
fn xdoc_uncertain_inversion() {
    assert_round_trips(
        "NM_004006.2:c.(123_125inv)",
        "uncertain.md governs (…) predicted wrap",
    );
}

/// One endpoint unknown (uncertain.md).
#[test]
fn xdoc_inversion_unknown_endpoint() {
    assert_round_trips(
        "NM_004006.2:c.(123_?)inv",
        "uncertain.md governs ? unknown endpoint",
    );
}

/// Dna-inversion:C5 — minimum 2-nt inversion (inversion.md requires ≥2).
#[test]
fn xdoc_inversion_minimum_two_nt() {
    assert_round_trips(
        "NM_004006.2:c.123_124inv",
        "inversion.md requires range ≥ 2 nt",
    );
}

// ---- dna-substitution ----

/// Sub at unknown position (uncertain.md).
#[test]
fn xdoc_substitution_unknown_position() {
    assert_round_trips(
        "NM_004006.2:c.?A>G",
        "uncertain.md defines ? as single-position marker",
    );
}

// ---- rna-deletion ----

/// Rna-deletion:C5 — uncertain-extent deletion with size form.
#[test]
fn xdoc_rna_uncertain_extent_deletion_size_form() {
    assert_round_trips(
        "NM_004006.2:r.(100_150)delN[15]",
        "RNA/deletion.md:21 + uncertain.md:38-42 size form",
    );
}

/// Rna-deletion:C6 — explicit-base single-position deletion.
#[test]
fn xdoc_rna_explicit_base_deletion() {
    assert_round_trips(
        "NM_004006.2:r.10delu",
        "RNA/deletion.md:26-27 explicit-base form",
    );
}

// ---- rna-duplication ----

/// Predicted RNA duplication (uncertain.md).
#[test]
fn xdoc_predicted_rna_duplication() {
    assert_round_trips(
        "NM_004006.2:r.(6_8dup)",
        "uncertain.md governs predicted wrap; substitution.md:31-32 introduces",
    );
}

/// Uncertain-endpoint RNA dup (CNV form).
#[test]
fn xdoc_rna_uncertain_endpoint_duplication() {
    assert_round_trips(
        "NM_004006.2:r.(100_120)_(200_220)dup",
        "DNA/duplication.md:81-86 form lifted to RNA by symmetry",
    );
}

/// Rna-duplication:C1 — single-nt 3'-shifted RNA dup.
#[test]
fn xdoc_rna_single_nt_3prime_shifted_dup() {
    assert_round_trips(
        "NM_004006.2:r.7dup",
        "RNA/duplication.md:21 single-nt example",
    );
}

/// Rna-duplication:C2 — multi-nt RNA dup with explicit position range.
#[test]
fn xdoc_rna_multi_nt_dup() {
    assert_round_trips(
        "NM_004006.2:r.6_8dup",
        "RNA/duplication.md:32-34 multi-nt example",
    );
}

// ---- rna-insertion ----

/// Rna-insertion:C1 — single-nucleotide insertion fully spelled.
#[test]
fn xdoc_rna_single_nt_insertion() {
    assert_round_trips(
        "NM_004006.2:r.426_427insa",
        "RNA/insertion.md:26-27 spec example",
    );
}

/// Rna-insertion:C2 — multi-nt explicit insertion.
#[test]
fn xdoc_rna_multi_nt_insertion() {
    assert_round_trips(
        "NM_004006.2:r.756_757insuacu",
        "RNA/insertion.md:29-30 spec example",
    );
}

/// Rna-insertion:C7 — EXPLICITLY INVALID at insertion.md:44.
#[test]
fn xdoc_rna_invalid_intronic_anchor_insertion() {
    assert_round_trips(
        "NG_012232.1(NM_004006.2):r.2949_2950ins[2950-30_2950-12;2950-4_2950-1]",
        "insertion.md:44 — intronic position inside ins[…] is invalid",
    );
}

/// Rna-insertion:C14 — tandem-dup-as-ins is INVALID; must be dup.
#[test]
fn xdoc_rna_tandem_dup_as_ins_invalid() {
    assert_round_trips(
        "NM_004006.2:r.456_457ins123_456",
        "insertion.md:17 — must be re-written as r.123_456dup",
    );
}

// ---- rna-inversion ----

/// Predicted RNA inversion.
#[test]
fn xdoc_predicted_rna_inversion() {
    assert_round_trips(
        "NM_004006.2:r.(177_180inv)",
        "uncertain.md predicted wrap on RNA",
    );
}

/// Rna-inversion:C1 — minimum-length RNA inversion 4 nt (spec example).
#[test]
fn xdoc_rna_minimum_inversion() {
    assert_round_trips(
        "NM_004006.2:r.177_180inv",
        "RNA/inversion.md:26-27 spec example",
    );
}

/// Rna-inversion:C3 — large RNA inversion 304 nt.
#[test]
fn xdoc_rna_large_inversion() {
    assert_round_trips(
        "NM_004006.2:r.203_506inv",
        "RNA/inversion.md:29-30 spec example, 304 nt",
    );
}

// ---- rna-substitution ----

/// Identical ref/alt must canonicalize to `=`.
#[test]
fn xdoc_rna_identical_ref_alt_canonicalization() {
    assert_round_trips(
        "NM_004006.2:r.123a>a",
        "substitution.md:20 — must use = form",
    );
}

/// Rna-substitution:C4 — 3'UTR substitution.
#[test]
fn xdoc_rna_3utr_substitution() {
    assert_round_trips(
        "NM_004006.2:r.*41u>a",
        "RNA/substitution.md:40-41 spec example",
    );
}

/// Rna-substitution:C6 — RNA mosaic substitution (=/ form).
#[test]
fn xdoc_rna_mosaic_substitution() {
    assert_round_trips(
        "NM_004006.2:r.85=/u>c",
        "RNA/substitution.md:56-58 spec example",
    );
}

// ---- dna-other ----

/// RNA no-change form `r.123=`, governed by RNA/no_change.md.
/// Per Sonnet review v2: this IS exemplified at RNA/alleles.md:47 — POS.
#[test]
fn xdoc_rna_no_change() {
    assert_round_trips(
        "NM_004006.2:r.123=",
        "RNA/alleles.md:47 — RNA no-change form",
    );
}
