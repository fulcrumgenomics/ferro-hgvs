//! Reject the RNA base `U`/`u` in a DNA-context edit (#486, `ENODNA`).
//!
//! DNA reference sequences (`g.`, `c.`, `n.`, `m.`, `o.`) use the alphabet
//! A/C/G/T; `U` is an RNA base and is invalid there (mutalyzer rejects with
//! `ENODNA`). `U`/`u` remains valid for RNA (`r.`) variants.

use ferro_hgvs::error::ErrorCode;
use ferro_hgvs::parse_hgvs;

fn assert_rejects(input: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_err(),
        "expected U-in-DNA {input:?} to be rejected; got: {:?}",
        result.map(|v| v.to_string())
    );
    let err = result.unwrap_err();
    assert_eq!(
        err.code(),
        Some(ErrorCode::InvalidEdit),
        "U-in-DNA rejection must carry a structured InvalidEdit code for {input:?}",
    );
}

fn assert_accepts(input: &str, expected: &str) {
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("{input:?} must parse: {e}"));
    assert_eq!(parsed.to_string(), expected, "round-trip mismatch");
}

// ---- Rejected: U/u in DNA edits across edit shapes ----

#[test]
fn rejects_substitution_alt_u() {
    assert_rejects("NM_000143.3:c.45A>U");
}

#[test]
fn rejects_delins_inserted_u() {
    assert_rejects("NM_000143.3:c.45_46delinsAUC");
}

#[test]
fn rejects_delins_deleted_u() {
    assert_rejects("NM_000143.3:c.45deluinsATC");
}

#[test]
fn rejects_ins_u() {
    assert_rejects("NM_000143.3:c.45_46insU");
}

#[test]
fn rejects_genome_substitution_u() {
    assert_rejects("NC_000001.11:g.45A>U");
}

#[test]
fn rejects_noncoding_ins_u() {
    assert_rejects("NR_004430.2:n.789_790insU");
}

#[test]
fn rejects_mito_substitution_u() {
    // `m.` (mitochondrial) is a DNA-context coordinate system.
    assert_rejects("NC_012920.1:m.45A>U");
}

#[test]
fn rejects_organellar_ins_u() {
    // `o.` (organellar) is a DNA-context coordinate system.
    assert_rejects("NC_000067.7:o.45_46insU");
}

#[test]
fn rejects_u_inside_cis_allele() {
    assert_rejects("NM_000143.3:c.[100A>G;45_46insU]");
}

// ---- Rejected: U inside tandem-repeat units (NaEdit::Repeat / MultiRepeat) ----

#[test]
fn rejects_repeat_unit_u() {
    // Repeat sequence unit carrying `U` (e.g. CUG[3]).
    assert_rejects("NM_000143.3:c.45_47CUG[3]");
}

#[test]
fn rejects_repeat_trailing_u() {
    // VEP-style trailing sequence after the count carrying `U` (e.g. [3]U).
    assert_rejects("NM_000143.3:c.45_47[3]U");
}

#[test]
fn rejects_multirepeat_unit_u() {
    // Complex multi-unit repeat whose first unit carries `U` (e.g. GU[2]GC[2]).
    assert_rejects("NM_000143.3:c.45_50GU[2]GC[2]");
}

// ---- The exact mutalyzer-corpus ENODNA cases (#326 errors axis) ----

#[test]
fn rejects_corpus_enodna_cases() {
    for input in [
        "NM_000143.3:c.45delU",
        "NM_000143.3:c.44_47delTCUG",
        "NM_000143.3:c.45deluinsATC",
        "NM_000143.3:c.44_47delTCuGinsATC",
        "NM_000143.3:c.45dupu",
        "NM_000143.3:c.44_47dupTCuG",
    ] {
        assert_rejects(input);
    }
}

// ---- Accepted: U/u valid for RNA; valid DNA bases accepted ----

#[test]
fn accepts_rna_substitution_u() {
    assert_accepts("NM_004006.2:r.456a>u", "NM_004006.2:r.456a>u");
}

#[test]
fn accepts_rna_ins_u() {
    assert_accepts("NM_004006.2:r.456_457insu", "NM_004006.2:r.456_457insu");
}

#[test]
fn accepts_dna_substitution_valid_base() {
    assert_accepts("NM_000143.3:c.45A>T", "NM_000143.3:c.45A>T");
}

#[test]
fn accepts_dna_ins_valid_base() {
    assert_accepts("NM_000143.3:c.45_46insT", "NM_000143.3:c.45_46insT");
}
