//! Audit for #81 § E1 — RNA lowercase enforcement coverage across edit
//! types and input casings.
//!
//! Per HGVS v21.0 RNA nomenclature, the nucleotide alphabet is
//! lowercase (`a`, `c`, `g`, `u`). `Display` on `RnaVariant` lowercases
//! the bases regardless of input casing (`src/hgvs/edit.rs:617+`,
//! `src/hgvs/variant.rs:1269+`). Prior coverage in #72 merge tests and
//! `tests/rna_spl_marker.rs` exercises a handful of edit shapes but
//! does not exhaustively pin the lattice of edit-type × input-casing.
//!
//! This file pins, for every confirmed RNA edit shape:
//!
//! 1. lowercase input round-trips through `parse` → `Display`,
//! 2. uppercase/mixed input is silently canonicalized to lowercase on
//!    `Display` (lenient input policy),
//! 3. `T`/`t` input is retained as lowercase `t` — the parser is
//!    lenient about thymine in `r.` variants and does **not**
//!    rewrite it to `u`. This pins observed behavior; whether to
//!    reject `t` or canonicalize to `u` is a parser-policy question
//!    out of scope for this audit.
//!
//! The predicted-wrapper form (`r.(<edit>)`) is intentionally NOT
//! covered here — probing shows the parens are dropped on `Display`
//! (and `r.(=)` fails to parse). That belongs to #81 § G3 (uncertain
//! edit round-trip) and gets its own audit.

use ferro_hgvs::parse_hgvs;

const ACC: &str = "NM_000088.3";

/// Parse `input`, format via `Display`, assert it equals `expected`.
#[track_caller]
fn assert_canonicalizes_to(input: &str, expected: &str) {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?} failed: {e}"));
    let out = format!("{}", v);
    assert_eq!(out, expected, "input={input:?}");
}

/// Round-trip helper: parse `s`, `Display`, assert equality.
#[track_caller]
fn assert_round_trips(s: &str) {
    assert_canonicalizes_to(s, s);
}

// =============================================================================
// SECTION 1 — Substitution (sub)
// =============================================================================
//
// `r.<pos><ref>><alt>` and `r.<pos>><alt>` (no-ref).

mod substitution {
    use super::*;

    #[test]
    fn lowercase_round_trips() {
        assert_round_trips(&format!("{ACC}:r.123a>g"));
        assert_round_trips(&format!("{ACC}:r.123c>u"));
    }

    #[test]
    fn uppercase_canonicalizes_to_lowercase() {
        assert_canonicalizes_to(&format!("{ACC}:r.123A>G"), &format!("{ACC}:r.123a>g"));
        assert_canonicalizes_to(&format!("{ACC}:r.123C>U"), &format!("{ACC}:r.123c>u"));
    }

    #[test]
    fn mixed_case_canonicalizes_to_lowercase() {
        assert_canonicalizes_to(&format!("{ACC}:r.123A>g"), &format!("{ACC}:r.123a>g"));
        assert_canonicalizes_to(&format!("{ACC}:r.123a>G"), &format!("{ACC}:r.123a>g"));
    }

    /// Thymine input is accepted and emitted as lowercase `t` — pins
    /// the lenient parser policy. The spec prefers `u`; a future
    /// parser tightening would reject `t` (or canonicalize to `u`).
    #[test]
    fn thymine_input_is_retained_as_lowercase_t() {
        assert_round_trips(&format!("{ACC}:r.123a>t"));
        assert_canonicalizes_to(&format!("{ACC}:r.123A>T"), &format!("{ACC}:r.123a>t"));
    }
}

// =============================================================================
// SECTION 2 — Deletion (del)
// =============================================================================
//
// `r.<start>_<end>del` (position-only) and
// `r.<start>_<end>del<seq>` (with stated reference bases).

mod deletion {
    use super::*;

    #[test]
    fn pos_only_round_trips() {
        assert_round_trips(&format!("{ACC}:r.123_125del"));
    }

    #[test]
    fn lowercase_stated_ref_round_trips() {
        assert_round_trips(&format!("{ACC}:r.123_125delaug"));
    }

    #[test]
    fn uppercase_stated_ref_canonicalizes_to_lowercase() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_125delAUG"),
            &format!("{ACC}:r.123_125delaug"),
        );
    }

    #[test]
    fn mixed_case_stated_ref_canonicalizes_to_lowercase() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_125delAuG"),
            &format!("{ACC}:r.123_125delaug"),
        );
    }

    #[test]
    fn thymine_in_stated_ref_is_retained_as_lowercase_t() {
        assert_round_trips(&format!("{ACC}:r.123_125delaut"));
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_125delAUT"),
            &format!("{ACC}:r.123_125delaut"),
        );
    }
}

// =============================================================================
// SECTION 3 — Insertion (ins)
// =============================================================================
//
// `r.<a>_<b>ins<seq>` — `<seq>` must be lowercase on emission.

mod insertion {
    use super::*;

    #[test]
    fn lowercase_round_trips() {
        assert_round_trips(&format!("{ACC}:r.123_124insauc"));
    }

    #[test]
    fn uppercase_inserted_seq_canonicalizes_to_lowercase() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_124insAUC"),
            &format!("{ACC}:r.123_124insauc"),
        );
    }

    #[test]
    fn mixed_case_inserted_seq_canonicalizes_to_lowercase() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_124insAuC"),
            &format!("{ACC}:r.123_124insauc"),
        );
    }
}

// =============================================================================
// SECTION 4 — Duplication (dup)
// =============================================================================
//
// `r.<start>_<end>dup` (position-only) and
// `r.<start>_<end>dup<seq>` (with stated reference bases).

mod duplication {
    use super::*;

    #[test]
    fn pos_only_round_trips() {
        assert_round_trips(&format!("{ACC}:r.123_125dup"));
    }

    #[test]
    fn lowercase_stated_ref_round_trips() {
        assert_round_trips(&format!("{ACC}:r.123_125dupaug"));
    }

    #[test]
    fn uppercase_stated_ref_canonicalizes_to_lowercase() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_125dupAUG"),
            &format!("{ACC}:r.123_125dupaug"),
        );
    }

    #[test]
    fn mixed_case_stated_ref_canonicalizes_to_lowercase() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_125dupAuG"),
            &format!("{ACC}:r.123_125dupaug"),
        );
    }
}

// =============================================================================
// SECTION 5 — Delins
// =============================================================================
//
// `r.<a>_<b>del[ref]ins<alt>` — both `ref` (if stated) and `alt`
// must be lowercase on emission.

mod delins {
    use super::*;

    #[test]
    fn lowercase_round_trips() {
        assert_round_trips(&format!("{ACC}:r.123_125delinsauc"));
    }

    #[test]
    fn uppercase_alt_canonicalizes_to_lowercase() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_125delinsAUC"),
            &format!("{ACC}:r.123_125delinsauc"),
        );
    }

    #[test]
    fn uppercase_stated_ref_and_lowercase_alt_canonicalizes() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_125delAUGinsauc"),
            &format!("{ACC}:r.123_125delauginsauc"),
        );
    }

    #[test]
    fn mixed_ref_and_mixed_alt_canonicalize_to_lowercase() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_125delAuGinsAuC"),
            &format!("{ACC}:r.123_125delauginsauc"),
        );
    }
}

// =============================================================================
// SECTION 6 — Inversion (inv)
// =============================================================================
//
// `r.<start>_<end>inv` (position-only) and
// `r.<start>_<end>inv<seq>` (with stated reference bases).

mod inversion {
    use super::*;

    #[test]
    fn pos_only_round_trips() {
        assert_round_trips(&format!("{ACC}:r.123_125inv"));
    }

    #[test]
    fn lowercase_stated_ref_round_trips() {
        assert_round_trips(&format!("{ACC}:r.123_125invaug"));
    }

    #[test]
    fn uppercase_stated_ref_canonicalizes_to_lowercase() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.123_125invAUG"),
            &format!("{ACC}:r.123_125invaug"),
        );
    }
}

// =============================================================================
// SECTION 7 — Repeat
// =============================================================================
//
// `r.<pos><unit>[N]` and multi-unit forms `r.<pos><u1>[N]<u2>[M]`.

mod repeat {
    use super::*;

    #[test]
    fn lowercase_unit_round_trips() {
        assert_round_trips(&format!("{ACC}:r.123aug[5]"));
    }

    #[test]
    fn uppercase_unit_canonicalizes_to_lowercase() {
        assert_canonicalizes_to(&format!("{ACC}:r.123AUG[5]"), &format!("{ACC}:r.123aug[5]"));
    }

    #[test]
    fn mixed_case_unit_canonicalizes_to_lowercase() {
        assert_canonicalizes_to(&format!("{ACC}:r.123Aug[5]"), &format!("{ACC}:r.123aug[5]"));
    }
}

// =============================================================================
// SECTION 8 — Whole-entity / sentinel forms
// =============================================================================
//
// `r.=`, `r.?`, `r.0`, and `r.spl?` are alphabet-free and must
// round-trip verbatim. Pin as regression guards so a fix to the
// uncertain-edit path (#81 G3) does not regress these.

mod whole_entity {
    use super::*;

    #[test]
    fn identity_round_trips() {
        assert_round_trips(&format!("{ACC}:r.="));
    }

    #[test]
    fn unknown_round_trips() {
        assert_round_trips(&format!("{ACC}:r.?"));
    }

    #[test]
    fn no_product_round_trips() {
        assert_round_trips(&format!("{ACC}:r.0"));
    }

    #[test]
    fn spl_unknown_round_trips() {
        assert_round_trips(&format!("{ACC}:r.spl?"));
    }
}

// =============================================================================
// SECTION 9 — Compound alleles (cis brackets)
// =============================================================================
//
// Lowercase enforcement inside `r.[...;...]` brackets — pins that
// each member edit's bases are lowercase regardless of input casing.

mod compound_cis {
    use super::*;

    #[test]
    fn two_sub_lowercase_round_trips() {
        assert_round_trips(&format!("{ACC}:r.[123a>g;125c>u]"));
    }

    #[test]
    fn two_sub_uppercase_canonicalizes_to_lowercase() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.[123A>G;125C>U]"),
            &format!("{ACC}:r.[123a>g;125c>u]"),
        );
    }

    #[test]
    fn mixed_edit_shapes_canonicalize_to_lowercase() {
        assert_canonicalizes_to(
            &format!("{ACC}:r.[123A>G;200_202delAUG;300_301insAUC]"),
            &format!("{ACC}:r.[123a>g;200_202delaug;300_301insauc]"),
        );
    }
}
