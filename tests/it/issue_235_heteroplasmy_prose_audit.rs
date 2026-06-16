//! Audit for the **remaining** work under #81 § F2 — heteroplasmy /
//! prose-shape acceptance and policy documentation.
//!
//! The F2 audit landed in #139 (`tests/mito_heteroplasmy_audit.rs`)
//! and the spec compact mosaic / chimeric form was implemented in
//! #153 (`tests/mosaic_chimeric_compact.rs`), closing #133. What is
//! **not yet pinned** in any audit:
//!
//! 1. Prose multi-allelic shapes with **non-slash separators**
//!    (`m.X>Y;T`, `m.X>Y,T`, `m.X>Y|T`).
//! 2. **Three-alt prose chains** (`m.3243A>G/T/C`).
//! 3. **IUPAC ambiguity codes on the alt** (`m.3243A>Y`, `A>R`, …) —
//!    silently accepted and round-trip cleanly. This is the closest
//!    spec-conformant shape to the prose "either-or" intent.
//! 4. **Range-edit prose** (`m.3243_3245delAGG/CTT`).
//! 5. **Mixed compact + bare-prose chains** (`m.3243=/A>G/T`).
//!
//! Policy (per #235 / #81 F2 closeout): the prose form
//! `m.<pos><ref>><alt>/<alt2>` is **rejected** and stays rejected.
//! Users have three spec-supported alternatives:
//!
//! - compound brackets `m.[3243A>G;3243A>T]`,
//! - dual fully-qualified slash `m.3243A>G/m.3243A>T`,
//! - spec compact mosaic / chimeric `m.3243=/A>G` (`=//A>G`).
//!
//! This file is **pure coverage** — no production change. Improving
//! the diagnostic message for the prose form (currently a generic
//! "Unknown variant type prefix" or "Unexpected trailing characters")
//! is a separate quality-of-error improvement and is not covered
//! here.

use ferro_hgvs::{parse_hgvs, HgvsVariant};

const ACC: &str = "NC_012920.1";

#[track_caller]
fn assert_rejects(input: &str) {
    let res = parse_hgvs(input);
    assert!(
        res.is_err(),
        "expected {input:?} to be rejected, but parsed: {:?}",
        res.ok().map(|v| format!("{}", v))
    );
}

#[track_caller]
fn assert_round_trips(input: &str) {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?} failed: {e}"));
    let out = format!("{}", v);
    assert_eq!(out, input, "round-trip mismatch for {input:?}");
}

// =============================================================================
// SECTION 1 — Prose with non-slash separators
// =============================================================================
//
// ClinVar / submitter prose sometimes uses `;`, `,`, or `|` to join
// alternative alleles. None of these are spec forms; all stay
// rejected. Pin so the rejection survives parser changes.

mod non_slash_separator_prose {
    use super::*;

    #[test]
    fn semicolon_bare_alt_is_rejected() {
        assert_rejects(&format!("{ACC}:m.3243A>G;T"));
    }

    #[test]
    fn comma_bare_alt_is_rejected() {
        assert_rejects(&format!("{ACC}:m.3243A>G,T"));
    }

    #[test]
    fn pipe_bare_alt_is_rejected() {
        assert_rejects(&format!("{ACC}:m.3243A>G|T"));
    }

    /// Trailing bare separator with no second alt — also rejected
    /// (guards against accidental "trailing comma" tolerance).
    #[test]
    fn trailing_separator_only_is_rejected() {
        assert_rejects(&format!("{ACC}:m.3243A>G;"));
        assert_rejects(&format!("{ACC}:m.3243A>G,"));
        assert_rejects(&format!("{ACC}:m.3243A>G|"));
    }
}

// =============================================================================
// SECTION 2 — Three-alt prose chains
// =============================================================================
//
// Single-alt prose `m.A>G/T` is already pinned as rejected by
// `tests/mito_heteroplasmy_audit.rs`. Pin the longer chain forms here.

mod prose_chain {
    use super::*;

    #[test]
    fn three_alts_with_slash_is_rejected() {
        assert_rejects(&format!("{ACC}:m.3243A>G/T/C"));
    }

    #[test]
    fn four_alts_with_slash_is_rejected() {
        assert_rejects(&format!("{ACC}:m.3243A>G/T/C/N"));
    }
}

// =============================================================================
// SECTION 3 — IUPAC ambiguity codes on the alt
// =============================================================================
//
// `Y` = C ∪ T, `R` = A ∪ G, `N` = any, `W` = A ∪ T, etc. Ferro
// silently accepts these on the alt and round-trips them. This is
// the closest spec-conformant shape to a prose "either-or" intent
// — the user can write `m.3243A>Y` for "C or T" without needing
// the mosaic form. Pin acceptance so this acceptance is not lost
// in a future parser tightening.

mod iupac_ambiguity_alt {
    use super::*;

    #[test]
    fn y_ambiguity_alt_round_trips() {
        assert_round_trips(&format!("{ACC}:m.3243A>Y"));
    }

    #[test]
    fn r_ambiguity_alt_round_trips() {
        assert_round_trips(&format!("{ACC}:m.3243A>R"));
    }

    #[test]
    fn n_ambiguity_alt_round_trips() {
        assert_round_trips(&format!("{ACC}:m.3243A>N"));
    }

    #[test]
    fn w_ambiguity_alt_round_trips() {
        assert_round_trips(&format!("{ACC}:m.3243A>W"));
    }

    /// Confirm the parse target is the mt variant type — the
    /// ambiguity code on the alt does not route to a different
    /// variant species.
    #[test]
    fn ambiguity_alt_parses_as_mt_variant() {
        let v = parse_hgvs(&format!("{ACC}:m.3243A>Y")).expect("parse");
        assert!(matches!(v, HgvsVariant::Mt(_)));
    }
}

// =============================================================================
// SECTION 4 — Range-edit prose
// =============================================================================
//
// Slash-tail prose on range edits (`del`, `dup`, `inv`) is rejected,
// matching the point-edit prose policy. The spec compact form
// `m.3243_3245=/del` IS accepted (#153); pin both so the policy
// is symmetric and explicit.

mod range_prose {
    use super::*;

    #[test]
    fn del_with_prose_alt_seq_is_rejected() {
        assert_rejects(&format!("{ACC}:m.3243_3245delAGG/CTT"));
    }

    #[test]
    fn dup_with_prose_alt_seq_is_rejected() {
        assert_rejects(&format!("{ACC}:m.3243_3245dupAGG/CTT"));
    }

    /// Companion: the spec compact form for del is accepted.
    /// Pin so the policy split is explicit.
    #[test]
    fn spec_compact_del_is_accepted() {
        assert_round_trips(&format!("{ACC}:m.3243_3245=/del"));
    }
}

// =============================================================================
// SECTION 5 — Mixed compact + bare-prose chains
// =============================================================================
//
// `m.3243=/A>G/T` — the LHS is the spec compact form, the RHS is
// a bare prose tail. Rejected (the tail is not a recognised
// alternative). Pin so accepting the spec compact form does not
// accidentally license the prose tail.

mod compact_plus_prose_chain {
    use super::*;

    #[test]
    fn compact_plus_bare_alt_is_rejected() {
        assert_rejects(&format!("{ACC}:m.3243=/A>G/T"));
    }

    #[test]
    fn chimeric_compact_plus_bare_alt_is_rejected() {
        assert_rejects(&format!("{ACC}:m.3243=//A>G/T"));
    }
}

// =============================================================================
// SECTION 6 — Supported alternatives to prose multi-allelic
// =============================================================================
//
// Document the three spec-supported alternatives users have for
// expressing multi-allelic mt sites. Pinned as regression guards
// — if any of these break, the F2 policy ("reject prose, here
// are the alternatives") becomes hollow.

mod supported_alternatives {
    use super::*;

    /// Cis compound brackets — the canonical way to express two
    /// alternative alleles at the same mt position.
    #[test]
    fn compound_brackets_cis_round_trips() {
        assert_round_trips(&format!("{ACC}:m.[3243A>G;3243A>T]"));
    }

    /// Dual fully-qualified slash — slash dispatcher with explicit
    /// accession on both halves.
    #[test]
    fn dual_fully_qualified_slash_round_trips() {
        assert_round_trips(&format!("{ACC}:m.3243A>G/{ACC}:m.3243A>T"));
    }

    /// Spec compact mosaic — single-alt at the same position
    /// against the reference (`=`).
    #[test]
    fn spec_compact_mosaic_round_trips() {
        assert_round_trips(&format!("{ACC}:m.3243=/A>G"));
    }

    /// Spec compact chimeric — single-alt at the same position
    /// against the reference (`=`), double slash.
    #[test]
    fn spec_compact_chimeric_round_trips() {
        assert_round_trips(&format!("{ACC}:m.3243=//A>G"));
    }
}
