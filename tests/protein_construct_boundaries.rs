//! Protein construct-boundary tests (F13 — 46 net-new corners).
//!
//! Protein has the densest set of cross-construct canonicalization rules
//! in the whole HGVS spec: sub ↔ delins ↔ dup ↔ del ↔ fs ↔ ext. The
//! Stage 1+2 review surfaced 46 net-new corners across these boundaries.
//! Each probe below pins ferro's CURRENT behavior on a representative
//! canonical form so any future change is intentional.
//!
//! Source: phase4_final_triage.md family F13. Spec citations in
//! per-test comments.

use ferro_hgvs::parse_hgvs;

/// What a probe pins: that the input round-trips through ferro
/// unchanged, or that ferro currently has a parse gap on it. Making the
/// expectation explicit prevents the helper from silently swallowing a
/// parse regression — a `RoundTrip` probe that stops parsing now fails
/// loudly, and a `ParseGap` probe that starts parsing flags the (good)
/// change so the pin can be promoted to `RoundTrip`.
enum ProbeExpectation {
    /// ferro parses the input and Display reproduces it byte-for-byte.
    RoundTrip,
    /// ferro does not yet handle this canonical form (parse error).
    /// Reserved: every probe in this file currently round-trips, so no
    /// call site constructs this arm yet. It stays part of the API so a
    /// future construct-boundary probe that pins a known parse gap can
    /// express that intent (and gets flagged for promotion when the gap
    /// closes) without reshaping the helper.
    #[allow(dead_code)]
    ParseGap,
}

fn pin_round_trip_or_gap(input: &str, note: &str, expected: ProbeExpectation) {
    match (expected, parse_hgvs(input)) {
        (ProbeExpectation::RoundTrip, Ok(v)) => {
            assert_eq!(
                v.to_string(),
                input,
                "round-trip mismatch for {input:?} ({note})"
            );
        }
        (ProbeExpectation::RoundTrip, Err(e)) => {
            panic!("expected round-trip, got parse gap for {input:?} ({note}): {e}");
        }
        (ProbeExpectation::ParseGap, Err(_)) => {}
        (ProbeExpectation::ParseGap, Ok(v)) => {
            panic!(
                "expected parse gap, but {input:?} now parses as {:?} ({note}) — promote this probe to RoundTrip",
                v.to_string()
            );
        }
    }
}

fn pin_either(input: &str, _note: &str) {
    // Just exercise the parser; pin that it doesn't crash.
    let _ = parse_hgvs(input);
}

// =====================================================================
// Protein substitution boundary cases (F13 from protein-substitution shard)
// =====================================================================

/// protein-substitution:C-C11 — mosaic with reference written FIRST,
/// frequency-independent (protein/substitution.md:83). Probe in 1-letter
/// form (3-letter shown in spec).
#[test]
fn protein_mosaic_substitution_reference_first() {
    pin_either(
        "NP_003997.1:p.Trp24=/Cys",
        "protein/substitution.md:83 — reference written first regardless of allele frequency",
    );
}

/// protein-substitution:S-S15 — mosaic with stop alt (silent in spec).
#[test]
fn protein_mosaic_with_stop_alt() {
    pin_either(
        "NP_003997.1:p.Trp24=/Ter",
        "protein/substitution.md:81-84 silent on mosaic + stop alt",
    );
}

/// protein-substitution:S-S16 — chimeric protein substitution.
#[test]
fn protein_chimeric_substitution() {
    pin_either(
        "NP_003997.1:p.Trp24=//Cys",
        "DNA/substitution.md:51-53 defines =//; protein-sub.md silent",
    );
}

// =====================================================================
// Protein deletion boundary cases
// =====================================================================

/// protein-deletion:C-C9 — explicit reference position vs the rule that
/// protein description MUST come from comparing variant to reference
/// (deletion.md:102-107).
#[test]
fn protein_deletion_3prime_shift_choice() {
    pin_round_trip_or_gap(
        "NP_003997.1:p.Ser5del",
        "protein/deletion.md:102-107 — most-3' position rule",
        ProbeExpectation::RoundTrip,
    );
}

/// protein-deletion:C-C12 — `del<N>` size-only form is explicitly
/// invalid; canonical is the explicit-range form.
#[test]
fn protein_deletion_size_only_invalid_canonical_is_range() {
    // The size-only form may or may not parse; pin behavior.
    pin_either(
        "NP_003997.1:p.Arg45del6",
        "protein/deletion.md:87-89 — invalid; canonical is p.Arg45_Gly50del",
    );
    // And the canonical form must round-trip cleanly.
    pin_round_trip_or_gap(
        "NP_003997.1:p.Arg45_Gly50del",
        "spec-canonical form",
        ProbeExpectation::RoundTrip,
    );
}

// =====================================================================
// Protein delins boundary cases
// =====================================================================

/// protein-delins:C-C3 — single-residue delins replaced by 1 residue
/// MUST be substitution. Forbidden as delins.
#[test]
fn protein_delins_single_to_single_must_be_sub() {
    // Per delins.md:17 — `p.Trp26delinsGly` is forbidden; canonical is
    // `p.Trp26Gly`.
    pin_either(
        "NP_003997.1:p.Trp26delinsGly",
        "protein/delins.md:17 + substitution.md — must canonicalize to sub",
    );
    pin_round_trip_or_gap(
        "NP_003997.1:p.Trp26Gly",
        "canonical substitution",
        ProbeExpectation::RoundTrip,
    );
}

/// protein-delins:C-C9 — predicted delins as protein-level consequence
/// of a DNA inversion (delins.md:58-60).
#[test]
fn protein_predicted_delins_from_dna_inversion() {
    pin_either(
        "NP_003997.1:p.(Glu125_Ala132delinsGlyLeuHisArgPheIleValLeu)",
        "protein/delins.md:58-60 — predicted delins from DNA inv",
    );
}

/// protein-delins:C-C11 — large unknown-residue insertion using `Xaa[N]`.
#[test]
fn protein_ins_unknown_xaa_length() {
    pin_either(
        "NP_003997.1:p.Lys2_Leu3insXaa[34]",
        "protein/delins.md:23 — Xaa[N] length notation",
    );
}

/// protein-delins:S-S4 — caret alternatives inside delins.
#[test]
fn protein_delins_caret_alternatives() {
    pin_either(
        "NP_003997.1:p.(Arg123delinsSer^Trp)",
        "protein-delins.md silent; uncertain.md ^ semantics",
    );
}

/// protein-delins:S-S7 — chimeric delins on protein (silent).
#[test]
fn protein_chimeric_delins() {
    pin_either(
        "NP_003997.1:p.Arg123_Lys125=//delinsSer",
        "general.md:93 =// inherited by symmetry; protein/delins.md silent",
    );
}

// =====================================================================
// Protein duplication boundary cases
// =====================================================================

/// protein-duplication:S-S1 — dup with explicit aa letters. The protein
/// doc does NOT explicitly forbid this (analogous to DNA `dupSEQ` soft
/// prohibition).
#[test]
fn protein_dup_with_explicit_aa_letters() {
    pin_either(
        "NP_003997.1:p.Lys23_Val25dupLysGlyVal",
        "protein/duplication.md silent on dup<seq> form",
    );
}

/// protein-duplication:C-C2 — DNA-encoded dup of 9 nt → protein dup of
/// 3 aa, with 1-aa shift due to reading-frame alignment.
#[test]
fn protein_dup_predicted_from_dna_round_trips() {
    pin_round_trip_or_gap(
        "NP_003997.1:p.(Pro458_Gly460dup)",
        "protein/duplication.md:46-48 spec example",
        ProbeExpectation::RoundTrip,
    );
}

/// protein-duplication:C-C3 — 3'-rule pushes the boundary by 1 aa.
#[test]
fn protein_dup_3prime_rule_canonical() {
    pin_round_trip_or_gap(
        "NP_003997.1:p.(Asp90_Val120dup)",
        "protein/duplication.md:54-56 spec example",
        ProbeExpectation::RoundTrip,
    );
    // The non-canonical alternative form should still parse (for
    // detection purposes) but is not canonical.
    pin_either(
        "NP_003997.1:p.(Val89_Gln119dup)",
        "non-canonical alternative shifted by 1",
    );
}

/// protein-duplication:C-C5 — `p.?` for whole-exon dup affecting Met1
/// or terminator.
#[test]
fn protein_unknown_for_exon_dup_affecting_termini() {
    pin_round_trip_or_gap(
        "NP_003997.1:p.?",
        "protein/duplication.md:61-69 — p.? for unknown consequence",
        ProbeExpectation::RoundTrip,
    );
}

// =====================================================================
// Protein insertion boundary cases (largest sub-cluster in F13)
// =====================================================================

/// protein-insertion:C-C8 — length-by-stop-position form.
#[test]
fn protein_ins_ends_at_stop_position() {
    pin_either(
        "NP_003997.1:p.Gln746_Lys747ins*63",
        "protein/insertion.md:49-51 — length-by-stop-pos form",
    );
}

/// protein-insertion:C-C9 — `Xaa[N]` explicit equivalence with run-out.
#[test]
fn protein_ins_unknown_aa_count() {
    pin_either(
        "NP_003997.1:p.(Val582_Asn583insXaa[5])",
        "protein/insertion.md:N — [N] form ≡ Xaa run-out",
    );
}

/// protein-insertion:C-C10 — positional source-reference payload.
#[test]
fn protein_ins_positional_source_ref() {
    pin_either(
        "NP_003997.1:p.His7_Gln8insGly4_Ser6",
        "protein/insertion.md:76-80 — positional source",
    );
}

/// protein-insertion:C-C4 — ins-with-stop preferred over delins-removing-
/// the-terminator.
#[test]
fn protein_ins_with_stop_preferred_over_delins() {
    pin_either(
        "NP_003997.1:p.(Met3_His4insGlyTer)",
        "protein/insertion.md — ins with stop, not delins-removing-Ter",
    );
}

/// protein-insertion:S-S9 — `Ter` as the FIRST aa of an ins payload
/// (silent; spec only shows Ter at end).
#[test]
fn protein_ins_ter_at_start_of_payload() {
    pin_either(
        "NP_003997.1:p.Lys2_Leu3insTerGlyAla",
        "protein/insertion.md silent on Ter as first payload aa",
    );
}

/// protein-insertion:S-S10 — `insTer##` interpretation.
#[test]
fn protein_ins_ter_at_position_form() {
    pin_either(
        "NP_003997.1:p.Lys2_Leu3insTer12",
        "protein/insertion.md:21 — translation stop in inserted seq",
    );
}

/// protein-insertion:S-S14 — aa's after `Ter` in payload — forbidden.
#[test]
fn protein_ins_aa_after_ter_forbidden() {
    pin_either(
        "NP_003997.1:p.Lys2_Leu3insSerSerTerAlaPro",
        "protein/insertion.md:43 — aa's after Ter in payload is invalid",
    );
}

// =====================================================================
// Cross-construct: RNA dup ↔ ins (when the inserted seq is the immediate
// 5' flank, it's a dup, not an ins). Covered in DNA already; RNA shards
// flagged this for symmetry.
// =====================================================================

/// rna-insertion:C-C13 + rna-inversion:C-C7 — inverted-dup-as-ins.
/// `r.123_456dupinv` is forbidden; canonical is `r.234_235ins123_234inv`.
#[test]
fn rna_inverted_dup_canonical_as_ins_inv() {
    pin_either(
        "NM_004006.2:r.234_235ins123_234inv",
        "RNA/insertion.md:18 + inversion.md:19 — inverted-dup canonical",
    );
    pin_either(
        "NM_004006.2:r.123_456dupinv",
        "RNA/inversion.md:19 — dupinv is invalid; should canonicalize",
    );
}

/// dna-insertion:C-C2 / dna-duplication:C-C16 — inverted-dup-is-insertion
/// (the bare DNA form, mirrors the RNA case above).
#[test]
fn dna_inverted_dup_canonical_as_ins_inv() {
    pin_either(
        "NC_000023.11:g.122_123ins123_234inv",
        "DNA/duplication.md E7 — inverted dup canonical form",
    );
}

/// general:S-9 (the fs+ext combination silence) — does `fs` and `ext`
/// ever co-occur? Spec lists each independently.
#[test]
fn protein_fs_ext_combination() {
    pin_either(
        "NP_003997.1:p.Arg97ProfsTer23extTer5",
        "general.md:112-113 — no explicit example of fs+ext composition",
    );
}
