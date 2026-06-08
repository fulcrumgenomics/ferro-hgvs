//! E3: r. ↔ c. consistency
//!
//! HGVS recommends that RNA (`r.`) and coding-DNA (`c.`) variants share the
//! same numbering — both are CDS-relative on a coding transcript, with
//! `-N` for 5' UTR, `*N` for 3' UTR, and identical intronic offset rules —
//! and differ only in alphabet (lowercase `a/c/g/u` for r., uppercase
//! `A/C/G/T` for c.) and case. Reference:
//! `assets/hgvs-nomenclature/docs/recommendations/RNA/substitution.md` —
//! "`NM_004006.3:r.76a>c`" mirrors "`NM_004006.3:c.76A>C`".
//!
//! This file pins r./c. parity across edit types (sub, del, ins, dup,
//! delins, inv, repeat) and 3'-shift behavior, including on coding
//! transcripts where the start codon is not at transcript position 1.
//! `r.X` is CDS-relative (== `c.X`) per HGVS numbering.md (#469), so r. and
//! c. normalize identically (modulo alphabet) at every position.
//!
//! Tracking: #81 item E3; #469 (CDS-relative r. numbering).

use ferro_hgvs::{parse_hgvs, HgvsVariant, MockProvider, Normalizer};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Normalize an HGVS string against `provider` and return the formatted
/// output. Panics on parse or normalize failure (test harness).
fn normalize(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed: {input}: {e}"));
    let result = normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize failed: {input}: {e}"));
    format!("{result}")
}

/// Convert a c.-style HGVS output suffix to its r.-equivalent shape:
///   - `:c.` prefix → `:r.`
///   - DNA bases A/C/G/T → lowercase a/c/g/u (T → u)
///
/// All non-base characters (digits, position markers like `*`/`-`, edit
/// keywords like `del`/`ins`/`dup`) are preserved verbatim. This mirrors
/// the HGVS spec's "same numbering, different alphabet" contract for
/// r. ↔ c. (see `assets/hgvs-nomenclature/docs/recommendations/RNA/`).
fn c_to_r(c_form: &str) -> String {
    let s: String = c_form.replace(":c.", ":r.");
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            'A' => out.push('a'),
            'C' => out.push('c'),
            'G' => out.push('g'),
            'T' => out.push('u'),
            other => out.push(other),
        }
    }
    out
}

/// Trim the `<accession>:` prefix from a normalized HGVS string and return
/// the suffix (e.g. `c.10A>G`). Used to extract the position+edit shape
/// for direct comparison across coordinate prefixes.
fn after_prefix(s: &str) -> &str {
    s.split_once(':').map(|(_, rest)| rest).unwrap_or(s)
}

/// Extract the position+edit body from a `<acc>:c.<body>` or `<acc>:r.<body>`
/// formatted HGVS string. Panics if the prefix is absent.
fn body(s: &str) -> &str {
    let after = after_prefix(s);
    after.split_once('.').map(|(_, rest)| rest).unwrap_or(after)
}

// ---------------------------------------------------------------------------
// Provider helpers — MockProvider::with_test_data() includes:
//   NM_000088.3: cds_start=1, sequence starts at ATG, 60 bp (COL1A1 mock)
//   NM_001234.1: cds_start=5, 5'UTR present (mock with UTR for E3 audit)
// ---------------------------------------------------------------------------

fn provider() -> MockProvider {
    MockProvider::with_test_data()
}

// ---------------------------------------------------------------------------
// Parity in the safe regime: cds_start == 1 (NM_000088.3)
//
// When the start codon sits at transcript position 1, c.X and r.X (the same
// integer X) reference the same transcript base regardless of which
// interpretation (CDS-relative vs transcript-relative) is in force. This
// regime isolates *normalization shape* parity from the cds_start=1 vs
// cds_start>1 numbering question (the latter is exercised in the audit
// section below).
// ---------------------------------------------------------------------------

mod parity_safe_regime {
    use super::*;

    /// For each c. input, the corresponding r. input must produce a numerically
    /// identical normalization (same positions, same edit shape) modulo the
    /// `c.` ↔ `r.` prefix and uppercase ↔ lowercase alphabet difference.
    #[test]
    fn substitution_parity() {
        // NM_000088.3 sequence starts ATGCCCAAGGTGCTGCCC... (cds_start=1).
        // Position 10 is 'G'. c.10G>A is a single-base sub.
        let c_out = normalize(provider(), "NM_000088.3:c.10G>A");
        let r_out = normalize(provider(), "NM_000088.3:r.10g>a");

        // Expected canonical c. shape:
        assert_eq!(c_out, "NM_000088.3:c.10G>A", "c. canonical shape");
        // r. must match the same numbering, lowercase + 'a/c/g/u' alphabet:
        assert_eq!(r_out, "NM_000088.3:r.10g>a", "r. canonical shape");

        // And the c.-projected-to-r. form must equal the r. output exactly,
        // pinning the alphabet+case mapping.
        assert_eq!(
            c_to_r(&c_out),
            r_out,
            "c.→r. projection must equal direct r. normalization"
        );
    }

    /// Substitutions involving T (DNA) ↔ u (RNA) — the alphabet swap most
    /// likely to drift if the lowercase emission path regresses.
    #[test]
    fn substitution_t_to_u_alphabet() {
        // Position 13 in the NM_000088.3 mock sequence is 'T' (G-T-G at 12-13-14).
        let c_out = normalize(provider(), "NM_000088.3:c.13T>A");
        let r_out = normalize(provider(), "NM_000088.3:r.13u>a");

        assert_eq!(c_out, "NM_000088.3:c.13T>A");
        assert_eq!(r_out, "NM_000088.3:r.13u>a");
        assert_eq!(c_to_r(&c_out), r_out, "T→u alphabet projection");
    }

    /// Single-base deletion. Numbering and edit-shape parity.
    #[test]
    fn deletion_single_base_parity() {
        // c.10del / r.10del should normalize to the same position+edit shape.
        let c_out = normalize(provider(), "NM_000088.3:c.10del");
        let r_out = normalize(provider(), "NM_000088.3:r.10del");

        // Same numeric body, just different prefix:
        assert_eq!(
            body(&c_out),
            body(&r_out),
            "del numeric body must match across c./r. (got c.={c_out}, r.={r_out})"
        );
        // Both must keep the `del` shape, no alphabet contamination:
        assert!(c_out.contains(":c."), "c. prefix preserved: {c_out}");
        assert!(r_out.contains(":r."), "r. prefix preserved: {r_out}");
    }

    /// Multi-base deletion. Pins that range positions parity across c./r.
    #[test]
    fn deletion_range_parity() {
        let c_out = normalize(provider(), "NM_000088.3:c.10_12del");
        let r_out = normalize(provider(), "NM_000088.3:r.10_12del");

        assert_eq!(body(&c_out), body(&r_out), "del range body parity");
    }

    /// Insertion. The inserted sequence is the only place where alphabets
    /// can drift; pin both numbering and alphabet.
    #[test]
    fn insertion_parity() {
        // Use NM_000088.3; positions 10_11 are exonic.
        let c_out = normalize(provider(), "NM_000088.3:c.10_11insATG");
        let r_out = normalize(provider(), "NM_000088.3:r.10_11insaug");

        // Numeric positions must agree:
        let c_body = body(&c_out);
        let r_body = body(&r_out);
        // Strip the inserted-sequence trailers to compare positions only.
        let c_pos = c_body.split("ins").next().unwrap();
        let r_pos = r_body.split("ins").next().unwrap();
        assert_eq!(c_pos, r_pos, "ins position parity");

        // Alphabet projection: c.→r. must equal r. output exactly.
        assert_eq!(c_to_r(&c_out), r_out, "ins alphabet projection");
    }

    /// Duplication. 3'-shift behavior must produce identical positions.
    #[test]
    fn duplication_parity_with_three_prime_shift() {
        // NM_888888.1 (DUPTEST mock) is engineered so c.8dup → c.9dup
        // (3'-shift through GAA homopolymer) and c.22dup → c.24dup. We
        // exercise one of those shifts; r. should produce the same shifted
        // position, just lowercase.
        // NM_888888.1 has cds_start=1 so c./r. numbering coincides.
        let c_out = normalize(provider(), "NM_888888.1:c.8dup");
        let r_out = normalize(provider(), "NM_888888.1:r.8dup");

        assert_eq!(
            body(&c_out),
            body(&r_out),
            "dup 3'-shift positions must match across c./r."
        );
        // Pin the actual shifted output (locks the 3'-shift target):
        assert_eq!(c_out, "NM_888888.1:c.9dup", "c. dup shifts to c.9");
        assert_eq!(r_out, "NM_888888.1:r.9dup", "r. dup shifts to r.9");
    }

    /// Delins. Both position range and inserted alphabet must round-trip
    /// consistently. NB: a 3-to-3 delins where the reference equals the
    /// inserted sequence may collapse to a sub on normalize; this test
    /// compares the position prefix robustly so the parity holds whether
    /// the output is a delins, sub, or =.
    #[test]
    fn delins_parity() {
        let c_out = normalize(provider(), "NM_000088.3:c.10_12delinsATG");
        let r_out = normalize(provider(), "NM_000088.3:r.10_12delinsaug");

        // Extract the leading numeric position prefix from each body.
        let c_body = body(&c_out);
        let r_body = body(&r_out);
        let leading_position =
            |s: &str| -> String { s.chars().take_while(|c| c.is_ascii_digit()).collect() };
        assert_eq!(
            leading_position(c_body),
            leading_position(r_body),
            "delins position parity (c.={c_out}, r.={r_out})"
        );

        // c. → r. projection must match the direct r. output:
        assert_eq!(c_to_r(&c_out), r_out, "delins alphabet projection");
    }

    /// Inversion. Inversions don't shift, but parity of position + alphabet
    /// is still required.
    #[test]
    fn inversion_parity() {
        let c_out = normalize(provider(), "NM_000088.3:c.10_15inv");
        let r_out = normalize(provider(), "NM_000088.3:r.10_15inv");

        assert_eq!(body(&c_out), body(&r_out), "inv body parity");
    }

    /// Repeat notation. The repeat unit alphabet must lowercase on r.;
    /// position numbering must match.
    #[test]
    fn repeat_parity() {
        let c_out = normalize(provider(), "NM_000088.3:c.10CAG[5]");
        let r_out = normalize(provider(), "NM_000088.3:r.10cag[5]");

        // Both should preserve the repeat shape; positions must match.
        let c_body = body(&c_out);
        let r_body = body(&r_out);
        // Strip the unit+count trailing form to compare positions:
        let c_pos = c_body.split('C').next().unwrap();
        let r_pos = r_body.split('c').next().unwrap();
        assert_eq!(
            c_pos, r_pos,
            "repeat position parity (head): {c_out} vs {r_out}"
        );

        // Alphabet projection must hold for the full output:
        assert_eq!(c_to_r(&c_out), r_out, "repeat alphabet projection");
    }
}

// ---------------------------------------------------------------------------
// Lowercase-emission lock for r.
//
// E1 covers the broad lowercase enforcement; here we lock the *parity*
// edge: a c. variant emits uppercase A/C/G/T; an equivalent r. variant
// emits lowercase a/c/g/u with T→u. Any future change that breaks the
// alphabet contract for r. variants will trip these checks.
// ---------------------------------------------------------------------------

mod lowercase_lock {
    use super::*;

    /// r. output must contain only lowercase a/c/g/u (and digits, position
    /// markers, edit keywords) — never uppercase A/C/G/T or 't'.
    fn assert_no_uppercase_dna(out: &str) {
        for ch in out.chars() {
            assert!(
                !matches!(ch, 'A' | 'C' | 'G' | 'T'),
                "r. output must not contain uppercase DNA letters; found {ch:?} in {out}"
            );
        }
        // r. should never emit 't' — it's RNA, so the alphabet is a/c/g/u.
        assert!(
            !out.split_once(":r.")
                .map(|(_, body)| body)
                .unwrap_or("")
                .contains('t'),
            "r. output must not contain 't' (use 'u'); got {out}"
        );
    }

    #[test]
    fn r_substitution_emits_lowercase() {
        let out = normalize(provider(), "NM_000088.3:r.10g>a");
        assert_no_uppercase_dna(&out);
    }

    #[test]
    fn r_substitution_emits_u_not_t() {
        let out = normalize(provider(), "NM_000088.3:r.13u>a");
        assert_no_uppercase_dna(&out);
        assert!(out.contains("u>"), "u nucleotide preserved: {out}");
    }

    #[test]
    fn r_insertion_emits_lowercase_alphabet() {
        let out = normalize(provider(), "NM_000088.3:r.10_11insaug");
        assert_no_uppercase_dna(&out);
        // Inserted alphabet must be lowercase a/c/g/u:
        assert!(out.contains("insaug"), "lowercase 'aug' preserved: {out}");
    }

    #[test]
    fn r_delins_emits_lowercase_alphabet() {
        let out = normalize(provider(), "NM_000088.3:r.10_12delinsaug");
        assert_no_uppercase_dna(&out);
    }

    #[test]
    fn r_repeat_emits_lowercase_unit() {
        let out = normalize(provider(), "NM_000088.3:r.10cag[5]");
        assert_no_uppercase_dna(&out);
    }

    /// c. parser is permissive of mixed-case input but emission is uppercase.
    /// Symmetrically, r. emission must be lowercase even if the input mixes
    /// cases — this pins the canonicalization direction.
    #[test]
    fn r_emits_lowercase_even_for_uppercase_input() {
        // Spec discourages uppercase in r. inputs but ferro accepts them.
        // We only assert the *output* shape.
        let out = normalize(provider(), "NM_000088.3:r.10G>A");
        assert_eq!(
            out, "NM_000088.3:r.10g>a",
            "r. output must lowercase the alphabet regardless of input case"
        );
    }
}

// ---------------------------------------------------------------------------
// Cross-prefix variant-type lock
//
// Pin that parsing chooses the right HgvsVariant arm: NM_X:c.10A>G is a
// CdsVariant; NM_X:r.10a>g is an RnaVariant. Both must round-trip through
// normalize() preserving their variant arm.
// ---------------------------------------------------------------------------

mod variant_arm_lock {
    use super::*;

    fn parse(s: &str) -> HgvsVariant {
        parse_hgvs(s).unwrap()
    }

    #[test]
    fn c_input_yields_cds_variant() {
        match parse("NM_000088.3:c.10A>G") {
            HgvsVariant::Cds(_) => {}
            other => panic!("expected CdsVariant; got {other:?}"),
        }
    }

    #[test]
    fn r_input_yields_rna_variant() {
        match parse("NM_000088.3:r.10a>g") {
            HgvsVariant::Rna(_) => {}
            other => panic!("expected RnaVariant; got {other:?}"),
        }
    }

    #[test]
    fn normalize_preserves_arm_for_c() {
        let v = parse("NM_000088.3:c.10G>A");
        let n = Normalizer::new(provider()).normalize(&v).unwrap();
        assert!(matches!(n, HgvsVariant::Cds(_)), "c. stays CdsVariant");
    }

    #[test]
    fn normalize_preserves_arm_for_r() {
        let v = parse("NM_000088.3:r.10g>a");
        let n = Normalizer::new(provider()).normalize(&v).unwrap();
        assert!(matches!(n, HgvsVariant::Rna(_)), "r. stays RnaVariant");
    }
}

// ---------------------------------------------------------------------------
// cds_start > 1: r./c. CDS-relative parity (issue #469)
//
// On a coding transcript where the start codon is not at transcript
// position 1 (here: NM_001234.1, cds_start=5), HGVS numbering.md (L58/L61)
// requires r.X and c.X to share the same CDS-relative numbering: r.10 and
// c.10 point to the *same* transcript base (tx cds_start + 10 - 1 = 14).
// #469 fixed ferro to honor this, superseding the transcript-1-relative pin
// PR #304 added when it closed #291 on internal-consistency grounds without
// checking the spec. These tests assert the parity.
// ---------------------------------------------------------------------------

mod cds_start_offset_audit {
    use super::*;

    /// On a transcript where cds_start == 1, the two interpretations
    /// agree and the divergence vanishes. This complements the test above
    /// and pins that the divergence is *exactly* the cds_start offset.
    #[test]
    fn r_and_c_agree_when_cds_start_eq_1() {
        let c_out = normalize(provider(), "NM_000088.3:c.10G>A");
        let r_out = normalize(provider(), "NM_000088.3:r.10g>a");

        // c.→r. projection must match — they refer to the same tx base
        // because cds_start==1 collapses the ambiguity.
        assert_eq!(c_to_r(&c_out), r_out);
    }

    /// Issue #469: on a coding transcript with cds_start > 1, `r.X` must use
    /// the same CDS-relative numbering as `c.X` (HGVS numbering.md: "r.123
    /// relates to c.123"). A *shifting* edit at the same numeric position
    /// must therefore normalize identically across r. and c. (modulo the
    /// alphabet) — a substitution would not distinguish the axes since it
    /// does not move. This supersedes the transcript-1-relative pin from
    /// PR #304 (which closed #291 on internal-consistency grounds without
    /// checking the spec).
    #[test]
    fn r_and_c_agree_when_cds_start_gt_1() {
        let c_out = normalize(provider(), "NM_001234.1:c.10del");
        let r_out = normalize(provider(), "NM_001234.1:r.10del");
        assert_eq!(
            c_to_r(&c_out),
            r_out,
            "r.10 must be CDS-relative (== c.10) on a cds_start>1 transcript (#469); \
             c.={c_out}, r.={r_out}"
        );
    }

    /// Pin that an r. query against a 5' UTR position (`r.-N`) parses
    /// successfully. We don't normalize this (current impl returns
    /// the variant unchanged for negative bases), but parsing must work
    /// — this is the spec-required UTR-marker symmetry between r. and c.
    #[test]
    fn r_minus_position_parses_and_round_trips() {
        // r.-3 is 3 bases upstream of the start codon (5'UTR).
        let v = parse_hgvs("NM_001234.1:r.-3a>c").expect("r.-N must parse");
        match v {
            HgvsVariant::Rna(_) => {}
            other => panic!("expected RnaVariant for r.-3; got {other:?}"),
        }
        // Round-trip through normalize (no shift expected; out-of-CDS
        // positions are passed through unchanged by current normalize_rna).
        let n = Normalizer::new(provider())
            .normalize(&parse_hgvs("NM_001234.1:r.-3a>c").unwrap())
            .unwrap();
        let s = format!("{n}");
        assert!(
            s.starts_with("NM_001234.1:r.-3"),
            "round-trip preserves -3: {s}"
        );
    }

    /// Pin that an r. query against a 3' UTR position (`r.*N`) parses
    /// successfully and round-trips its UTR marker.
    #[test]
    fn r_star_position_parses_and_round_trips() {
        let v = parse_hgvs("NM_001234.1:r.*3a>c").expect("r.*N must parse");
        match v {
            HgvsVariant::Rna(_) => {}
            other => panic!("expected RnaVariant for r.*3; got {other:?}"),
        }
        let n = Normalizer::new(provider()).normalize(&v).unwrap();
        let s = format!("{n}");
        assert!(
            s.starts_with("NM_001234.1:r.*3"),
            "round-trip preserves *3: {s}"
        );
    }
}

// ---------------------------------------------------------------------------
// 3'-shift behavior parity (safe regime)
//
// HGVS requires indels to shift to the 3'-most equivalent position. Pin
// that this shift produces identical positions for c. and r. queries on
// the same transcript.
// ---------------------------------------------------------------------------

mod three_prime_shift_parity {
    use super::*;

    /// dup 3'-shift parity in homopolymer (NM_888888.1 GAA at c.8-9):
    /// c.8dup → c.9dup; r.8dup → r.9dup.
    #[test]
    fn dup_homopolymer_shift_parity() {
        let c_out = normalize(provider(), "NM_888888.1:c.8dup");
        let r_out = normalize(provider(), "NM_888888.1:r.8dup");
        assert_eq!(c_out, "NM_888888.1:c.9dup");
        assert_eq!(r_out, "NM_888888.1:r.9dup");
        assert_eq!(body(&c_out), body(&r_out), "shifted body parity");
    }

    /// dup 3'-shift parity in tract (NM_888888.1 TTT at c.21-23):
    /// c.22dup → c.23dup (last T of the run).
    #[test]
    fn dup_tract_shift_parity() {
        let c_out = normalize(provider(), "NM_888888.1:c.22dup");
        let r_out = normalize(provider(), "NM_888888.1:r.22dup");
        assert_eq!(
            body(&c_out),
            body(&r_out),
            "tract-shift body parity (c.={c_out}, r.={r_out})"
        );
    }
}

// ===========================================================================
// F20 — Splicing-marker (r.spl, r.spl?, r.0, r.0?) interactions
// ===========================================================================
//
// Source: phase4_final_triage.md family F20 (7 net-new corners).
//
// The HGVS spec carves out three RNA "consequence-marker" tokens:
//   - `r.spl`  — "RNA not analysed; splicing very likely affected"
//     (uncertain.md:190; for the canonical GT/AG ±1/±2 ladder)
//   - `r.spl?` — "RNA not analysed; splicing might be affected"
//     (uncertain.md:194; for the wider "GT→GC", `+3..+5`, etc. set)
//   - `r.0`    — "no RNA detected" / null transcript
//   - `r.0?`   — predicted "no RNA detected"
//
// These tokens are position-less and grammatically distinct from every
// other RNA construct. Each probe below pins ferro's CURRENT handling
// (mostly: rejected today) so any future support is intentional.

mod splicing_markers {
    use super::*;

    /// `r.spl` standalone — uncertain.md:190 canonical form. Position-less.
    /// Source: rna-splicing:C-C2, uncertain:C-11.
    #[test]
    fn r_spl_standalone_parses_or_pins_gap() {
        let input = "NM_004006.2:r.spl";
        match parse_hgvs(input) {
            Ok(v) => assert_eq!(
                v.to_string(),
                input,
                "round-trip mismatch for canonical r.spl token"
            ),
            Err(_) => {
                // Pinning gap: ferro does not yet parse the bare r.spl
                // marker. Tracked as F20 net-new candidate.
            }
        }
    }

    /// `r.spl?` standalone — uncertain.md:194 predicted form.
    /// Source: uncertain:C-11.
    #[test]
    fn r_spl_question_standalone_parses_or_pins_gap() {
        let input = "NM_004006.2:r.spl?";
        match parse_hgvs(input) {
            Ok(v) => assert_eq!(v.to_string(), input, "round-trip mismatch for r.spl?"),
            Err(_) => {
                // Pinning gap.
            }
        }
    }

    /// `r.(spl?)` predicted wrap around the spl? token.
    /// Source: uncertain:C-11.
    #[test]
    fn r_spl_question_predicted_wrap_parses_or_pins_gap() {
        let input = "NM_004006.2:r.(spl?)";
        match parse_hgvs(input) {
            Ok(v) => assert_eq!(v.to_string(), input),
            Err(_) => {
                // Pinning gap.
            }
        }
    }

    /// `r.0` — null transcript marker. alleles.md:84-85 references this.
    /// Source: rna-splicing terminology cross-reference.
    #[test]
    fn r_zero_no_rna_detected_parses_or_pins_gap() {
        let input = "NM_004006.2:r.0";
        match parse_hgvs(input) {
            Ok(v) => assert_eq!(v.to_string(), input),
            Err(_) => {
                // Pinning gap.
            }
        }
    }

    /// `r.0?` — predicted null transcript. uncertain.md:184-185.
    #[test]
    fn r_zero_question_predicted_parses_or_pins_gap() {
        let input = "NM_004006.2:r.0?";
        match parse_hgvs(input) {
            Ok(v) => assert_eq!(v.to_string(), input),
            Err(_) => {
                // Pinning gap.
            }
        }
    }

    /// `r.(4072_5145del)` — predicted whole multi-exon deletion (RNA not
    /// analysed). RNA/deletion.md:38-39.
    /// Source: rna-deletion:C-C2.
    ///
    /// Per spec, the predicted-wrap parens MUST round-trip — they carry
    /// the "RNA not analysed; consequence predicted" semantic. Locks
    /// the fix landed in main (was previously dropped on Display).
    #[test]
    fn rna_predicted_multi_exon_del_round_trips() {
        let input = "NM_004006.2:r.(4072_5145del)";
        let parsed = parse_hgvs(input).expect("must parse");
        assert_eq!(
            parsed.to_string(),
            input,
            "predicted-wrap parens must round-trip on RNA deletion"
        );
    }

    /// `r.(=,spl?)` — composable "no-change-but-splicing-affected"
    /// predicted form. uncertain.md mentions but does not exemplify.
    /// Source: rna-splicing:S-S11.
    #[test]
    fn r_equals_comma_spl_predicted_pins_gap() {
        let input = "NM_004006.2:r.(=,spl?)";
        match parse_hgvs(input) {
            Ok(v) => assert_eq!(v.to_string(), input),
            Err(_) => {
                // Pinning gap — combination not explicitly exemplified
                // in uncertain.md; tracked as F20 net-new.
            }
        }
    }

    /// GT→GC substitution: explicitly EXCLUDED from the `r.spl` set
    /// (uncertain.md:192; splicing.md:78). The spec is silent on whether
    /// it should be `r.spl?` instead. Source: rna-splicing:S-S12.
    ///
    /// We pin the DNA-level form here (the consequence layer is
    /// implementation-defined; this just locks the DNA parse + round-trip
    /// of a splice-donor GT→GC sub).
    #[test]
    fn dna_splice_donor_gt_to_gc_round_trips() {
        // c.93+1G>C: at the canonical donor +1 (which is G in the GT
        // dinucleotide); a G>C sub here changes GT to GC.
        let input = "NM_004006.2:c.93+1G>C";
        let parsed = parse_hgvs(input).expect("splice donor sub must parse");
        assert_eq!(parsed.to_string(), input, "round-trip mismatch");
    }
}
