//! Miscellaneous spec-coverage corners across F6, F7, F10, F11, F14, F15,
//! F16, F17, F18 (77 net-new corners; ~25 representative probes selected).
//!
//! Each probe pins ferro's CURRENT behavior on a corner the construct shard
//! identified. Source: phase4_final_triage.md families F6 (mt origin
//! non-deletion), F7 (`=` no-change combinations), F10 (refseq/accession),
//! F11 (numbering corners), F14 (repeat canonicalization), F15 (adjoined /
//! `::` / complex), F16 (uncertain/predicted parens), F17 (IUPAC), F18
//! (grammar meta-glyphs).

use ferro_hgvs::parse_hgvs;

/// Pin that ferro currently parses `input`. A regression that stops
/// parsing an accepted form fails loudly. Some probes pin a current
/// over-acceptance (a spec-invalid form ferro still parses) — the note
/// flags those so they can be promoted to `pin_reject` once tightened.
fn pin_accept(input: &str, note: &str) {
    assert!(
        parse_hgvs(input).is_ok(),
        "expected parse success for {input:?} ({note})",
    );
}

/// Pin that ferro currently rejects `input` — invalid / spec-forbidden
/// forms, or valid forms ferro does not yet handle. When the behavior
/// flips, this fails loudly so the probe can be re-classified.
fn pin_reject(input: &str, note: &str) {
    assert!(
        parse_hgvs(input).is_err(),
        "expected parse failure for {input:?} ({note})",
    );
}

fn pin_round_trip(input: &str, expected: &str, note: &str) {
    let parsed = parse_hgvs(input)
        .unwrap_or_else(|e| panic!("expected parse of {input:?} ({note}); got: {e}"));
    assert_eq!(parsed.to_string(), expected, "for {input:?} ({note})");
}

// =====================================================================
// F6 — Origin/wraparound for non-deletion edits on circular contigs
// =====================================================================

/// dna-duplication:S-S6 — Duplication at `m.1` (mt origin).
#[test]
fn f6_mt_origin_duplication() {
    pin_accept(
        "NC_012920.1:m.1dup",
        "DNA/duplication.md silent on mt origin; SVD-WG006 deletion-only carve-out",
    );
}

/// dna-duplication:S-S7 — Duplication at `o.1` (chloroplast/plasmid).
#[test]
fn f6_o_origin_duplication() {
    pin_accept(
        "NC_000067.7:o.1dup",
        "DNA/duplication.md silent on o. origin",
    );
}

// =====================================================================
// F7 — `=` no-change combinations
// =====================================================================

/// dna-other:S-S13 — `=` confirmation that is itself uncertain.
#[test]
fn f7_uncertain_no_change_form() {
    pin_accept(
        "NM_004006.2:c.(123=)",
        "DNA/other.md silent; uncertain.md governs (…) wrap",
    );
}

/// dna-other:S-S14 — Phase-unknown two `=` confirmations.
#[test]
fn f7_phase_unknown_two_no_change() {
    pin_accept(
        "NM_004006.2:c.123=(;)456=",
        "DNA/other.md silent on phase-unknown two = confirmations",
    );
}

/// general:C-12 — Methylation `|met=` combines modification pipe with `=`.
#[test]
fn f7_methylation_with_no_change() {
    pin_accept(
        "NC_000023.11:g.123|met=",
        "general.md:125 + DNA/other.md:46-47",
    );
}

// =====================================================================
// F10 — Refseq / accession compliance
// =====================================================================

/// checklist:C-11 — `chr1:` is invalid as primary accession.
#[test]
fn f10_chr_prefix_invalid() {
    // Over-acceptance: `chr1:` is not a valid primary ref seq per
    // checklist.md, but ferro currently parses it. Promote to
    // `pin_reject` once accession validation is enforced.
    pin_accept(
        "chr1:g.1000_1005del",
        "checklist.md + style-examples.csv:7 — chr1 not a valid ref seq; currently parsed",
    );
}

/// refseq:S-S7 — `MT-XXX` as primary identifier.
#[test]
fn f10_mt_gene_as_primary_identifier() {
    // Over-acceptance: refseq.md:41 forbids a gene symbol as primary
    // identifier, but ferro currently parses `MT-TL1:`. Promote to
    // `pin_reject` once gene-as-primary is refused.
    pin_accept(
        "MT-TL1:m.3243A>G",
        "refseq.md:41 forbids gene-as-primary; mt parens-only exception; currently parsed",
    );
}

/// refseq:S-S6 — LRG with multi-version sub-LRG (LRG_199t3).
#[test]
fn f10_lrg_multi_version() {
    pin_accept(
        "LRG_199t3:c.76A>G",
        "refseq.md:159 mentions LRG_199t3 exists",
    );
}

// =====================================================================
// F11 — Numbering corner cases
// =====================================================================

/// checklist:C-2 — `c.123-65_-50` ambiguous endpoint.
#[test]
fn f11_ambiguous_intronic_endpoint() {
    pin_accept(
        "NM_004006.2:c.123-65_-50del",
        "checklist.md:26 — -50 endpoint ambiguous",
    );
}

/// dna-delins:C-C1 — `g.pter` position token.
#[test]
fn f11_pter_position_token() {
    pin_accept(
        "NC_000001.11:g.pter_1000del",
        "DNA/delins.md C1 — pter only in translocation example",
    );
}

/// protein-substitution:C-C5 — `p.0?` vs `p.0`.
#[test]
fn f11_protein_zero_with_predicted() {
    pin_accept(
        "NP_003997.1:p.0?",
        "protein/substitution.md:46-48 — predicted no-protein",
    );
    pin_accept(
        "NP_003997.1:p.0",
        "protein/substitution.md:46-48 — observed no-protein",
    );
}

/// dna-repeated:C-C14 — `c.0` is not a valid position.
#[test]
fn f11_c_zero_not_valid_position() {
    pin_reject(
        "NM_004006.2:c.0A>G",
        "numbering.md:31 — c.0 does not exist; should reject",
    );
}

// =====================================================================
// F14 — Repeat-tract canonicalization
// =====================================================================

/// rna-repeated:C-C3 — Redundant form invalid.
#[test]
fn f14_redundant_repeat_form_invalid() {
    pin_accept(
        "NM_004006.2:r.-125_-123cug[4]",
        "RNA/repeated.md:22 — redundant info form is invalid",
    );
}

/// dna-repeated:C-C4 — `c.-129CGG[79]` invalid for FMR1 locus.
#[test]
fn f14_fmr1_invalid_explicit_unit() {
    pin_accept(
        "NM_002024.5:c.-129CGG[79]",
        "DNA/repeated.md:63-65 — locus-specific invalidation",
    );
}

/// rna-duplication:S-S9 — Triplication `[3]` form on RNA.
#[test]
fn f14_rna_triplication_form() {
    // Valid per DNA/duplication.md:94-95 lifted to RNA, but ferro does
    // not yet parse this uncertain-breakpoint repeat form — pin the gap.
    pin_reject(
        "NM_004006.2:r.(low5_high5)_(low3_high3)[3]",
        "DNA/duplication.md:94-95 lifted to RNA by symmetry — currently a parse gap",
    );
}

// =====================================================================
// F15 — Adjoined transcript / fusion / `::` / complex structural
// =====================================================================

/// rna-adjoined_transcript canonical form.
#[test]
fn f15_adjoined_transcript_two_partner() {
    pin_accept(
        "NM_004006.2:r.1_100::NM_000088.3:r.200_300",
        "RNA/adjoined_transcript.md:17 — two-partner form",
    );
}

/// dna-complex:S-S4 — `::` use on DNA-level complex.
#[test]
fn f15_double_colon_on_dna() {
    pin_reject(
        "NC_000023.11:g.1000::NC_000023.11:g.5000",
        "DNA/complex silent; general.md:65 :: is fusion-only",
    );
}

/// Three-partner adjoined transcript — explicitly INVALID
/// (adjoined_transcript.md:17 two-partner only).
#[test]
fn f15_three_partner_adjoined_invalid() {
    pin_reject(
        "NM_004006.2:r.1_100::NM_000088.3:r.200_300::NM_111111.1:r.50_150",
        "adjoined_transcript.md:17 — exactly two partners",
    );
}

// =====================================================================
// F16 — Uncertain / predicted forms (parens placement, `^`, `?`)
// =====================================================================

/// general:C-8 — `(...)_(...)del` uncertain breakpoint form.
#[test]
fn f16_uncertain_breakpoint_form_round_trips() {
    pin_round_trip(
        "NC_000023.9:g.(123456_234567)_(345678_456789)del",
        "NC_000023.9:g.(123456_234567)_(345678_456789)del",
        "general.md:86 + uncertain.md — canonical uncertain breakpoint form",
    );
}

/// general:C-9 — `?` inside uncertain breakpoint for open-side.
#[test]
fn f16_open_side_uncertain_breakpoint() {
    pin_round_trip(
        "NC_000023.9:g.(?_234567)_(345678_?)del",
        "NC_000023.9:g.(?_234567)_(345678_?)del",
        "general.md + uncertain.md — open-on-one-side form",
    );
}

/// protein-alleles:S-S8 — Predicted-each-cis form `p.[(a);(b)]`.
#[test]
fn f16_protein_predicted_each_cis() {
    pin_reject(
        "NP_003997.1:p.[(Trp24Cys);(Lys25Arg)]",
        "protein/alleles.md silent on predicted-each-cis combination; ferro currently rejects",
    );
}

/// protein-deletion:S-S7 — Caret alternatives on a deletion.
#[test]
fn f16_protein_deletion_caret_alternatives() {
    pin_reject(
        "NP_003997.1:p.(Trp26^Trp27)del",
        "protein/deletion.md silent; uncertain.md ^ governs; ferro currently rejects",
    );
}

// =====================================================================
// F17 — IUPAC / ambiguity codes
// =====================================================================

/// general:S-8 — protein stop-codon glyph choice (Ter vs * vs X).
/// general.md:54 — Ter ≡ *; X has been removed.
#[test]
fn f17_protein_stop_glyph_x_removed() {
    // Over-acceptance: `X` as a stop glyph was removed (general.md:54),
    // but ferro currently parses `X` as the unknown-residue `Xaa` rather
    // than rejecting. Pins that current reinterpretation.
    pin_accept(
        "NP_003997.1:p.Trp24X",
        "general.md:54 — X removed; ferro currently parses X as Xaa",
    );
    pin_round_trip(
        "NP_003997.1:p.Trp24Ter",
        "NP_003997.1:p.Trp24Ter",
        "Ter canonical form",
    );
}

/// dna-repeated:C-C10 — Mixed-repeat with IUPAC-collapsed unit `GGM[108]`.
#[test]
fn f17_iupac_collapsed_repeat_unit() {
    pin_accept(
        "NM_004006.2:c.123GGM[108]",
        "DNA/repeated.md — IUPAC-collapsed mixed-repeat unit",
    );
}

/// rna-adjoined_transcript:C-C8 — `c.`/`n.` prefixes FORBIDDEN inside
/// adjoined transcripts.
#[test]
fn f17_adjoined_c_prefix_forbidden() {
    pin_reject(
        "NM_004006.2:c.1_100::NM_000088.3:c.200_300",
        "adjoined_transcript.md:18 — only r. allowed inside ::",
    );
}

// =====================================================================
// F18 — Grammar / EBNF meta-glyph collisions
// =====================================================================

/// general:C-1 — `p.Trp41*` parses but ferro canonicalizes to `Ter` on
/// Display. general.md:54 says both glyphs are acceptable; ferro picks
/// one canonical (`Ter`).
#[test]
fn f18_star_glyph_normalizes_to_ter_on_protein() {
    pin_round_trip(
        "NP_003997.1:p.Trp41*",
        "NP_003997.1:p.Trp41Ter",
        "ferro canonicalizes * → Ter on protein output (general.md:54 allows both)",
    );
    pin_round_trip(
        "NP_003997.1:p.Trp41Ter",
        "NP_003997.1:p.Trp41Ter",
        "Ter form is canonical",
    );
}

/// general:C-1 — `*N` on `c.` is a 3'-UTR-relative position, disambiguated
/// from `*` stop by prefix.
#[test]
fn f18_star_as_3utr_position_on_c() {
    // The point of this case is the `*41` 3'-UTR-relative position; the
    // reference base must be a DNA base (`U` is RNA-only and rejected on a
    // `c.` description — #486).
    pin_round_trip(
        "NM_004006.2:c.*41T>A",
        "NM_004006.2:c.*41T>A",
        "general.md:77 — * as 3'-UTR-relative position on c.",
    );
}
