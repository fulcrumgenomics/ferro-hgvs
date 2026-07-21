//! Allele-bracket grammar corner cases (F12 — 20 net-new corners).
//!
//! F12 was Stage 1's noisiest family per the Stage 2 review: Stage 1
//! conflated "3 variants in 1 arm" (POS, exemplified at
//! DNA/alleles.md:56) with "3 arms total" (POL, silent). This file pins
//! each corner explicitly so the distinction is preserved.
//!
//! Source: phase4_final_triage.md family F12.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Pin that ferro currently parses `input`. Use for forms that are valid
/// (or that ferro intentionally accepts) — a regression that stops
/// parsing them fails loudly.
fn pin_accept(input: &str, note: &str) {
    assert!(
        parse_hgvs(input).is_ok(),
        "expected parse success for {input:?} ({note})",
    );
}

/// Pin that ferro currently rejects `input`. Use for invalid / retracted
/// forms (and for valid forms ferro does not yet handle). When ferro's
/// behavior flips — an invalid form starts parsing, or a parse-gap on a
/// valid form closes — this fails loudly so the probe can be re-classified.
fn pin_reject(input: &str, note: &str) {
    assert!(
        parse_hgvs(input).is_err(),
        "expected parse failure for {input:?} ({note})",
    );
}

/// Pin only that ferro doesn't crash on `input`. Reserved for genuinely
/// spec-silent forms where neither accept nor reject is the
/// spec-mandated outcome.
fn pin_either(input: &str, _note: &str) {
    let _ = parse_hgvs(input);
}

fn pin_round_trip(input: &str, expected: &str, note: &str) {
    let parsed = parse_hgvs(input)
        .unwrap_or_else(|e| panic!("expected parse of {input:?} ({note}); got error: {e}"));
    assert_eq!(parsed.to_string(), expected, "for {input:?} ({note})");
}

/// Pin the canonical form ferro's normalizer produces for `input`.
/// Unlike `pin_round_trip` (which only checks parse + Display), this
/// runs the full normalize pass so adjacency-driven canonicalization
/// (e.g. an adjacent-variant allele collapsing to a `delins`) is
/// actually exercised. Uses an empty `MockProvider` — these
/// adjacency rewrites are sequence-independent.
fn pin_canonicalizes_to(input: &str, expected: &str, note: &str) {
    let parsed = parse_hgvs(input)
        .unwrap_or_else(|e| panic!("expected parse of {input:?} ({note}); got error: {e}"));
    let normalized = Normalizer::new(MockProvider::new())
        .normalize(&parsed)
        .unwrap_or_else(|e| panic!("expected normalize of {input:?} ({note}); got error: {e}"));
    assert_eq!(normalized.to_string(), expected, "for {input:?} ({note})");
}

// =====================================================================
// Three variants in one arm — VALID (spec example) vs three arms — SILENT
// =====================================================================

/// DNA/alleles.md:56 explicitly shows three variants in a cis allele.
/// dna-alleles:S-S9 (and related).
#[test]
fn allele_three_variants_one_arm_dna_round_trips() {
    let input = "NM_004006.2:c.[76A>C;87G>T;112T>G]";
    pin_round_trip(
        input,
        input,
        "DNA/alleles.md:56 three-variants-in-one-arm spec example",
    );
}

/// Three-arms-total form is spec-silent (would represent polyploidy).
#[test]
fn allele_three_arms_total_polyploidy_currently_handled() {
    pin_either(
        "NM_004006.2:c.[76A>C];[87G>T];[112T>G]",
        "polyploidy: 3 arms (spec-silent; only 2-arm trans shown)",
    );
}

// =====================================================================
// Predicted variant inside allele — protein-style placement
// =====================================================================

/// dna-alleles:S-S9 / protein-alleles:C-C25 — predicted variant inside
/// allele bracket. protein/alleles.md:91 shows `p.[Phe233Leu;(Cys690Trp)]`.
#[test]
fn allele_with_predicted_variant_inside() {
    pin_accept(
        "NM_004006.2:c.[76A>C;(87G>T)]",
        "predicted-inside-allele form analogous to protein/alleles.md:91",
    );
}

#[test]
fn protein_allele_with_predicted_variant_inside() {
    // protein/alleles.md:91 spec example: a predicted member inside a
    // protein allele bracket now parses and round-trips (#544).
    pin_round_trip(
        "NP_003997.1:p.[Phe233Leu;(Cys690Trp)]",
        "NP_003997.1:p.[Phe233Leu;(Cys690Trp)]",
        "protein/alleles.md:91 spec example — predicted member inside allele bracket",
    );
}

// =====================================================================
// Bracket-shortcut FORBIDDEN forms (DNA + RNA + protein)
// =====================================================================

/// dna-alleles:C-C15 / rna-alleles:C-C2 / protein-alleles:C-C2 — the
/// factored-prefix and empty-arm forms are explicitly invalid per
/// DNA/alleles.md:54 and inherited by RNA + protein.
#[test]
fn allele_factored_prefix_form_invalid() {
    pin_reject(
        "NM_004006.2:c.2376[G>C];[G>C]",
        "DNA/alleles.md:54 — factored-prefix form is invalid",
    );
}

#[test]
fn allele_empty_arm_form_invalid() {
    pin_reject(
        "NM_004006.2:c.2376G>C[];[]",
        "DNA/alleles.md:54 — empty-arm form is invalid",
    );
}

#[test]
fn protein_allele_factored_prefix_form_invalid() {
    pin_reject(
        "NP_003997.1:p.Ser68[Arg];[Arg]",
        "protein/alleles.md:C2 — inherited from DNA/alleles.md:54",
    );
}

#[test]
fn rna_allele_factored_prefix_form_invalid() {
    pin_reject(
        "NM_004006.2:r.76[a>u];[a>u]",
        "RNA/alleles.md:C2 — inherited from DNA/alleles.md:54",
    );
}

// =====================================================================
// Adjacent two-variant cis MUST be delins
// =====================================================================

/// dna-alleles:C-C2 — adjacent positions (≤1 nt apart) re-encoded as
/// delins per general.md:34-39.
#[test]
fn allele_adjacent_two_variants_must_be_delins() {
    // The adjacent-position allele canonicalizes to the delins under
    // normalize (general.md:34-39) — verify the rewrite, not just that
    // both forms parse independently.
    pin_canonicalizes_to(
        "NM_004006.2:c.[2076G>T;2077C>G]",
        "NM_004006.2:c.2076_2077delinsTG",
        "general.md:34-39 — adjacent positions canonicalize to delins",
    );
    // The canonical delins form round-trips unchanged.
    pin_round_trip(
        "NM_004006.2:c.2076_2077delinsTG",
        "NM_004006.2:c.2076_2077delinsTG",
        "canonical delins of adjacent positions",
    );
}

/// protein-delins:C-C5 — two consecutive aa changes MUST be delins per
/// spec, but ferro does NOT yet canonicalize the protein allele form to
/// delins (the DNA-axis adjacency rewrite has no protein counterpart).
/// Pin the current divergence: the allele normalizes to itself. Promote
/// to a `pin_canonicalizes_to(..., delins, ...)` once protein adjacency
/// canonicalization is implemented.
#[test]
fn protein_consecutive_aa_changes_must_be_delins() {
    pin_canonicalizes_to(
        "NP_003997.1:p.[Arg76Ser;Cys77Trp]",
        "NP_003997.1:p.Arg76_Cys77delinsSerTrp",
        "protein/substitution.md:23 + delins.md:62-64 — two consecutive-residue cis \
         substitutions canonicalize to a single delins",
    );
    // The spec-canonical delins form itself round-trips unchanged.
    pin_round_trip(
        "NP_003997.1:p.Arg76_Cys77delinsSerTrp",
        "NP_003997.1:p.Arg76_Cys77delinsSerTrp",
        "canonical delins form for consecutive aa changes",
    );
}

/// protein-delins:C-C10 — two changes SEPARATED by ≥1 unchanged
/// residue stay as allele, NOT delins (delins.md:62-64).
#[test]
fn protein_separated_changes_stay_as_allele() {
    pin_accept(
        "NP_003997.1:p.[Ser44Arg;Trp46Arg]",
        "protein/delins.md:62-64 — stays as allele (separated by unchanged res)",
    );
}

/// A predicted (`( )`) cis pair of consecutive-residue substitutions
/// canonicalizes to a predicted delins — the wrapper is carried onto the
/// combined form (protein/substitution.md:23).
#[test]
fn protein_predicted_consecutive_changes_canonicalize_to_predicted_delins() {
    pin_canonicalizes_to(
        "NP_003997.1:p.[(Arg76Ser);(Cys77Trp)]",
        "NP_003997.1:p.(Arg76_Cys77delinsSerTrp)",
        "protein/substitution.md:23 — predicted consecutive subs → predicted delins",
    );
}

/// A run of three or more consecutive-residue substitutions merges into one
/// delins spanning the whole run.
#[test]
fn protein_three_consecutive_changes_merge_into_one_delins() {
    pin_canonicalizes_to(
        "NP_003997.1:p.[Arg76Ser;Cys77Trp;Gly78Asp]",
        "NP_003997.1:p.Arg76_Gly78delinsSerTrpAsp",
        "protein/substitution.md:23 — a maximal adjacent run collapses to one delins",
    );
}

/// A partial allele — one adjacent run plus a separated member — merges only
/// the run and keeps the distant member as its own allele arm.
#[test]
fn protein_partial_run_merges_only_the_adjacent_pair() {
    pin_canonicalizes_to(
        "NP_003997.1:p.[Arg76Ser;Cys77Trp;Gly90Asp]",
        "NP_003997.1:p.[Arg76_Cys77delinsSerTrp;Gly90Asp]",
        "protein/substitution.md:23 — only the strictly-adjacent 76_77 run coalesces",
    );
}

/// Separated changes survive **normalization** unchanged — the coalescing
/// pass must not reach across the unchanged residue at 45 (delins.md:63).
#[test]
fn protein_separated_changes_survive_normalization() {
    pin_canonicalizes_to(
        "NP_003997.1:p.[Ser44Arg;Trp46Arg]",
        "NP_003997.1:p.[Ser44Arg;Trp46Arg]",
        "protein/delins.md:63 — a one-residue gap keeps the members separate",
    );
}

/// Two substitutions at the **same** residue are contradictory (one residue
/// cannot be two things in cis), not a consecutive-residue delins — so they are
/// not coalesced. #1098's cis-member sort still reorders them into a
/// deterministic (descriptor tie-break) order so the canonical string is
/// input-order-independent; a same-residue two-member bracket stays spec-valid.
#[test]
fn protein_same_residue_changes_are_not_coalesced() {
    pin_canonicalizes_to(
        "NP_003997.1:p.[Asp2His;Asp2Glu]",
        "NP_003997.1:p.[Asp2Glu;Asp2His]",
        "protein/substitution.md:23 — same-residue pair is not a delins; #1098 orders members",
    );
}

/// A `trans` allele (`];[`) describes two different physical alleles, so its
/// members are never coalesced into one cis delins.
#[test]
fn protein_trans_consecutive_changes_stay_separate() {
    pin_canonicalizes_to(
        "NP_003997.1:p.[Arg76Ser];[Cys77Trp]",
        "NP_003997.1:p.[Arg76Ser];[Cys77Trp]",
        "protein/substitution.md:23 — trans members are independent, never a cis delins",
    );
}

/// A run containing a nonsense (to-`Ter`) substitution must NOT coalesce: the
/// spec forbids listing residues after the stop, so `delins…Ter…` would be
/// spec-invalid AND unparseable (a round-trip break). The allele is left as
/// authored. Asserting the output equals the (parseable) input also pins that
/// normalize does not emit a string it cannot re-read.
#[test]
fn protein_to_ter_run_is_not_coalesced() {
    pin_canonicalizes_to(
        "NP_003997.1:p.[Arg76Ter;Cys77Trp]",
        "NP_003997.1:p.[Arg76Ter;Cys77Trp]",
        "protein/substitution.md:20, protein/delins.md:45 — no residues after Ter; do not merge",
    );
}

/// Members authored in descending residue order are reordered to ascending
/// residue order by #1098's cis-member sort.
///
/// KNOWN LIMITATION (protein axis): the descending input reorders to the bracket
/// `p.[Arg76Ser;Cys77Trp]` but is *not* re-merged into the delins
/// `protein/substitution.md:23` requires (`p.Arg76_Cys77delinsSerTrp`). #1103/#1106
/// makes the *nucleotide* merge input-order-independent by sorting before the
/// `merge_consecutive_edits` pipeline, but the protein delins-canonicalization
/// (`coalesce_protein_adjacent_substitutions`) is a separate pass that still runs
/// before the sort, so this protein order-dependence persists after #1106 — a
/// candidate follow-up (re-run the protein coalesce after the sort). The sort is
/// nonetheless the machinery that makes the ascending delins reachable at all.
#[test]
fn protein_descending_order_members_reorder_to_residue_order() {
    pin_canonicalizes_to(
        "NP_003997.1:p.[Cys77Trp;Arg76Ser]",
        "NP_003997.1:p.[Arg76Ser;Cys77Trp]",
        "#1098 orders cis members; protein delins-canonicalization does not re-run after the sort",
    );
}

// =====================================================================
// Single-bracket allele description (misleading per spec)
// =====================================================================

/// dna-alleles + protein-alleles:S-S19 / C-C14 — single-bracket forms
/// like `p.[Ser73Arg]` are "misleading" per alleles.md:100 but not
/// explicitly invalid.
#[test]
fn protein_single_bracket_allele_misleading() {
    pin_accept(
        "NP_003997.1:p.[Ser73Arg]",
        "alleles.md:100 — misleading but not explicitly invalid",
    );
}

// =====================================================================
// Allele × inversion-with-substitution (complex compositions)
// =====================================================================

/// dna-alleles:C-C17 / dna-complex:C2/C5 — allele combining inv + sub
/// at breakpoint (complex.md:113).
#[test]
fn allele_inv_with_substitution_at_breakpoint() {
    pin_accept(
        "NC_000023.11:g.[776788_93191545inv;93191546T>C]",
        "complex.md:113 — sub at inv breakpoint",
    );
}

/// dna-duplication:C-C14 — allele combining dup + ins.
#[test]
fn allele_dup_with_ins() {
    pin_accept(
        "NM_004006.2:c.[20_21dup;50_51insATG]",
        "DNA/duplication.md C14 — dup+ins composition",
    );
}

// =====================================================================
// Allele with two deletions on one arm
// =====================================================================

/// protein-deletion:S-S11 — allele with two protein deletions on one
/// arm (deletion.md silent; alleles.md:29 shows for sub by analogy).
#[test]
fn protein_allele_two_deletions_one_arm() {
    pin_accept(
        "NP_003997.1:p.[Trp4del;Lys23del]",
        "deletion.md silent; alleles.md:29 by analogy",
    );
}

// =====================================================================
// , (comma) separator on RNA + protein — multiple transcripts/products
// =====================================================================

/// rna-alleles:S-S17 — `,` × phase-unknown `(;)` combination is silent.
#[test]
fn rna_comma_separator_with_phase_unknown() {
    pin_either(
        "NM_004006.2:r.[897u>g(,)832_960del]",
        "RNA/alleles.md — (,) combination silent",
    );
}

/// protein-alleles:S-S15 — three-or-more outputs in `,` form.
#[test]
fn protein_comma_separator_three_outputs() {
    pin_either(
        "NP_003997.1:p.[(Lys31Asn,Val25_Lys31del,Ser40Arg)]",
        "protein/alleles.md:77 + substitution.md:74 show 2-output; 3+ silent",
    );
}

// =====================================================================
// Retracted forms (historical)
// =====================================================================

/// dna-repeated:C-C13 — `c.IVS#` form is retracted (numbering.md:107-113;
/// checklist.md:24).
#[test]
fn ivs_form_retracted() {
    pin_reject(
        "NM_004006.2:c.IVS4-2A>G",
        "numbering.md:107-113 — c.IVS# is retracted; should reject",
    );
}

/// rna-inversion:C-C6 — historic `o` opposite-strand insertion marker
/// retracted; canonical is `inv` payload form.
#[test]
fn rna_o_marker_retracted() {
    // The historic `o` opposite-strand marker is retracted, so ferro
    // should ultimately reject this. It currently still parses (the
    // marker is absorbed into the delins payload) — pin that over-
    // acceptance and promote to `pin_reject` once the marker is refused.
    pin_accept(
        "NM_004006.2:r.124_500delinsoAB053210.2:r.1289-365_1289-73",
        "RNA/inversion.md — historic 'o' marker retracted; currently still parsed",
    );
}
