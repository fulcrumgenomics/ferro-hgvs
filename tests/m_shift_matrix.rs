//! Mitochondrial (`m.`) 3' shift coverage matrix — issue #210.
//!
//! Mirrors `tests/del_shift_matrix.rs` (#81 A5), `tests/dup_shift_matrix.rs`
//! (#81 A6), and `tests/ins_shift_matrix.rs` (#81 A1/A7) but exercises
//! the `m.` (mitochondrial) coordinate system.
//!
//! The mito genome is genomic-style and plus-strand-only — the
//! `normalize_mt` path now runs the same window-based shuffle +
//! repeat-notation canonicalization as `normalize_genome`. Mito is
//! exempt from the codon-frame `unit_len % 3 == 0` restriction (the
//! HGVS spec's mito chapter does not carry the c.-context clause), so
//! homopolymer dup/del/ins emit `A[N]` notation just like `g.`.
//!
//! Out of scope here: origin-crossing (wraparound) variants like
//! `m.16569_1del`. Those are rejected at parse time today and are
//! tracked under #129 (F1) — see `tests/mito_circular_audit.rs`.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};
use rstest::rstest;

const MT_ACCESSION: &str = "NC_012920.1";

/// 1-based HGVS position of the first base of the core region. Mirrors
/// `common::synthetic::PAD_OFFSET + 1 = 257` but is redefined here so
/// the test file is self-contained (and so `add_genomic_sequence` is
/// called with the actual mitochondrial accession).
const C0: u64 = 257;

/// 256 bp of `ACGT` padding on each side keeps the normalizer's 100 bp
/// shuffle window in range on either side of the core tract.
const PAD: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
     ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
     ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\
     ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

/// Build a `MockProvider` containing the mitochondrial accession with
/// the given core sequence sandwiched between `PAD` on each side.
fn provider(core: &str) -> MockProvider {
    let seq = format!("{}{}{}", PAD, core, PAD);
    let mut p = MockProvider::new();
    p.add_genomic_sequence(MT_ACCESSION, seq);
    p
}

/// Substitute 1-based core positions into an HGVS template — `{N}` is
/// replaced with `(C0 + N - 1)`. Matches the helper used in
/// `del_shift_matrix.rs::genomic`.
fn hgvs(template: &str, args: &[u64]) -> String {
    let mut s = template.to_string();
    for (i, p) in args.iter().enumerate() {
        s = s.replace(&format!("{{{}}}", i), &(C0 + p - 1).to_string());
    }
    s
}

/// Normalize one variant string and return the formatted output. Panics
/// on parse / normalize failure with the input for context.
fn normalize_to_string(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let variant =
        parse_hgvs(input).unwrap_or_else(|e| panic!("Failed to parse {:?}: {}", input, e));
    let normalized = normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("Normalize failed for {:?}: {}", input, e));
    format!("{}", normalized)
}

// =============================================================================
// Substitution: passes through unchanged (no shuffle invariant for SNVs).
// =============================================================================

#[test]
fn sub_passes_through_unchanged() {
    // Plain SNV in a non-repeat region — normalize is a no-op.
    let p = provider("ACGTACGTACGTACGT");
    let result = normalize_to_string(p, &format!("{}:{}", MT_ACCESSION, hgvs("m.{0}A>G", &[5])));
    assert_eq!(
        result,
        format!("{}:{}", MT_ACCESSION, hgvs("m.{0}A>G", &[5]))
    );
}

// =============================================================================
// Deletion: 3' shift and `unit[N-k]` canonicalization.
// =============================================================================

#[test]
fn single_base_del_no_shift_when_neighbors_differ() {
    // Core: ACGTACGT. Delete 'G' at core 3. ref[del_start]='G',
    // ref[del_end]='T' → mismatch, no shift, stays as del.
    let p = provider("ACGTACGT");
    let result = normalize_to_string(p, &format!("{}:{}", MT_ACCESSION, hgvs("m.{0}del", &[3])));
    assert_eq!(
        result,
        format!("{}:{}", MT_ACCESSION, hgvs("m.{0}del", &[3]))
    );
}

#[test]
fn single_base_del_in_homopolymer_shifts_three_prime() {
    // Core: ACAAAACG. A-tract at core 3..6. Delete A at core 3.
    // Shifts to 3' end (core 6); k=1 stays as del per the existing
    // tests/del_shift_matrix.rs::genomic::single_base_del_in_homopolymer.
    let p = provider("ACAAAACG");
    let result = normalize_to_string(p, &format!("{}:{}", MT_ACCESSION, hgvs("m.{0}del", &[3])));
    assert_eq!(
        result,
        format!("{}:{}", MT_ACCESSION, hgvs("m.{0}del", &[6]))
    );
}

#[rstest]
// Case 1: AC tandem (ACAC at core 3..6). Del 1 AC at core 3-4 → del at core 5-6.
#[case("ACACACGT", "m.{0}_{1}del", &[3u64, 4u64], "m.{0}_{1}del", &[5u64, 6u64])]
// Case 2: GCA tandem (GCAGCAGCA at core 3..11). Del 1 GCA at core 3-5 → del at core 9-11.
#[case("ACGCAGCAGCATG", "m.{0}_{1}del", &[3u64, 5u64], "m.{0}_{1}del", &[9u64, 11u64])]
fn multi_base_del_in_tandem_one_unit(
    #[case] core: &str,
    #[case] in_template: &str,
    #[case] in_args: &[u64],
    #[case] out_template: &str,
    #[case] out_args: &[u64],
) {
    let p = provider(core);
    let input = format!("{}:{}", MT_ACCESSION, hgvs(in_template, in_args));
    let expected = format!("{}:{}", MT_ACCESSION, hgvs(out_template, out_args));
    assert_eq!(normalize_to_string(p, &input), expected);
}

#[rstest]
// Case 1: 2 A's from 5-A homopolymer → A[3].
#[case("ACAAAAACG", "m.{0}_{1}del", &[3u64, 4u64], "m.{0}_{1}A[3]", &[3u64, 7u64])]
// Case 2: 2 GCAs from 3-GCA tract → GCA[1].
#[case("ACGCAGCAGCATG", "m.{0}_{1}del", &[3u64, 8u64], "m.{0}_{1}GCA[1]", &[3u64, 11u64])]
// Case 3: cyclic rotation — del CAGCAG (r=1 of GCA) → GCA[1] via
// shuffle phase-alignment lemma.
#[case("TTGCAGCAGCATT", "m.{0}_{1}del", &[4u64, 9u64], "m.{0}_{1}GCA[1]", &[3u64, 11u64])]
fn multi_base_del_in_tandem_multiple_units_emits_repeat(
    #[case] core: &str,
    #[case] in_template: &str,
    #[case] in_args: &[u64],
    #[case] out_template: &str,
    #[case] out_args: &[u64],
) {
    let p = provider(core);
    let input = format!("{}:{}", MT_ACCESSION, hgvs(in_template, in_args));
    let expected = format!("{}:{}", MT_ACCESSION, hgvs(out_template, out_args));
    assert_eq!(normalize_to_string(p, &input), expected);
}

// =============================================================================
// Duplication: 3' shift, `unit[N+k]` canonicalization, single-copy stays as dup.
// =============================================================================

#[test]
fn single_base_dup_no_shift_when_neighbors_differ() {
    let p = provider("ACGTACGT");
    let result = normalize_to_string(p, &format!("{}:{}", MT_ACCESSION, hgvs("m.{0}dup", &[3])));
    assert_eq!(
        result,
        format!("{}:{}", MT_ACCESSION, hgvs("m.{0}dup", &[3]))
    );
}

#[test]
fn single_base_dup_in_homopolymer_shifts_three_prime() {
    // A-tract at core 3..6. Dup A at core 3 → shifts to 3' end of
    // tract (core 6); dup_len=1 stays as dup (matches the existing
    // dup_shift_matrix.rs convention for single-base dup).
    let p = provider("ACAAAACG");
    let result = normalize_to_string(p, &format!("{}:{}", MT_ACCESSION, hgvs("m.{0}dup", &[3])));
    assert_eq!(
        result,
        format!("{}:{}", MT_ACCESSION, hgvs("m.{0}dup", &[6]))
    );
}

#[rstest]
// Case 1: 2 A's dup in 5-A homopolymer → A[7].
#[case("ACAAAAACG", "m.{0}_{1}dup", &[3u64, 4u64], "m.{0}_{1}A[7]", &[3u64, 7u64])]
// Case 2: 2 GCAs dup in 3-GCA tract → GCA[5].
#[case("ACGCAGCAGCATG", "m.{0}_{1}dup", &[3u64, 8u64], "m.{0}_{1}GCA[5]", &[3u64, 11u64])]
// Case 3: cyclic rotation — dup CAGCAG (r=1 of GCA, 2 copies) → GCA[5].
#[case("TTGCAGCAGCATT", "m.{0}_{1}dup", &[4u64, 9u64], "m.{0}_{1}GCA[5]", &[3u64, 11u64])]
fn multi_base_dup_in_tandem_multiple_units_emits_repeat(
    #[case] core: &str,
    #[case] in_template: &str,
    #[case] in_args: &[u64],
    #[case] out_template: &str,
    #[case] out_args: &[u64],
) {
    let p = provider(core);
    let input = format!("{}:{}", MT_ACCESSION, hgvs(in_template, in_args));
    let expected = format!("{}:{}", MT_ACCESSION, hgvs(out_template, out_args));
    assert_eq!(normalize_to_string(p, &input), expected);
}

// =============================================================================
// Insertion: ins matching adjacent repeat unit → `[N+1]` increment
// (B1); cyclic-rotation single-copy ins → dup at most-3' position
// (#155 / cyclic_rotation_ins_one_unit pattern).
// =============================================================================

#[test]
fn single_base_ins_matching_homopolymer_becomes_dup() {
    // A-tract A[4] at core 3..6. Insert one A at core 2_3 (just
    // before the tract). Per `insertion_to_repeat`, a single-copy
    // ins of a unit already adjacent canonicalizes to a `dup` at
    // the most-3' matching position, not to `A[N+1]` repeat
    // notation (the [N+k] form requires k >= 2 added unit copies).
    // Same convention as the existing `tests/ins_shift_matrix.rs`
    // single-base ins cases.
    let p = provider("ACAAAACG");
    let input = format!("{}:{}", MT_ACCESSION, hgvs("m.{0}_{1}insA", &[2u64, 3u64]));
    let expected = format!("{}:{}", MT_ACCESSION, hgvs("m.{0}dup", &[6u64]));
    assert_eq!(normalize_to_string(p, &input), expected);
}

#[test]
fn ins_one_unit_aligned_becomes_dup() {
    // AC tandem AC[2] at core 3..6. Insert one AC at core 2_3
    // (matches adjacent unit). After 3'-shift it canonicalizes to a
    // dup at the most-3' AC. Mirrors the B1 / A1 pattern for
    // unit_len=2.
    let p = provider("ACACACGT");
    let input = format!("{}:{}", MT_ACCESSION, hgvs("m.{0}_{1}insAC", &[2u64, 3u64]));
    let expected = format!("{}:{}", MT_ACCESSION, hgvs("m.{0}_{1}dup", &[5u64, 6u64]));
    assert_eq!(normalize_to_string(p, &input), expected);
}

#[test]
fn ins_two_units_aligned_becomes_repeat_increment() {
    // AC tandem AC[4] at core 1..8 (the entire core minus the
    // trailing GT). Insert two ACs (`insACAC`) just before the
    // tract; after 3'-shift `insertion_to_repeat` emits AC[6] over
    // the reference-tract extent (the full 4-AC run at core 1..8).
    let p = provider("ACACACACGT");
    let input = format!(
        "{}:{}",
        MT_ACCESSION,
        hgvs("m.{0}_{1}insACAC", &[2u64, 3u64])
    );
    let expected = format!("{}:{}", MT_ACCESSION, hgvs("m.{0}_{1}AC[6]", &[1u64, 8u64]));
    assert_eq!(normalize_to_string(p, &input), expected);
}

// =============================================================================
// Round-trip equivalence: dup-form and repeat-form of the same edit
// must canonicalize identically. Sanity check that the m. path matches
// the g. path's behavior on equivalent input.
// =============================================================================

#[test]
fn dup_form_and_repeat_form_canonicalize_identically() {
    // 4-A homopolymer, 2 A's added. `m.{0}_{1}dup` and `m.{0}_{1}A[6]`
    // describe the same allele and must normalize to the same canonical
    // form.
    let p1 = provider("ACAAAACG");
    let p2 = provider("ACAAAACG");
    let r1 = normalize_to_string(
        p1,
        &format!("{}:{}", MT_ACCESSION, hgvs("m.{0}_{1}dup", &[3u64, 4u64])),
    );
    let r2 = normalize_to_string(
        p2,
        &format!("{}:{}", MT_ACCESSION, hgvs("m.{0}_{1}A[6]", &[3u64, 6u64])),
    );
    assert_eq!(r1, r2, "dup-form and repeat-form must agree");
}

// =============================================================================
// Fallback path: no provider data → minimal-notation canonicalization,
// not a panic. Preserves the existing behavior pinned by
// `mito_circular_audit.rs::audit_mt_non_wrapping_sub_normalizes_to_self`
// and friends.
// =============================================================================

#[test]
fn no_provider_data_falls_back_to_minimal_notation() {
    // Empty provider — `normalize_mt` falls back to
    // `canonicalize_mt_variant` (the genome-variant analog). A `dup`
    // with explicit sequence has the sequence stripped to minimal
    // notation but otherwise round-trips.
    let normalizer = Normalizer::new(MockProvider::new());
    let variant = parse_hgvs("NC_012920.1:m.100_102dupACG").expect("parse");
    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize should not error without provider data");
    // Minimal notation drops the explicit sequence in `dup`.
    assert_eq!(format!("{}", normalized), "NC_012920.1:m.100_102dup");
}

#[test]
fn no_provider_data_substitution_passes_through() {
    // Substitutions short-circuit via `needs_normalization` before the
    // window-fetch / fallback paths, so they round-trip even with an
    // empty provider.
    let normalizer = Normalizer::new(MockProvider::new());
    let result = parse_hgvs("NC_012920.1:m.100A>G")
        .and_then(|v| normalizer.normalize(&v))
        .expect("substitution passes through normalize");
    assert_eq!(format!("{}", result), "NC_012920.1:m.100A>G");
}

#[test]
fn decorated_position_with_offset_routes_to_canonicalize_fallback() {
    // `m.` positions can carry intronic-style offsets in HGVS (e.g.
    // `m.100+1`). Window-based normalization remaps from `pos.base`
    // and reconstructs with `GenomePos::new`, which would silently
    // drop the offset. `normalize_mt` must detect decorated positions
    // and route to the canonicalize-only fallback so the offset is
    // preserved verbatim — *even when reference data is present*, so
    // the guard fires before `get_sequence`. Use an indel form
    // (substitutions short-circuit before the guard).
    let normalizer = Normalizer::new(provider("ACGT"));
    let variant =
        parse_hgvs(&format!("{}:m.{}+1_{}+1delAC", MT_ACCESSION, C0, C0 + 1)).expect("parse");
    let normalized = normalizer
        .normalize(&variant)
        .expect("decorated m. positions must not panic in normalize");
    // Minimal-notation cleanup strips the explicit `del` sequence but
    // preserves the offsets — they are not collapsed to plain bases.
    assert_eq!(
        format!("{}", normalized),
        format!("{}:m.{}+1_{}+1del", MT_ACCESSION, C0, C0 + 1)
    );
}

// =============================================================================
// Delins: confirms the window-shuffle path fires for `delins` too,
// not just `del`/`dup`/`ins`. Covers the edit-priority rewrites that
// `canonicalize_delins` applies (e.g. same-ref `delins` → `=`, ins-only
// `delins` → `ins`).
// =============================================================================

#[test]
fn delins_same_ref_collapses_to_identity() {
    // Core "ACGT" at core 1..4. `m.<C0+1>_<C0+2>delinsCG` replaces
    // ref "CG" with "CG" — same bases → identity (`=`). Mirrors the
    // genomic edit-type priority rule (sub > del > inv > dup > ins);
    // confirms the same path runs for `m.`.
    let p = provider("ACGTACGT");
    let input = format!(
        "{}:{}",
        MT_ACCESSION,
        hgvs("m.{0}_{1}delinsCG", &[2u64, 3u64])
    );
    let expected = format!("{}:{}", MT_ACCESSION, hgvs("m.{0}_{1}=", &[2u64, 3u64]));
    assert_eq!(normalize_to_string(p, &input), expected);
}

#[test]
fn delins_rewrites_to_insertion_then_canonicalizes_to_dup() {
    // A-tract A[4] at core 3..6. `delinsAAA` over core 3-4 replaces
    // "AA" with "AAA"; shared-affix trimming reduces to a 1-base
    // insertion of "A". Per HGVS spec the canonical short form for an
    // `A` insertion into an existing A-tract is `dup` (the inserted
    // base duplicates the adjacent reference base under the 3'-rule).
    // Closes-after #356: the delins → ins canonicalization now
    // recurses into `normalize_na_edit` so the full Insertion pipeline
    // (3'-shuffle + `insertion_to_duplication` recognizer) runs.
    let p = provider("ACAAAACG");
    let input = format!(
        "{}:{}",
        MT_ACCESSION,
        hgvs("m.{0}_{1}delinsAAA", &[3u64, 4u64])
    );
    // 3'-shuffle lands the resulting `insA` at the 3'-most position of
    // the A-tract (core position 6 = `C0 + 6 - 1 = 262`), where it
    // duplicates the preceding `A`. Canonical form is `m.262dup`.
    let expected = format!("{}:{}", MT_ACCESSION, hgvs("m.{0}dup", &[6u64]));
    assert_eq!(normalize_to_string(p, &input), expected);
}
