//! Codon-frame exception extensions (issue #275, follow-up to #79 / #104).
//!
//! PR #104 implemented the HGVS codon-frame exception for `c.` SNV pairs
//! separated by one unchanged nucleotide. Three sibling cases were flagged
//! as out of scope at the time:
//!
//! 1. `r.` coding regions — RNA coordinates use the same CDS-relative axis
//!    as `c.`, and the same "two variants together affecting one amino acid"
//!    rule applies. The merge pass should accept `Region::Rna` alongside
//!    `Region::Cds`.
//! 2. Chains involving 3+ SNVs across a codon — the original pair-detector
//!    required `prev_a` to still be a single-base SUB anchor. A chain of
//!    `sub`+`sub` (strict adjacency) followed by a gap-of-1 sub in the same
//!    codon should still trigger the carve-out, even though `prev_a` has
//!    already grown into a 2-base delins via the first merge.
//! 3. `sub`+`del` pairs separated by one unchanged nucleotide — the pair
//!    rule applied only when both endpoints were single-base subs. A del
//!    on either side of a gap-of-1 pair within one codon still satisfies
//!    the spec rationale ("together affecting one amino acid").
//!
//! Most tests use the same shared `simple_transcript` fixture as the existing
//! `merge_consecutive_edits_tests.rs`, with `cds_start = 1` so `c.N` and
//! `r.N` index the same transcript byte. Reference bases (1-indexed): the
//! transcript sequence is `ATGCAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGT...`,
//! i.e. `c.1=A`, `c.2=T`, `c.3=G`, `c.4=C`, `c.5=A`, `c.6=A`, `c.7=A`,
//! `c.8=A`, `c.9=A`, `c.10=C`, `c.11=C`, `c.12=C`, `c.13=C`, `c.14=C`,
//! `c.15=G`, ...
//!
//! The minus-strand test uses its own ad-hoc provider so the codon-frame
//! relaxation can be pinned on a `Strand::Minus` fixture — see
//! `provider_simple_minus_strand` for the rationale.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Shared plus-strand fixture (`cds_start = 1`) used by almost every test
/// in this file. Naming follows the `provider_simple()` / `normalize_with()`
/// pattern from `tests/issue_165_delins_sub_only_decompose.rs`.
fn provider_simple() -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence: String =
        "ATGCAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGT".to_string();
    let len = sequence.len() as u64;
    let exons = vec![Exon::new(1, 1, len)];
    let transcript = Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        sequence,
        Some(1),
        Some(60),
        exons,
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    provider
}

/// Minus-strand fixture with the same transcript-orientation sequence as
/// `provider_simple()`. `lookup_codon_middle_ref` reads `tx.sequence`
/// directly (transcript 5'→3' orientation) and never goes through the
/// genome, so the codon-frame merge fires regardless of strand. This test
/// pins that property — a parallel sibling to the plus-strand tests below.
fn provider_simple_minus_strand() -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence: String =
        "ATGCAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGT".to_string();
    let len = sequence.len() as u64;
    let exons = vec![Exon::new(1, 1, len)];
    let transcript = Transcript::new(
        "NM_MINUS.1".to_string(),
        Some("MINUS".to_string()),
        Strand::Minus,
        sequence,
        Some(1),
        Some(60),
        exons,
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    provider
}

/// Fixture with `cds_start = 10` so `c.-N` (5'UTR) and `c.*N` (3'UTR)
/// positions are well-defined. Transcript-byte layout:
///   tx pos:    1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 ...
///   c. axis:  -9 -8 -7 -6 -5 -4 -3 -2 -1 1 2 3 4 5 ...
///   base:      A T G C A A A A A  C  C  C  C  C ...
/// CDS span = c.1..c.51 (cds_start=10, cds_end=60). 3'UTR begins at `c.*1`.
fn provider_with_utr() -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence: String =
        "ATGCAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGT".to_string();
    let len = sequence.len() as u64;
    let exons = vec![Exon::new(1, 1, len)];
    let transcript = Transcript::new(
        "NM_UTR.1".to_string(),
        Some("UTR".to_string()),
        Strand::Plus,
        sequence,
        Some(10),
        Some(60),
        exons,
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    provider
}

fn normalize_with(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs(input).expect("parse failed");
    let normalized = normalizer.normalize(&variant).expect("normalize failed");
    format!("{}", normalized)
}

// =====================================================================
// Item 1: codon-frame exception for `r.` coding regions.
//
// `r.` positions are CDS-relative (same axis as `c.`), so the same
// codon-frame exception applies. The merged delins must preserve the
// RNA-display convention (lowercase nucleotides).
// =====================================================================

#[test]
fn codon_frame_rna_two_subs_one_codon() {
    // r.10 = C, r.11 = C (middle ref), r.12 = C. Codon = 4 (positions
    // 10-12). User-declared alts `g` and `a` are preserved. The middle
    // ref base comes from the transcript sequence; r. display rules
    // emit it lowercase.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:r.[10c>g;12c>a]"),
        "NM_TEST.1:r.10_12delinsgca",
    );
}

#[test]
fn no_codon_frame_rna_pair_straddles_codon_boundary() {
    // r.3 sits in codon 1 (positions 1-3); r.5 sits in codon 2
    // (positions 4-6). Different codons → no merge.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:r.[3g>u;5a>c]"),
        "NM_TEST.1:r.[3g>u;5a>c]",
    );
}

#[test]
fn no_codon_frame_rna_gap_of_two() {
    // Gap of 2 nt is beyond the spec carve-out, which is exactly one.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:r.[10c>g;13c>a]"),
        "NM_TEST.1:r.[10c>g;13c>a]",
    );
}

// =====================================================================
// Item 2: chains of 3+ SNVs across a codon.
//
// Original guard required `prev_a.alt.len() == 1`, blocking codon-frame
// after a strict-adjacency merge has already grown `prev_a` into a 2+
// base delins. After relaxation, a chain like `[sub at p; sub at p+1;
// sub at p+3]` (strict 3 then gap-of-1 to same-codon) still merges into
// a single delins covering the full span, then the post-merge decompose
// pass splits it back per the codon-frame and `general.md:34` rules.
// =====================================================================

#[test]
fn codon_frame_chain_strict_then_gap_in_next_codon() {
    // Positions 3 (codon 1), 4 (codon 2), 6 (codon 2). 3+4 strict, 4+6
    // gap-of-1 same-codon. Refs: c.3=G, c.4=C, c.5=A (middle), c.6=A.
    //
    // Trace:
    //   Strict merge 3+4  →  delins[3,4] alt=TG.
    //   Codon-frame [3,4] + 6 (relaxed)  →  delins[3,6] alt=TGAT
    //     (middle ref=A injected at position 5).
    //   Post-decompose: ref=GCAA vs alt=TGAT  →  Sub@3, Sub@4, Id@5, Sub@6.
    //   build_split_variants (CDS, codon-frame aware) keeps Sub@3 alone
    //     (no Sub-Id-Sub triplet starting at i=0); at i=1 the
    //     [Sub@4; Id@5; Sub@6] triplet matches and same_codon(4,6) = codon 2,
    //     so it emits as a single 3-base delins. The trailing Sub@3 falls
    //     into the substitution run and flushes as `c.3G>T`.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:c.[3G>T;4C>G;6A>T]"),
        "NM_TEST.1:c.[3G>T;4_6delinsGAT]",
    );
}

#[test]
fn codon_frame_chain_rna_strict_then_gap_in_next_codon() {
    // RNA analogue of the c. chain test above. Same positions, lowercase
    // bases for r. display.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:r.[3g>u;4c>g;6a>u]"),
        "NM_TEST.1:r.[3g>u;4_6delinsgau]",
    );
}

#[test]
fn no_codon_frame_chain_strict_then_gap_crosses_codon() {
    // Positions 4, 5 (codon 2), 7 (codon 3). 4+5 strict, but 5+7 spans
    // codons 2 and 3 — the codon-frame exception explicitly does NOT
    // extend across codon boundaries. Result: a strict 2-base delins at
    // [4,5] plus an unmerged sub at 7.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:c.[4C>G;5A>C;7A>T]"),
        "NM_TEST.1:c.[4_5delinsGC;7A>T]",
    );
}

// =====================================================================
// Item 3: sub+del pairs separated by one unchanged nucleotide.
//
// The exception fired only for sub+sub. With one of the two edits a
// single-base del, the same spec rationale applies: two edits in one
// codon together affecting one amino acid become a delins. The merge
// extends a 1-base del (alt.len() == 0) the same way it extends a
// 1-base sub.
// =====================================================================

#[test]
fn codon_frame_sub_then_del_one_codon() {
    // c.1 = A, c.2 = T (middle), c.3 = G. Codon 1. SUB at 1 replaces A
    // with C; DEL at 3 removes the G. Merged span [1,3] with alt
    // [C(sub), T(middle ref), <empty from del>] = "CT".
    // Canonicalize: no shared affix between ATG and CT (prefix C vs A,
    // suffix T vs G), so the form stays as a 3-base→2-base delins.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:c.[1A>C;3del]"),
        "NM_TEST.1:c.1_3delinsCT",
    );
}

#[test]
fn codon_frame_del_then_sub_one_codon() {
    // DEL at 1, SUB C>A at 3 (ref G). Merged alt = [<empty>, T(middle
    // ref), A(sub)] = "TA". Ref ATG vs alt TA: prefix T vs A no, suffix
    // A vs G no — no trim. Final delins [1,3] alt "TA".
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:c.[1del;3G>A]"),
        "NM_TEST.1:c.1_3delinsTA",
    );
}

#[test]
fn codon_frame_sub_then_del_rna() {
    // RNA analogue of sub+del: r. axis matches c. axis, so the same
    // codon-frame rule applies. Lowercase RNA bases per r. display —
    // the RNA alphabet is `a/c/g/u` per HGVS spec
    // (recommendations/RNA/{substitution,delins,repeated}.md, e.g.
    // `r.142_144delinsugg`), and PR #293 wired the display-side T→u
    // canonicalization so the merged delins emits `cu`, not `ct`.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:r.[1a>c;3del]"),
        "NM_TEST.1:r.1_3delinscu",
    );
}

#[test]
fn no_codon_frame_sub_then_del_pair_straddles_codon_boundary() {
    // SUB at 12 (codon 4), DEL at 14 (codon 5) → different codons, no
    // merge. Position 14 is the last C in the CCCCC run at positions
    // 10-14, next base is G at 15, so the del does not right-shift.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:c.[12C>G;14del]"),
        "NM_TEST.1:c.[12C>G;14del]",
    );
}

#[test]
fn no_codon_frame_sub_then_del_gap_of_two() {
    // Gap of 2 nt (positions 3 apart) exceeds the spec carve-out, which
    // is exactly one unchanged base. SUB at 1 + DEL at 4 (stable: next
    // base at 5 is A). The pair passes through unchanged.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:c.[1A>C;4del]"),
        "NM_TEST.1:c.[1A>C;4del]",
    );
}

// =====================================================================
// Extra coverage for cases the relaxation now admits.
//
// These tests pin the behavior of inputs that the merge-eligibility
// check accepts but the original issue scope did not enumerate. For
// each, the assertion reflects the spec rationale ("two variants in
// one codon together affecting one amino acid → emit as a delins").
// =====================================================================

#[test]
fn codon_frame_del_then_del_one_codon() {
    // Both endpoints are single-base deletions; codon-frame eligibility
    // requires `prev_a.start <= prev_a.end` and `next.alt.len() <= 1`,
    // which a 1-base del (start == end, alt empty) satisfies on both
    // sides. c.1=A, c.2=T (middle), c.3=G — codon 1. Two dels in one
    // codon leave only the middle T, an unambiguous "two edits together
    // affecting one amino acid" case, so the spec rationale supports
    // merging. Merged alt = [<empty>, T, <empty>] = "T". Ref ATG vs alt
    // T: prefix A vs T no, suffix G vs T no — final delins [1,3] alt
    // "T". The post-merge canonicalization (issue #165 path) decomposes
    // the resulting delins only when interior identities show up between
    // two SUB edges, so a sub-only split does not fire here.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:c.[1del;3del]"),
        "NM_TEST.1:c.1_3delinsT",
    );
}

#[test]
fn codon_frame_multibase_delins_head_then_sub_one_codon() {
    // User-supplied multi-base delins as the chain head. The merge pass
    // alone would not eat it (anchor_for_variant would yield
    // `start=3, end=4`, not the 1-base shape codon-frame requires), but
    // the post-canonicalization split decomposes `c.3_4delinsCG` into
    // its component subs first — `c.3G>C` and `c.4C>G` — and the merge
    // pass then chains those with the trailing sub at c.6 (codon 2). The
    // 4+6 gap-of-1 leg is codon-frame eligible and fires. End state:
    //   strict merge 3+4 → delins[3,4] alt=CG.
    //   codon-frame [3,4] + 6 → delins[3,6] alt=CGAT.
    //   post-decompose: ref GCAA vs alt CGAT → Sub@3, Sub@4, Id@5, Sub@6.
    //   build_split_variants emits Sub@3 (`c.3G>C`), then the
    //   [Sub@4; Id@5; Sub@6] triplet in codon 2 as a 3-base delins.
    // The relaxation's spec rationale applies cleanly here: codons 1
    // (c.3 alone) and 2 (c.4..6 together) each emit in canonical form.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:c.[3_4delinsCG;6A>T]"),
        "NM_TEST.1:c.[3G>C;4_6delinsGAT]",
    );
}

// =====================================================================
// Negative coverage: intronic / UTR / wrong-region inputs must NOT
// trigger the relaxed codon-frame path.
//
// The codon-frame eligibility check in `merge_consecutive_edits` requires
// `Region::Cds` or `Region::Rna`. `simple_cds_pos` / `simple_rna_pos`
// reject intronic offsets outright (return None), and 5'UTR / 3'UTR
// positions map to `Region::FivePrimeUtr` / `Region::ThreePrimeUtr` and
// therefore fail the region match. Both rejections leave the chain
// un-merged.
// =====================================================================

#[test]
fn no_codon_frame_intronic_offset_pair() {
    // c.1+1A>T sits on an intronic offset (`base=1, offset=+1`), which
    // `simple_cds_pos` rejects via `pos.is_intronic()`. The pair never
    // becomes a merge candidate, so it passes through unchanged.
    assert_eq!(
        normalize_with(provider_simple(), "NM_TEST.1:c.[1+1A>T;3G>A]"),
        "NM_TEST.1:c.[1+1A>T;3G>A]",
    );
}

#[test]
fn no_codon_frame_five_prime_utr_pair() {
    // c.-3 and c.-1 both live in the 5'UTR. `simple_cds_pos` returns
    // `Region::FivePrimeUtr` rather than `Region::Cds`, and the
    // codon-frame eligibility check explicitly fails the
    // `Region::Cds | Region::Rna` filter. UTR positions are not part of
    // the CDS-relative codon grid, so the relaxation must not fire.
    assert_eq!(
        normalize_with(provider_with_utr(), "NM_UTR.1:c.[-3C>G;-1A>T]"),
        "NM_UTR.1:c.[-3C>G;-1A>T]",
    );
}

#[test]
fn no_codon_frame_three_prime_utr_pair() {
    // c.*1 and c.*3 are 3'UTR (`utr3=true`), which `simple_cds_pos`
    // tags as `Region::ThreePrimeUtr`. Same rationale as the 5'UTR
    // case: outside the CDS codon grid, no relaxation. Refs at the
    // matching transcript positions follow the test fixture
    // (`provider_with_utr` uses `cds_start = 10`, `cds_end = 60`, so
    // c.*1 maps to transcript index 61). Out-of-fixture refs are
    // unimportant — the assertion only requires that the codon-frame
    // path does not fire.
    assert_eq!(
        normalize_with(provider_with_utr(), "NM_UTR.1:c.[*1A>G;*3C>T]"),
        "NM_UTR.1:c.[*1A>G;*3C>T]",
    );
}

// =====================================================================
// Minus-strand transcript: the codon-frame relaxation must fire on
// `Strand::Minus` fixtures too. `lookup_codon_middle_ref` reads
// `tx.sequence` directly (transcript 5'→3' orientation), which is
// strand-independent — the merge pass does not consult genomic
// orientation at all.
// =====================================================================

#[test]
fn codon_frame_minus_strand_two_subs_one_codon() {
    // Mirror of `codon_frame_rna_two_subs_one_codon` but on a minus-
    // strand fixture. r.10 = C, r.11 = C (middle ref), r.12 = C in
    // transcript orientation. The codon-frame merge fires identically.
    assert_eq!(
        normalize_with(provider_simple_minus_strand(), "NM_MINUS.1:r.[10c>g;12c>a]",),
        "NM_MINUS.1:r.10_12delinsgca",
    );
}
