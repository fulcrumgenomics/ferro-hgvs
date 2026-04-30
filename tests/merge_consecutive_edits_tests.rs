//! Tests for HGVS-spec consecutive-edit merging in alleles.
//! See docs/superpowers/specs/2026-04-30-merge-consecutive-allele-edits-design.md.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

fn normalize_to_string(input: &str) -> String {
    let normalizer = Normalizer::new(MockProvider::new());
    let variant = parse_hgvs(input).expect("parse failed");
    let normalized = normalizer.normalize(&variant).expect("normalize failed");
    format!("{}", normalized)
}

#[test]
fn test_merge_consecutive_subs_genome() {
    // Issue #72 example: two adjacent SNVs collapse to a single delins.
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000G>A;1001A>C]"),
        "NC_000001.11:g.1000_1001delinsAC",
    );
}

#[test]
fn test_merge_consecutive_dels_genome() {
    // Issue #72 example: two adjacent single-nt deletions become a single ranged del.
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000del;1001del]"),
        "NC_000001.11:g.1000_1001del",
    );
}

#[test]
fn test_merge_sub_then_del() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000G>A;1001del]"),
        "NC_000001.11:g.1000_1001delinsA",
    );
}

#[test]
fn test_merge_del_then_sub() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000del;1001A>C]"),
        "NC_000001.11:g.1000_1001delinsC",
    );
}

#[test]
fn test_merge_dels_drops_explicit_sequence() {
    // Per design doc: del+del with explicit ref sequences emits the no-sequence form.
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000delA;1001delC]"),
        "NC_000001.11:g.1000_1001del",
    );
}

#[test]
fn test_merge_dels_drops_length() {
    // Per design doc: del+del with length specifiers emits no-sequence/no-length form.
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000_1002del3;1003_1004del2]"),
        "NC_000001.11:g.1000_1004del",
    );
}

#[test]
fn test_merge_sub_then_delins() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000G>A;1001_1002delinsTT]"),
        "NC_000001.11:g.1000_1002delinsATT",
    );
}

#[test]
fn test_merge_delins_then_sub() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000_1001delinsTT;1002A>C]"),
        "NC_000001.11:g.1000_1002delinsTTC",
    );
}

#[test]
fn test_merge_delins_then_delins() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000_1001delinsTT;1002_1003delinsAA]"),
        "NC_000001.11:g.1000_1003delinsTTAA",
    );
}

#[test]
fn test_merge_skips_non_literal_delins() {
    // delins with a non-Literal payload (e.g., delins10) is not safely
    // concatenable; the design doc requires the pair to pass through.
    let input = "NC_000001.11:g.[1000G>A;1001_1010delins10]";
    let result = normalize_to_string(input);
    // The output must still contain both edits separately (unchanged).
    assert!(result.contains("1000G>A"), "expected 1000G>A in {}", result);
    assert!(
        result.contains("1001_1010delins"),
        "expected 1001_1010delins in {}",
        result
    );
    assert!(result.contains(';'), "expected separator in {}", result);
}

#[test]
fn test_merge_sub_then_ins() {
    // Spec FAQ analogue (in g. context):
    // sub at 100 + ins between 100 and 101 -> delins at 100 with alt = sub.alt + ins.bases.
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[100G>A;100_101insTA]"),
        "NC_000001.11:g.100delinsATA",
    );
}

#[test]
fn test_merge_ins_then_sub() {
    // Mirror: ins between 100 and 101 + sub at 101 -> delins at 101 with alt = ins.bases + sub.alt.
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[100_101insTA;101G>A]"),
        "NC_000001.11:g.101delinsTAA",
    );
}

#[test]
fn test_merge_del_then_ins() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[100del;100_101insTA]"),
        "NC_000001.11:g.100delinsTA",
    );
}

#[test]
fn test_merge_ins_then_del() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[100_101insTA;101del]"),
        "NC_000001.11:g.101delinsTA",
    );
}

#[test]
fn test_merge_two_ins_same_boundary_preserves_input_order() {
    // Two ins at the same boundary p|p+1 collapse into a single ins; bases
    // are concatenated in input order (T then A -> TA).
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[100_101insT;100_101insA]"),
        "NC_000001.11:g.100_101insTA",
    );
}

#[test]
fn test_merge_skips_non_literal_ins() {
    // Ins with a non-Literal payload (e.g., ins10) is not safely
    // concatenable; the pair passes through unchanged.
    let input = "NC_000001.11:g.[100G>A;100_101ins10]";
    let result = normalize_to_string(input);
    assert!(result.contains("100G>A"), "expected 100G>A in {}", result);
    assert!(
        result.contains("100_101ins10"),
        "expected 100_101ins10 in {}",
        result
    );
    assert!(result.contains(';'), "expected separator in {}", result);
}

fn provider_with_simple_transcript() -> MockProvider {
    use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
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

fn normalize_with_provider(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs(input).expect("parse failed");
    let normalized = normalizer.normalize(&variant).expect("normalize failed");
    format!("{}", normalized)
}

#[test]
fn test_merge_cds_consecutive_subs() {
    assert_eq!(
        normalize_with_provider(
            provider_with_simple_transcript(),
            "NM_TEST.1:c.[10A>G;11A>C]",
        ),
        "NM_TEST.1:c.10_11delinsGC",
    );
}

#[test]
fn test_merge_tx_consecutive_subs() {
    // n. compact form isn't accepted by the parser; use expanded form.
    assert_eq!(
        normalize_with_provider(
            provider_with_simple_transcript(),
            "[NM_TEST.1:n.10A>G;NM_TEST.1:n.11A>C]",
        ),
        "NM_TEST.1:n.10_11delinsGC",
    );
}

#[test]
fn test_merge_rna_consecutive_subs_lowercase() {
    // RNA uses lowercase nucleotides per HGVS spec; merged alt must preserve case.
    assert_eq!(
        normalize_with_provider(
            provider_with_simple_transcript(),
            "NM_TEST.1:r.[10a>g;11a>c]",
        ),
        "NM_TEST.1:r.10_11delinsgc",
    );
}

#[test]
fn test_merge_mt_consecutive_subs() {
    // m. compact form isn't accepted by the parser; use expanded form.
    assert_eq!(
        normalize_to_string("[NC_012920.1:m.100G>A;NC_012920.1:m.101A>C]"),
        "NC_012920.1:m.100_101delinsAC",
    );
}

// =====================================================================
// Negative cases — must round-trip unchanged.
// =====================================================================

#[test]
fn test_no_merge_one_nt_gap() {
    // One unchanged nucleotide between variants -> spec keeps them separate.
    let result = normalize_to_string("NC_000001.11:g.[100G>A;102C>T]");
    assert!(result.contains("100G>A"), "got {}", result);
    assert!(result.contains("102C>T"), "got {}", result);
    assert!(result.contains(';'), "got {}", result);
}

#[test]
fn test_no_merge_different_accessions() {
    let result = normalize_to_string("[NC_000001.11:g.100G>A;NC_000002.11:g.101A>C]");
    assert!(result.contains("NC_000001.11"), "got {}", result);
    assert!(result.contains("NC_000002.11"), "got {}", result);
}

#[test]
fn test_no_merge_different_variant_types() {
    // Genome and Cds in the same allele bracket -> not mergeable.
    let result = normalize_with_provider(
        provider_with_simple_transcript(),
        "[NC_000001.11:g.100G>A;NM_TEST.1:c.10A>C]",
    );
    assert!(result.contains("g.100G>A"), "got {}", result);
    assert!(result.contains("c.10A>C"), "got {}", result);
}

#[test]
fn test_no_merge_trans_phase() {
    // [a];[b] (semicolon between bracket pairs) is trans phase per HGVS.
    // The compact form g.[a];[b] isn't accepted by the parser; use expanded form.
    let result = normalize_to_string("[NC_000001.11:g.100G>A];[NC_000001.11:g.101A>C]");
    assert!(result.contains("100G>A"), "got {}", result);
    assert!(result.contains("101A>C"), "got {}", result);
    // No delins should appear.
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_no_merge_intronic_position() {
    // Intronic positions (non-zero offset) are excluded from the merge pass.
    let result = normalize_with_provider(
        provider_with_simple_transcript(),
        "NM_TEST.1:c.[10+1A>G;10+2T>G]",
    );
    assert!(result.contains("10+1A>G"), "got {}", result);
    assert!(result.contains("10+2T>G"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_no_merge_utr_boundary() {
    // c.-1 and c.1 are physically adjacent but no valid HGVS range syntax
    // spans the 5'UTR / CDS boundary (c.-1_1 doesn't exist).
    let result = normalize_with_provider(
        provider_with_simple_transcript(),
        "NM_TEST.1:c.[-1A>G;1A>T]",
    );
    assert!(result.contains("-1A>G"), "got {}", result);
    assert!(result.contains("1A>T"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_no_merge_uncertain_edit() {
    // Mu::Uncertain (paren-wrapped edit) is not mergeable.
    // The g. uncertain-edit syntax isn't supported by the parser, but c.
    // accepts it via the expanded allele form.
    let result = normalize_with_provider(
        provider_with_simple_transcript(),
        "[NM_TEST.1:c.(10A>G);NM_TEST.1:c.11A>C]",
    );
    assert!(result.contains("10(A>G)"), "got {}", result);
    assert!(result.contains("11A>C"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_no_merge_duplication_adjacent_to_sub() {
    let result = normalize_to_string("NC_000001.11:g.[100dup;101A>C]");
    assert!(result.contains("100dup"), "got {}", result);
    assert!(result.contains("101A>C"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_no_merge_inversion_adjacent_to_sub() {
    let result = normalize_to_string("NC_000001.11:g.[100_102inv;103A>C]");
    assert!(result.contains("100_102inv"), "got {}", result);
    assert!(result.contains("103A>C"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_no_merge_two_ins_different_boundaries() {
    // Two ins separated by an unchanged nucleotide at position 101 -> no merge.
    let result = normalize_to_string("NC_000001.11:g.[100_101insT;101_102insA]");
    assert!(result.contains("100_101insT"), "got {}", result);
    assert!(result.contains("101_102insA"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_no_merge_reverse_input_order() {
    // Input listed in reverse position order; the walk preserves input order
    // and checks input-order adjacency (1001's anchor end + 1 != 1000's start),
    // so no merge happens. Pins the no-resort decision.
    let result = normalize_to_string("NC_000001.11:g.[1001A>C;1000G>A]");
    assert!(result.contains("1001A>C"), "got {}", result);
    assert!(result.contains("1000G>A"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

// =====================================================================
// Chains and mixed-edit-type sequences.
// =====================================================================

#[test]
fn test_merge_chain_three_consecutive_subs() {
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000G>A;1001A>C;1002C>T]"),
        "NC_000001.11:g.1000_1002delinsACT",
    );
}

#[test]
fn test_merge_chain_sub_ins_sub_at_shared_boundary() {
    // Chain across one shared boundary: sub at 100 + ins between 100 and 101 + sub at 101.
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[100G>A;100_101insT;101A>C]"),
        "NC_000001.11:g.100_101delinsATC",
    );
}

#[test]
fn test_merge_long_chain_mixed_types() {
    // 4-element chain spanning sub, del, sub, sub at consecutive positions.
    // Position 100 sub, 101 del, 102 sub, 103 sub -> single delins 100..=103 alt = A_AA = AAA.
    // (The del at 101 contributes no alt, so the merged alt is "A" + "" + "A" + "A" = "AAA".)
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[100G>A;101del;102C>A;103T>A]"),
        "NC_000001.11:g.100_103delinsAAA",
    );
}

#[test]
fn test_merge_same_position_twice_no_merge() {
    // Two edits at the SAME position (zero gap, but overlap not adjacency)
    // must not collapse into a single delins.
    let result = normalize_to_string("NC_000001.11:g.[100G>A;100A>C]");
    // The output should still contain both — we don't synthesize a "double mutation" form.
    assert!(result.contains("100G>A"), "got {}", result);
    assert!(result.contains("100A>C"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
#[ignore = "Codon-frame exception for c. — tracked in issue #79"]
fn test_codon_frame_exception_deferred() {
    // Per HGVS spec: c.[145C>T;147C>G] (positions 145-147 share codon 49)
    // must be expressed as c.145_147delinsTGG. Implementing this needs
    // reference-provider access for the unchanged middle nucleotide and
    // codon-frame logic; deferred to issue #79.
    let provider = provider_with_simple_transcript();
    assert_eq!(
        normalize_with_provider(provider, "NM_TEST.1:c.[10A>G;12A>C]"),
        "NM_TEST.1:c.10_12delinsGAC",
    );
}
