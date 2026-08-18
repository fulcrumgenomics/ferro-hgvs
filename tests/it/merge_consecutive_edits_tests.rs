//! Tests for HGVS-spec consecutive-edit merging in alleles.
//! See docs/superpowers/specs/2026-04-30-merge-consecutive-allele-edits-design.md.

use ferro_hgvs::{hgvs_to_spdi_simple, parse_hgvs, MockProvider, NormalizeConfig, Normalizer};

fn normalize_to_string(input: &str) -> String {
    let normalizer = Normalizer::new(MockProvider::new());
    let variant = parse_hgvs(input).expect("parse failed");
    let normalized = normalizer.normalize(&variant).expect("normalize failed");
    format!("{}", normalized)
}

fn normalize_to_string_and_warning_codes(input: &str) -> (String, Vec<String>) {
    let normalizer = Normalizer::new(MockProvider::new());
    let parsed = parse_hgvs(input).expect("parse");
    let r = normalizer
        .normalize_with_diagnostics(&parsed)
        .expect("normalize");
    let codes: Vec<String> = r.warnings.iter().map(|w| w.code().to_string()).collect();
    (r.result.to_string(), codes)
}

fn normalize_with_provider_and_warning_codes(
    provider: MockProvider,
    input: &str,
) -> (String, Vec<String>) {
    let normalizer = Normalizer::new(provider);
    let parsed = parse_hgvs(input).expect("parse");
    let r = normalizer
        .normalize_with_diagnostics(&parsed)
        .expect("normalize");
    let codes: Vec<String> = r.warnings.iter().map(|w| w.code().to_string()).collect();
    (r.result.to_string(), codes)
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
fn test_two_ins_same_boundary_preserved_as_overlap() {
    // Two separate insertions at the same boundary p|p+1 are an order-ambiguous
    // overlap conflict, not a merge: there is no canonical order for the two
    // inserted sequences (HGVS expresses "insert both" as a single ordered
    // compound payload `ins[T;A]`, not two members — general.md:79). mutalyzer
    // rejects this as `EOVERLAP` and VariantValidator as `AlleleSyntaxError`;
    // ferro's strict mode rejects it (W5002, #486) and the non-strict path
    // warn-and-preserves it as authored rather than fabricating an
    // order-dependent, non-idempotent `insTA` (#1004).
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[100_101insT;100_101insA]"),
        "NC_000001.11:g.[100_101insT;100_101insA]",
    );
}

#[test]
fn test_merge_skips_non_literal_ins() {
    // Ins with a non-Literal payload (e.g., insN[10]) is not safely
    // concatenable; the pair passes through unchanged.
    //
    // The payload used to be spelled `ins10`. #1789 refuses that at normalize
    // in every mode (`W3029`, `checklist.md:33`), so it can no longer stand in
    // for "a non-literal payload" here — but `insN[10]` is the spelling
    // `DNA/insertion.md:77` gives for the same claim, is equally non-literal,
    // and reaches the same `merge` arm.
    let input = "NC_000001.11:g.[100G>A;100_101insN[10]]";
    let result = normalize_to_string(input);
    assert!(result.contains("100G>A"), "expected 100G>A in {}", result);
    assert!(
        result.contains("100_101insN[10]"),
        "expected 100_101insN[10] in {}",
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
    // Positions 10-11 in the simple transcript are "CC". Use ref-matching
    // SNV inputs so the merged delins reflects the actual ref bases (under
    // issue #160 the merged delins is re-canonicalized post-merge, which
    // would otherwise rewrite a ref-mismatched delins to its actual-ref
    // canonical form).
    assert_eq!(
        normalize_with_provider(
            provider_with_simple_transcript(),
            "NM_TEST.1:c.[10C>G;11C>A]",
        ),
        "NM_TEST.1:c.10_11delinsGA",
    );
}

#[test]
fn test_merge_tx_consecutive_subs() {
    // n. compact form isn't accepted by the parser; use expanded form.
    // Same ref-matching rationale as test_merge_cds_consecutive_subs.
    assert_eq!(
        normalize_with_provider(
            provider_with_simple_transcript(),
            "[NM_TEST.1:n.10C>G;NM_TEST.1:n.11C>A]",
        ),
        "NM_TEST.1:n.10_11delinsGA",
    );
}

#[test]
fn test_merge_rna_consecutive_subs_lowercase() {
    // RNA uses lowercase nucleotides per HGVS spec; merged alt must preserve case.
    // Same ref-matching rationale as test_merge_cds_consecutive_subs.
    assert_eq!(
        normalize_with_provider(
            provider_with_simple_transcript(),
            "NM_TEST.1:r.[10c>g;11c>a]",
        ),
        "NM_TEST.1:r.10_11delinsga",
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
    // Per #241, the canonical Display wraps position+edit together,
    // mirroring the protein predicted form (`p.(Arg248Gln)`).
    assert!(result.contains("(10A>G)"), "got {}", result);
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
fn test_merge_reverse_input_order_is_order_independent() {
    // #1103: cis members are sorted into genomic order *before* the merge, so an
    // adjacent pair listed in reverse position order fires the same merge as the
    // forward order. (Before #1103 the walk preserved input order and checked
    // input-order adjacency, so `[1001A>C;1000G>A]` was left unmerged — the exact
    // order-dependence #1103 removes.)
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1001A>C;1000G>A]"),
        "NC_000001.11:g.1000_1001delinsAC",
    );
    // The forward order yields the identical canonical string.
    assert_eq!(
        normalize_to_string("NC_000001.11:g.[1000G>A;1001A>C]"),
        "NC_000001.11:g.1000_1001delinsAC",
    );
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
    // Two SUBs at the SAME position (zero gap, but overlap not adjacency)
    // must not collapse into a single delins. Emit OVERLAP_CONFLICTING_EDITS.
    let (text, codes) = normalize_to_string_and_warning_codes("NC_000001.11:g.[100G>A;100A>C]");
    assert!(text.contains("100G>A"), "got {}", text);
    assert!(text.contains("100A>C"), "got {}", text);
    assert!(!text.contains("delins"), "got {}", text);
    assert!(
        codes.contains(&"OVERLAP_CONFLICTING_EDITS".to_string()),
        "expected OVERLAP_CONFLICTING_EDITS, got {:?}",
        codes
    );
}

#[test]
fn test_merge_chain_then_barrier_then_chain() {
    // Three consecutive subs, then a duplication (a non-mergeable barrier),
    // then two more consecutive subs. The expected output is two separate
    // merged delins surrounding the dup. Pins the chain-flush-on-barrier
    // and chain-restart-after-barrier paths used by deferred-materialization.
    let result = normalize_to_string("NC_000001.11:g.[100G>A;101A>C;102C>T;103dup;104A>G;105G>T]");
    assert!(result.contains("100_102delinsACT"), "got {}", result);
    assert!(result.contains("103dup"), "got {}", result);
    assert!(result.contains("104_105delinsGT"), "got {}", result);
    // No spurious extra merges.
    assert_eq!(
        result.matches("delins").count(),
        2,
        "expected exactly two delins in {}",
        result
    );
}

#[test]
fn test_merge_long_consecutive_chain_correctness() {
    // 10-element chain of consecutive substitutions. The merged alt must be
    // the concatenation of all alt bases in input order. Pins correctness
    // for chains long enough that the deferred-materialization path
    // dominates the work.
    let result =
        normalize_to_string("NC_000001.11:g.[1000G>A;1001A>C;1002C>T;1003T>G;1004G>A;1005A>C;1006C>T;1007T>G;1008G>A;1009A>C]");
    assert_eq!(result, "NC_000001.11:g.1000_1009delinsACTGACTGAC");
}

#[test]
fn test_merge_chain_at_input_end() {
    // Chain spans the last variants in the allele — the loop ends mid-chain
    // and the deferred-materialization version must flush at loop end, not
    // leave the last chain unemitted.
    let result = normalize_to_string("NC_000001.11:g.[100dup;101G>A;102A>C]");
    assert!(result.contains("100dup"), "got {}", result);
    assert!(result.contains("101_102delinsAC"), "got {}", result);
}

#[test]
fn test_unmerged_pending_passes_through_unchanged() {
    // A merge-eligible variant followed by a non-mergeable barrier with no
    // mergeable predecessor: the eligible variant must emit AS-IS (a
    // substitution, not a delins). Pins the was_merged=false flush path.
    let result = normalize_to_string("NC_000001.11:g.[100G>A;200dup]");
    assert!(result.contains("100G>A"), "got {}", result);
    assert!(result.contains("200dup"), "got {}", result);
    // CRITICAL: the lone sub must NOT have been turned into a delins by
    // building from its anchor — that's the pitfall the was_merged flag
    // exists to prevent.
    assert!(!result.contains("delins"), "got {}", result);
}

// =====================================================================
// SPDI/VCF interaction with merge.
//
// Audit findings (recorded for future readers, not just tests):
//   * Neither hgvs_to_spdi_simple nor any vcf::to_hgvs path *produces*
//     HgvsVariant::Allele, so the merge pass cannot affect VCF/SPDI -> HGVS.
//   * hgvs_to_spdi_simple accepts only HgvsVariant::Genome. Merge collapses
//     a cis allele into a Genome (or Cds/etc.), changing which branch is
//     taken downstream. For a sub-sub merge the result is a Delins, which
//     SPDI cannot encode without the deleted reference sequence — that
//     limitation is independent of merge and the test below pins it.
//   * vcf::from_hgvs::HgvsToVcfConverter erred on multi-variant Cis
//     alleles. Post-merge, those alleles become a single variant and now
//     follow the single-variant convert path (success there is gated on
//     reference availability for delins, which MockProvider lacks).
// =====================================================================

#[test]
fn test_merged_allele_changes_variant_kind() {
    // Before merge: HgvsVariant::Allele wrapping two Genome subs.
    // After merge: a single Genome variant (no Allele wrapper) — so
    // downstream code dispatching on variant kind hits the Genome branch.
    use ferro_hgvs::HgvsVariant;
    let normalizer = Normalizer::new(MockProvider::new());
    let parsed = parse_hgvs("NC_000001.11:g.[1000G>A;1001A>C]").expect("parse failed");
    assert!(matches!(parsed, HgvsVariant::Allele(_)));
    let normalized = normalizer.normalize(&parsed).expect("normalize failed");
    assert!(
        matches!(normalized, HgvsVariant::Genome(_)),
        "expected Genome after merge, got {:?}",
        normalized
    );
}

#[test]
fn test_singleton_cis_allele_preserves_wrapper() {
    // The unwrap in normalize_allele must only fire when a merge actually
    // collapsed multiple sub-variants — not for pre-existing singletons.
    // Parsing `[1000G>A]` yields an AlleleVariant with one sub-variant; with
    // no merge happening, the wrapper must round-trip so programmatic callers
    // still see HgvsVariant::Allele. (Display already renders the bracket-free
    // form via AlleleVariant::Display, so user-visible output is unchanged.)
    use ferro_hgvs::HgvsVariant;
    let normalizer = Normalizer::new(MockProvider::new());
    let parsed = parse_hgvs("NC_000001.11:g.[1000G>A]").expect("parse failed");
    assert!(matches!(parsed, HgvsVariant::Allele(_)));
    let normalized = normalizer.normalize(&parsed).expect("normalize failed");
    assert!(
        matches!(normalized, HgvsVariant::Allele(_)),
        "expected Allele preserved (no merge happened), got {:?}",
        normalized
    );
}

#[test]
fn test_merged_delins_spdi_needs_reference() {
    // Documents that hgvs_to_spdi_simple cannot encode a merged delins
    // without reference data: the merged form drops per-base ref info from
    // the input subs, and SPDI requires the deleted bases. Pre-merge this
    // call would have failed with UnsupportedVariantType (Allele); post-merge
    // it fails with MissingReferenceData (Delins). The capability hasn't
    // changed — just the failure mode.
    let normalizer = Normalizer::new(MockProvider::new());
    let parsed = parse_hgvs("NC_000001.11:g.[1000G>A;1001A>C]").expect("parse failed");
    let normalized = normalizer.normalize(&parsed).expect("normalize failed");
    let err = hgvs_to_spdi_simple(&normalized).expect_err("delins requires ref data");
    let msg = format!("{:?}", err);
    assert!(
        msg.contains("MissingReferenceData") || msg.contains("reference"),
        "expected ref-data failure, got {}",
        msg
    );
}

#[test]
fn test_unmerged_allele_still_unsupported_by_spdi() {
    // Sanity check the non-merge path: hgvs_to_spdi_simple still refuses
    // an Allele that didn't collapse (variants stay as separate sub-variants).
    let normalizer = Normalizer::new(MockProvider::new());
    // One nt gap -> no merge, still an Allele after normalize.
    let parsed = parse_hgvs("NC_000001.11:g.[100G>A;102C>T]").expect("parse failed");
    let normalized = normalizer.normalize(&parsed).expect("normalize failed");
    let result = hgvs_to_spdi_simple(&normalized);
    assert!(
        result.is_err(),
        "unmerged Allele should not convert to SPDI, got {:?}",
        result
    );
}

// =====================================================================
// Merge runs even when overlap prevention is disabled.
//
// The merge pass uses a strict `a.end + 1 == b.start` adjacency check,
// so an overlapping pair never merges regardless of the prevent_overlap
// flag — but adjacent pairs still must merge.
// =====================================================================

fn normalize_to_string_no_overlap_check(input: &str) -> String {
    let config = NormalizeConfig::default().with_overlap_prevention(false);
    let normalizer = Normalizer::with_config(MockProvider::new(), config);
    let variant = parse_hgvs(input).expect("parse failed");
    let normalized = normalizer.normalize(&variant).expect("normalize failed");
    format!("{}", normalized)
}

#[test]
fn test_merge_still_runs_with_prevent_overlap_disabled() {
    // Adjacent pair must still collapse to a delins; merge is independent
    // of the overlap-prevention pre-pass.
    assert_eq!(
        normalize_to_string_no_overlap_check("NC_000001.11:g.[1000G>A;1001A>C]"),
        "NC_000001.11:g.1000_1001delinsAC",
    );
}

#[test]
fn test_overlap_pair_does_not_merge_with_prevent_overlap_disabled() {
    // Same-position pair is an overlap, not adjacency. Even with the
    // overlap-prevention pre-pass disabled, the strict adjacency check
    // in merge_consecutive_edits refuses to combine them.
    let result = normalize_to_string_no_overlap_check("NC_000001.11:g.[100G>A;100A>C]");
    assert!(result.contains("100G>A"), "got {}", result);
    assert!(result.contains("100A>C"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

// =====================================================================
// Codon-frame exception (issue #79).
//
// HGVS spec carves out a single exception to the strictly-consecutive
// rule: two `c.` exonic SNVs separated by exactly one nucleotide that
// fall within the same codon must merge into a single delins, with
// the unchanged middle nucleotide preserved verbatim. Spec example:
//
//   `c.[145C>T;147C>G]`  →  `c.145_147delinsTGG`
//
// Applies only to the CDS proper. UTR / `n.` / `r.` / `g.` have no
// codon frame and continue to round-trip unchanged for gap-of-1 inputs.
// =====================================================================

#[test]
fn test_codon_frame_two_subs_one_codon() {
    // Tx sequence: ATGCAAAAACCCCCGGGGG... (cds_start = 1).
    //   c.10 = C, c.11 = C, c.12 = C.
    // Codon for c.10 is `(10-1)/3 = 3`; codon for c.12 is `(12-1)/3 = 3`.
    // Same codon. The unchanged middle base is the actual reference at
    // c.11 = C, so the merged delins alt is
    // `G(c.10) + C(ref c.11) + A(c.12) = "GCA"`.
    //
    // Inputs use ref-matching SNVs (`c.10C>G;c.12C>A`) so the codon-frame
    // merge result is a real 3-base delins; under issue #160 the merged
    // delins is re-canonicalized post-merge, which would otherwise trim
    // the codon-frame-merged form down to a single sub if the second
    // SNV's stated alt happened to equal the actual middle reference.
    assert_eq!(
        normalize_with_provider(
            provider_with_simple_transcript(),
            "NM_TEST.1:c.[10C>G;12C>A]",
        ),
        "NM_TEST.1:c.10_12delinsGCA",
    );
}

#[test]
fn test_no_codon_frame_pair_straddles_codon_boundary() {
    // c.3 sits in codon 0 (positions 1-3); c.5 sits in codon 1
    // (positions 4-6). Different codons → no merge.
    let result =
        normalize_with_provider(provider_with_simple_transcript(), "NM_TEST.1:c.[3G>T;5A>C]");
    assert!(result.contains("3G>T"), "got {}", result);
    assert!(result.contains("5A>C"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_no_codon_frame_gap_of_two() {
    // Gap of 2 nt is beyond the spec carve-out, which is exactly one.
    let result = normalize_with_provider(
        provider_with_simple_transcript(),
        "NM_TEST.1:c.[10A>G;13A>C]",
    );
    assert!(result.contains("10A>G"), "got {}", result);
    assert!(result.contains("13A>C"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_no_codon_frame_genomic() {
    // `g.` has no codon frame — gap-of-1 must round-trip unchanged.
    let result = normalize_to_string("NC_000001.11:g.[145G>A;147C>T]");
    assert!(result.contains("145G>A"), "got {}", result);
    assert!(result.contains("147C>T"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

#[test]
fn test_codon_frame_then_strict_chain() {
    // The merge stage from issue #79 still extends a codon-frame pair
    // with a strictly-adjacent third SNV — c.10A>G + c.12A>C codon-frame
    // merge into `c.10_12delinsGCC`, then c.13A>T extends to
    // `c.10_13delinsGCCT`. After A10 (issue #165) the post-merge
    // `decompose_delins` pass runs over the 4-base result.
    //
    // Sequence: ATGCAAAAACCCCC... → c.10..c.13 = `CCCC` (the user's
    // declared "A" refs are ignored by the merge and recovered from the
    // provider during decomposition). Pair-by-pair canonicalization:
    //   pos 10: ref C, alt G (sub)
    //   pos 11: ref C, alt C (identity — codon-frame synthesized middle)
    //   pos 12: ref C, alt C (identity — user-declared `12A>C` is a
    //                         no-op against the actual reference)
    //   pos 13: ref C, alt T (sub)
    //
    // Per `general.md:34` the canonical form is two separate subs at
    // c.10 and c.13 — they are >1 unchanged nucleotide apart, so neither
    // adjacency (`substitution.md`) nor the codon-frame exception
    // (`general.md:35-38`) applies. The intermediate 4-base delins
    // produced by #79's chain extension was a ferro-internal stop on the
    // way to this canonical form; A10 collapses it.
    //
    // Spec basis for the chain undoing itself: the codon-frame exception
    // is defined for a *pair* of variants ("two variants … together
    // affecting one amino acid"), not for a chain. c.10 and c.13 sit in
    // codons 4 and 5 respectively, so the chain spans two amino acids
    // and the exception does not extend to it.
    assert_eq!(
        normalize_with_provider(
            provider_with_simple_transcript(),
            "NM_TEST.1:c.[10A>G;12A>C;13A>T]",
        ),
        "NM_TEST.1:c.[10C>G;13C>T]",
    );
}

#[test]
fn test_no_codon_frame_without_provider_transcript() {
    // Without a transcript registered on the provider, the middle-base
    // reference lookup fails and the merge passes through unchanged
    // (no error). Demonstrates the graceful-fallback contract.
    let result = normalize_to_string("NM_NOSEQ.99:c.[10A>G;12A>C]");
    assert!(result.contains("10A>G"), "got {}", result);
    assert!(result.contains("12A>C"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

// =====================================================================
// Same-region UTR adjacency (issue #89).
//
// PR #80 (issue #72) implements consecutive-edit merging in cis alleles
// but rejects every UTR position outright. Per the design doc and the
// HGVS spec, same-region adjacency *within* a UTR is valid and should
// merge — only crossings (5'UTR↔CDS, CDS↔3'UTR) remain barriers, since
// no HGVS range syntax spans those (e.g., `c.-1_1` does not exist).
// The corresponding crossing tests above (`test_no_merge_utr_boundary`)
// must continue to pass after this change.
// =====================================================================

/// 60-base transcript with 5'UTR (10 bases) + CDS (40 bases) + 3'UTR (10
/// bases). Sequence chosen so that the substitution ref-bases used in
/// the tests below match: `c.-2 = A`, `c.-1 = C`, `c.*1 = A`, `c.*2 = C`.
fn provider_with_utr_transcript() -> MockProvider {
    use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
    let mut provider = MockProvider::new();
    // 5'UTR (positions 1-10):     A A A A A A A A A C   → c.-10..c.-1
    // CDS    (positions 11-50):   ATG + 12 × ATGC      → c.1..c.40
    // 3'UTR  (positions 51-60):   A C T T T T T T T T   → c.*1..c.*10
    let sequence: String =
        "AAAAAAAAACATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCACTTTTTTTT".to_string();
    assert_eq!(sequence.len(), 60, "fixture length must be 60");
    let len = sequence.len() as u64;
    let exons = vec![Exon::new(1, 1, len)];
    let transcript = Transcript::new(
        "NM_TESTUTR.1".to_string(),
        Some("UTR_TEST".to_string()),
        Strand::Plus,
        sequence,
        Some(11),
        Some(50),
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

#[test]
fn test_merge_5utr_consecutive_subs_cds() {
    // `c.-2A>G;c.-1C>T` are physically adjacent within the 5'UTR.
    // No HGVS range crosses the 5'UTR/CDS boundary, but a range *within*
    // the 5'UTR is valid (`c.-2_-1`); the merged form is delinsGT, which
    // is the reverse complement of the deleted AC → canonicalized to inv.
    assert_eq!(
        normalize_with_provider(
            provider_with_utr_transcript(),
            "NM_TESTUTR.1:c.[-2A>G;-1C>T]",
        ),
        "NM_TESTUTR.1:c.-2_-1inv",
    );
}

#[test]
fn test_merge_3utr_consecutive_subs_cds() {
    // `c.*1A>G;c.*2C>T` are physically adjacent within the 3'UTR. Merged
    // form `delinsGT` over deleted `AC` canonicalizes to inv (rev-comp).
    assert_eq!(
        normalize_with_provider(
            provider_with_utr_transcript(),
            "NM_TESTUTR.1:c.[*1A>G;*2C>T]",
        ),
        "NM_TESTUTR.1:c.*1_*2inv",
    );
}

#[test]
fn test_merge_5utr_consecutive_subs_rna() {
    // r. mirrors c. for 5'UTR. The merged edit deletes the 5'UTR reference
    // `ac` and inserts `gu`, which is its reverse complement — so the
    // canonical form is an inversion, exactly as the DNA-axis analog
    // (`c.-2_-1delinsGT` -> `c.-2_-1inv`). Before #736 the RNA path kept this
    // as `delinsgu` because the inserted `u` (`Base::U`) was compared against
    // the DNA reference `T` and never matched the reverse complement; with the
    // U->T normalization the inversion is now recognized on r. too.
    assert_eq!(
        normalize_with_provider(
            provider_with_utr_transcript(),
            "NM_TESTUTR.1:r.[-2a>g;-1c>u]",
        ),
        "NM_TESTUTR.1:r.-2_-1inv",
    );
}

#[test]
fn test_merge_3utr_consecutive_subs_rna() {
    // 3'UTR analog of the above: ref `ac` -> `gu` (its reverse complement) is
    // canonicalized to `inv`, matching the DNA axes (#736).
    assert_eq!(
        normalize_with_provider(
            provider_with_utr_transcript(),
            "NM_TESTUTR.1:r.[*1a>g;*2c>u]",
        ),
        "NM_TESTUTR.1:r.*1_*2inv",
    );
}

#[test]
fn test_no_merge_3utr_cds_boundary() {
    // The c.40 (last CDS base) ↔ c.*1 (first 3'UTR base) crossing must
    // NOT merge — there is no HGVS range syntax spanning the CDS/3'UTR
    // boundary (`c.40_*1` does not exist), even though the positions are
    // physically adjacent in the transcript. Companion to the existing
    // `test_no_merge_utr_boundary` (5'UTR↔CDS).
    let result = normalize_with_provider(
        provider_with_utr_transcript(),
        "NM_TESTUTR.1:c.[40C>T;*1A>G]",
    );
    assert!(result.contains("40C>T"), "got {}", result);
    assert!(result.contains("*1A>G"), "got {}", result);
    assert!(!result.contains("delins"), "got {}", result);
}

// =====================================================================
// Scope-C: Coincident-bounds overlap warning tests.
//
// These tests verify that edits with identical boundaries (or same
// single position) emit OVERLAP_CONFLICTING_EDITS when they are
// conflicting cis edits on the same allele.
// =====================================================================

#[test]
fn test_overlap_sub_then_del_at_same_base() {
    let (text, codes) = normalize_to_string_and_warning_codes("NC_000001.11:g.[100A>C;100del]");
    assert!(text.contains("100A>C"), "got {}", text);
    assert!(text.contains("100del"), "got {}", text);
    assert!(
        codes.contains(&"OVERLAP_CONFLICTING_EDITS".to_string()),
        "got {:?}",
        codes
    );
}

#[test]
fn test_overlap_del_then_sub_at_same_base() {
    let (text, codes) = normalize_to_string_and_warning_codes("NC_000001.11:g.[100del;100A>C]");
    assert!(text.contains("100del"), "got {}", text);
    assert!(text.contains("100A>C"), "got {}", text);
    assert!(
        codes.contains(&"OVERLAP_CONFLICTING_EDITS".to_string()),
        "got {:?}",
        codes
    );
}

#[test]
fn test_overlap_two_dels_at_same_base() {
    let (text, codes) = normalize_to_string_and_warning_codes("NC_000001.11:g.[100del;100del]");
    assert_eq!(
        text.matches("100del").count(),
        2,
        "expected both dels preserved, got {}",
        text
    );
    assert!(
        codes.contains(&"OVERLAP_CONFLICTING_EDITS".to_string()),
        "got {:?}",
        codes
    );
}

#[test]
fn test_overlap_coincident_range_del_inv() {
    let (text, codes) = normalize_with_provider_and_warning_codes(
        provider_with_simple_transcript(),
        "NM_TEST.1:c.[100_103del;100_103inv]",
    );
    assert!(text.contains("100_103del"), "got {}", text);
    assert!(text.contains("100_103inv"), "got {}", text);
    assert!(
        codes.contains(&"OVERLAP_CONFLICTING_EDITS".to_string()),
        "got {:?}",
        codes
    );
}

#[test]
fn test_overlap_coincident_range_two_dels() {
    let (text, codes) =
        normalize_to_string_and_warning_codes("NC_000001.11:g.[100_103del;100_103del]");
    assert_eq!(text.matches("100_103del").count(), 2, "got {}", text);
    assert!(
        codes.contains(&"OVERLAP_CONFLICTING_EDITS".to_string()),
        "got {:?}",
        codes
    );
}

#[test]
fn test_no_overlap_warning_for_trans_alleles() {
    let (_text, codes) =
        normalize_to_string_and_warning_codes("[NC_000001.11:g.100A>C];[NC_000001.11:g.100A>G]");
    assert!(
        !codes.contains(&"OVERLAP_CONFLICTING_EDITS".to_string()),
        "trans alleles must not emit OVERLAP_CONFLICTING_EDITS, got {:?}",
        codes
    );
}

#[test]
fn test_no_overlap_warning_for_multi_accession_cis() {
    let (_text, codes) =
        normalize_to_string_and_warning_codes("[NC_000001.11:g.100A>C;NC_000002.11:g.100A>G]");
    assert!(
        !codes.contains(&"OVERLAP_CONFLICTING_EDITS".to_string()),
        "multi-accession cis must not emit OVERLAP_CONFLICTING_EDITS, got {:?}",
        codes
    );
}

#[test]
fn test_no_overlap_warning_for_adjacency() {
    let (_text, codes) = normalize_to_string_and_warning_codes("NC_000001.11:g.[100A>C;101A>G]");
    assert!(
        !codes.contains(&"OVERLAP_CONFLICTING_EDITS".to_string()),
        "adjacency must not emit OVERLAP_CONFLICTING_EDITS, got {:?}",
        codes
    );
}

// =====================================================================
// Codon-frame exception on a LENGTH-CHANGING block (#1235).
//
// `general.md:35`'s exception is not only about two substitutions. Once
// #1235 widened the sequence-first axis gate to `c.`/`n.`/`r.`, a lone
// coding `delins` whose payload coincidentally aligns against one
// interior base started splitting into an allele — the shape
// `delins.md:44-47` keeps whole. `apply_coding_codon_exception` cannot
// reach it: that pass matches `[Sub@p; Identity@p+1; Sub@p+2]` exactly,
// and here neither member is a substitution.
//
// The axis is the whole discriminator, so both halves are asserted on
// the same block and the same reference. `NM_TEST.1` has `cds_start = 1`,
// so `c.N` and `n.N` address the identical base.
// =====================================================================

/// `c.2_5` is `TGCA` and the payload `AAC` aligns only at its last column
/// (`C`), leaving a one-base gap between the two derived runs — so the
/// *distance* conjunct of `general.md:35` is met. The second conjunct is not:
/// `c.2` is in codon 1 and `c.5` in codon 2, and no codon holds four
/// consecutive positions, so the pair cannot together affect one amino acid at
/// any phase. `:35` is an exception to `:34`, not a replacement for it, so `:34`
/// governs and the block is described individually — the same answer the `n.`
/// axis reaches, because on this block the reading frame changes nothing.
///
/// # This expectation was flipped, and the old rationale is the reason
///
/// It previously read `c.2_5delinsAAC`, on the argument that a merge wider than
/// the spec's substitution triplet "is governed by the distance rule alone — the
/// reading the SVD-WG states it is moving to outright (`general.md:39`)". That
/// citation is dead text. `general.md:36-39` is a future-tense description of
/// proposal **SVD-WG010**, which `docs/consultation/SVD-WG010.md:5` records as
/// **rejected**, and the decided ruling `delins-codon-carve-out-gap-one` says so
/// in terms: "Do not cite it, and do not cite SVD-WG010 as a live rule either.
/// The codon-conditioned exception at `:18` / `general.md:35` remains the law."
/// That ruling also states the scope this file's expectation had dropped — "the
/// carve-out is scoped to the coding sequence **and to one amino acid**".
///
/// So the flip is not a new policy. It is this test catching up with a
/// correction already recorded against the clause it cites.
///
/// The `c.` assertion is now `input != output`, so unlike the old expectation it
/// cannot be satisfied by an unrelated decline — a window clamp, the weight
/// bound, the round-trip guard — and the `n.` control below is what pins that
/// the derivation reaches the block at all.
///
/// The discriminating half stays in
/// `merge::tests::coalesce_coding_frame_separation_declines_a_span_wider_than_a_codon`,
/// which asserts against the pass directly and at every codon phase.
///
/// # RE-PINNED BY #1610, THEN RE-PINNED AGAIN BY #2155 — AND THE TWO AXES NOW
/// AGREE ONCE MORE
///
/// This pair used to rest on the axes agreeing: `NM_TEST.1` has `cds_start = 1`,
/// so `c.p` and `n.p` name the same base, and nothing distinguished them. #1610
/// broke that agreement — see below — and #2155 restores it, for a different
/// reason than the original one.
///
/// The block is `TGCA -> AAC`, four reference bases against three, whose only
/// interior match is the payload's own `C`.
/// `merge::split_is_a_placed_gap_coincidence` (#1610) keeps such a block whole
/// (`DNA/delins.md:46`/`:47`) — gated on
/// `CoincidenceCarveOut::may_disbelieve_a_separation`, per
/// `rulings[delins-payload-coincidence-carve-out-is-coding-dna-scoped]`. #1610
/// shipped that gate scoped to the coding DNA axis only, so `c.` returned the
/// input whole while `n.` still split. **#2155 supersedes that scope to every
/// DNA axis** (`c./g./m./n.`, `r.` excluded), so `n.` is now in reach too and
/// both axes return the block whole.
///
/// The codon conjunct is no longer discriminated here on either axis, because
/// there are no two derived members left to ask about on either. Its guard is
/// `merge::tests::coalesce_coding_frame_separation_declines_a_span_wider_than_a_codon`,
/// which builds its two-member fixture by hand for exactly this reason, and
/// `a_length_changing_block_the_exception_still_cannot_reach` below, which is an
/// end-to-end case that still splits on both axes (its split members carry a
/// lone substitution, condition 4 of `split_is_a_placed_gap_coincidence`, which
/// the axis widen does not touch).
#[test]
fn coding_length_changing_block_across_codons_is_described_individually() {
    assert_eq!(
        normalize_with_provider(
            provider_with_simple_transcript(),
            "NM_TEST.1:n.2_5delinsAAC"
        ),
        "NM_TEST.1:n.2_5delinsAAC",
        "#2155 widens `delins.md:47`'s carve-out to every DNA axis, so `n.` is \
         now in reach and the block merges just as `c.` does",
    );
    assert_eq!(
        normalize_with_provider(
            provider_with_simple_transcript(),
            "NM_TEST.1:c.2_5delinsAAC"
        ),
        "NM_TEST.1:c.2_5delinsAAC",
        "#1610: on the coding DNA axis the block is kept whole by the \
         partitioner, so no member pair reaches the codon exception. Same \
         reference, same block, same coordinates as the `n.` assertion above — \
         and, as of #2155, the same output too",
    );
}

/// A length-changing block that **does** still split across a codon boundary,
/// end to end.
///
/// Added with #1610, which took the discriminating power out of the pair above.
/// `c.2_4` is `TGC` and the payload `CG` is two bases, so the block is a net
/// deletion like the one above — but its derived split carries a lone
/// **substitution**, which is a rank-1 type the split genuinely buys, so
/// `split_is_a_placed_gap_coincidence` declines and the two members survive to
/// be asked the codon question. `c.2` is in codon 1 and `c.4` in codon 2, so
/// `general.md:35`'s second conjunct is unmet and `:34` governs.
///
/// Both axes are asserted, on one block and one reference, so the reading frame
/// is the only thing that varies — the property the pair above used to carry.
#[test]
fn a_length_changing_block_the_exception_still_cannot_reach() {
    assert_eq!(
        normalize_with_provider(provider_with_simple_transcript(), "NM_TEST.1:n.2_4delinsCG"),
        "NM_TEST.1:n.[2T>C;4del]",
        "an axis with no reading frame gets `general.md:34`'s plain rule",
    );
    assert_eq!(
        normalize_with_provider(provider_with_simple_transcript(), "NM_TEST.1:c.2_4delinsCG"),
        "NM_TEST.1:c.[2T>C;4del]",
        "`c.2` and `c.4` sit in different codons, so `general.md:35`'s \
         one-amino-acid conjunct is unmet and the members stay individual",
    );
}

/// The same block on `n.`, asserted on its own so the axis half of the rule
/// keeps a test of its own rather than living only as the control above. The
/// deletion 3'-shifts down the `A` run at `n.5-9`, which is why the coding
/// exception has to be evaluated *before* the shift: afterwards the two runs
/// are four bases apart and no one-base rule could rejoin them even where the
/// codon test allows it.
///
/// # UNMOVED BY #1610, RE-PINNED BY #2155
///
/// An earlier revision of this branch (pre-#1610) pinned this to the whole
/// `n.2_5delinsAAC`, then #1610 shipped `split_is_a_placed_gap_coincidence`
/// scoped to the coding DNA axis only, moving `n.` to the split pinned above
/// for a time. **#2155 supersedes that scope to every DNA axis**
/// (`rulings[delins-payload-coincidence-carve-out-is-coding-dna-scoped]`), so
/// `n.` is back in reach and this row returns to the whole form — the same
/// string the pre-#1610 revision carried, though for the current, wider-scope
/// reason rather than the original axis-blind one.
///
/// This is the axis **#1610 was filed on**, so as of #2155 that issue's own
/// reproduction is closed on this axis too. See the module docs of
/// `issue_1610_lone_unequal_length_delins`.
///
/// NOTE: the function name is kept for continuity with the pre-#2155 history
/// above; despite reading "keeps_the_split", this test now asserts the block
/// **merges** on `n.` (see the assertion below). The name is left stable so
/// cross-references do not churn — the doc comment is the source of truth.
#[test]
fn non_coding_length_changing_block_keeps_the_split() {
    let result = normalize_with_provider(
        provider_with_simple_transcript(),
        "NM_TEST.1:n.2_5delinsAAC",
    );
    assert_eq!(
        result, "NM_TEST.1:n.2_5delinsAAC",
        "#2155 widens `delins.md:47`'s carve-out to every DNA axis, so `n.` \
         merges the block whole just as `c.` does",
    );
}

// =====================================================================
// The amino-acid half of `DNA/delins.md:18`, on a codon-designed core.
//
// `:17` says two variants separated by one or more nucleotides are
// described individually. `:18` excepts the pair that is separated by
// *one* nucleotide **and** together affects **one amino acid** — two
// conjuncts, not one. The tests below vary only the second, on one
// reference, one block shape and one pair of edit types, so the codon
// boundary is the whole of the difference between them.
//
// The precondition is stated identically on the length-changing pages —
// `DNA/deletion.md:19`, `DNA/insertion.md:20`, `DNA/duplication.md:23`,
// `DNA/inversion.md:21` — so a deletion paired with a substitution gets
// no relaxation for being length-changing. `DNA/delins.md:81` restates
// it as an exclusion: two variants separated by one or more nucleotides
// are described individually "(unless they together affect one amino
// acid)".
// =====================================================================

/// A 30 nt coding transcript carrying the same `ACG` triplet at two codon
/// phases, and nothing else that can move.
///
/// `CDS_START` is 1, so `c.p` is transcript position `p` and codon `k`
/// spans `c.3k-2 ..= c.3k`.
///
/// * `c.10_12` is `ACG` and is codon 4 exactly — an aligned triplet.
/// * `c.14_16` is `ACG` and straddles codon 5 (`c.13_15`) and codon 6
///   (`c.16_18`).
///
/// The flanks are `GCT` repeats so that neither `A` sits in a run: `c.9`
/// and `c.13` are `T`, `c.11` and `c.15` are `C`, so the deletion member
/// has nowhere to 3'-shift and the two cases cannot differ by shuffling.
fn provider_with_codon_designed_core() -> MockProvider {
    use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
    //                    c.1      c.10 c.13 c.17
    //                    |        |    |    |
    let sequence = "GCTGCTGCT".to_string() + "ACG" + "T" + "ACG" + "T" + "GCTGCTGCTGCTG";
    let len = sequence.len() as u64;
    assert_eq!(len, 30, "the codon arithmetic in this module assumes 30 nt");
    assert_eq!(
        &sequence[9..12],
        "ACG",
        "c.10_12 must be the aligned triplet"
    );
    assert_eq!(
        &sequence[13..16],
        "ACG",
        "c.14_16 must be the straddling triplet"
    );
    let transcript = Transcript::new(
        "NM_CODON.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        sequence,
        Some(1),
        Some(len),
        vec![Exon::new(1, 1, len)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    let mut provider = MockProvider::new();
    provider.add_transcript(transcript);
    provider
}

/// **The exception applies.** `c.10del` and `c.12G>T` are separated by the
/// unchanged `c.11`, and `c.10`, `c.11`, `c.12` are codon 4 — so they
/// together affect one amino acid and `delins.md:18` asks for a `delins`.
///
/// The pair is length-changing (a deletion and a substitution, net -1),
/// which is exactly the point: the clause names no edit type, and the
/// length-changing pages state the same precondition without weakening it.
#[test]
fn one_base_gap_within_one_codon_merges() {
    assert_eq!(
        normalize_with_provider(
            provider_with_codon_designed_core(),
            "NM_CODON.1:c.[10del;12G>T]"
        ),
        "NM_CODON.1:c.10_12delinsCT",
    );
}

/// **The exception does not apply.** `c.14del` and `c.16G>T` are separated
/// by the unchanged `c.15` — the *distance* conjunct is met, identically to
/// the test above — but `c.14` is in codon 5 and `c.16` in codon 6, so the
/// pair does not together affect one amino acid. `:18` is an exception to
/// `:17`, not a replacement for it, so `:17` governs and the two variants
/// are described individually.
///
/// The only difference from the merging case is the codon phase: same
/// reference triplet, same two edit types, same one-base separation, same
/// axis. That is what makes this pair a test of the precondition rather
/// than of the shape.
///
/// # What this assertion does not establish
///
/// It is `input == output`, which — as the doc on
/// `coding_length_changing_block_across_codons_is_described_individually`
/// says of the expectation it replaced — *any* upstream decline satisfies as
/// readily as the codon test: a window clamp, the weight bound, the
/// round-trip guard. The `n.` control below is `input == output` at this
/// position too, so the pair cannot discriminate at `c.14_16` on its own.
/// What it does discriminate is the *reading-frame* half at the other
/// position, where `one_base_gap_within_one_codon_merges` yields a string
/// that differs from its input while the `n.` row does not.
///
/// The codon-phase half is pinned directly instead, by
/// `merge::tests::coalesce_coding_frame_separation_merges_only_within_one_codon`,
/// which calls the pass on one set of pieces at `w_lo` 1 and 2 — same
/// payloads, same flags, one integer apart — so only the phase can account
/// for the difference in its two answers.
#[test]
fn one_base_gap_across_a_codon_boundary_stays_split() {
    assert_eq!(
        normalize_with_provider(
            provider_with_codon_designed_core(),
            "NM_CODON.1:c.[14del;16G>T]"
        ),
        "NM_CODON.1:c.[14del;16G>T]",
        "`c.14` is in codon 5 and `c.16` in codon 6, so `delins.md:18`'s \
         second conjunct is unmet and `:17` governs — the two variants are \
         described individually",
    );
}

/// The control for the pair above: on `n.`, where there is no reading
/// frame, *neither* position merges. Without it the split at `c.14_16`
/// could be produced by any axis-blind decline — a window clamp, the
/// weight bound, the round-trip guard — rather than by the codon test.
///
/// `NM_CODON.1` has `cds_start = 1`, so `n.p` and `c.p` address the same
/// base and the four descriptions differ in nothing but their axis letter.
#[test]
fn without_a_reading_frame_neither_position_merges() {
    for position in [10, 14] {
        let input = format!("NM_CODON.1:n.[{position}del;{}G>T]", position + 2);
        assert_eq!(
            normalize_with_provider(provider_with_codon_designed_core(), &input),
            input,
            "`general.md:34`'s plain rule governs an axis with no amino acid to affect",
        );
    }
}
