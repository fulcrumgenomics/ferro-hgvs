//! Regression corpus for the tests that guard `msto`'s ferro-hgvs issues.
//!
//! The cis-allele normalizer is being rewritten to derive member boundaries
//! from the denoted sequence (`FERRO_SEQFIRST_SHADOW=1` reports where the new
//! splitter would differ from the live one). A prior review found that
//! driving live output from the new splitter changes the rendered string for
//! several tests tied to issues `msto` filed, without changing which
//! sequence is denoted. That review is not committed to this repository (it
//! is a working artifact); this module is the durable, checked artifact that
//! survives it.
//!
//! [`MSTO_ISSUES`] enumerates every issue `msto` has filed against
//! `fulcrumgenomics/ferro-hgvs` (`gh issue list --author msto --state all`)
//! and, for each, the test functions that reference it by number
//! (`issue_<N>`, `#<N>`, or `issue #<N>`) — in the function's own name, in a
//! comment immediately above it, or, for a module that is one issue family end
//! to end, in that module's doc. The two sections below record why the
//! module-doc form was admitted and where it is *not* claimed. Where no such
//! reference exists the `tests` slice is empty — that is itself the finding for
//! a `high`-labeled issue, per
//! [`no_high_priority_issue_is_unguarded_without_being_named`].
//!
//! ## Python guards count, and used not to
//!
//! The catalog is a report about what is *unguarded*, so a guard it cannot see
//! is a false negative in the one direction that matters. Until now it scanned
//! only Rust — `tests/` and `src/` — and the bindings' guards live in
//! `tests/python/`. Four issues were listed as unguarded while a dedicated
//! pytest module guarded each: **#1050, #1159, #1245, #1246**.
//!
//! Python attribution works differently and the rule above is adjusted for it
//! rather than bent. A Rust suite puts many unrelated tests in one file, so the
//! reference has to be per-function. The bindings' guards are instead written one
//! module per issue, and the issue number is named in the **module docstring** —
//! so for a Python entry the reference may be in that docstring rather than
//! beside each function, and [`source_defines_fn`] accepts `def` as well as `fn`.
//!
//! Three of the four carry the number in the filename too
//! (`test_issue_<N>_*.py`); #1159's module is named for its behaviour
//! (`test_reference_aware_spdi.py`) and names the issue in its first docstring
//! line. The docstring, not the filename, is what the attribution rests on — a
//! filename convention that three of four follow is not one to state as the rule.
//!
//! Only issues whose slice was **empty** were audited this way, because only
//! those produced a false "unguarded" claim. An issue that already lists a Rust
//! guard may still have an unlisted Python one; that understates coverage
//! without misreporting the issue as unguarded.
//!
//! ## Table-driven guards count, and used not to
//!
//! The same false negative, in a second shape. `tests/it/reported_confluence_pairs.rs`
//! and `tests/it/reported_partition_verdicts.rs` guard **#1419, #1420 and
//! #1421** — every row in both is one of those issues' reported spellings — but
//! the issue numbers live in the tables' *data rows* (`"1419-r1/cis"`,
//! `"1421-n3/span"`) and in the module docs, never in a function name or in the
//! comment above one. So the per-function rule could not see them, and all three
//! were listed as unguarded while eighteen pinned expectations and four
//! whole-corpus properties covered them.
//!
//! The accommodation is the Python one generalized rather than a new exception:
//! **where a module is one issue family end to end and says so in its module
//! doc, the module doc is the attribution.** That is sound for the same reason it
//! was sound for pytest — the unit of authorship is the module, not the function
//! — and it is *only* claimed for modules where it is true. A grab-bag suite
//! that happens to mention an issue in its header still needs the per-function
//! reference, because there the module doc says nothing about any individual
//! test.
//!
//! Note what this does not do: it does not lower the bar for **#200**, which
//! stays in [`no_high_priority_issue_is_unguarded_without_being_named`]. No
//! module anywhere is about #200, so there is nothing for a module-doc rule to
//! find. The rule widened to admit guards that exist; it did not widen until the
//! finding disappeared.
//!
//! ## An issue whose satisfaction is a property, not an output
//!
//! **#1430** is the one entry here that could never have a conventional
//! regression test, and its slice is populated on a different basis from every
//! other row. It is a design proposal — coerce a description to its denoted
//! sequence, derive a representative variant from that sequence, then optionally
//! apply the normalizer for spec compliance — and it states no input/expected
//! pair at all. There is no string to pin.
//!
//! What it can be measured by is whether each of its three stages holds as a
//! property, so its slice names one *or more* guards per stage rather than a
//! reproduction — currently two, three and two, plus one whole-proposal census.
//! See [`issue_1430_is_measured_by_a_property_not_by_an_expected_output`], which
//! pins that reading so it cannot quietly be replaced by an ordinary
//! input/expected row. What it requires is that no stage is left uncovered, not
//! that any stage carry a particular number of guards.
//!
//! ## A citation is not a guarantee that the test guards the issue
//!
//! The attribution rule is mechanical — `#<N>` in the function's name, in the
//! comment above it, or in the doc of a module that is one issue family end to
//! end — and that is deliberate, because a judgement-based rule rots faster
//! than the thing it describes. But mechanical means it also catches
//! *cross-references*: a test whose comment names an issue to explain why the
//! fixture was built a certain way, or as a historical aside, reads identically
//! to a test that reproduces it.
//!
//! An audit of every cataloged pair found two such rows, and they are recorded
//! here rather than removed, because both satisfy the rule as written and
//! deleting them would leave the table no longer a faithful application of its
//! own stated rule:
//!
//! * **#1034** lists
//!   `tests/it/issue_291_rna_axis_convention.rs::rna_canonical_split_fetch_is_cds_relative`.
//!   That test pins that the r. reference window is fetched CDS-relative
//!   (a #469 claim). #1034 appears only as the reason a *full-run* inversion was
//!   used to build the case rather than a sub-run.
//! * **#310** lists
//!   `src/reference/transcript.rs::test_transcript_default_is_minimal_and_non_coding`.
//!   That test pins `Transcript::default()`. #310 appears only as an anecdote
//!   about a past field addition.
//!
//! So read a populated `tests` slice as *"a test names this issue"*, which is
//! what it measures, and not as *"this issue is covered"*. The empty slices are
//! the load-bearing half and they are unaffected: a cross-reference can inflate
//! coverage, never hide a gap.
//!
//! ## What this module does NOT claim
//!
//! It does not re-run the shadow/live comparison, so it makes no claim about
//! which tests currently pass or fail under the rewrite. The `note` field on
//! the eleven pinned at-risk issues below summarizes the prior review's own
//! grouping (three mechanism clusters) as a risk flag for what to re-check
//! before the rewrite's default is flipped — it is not something a test in
//! this module executes or verifies.
//!
//! ## Self-checks
//!
//! A mapping nobody validates rots silently (a test gets renamed or deleted
//! and the table quietly stops guarding anything). So this module checks
//! itself:
//!
//! - [`every_cataloged_test_exists_in_its_source_file`] greps each
//!   `(file, fn)` pair's source file for a matching `fn` definition.
//! - [`issue_and_test_totals_match_their_pinned_counts`] pins the issue
//!   count, the `high`-label count, and the total number of distinct
//!   cataloged tests, so any addition/removal of an issue or a test
//!   reference is visible in a diff instead of buried in a 9,000-test
//!   summary.
//! - [`pinned_at_risk_issues_carry_a_mechanism_note`] checks that exactly
//!   the eleven issues named in the migration review carry a `note`, and
//!   that no other issue does.

use regex::Regex;
use std::collections::HashSet;
use std::path::Path;

/// Issue state as of the last `gh issue list` enumeration.
///
/// An enum rather than a `&'static str`, so a typo such as `"Closed"` is a
/// compile error instead of something a runtime assertion has to catch. The same
/// pattern is used for `Pinning` in
/// `src/normalize/merge/splitter_reproducer_corpus.rs`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum State {
    Open,
    Closed,
}

/// One `msto`-filed issue and the tests (if any) that reference it by number.
struct MstoIssue {
    /// GitHub issue number in `fulcrumgenomics/ferro-hgvs`.
    number: u32,
    /// Whether the issue carries the `high` label.
    high: bool,
    /// Issue state as of the last `gh issue list` enumeration.
    state: State,
    /// `(repo-relative file path, test function name)` pairs that reference
    /// this issue number in the function's own name, in a comment immediately
    /// above it, or — where the whole module is one issue family — in that
    /// module's doc. Empty means no such reference was found.
    tests: &'static [(&'static str, &'static str)],
    /// One-line mechanism note for the eleven issues the sequence-first
    /// migration review flagged as at risk. `None` for every other issue.
    note: Option<&'static str>,
}

/// Every issue `msto` has filed against `fulcrumgenomics/ferro-hgvs`
/// (`gh issue list -R fulcrumgenomics/ferro-hgvs --state all --author msto`),
/// mapped to the tests that reference it by number.
///
/// Assembled by searching `tests/` (Rust and Python) and `src/` for
/// `issue_<N>`, `#<N>`, and `issue #<N>` occurring in a test function's own
/// name, in a comment directly above it, or in the doc of a module that is one
/// issue family end to end. This is a *citation*-based search: a test that
/// exercises an issue's behavior without the number appearing in any of those
/// three places is not in this table.
///
/// The module-doc form is not a loophole — it is claimed only where the module
/// really is about one issue family, which is how the bindings' guards and the
/// two `reported_*` tables are written. See the module doc for both cases and
/// for the one it deliberately does not reach (#200).
#[rustfmt::skip]
static MSTO_ISSUES: &[MstoIssue] = &[
    MstoIssue { number: 30, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 31, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 33, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 35, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 36, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 40, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 45, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 46, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 47, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 51, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 54, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 59, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 60, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 61, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 62, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 67, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 69, high: false, state: State::Closed, tests: &[("src/hgvs/parser/accession.rs", "test_parse_simple_accession_with_gene_selector")], note: None },
    MstoIssue { number: 72, high: false, state: State::Closed, tests: &[("tests/it/allele_trans_phase.rs", "test_trans_does_not_trigger_consecutive_edit_merge")], note: None },
    MstoIssue { number: 73, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 74, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 75, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 79, high: false, state: State::Closed, tests: &[("tests/it/issue_165_delins_sub_only_decompose.rs", "cds_codon_frame_pair_preserved_as_delins"), ("tests/it/merge_consecutive_edits_tests.rs", "test_codon_frame_then_strict_chain"), ("tests/it/merge_consecutive_edits_tests.rs", "test_codon_frame_two_subs_one_codon")], note: None },
    MstoIssue { number: 82, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 83, high: false, state: State::Closed, tests: &[("tests/it/hgvs_spec_normalization_tests.rs", "harvester_yields_only_real_variants"), ("tests/it/hgvs_spec_normalization_tests.rs", "pinned_v21_normalization_behavior"), ("tests/it/protein_unknown_roundtrip.rs", "protein_unknown_trans_compound_with_position_unknown_arm_round_trips"), ("tests/it/protein_unknown_roundtrip.rs", "protein_unknown_trans_compound_with_predicted_arm_round_trips")], note: None },
    MstoIssue { number: 84, high: false, state: State::Closed, tests: &[("tests/it/hgvs_spec_normalization_tests.rs", "harvester_yields_only_real_variants"), ("tests/it/hgvs_spec_normalization_tests.rs", "pinned_v21_normalization_behavior")], note: None },
    MstoIssue { number: 87, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 160, high: false, state: State::Closed, tests: &[("tests/it/normalize_tests.rs", "cds_full_span_revcomp_in_cis_merges_to_inv"), ("tests/it/normalize_tests.rs", "cds_sub_span_revcomp_stays_single_delins"), ("tests/it/normalize_tests.rs", "cds_user_typed_delins_with_revcomp_subspan_stays_delins_symmetric"), ("tests/it/normalize_tests.rs", "full_span_revcomp_in_cis_merges_to_inv"), ("tests/it/normalize_tests.rs", "mt_user_typed_full_span_revcomp_delins_becomes_inv"), ("tests/it/normalize_tests.rs", "near_miss_complement_only_stays_delins"), ("tests/it/normalize_tests.rs", "near_miss_reverse_only_stays_delins"), ("tests/it/normalize_tests.rs", "rna_full_span_revcomp_in_cis_merges_to_inv"), ("tests/it/normalize_tests.rs", "rna_sub_span_revcomp_stays_single_delins"), ("tests/it/normalize_tests.rs", "single_nt_complement_remains_substitution"), ("tests/it/normalize_tests.rs", "split_drops_identity_subedit_no_eq_emitted"), ("tests/it/normalize_tests.rs", "sub_span_revcomp_stays_single_delins"), ("tests/it/normalize_tests.rs", "three_nt_revcomp_subrun_in_contiguous_change_stays_delins"), ("tests/it/normalize_tests.rs", "tx_user_typed_delins_with_full_span_revcomp"), ("tests/it/normalize_tests.rs", "user_typed_delins_with_revcomp_subspan_stays_delins_symmetric")], note: Some("Cluster A: the new splitter allows compensating gaps and shreds a contiguous replacement into pieces.") },
    MstoIssue { number: 161, high: false, state: State::Closed, tests: &[("tests/it/del_shift_matrix.rs", "issue_161_bracket_cis_allele_2bp_del_shifts_one_rotation"), ("tests/it/del_shift_matrix.rs", "issue_161_bracket_cis_allele_4bp_del_shifts_one_rotation"), ("tests/it/del_shift_matrix.rs", "issue_161_canonical_form_4bp_del_matches_bracket_form")], note: None },
    MstoIssue { number: 179, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 180, high: false, state: State::Closed, tests: &[("tests/it/issue_180_allele_3prime_shift.rs", "allele_ins_canonicalized_to_dup_then_shifted_3prime"), ("tests/it/issue_180_allele_3prime_shift.rs", "allele_norm_3prime_shift_is_idempotent"), ("tests/it/issue_180_allele_3prime_shift.rs", "allele_seven_single_dels_spanning_run_boundary_shift_3prime"), ("tests/it/issue_180_allele_3prime_shift.rs", "allele_single_del_in_bracket_shifts_to_3prime_end_of_run"), ("tests/it/issue_180_allele_3prime_shift.rs", "allele_two_single_dels_shift_to_3prime_end_of_homopolymer"), ("tests/it/issue_180_allele_3prime_shift.rs", "issue_181_example_1_adjacent_dels_no_duplicate_position"), ("tests/it/issue_180_allele_3prime_shift.rs", "issue_181_example_2_eight_adjacent_dels_no_overlap"), ("tests/it/issue_180_allele_3prime_shift.rs", "issue_181_example_3_mixed_dels_and_ins_no_duplicate_position")], note: None },
    MstoIssue { number: 181, high: false, state: State::Closed, tests: &[("tests/it/issue_180_allele_3prime_shift.rs", "issue_181_example_1_adjacent_dels_no_duplicate_position"), ("tests/it/issue_180_allele_3prime_shift.rs", "issue_181_example_2_eight_adjacent_dels_no_overlap"), ("tests/it/issue_180_allele_3prime_shift.rs", "issue_181_example_3_mixed_dels_and_ins_no_duplicate_position")], note: None },
    MstoIssue { number: 182, high: false, state: State::Closed, tests: &[("tests/it/issue_165_delins_sub_only_decompose.rs", "cds_adjacent_sub_pair_stays_as_delins"), ("tests/it/issue_182_postcanon_adjacency.rs", "identity_at_between_subs_prevents_grouping"), ("tests/it/issue_182_postcanon_adjacency.rs", "inv_flanked_by_two_sub_runs_each_side_groups_both_flanks"), ("tests/it/issue_182_postcanon_adjacency.rs", "single_sub_flank_stays_as_substitution"), ("tests/it/issue_182_postcanon_adjacency.rs", "two_subs_leading_an_inv_merge_to_delins"), ("tests/it/issue_182_postcanon_adjacency.rs", "two_subs_trailing_an_inv_merge_to_delins")], note: None },
    MstoIssue { number: 183, high: true, state: State::Closed, tests: &[("src/reference/annotation/builder.rs", "ladder_step3_cds_as_exon_for_issue_183"), ("src/reference/annotation/builder.rs", "ladder_step3_preserves_utr_when_cds_is_subset_of_mrna"), ("src/reference/annotation/mod.rs", "load_gff3_issue_183_single_exon_no_exon_line")], note: None },
    MstoIssue { number: 184, high: false, state: State::Closed, tests: &[("tests/it/cli_build_transcript.rs", "build_transcript_emit_genomic_sequences_makes_reference_genome_capable"), ("tests/it/cli_build_transcript.rs", "build_transcript_emit_genomic_sequences_stores_forward_bytes_on_minus_strand"), ("tests/it/cli_build_transcript.rs", "build_transcript_emits_valid_json_for_single_exon"), ("tests/it/cli_build_transcript.rs", "build_transcript_rejects_cds_beyond_contig"), ("tests/it/cli_build_transcript.rs", "build_transcript_supports_minus_strand_and_custom_id"), ("tests/it/cli_build_transcript.rs", "build_transcript_without_flag_omits_genomic_sequences")], note: None },
    MstoIssue { number: 185, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 200, high: true, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 310, high: true, state: State::Closed, tests: &[("src/project/protein/identity.rs", "whole_molecule_identity_anchors_on_the_reference_initiator"), ("src/reference/transcript.rs", "test_transcript_default_is_minimal_and_non_coding"), ("tests/it/gene_selector_roundtrip.rs", "display_drops_gene_selector_refseq_nm"), ("tests/it/issue_1099_all_silent_protein.rs", "gene_symbol_rendering_matches_the_changed_residue_sibling"), ("tests/it/issue_310_projector_non_refseq.rs", "non_refseq_id_with_explicit_protein_id_uses_protein_id"), ("tests/it/issue_310_projector_non_refseq.rs", "non_refseq_id_without_protein_id_falls_back_to_transcript_id"), ("tests/it/issue_310_projector_non_refseq.rs", "parser_still_accepts_legacy_np_with_selector"), ("tests/it/issue_310_projector_non_refseq.rs", "refseq_nm_prefix_without_protein_falls_back_to_transcript_id"), ("tests/it/issue_310_projector_non_refseq.rs", "refseq_xm_prefix_without_protein_falls_back_to_transcript_id")], note: None },
    MstoIssue { number: 311, high: true, state: State::Closed, tests: &[("tests/it/issue_311_transcript_prefix_chromosome.rs", "allele_normalizes_when_tx_id_differs_only_in_case"), ("tests/it/issue_311_transcript_prefix_chromosome.rs", "allele_normalizes_when_tx_id_does_not_collide"), ("tests/it/issue_311_transcript_prefix_chromosome.rs", "allele_normalizes_when_tx_id_equals_chromosome"), ("tests/it/issue_311_transcript_prefix_chromosome.rs", "allele_normalizes_when_tx_id_starts_with_chromosome"), ("tests/it/issue_311_transcript_prefix_chromosome.rs", "delins_roundtrips_when_tx_id_starts_with_chromosome"), ("tests/it/issue_311_transcript_prefix_chromosome.rs", "get_protein_sequence_versioning_contract"), ("tests/it/issue_311_transcript_prefix_chromosome.rs", "get_transcript_does_not_collide_with_chromosome_prefix"), ("tests/it/issue_311_transcript_prefix_chromosome.rs", "get_transcript_rejects_non_version_boundary_prefix"), ("tests/it/issue_311_transcript_prefix_chromosome.rs", "get_transcript_rejects_versioned_query_against_other_version"), ("tests/it/issue_311_transcript_prefix_chromosome.rs", "get_transcript_resolves_exact_id"), ("tests/it/issue_311_transcript_prefix_chromosome.rs", "get_transcript_resolves_unversioned_query_to_versioned_key")], note: None },
    MstoIssue { number: 1026, high: false, state: State::Closed, tests: &[("src/reference/mock.rs", "from_json_accepts_consistent_genomic_reference"), ("src/reference/mock.rs", "from_json_skips_backing_check_without_genomic_sequences"), ("tests/it/issue_1026_genome_capable_json.rs", "genome_capable_reference_normalizes_intronic"), ("tests/it/issue_1026_genome_capable_json.rs", "transcripts_only_reference_cannot_normalize_intronic")], note: None },
    MstoIssue { number: 1027, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 1028, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 1034, high: true, state: State::Closed, tests: &[("tests/it/issue_1034_inv_subspan_delins.rs", "contiguous_delins_with_internal_revcomp_subrun_stays_delins"), ("tests/it/issue_1034_inv_subspan_delins.rs", "contiguous_delins_with_leading_revcomp_subrun_stays_delins"), ("tests/it/issue_1034_inv_subspan_delins.rs", "full_run_reverse_complement_still_emits_inv"), ("tests/it/issue_1034_inv_subspan_delins.rs", "two_base_full_run_reverse_complement_still_emits_inv"), ("tests/it/issue_1040_inv_overrecognition_probes.rs", "issue_1034_minimal_contiguous_run_stays_delins"), ("tests/it/issue_291_rna_axis_convention.rs", "rna_canonical_split_fetch_is_cds_relative"), ("tests/it/normalize_tests.rs", "cds_full_span_revcomp_in_cis_merges_to_inv")], note: Some("Cluster A: the new splitter allows compensating gaps and shreds a contiguous replacement into pieces.") },
    MstoIssue { number: 1035, high: false, state: State::Closed, tests: &[("tests/it/build_transcript_library_parity.rs", "library_matches_cli_minus_strand_custom_id_gene"), ("tests/it/build_transcript_library_parity.rs", "library_matches_cli_plus_strand"), ("tests/it/build_transcript_library_parity.rs", "library_matches_cli_with_emit_genomic_sequences"), ("tests/it/convert_gff_library_parity.rs", "library_matches_cli_transcript_only"), ("tests/it/convert_gff_library_parity.rs", "library_matches_cli_with_emit_genomic_sequences"), ("tests/it/convert_gff_library_parity.rs", "library_matches_cli_with_gene_filter")], note: None },
    MstoIssue { number: 1040, high: true, state: State::Closed, tests: &[("tests/it/issue_1040_inv_overrecognition_probes.rs", "issue_1034_minimal_contiguous_run_stays_delins"), ("tests/it/issue_1040_inv_overrecognition_probes.rs", "issue_1040_compound_allele_near_3prime_end_stays_delins"), ("tests/it/issue_1040_inv_overrecognition_probes.rs", "issue_1040_contiguous_run_near_3prime_end_stays_delins"), ("tests/it/issue_1040_inv_overrecognition_probes.rs", "issue_1040_literal_ten_nt_contiguous_run_stays_delins"), ("tests/it/issue_1040_inv_overrecognition_probes.rs", "issue_1040_same_run_as_delins_input_stays_delins"), ("tests/it/issue_1040_inv_overrecognition_probes.rs", "length_changing_revcomp_delins_stays_delins"), ("tests/it/issue_1040_inv_overrecognition_probes.rs", "revcomp_runs_in_distinct_codons_are_individual_on_cds_axis"), ("tests/it/issue_1040_inv_overrecognition_probes.rs", "revcomp_runs_separated_by_identity_are_individual_on_genomic_axis"), ("tests/it/issue_1040_inv_overrecognition_probes.rs", "shared_affix_trimming_reveals_inner_inv"), ("tests/it/issue_1040_inv_overrecognition_probes.rs", "two_nt_reverse_complement_emits_inv"), ("tests/it/issue_1040_inv_overrecognition_probes.rs", "whole_run_reverse_complement_emits_inv")], note: Some("Cluster A: the new splitter allows compensating gaps and shreds a contiguous replacement into pieces.") },
    MstoIssue { number: 1041, high: true, state: State::Closed, tests: &[("tests/it/issue_1041_repro.rs", "issue_1041_compound_subs_near_3prime_end_merge_to_inv"), ("tests/it/issue_1041_repro.rs", "issue_1041_root_cause_also_restores_3prime_shift"), ("tests/it/issue_1041_repro.rs", "issue_1041_span_past_contig_end_passes_through_unchanged"), ("tests/it/issue_1041_repro.rs", "issue_1041_whole_run_revcomp_is_inv_regardless_of_3prime_proximity")], note: Some("Cluster A: the new splitter allows compensating gaps and shreds a contiguous replacement into pieces.") },
    MstoIssue { number: 1050, high: false, state: State::Closed, tests: &[("tests/python/test_issue_1050_protein_render_style.py", "test_default_is_three_letter_ter"), ("tests/python/test_issue_1050_protein_render_style.py", "test_constructor_protein_stop_star_switches_stop_token"), ("tests/python/test_issue_1050_protein_render_style.py", "test_per_call_option_b_overrides_default")], note: None },
    MstoIssue { number: 1051, high: false, state: State::Closed, tests: &[("tests/it/gene_selector_display_preserve.rs", "display_drops_redundant_gene_selector_with_compound_ref"), ("tests/it/gene_selector_roundtrip.rs", "display_drops_gene_selector_in_compact_cis_allele"), ("tests/it/gene_selector_roundtrip.rs", "display_drops_gene_selector_refseq_nm"), ("tests/it/gene_selector_roundtrip.rs", "display_drops_redundant_gene_selector_with_compound_ref"), ("tests/it/gene_selector_roundtrip.rs", "parse_captures_gene_symbol_refseq_nm"), ("tests/it/gene_selector_roundtrip.rs", "reparse_preserves_gene_symbol_on_nontranscript_reference")], note: None },
    MstoIssue { number: 1052, high: true, state: State::Closed, tests: &[("tests/it/issue_1052_substitution_refseq.rs", "correct_ref_cds_sub_no_warning"), ("tests/it/issue_1052_substitution_refseq.rs", "correct_ref_exonic_sub_at_transcript_boundary_no_capability_warning"), ("tests/it/issue_1052_substitution_refseq.rs", "correct_ref_exonic_sub_at_transcript_boundary_passes_strict"), ("tests/it/issue_1052_substitution_refseq.rs", "correct_ref_genomic_sub_no_warning"), ("tests/it/issue_1052_substitution_refseq.rs", "correct_ref_mito_sub_no_warning"), ("tests/it/issue_1052_substitution_refseq.rs", "genomic_sub_without_reference_passes_through"), ("tests/it/issue_1052_substitution_refseq.rs", "intronic_sub_bare_nm_stays_on_pre_existing_eintronic_warning_only"), ("tests/it/issue_1052_substitution_refseq.rs", "intronic_sub_without_genomic_sequence_is_silent_passthrough"), ("tests/it/issue_1052_substitution_refseq.rs", "uncertain_correct_ref_genomic_sub_preserves_wrapper"), ("tests/it/issue_1052_substitution_refseq.rs", "uncertain_wrong_ref_cds_sub_stays_silent_in_lenient"), ("tests/it/issue_1052_substitution_refseq.rs", "uncertain_wrong_ref_cds_sub_stays_silent_in_strict"), ("tests/it/issue_1052_substitution_refseq.rs", "uncertain_wrong_ref_genomic_sub_stays_silent_in_lenient_and_strict"), ("tests/it/issue_1052_substitution_refseq.rs", "uncertain_wrong_ref_mito_sub_stays_silent"), ("tests/it/issue_1052_substitution_refseq.rs", "wrong_ref_cds_sub_warns"), ("tests/it/issue_1052_substitution_refseq.rs", "wrong_ref_genomic_sub_rejects_in_strict"), ("tests/it/issue_1052_substitution_refseq.rs", "wrong_ref_genomic_sub_warns_in_lenient"), ("tests/it/issue_1052_substitution_refseq.rs", "wrong_ref_mito_sub_rejects_in_strict"), ("tests/it/issue_1052_substitution_refseq.rs", "wrong_ref_mito_sub_warns_in_lenient")], note: None },
    MstoIssue { number: 1059, high: true, state: State::Closed, tests: &[("src/hgvs/variant.rs", "test_coordinate_axis_allele_shares_member_axis"), ("src/hgvs/variant.rs", "test_coordinate_axis_enum_groupings"), ("src/hgvs/variant.rs", "test_coordinate_axis_genome_ring_is_genomic"), ("src/hgvs/variant.rs", "test_coordinate_axis_leaf_kinds"), ("src/hgvs/variant.rs", "test_coordinate_axis_markers_contribute_no_axis"), ("src/hgvs/variant.rs", "test_coordinate_axis_mixed_allele_is_none"), ("src/hgvs/variant.rs", "test_coordinate_axis_nested_allele")], note: None },
    MstoIssue { number: 1070, high: true, state: State::Closed, tests: &[("src/project/projector.rs", "classify_cis_member_intronic_span"), ("src/project/projector.rs", "project_cis_exon_boundary_member_forces_fallback"), ("src/project/projector.rs", "project_cis_inframe_with_pure_intron_member"), ("src/project/projector.rs", "project_cis_init_member_wins_over_combined"), ("src/project/projector.rs", "project_cis_insertion_member_nets_in_frame"), ("src/project/projector.rs", "project_cis_net_in_frame_allele_is_not_frameshift"), ("src/project/projector.rs", "project_cis_net_in_frame_but_mid_shift_stop_takes_the_frameshift_arm"), ("src/project/projector.rs", "project_cis_overlapping_frameshift_members_fall_back"), ("src/project/projector.rs", "project_cis_separated_inframe_deletions_stay_individual"), ("src/project/projector.rs", "project_trans_net_shifting_allele_is_frameshift"), ("src/project/protein/indel.rs", "cis_combined_rejects_out_of_bounds_span_without_panic")], note: None },
    MstoIssue { number: 1076, high: true, state: State::Closed, tests: &[("src/project/projector.rs", "combine_by_codon_adjacent_codons_form_delins"), ("src/project/projector.rs", "combine_by_codon_same_codon_missense"), ("src/project/projector.rs", "combine_by_codon_same_codon_nonsense"), ("src/project/projector.rs", "combine_by_codon_same_codon_synonymous"), ("src/project/projector.rs", "combine_by_codon_same_position_overlap_returns_none"), ("src/project/projector.rs", "combine_by_codon_separated_codons_stay_bracketed"), ("src/project/projector.rs", "combine_by_codon_single_member_returns_none"), ("src/project/projector.rs", "combine_by_codon_synonymous_group_drops_leaving_single"), ("src/project/projector.rs", "project_cis_intron_split_codon_combines_to_single_missense"), ("src/project/projector.rs", "project_cis_same_codon_all_three_bases_combine_to_single_missense"), ("src/project/projector.rs", "project_cis_same_codon_bases_1_and_2_combine_to_single_missense"), ("src/project/projector.rs", "project_cis_same_codon_bases_1_and_3_combine_to_single_missense"), ("src/project/projector.rs", "project_cis_same_codon_combines_to_nonsense"), ("src/project/projector.rs", "project_cis_same_codon_minus_strand_combines_to_single_missense"), ("src/project/projector.rs", "project_cis_same_position_overlap_stays_bracketed"), ("src/project/projector.rs", "project_cis_terminal_stop_member_keeps_extension_form")], note: None },
    MstoIssue { number: 1097, high: true, state: State::Closed, tests: &[("tests/it/issue_1097_range_sub_refseq_reject.rs", "correct_ref_range_sub_normalizes_cleanly"), ("tests/it/issue_1097_range_sub_refseq_reject.rs", "per_base_allele_subs_wrong_ref_rejects"), ("tests/it/issue_1097_range_sub_refseq_reject.rs", "single_base_sub_wrong_ref_rejects"), ("tests/it/issue_1097_range_sub_refseq_reject.rs", "three_base_range_sub_wrong_ref_rejects"), ("tests/it/issue_1097_range_sub_refseq_reject.rs", "two_base_range_sub_wrong_ref_rejects")], note: None },
    MstoIssue { number: 1098, high: true, state: State::Closed, tests: &[("tests/it/allele_grammar_corners.rs", "protein_descending_order_members_coalesce_to_the_same_delins"), ("tests/it/allele_grammar_corners.rs", "protein_same_residue_changes_are_not_coalesced"), ("tests/it/cis_allele_confluence_proptest.rs", "an_indel_haplotype_is_authored_order_independent"), ("tests/it/issue_1116_protein_coalesce_order_independence.rs", "a_gapped_run_stays_a_bracket_in_both_orders"), ("tests/it/issue_221_cis_three_plus_no_merge.rs", "all_permutations_normalize_to_same_canonical_string"), ("tests/it/issue_221_cis_three_plus_no_merge.rs", "descending_position_order_sorts_to_ascending")], note: None },
    MstoIssue { number: 1099, high: false, state: State::Closed, tests: &[("src/project/projector.rs", "combine_by_codon_consecutive_silent_codons_form_a_range"), ("src/project/projector.rs", "combine_by_codon_ignores_a_member_that_rewrites_nothing"), ("src/project/projector.rs", "combine_by_codon_separated_silent_codons_are_bracketed"), ("src/project/projector.rs", "combine_by_residue_all_silent_names_every_rewritten_codon"), ("src/project/projector.rs", "combine_by_residue_ignores_a_member_that_rewrites_nothing"), ("tests/it/issue_1099_all_silent_protein.rs", "all_silent_cis_allele_names_every_rewritten_codon"), ("tests/it/issue_1099_all_silent_protein.rs", "all_silent_members_in_one_codon_stay_a_single_residue"), ("tests/it/issue_1099_all_silent_protein.rs", "an_unchanged_interior_codon_is_not_named"), ("tests/it/issue_1099_all_silent_protein.rs", "emitted_identities_reparse"), ("tests/it/issue_1099_all_silent_protein.rs", "emitted_identities_survive_normalization_unchanged"), ("tests/it/issue_1099_all_silent_protein.rs", "gene_symbol_rendering_matches_the_changed_residue_sibling"), ("tests/it/issue_1099_all_silent_protein.rs", "output_is_independent_of_member_order"), ("tests/it/issue_1099_all_silent_protein.rs", "separated_silent_codons_render_as_a_bracket"), ("tests/it/issue_1099_all_silent_protein.rs", "silent_delins_names_its_codons_instead_of_whole_protein_identity"), ("tests/it/issue_1099_all_silent_protein.rs", "single_silent_substitution_is_unchanged"), ("tests/it/issue_1099_all_silent_protein.rs", "the_codon_level_combiner_also_brackets_separated_codons"), ("tests/it/issue_1099_all_silent_protein.rs", "three_consecutive_silent_codons_render_as_one_range"), ("tests/it/protein_render_policy.rs", "a_silent_protein_allele_honours_the_render_style")], note: None },
    MstoIssue { number: 1102, high: false, state: State::Closed, tests: &[], note: None },
    MstoIssue { number: 1139, high: false, state: State::Closed, tests: &[("tests/it/gene_selector_roundtrip.rs", "assembly_ref_keeps_its_gene_selector"), ("tests/it/issue_1139_gene_symbol_specification.rs", "a_reference_without_a_gene_symbol_is_unchanged"), ("tests/it/issue_1139_gene_symbol_specification.rs", "gene_symbol_is_not_stacked_as_a_second_specification_on_the_coding_axis"), ("tests/it/issue_1139_gene_symbol_specification.rs", "gene_symbol_is_not_stacked_as_a_second_specification_on_the_noncoding_axis"), ("tests/it/issue_1139_gene_symbol_specification.rs", "projected_coding_name_round_trips_through_the_parser"), ("tests/it/issue_1139_gene_symbol_specification.rs", "projected_noncoding_name_round_trips_through_the_parser"), ("tests/it/issue_1139_gene_symbol_specification.rs", "the_protein_axis_stays_gene_symbol_free")], note: None },
    MstoIssue { number: 1157, high: true, state: State::Closed, tests: &[("tests/it/issue_1157_delins_reduction_shift.rs", "decomposed_cis_allele_is_idempotent"), ("tests/it/issue_1157_delins_reduction_shift.rs", "decomposed_cis_allele_is_not_collapsed_into_a_spanning_delins"), ("tests/it/issue_1157_delins_reduction_shift.rs", "delins_reducing_to_deletion_is_idempotent"), ("tests/it/issue_1157_delins_reduction_shift.rs", "delins_reducing_to_deletion_is_three_prime_shifted"), ("tests/it/issue_1157_delins_reduction_shift.rs", "delins_reducing_to_deletion_matches_direct_deletion"), ("tests/it/issue_1157_delins_reduction_shift.rs", "delins_reducing_to_duplication_is_idempotent"), ("tests/it/issue_1157_delins_reduction_shift.rs", "delins_reducing_to_duplication_is_three_prime_shifted"), ("tests/it/issue_1157_delins_reduction_shift.rs", "delins_reducing_to_duplication_matches_direct_duplication"), ("tests/it/issue_1157_delins_reduction_shift.rs", "sequence_identical_delins_and_allele_normalize_equal"), ("tests/it/issue_1157_delins_reduction_shift.rs", "single_delins_is_decomposed_at_its_unchanged_bases"), ("tests/it/issue_1157_five_prime_insertion_rotation.rs", "five_prime_multibase_insertion_rotates_and_is_idempotent"), ("tests/it/issue_1157_five_prime_insertion_rotation.rs", "five_prime_rotated_insertion_is_a_fixed_point"), ("tests/it/issue_1158_equivalence_resulting_sequence.rs", "delins_reducing_to_deletion_matches_direct_deletion")], note: Some("Cluster B: a non-unique minimal alignment leaves no dominator, so the block stays one spanning member.") },
    MstoIssue { number: 1158, high: true, state: State::Closed, tests: &[("tests/it/issue_1158_equivalence_resulting_sequence.rs", "decomposed_allele_with_different_sequence_stays_not_equivalent"), ("tests/it/issue_1158_equivalence_resulting_sequence.rs", "delins_and_decomposed_allele_with_same_sequence_are_equivalent"), ("tests/it/issue_1158_equivalence_resulting_sequence.rs", "delins_reducing_to_deletion_matches_direct_deletion"), ("tests/it/issue_1158_equivalence_resulting_sequence.rs", "delins_with_different_resulting_sequence_stays_not_equivalent"), ("tests/it/issue_1158_equivalence_resulting_sequence.rs", "identical_variants_still_identical"), ("tests/it/issue_1158_equivalence_resulting_sequence.rs", "sequence_match_is_equivalent"), ("tests/it/issue_1158_equivalence_resulting_sequence.rs", "substitution_written_two_ways_still_normalized_match")], note: Some("Cluster B: confluence diverges rather than relocates (NormalizedMatch degrades to SequenceMatch).") },
    MstoIssue { number: 1159, high: false, state: State::Closed, tests: &[("tests/python/test_reference_aware_spdi.py", "test_the_module_level_conversion_still_declines"), ("tests/python/test_reference_aware_spdi.py", "test_canonical_spdi_is_the_same_for_two_encodings")], note: None },
    MstoIssue { number: 1229, high: true, state: State::Closed, tests: &[("tests/it/cis_allele_confluence_proptest.rs", "members_are_disjoint_and_ascending"), ("tests/it/issue_1235_cis_allele_confluence.rs", "both_public_exits_emit_the_same_canonical_variant"), ("tests/it/issue_1235_cis_allele_confluence.rs", "issue_1229_adjacent_inv_member_coalesces"), ("tests/it/issue_1249_inv_one_base_residue.rs", "every_spelling_of_the_variant_normalizes_alike")], note: None },
    MstoIssue { number: 1230, high: true, state: State::Closed, tests: &[("tests/it/issue_1235_cis_allele_confluence.rs", "issue_1230_inv_with_unchanged_interior_stays_inv")], note: None },
    MstoIssue { number: 1231, high: true, state: State::Closed, tests: &[("tests/it/issue_1235_cis_allele_confluence.rs", "issue_1231_dup_del_reduces_to_substitutions")], note: Some("Cluster B: a non-unique minimal alignment leaves no dominator, so the block stays one spanning member.") },
    MstoIssue { number: 1232, high: true, state: State::Closed, tests: &[("src/normalize/merge.rs", "shadow_merges_across_one_unchanged_base_where_live_splits"), ("tests/it/issue_1235_cis_allele_confluence.rs", "issue_1232_spanning_delins_splits_at_unchanged_base"), ("tests/it/issue_1235_cis_allele_confluence.rs", "soft_masked_reference_yields_the_same_canonical_form")], note: Some("Cluster B: a non-unique minimal alignment leaves no dominator, so the block stays one spanning member.") },
    MstoIssue { number: 1233, high: true, state: State::Closed, tests: &[("tests/it/issue_1235_cis_allele_confluence.rs", "issue_1233_ins_del_reduces_to_substitutions")], note: Some("Cluster B: a non-unique minimal alignment leaves no dominator, so the block stays one spanning member.") },
    MstoIssue { number: 1234, high: true, state: State::Closed, tests: &[("tests/it/cis_allele_confluence_proptest.rs", "members_are_disjoint_and_ascending"), ("tests/it/issue_1234_sibling_clamped_shift.rs", "a_distant_sibling_does_not_clamp_the_shift"), ("tests/it/issue_1234_sibling_clamped_shift.rs", "a_lone_deletion_still_shifts_to_its_three_prime_most_position"), ("tests/it/issue_1234_sibling_clamped_shift.rs", "an_uncertain_allele_keeps_its_predicted_marker"), ("tests/it/issue_1234_sibling_clamped_shift.rs", "deletion_shift_stops_before_a_sibling_substitution"), ("tests/it/issue_1234_sibling_clamped_shift.rs", "normalization_preserves_the_resulting_sequence"), ("tests/it/issue_1234_sibling_clamped_shift.rs", "normalized_cis_members_never_overlap"), ("tests/it/issue_1234_sibling_clamped_shift.rs", "the_certain_spelling_of_that_allele_still_clamps"), ("tests/it/issue_1234_sibling_clamped_shift.rs", "the_clamped_spelling_reaches_the_same_string"), ("tests/it/issue_1235_cis_allele_confluence.rs", "issue_1234_lone_deletion_still_shifts_fully"), ("tests/it/issue_1235_cis_allele_confluence.rs", "normalization_preserves_the_resulting_sequence")], note: Some("Cluster C: the fix lives in a repair pass (clamp_sibling_crossing_shifts) the rewrite's splitter never reaches.") },
    MstoIssue { number: 1235, high: true, state: State::Closed, tests: &[("src/normalize/merge.rs", "shadow_merges_across_one_unchanged_base_where_live_splits"), ("tests/it/cis_allele_confluence_proptest.rs", "members_are_disjoint_and_ascending"), ("tests/it/idempotency_tests.rs", "test_overlap_conflicting_allele_is_idempotent_across_respellings"), ("tests/it/issue_1235_cis_allele_confluence.rs", "both_public_exits_emit_the_same_canonical_variant"), ("tests/it/issue_1235_cis_allele_confluence.rs", "issue_1229_adjacent_inv_member_coalesces"), ("tests/it/issue_1235_cis_allele_confluence.rs", "issue_1230_inv_with_unchanged_interior_stays_inv"), ("tests/it/issue_1235_cis_allele_confluence.rs", "issue_1231_dup_del_reduces_to_substitutions"), ("tests/it/issue_1235_cis_allele_confluence.rs", "issue_1232_spanning_delins_splits_at_unchanged_base"), ("tests/it/issue_1235_cis_allele_confluence.rs", "issue_1233_ins_del_reduces_to_substitutions"), ("tests/it/issue_1235_cis_allele_confluence.rs", "issue_1234_lone_deletion_still_shifts_fully"), ("tests/it/issue_1235_cis_allele_confluence.rs", "normalization_is_idempotent"), ("tests/it/issue_1235_cis_allele_confluence.rs", "normalization_preserves_the_resulting_sequence"), ("tests/it/issue_1235_cis_allele_confluence.rs", "normalized_cis_members_are_disjoint_and_ordered"), ("tests/it/issue_1235_cis_allele_confluence.rs", "shifted_insertion_piece_rotates_its_payload"), ("tests/it/issue_1235_cis_allele_confluence.rs", "soft_masked_reference_yields_the_same_canonical_form"), ("tests/it/issue_1235_transcript_axes.rs", "cds_axis_converges_on_separated_changes"), ("tests/it/issue_1235_transcript_axes.rs", "cds_axis_reduces_a_mixed_type_allele_to_substitutions"), ("tests/it/issue_1235_transcript_axes.rs", "cds_axis_splits_across_a_codon_boundary"), ("tests/it/issue_1235_transcript_axes.rs", "coding_rna_axis_keeps_its_reading_frame"), ("tests/it/issue_1235_transcript_axes.rs", "noncoding_axis_converges_on_separated_changes"), ("tests/it/issue_1235_transcript_axes.rs", "noncoding_rna_axis_converges_without_a_reading_frame"), ("tests/it/issue_1235_transcript_axes.rs", "overlap_conflicting_allele_is_not_canonicalized"), ("tests/it/issue_1235_transcript_axes.rs", "rna_axis_converges_and_keeps_the_u_alphabet"), ("tests/it/issue_1235_transcript_axes.rs", "rna_axis_reduces_a_mixed_type_allele_and_keeps_the_u_alphabet"), ("tests/it/issue_1249_inv_one_base_residue.rs", "every_spelling_of_the_variant_normalizes_alike"), ("tests/it/issue_1325_repeat_growth_swallows_junction.rs", "a_repeat_growing_past_its_tract_releases_a_swallowed_junction"), ("tests/it/issue_1394_repeat_growth_swallows_deletion.rs", "a_repeat_growing_past_its_tract_releases_a_swallowed_deletion")], note: Some("Cluster B (c./r. axes): a non-unique minimal alignment leaves no dominator, so the block stays spanning.") },
    MstoIssue { number: 1244, high: true, state: State::Closed, tests: &[("tests/it/issue_1244_equivalence_overlap_panic.rs", "a_disjoint_cis_allele_is_not_declined"), ("tests/it/issue_1244_equivalence_overlap_panic.rs", "an_overlapping_allele_compared_with_itself_is_still_decided"), ("tests/it/issue_1244_equivalence_overlap_panic.rs", "an_overlapping_cis_allele_returns_a_verdict_instead_of_panicking"), ("tests/it/issue_1244_equivalence_overlap_panic.rs", "members_that_merely_abut_are_not_treated_as_overlapping"), ("tests/it/issue_1244_equivalence_overlap_panic.rs", "normalizing_the_same_overlapping_allele_does_not_panic"), ("tests/it/issue_1244_equivalence_overlap_panic.rs", "pure_insertions_at_one_position_are_not_treated_as_overlapping"), ("tests/it/issue_1244_equivalence_overlap_panic.rs", "the_crash_is_symmetric_in_argument_order"), ("tests/it/issue_1244_equivalence_overlap_panic.rs", "three_insertions_at_one_position_are_still_not_overlapping")], note: None },
    MstoIssue { number: 1245, high: false, state: State::Closed, tests: &[("tests/python/test_issue_1245_enum_contract.py", "test_members_are_hashable"), ("tests/python/test_issue_1245_enum_contract.py", "test_members_work_as_set_and_dict_keys"), ("tests/python/test_issue_1245_enum_contract.py", "test_members_are_orderable")], note: None },
    MstoIssue { number: 1246, high: false, state: State::Closed, tests: &[("tests/python/test_issue_1246_changed_positions_list.py", "test_changed_positions_is_a_list"), ("tests/python/test_issue_1246_changed_positions_list.py", "test_changed_positions_holds_ints"), ("tests/python/test_issue_1246_changed_positions_list.py", "test_changed_positions_is_json_serializable")], note: None },
    // Filed 2026-08-05. #1419–#1421 are the surviving non-confluence follow-ups
    // to the #1229–#1234 family; #1430 proposes coercing a description to its
    // denoted sequence before normalizing, which is the same seam the
    // sequence-first splitter sits on.
    //
    // All four were listed as unguarded on the per-function attribution rule.
    // Three of them were not: `reported_confluence_pairs` and
    // `reported_partition_verdicts` are table-driven and carry the issue numbers
    // in their rows and module docs, which is the accommodation the module doc
    // above sets out. The shared guards are listed against each of the three
    // because each module's assertions sweep every row; the two per-issue
    // extras (`reported_spans_change_the_columns_the_reports_state` for #1420,
    // `the_1421_spans_separate_by_two_nucleotides_not_one` for #1421) test only
    // that issue's rows and are listed only there.
    //
    // #1802 renamed one shared guard and added one, so all three lists move
    // together: `each_pair_reaches_its_canonical_form_from_exactly_one_spelling`
    // became `each_pair_reaches_its_wanted_form_from_the_spellings_its_state_names`
    // when the pair model became a closed three-state enumeration, and
    // `the_pair_state_census_holds` is the new census over that enumeration. Both
    // sweep every row, so both are shared. That is the intentional addition
    // PINNED_UNIQUE_TEST_COUNT moved 289 -> 290 for.
    MstoIssue { number: 1419, high: true, state: State::Open, tests: &[("tests/it/reported_confluence_pairs.rs", "every_reported_output_is_a_fixed_point"), ("tests/it/reported_confluence_pairs.rs", "every_reported_pair_denotes_one_sequence"), ("tests/it/reported_confluence_pairs.rs", "no_reported_pair_normalizes_to_a_different_sequence"), ("tests/it/reported_confluence_pairs.rs", "the_reported_pair_census_is_unchanged"), ("tests/it/reported_partition_verdicts.rs", "each_pair_reaches_its_wanted_form_from_the_spellings_its_state_names"), ("tests/it/reported_partition_verdicts.rs", "every_canonical_row_already_prints_its_wanted_form"), ("tests/it/reported_partition_verdicts.rs", "every_gap_row_is_returned_exactly_as_authored"), ("tests/it/reported_partition_verdicts.rs", "every_reported_pair_is_still_one_variant_by_equivalence"), ("tests/it/reported_partition_verdicts.rs", "every_reported_spelling_normalizes_as_recorded_under_five_prime"), ("tests/it/reported_partition_verdicts.rs", "every_reported_spelling_still_normalizes_as_the_reports_recorded"), ("tests/it/reported_partition_verdicts.rs", "only_the_named_rows_answer_differently_in_the_two_directions"), ("tests/it/reported_partition_verdicts.rs", "reference_bases_are_what_the_reports_state"), ("tests/it/reported_partition_verdicts.rs", "the_open_gap_census_holds"), ("tests/it/reported_partition_verdicts.rs", "the_pair_state_census_holds"), ("tests/it/reported_partition_verdicts.rs", "the_spec_authority_census_holds"), ("tests/it/reported_partition_verdicts.rs", "the_two_reported_modules_describe_the_same_pairs"), ("tests/it/reported_partition_verdicts.rs", "the_wanted_form_of_every_gap_row_is_its_siblings_output")], note: None },
    // #1420's `the_1420_v2_pair_does_not_converge_by_re_derivation` was DELETED by
    // operator ruling of 2026-08-13 — it read separation off the input's spelling,
    // which `rulings[separation-is-a-property-of-the-spelling-not-of-the-variant]`
    // forbids. The issue keeps its slice through the remaining table-driven
    // guards; `the_reported_pair_census_is_unchanged` now carries the row, whose
    // pair converges at 9 of 9.
    MstoIssue { number: 1420, high: true, state: State::Open, tests: &[("tests/it/reported_confluence_pairs.rs", "every_reported_output_is_a_fixed_point"), ("tests/it/reported_confluence_pairs.rs", "every_reported_pair_denotes_one_sequence"), ("tests/it/reported_confluence_pairs.rs", "no_reported_pair_normalizes_to_a_different_sequence"), ("tests/it/reported_confluence_pairs.rs", "the_reported_pair_census_is_unchanged"), ("tests/it/reported_partition_verdicts.rs", "each_pair_reaches_its_wanted_form_from_the_spellings_its_state_names"), ("tests/it/reported_partition_verdicts.rs", "every_canonical_row_already_prints_its_wanted_form"), ("tests/it/reported_partition_verdicts.rs", "every_gap_row_is_returned_exactly_as_authored"), ("tests/it/reported_partition_verdicts.rs", "every_reported_pair_is_still_one_variant_by_equivalence"), ("tests/it/reported_partition_verdicts.rs", "every_reported_spelling_normalizes_as_recorded_under_five_prime"), ("tests/it/reported_partition_verdicts.rs", "every_reported_spelling_still_normalizes_as_the_reports_recorded"), ("tests/it/reported_partition_verdicts.rs", "only_the_named_rows_answer_differently_in_the_two_directions"), ("tests/it/reported_partition_verdicts.rs", "reference_bases_are_what_the_reports_state"), ("tests/it/reported_partition_verdicts.rs", "reported_spans_change_the_columns_the_reports_state"), ("tests/it/reported_partition_verdicts.rs", "the_open_gap_census_holds"), ("tests/it/reported_partition_verdicts.rs", "the_pair_state_census_holds"), ("tests/it/reported_partition_verdicts.rs", "the_spec_authority_census_holds"), ("tests/it/reported_partition_verdicts.rs", "the_two_reported_modules_describe_the_same_pairs"), ("tests/it/reported_partition_verdicts.rs", "the_wanted_form_of_every_gap_row_is_its_siblings_output")], note: None },
    MstoIssue { number: 1421, high: true, state: State::Open, tests: &[("tests/it/reported_confluence_pairs.rs", "every_reported_output_is_a_fixed_point"), ("tests/it/reported_confluence_pairs.rs", "every_reported_pair_denotes_one_sequence"), ("tests/it/reported_confluence_pairs.rs", "no_reported_pair_normalizes_to_a_different_sequence"), ("tests/it/reported_confluence_pairs.rs", "the_reported_pair_census_is_unchanged"), ("tests/it/reported_partition_verdicts.rs", "each_pair_reaches_its_wanted_form_from_the_spellings_its_state_names"), ("tests/it/reported_partition_verdicts.rs", "every_canonical_row_already_prints_its_wanted_form"), ("tests/it/reported_partition_verdicts.rs", "every_gap_row_is_returned_exactly_as_authored"), ("tests/it/reported_partition_verdicts.rs", "every_reported_pair_is_still_one_variant_by_equivalence"), ("tests/it/reported_partition_verdicts.rs", "every_reported_spelling_normalizes_as_recorded_under_five_prime"), ("tests/it/reported_partition_verdicts.rs", "every_reported_spelling_still_normalizes_as_the_reports_recorded"), ("tests/it/reported_partition_verdicts.rs", "only_the_named_rows_answer_differently_in_the_two_directions"), ("tests/it/reported_partition_verdicts.rs", "reference_bases_are_what_the_reports_state"), ("tests/it/reported_partition_verdicts.rs", "the_1421_spans_separate_by_two_nucleotides_not_one"), ("tests/it/reported_partition_verdicts.rs", "the_open_gap_census_holds"), ("tests/it/reported_partition_verdicts.rs", "the_pair_state_census_holds"), ("tests/it/reported_partition_verdicts.rs", "the_spec_authority_census_holds"), ("tests/it/reported_partition_verdicts.rs", "the_two_reported_modules_describe_the_same_pairs"), ("tests/it/reported_partition_verdicts.rs", "the_wanted_form_of_every_gap_row_is_its_siblings_output")], note: None },
    // #1430's slice is one or more guards per proposed stage, not a
    // reproduction — see
    // `issue_1430_is_measured_by_a_property_not_by_an_expected_output`.
    MstoIssue { number: 1430, high: true, state: State::Open, tests: &[("src/normalize/merge.rs", "an_all_survivor_payload_is_not_an_alignment_artifact"), ("src/normalize/merge.rs", "the_coalesce_pass_reaches_the_spec_worked_example"), ("src/normalize/seqfirst/partition.rs", "canonical_members_claim_exactly_the_blocks_edit_distance"), ("src/normalize/seqfirst/partition.rs", "canonical_members_rebuild_the_alternate_block"), ("src/normalize/seqfirst/partition.rs", "round_trips_exhaustively_over_a_small_alphabet"), ("tests/it/reported_confluence_pairs.rs", "every_reported_pair_denotes_one_sequence"), ("tests/it/reported_confluence_pairs.rs", "no_reported_pair_normalizes_to_a_different_sequence"), ("tests/it/reported_confluence_pairs.rs", "the_reported_pair_census_is_unchanged")], note: None },
];

/// #1430's slice, split by the stage of its proposal each guard measures.
///
/// Pinned so the entry cannot be quietly converted into an ordinary
/// input/expected row: #1430 states no expected output, so anything that looked
/// like one would be somebody's invention attributed to the issue.
const ISSUE_1430_STAGE_GUARDS: &[(u8, &str, &str)] = &[
    // Stage 1 — coerce the description to its denoted sequence, deterministically
    // and confluently. The applier in `common::cis_apply_oracle` *is* that
    // coercion, and it is independent of the normalizer, so these two say the
    // sequence a description denotes is well defined and survives normalization.
    (
        1,
        "tests/it/reported_confluence_pairs.rs",
        "every_reported_pair_denotes_one_sequence",
    ),
    (
        1,
        "tests/it/reported_confluence_pairs.rs",
        "no_reported_pair_normalizes_to_a_different_sequence",
    ),
    // Stage 2 — derive a representative variant from that sequence,
    // deterministically. `seqfirst::partition` is the implementation of stage 2;
    // these pin that its output rebuilds the alternate block, costs exactly the
    // block's edit distance, and round-trips exhaustively over a small alphabet.
    (
        2,
        "src/normalize/seqfirst/partition.rs",
        "canonical_members_rebuild_the_alternate_block",
    ),
    (
        2,
        "src/normalize/seqfirst/partition.rs",
        "canonical_members_claim_exactly_the_blocks_edit_distance",
    ),
    (
        2,
        "src/normalize/seqfirst/partition.rs",
        "round_trips_exhaustively_over_a_small_alphabet",
    ),
    // Stage 3 — apply spec compliance on top of the derived variant. The
    // coalesce pass is the step-3 pass `merge::coalesce` cites #1430 for, and
    // `the_coalesce_pass_reaches_the_spec_worked_example` runs it against
    // `delins.md:44-47`'s own worked example.
    (
        3,
        "src/normalize/merge.rs",
        "the_coalesce_pass_reaches_the_spec_worked_example",
    ),
    (
        3,
        "src/normalize/merge.rs",
        "an_all_survivor_payload_is_not_an_alignment_artifact",
    ),
    // The measure of the whole proposal: does the reported family converge? A
    // census rather than a pin, because #1430's success criterion is confluence
    // and not any particular string.
    (
        0,
        "tests/it/reported_confluence_pairs.rs",
        "the_reported_pair_census_is_unchanged",
    ),
];

/// The eleven issues the sequence-first migration review flagged as at risk,
/// pinned so drift in [`MSTO_ISSUES`] is caught by
/// [`pinned_at_risk_issues_carry_a_mechanism_note`].
const PINNED_AT_RISK_ISSUES: &[u32] = &[
    160, 1034, 1040, 1041, 1157, 1158, 1231, 1232, 1233, 1234, 1235,
];

/// Pinned totals, checked by [`issue_and_test_totals_match_their_pinned_counts`].
/// A change to any of these numbers means [`MSTO_ISSUES`] changed — update the
/// constant only after confirming the new count is intentional (a `gh issue
/// list` re-enumeration, a newly-tagged test, or a rename/deletion this table
/// should reflect), not to silence a failing assertion.
const PINNED_ISSUE_COUNT: usize = 73;
const PINNED_HIGH_ISSUE_COUNT: usize = 27;
/// Raised from 265 when #1419/#1420/#1421/#1430 gained their slices: the guards
/// already existed, the catalog could not see them (see the module doc's
/// "Table-driven guards count" and "An issue whose satisfaction is a property"
/// sections). No test was written to move this number.
///
/// Raised again to 289 by
/// `reported_confluence_pairs::the_1420_v2_pair_does_not_converge_by_re_derivation`
/// — a genuinely new test, and the one exception to the sentence above. It names
/// the string `1420-v2` converges on when the merge-across-unchanged-bases veto
/// stops firing, so the pair cannot be recorded as converging on #1420's own
/// wanted form by way of a `general.md:34` violation.
///
/// Raised to 290 by
/// `reported_partition_verdicts::the_pair_state_census_holds` (#1802), the
/// census over the closed three-state pair model. The same change *renamed*
/// `each_pair_reaches_its_canonical_form_from_exactly_one_spelling` to
/// `each_pair_reaches_its_wanted_form_from_the_spellings_its_state_names`,
/// which moves no count — see the note above the `#1419`/`#1420`/`#1421` rows,
/// where all three lists carry both edits.
///
/// Lowered to 289 by #1616: the test two paragraphs above,
/// `the_1420_v2_pair_does_not_converge_by_re_derivation`, was DELETED by
/// operator ruling of 2026-08-13. It read separation off the input's spelling,
/// which `rulings[separation-is-a-property-of-the-spelling-not-of-the-variant]`
/// forbids. An intentional removal, not corpus rot: #1420 keeps its slice
/// through the remaining table-driven guards.
///
/// **Both sides of this constant were wrong at the merge, which is why it is
/// re-derived here rather than resolved.** #1616 was written when the pin read
/// 289 and recorded `289 -> 288`; #1802 then raised the base to 290, so the same
/// deletion lands at 289. Two branches carrying different digits for one event
/// is exactly the shape git merges cleanly into a wrong number.
const PINNED_UNIQUE_TEST_COUNT: usize = 289;

/// Reads `relative_path` (relative to the crate root) and returns whether it
/// contains a function definition named `fn_name` — i.e. a line containing
/// `fn <fn_name>` or `def <fn_name>` immediately followed by `(` or `<`
/// (accounting for generics), possibly preceded by `pub`/`async`/whitespace.
fn source_defines_fn(relative_path: &str, fn_name: &str) -> bool {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join(relative_path);
    let content = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("could not read {relative_path}: {e}"));
    // `def` as well as `fn`: the bindings' guards are pytest modules, and a
    // catalog that cannot see them reports guarded issues as unguarded.
    let pattern = format!(r"(?:fn|def)\s+{}\s*[(<]", regex::escape(fn_name));
    Regex::new(&pattern)
        .expect("pattern built from an escaped literal is always valid")
        .is_match(&content)
}

/// Every `(file, fn)` pair cataloged in [`MSTO_ISSUES`] must still name a real
/// function in its source file. A pair that fails here means the guarding
/// test was renamed or deleted without updating this table — the exact
/// silent-rot failure mode this module exists to catch.
#[test]
fn every_cataloged_test_exists_in_its_source_file() {
    let mut missing = Vec::new();
    for issue in MSTO_ISSUES {
        for &(file, fn_name) in issue.tests {
            if !source_defines_fn(file, fn_name) {
                missing.push(format!("#{}: {file}::{fn_name}", issue.number));
            }
        }
    }
    assert!(
        missing.is_empty(),
        "cataloged tests no longer found in their source files (renamed, moved, or deleted; \
         update MSTO_ISSUES in tests/it/msto_regression_corpus.rs):\n{}",
        missing.join("\n")
    );
}

/// Pins the issue count, the `high`-label count, and the number of distinct
/// tests cataloged across every issue, so an edit to [`MSTO_ISSUES`] that
/// silently changes one of these totals shows up as a failing assertion
/// instead of only as a diff someone has to notice.
#[test]
fn issue_and_test_totals_match_their_pinned_counts() {
    assert_eq!(
        MSTO_ISSUES.len(),
        PINNED_ISSUE_COUNT,
        "msto issue count changed \u{2014} re-run `gh issue list -R fulcrumgenomics/ferro-hgvs \
         --state all --author msto` and update MSTO_ISSUES and PINNED_ISSUE_COUNT"
    );

    let high_count = MSTO_ISSUES.iter().filter(|i| i.high).count();
    assert_eq!(
        high_count, PINNED_HIGH_ISSUE_COUNT,
        "count of `high`-labeled msto issues changed \u{2014} update PINNED_HIGH_ISSUE_COUNT \
         after confirming the new count against GitHub"
    );

    let unique_tests: HashSet<(&str, &str)> = MSTO_ISSUES
        .iter()
        .flat_map(|issue| issue.tests.iter().copied())
        .collect();
    assert_eq!(
        unique_tests.len(),
        PINNED_UNIQUE_TEST_COUNT,
        "total distinct cataloged tests changed \u{2014} update PINNED_UNIQUE_TEST_COUNT after \
         confirming the new count is an intentional addition/removal, not corpus rot"
    );

    // No duplicate issue numbers, and no test pair claimed by an issue twice.
    let mut seen_numbers = HashSet::new();
    for issue in MSTO_ISSUES {
        assert!(
            seen_numbers.insert(issue.number),
            "issue #{} appears more than once in MSTO_ISSUES",
            issue.number
        );
        let mut seen_pairs = HashSet::new();
        for &pair in issue.tests {
            assert!(
                seen_pairs.insert(pair),
                "issue #{} lists {}::{} more than once",
                issue.number,
                pair.0,
                pair.1
            );
        }
    }
}

/// Exactly the eleven issues the sequence-first migration review named as at
/// risk carry a `note`; no other issue does. Catches both directions of
/// drift: an at-risk issue silently losing its note, and an unrelated issue
/// picking one up by copy-paste.
#[test]
fn pinned_at_risk_issues_carry_a_mechanism_note() {
    assert_eq!(PINNED_AT_RISK_ISSUES.len(), 11);

    let noted: HashSet<u32> = MSTO_ISSUES
        .iter()
        .filter(|i| i.note.is_some())
        .map(|i| i.number)
        .collect();
    let pinned: HashSet<u32> = PINNED_AT_RISK_ISSUES.iter().copied().collect();
    assert_eq!(
        noted, pinned,
        "the set of issues carrying a mechanism `note` no longer matches PINNED_AT_RISK_ISSUES"
    );

    // Ten of the eleven pinned issues are `high`; #160 is the review's sole
    // non-`high` exception (it surfaced from a full-suite run, not the
    // `high`-labeled backlog).
    let high_count = PINNED_AT_RISK_ISSUES
        .iter()
        .filter(|&&number| {
            MSTO_ISSUES
                .iter()
                .find(|i| i.number == number)
                .unwrap_or_else(|| panic!("pinned at-risk issue #{number} is not in MSTO_ISSUES"))
                .high
        })
        .count();
    assert_eq!(
        high_count, 10,
        "expected exactly ten of the eleven pinned at-risk issues to be labeled `high` \
         (#160 is the sole non-high exception)"
    );
}

/// Exactly the four issues known to be open are marked as such:
/// #1419/#1420/#1421 (the non-confluence follow-ups to the #1229–#1234 family)
/// and #1430 (the sequence-coercion proposal).
///
/// Two issues left the open set in this re-enumeration. #1235, the rewrite's own
/// tracking issue, closed as completed. #1159, the SPDI feature request, closed
/// 2026-08-06 — after the first pass of this catalog was written, which is why
/// the accompanying re-enumeration was re-run against GitHub before merge rather
/// than trusted from the earlier `gh issue list`.
///
/// The "is the state one of the two legal values" half of this check is gone,
/// because `State` makes an illegal value unrepresentable — the compiler owns it
/// now.
#[test]
fn issue_state_matches_the_known_open_set() {
    let open: HashSet<u32> = MSTO_ISSUES
        .iter()
        .filter(|i| i.state == State::Open)
        .map(|i| i.number)
        .collect();
    assert_eq!(
        open,
        HashSet::from([1419, 1420, 1421, 1430]),
        "the set of open msto issues changed \u{2014} re-check against GitHub and update \
         MSTO_ISSUES"
    );
}

/// A `high`-labeled issue with no cataloged test is a finding, not an error
/// in this table: it means no test in the corpus can be traced to that issue
/// by number. Pinned to the set actually found so that closing the gap (or
/// opening a new one) is visible rather than silently absorbed.
#[test]
fn no_high_priority_issue_is_unguarded_without_being_named() {
    let unguarded: Vec<u32> = MSTO_ISSUES
        .iter()
        .filter(|i| i.high && i.tests.is_empty())
        .map(|i| i.number)
        .collect();
    assert_eq!(
        unguarded,
        vec![200],
        "the set of unguarded high-priority msto issues changed; if this is a newly-found gap, \
         leave it here as a finding rather than silently dropping it"
    );
}

/// #1430 is satisfied by a property holding, not by an output matching, and its
/// slice records that rather than pretending otherwise.
///
/// Every other entry in [`MSTO_ISSUES`] points at tests that reproduce a
/// reported input and check what comes out. #1430 has no reported input: it
/// proposes a three-stage pipeline (coerce to sequence, derive a representative
/// variant, apply spec compliance) and argues for it on determinism and
/// confluence. So the slice names one or more guards **per stage**, and this
/// test pins that structure — that all three stages are covered, that the two
/// tables agree, and that the whole-proposal measure is a confluence census
/// rather than a string.
///
/// Without this, the natural next edit is for somebody to "fix" the odd-looking
/// entry by writing an input/expected row for #1430. That row would be an
/// invention attributed to an issue that states no expected output, which is a
/// worse failure than the empty slice this replaced.
#[test]
fn issue_1430_is_measured_by_a_property_not_by_an_expected_output() {
    let issue = MSTO_ISSUES
        .iter()
        .find(|i| i.number == 1430)
        .expect("#1430 is in MSTO_ISSUES");

    let from_stages: HashSet<(&str, &str)> = ISSUE_1430_STAGE_GUARDS
        .iter()
        .map(|&(_, file, fn_name)| (file, fn_name))
        .collect();
    let from_slice: HashSet<(&str, &str)> = issue.tests.iter().copied().collect();
    assert_eq!(
        from_slice, from_stages,
        "#1430's cataloged tests and ISSUE_1430_STAGE_GUARDS disagree; every guard \
         claimed for #1430 must be attributable to a stage of its proposal"
    );

    for stage in [1u8, 2, 3] {
        assert!(
            ISSUE_1430_STAGE_GUARDS.iter().any(|&(s, _, _)| s == stage),
            "no guard is claimed for stage {stage} of #1430's proposal; an \
             uncovered stage is a finding to leave visible, not one to absorb"
        );
    }

    // Stage 0 is the whole-proposal measure, and it must be the confluence
    // census — the only assertion in the family that reports whether two
    // spellings reach one string, which is #1430's actual success criterion.
    let overall: Vec<&str> = ISSUE_1430_STAGE_GUARDS
        .iter()
        .filter(|&&(stage, _, _)| stage == 0)
        .map(|&(_, _, fn_name)| fn_name)
        .collect();
    assert_eq!(
        overall,
        vec!["the_reported_pair_census_is_unchanged"],
        "#1430's overall measure must stay the confluence census; it has no \
         expected output to pin instead"
    );
}
