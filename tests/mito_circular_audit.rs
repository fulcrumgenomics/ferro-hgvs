//! Audit and behavior-pinning tests for mitochondrial (`m.`) circular
//! reference handling — issue #81 item F1 (wraparound at the origin).
//!
//! Tracking issue: **#129** ("Mitochondrial circular reference
//! handling (m. coordinate system) (#81 F1)"). PR #106 explicitly
//! excluded the `m.` dimension from the `del` 3'-shift coverage matrix
//! because "circular wraparound is its own design"; this file pins
//! what ferro does today so any future design lands as an intentional
//! behavior change with a visible diff.
//!
//! ### Background — the HGVS spec
//!
//! - The mitochondrial genome (rCRS, NC_012920.1) is **circular**, length 16569.
//!   Position 16569 is followed by position 1 with no gap. See
//!   `assets/hgvs-nomenclature/docs/background/numbering.md` lines 6–8.
//! - The HGVS recommendations therefore allow ranges that **cross the
//!   origin**, written as `m.<high>_<low>` where `high > low`. The
//!   exception is called out only on the **deletion** page:
//!   `assets/hgvs-nomenclature/docs/recommendations/DNA/deletion.md`
//!   line 18 says "when a circular genomic reference sequence is used
//!   (`o` and `m` prefix) nucleotide positions may be listed from 3' to
//!   5' when the deletion includes both the last and first nucleotides
//!   of the reference sequence". The corresponding text is **not**
//!   present in the `dup` / `ins` / `inv` / `delins` pages — that
//!   asymmetry is one of the open ambiguities F1 has to resolve.
//! - Normalization on a circular reference also has its own rules:
//!   3'-shifting may have to wrap, and "leftmost-shifted" / "rightmost-
//!   shifted" lose their absolute meaning — they are defined modulo
//!   the contig length.
//!
//! ### What the spec corpus exercises today
//!
//! `tests/fixtures/grammar/hgvs_spec_normalization.json` carries rows
//! for `NC_012920.1:m.3243A>G`, `NC_012920.1:m.3460G>A`, and
//! `NC_012920.1(MT-ND1):m.3460del` (with and without gene symbol).
//! **Zero** rows in the spec corpus exercise wraparound, full-mtDNA
//! span, or any other circular-topology semantic. Anything F1
//! introduces will need a corresponding extension of the spec-corpus
//! fixtures — that is also out of scope here.
//!
//! ### What ferro does today (May 2026, post-#129 Path 1)
//!
//! Partial handling of circular references:
//!
//! - `src/hgvs/parser/variant.rs::parse_genome_interval_for_circular`
//!   plus the `check_circular_reversed_range` post-parse validator
//!   accept wraparound `del`, `delins`, and `dup` ranges on `m.`/`o.`
//!   per SVD-WG006 (and `deletion.md:17`). `ins` and `inv` wraparound
//!   forms are rejected by `check_circular_reversed_range` — the spec
//!   is silent on them, and ferro prefers a clear parse error to the
//!   prior silent-accept with mis-computed indel length. (#129 Path 1;
//!   see the file-level assertions for each edit kind.)
//! - `src/normalize/mod.rs::normalize_mt` now mirrors `normalize_genome`
//!   for non-origin-crossing variants (issue #210): window-based 3'
//!   shuffle and repeat-notation canonicalization, with
//!   `is_coding=false`. Wraparound variants either fail to parse
//!   (`del`/`delins`/`sub`) or fall through to the no-reference
//!   `canonicalize_mt_variant` fallback (`dup`/`ins`/`inv`), preserving
//!   the pinned behavior below until F1/#129 introduces circular-aware
//!   semantics. The audit assertions in this file pin the wraparound
//!   parse / pass-through surface — they continue to hold because the
//!   tests use `MockProvider::new()` (no reference data) and bare
//!   substitutions short-circuit before normalization.
//! - `src/vcf/from_hgvs.rs` synthesizes a `GenomeVariant` from the
//!   `MtVariant` and routes it through the linear genomic conversion
//!   path; nothing knows the contig is circular.
//! - `src/spdi/convert.rs::hgvs_to_spdi_simple` rejects all
//!   non-`Genome` variants outright — including `Mt`. Section 13 pins
//!   that policy.
//! - `src/python_helpers.rs::compute_span` has a defensive `None` for
//!   origin-crossing **circular (`o.`)** intervals where `end < start`,
//!   but the same guard does not exist for `m.`. Section 6 pins the
//!   resulting silent miscompute.
//! - There is no contig-length validation: `m.16570A>G` (one past the
//!   end of NC_012920.1) parses cleanly today. Section 7 pins this.
//!
//! ### Cross-cuts
//!
//! - **F2 (heteroplasmy, PR #139, tracking #133).** Heteroplasmy uses
//!   `=/` (mosaic) and `=//` (chimeric) at the variant level, not as
//!   a circular concern. The only cross-cut is via compound alleles
//!   (Section 8 below): all `m.[…]` allele forms are rejected today,
//!   so there is no path on which F1 and F2 can collide here.
//! - **PR #106 (`del` 3'-shift coverage matrix).** Excluded the `m.`
//!   axis with the explicit note "circular wraparound is its own
//!   design"; this audit and #129 are the bookkeeping for that note.
//!
//! ### What this file pins
//!
//! Each test parses a string and asserts on the **observed outcome**
//! (Ok with a specific render, or Err). When F1 lands, the assertions
//! that pin "this fails today" must flip to "this succeeds with the
//! expected canonical form" — that diff is the audit trail.

use ferro_hgvs::hgvs::AllelePhase;
use ferro_hgvs::python_helpers::get_indel_length;
use ferro_hgvs::{hgvs_to_spdi_simple, parse_hgvs, HgvsVariant, MockProvider, Normalizer};

// -----------------------------------------------------------------------------
// Scenario 1 — non-wrapping substitution.
//
// Sanity check: a plain `m.` substitution that does not cross the
// origin parses, round-trips, and passes through `normalize` unchanged.
// This is the baseline that any F1 implementation must preserve.
// -----------------------------------------------------------------------------

/// Non-wrapping `m.` substitution parses and round-trips today.
#[test]
fn audit_mt_non_wrapping_sub_parses_and_round_trips() {
    let input = "NC_012920.1:m.100A>G";
    let variant =
        parse_hgvs(input).unwrap_or_else(|e| panic!("expected parse to succeed for {input}: {e}"));
    assert_eq!(format!("{}", variant), input);
}

/// Non-wrapping `m.` substitution passes through `normalize` unchanged
/// (substitutions short-circuit before normalization regardless of the
/// `normalize_mt` body).
#[test]
fn audit_mt_non_wrapping_sub_normalizes_to_self() {
    let input = "NC_012920.1:m.100A>G";
    let variant = parse_hgvs(input).expect("parse should succeed");
    let normalizer = Normalizer::new(MockProvider::new());
    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize should succeed for non-wrapping m. substitution");
    assert_eq!(format!("{}", normalized), input);
}

// -----------------------------------------------------------------------------
// Scenario 2 — wrapping deletion at the origin: `m.16569_1del`.
//
// Spec-canonical wraparound: a 2-nt deletion across positions 16569 and
// 1 (in that circular order). Today the parser rejects this at the
// inverted-range check. When F1 lands, this should parse to a
// circular-aware interval and normalize accordingly.
// -----------------------------------------------------------------------------

/// `m.16569_1del` is the SVD-WG006 example — spec-authorised wraparound
/// del. Per HGVS v19.01 + `deletion.md:17`, must parse and round-trip.
/// (#129 flipped this from the prior parser-rejection assertion.)
#[test]
fn audit_mt_wrapping_del_at_origin_parses() {
    let input = "NC_012920.1:m.16569_1del";
    let variant =
        parse_hgvs(input).expect("spec-authorised wraparound del must parse (#129, SVD-WG006)");
    assert_eq!(
        format!("{}", variant),
        input,
        "wraparound del must round-trip 3'→5' on m.",
    );
}

// -----------------------------------------------------------------------------
// Scenario 3 — longer wrapping deletion: `m.16569_5del`.
//
// 6-nt deletion spanning the origin (positions 16569, 1, 2, 3, 4, 5).
// Same parser rejection as scenario 2; pinned separately because a
// future F1 implementation may distinguish the 2-nt edge case from the
// general N-nt wraparound case.
// -----------------------------------------------------------------------------

/// Longer wraparound deletion also parses per the same SVD-WG006
/// exception. (#129 flipped this from the prior rejection assertion.)
#[test]
fn audit_mt_wrapping_del_longer_parses() {
    let input = "NC_012920.1:m.16569_5del";
    let variant =
        parse_hgvs(input).expect("longer wraparound del must parse on m. (#129, SVD-WG006)");
    assert_eq!(format!("{}", variant), input);
}

// -----------------------------------------------------------------------------
// Scenario 4 — full mtDNA span as a duplication: `m.1_16569dup`.
//
// Not strictly a wraparound (1 ≤ 16569), but represents the entire
// circular contig as a single duplication. It exercises the same
// "circular contig length" semantics: in a linear view the result is a
// 33138-nt sequence, but on a circular reference it is degenerate
// (every base duplicated relative to itself). F1 must decide how to
// handle this; today it parses (because `dup` is exempt from the
// inverted-range check, and 1 ≤ 16569 anyway). With #210 in place,
// `normalize_mt` would try the window fetch, but `MockProvider::new()`
// has no NC_012920.1 sequence so it falls back to
// `canonicalize_mt_variant` (minimal-notation cleanup, no shift).
// -----------------------------------------------------------------------------

/// Full mtDNA span as a duplication parses today; pin the round-trip
/// to expose any future change in canonical form.
#[test]
fn audit_mt_full_span_dup_round_trips_today() {
    let input = "NC_012920.1:m.1_16569dup";
    let variant = parse_hgvs(input)
        .unwrap_or_else(|e| panic!("PINNED: expected parse to succeed today: {e}"));
    assert_eq!(
        format!("{}", variant),
        input,
        "PINNED: ferro currently round-trips full-mtDNA dup verbatim; \
         F1 may decide to emit a circular-aware canonical form (e.g. \
         reject, warn, or rewrite as `=`)."
    );
}

/// Full mtDNA span passes through `normalize` unchanged today.
#[test]
fn audit_mt_full_span_dup_normalizes_to_self_today() {
    let input = "NC_012920.1:m.1_16569dup";
    let variant = parse_hgvs(input).expect("parse should succeed");
    let normalizer = Normalizer::new(MockProvider::new());
    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize should succeed (no-reference fallback)");
    assert_eq!(
        format!("{}", normalized),
        input,
        "PINNED: without provider data `normalize_mt` falls back to \
         minimal-notation cleanup, which preserves a count-less span dup \
         verbatim. F1 must define what circular-aware normalization \
         does for a full-mtDNA span."
    );
}

// -----------------------------------------------------------------------------
// Scenario 5 — wrapping non-`del` edits at the origin.
//
// The spec text in `recommendations/DNA/deletion.md` line 18 calls out
// 3'→5' position order as a deletion-specific exception for circular
// references. The `dup`, `ins`, `inv`, `delins` pages do **not** repeat
// the exception. Ferro's parser today follows the inverse asymmetry:
//   - `del` and `delins` reject inverted ranges (Verify error)
//   - `dup`, `ins`, `inv` accept inverted ranges (parser carve-out for
//     "structural" edits, unrelated to circular semantics)
// F1 must close this gap one way or the other: either lift the `del` /
// `delins` rejection so wraparound deletions parse, or tighten the
// `dup` / `ins` / `inv` carve-out so wraparound forms are not silently
// accepted without circular-aware normalize support.
// -----------------------------------------------------------------------------

/// Wrapping `delins` at the origin parses per the same SVD-WG006
/// exception (delins inherits via del+ins composition). (#129 flipped
/// this from the prior rejection assertion.)
#[test]
fn audit_mt_wrapping_delins_parses() {
    let input = "NC_012920.1:m.16569_1delinsT";
    let variant = parse_hgvs(input).expect("wraparound delins must parse on m. (#129, SVD-WG006)");
    assert_eq!(format!("{}", variant), input);
}

/// Wrapping `dup` at the origin parses per SVD-WG006:23 (explicit
/// example `J01749.1:o.4344_197dup`) and SVD-WG006:33 ("simplifies
/// the description of deletions/duplications"). (#129 made this
/// official-rather-than-silently-accepted; the round-trip and
/// no-reference-fallback behaviour are unchanged.)
#[test]
fn audit_mt_wrapping_dup_parses() {
    let input = "NC_012920.1:m.16560_5dup";
    let variant =
        parse_hgvs(input).expect("spec-authorised wraparound dup must parse (#129, SVD-WG006)");
    assert_eq!(format!("{}", variant), input);
    let normalizer = Normalizer::new(MockProvider::new());
    let normalized = normalizer
        .normalize(&variant)
        .expect("no-reference fallback normalize should succeed");
    assert_eq!(
        format!("{}", normalized),
        input,
        "PINNED: with no provider data, wraparound dup falls back to \
         minimal-notation cleanup and is preserved unchanged. \
         Wraparound 3'-shift is still a no-op (matches mutalyzer + \
         biocommons + strict spec reading).",
    );
}

/// Wrapping `inv` is rejected per #129 (spec-silent on inv wraparound;
/// SVD-WG006 names only del/dup in its rationale and examples).
/// Was previously silently accepted with mis-computed length.
#[test]
fn audit_mt_wrapping_inv_rejected() {
    let input = "NC_012920.1:m.16500_500inv";
    let result = parse_hgvs(input);
    assert!(
        result.is_err(),
        "wraparound inv has no spec exception (SVD-WG006 names only \
         del/dup); must be rejected. Got {result:?}",
    );
}

/// Wrapping `ins` is rejected per #129 (spec-silent on ins
/// wraparound). Was previously silently accepted.
#[test]
fn audit_mt_wrapping_ins_rejected() {
    let input = "NC_012920.1:m.16569_1insT";
    let result = parse_hgvs(input);
    assert!(
        result.is_err(),
        "wraparound ins has no spec exception (SVD-WG006 names only \
         del/dup); must be rejected. Got {result:?}",
    );
}

// -----------------------------------------------------------------------------
// Scenario 6 — silent miscompute on wraparound `dup`.
//
// `python_helpers::compute_span` has a `Circular(_) && end < start →
// None` guard but no `Mt` guard. Because wraparound `dup` parses today
// (Scenario 5), a downstream caller asking `get_indel_length` for a
// wrapping m. dup will get a *negative or otherwise wrong* number,
// silently. F1 must either extend the guard to `Mt` or compute the
// real circular span using contig length. Pin today's value so any
// fix surfaces.
// -----------------------------------------------------------------------------

/// Wraparound `dup` returns `None` from the no-provider
/// `get_indel_length` since `compute_span` now guards `Mt(_)`. Pin so
/// any regression (silent wrong value) surfaces. The spec-correct
/// value (+15) is asserted via `get_indel_length_with_provider` in
/// `tests/issue_399_mt_circular_followup.rs`.
#[test]
fn audit_mt_wrapping_dup_indel_length_is_none_without_provider() {
    let input = "NC_012920.1:m.16560_5dup";
    let variant = parse_hgvs(input).expect("parse should succeed per #129 / SVD-WG006");
    assert_eq!(
        get_indel_length(&variant),
        None,
        "PINNED: wraparound m. dup yields None without a provider; \
         the provider-backed variant in tests/issue_399_mt_circular_followup.rs \
         asserts the spec-correct +15."
    );
}

/// Wraparound `del` returns `None` from the no-provider `get_indel_length`
/// since `compute_span` now guards `Mt(_)`. Pin so any regression (silent
/// wrong value) surfaces. The spec-correct value (-2) is asserted via
/// `get_indel_length_with_provider` in `tests/issue_399_mt_circular_followup.rs`.
#[test]
fn audit_mt_wrapping_del_indel_length_is_none_without_provider() {
    let input = "NC_012920.1:m.16569_1del";
    let variant = parse_hgvs(input).expect("parse should succeed per #129 / SVD-WG006");
    assert_eq!(
        get_indel_length(&variant),
        None,
        "PINNED: wraparound m. del yields None without a provider; \
         the provider-backed variant in tests/issue_399_mt_circular_followup.rs \
         asserts the spec-correct -2."
    );
}

/// Wraparound `delins` returns `None` from the no-provider `get_indel_length`
/// since `compute_span` now guards `Mt(_)`. Pin so any regression (silent
/// wrong value) surfaces. The spec-correct value (-1) is asserted via
/// `get_indel_length_with_provider` in `tests/issue_399_mt_circular_followup.rs`.
#[test]
fn audit_mt_wrapping_delins_indel_length_is_none_without_provider() {
    let input = "NC_012920.1:m.16569_1delinsT";
    let variant = parse_hgvs(input).expect("parse should succeed per #129 / SVD-WG006");
    assert_eq!(
        get_indel_length(&variant),
        None,
        "PINNED: wraparound m. delins yields None without a provider; \
         the provider-backed variant in tests/issue_399_mt_circular_followup.rs \
         asserts the spec-correct -1."
    );
}

// -----------------------------------------------------------------------------
// Scenario 7 — out-of-bounds and zero/negative positions.
//
// Spec: numbering is 1..=L (L = 16569 for NC_012920.1); there is no
// position 0 (`numbering.md` line 20 documents the convention for
// coding sequences; the same 1-based no-zero rule applies to `m.`).
// Ferro today rejects 0 and negative but **accepts** positions past
// the contig end. F1 must decide whether to add length validation —
// it will need contig length anyway for wraparound math.
// -----------------------------------------------------------------------------

/// `m.16570A>G` (one past the end of NC_012920.1) parses today; ferro
/// has no contig-length validation.
#[test]
fn audit_mt_position_past_end_silently_accepted_today() {
    let input = "NC_012920.1:m.16570A>G";
    let variant = parse_hgvs(input).unwrap_or_else(|e| {
        panic!(
            "PINNED: ferro has no length validation for m. accessions; \
             m.16570A>G parses today. Got Err({e})."
        )
    });
    assert_eq!(format!("{}", variant), input);
}

/// `m.0A>G` is rejected today (1-based numbering, no zero).
/// Baseline that any F1 implementation must preserve.
#[test]
fn audit_mt_position_zero_rejected_today() {
    let input = "NC_012920.1:m.0A>G";
    assert!(
        parse_hgvs(input).is_err(),
        "PINNED: m. position 0 is rejected; preserve under F1."
    );
}

/// `m.-1A>G` is rejected today (negative positions not legal on `m.`).
#[test]
fn audit_mt_negative_position_rejected_today() {
    let input = "NC_012920.1:m.-1A>G";
    assert!(
        parse_hgvs(input).is_err(),
        "PINNED: m. negative position is rejected; preserve under F1."
    );
}

// -----------------------------------------------------------------------------
// Scenario 8 — compound alleles on `m.`.
//
// Spec permits compound HGVS allele descriptions in cis (`[a;b]`),
// trans (`[a];[b]`), and the homozygous shorthand `[a](;)`. Ferro's
// parser accepts trans (`[a];[b]`) on `m.` (PR #146 / sibling C1) and
// cis (`[a;b]`) plus unknown-phase (`[a(;)b]`) on `m.` (PR #148 /
// sibling C3, issue #123). The bracket-suffix homozygous shorthand
// (`[a](;)`) is still rejected — it is a separate construct from the
// bracketless unknown-phase form covered by #123.
// -----------------------------------------------------------------------------

/// Compound `m.[…;…]` (cis) parses today (PR #148 / issue #123) as an
/// `AlleleVariant` in `Cis` phase.
#[test]
fn audit_mt_compound_cis_allele_accepted_today() {
    let input = "NC_012920.1:m.[100A>G;200T>C]";
    let variant = parse_hgvs(input).expect("compound m. cis allele should parse");
    let HgvsVariant::Allele(allele) = &variant else {
        panic!("expected HgvsVariant::Allele, got {variant:?}");
    };
    assert_eq!(allele.phase, AllelePhase::Cis);
    assert_eq!(allele.variants.len(), 2);
    assert!(matches!(allele.variants[0], HgvsVariant::Mt(_)));
    assert!(matches!(allele.variants[1], HgvsVariant::Mt(_)));
}

/// Compound `m.[…];[…]` (trans) parses today as an `AlleleVariant` in
/// `Trans` phase. PR #146 (sibling C1) added the compact-prefix trans
/// shorthand for `m.` to parallel `g.`/`n.`. Heteroplasmy is the
/// canonical use case (different alleles on different mtDNA copies).
#[test]
fn audit_mt_compound_trans_allele_accepted_today() {
    let input = "NC_012920.1:m.[100A>G];[200T>C]";
    let variant = parse_hgvs(input).expect("compound m. trans allele should parse");
    let HgvsVariant::Allele(allele) = &variant else {
        panic!("expected HgvsVariant::Allele, got {variant:?}");
    };
    assert_eq!(allele.phase, AllelePhase::Trans);
    assert_eq!(allele.variants.len(), 2);
    assert!(matches!(allele.variants[0], HgvsVariant::Mt(_)));
    assert!(matches!(allele.variants[1], HgvsVariant::Mt(_)));
}

/// `m.[…](;)` (homozygous shorthand) is rejected today.
#[test]
fn audit_mt_homozygous_shorthand_allele_rejected_today() {
    let input = "NC_012920.1:m.[100A>G](;)";
    assert!(
        parse_hgvs(input).is_err(),
        "PINNED: m. homozygous-shorthand allele rejected today."
    );
}

// -----------------------------------------------------------------------------
// Scenario 9 — fully-qualified cross-coord allele.
//
// Spec describes "fully qualified" alleles where each variant carries
// its own `accession:type.` prefix; ferro accepts mixing `m.` and `g.`
// in such an allele today. Pinned to expose any future tightening
// (e.g. rejecting cross-coord alleles entirely, or normalizing them).
// -----------------------------------------------------------------------------

/// Cross-coord allele `[NC_012920.1:m.…;NC_000001.11:g.…]` parses and
/// round-trips today.
#[test]
fn audit_mt_cross_coord_allele_round_trips_today() {
    let input = "[NC_012920.1:m.100A>G;NC_000001.11:g.100A>G]";
    let variant = parse_hgvs(input)
        .unwrap_or_else(|e| panic!("PINNED: cross-coord m./g. allele parses today: Err({e})"));
    assert_eq!(format!("{}", variant), input);
}

// -----------------------------------------------------------------------------
// Scenario 10 — non-canonical mtDNA accession aliases.
//
// Spec prefers `NC_012920.1` (`numbering.md` line 8). Ferro is
// permissive: `chrM`, `chrMT`, and `MT` are all accepted as
// accessions. Pinned because F1 may want to whitelist mtDNA accessions
// (the implementation needs to know "is this contig circular?" anyway,
// which is a per-accession decision).
// -----------------------------------------------------------------------------

/// `chrM:m.100A>G` parses today as `HgvsVariant::Mt`.
#[test]
fn audit_mt_accession_chrm_lower_round_trips_today() {
    let input = "chrM:m.100A>G";
    let variant =
        parse_hgvs(input).unwrap_or_else(|e| panic!("PINNED: chrM accession parses today: {e}"));
    assert!(matches!(variant, HgvsVariant::Mt(_)));
    assert_eq!(format!("{}", variant), input);
}

/// `chrMT:m.100A>G` parses today as `HgvsVariant::Mt`.
#[test]
fn audit_mt_accession_chrmt_round_trips_today() {
    let input = "chrMT:m.100A>G";
    let variant =
        parse_hgvs(input).unwrap_or_else(|e| panic!("PINNED: chrMT accession parses today: {e}"));
    assert!(matches!(variant, HgvsVariant::Mt(_)));
    assert_eq!(format!("{}", variant), input);
}

/// `MT:m.100A>G` parses today as `HgvsVariant::Mt`.
#[test]
fn audit_mt_accession_mt_round_trips_today() {
    let input = "MT:m.100A>G";
    let variant =
        parse_hgvs(input).unwrap_or_else(|e| panic!("PINNED: MT accession parses today: {e}"));
    assert!(matches!(variant, HgvsVariant::Mt(_)));
    assert_eq!(format!("{}", variant), input);
}

// -----------------------------------------------------------------------------
// Scenario 11 — coding-style position decorations on `m.`.
//
// mtDNA has no introns conventionally and no UTRs (`m.` has no `+`,
// `-`, or `*` per `numbering.md` line 7). Ferro today nonetheless
// accepts `m.100+1A>G` and `m.100-1A>G` (intronic-style offsets),
// because the genome-position parser is shared with `g.` and does not
// know `m.` should reject them. `m.*100A>G` (UTR-style) is rejected,
// fortuitously, because `*` is not in the `m.`-position grammar.
// F1 may want to reject `+`/`-` offsets explicitly.
// -----------------------------------------------------------------------------

/// `m.100+1A>G` (intronic-style positive offset) parses today.
#[test]
fn audit_mt_intronic_plus_offset_silently_accepted_today() {
    let input = "NC_012920.1:m.100+1A>G";
    let variant = parse_hgvs(input).unwrap_or_else(|e| {
        panic!(
            "PINNED: m. intronic-style `+` offset parses today even though \
             mtDNA has no introns. Got Err({e})."
        )
    });
    assert_eq!(format!("{}", variant), input);
}

/// `m.100-1A>G` (intronic-style negative offset) parses today.
#[test]
fn audit_mt_intronic_minus_offset_silently_accepted_today() {
    let input = "NC_012920.1:m.100-1A>G";
    let variant = parse_hgvs(input).unwrap_or_else(|e| {
        panic!(
            "PINNED: m. intronic-style `-` offset parses today even though \
             mtDNA has no introns. Got Err({e})."
        )
    });
    assert_eq!(format!("{}", variant), input);
}

/// `m.*100A>G` (UTR-style asterisk position) is rejected today;
/// preserve as the correct behavior under F1.
#[test]
fn audit_mt_utr_asterisk_position_rejected_today() {
    let input = "NC_012920.1:m.*100A>G";
    assert!(
        parse_hgvs(input).is_err(),
        "PINNED: m. UTR-style `*` position rejected; preserve under F1."
    );
}

// -----------------------------------------------------------------------------
// Scenario 12 — gene-symbol preserved on round-trip.
//
// `NC_012920.1(MT-ND1):m.3460G>A` (a row that exists in the spec
// normalization fixture) parses today and the rendered output
// preserves the `(MT-ND1)` qualifier byte-for-byte. This is unrelated
// to wraparound but is surfaced by the same input class — pinned here
// so any future regression in `m.`-variant Display surfaces as a
// visible diff.
// -----------------------------------------------------------------------------

/// Gene-symbol qualifier on `m.` is preserved on round-trip.
#[test]
fn audit_mt_gene_symbol_preserved_on_round_trip() {
    let input = "NC_012920.1(MT-ND1):m.3460G>A";
    let variant = parse_hgvs(input).expect("parse should succeed");
    assert_eq!(
        format!("{}", variant),
        input,
        "PINNED: the `(MT-ND1)` gene-symbol qualifier is preserved on \
         render for m. variants."
    );
}

// -----------------------------------------------------------------------------
// Scenario 13 — SPDI conversion on `m.`.
//
// `hgvs_to_spdi_simple` now accepts plain non-wrapping `m.` substitutions
// and emits SPDI on the mitochondrial accession directly (closed by #116,
// matching VCF's treat-as-linear policy from `src/vcf/from_hgvs.rs`).
// Wraparound m. variants where linear semantics break down remain a
// separate F1 question.
// -----------------------------------------------------------------------------

/// `hgvs_to_spdi_simple` accepts a plain non-wrapping `m.` substitution
/// post-#116 and emits SPDI on the mt accession directly. Pinned so any
/// regression that re-rejects m. → SPDI surfaces immediately.
#[test]
fn audit_mt_spdi_conversion_accepted_after_issue_116() {
    let variant = parse_hgvs("NC_012920.1:m.100A>G").expect("parse should succeed");
    let spdi = hgvs_to_spdi_simple(&variant).expect("m. → SPDI must succeed post-#116");
    assert_eq!(spdi.sequence, "NC_012920.1");
    assert_eq!(spdi.position, 99); // 1-based 100 → 0-based 99
    assert_eq!(spdi.deletion, "A");
    assert_eq!(spdi.insertion, "G");
}
