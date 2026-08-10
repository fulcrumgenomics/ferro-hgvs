//! #1578: every post-parse semantic validator must reach a ring's `::`-joined
//! segments and a supernumerary marker's inner variant.
//!
//! A ring's segments are `LocEdit<GenomeInterval, NaEdit>` hanging off
//! `GenomeRing` rather than nested `HgvsVariant`s, so no validator reaches
//! them by matching on variant kind — they are validated only because the
//! parse driver's `for_each_leaf` walk re-wraps each segment as an ordinary
//! `g.` leaf (and unwraps the `sup` marker) before running the leaf rules.
//! These tests pin that delivery for each rule that had a measured escape
//! (#1578; the insertion-anchor rules are pinned alongside #1264's tests).
//!
//! Parser tests are the only guard that can: an accepted escape round-trips —
//! ferro parses it *and* renders it back — so none of the three normalization
//! oracles can see this class.
//!
//! Both ring spellings are covered throughout: the bare `::`-join and the
//! supernumerary-wrapped `g.[..::..]sup`.
//!
//! Each rejection test asserts error **parity** with the standalone spelling
//! of the same violation — not just that the ring is rejected, but that it is
//! rejected by the *same* validator the standalone input hits. That is the
//! invariant `for_each_leaf` guarantees: a ring segment is validated exactly
//! like a standalone description. Where parity deliberately does not hold
//! byte-for-byte, the test pins the divergence and says why.

use ferro_hgvs::parse_hgvs;

/// Assert that both ring spellings of a violation fail with byte-identical
/// diagnostics to the standalone spelling — the invariant `for_each_leaf`
/// guarantees: a ring segment is validated exactly like a standalone
/// description.
fn assert_ring_error_parity(standalone: &str, ring_inputs: [&str; 2]) {
    let standalone_err = parse_hgvs(standalone)
        .expect_err("the standalone spelling is pinned as rejected")
        .to_string();
    for input in ring_inputs {
        let ring_err = parse_hgvs(input)
            .expect_err(&format!(
                "`{input}` must be rejected like its standalone spelling"
            ))
            .to_string();
        assert_eq!(
            ring_err, standalone_err,
            "`{input}` must fail with the same diagnostic as `{standalone}` — \
             the walker validates ring segments exactly like standalone descriptions",
        );
    }
}

/// `dupins` is "a format not used in HGVS nomenclature"
/// (`DNA/duplication.md:92` — `<code class="invalid">` markup, the spec's
/// strongest prohibition register). The bare spelling is rejected; a ring
/// segment spelling it must be rejected too. Adjudicated-correct.
#[test]
fn a_ring_segment_dupins_is_rejected() {
    assert_ring_error_parity(
        "NC_000022.11:g.100_101dupinsATG",
        [
            "NC_000022.11:g.100_101dupinsATG::200_201insG",
            "NC_000022.11:g.[100_101dupinsATG::200_201insG]sup",
        ],
    );
}

/// A size-count suffix on a single-position anchor (`g.123del3`, `g.123dup3`)
/// is not allowed; the canonical form names both endpoints
/// (`checklist.md:49`, `DNA/deletion.md:117`, `DNA/duplication.md:140`).
/// Adjudicated-correct.
///
/// The `del` and `dup` sub-cases behave differently here and are checked
/// separately:
///
/// - `dup` carries a fully static message (the spec declines to disambiguate
///   where the duplication starts, so there is no canonical rewrite to name —
///   see `validate_no_point_size_suffix`'s `SizedEdit::Duplication` arm), so
///   ring and standalone match byte-for-byte.
/// - `del` names a canonical rewrite computed by
///   `suggest_point_del_range_form`, which scans the *whole* original input
///   for a `del<size>` token terminated by `)`/`;`/`]`/whitespace/end-of-string
///   (`corrections::point_size_token`). Inside a ring the token is instead
///   followed by `::`, which is not in that terminator set, so the scan finds
///   no rewritable token and the message falls back to the generic
///   `` `<start>_<end>del`, not `<start>del<size>` `` text — legitimately
///   different from the standalone message, which does compute a concrete
///   rewrite. Both ring spellings produce that identical fallback text, so
///   the `del` assertions below compare ring against sup-ring exactly and
///   fall back to a shared substring against standalone.
#[test]
fn a_ring_segment_size_count_suffix_is_rejected() {
    assert_ring_error_parity(
        "NC_000022.11:g.123dup3",
        [
            "NC_000022.11:g.123dup3::200_201insGG",
            "NC_000022.11:g.[123dup3::200_201insGG]sup",
        ],
    );

    const STANDALONE_DEL: &str = "NC_000022.11:g.123del3";
    let standalone_del_err = parse_hgvs(STANDALONE_DEL)
        .expect_err("the standalone del3 spelling is pinned as rejected")
        .to_string();
    // A distinctive substring shared by every del3 message regardless of
    // which canonical-rewrite tail it names — see the fn doc comment.
    const NAMES_ONLY_ONE_ENDPOINT: &str =
        "a deletion sized by a suffix on a single-position anchor names only one endpoint";
    assert!(
        standalone_del_err.contains(NAMES_ONLY_ONE_ENDPOINT),
        "sanity: the standalone del3 message must itself carry the shared substring; \
         got: {standalone_del_err}"
    );

    let ring_del_err = parse_hgvs("NC_000022.11:g.123del3::200_201insGG")
        .expect_err("the ring del3 segment must be rejected like its standalone spelling")
        .to_string();
    let supring_del_err = parse_hgvs("NC_000022.11:g.[123del3::200_201insGG]sup")
        .expect_err("the sup-wrapped ring del3 segment must be rejected too")
        .to_string();
    for (label, err) in [("ring", &ring_del_err), ("sup-ring", &supring_del_err)] {
        assert!(
            err.contains(NAMES_ONLY_ONE_ENDPOINT),
            "{label} del3 must fail on the same size-count rule as `{STANDALONE_DEL}` \
             (substring fallback: the whole-input canonical rewrite legitimately differs \
             inside a `::`-join, see the fn doc comment); got: {err}"
        );
        // Pin the fallback tail too: if the ring path ever regains a concrete
        // whole-input rewrite (e.g. the terminator-set follow-up lands), this
        // divergence is over and the test should be collapsed into
        // `assert_ring_error_parity`, not left half-pinned.
        assert!(
            err.contains("write `<start>_<end>del`, not `<start>del<size>`"),
            "{label} del3 currently falls back to the generic (non-rewritten) repair \
             text; got: {err}"
        );
    }
    // The two ring spellings do agree byte-for-byte with each other, which
    // the substring fallback above does not, by itself, pin.
    assert_eq!(
        ring_del_err, supring_del_err,
        "the bare ring and the sup-wrapped ring must fail identically for `del3` — both \
         fall back to the generic (non-rewritten) message"
    );
}

/// An inversion covers more than one nucleotide by definition; a
/// one-nucleotide inversion is described as a substitution to the
/// complementary base (`DNA/inversion.md:16`). Adjudicated-correct.
#[test]
fn a_ring_segment_single_nucleotide_inversion_is_rejected() {
    assert_ring_error_parity(
        "NC_000022.11:g.100inv",
        [
            "NC_000022.11:g.100inv::200_201insGG",
            "NC_000022.11:g.[100inv::200_201insGG]sup",
        ],
    );
}

/// A substitution replaces one nucleotide by one other
/// (`DNA/substitution.md:30`, `DNA/delins.md:73`); a ring segment spelling a
/// multi-base substitution must be rejected exactly like the standalone
/// spelling. Adjudicated-correct.
///
/// This is the one rejection whose message deliberately differs between the
/// whole-description and member/segment spellings: a leaf that IS the whole
/// description gets a full replacement to paste, but an allele member or ring
/// segment gets only its edit respelled — offering the repaired leaf as a
/// whole-variant replacement would tell the reader to silently drop the
/// sibling members or ring segments (and any `sup` marker). Pinned below on
/// both ring spellings and on an allele member.
#[test]
fn a_ring_segment_multibase_substitution_is_rejected() {
    const RULE: &str = "naming several reference bases is not allowed";

    let standalone_err = parse_hgvs("NC_000022.11:g.100_101AA>TT")
        .expect_err("the standalone spelling is pinned as rejected")
        .to_string();
    assert!(
        standalone_err.contains(RULE)
            && standalone_err.contains("write `NC_000022.11:g.100_101delinsTT` instead"),
        "the standalone rejection names the full replacement description; got: {standalone_err}"
    );

    let ring_err = parse_hgvs("NC_000022.11:g.100_101AA>TT::200_201insGG")
        .expect_err("the ring segment must be rejected like its standalone spelling")
        .to_string();
    let supring_err = parse_hgvs("NC_000022.11:g.[100_101AA>TT::200_201insGG]sup")
        .expect_err("the sup-wrapped ring segment must be rejected too")
        .to_string();
    for (label, err) in [("ring", &ring_err), ("sup-ring", &supring_err)] {
        assert!(
            err.contains(RULE),
            "{label}: same rule as standalone must fire; got: {err}"
        );
        assert!(
            err.contains("rewrite that edit as `g.100_101delinsTT`"),
            "{label}: the suggestion must respell only the offending segment; got: {err}"
        );
        assert!(
            !err.contains("write `NC_000022.11:g.100_101delinsTT` instead"),
            "{label}: the suggestion must NOT offer a whole-variant replacement that \
             drops the ring's other segments; got: {err}"
        );
    }
    assert_eq!(
        ring_err, supring_err,
        "the bare ring and the sup-wrapped ring must fail identically"
    );

    // The same scoping holds for an allele member: the suggestion respells
    // the offending edit rather than offering a single-member replacement
    // that would drop `;100del`.
    let member_err = parse_hgvs("NM_000088.3:c.[79_80GC>TT;100del]")
        .expect_err("a multibase substitution in an allele member is rejected")
        .to_string();
    assert!(
        member_err.contains(RULE) && member_err.contains("rewrite that edit as `c.79_80delinsTT`"),
        "the member suggestion must respell only the offending edit; got: {member_err}"
    );
    assert!(
        !member_err.contains("write `NM_000088.3:c.79_80delinsTT` instead"),
        "the member suggestion must NOT offer a whole-variant replacement that drops \
         the sibling member; got: {member_err}"
    );
}

/// `U` is an RNA base; DNA reference sequences use A/C/G/T (mutalyzer's
/// `ENODNA`). Adjudicated-correct against the bare form's pinned rejection.
#[test]
fn a_ring_segment_u_in_dna_is_rejected() {
    assert_ring_error_parity(
        "NC_000022.11:g.100_101delinsUU",
        [
            "NC_000022.11:g.100_101delinsUU::200_201insGG",
            "NC_000022.11:g.[100_101delinsUU::200_201insGG]sup",
        ],
    );
}

/// A non-coding RNA reference (NR_/XR_) has no CDS and is not genomic, so a
/// `g.` description on it is a coordinate-system mismatch (mutalyzer's
/// `ECOORDINATESYSTEMMISMATCH`). This is the accession-level check: an
/// entire ring on an NR_ accession escaped. Adjudicated-correct.
#[test]
fn a_ring_on_a_noncoding_rna_accession_is_rejected() {
    assert_ring_error_parity(
        "NR_000022.1:g.100del",
        [
            "NR_000022.1:g.100del::200_201insGG",
            "NR_000022.1:g.[100del::200_201insGG]sup",
        ],
    );
}

/// The editless sup-marked spelling `NR_…:g.pter_qtersup` must be rejected
/// exactly like its bare spelling: `sup` marks the presence of an additional
/// copy, not a different reference, so the coordinate-system rule applies
/// unchanged. Parity with the bare spelling — a pinned coordinate-system
/// mismatch — is the adjudication. (The `sup` inner is only reachable via
/// `for_each_leaf`; this acceptance change is disclosed on #1578.)
#[test]
fn a_supernumerary_marker_on_a_noncoding_rna_accession_is_rejected() {
    let bare_err = parse_hgvs("NR_000022.1:g.pter_qter")
        .expect_err("the bare editless NR_ spelling is pinned as rejected")
        .to_string();
    let sup_err = parse_hgvs("NR_000022.1:g.pter_qtersup")
        .expect_err("the sup-marked spelling must be rejected like the bare one")
        .to_string();
    assert_eq!(
        sup_err, bare_err,
        "`sup` marks presence of the sequence, not a different reference — the \
         coordinate-system rejection must be identical",
    );
}

/// Pins the driver's leaf-major diagnostic precedence (see `for_each_leaf`'s
/// doc): every leaf runs all rules before the next leaf is visited, so for an
/// allele whose members violate *different* rules, member order decides which
/// diagnostic surfaces. Deliberate — the contrasting rule-major order (each
/// rule sweeping every member before the next rule runs) would report the
/// `dupins` error for both orders below, so this test distinguishes the two.
#[test]
fn a_multi_member_allele_reports_the_first_member_in_error() {
    let dupins_first = parse_hgvs("NC_000022.11:g.[100_101dupinsAT;200insG]")
        .expect_err("both members are invalid")
        .to_string();
    assert!(
        dupins_first.contains("dupins"),
        "member 1's dupins violation must surface first; got: {dupins_first}"
    );

    let insertion_first = parse_hgvs("NC_000022.11:g.[200insG;100_101dupinsAT]")
        .expect_err("both members are invalid")
        .to_string();
    assert!(
        insertion_first.contains("single-position insertion"),
        "member 1's insertion violation must surface first; got: {insertion_first}"
    );
}

/// The control: a ring of legal edits stays valid and round-trips, in both
/// spellings. Without this, tightening the ring path could pass by
/// rejecting every ring.
#[test]
fn a_ring_of_legal_edits_still_parses() {
    for input in [
        "NC_000022.11:g.100_101dup::200_201insG",
        "NC_000022.11:g.[100_101dup::200_201insG]sup",
    ] {
        let variant =
            parse_hgvs(input).unwrap_or_else(|e| panic!("`{input}` must stay valid: {e}"));
        assert_eq!(variant.to_string(), input, "`{input}` must round-trip");
    }
}
