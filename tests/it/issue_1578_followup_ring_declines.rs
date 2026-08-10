//! #1578 follow-up (item 4): the three modules that hand-roll `Allele`
//! recursion and were **already** ring-safe, pinned so they stay that way.
//!
//! The follow-up asked whether the structural blindness #1578 fixed in the
//! parser is live in `src/spdi/apply.rs`, `src/project/projector.rs`,
//! `src/vcf/from_hgvs.rs` and `src/equivalence/checker.rs` — each of which
//! matches on `HgvsVariant` by hand and recurses into `Allele` members. It was
//! live in exactly one of them: `equivalence/checker.rs` returned a *positive
//! wrong verdict*, fixed and pinned in
//! `issue_1578_followup_equivalence_rungs.rs`. The other three decline, which
//! is the correct behaviour and is what this file guards.
//!
//! **A passing test here is a claim about today, so it is worth saying what
//! makes it non-vacuous.** Each module reaches its ring decline through a
//! *different* mechanism, and none of the three is a catch-all that would keep
//! working if the ring arm were dropped:
//!
//! | module | how a ring is declined |
//! |---|---|
//! | `spdi/apply.rs` | `variant_edit_triples`' `single => vec![single]` catch-all hands the ring to `hgvs_to_spdi`, which has no ring arm and returns `UnsupportedVariantType`; `.ok()?` turns that into a whole-variant decline |
//! | `vcf/from_hgvs.rs` | explicit `GenomeRing` / `Supernumerary` arms (`convert`), and `convert_all` reaches them per member |
//! | `project/projector.rs` | explicit `GenomeRing` / `Supernumerary` arms in both `extract_contig_and_pos` and `project_to_genomic_nc` |
//!
//! Three shapes are exercised against each: the bare `::` join, the
//! `sup`-wrapped ring, and a ring as a **cis-allele member**. That last one is
//! reachable only by construction — measured at this commit, the parser rejects
//! `g.[100_101del::200_201del;300A>T]`, `g.[300A>T;100_101del::200_201del]` and
//! `g.[[100_101del::200_201del];300A>T]` alike — so it is built by hand here.
//! It is the shape that matters for this follow-up: it is the one that walks a
//! ring *through* the hand-rolled `Allele` recursion rather than meeting it at
//! the top level, and therefore the one a kind-matching arm can miss.

use ferro_hgvs::data::cdot::CdotMapper;
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::hgvs::variant::{AllelePhase, AlleleVariant, HgvsVariant};
use ferro_hgvs::parse_hgvs;
use ferro_hgvs::reference::mock::JsonProvider;
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::VariantProjector;

/// The bare ring, the `sup`-wrapped ring, and the ring-as-cis-allele-member.
///
/// The allele pairs a *legal* genomic member with the ring, so a decline can only
/// come from walking into the second member — the same construction #1587 uses
/// to prove an allele-recursion test is not passing for an unrelated reason.
fn ring_shapes() -> Vec<(&'static str, HgvsVariant)> {
    let ring = parse_hgvs("NC_000022.11:g.100_101del::200_201del").expect("ring parses");
    let sup =
        parse_hgvs("NC_000022.11:g.[100_101del::200_201del]sup").expect("sup-wrapped ring parses");
    let legal_member = parse_hgvs("NC_000022.11:g.500A>T").expect("plain genomic member parses");
    let ring_in_allele = HgvsVariant::Allele(AlleleVariant::new(
        vec![legal_member, ring.clone()],
        AllelePhase::Cis,
    ));
    vec![
        ("bare ring", ring),
        ("sup-wrapped ring", sup),
        ("ring as a cis-allele member", ring_in_allele),
    ]
}

/// The legal genomic member the ring shapes are paired with, used on its own as
/// each test's positive control.
fn legal_member() -> HgvsVariant {
    parse_hgvs("NC_000022.11:g.500A>T").expect("plain genomic member parses")
}

/// The allele-shaped positive control: two *legal* genomic members in cis. The
/// ring shapes below are also alleles, so without this the declines could come
/// from the allele wrapper rather than from the ring inside it.
///
/// Both members state `A` as their reference base, for the reason `provider()`
/// documents: the contig is all-`A`, and a member whose stated base disagrees
/// makes the control fail for a reason that has nothing to do with rings.
fn legal_allele() -> HgvsVariant {
    HgvsVariant::Allele(AlleleVariant::new(
        vec![
            parse_hgvs("NC_000022.11:g.500A>T").expect("first legal member parses"),
            parse_hgvs("NC_000022.11:g.600A>G").expect("second legal member parses"),
        ],
        AllelePhase::Cis,
    ))
}

/// A contig long enough to serve every position the ring shapes name, so a
/// decline is never merely "the provider had no bases".
///
/// All-`A` **on purpose**: it makes `g.500A>T`'s stated reference base agree with
/// the contig. An earlier revision used a repeating `ACGT` filler, which put `T`
/// at position 500 and made the control *fail* — with a message that also
/// contains "cannot apply to reference", so `spdi_apply_declines_every_ring_shape`
/// passed while proving nothing about rings. That is why each test below asserts
/// the **ring-specific** reason and runs the control first.
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NC_000022.11", "A".repeat(1000));
    provider
}

/// `apply_to_reference` must decline a ring rather than apply some subset of its
/// segments. This is the module whose decline is *indirect* — it falls out of
/// `hgvs_to_spdi` having no ring arm — so it is the one most exposed to a future
/// change that gives SPDI a ring representation without revisiting what applying
/// a ring to a linear reference would even mean.
#[test]
fn spdi_apply_declines_every_ring_shape() {
    let provider = provider();
    // Control: the legal member alone applies cleanly, so a decline below is
    // attributable to the ring rather than to the fixture.
    ferro_hgvs::spdi::apply::apply_to_reference(&legal_member(), &provider)
        .expect("control: the legal member alone must apply");

    for (label, variant) in ring_shapes() {
        let error = ferro_hgvs::spdi::apply::apply_to_reference(&variant, &provider)
            .expect_err(&format!("{label}: applying a ring must decline"));
        let rendered = error.to_string();
        // The *ring-specific* reason, not the generic prefix: `apply_to_reference`
        // has a second decline ("its members overlap, or a stated reference base
        // disagrees") that shares the prefix, and matching the prefix alone would
        // accept a decline that had nothing to do with the ring.
        assert!(
            rendered.contains("no single resulting sequence is defined"),
            "{label}: the decline must be the SPDI-cannot-represent-it one, got: {rendered}"
        );
    }
}

/// The VCF converter must decline a ring on both entry points. `convert_all` is
/// the decomposing one, and for the allele shape it must fail the whole call
/// naming the offending member rather than silently emitting one record for the
/// legal member and dropping the ring — a partial success is the failure mode
/// worth guarding, because it looks like a successful conversion.
#[test]
fn vcf_conversion_declines_every_ring_shape() {
    let provider = MockProvider::with_test_data();
    let transcript =
        ferro_hgvs::reference::ReferenceProvider::get_transcript(&provider, "NM_000088.3")
            .expect("the shared test fixture carries NM_000088.3");
    let converter = ferro_hgvs::vcf::HgvsToVcfConverter::new(&transcript, &provider);

    // Control: the legal member converts, and a two-member legal allele
    // decomposes to two records — so `convert_all` is genuinely walking members.
    converter
        .convert(&legal_member())
        .expect("control: the legal member alone must convert");
    let legal_allele = HgvsVariant::Allele(AlleleVariant::new(
        vec![
            legal_member(),
            parse_hgvs("NC_000022.11:g.600C>G").expect("second legal member parses"),
        ],
        AllelePhase::Cis,
    ));
    assert_eq!(
        converter
            .convert_all(&legal_allele)
            .expect("control: a legal allele must convert")
            .len(),
        2,
        "control: `convert_all` must emit one record per member",
    );

    for (label, variant) in ring_shapes() {
        converter
            .convert(&variant)
            .err()
            .unwrap_or_else(|| panic!("{label}: `convert` must decline a ring"));
        let error = converter
            .convert_all(&variant)
            .err()
            .unwrap_or_else(|| panic!("{label}: `convert_all` must decline a ring"));
        let rendered = error.to_string();
        assert!(
            rendered.contains("not supported"),
            "{label}: the decline must say the shape is unsupported, got: {rendered}"
        );
    }
}

/// `project_to_genomic` must decline a ring. An **empty** cdot is deliberate: a
/// ring is refused on variant kind before any transcript lookup, so a decline
/// that depended on cdot being populated would be evidence of the wrong thing.
#[test]
fn projection_declines_every_ring_shape() {
    let projector = VariantProjector::new(Projector::new(CdotMapper::new()), JsonProvider::new());
    // Control: a plain genomic member, and a legal genomic allele, both project
    // through this empty-cdot harness — so the declines below are about the ring.
    assert_eq!(
        projector
            .project_to_genomic(&legal_member())
            .expect("control: a plain genomic variant must project")
            .to_string(),
        "NC_000022.11:g.500A>T",
    );
    // The allele half of that control, which the comment above promised and the
    // code did not make. Without it the ring-in-allele decline below is not shown
    // to be about the *ring*: a projector that declined every `Allele` would pass
    // this test just as well.
    projector
        .project_to_genomic(&legal_allele())
        .expect("control: a cis allele of two legal genomic members must project");

    for (label, variant) in ring_shapes() {
        let error = projector
            .project_to_genomic(&variant)
            .expect_err(&format!("{label}: projecting a ring must decline"));
        let rendered = error.to_string();
        assert!(
            rendered.contains("ring") || rendered.contains("supernumerary"),
            "{label}: the decline must name the shape it refused, got: {rendered}"
        );
    }
}
