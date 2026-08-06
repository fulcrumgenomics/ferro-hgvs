//! Integration tests pinning LRG accession behavior across all spec-derived forms.
//!
//! The HGVS specification (v21.0, see `assets/hgvs-nomenclature/docs/background/refseq.md`)
//! defines three LRG accession forms:
//! - `LRG_<N>`        — the LRG genomic record itself (uses `g.`)
//! - `LRG_<N>t<M>`    — transcript M of LRG_N (uses `c.`, `n.`, or `r.`)
//! - `LRG_<N>p<M>`    — protein M of LRG_N (uses `p.`)
//!
//! `Accession::inferred_variant_type` returns the *primary* type for each
//! class — coding (`c`) for transcripts, protein (`p`) for proteins, genomic
//! (`g`) for the bare record. It is informational, not a validator: a
//! `LRG_199t1:n.5C>T` parses successfully even though the inferred primary
//! type is `c`.

use std::collections::HashMap;

use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::{parse_hgvs_lenient, parse_hgvs_with_config};
use ferro_hgvs::hgvs::variant::Accession;
use ferro_hgvs::parse_hgvs;
use rstest::rstest;

/// `LRG_<N>t<M>` returns "c" regardless of the variant block's coordinate type.
/// The function answers "what is the primary coordinate type for this accession?"
/// not "what coordinate types are valid with it?".
#[test]
fn lrg_transcript_inferred_primary_is_c_for_all_variant_types() {
    let acc = Accession::new("LRG", "199t1", None);
    assert_eq!(acc.inferred_variant_type(), Some("c"));

    // The variant block can pick any of c./n./r. — they all parse:
    for input in &[
        "LRG_199t1:c.357+1G>A",
        "LRG_199t1:n.5C>T",
        "LRG_199t1:r.426_427insa",
    ] {
        let v = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} failed: {e}"));
        assert_eq!(v.to_string(), *input, "round-trip mismatch for {input}");
    }
}

/// `LRG_<N>p<M>` always maps to "p".
#[test]
fn lrg_protein_inferred_is_p() {
    let acc = Accession::new("LRG", "199p1", None);
    assert_eq!(acc.inferred_variant_type(), Some("p"));
}

/// Bare `LRG_<N>` always maps to "g".
#[test]
fn lrg_bare_inferred_is_g() {
    let acc = Accession::new("LRG", "199", None);
    assert_eq!(acc.inferred_variant_type(), Some("g"));
}

/// `inferred_variant_type` is a pure function over `(prefix, number)` —
/// idempotent and reproducible across calls.
#[test]
fn lrg_inferred_type_is_idempotent() {
    let acc = Accession::new("LRG", "1t1", None);
    let first = acc.inferred_variant_type();
    let second = acc.inferred_variant_type();
    let third = acc.inferred_variant_type();
    assert_eq!(first, second);
    assert_eq!(second, third);
    assert_eq!(first, Some("c"));
}

/// Two accessions with identical `(prefix, number, version)` must hash
/// and compare equal regardless of construction path. Insert one via
/// `Accession::new` and look up via a parsed accession — both should hit
/// the same hash bucket.
#[test]
fn lrg_accessions_hash_and_eq_consistently() {
    let constructed = Accession::new("LRG", "1t1", None);
    let parsed = match parse_hgvs("LRG_1t1:c.100A>G").unwrap() {
        ferro_hgvs::HgvsVariant::Cds(c) => c.accession,
        other => panic!("expected Cds variant, got {other:?}"),
    };

    assert_eq!(
        constructed, parsed,
        "Eq must hold across construction paths"
    );

    let mut map: HashMap<Accession, &'static str> = HashMap::new();
    map.insert(constructed, "marker");
    assert_eq!(
        map.get(&parsed),
        Some(&"marker"),
        "Hash must match across construction paths"
    );
}

/// Malformed LRG numbers (no digits before/after the discriminator, alpha
/// after, etc.) fall back to "g" — the same conservative answer the
/// pre-#141 code returned for *any* LRG. This preserves backward
/// compatibility with inputs that ferro accepts leniently.
#[test]
fn lrg_malformed_number_falls_back_to_g() {
    // No digits before `t` (e.g., `LRG_t1`).
    assert_eq!(
        Accession::new("LRG", "t1", None).inferred_variant_type(),
        Some("g")
    );
    // No digits after `t`.
    assert_eq!(
        Accession::new("LRG", "1t", None).inferred_variant_type(),
        Some("g")
    );
    // Alpha discriminator suffix.
    assert_eq!(
        Accession::new("LRG", "1tabc", None).inferred_variant_type(),
        Some("g")
    );
    // Two discriminator letters in a row.
    assert_eq!(
        Accession::new("LRG", "1tt1", None).inferred_variant_type(),
        Some("g")
    );
    // Discriminator chained (`1t1p1`) — first non-digit is `t`, but the
    // span after contains a non-digit `p`, so falls back to `g`.
    assert_eq!(
        Accession::new("LRG", "1t1p1", None).inferred_variant_type(),
        Some("g")
    );
    // Uppercase discriminator (spec mandates lowercase). Lenient: returns g.
    assert_eq!(
        Accession::new("LRG", "1T1", None).inferred_variant_type(),
        Some("g")
    );
    assert_eq!(
        Accession::new("LRG", "1P1", None).inferred_variant_type(),
        Some("g")
    );
}

/// `LRG_1(NM_000088.3):c.459A>G` is a valid compound reference: LRG outer
/// (genomic context) + NM transcript inner. The parsed Accession's primary
/// `prefix` is `NM`, and `inferred_variant_type` returns `c`.
#[test]
fn lrg_outer_compound_with_nm_inner() {
    let v = parse_hgvs("LRG_1(NM_000088.3):c.459A>G").unwrap();
    assert_eq!(v.to_string(), "LRG_1(NM_000088.3):c.459A>G");

    let acc = match v {
        ferro_hgvs::HgvsVariant::Cds(c) => c.accession,
        other => panic!("expected Cds, got {other:?}"),
    };
    assert_eq!(&*acc.prefix, "NM");
    assert_eq!(acc.inferred_variant_type(), Some("c"));
    let ctx = acc
        .genomic_context
        .as_ref()
        .expect("must have genomic context");
    assert_eq!(&*ctx.prefix, "LRG");
    assert!(ctx.is_lrg());
}

/// LRG accessions never carry a version per HGVS spec (SVD-WG008). Ferro
/// parses `LRG_1.2:g.…` leniently but drops the version so Display emits
/// the canonical `LRG_1:g.…`.
#[test]
fn lrg_version_is_dropped_on_construction() {
    // Direct construction.
    let acc_with_ver = Accession::new("LRG", "1", Some(2));
    assert_eq!(
        acc_with_ver.version, None,
        "Accession::new must drop LRG version"
    );
    assert_eq!(acc_with_ver.full(), "LRG_1");

    let acc_styled = Accession::with_style("LRG", "1t1", Some(7), false);
    assert_eq!(acc_styled.version, None);
    assert_eq!(acc_styled.full(), "LRG_1t1");
}

/// Round-trip: parser may accept a stray `.N` after LRG, but Display drops it.
#[test]
fn lrg_version_round_trip_drops_illegal_version() {
    let v = parse_hgvs("LRG_1.2:g.100A>G").unwrap();
    assert_eq!(
        v.to_string(),
        "LRG_1:g.100A>G",
        "Display must canonicalize away the illegal LRG version"
    );
}

/// `NC_…(LRG_<N>t<M>)` compound reference: NC outer (chromosome) +
/// LRG transcript inner. Today this fails because the inner-accession
/// dispatcher in `parse_compound_inner` doesn't recognize the `LRG_`
/// 3-letter prefix; this test pins the fix.
#[test]
fn nc_outer_compound_with_lrg_transcript_inner() {
    let inputs_and_expected_inferred = [
        ("NC_000023.11(LRG_199t1):c.357+1G>A", "c"),
        ("NC_000023.11(LRG_199):g.1A>G", "g"),
        ("NC_000023.11(LRG_199p1):p.Trp24Cys", "p"),
    ];
    for (input, expected) in inputs_and_expected_inferred {
        let v = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} failed: {e}"));
        assert_eq!(v.to_string(), input, "round-trip for {input}");
        let acc = match v {
            ferro_hgvs::HgvsVariant::Cds(c) => c.accession,
            ferro_hgvs::HgvsVariant::Genome(g) => g.accession,
            ferro_hgvs::HgvsVariant::Protein(p) => p.accession,
            other => panic!("unexpected variant kind for {input}: {other:?}"),
        };
        assert!(acc.is_lrg(), "primary accession must be LRG for {input}");
        assert_eq!(
            acc.inferred_variant_type(),
            Some(expected),
            "inferred type mismatch for {input}",
        );
        let ctx = acc.genomic_context.as_ref().expect("must have NC context");
        assert_eq!(&*ctx.prefix, "NC");
    }
}

/// Cross-LRG compound allele in cis. Both inner accessions are LRG
/// transcripts and each must report `c` as the inferred primary type.
#[test]
fn cross_lrg_compound_allele_in_cis() {
    let input = "[LRG_1t1:c.1A>G;LRG_2t1:c.2A>G]";
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} failed: {e}"));
    assert_eq!(v.to_string(), input);
}

/// Same LRG transcript, two variants in cis.
#[test]
fn lrg_compound_allele_same_transcript() {
    let input = "LRG_199t1:c.[633A>G;635T>C]";
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} failed: {e}"));
    assert_eq!(v.to_string(), input);
}

/// LRG transcript heterozygous trans-allele.
#[test]
fn lrg_compound_allele_in_trans() {
    let input = "LRG_199t1:c.[76A>G];[103del]";
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} failed: {e}"));
    assert_eq!(v.to_string(), input);
}

/// Spec corpus round-trips. Examples lifted from
/// `assets/hgvs-nomenclature/docs/recommendations/` and
/// `assets/hgvs-nomenclature/docs/background/refseq.md`.
#[rstest]
// Bare LRG genomic
#[case("LRG_199:g.954966C>T", "g")]
#[case("LRG_476:g.4950_39800=", "g")]
// LRG transcript with c. (background/refseq.md, recommendations/DNA/delins.md)
#[case("LRG_199t1:c.357+1G>A", "c")]
#[case("LRG_199t1:c.4661delinsTC", "c")]
#[case("LRG_199t1:c.145_147delinsTGG", "c")]
// LRG transcript with r. (recommendations/RNA/insertion.md, deletion.md)
#[case("LRG_199t1:r.426_427insa", "c")] // primary inferred = c, real = r
#[case("LRG_2t1:r.1034_1036del", "c")]
// LRG protein (recommendations/protein/substitution.md)
#[case("LRG_199p1:p.Trp24Cys", "p")]
fn lrg_spec_corpus_round_trip(#[case] input: &str, #[case] expected_primary: &str) {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} failed: {e}"));
    assert_eq!(v.to_string(), input, "round-trip mismatch for {input}");

    let acc = v
        .accession()
        .unwrap_or_else(|| panic!("{input} should expose an accession"));
    assert!(acc.is_lrg(), "{input} should have an LRG accession");
    assert_eq!(
        acc.inferred_variant_type(),
        Some(expected_primary),
        "primary inferred type mismatch for {input}",
    );
}

/// Edge values that ARE valid digit-shaped — pinned as `c`/`p` because the
/// validator only checks "digits-only either side of the discriminator".
#[test]
fn lrg_zero_and_leading_zero_discriminators_are_accepted() {
    // `LRG_0` exists nowhere in the wild but ferro accepts it leniently.
    assert_eq!(
        Accession::new("LRG", "0", None).inferred_variant_type(),
        Some("g")
    );
    // Leading zero in M (`t01`): digits-only, accepted as `c`.
    assert_eq!(
        Accession::new("LRG", "1t01", None).inferred_variant_type(),
        Some("c")
    );
    // Massive numbers — digits-only check is independent of magnitude.
    assert_eq!(
        Accession::new("LRG", "999999t999999", None).inferred_variant_type(),
        Some("c")
    );
}

// ---------------------------------------------------------------------------
// Strict-mode acceptance of the versionless LRG spelling (#1498)
// ---------------------------------------------------------------------------

/// Strict mode must accept a bare `LRG_<N>` — it is the only legal spelling.
///
/// W3001 (`missing accession version`) is scoped by the spec to databases that
/// version at all: `background/refseq.md` L23 requires a version "only when the
/// reference sequence databases use versioning to distinguish between unique
/// sequences", naming RefSeq and Ensembl. LRG is listed separately (L66-69) as
/// `LRG_199` / `LRG_199t1` / `LRG_199p1`, none versioned.
///
/// Before #1498, strict mode — which is `ferro normalize`'s **default**
/// `--error-mode` — rejected every one of these with "Accession 'LRG_292' is
/// missing a version suffix", making LRG genomic variants unusable by default.
#[test]
fn strict_mode_accepts_versionless_lrg_accessions() {
    for input in [
        "LRG_292:g.160875del",
        "LRG_1:g.50dup",
        "LRG_199t1:c.79del",
        "LRG_199p1:p.Ser68Arg",
        "LRG_163t1:n.5C>T",
    ] {
        let parsed = parse_hgvs_with_config(input, ErrorConfig::strict());
        assert!(
            parsed.is_ok(),
            "strict mode must accept the versionless LRG form {input}: {:?}",
            parsed.err(),
        );
    }
}

/// The discriminating half: exempting LRG must not exempt RefSeq or Ensembl,
/// which the spec explicitly does require to be versioned.
#[test]
fn strict_mode_still_rejects_versionless_refseq_and_ensembl() {
    for (input, accession) in [
        ("NM_000088:c.100A>G", "NM_000088"),
        ("NG_012232:g.100del", "NG_012232"),
        ("NC_000023:g.100del", "NC_000023"),
        ("ENST00000380152:c.100A>G", "ENST00000380152"),
    ] {
        // Assert *which* error, not merely that there is one. A bare `is_err()`
        // is satisfied by any unrelated parse failure, and this is the control
        // that rules out the over-broad fix (marking every family
        // version-optional, or dropping the check) — so it has to fail for the
        // specific reason it names, or it stops discriminating.
        let err = match parse_hgvs_with_config(input, ErrorConfig::strict()) {
            Ok(v) => panic!("strict mode must still reject {input}, got {v:?}"),
            Err(e) => e.to_string(),
        };
        assert!(
            err.contains("missing a version suffix"),
            "{input} must be rejected for its missing version, got: {err}"
        );
        assert!(
            err.contains(accession),
            "the diagnostic for {input} must name `{accession}`, got: {err}"
        );
    }
}

/// Lenient mode emits no W3001 for a versionless LRG either — the warning was
/// as wrong as the rejection, it was just survivable.
#[test]
fn lenient_mode_emits_no_missing_version_warning_for_lrg() {
    let parsed = parse_hgvs_lenient("LRG_292:g.160875del").unwrap();
    assert!(
        !parsed
            .warnings
            .iter()
            .any(|w| w.error_type.code() == "W3001"),
        "got warnings: {:?}",
        parsed.warnings,
    );
    // Still fires for a genuinely versionless RefSeq accession.
    let refseq = parse_hgvs_lenient("NM_000088:c.100A>G").unwrap();
    assert!(refseq
        .warnings
        .iter()
        .any(|w| w.error_type.code() == "W3001"));
}
