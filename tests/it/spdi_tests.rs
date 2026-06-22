//! SPDI Roundtrip Tests
//!
//! These tests validate HGVS ↔ SPDI consistency using an offline, hand-authored
//! synthetic fixture (`tests/fixtures/validation/ncbi_variation.json`) rather than
//! a live NCBI Variation Services fetch — the fixture is committed and curated, so
//! the tests run deterministically without network access. SPDI (Sequence Position
//! Deletion Insertion) is NCBI's canonical representation for sequence variants.
//!
//! Tests verify that:
//! 1. HGVS expressions parse correctly
//! 2. When the fixture records an SPDI, we can parse the corresponding HGVS
//! 3. Roundtrip conversions maintain semantic equivalence

use ferro_hgvs::spdi::convert::{hgvs_to_spdi_simple, spdi_to_hgvs};
use ferro_hgvs::spdi::SpdiVariant;
use ferro_hgvs::{parse_hgvs, HgvsVariant};
use serde::Deserialize;
use std::fs;
use std::path::Path;

// ============================================================================
// NCBI Variation Services Fixture Types
// ============================================================================

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct NcbiVariationFixture {
    source: String,
    api_base: String,
    generated: String,
    hgvs_conversions: HgvsConversions,
    rsid_lookups: RsidLookups,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct HgvsConversions {
    total: usize,
    successful: usize,
    variants: Vec<HgvsConversion>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct HgvsConversion {
    input_hgvs: String,
    spdi_result: Option<SpdiResult>,
    roundtrip_hgvs: Option<RoundtripResult>,
    vcf: Option<serde_json::Value>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct SpdiResult {
    #[serde(default)]
    data: Option<SpdiData>,
    #[serde(default)]
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct SpdiData {
    #[serde(default)]
    spdis: Vec<Spdi>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct Spdi {
    seq_id: String,
    position: i64,
    deleted_sequence: String,
    inserted_sequence: String,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct RoundtripResult {
    #[serde(default)]
    data: Option<RoundtripData>,
    #[serde(default)]
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct RoundtripData {
    #[serde(default)]
    hgvs_list: Vec<String>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct RsidLookups {
    total: usize,
    successful: usize,
    variants: Vec<RsidLookup>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct RsidLookup {
    rsid: String,
    info: Option<serde_json::Value>,
    #[serde(default)]
    hgvs_expressions: Vec<String>,
    spdi: Option<serde_json::Value>,
}

// ============================================================================
// Test Report
// ============================================================================

#[derive(Default)]
struct SpdiTestReport {
    total_hgvs: usize,
    hgvs_parsed: usize,
    hgvs_with_spdi: usize,
    roundtrip_available: usize,
    roundtrip_parsed: usize,
    roundtrip_verified: usize,
    parse_failures: Vec<String>,
}

impl SpdiTestReport {
    fn summary(&self) -> String {
        format!(
            "HGVS total: {}, parsed: {}, with SPDI: {}, SPDI roundtrips verified: {}, roundtrip available: {}, roundtrip parsed: {}",
            self.total_hgvs,
            self.hgvs_parsed,
            self.hgvs_with_spdi,
            self.roundtrip_verified,
            self.roundtrip_available,
            self.roundtrip_parsed
        )
    }
}

// ============================================================================
// HGVS to SPDI Tests
// ============================================================================

#[test]
fn test_hgvs_to_spdi_parsing() {
    let fixture = load_fixture();

    assert!(
        !fixture.hgvs_conversions.variants.is_empty(),
        "fixture must contain HGVS conversions; otherwise this test passes vacuously"
    );

    // The fixture's aggregate header counts must agree with the actual data so
    // they cannot silently drift: `total` is the number of recorded conversions,
    // and `successful` is the number that carry an SPDI result.
    assert_eq!(
        fixture.hgvs_conversions.total,
        fixture.hgvs_conversions.variants.len(),
        "hgvs_conversions.total must equal the number of recorded conversions"
    );
    let hgvs_with_spdi_result = fixture
        .hgvs_conversions
        .variants
        .iter()
        .filter(|conversion| conversion.spdi_result.is_some())
        .count();
    assert_eq!(
        fixture.hgvs_conversions.successful, hgvs_with_spdi_result,
        "hgvs_conversions.successful must equal the number of conversions with an SPDI result"
    );

    let mut report = SpdiTestReport {
        total_hgvs: fixture.hgvs_conversions.variants.len(),
        ..Default::default()
    };

    for conversion in &fixture.hgvs_conversions.variants {
        // Every input HGVS in the fixture must parse.
        let parse_result = parse_hgvs(&conversion.input_hgvs);
        match parse_result {
            Ok(_) => report.hgvs_parsed += 1,
            Err(_) => report.parse_failures.push(conversion.input_hgvs.clone()),
        }

        // When the fixture carries an SPDI for this variant, verify the full
        // HGVS -> SPDI -> HGVS roundtrip is exact. We only assert SPDI equality
        // for variants whose HGVS -> SPDI conversion is computable offline (i.e.
        // without a reference provider): genomic substitutions. Such variants
        // need no reference bases, so `hgvs_to_spdi_simple` resolves them fully.
        // The asserted SPDI values therefore come from ferro's own converter,
        // never a hand-guessed relationship.
        if let Some(ref spdi_result) = conversion.spdi_result {
            if let Some(ref data) = spdi_result.data {
                if let Some(expected) = data.spdis.first() {
                    report.hgvs_with_spdi += 1;
                    verify_spdi_roundtrip(&conversion.input_hgvs, expected, &mut report);
                }
            }
        }

        // Any roundtrip HGVS strings the fixture records must parse AND be
        // semantically equal to the input variant — not merely syntactically
        // valid. For these genomic substitutions NCBI's SPDI->HGVS roundtrip
        // reproduces the input exactly, so we assert full variant equality
        // (catches silent position/allele drift, which an is_ok() check would
        // miss).
        if let Some(ref roundtrip) = conversion.roundtrip_hgvs {
            if let Some(ref data) = roundtrip.data {
                if !data.hgvs_list.is_empty() {
                    report.roundtrip_available += 1;
                    let input_variant = parse_hgvs(&conversion.input_hgvs).unwrap_or_else(|e| {
                        panic!("input HGVS should parse: {}: {e:?}", conversion.input_hgvs)
                    });
                    for hgvs in &data.hgvs_list {
                        let round = parse_hgvs(hgvs).unwrap_or_else(|e| {
                            panic!("roundtrip HGVS should parse: {hgvs}: {e:?}")
                        });
                        assert_eq!(
                            round, input_variant,
                            "roundtrip HGVS {hgvs} should be semantically equal to input {}",
                            conversion.input_hgvs
                        );
                    }
                    report.roundtrip_parsed += 1;
                }
            }
        }
    }

    println!("\n=== SPDI Roundtrip Test Report ===");
    println!("{}", report.summary());

    // All input HGVS in the synthetic fixture are valid and must parse.
    assert!(
        report.parse_failures.is_empty(),
        "all input HGVS should parse; failures: {:?}",
        report.parse_failures
    );
    assert_eq!(report.hgvs_parsed, report.total_hgvs);

    // The fixture carries genomic-substitution SPDIs; each must roundtrip
    // exactly. Guard against the assertions becoming a no-op.
    assert!(
        report.hgvs_with_spdi >= 4,
        "expected at least 4 variants with verifiable SPDI, found {}",
        report.hgvs_with_spdi
    );
    assert_eq!(
        report.roundtrip_verified, report.hgvs_with_spdi,
        "every variant with an SPDI should pass the exact HGVS<->SPDI roundtrip"
    );
}

/// Verify that an input genomic-substitution HGVS converts to exactly the
/// SPDI recorded in the fixture, and that converting that SPDI back yields an
/// HGVS expression that parses and re-expresses the original variant.
///
/// This proves the asserted HGVS<->SPDI relationship genuinely holds rather
/// than trusting a hand-authored value. Only genomic substitutions are
/// asserted here because their HGVS->SPDI conversion needs no reference
/// provider (`hgvs_to_spdi_simple`); transcript/`c.` variants are exercised
/// on the parse path only (their conversion requires CDS metadata / reference
/// bases that a unit fixture cannot supply offline).
fn verify_spdi_roundtrip(input_hgvs: &str, expected: &Spdi, report: &mut SpdiTestReport) {
    let variant = parse_hgvs(input_hgvs)
        .unwrap_or_else(|e| panic!("input HGVS should parse: {input_hgvs}: {e:?}"));

    let expected_pos = u64::try_from(expected.position)
        .unwrap_or_else(|_| panic!("SPDI position must be non-negative: {}", expected.position));

    // Independent oracle (not ferro): the fixture's SPDI must satisfy the SPDI
    // spec's own arithmetic, derived directly from the input HGVS string by a
    // standalone parse. For a genomic substitution `<acc>:g.<pos><ref>><alt>`
    // the canonical SPDI is `<acc>:<pos-1>:<ref>:<alt>` (0-based interbase
    // coordinate, deleted = ref base, inserted = alt base). This check does not
    // call any ferro converter, so passing it is not circular with the forward
    // conversion asserted below.
    let oracle = expected_spdi_from_genomic_substitution(input_hgvs);
    assert_eq!(
        (
            expected.seq_id.as_str(),
            expected_pos,
            expected.deleted_sequence.as_str(),
            expected.inserted_sequence.as_str(),
        ),
        (
            oracle.sequence.as_str(),
            oracle.position,
            oracle.deletion.as_str(),
            oracle.insertion.as_str(),
        ),
        "fixture SPDI for {input_hgvs} must match the independently derived SPDI"
    );

    // Forward: ferro's HGVS -> SPDI must equal the (independently validated)
    // fixture SPDI.
    let computed = hgvs_to_spdi_simple(&variant)
        .unwrap_or_else(|e| panic!("genomic HGVS should convert to SPDI: {input_hgvs}: {e:?}"));
    assert_eq!(
        computed.sequence, expected.seq_id,
        "SPDI seq_id for {input_hgvs}"
    );
    assert_eq!(
        computed.position, expected_pos,
        "SPDI position for {input_hgvs}"
    );
    assert_eq!(
        computed.deletion, expected.deleted_sequence,
        "SPDI del for {input_hgvs}"
    );
    assert_eq!(
        computed.insertion, expected.inserted_sequence,
        "SPDI ins for {input_hgvs}"
    );

    // Backward: rebuild the SPDI from the fixture fields, convert to HGVS, and
    // assert the result is semantically equal to the original variant (full
    // variant equality, not just a string/parse check).
    let spdi = SpdiVariant::new(
        &expected.seq_id,
        expected_pos,
        &expected.deleted_sequence,
        &expected.inserted_sequence,
    );
    let back = spdi_to_hgvs(&spdi)
        .unwrap_or_else(|e| panic!("SPDI should convert back to HGVS: {spdi}: {e:?}"));
    assert_eq!(
        back, variant,
        "SPDI {spdi} should round-trip to the original variant {input_hgvs}"
    );

    report.roundtrip_verified += 1;
}

/// Independently derive the canonical SPDI for a genomic-substitution HGVS
/// string of the form `<accession>:g.<pos><ref>><alt>`, using only the SPDI
/// specification's coordinate arithmetic (0-based interbase position =
/// 1-based HGVS position − 1; deleted = ref base; inserted = alt base).
///
/// This deliberately does not call any ferro conversion routine so it can
/// serve as an external oracle for the fixture's recorded SPDI values, rather
/// than re-deriving them from the same code under test.
fn expected_spdi_from_genomic_substitution(input_hgvs: &str) -> SpdiVariant {
    let (accession, rest) = input_hgvs
        .split_once(":g.")
        .unwrap_or_else(|| panic!("expected a genomic (g.) substitution: {input_hgvs}"));

    // rest looks like "<pos><ref>><alt>", e.g. "7674220C>T".
    let (lhs, alt) = rest
        .split_once('>')
        .unwrap_or_else(|| panic!("expected a substitution (`>`): {input_hgvs}"));
    let ref_base = lhs
        .chars()
        .next_back()
        .unwrap_or_else(|| panic!("missing reference base: {input_hgvs}"));
    let pos_str = &lhs[..lhs.len() - ref_base.len_utf8()];
    let one_based: u64 = pos_str
        .parse()
        .unwrap_or_else(|_| panic!("invalid position {pos_str:?} in {input_hgvs}"));
    assert!(
        one_based >= 1,
        "1-based position must be >= 1: {input_hgvs}"
    );

    SpdiVariant::substitution(
        accession,
        one_based - 1,
        ref_base.to_string(),
        alt.to_string(),
    )
}

/// Load the synthetic, committed NCBI variation fixture. The fixture is a
/// checked-in test asset (not a network fetch), so its absence is a real test
/// failure rather than a reason to skip.
fn load_fixture() -> NcbiVariationFixture {
    let fixture_path = Path::new("tests/fixtures/validation/ncbi_variation.json");
    let content = fs::read_to_string(fixture_path).unwrap_or_else(|e| {
        panic!(
            "failed to read committed fixture {}: {e}",
            fixture_path.display()
        )
    });
    serde_json::from_str(&content).expect("failed to parse ncbi_variation.json")
}

// ============================================================================
// rsID Lookup Tests
// ============================================================================

#[test]
fn test_rsid_hgvs_parsing() {
    let fixture = load_fixture();

    assert!(
        !fixture.rsid_lookups.variants.is_empty(),
        "fixture must contain rsID lookups; otherwise this test passes vacuously"
    );

    // The fixture's aggregate header counts must agree with the actual data so
    // they cannot silently drift: `total` is the number of recorded lookups, and
    // `successful` is the number that resolved to at least one HGVS expression.
    assert_eq!(
        fixture.rsid_lookups.total,
        fixture.rsid_lookups.variants.len(),
        "rsid_lookups.total must equal the number of recorded lookups"
    );
    let rsids_resolved = fixture
        .rsid_lookups
        .variants
        .iter()
        .filter(|lookup| !lookup.hgvs_expressions.is_empty())
        .count();
    assert_eq!(
        fixture.rsid_lookups.successful, rsids_resolved,
        "rsid_lookups.successful must equal the number of lookups with HGVS expressions"
    );

    let mut total_hgvs = 0;
    let mut rsids_with_hgvs = 0;

    for lookup in &fixture.rsid_lookups.variants {
        if !lookup.hgvs_expressions.is_empty() {
            rsids_with_hgvs += 1;

            for hgvs in &lookup.hgvs_expressions {
                total_hgvs += 1;
                // Every HGVS expression associated with an rsID must parse, and
                // re-parsing the variant's own rendering must yield an equal
                // variant (display/parse idempotency). This verifies semantic
                // stability, not merely that the string is syntactically valid.
                let variant = parse_hgvs(hgvs).unwrap_or_else(|e| {
                    panic!("HGVS for {} should parse: {hgvs}: {e:?}", lookup.rsid)
                });
                // Independent oracle (a different internal path than the
                // render->reparse idempotency check below): the parser must store
                // the accession exactly as written. The accession is the literal
                // substring before the first ':'; comparing it against the parsed
                // variant's own `accession()` field exercises field extraction, not
                // the Display+parse cycle, so the two assertions cannot mask a
                // shared parser/renderer defect.
                let expected_accession = hgvs
                    .split_once(':')
                    .map(|(accession, _)| accession)
                    .unwrap_or_else(|| panic!("HGVS for {} must contain ':': {hgvs}", lookup.rsid));
                let parsed_accession = variant
                    .accession()
                    .unwrap_or_else(|| {
                        panic!(
                            "variant for {} should carry an accession: {hgvs}",
                            lookup.rsid
                        )
                    })
                    .to_string();
                assert_eq!(
                    parsed_accession, expected_accession,
                    "parsed accession for {} should match the input literal: {hgvs}",
                    lookup.rsid
                );
                let reparsed = parse_hgvs(&variant.to_string()).unwrap_or_else(|e| {
                    panic!(
                        "rendered HGVS for {} should re-parse: {variant}: {e:?}",
                        lookup.rsid
                    )
                });
                assert_eq!(
                    reparsed, variant,
                    "HGVS for {} should be display/parse idempotent: {hgvs}",
                    lookup.rsid
                );
            }
        }
    }

    println!("\n=== rsID HGVS Parsing Report ===");
    println!(
        "rsIDs with HGVS: {}/{}",
        rsids_with_hgvs,
        fixture.rsid_lookups.variants.len()
    );
    println!("Total HGVS expressions: {total_hgvs}");

    // Guard against the per-expression assertions never running.
    assert!(
        rsids_with_hgvs >= 1,
        "at least one rsID should carry HGVS expressions"
    );
    assert!(
        total_hgvs > 0,
        "expected at least one HGVS expression to verify"
    );
}

// ============================================================================
// SPDI Format Validation Tests
// ============================================================================

#[test]
fn test_spdi_format_understanding() {
    // Test that we understand SPDI format correctly
    // SPDI: Sequence:Position:Deletion:Insertion

    // Example: NM_000518.4:76:G:GG represents a G duplication
    // This corresponds to NM_000518.4:c.27dupG

    // Verify we can parse the HGVS form
    let result = parse_hgvs("NM_000518.4:c.27dupG");
    assert!(
        result.is_ok(),
        "Should parse canonical duplication notation"
    );

    // Verify we can parse the equivalent insertion form
    let result2 = parse_hgvs("NM_000518.4:c.27_28insG");
    assert!(result2.is_ok(), "Should parse insertion notation");
}

#[test]
fn test_spdi_deletion_equivalence() {
    // SPDI deletion: NM_000492.4:1653:CTT: (empty insertion = deletion)
    // HGVS: NM_000492.4:c.1521_1523del

    let result = parse_hgvs("NM_000492.4:c.1521_1523del");
    assert!(result.is_ok(), "Should parse deletion");

    // With explicit deleted sequence
    let result2 = parse_hgvs("NM_000492.4:c.1521_1523delCTT");
    assert!(result2.is_ok(), "Should parse deletion with sequence");
}

#[test]
fn test_spdi_substitution_equivalence() {
    // SPDI substitution: NC_000017.11:7674220:C:T
    // HGVS: NC_000017.11:g.7674220C>T

    let result = parse_hgvs("NC_000017.11:g.7674220C>T");
    assert!(result.is_ok(), "Should parse genomic substitution");

    if let Ok(variant) = result {
        match variant {
            HgvsVariant::Genome(genome) => {
                assert_eq!(format!("{}", genome.accession), "NC_000017.11");
            }
            _ => panic!("Expected genome variant"),
        }
    }
}

#[test]
fn test_spdi_indel_equivalence() {
    // SPDI complex: NM_000546.6:100:ACG:TTT
    // HGVS: NM_000546.6:c.100_102delinsTTT

    let result = parse_hgvs("NM_000546.6:c.100_102delinsTTT");
    assert!(result.is_ok(), "Should parse delins");

    if let Ok(variant) = result {
        match variant {
            HgvsVariant::Cds(_) => {} // Expected - delins is a coding variant
            _ => panic!("Expected Cds variant for delins"),
        }
    }
}

// ============================================================================
// Coordinate System Tests
// ============================================================================

#[test]
fn test_spdi_zero_based_understanding() {
    // SPDI uses 0-based coordinates
    // HGVS uses 1-based coordinates

    // SPDI: NM_000546.6:214:C:G (0-based position 214)
    // HGVS: NM_000546.6:c.215C>G (1-based position 215)

    let result = parse_hgvs("NM_000546.6:c.215C>G");
    assert!(result.is_ok(), "Should parse with 1-based coordinate");
}

#[test]
fn test_spdi_half_open_intervals() {
    // SPDI uses half-open intervals [start, end)
    // HGVS uses closed intervals [start, end]

    // For a 3bp deletion:
    // SPDI: position=1520, deleted_sequence=CTT (deletes positions 1520,1521,1522)
    // HGVS: c.1521_1523del (deletes positions 1521,1522,1523 in 1-based)

    let result = parse_hgvs("NM_000492.4:c.1521_1523del");
    assert!(result.is_ok());
}
