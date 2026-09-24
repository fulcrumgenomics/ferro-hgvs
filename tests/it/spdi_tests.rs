//! SPDI Roundtrip Tests
//!
//! These tests validate HGVS ↔ SPDI consistency using an offline, hand-authored
//! synthetic fixture (`tests/fixtures/validation/ncbi_variation.json`) rather than
//! a live NCBI Variation Services fetch — the fixture is committed and curated, so
//! the tests run deterministically without network access. SPDI (Sequence Position
//! Deletion Insertion) is NCBI's canonical representation for sequence variants.
//!
//! # Live data (#2256)
//!
//! The weekly External API Validation workflow deliberately writes a live NCBI
//! fetch over the fixture and runs these tests against it, with
//! `FERRO_SPDI_LIVE_FIXTURE=1` set. That is the only place ferro's genomic
//! HGVS -> SPDI conversion is compared with NCBI's own answers. Everywhere else
//! (PR CI, a local run) the committed curated fixture is read. [`FixtureKind`]
//! requires the marker and the fixture's `source` to agree, so a live fetch
//! committed over the fixture fails PR CI, and a weekly run that silently fell
//! back to the curated file fails too.
//!
//! What differs on live data, and how the tests handle it:
//!
//! - NCBI returns SPDIs for transcript (`c.`) inputs as well, in transcript-
//!   sequence coordinates (e.g. `NM_000546.6:356:C:G`). The round-trip checks
//!   apply only to single-base genomic substitutions, the one shape the offline
//!   oracle and converter handle, so those SPDIs are counted and skipped, not
//!   verified.
//! - NCBI's SPDI -> HGVS endpoint returns one `hgvs` string. The curated fixture
//!   uses the same shape.
//! - NCBI can be partly down (22 of 25 conversions returned HTTP 502 on
//!   2026-08-30). A genomic input may lack an SPDI only when its error is
//!   transient (HTTP 5xx, HTTP 429, or a network exception). A permanent error
//!   such as HTTP 400, which is what a wrong reference base in the input gets,
//!   fails the test. At least one genomic SPDI must be verified, so a total
//!   outage fails rather than passing having checked nothing.
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
    /// The HGVS the SPDI converts back to. NCBI's `/spdi/{spdi}/hgvs` endpoint
    /// returns a single string, and the curated fixture uses the same shape.
    #[serde(default)]
    hgvs: Option<String>,
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
    genomic_substitutions: usize,
    hgvs_with_spdi: usize,
    roundtrip_parsed: usize,
    roundtrip_verified: usize,
    transient_failures: usize,
    spdi_skipped_not_genomic_substitution: usize,
    parse_failures: Vec<String>,
}

impl SpdiTestReport {
    fn summary(&self) -> String {
        format!(
            "HGVS total: {}, parsed: {}, genomic substitutions: {}, with SPDI: {}, SPDI roundtrips verified: {}, roundtrip HGVS checked: {}, transient fetch failures: {}, SPDIs skipped (not a genomic substitution): {}",
            self.total_hgvs,
            self.hgvs_parsed,
            self.genomic_substitutions,
            self.hgvs_with_spdi,
            self.roundtrip_verified,
            self.roundtrip_parsed,
            self.transient_failures,
            self.spdi_skipped_not_genomic_substitution
        )
    }
}

// ============================================================================
// HGVS to SPDI Tests
// ============================================================================

#[test]
fn test_hgvs_to_spdi_parsing() {
    let fixture = load_fixture();
    let expect_live = std::env::var(LIVE_FIXTURE_ENV).is_ok_and(|value| value == "1");
    let kind = FixtureKind::from_source(&fixture.source, expect_live);

    assert!(
        !fixture.hgvs_conversions.variants.is_empty(),
        "fixture must contain HGVS conversions; otherwise this test passes vacuously"
    );

    // The fixture's aggregate header counts must agree with the actual data so
    // they cannot silently drift: `total` is the number of recorded conversions,
    // and `successful` is the number that carry a *successful* SPDI result (an
    // `spdi_result` that is present and reports no error).
    assert_eq!(
        fixture.hgvs_conversions.total,
        fixture.hgvs_conversions.variants.len(),
        "hgvs_conversions.total must equal the number of recorded conversions"
    );
    let hgvs_with_spdi_result = count_successful_conversions(&fixture.hgvs_conversions.variants);
    assert_eq!(
        fixture.hgvs_conversions.successful, hgvs_with_spdi_result,
        "hgvs_conversions.successful must equal the number of conversions with a successful SPDI result"
    );

    println!("\n=== SPDI Roundtrip Test Report ===");
    println!(
        "fixture: {:?} ({kind:?}), generated {}",
        fixture.source, fixture.generated
    );
    let report = check_hgvs_conversions(&fixture.hgvs_conversions.variants, kind);
    println!("{}", report.summary());

    // The curated fixture must also exercise the live shape the skip handles: a
    // conversion that carries an SPDI but is not a genomic substitution.
    if kind == FixtureKind::Curated {
        assert!(
            report.spdi_skipped_not_genomic_substitution >= 1,
            "the curated fixture must include a conversion that carries an SPDI but is not a genomic substitution"
        );
    }
}

/// Which NCBI fixture a run is reading. See the module docs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FixtureKind {
    /// The committed, hand-authored fixture.
    Curated,
    /// The live fetch the weekly workflow writes over it.
    Live,
}

impl FixtureKind {
    /// Decide the kind from the fixture's `source` and whether the run expects
    /// live data (`FERRO_SPDI_LIVE_FIXTURE=1`). The two must agree: a live
    /// source without the marker is a live fetch committed over the curated
    /// fixture, and the marker without a live source is a weekly run whose fetch
    /// never replaced the file.
    fn from_source(source: &str, expect_live: bool) -> Self {
        let is_live = source == LIVE_FETCH_SOURCE;
        match (expect_live, is_live) {
            (true, true) => FixtureKind::Live,
            (false, false) => FixtureKind::Curated,
            (true, false) => panic!(
                "{LIVE_FIXTURE_ENV}=1 but the fixture's source is {source:?}: the live NCBI \
                 fetch did not replace the committed fixture, so live data is not being tested"
            ),
            (false, true) => panic!(
                "the fixture's source is the live fetch's ({LIVE_FETCH_SOURCE:?}) but \
                 {LIVE_FIXTURE_ENV} is unset: a live fetch must never be committed over the \
                 curated fixture"
            ),
        }
    }
}

/// Whether a fetch error recorded by `scripts/fetch_ncbi_variation.py` may be
/// transient: an HTTP 5xx or 429, or a network exception (anything that is not
/// an `HTTP <status>` string). Any other HTTP status, typically a 400 for an
/// input NCBI rejects, is permanent.
fn is_transient_fetch_error(error: &str) -> bool {
    match error
        .strip_prefix("HTTP ")
        .and_then(|status| status.parse::<u16>().ok())
    {
        Some(status) => status == 429 || (500..600).contains(&status),
        None => true,
    }
}

/// Check every conversion and return the report. Per genomic substitution:
///
/// - with an SPDI, the HGVS -> SPDI -> HGVS round trip must hold exactly, and
///   any SPDI -> HGVS string recorded must denote the input variant;
/// - on the curated fixture, the SPDI and the SPDI -> HGVS string are required;
/// - on live data, either may be missing only for a transient fetch error.
///
/// Conversions that are not genomic substitutions are skipped (counted when
/// they carry an SPDI). At least one genomic SPDI must be verified.
fn check_hgvs_conversions(conversions: &[HgvsConversion], kind: FixtureKind) -> SpdiTestReport {
    let mut report = SpdiTestReport {
        total_hgvs: conversions.len(),
        ..Default::default()
    };

    for conversion in conversions {
        let input_hgvs = conversion.input_hgvs.as_str();
        // Every input HGVS in the fixture must parse.
        match parse_hgvs(input_hgvs) {
            Ok(_) => report.hgvs_parsed += 1,
            Err(_) => report.parse_failures.push(input_hgvs.to_string()),
        }

        // Only genomic substitutions convert to SPDI without a reference, so
        // only they are round-tripped; the oracle decides the shape.
        let Some(oracle) = expected_spdi_from_genomic_substitution(input_hgvs) else {
            let carries_spdi = conversion
                .spdi_result
                .as_ref()
                .and_then(|result| result.data.as_ref())
                .is_some_and(|data| !data.spdis.is_empty());
            if carries_spdi {
                report.spdi_skipped_not_genomic_substitution += 1;
            }
            continue;
        };
        report.genomic_substitutions += 1;

        let spdi_result = conversion.spdi_result.as_ref();
        let expected = spdi_result
            .and_then(|result| result.data.as_ref())
            .and_then(|data| data.spdis.first());
        let Some(expected) = expected else {
            require_transient_failure(
                kind,
                input_hgvs,
                "SPDI",
                spdi_result.and_then(|result| result.error.as_deref()),
                &mut report,
            );
            continue;
        };
        report.hgvs_with_spdi += 1;
        verify_spdi_roundtrip(input_hgvs, expected, &oracle, &mut report);

        // The SPDI -> HGVS string NCBI recorded must parse AND denote the input
        // variant, not merely be syntactically valid: full variant equality
        // catches silent position or allele drift.
        let roundtrip = conversion.roundtrip_hgvs.as_ref();
        let Some(hgvs) = roundtrip
            .and_then(|result| result.data.as_ref())
            .and_then(|data| data.hgvs.as_ref())
        else {
            require_transient_failure(
                kind,
                input_hgvs,
                "SPDI -> HGVS roundtrip",
                roundtrip.and_then(|result| result.error.as_deref()),
                &mut report,
            );
            continue;
        };
        let input_variant = parse_hgvs(input_hgvs)
            .unwrap_or_else(|e| panic!("input HGVS should parse: {input_hgvs}: {e:?}"));
        let round = parse_hgvs(hgvs)
            .unwrap_or_else(|e| panic!("roundtrip HGVS should parse: {hgvs}: {e:?}"));
        assert_eq!(
            round, input_variant,
            "roundtrip HGVS {hgvs} should be semantically equal to input {input_hgvs}"
        );
        report.roundtrip_parsed += 1;
    }

    assert!(
        report.parse_failures.is_empty(),
        "all input HGVS should parse; failures: {:?}",
        report.parse_failures
    );
    assert!(
        report.genomic_substitutions >= 1,
        "the fixture holds no genomic substitution, so nothing can be round-tripped; \
         a live fetch run with a small `--limit` requests none"
    );
    assert!(
        report.hgvs_with_spdi >= 1,
        "none of the {} genomic substitutions has an SPDI, so nothing was verified \
         (a total upstream outage on a live fetch)",
        report.genomic_substitutions
    );
    report
}

/// A genomic substitution lacks a result: fail unless this is live data and the
/// recorded error may be transient. A missing record (`null`) is never
/// transient; the fetch script records an attempted call as data or an error.
fn require_transient_failure(
    kind: FixtureKind,
    input_hgvs: &str,
    what: &str,
    error: Option<&str>,
    report: &mut SpdiTestReport,
) {
    match (kind, error) {
        (FixtureKind::Live, Some(error)) if is_transient_fetch_error(error) => {
            report.transient_failures += 1;
        }
        _ => panic!(
            "{input_hgvs}: no {what} in the {kind:?} fixture (error: {error:?}); \
             only a live fetch may lack one, and only for a transient error"
        ),
    }
}

/// Regression guard for issue #2017: a failed HGVS->SPDI conversion is recorded
/// by `scripts/fetch_ncbi_variation.py` as a *non-null* `spdi_result` carrying
/// an `error` field (the API's error body is stored inline), not as null. The
/// `successful` header counts only genuine successes, so the fixture
/// self-consistency check must too. Counting `spdi_result.is_some()` treats the
/// errored row as successful and makes the assertion fail on any nightly whose
/// live fetch hits an upstream API error — exactly the reported failure
/// (`successful` = 19 vs 25 `spdi_result` present). This pins the fetch harness'
/// two failure representations (null and error-object) against the counter.
#[test]
fn errored_conversions_are_excluded_from_the_successful_count() {
    // Shape emitted by a live fetch: one genuine success, one API error stored
    // as a non-null `spdi_result` (the case that broke the nightly), and one
    // failure stored as null (the committed fixture's convention). Only the
    // first is a success, so `successful` is 1.
    let json = r#"{
        "total": 3,
        "successful": 1,
        "variants": [
            {
                "input_hgvs": "NC_000017.11:g.7674220C>T",
                "spdi_result": {"data": {"spdis": [
                    {"seq_id": "NC_000017.11", "position": 7674219,
                     "deleted_sequence": "C", "inserted_sequence": "T"}
                ]}}
            },
            {
                "input_hgvs": "NM_000546.6:c.215C>G",
                "spdi_result": {"error": "HTTP 500", "input": "NM_000546.6:c.215C>G"}
            },
            {
                "input_hgvs": "NM_000492.4:c.1521_1523del",
                "spdi_result": null
            }
        ]
    }"#;
    let conversions: HgvsConversions =
        serde_json::from_str(json).expect("regression fixture JSON should deserialize");

    // Both the errored row and the null row are excluded; only the error-free
    // result counts, matching the `successful` header.
    assert_eq!(count_successful_conversions(&conversions.variants), 1);
    assert_eq!(
        count_successful_conversions(&conversions.variants),
        conversions.successful,
        "successful header must equal the count of error-free spdi_result objects"
    );

    // Guard against the guard passing vacuously: the errored row really is
    // present and non-null, so a naive is_some() count would (wrongly) be 2.
    let non_null = conversions
        .variants
        .iter()
        .filter(|c| c.spdi_result.is_some())
        .count();
    assert_eq!(
        non_null, 2,
        "the errored conversion must be stored as a non-null spdi_result"
    );
}

/// Parse an inline `hgvs_conversions.variants` array for the checker tests.
fn conversions_from(json: &str) -> Vec<HgvsConversion> {
    serde_json::from_str(json).expect("checker test JSON should deserialize")
}

/// One genomic substitution NCBI answered, in the live shape.
const GENOMIC_SUCCESS: &str = r#"{
    "input_hgvs": "NC_000017.11:g.7674220C>T",
    "spdi_result": {"data": {"spdis": [{"seq_id": "NC_000017.11", "position": 7674219,
        "deleted_sequence": "C", "inserted_sequence": "T"}]}},
    "roundtrip_hgvs": {"data": {"hgvs": "NC_000017.11:g.7674220C>T"}}
}"#;

#[test]
fn a_live_partial_outage_still_verifies_what_came_back() {
    let conversions = conversions_from(&format!(
        r#"[{GENOMIC_SUCCESS},
            {{"input_hgvs": "NC_000007.14:g.140753336A>T",
              "spdi_result": {{"error": "HTTP 502", "input": "NC_000007.14:g.140753336A>T"}}}},
            {{"input_hgvs": "NC_000013.11:g.32316461A>T",
              "spdi_result": {{"error": "HTTPSConnectionPool: Read timed out."}}}},
            {{"input_hgvs": "NM_000546.6:c.215C>G",
              "spdi_result": {{"data": {{"spdis": [{{"seq_id": "NM_000546.6", "position": 356,
                  "deleted_sequence": "C", "inserted_sequence": "G"}}]}}}},
              "roundtrip_hgvs": {{"data": {{"hgvs": "NM_000546.6:c.215C>G"}}}}}},
            {{"input_hgvs": "NM_001126112.3:c.35G>T",
              "spdi_result": {{"error": "HTTP 400"}}}}]"#
    ));
    let report = check_hgvs_conversions(&conversions, FixtureKind::Live);
    assert_eq!(report.roundtrip_verified, 1);
    assert_eq!(report.roundtrip_parsed, 1);
    assert_eq!(report.transient_failures, 2);
    assert_eq!(report.spdi_skipped_not_genomic_substitution, 1);
}

#[test]
#[should_panic(expected = "only a live fetch may lack one, and only for a transient error")]
fn a_permanent_rejection_of_a_genomic_input_fails_on_live_data() {
    // A wrong reference base gets HTTP 400; that must not pass as an outage.
    let conversions = conversions_from(&format!(
        r#"[{GENOMIC_SUCCESS},
            {{"input_hgvs": "NC_000017.11:g.7673802G>A",
              "spdi_result": {{"error": "HTTP 400"}}}}]"#
    ));
    check_hgvs_conversions(&conversions, FixtureKind::Live);
}

#[test]
#[should_panic(expected = "so nothing was verified")]
fn a_total_live_outage_fails() {
    let conversions = conversions_from(
        r#"[{"input_hgvs": "NC_000017.11:g.7674220C>T",
             "spdi_result": {"error": "HTTP 502"}},
            {"input_hgvs": "NC_000007.14:g.140753336A>T",
             "spdi_result": {"error": "HTTP 503"}}]"#,
    );
    check_hgvs_conversions(&conversions, FixtureKind::Live);
}

#[test]
#[should_panic(expected = "holds no genomic substitution")]
fn a_fetch_that_requested_no_genomic_substitution_fails_with_its_own_message() {
    let conversions = conversions_from(
        r#"[{"input_hgvs": "NM_000546.6:c.215C>G",
             "spdi_result": {"data": {"spdis": [{"seq_id": "NM_000546.6", "position": 356,
                 "deleted_sequence": "C", "inserted_sequence": "G"}]}}}]"#,
    );
    check_hgvs_conversions(&conversions, FixtureKind::Live);
}

#[test]
#[should_panic(expected = "only a live fetch may lack one")]
fn the_curated_fixture_may_not_lack_a_genomic_spdi_even_for_a_transient_error() {
    let conversions = conversions_from(&format!(
        r#"[{GENOMIC_SUCCESS},
            {{"input_hgvs": "NC_000007.14:g.140753336A>T",
              "spdi_result": {{"error": "HTTP 502"}}}}]"#
    ));
    check_hgvs_conversions(&conversions, FixtureKind::Curated);
}

#[test]
#[should_panic(expected = "no SPDI -> HGVS roundtrip")]
fn the_curated_fixture_may_not_lack_a_roundtrip_hgvs() {
    let conversions = conversions_from(
        r#"[{"input_hgvs": "NC_000017.11:g.7674220C>T",
             "spdi_result": {"data": {"spdis": [{"seq_id": "NC_000017.11", "position": 7674219,
                 "deleted_sequence": "C", "inserted_sequence": "T"}]}}}]"#,
    );
    check_hgvs_conversions(&conversions, FixtureKind::Curated);
}

#[test]
fn only_http_5xx_429_and_network_errors_are_transient() {
    for transient in [
        "HTTP 500",
        "HTTP 502",
        "HTTP 503",
        "HTTP 429",
        "Read timed out.",
    ] {
        assert!(is_transient_fetch_error(transient), "{transient}");
    }
    for permanent in ["HTTP 400", "HTTP 404", "HTTP 422"] {
        assert!(!is_transient_fetch_error(permanent), "{permanent}");
    }
}

#[test]
fn the_fixture_kind_requires_the_marker_and_the_source_to_agree() {
    assert_eq!(
        FixtureKind::from_source(LIVE_FETCH_SOURCE, true),
        FixtureKind::Live
    );
    assert_eq!(
        FixtureKind::from_source("synthetic (hand-authored)", false),
        FixtureKind::Curated
    );
    let committed_live = std::panic::catch_unwind(|| {
        FixtureKind::from_source(LIVE_FETCH_SOURCE, false);
    });
    assert!(committed_live.is_err(), "a committed live fetch must fail");
    let silent_fallback = std::panic::catch_unwind(|| {
        FixtureKind::from_source("synthetic (hand-authored)", true);
    });
    assert!(
        silent_fallback.is_err(),
        "a weekly run on curated data must fail"
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
fn verify_spdi_roundtrip(
    input_hgvs: &str,
    expected: &Spdi,
    oracle: &SpdiVariant,
    report: &mut SpdiTestReport,
) {
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
/// Returns `None` for any other shape: another axis, a gene selector, a
/// non-substitution, a multi-base or non-`ACGT` allele, or a position that is
/// not a plain number of at least 1. This is also the definition of which
/// conversions the round-trip checks cover.
///
/// This deliberately does not call any ferro conversion routine so it can
/// serve as an external oracle for the fixture's recorded SPDI values, rather
/// than re-deriving them from the same code under test.
fn expected_spdi_from_genomic_substitution(input_hgvs: &str) -> Option<SpdiVariant> {
    let (accession, rest) = input_hgvs.split_once(":g.")?;
    if accession.is_empty() || accession.contains('(') {
        return None;
    }
    // rest looks like "<pos><ref>><alt>", e.g. "7674220C>T".
    let (lhs, alt) = rest.split_once('>')?;
    let is_base = |c: &char| matches!(c, 'A' | 'C' | 'G' | 'T');
    let mut alt_chars = alt.chars();
    let alt_base = alt_chars.next().filter(is_base)?;
    if alt_chars.next().is_some() {
        return None;
    }
    let ref_base = lhs.chars().next_back().filter(is_base)?;
    let pos_str = &lhs[..lhs.len() - ref_base.len_utf8()];
    if pos_str.is_empty() || !pos_str.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let one_based: u64 = pos_str.parse().ok()?;
    let zero_based = one_based.checked_sub(1)?;
    Some(SpdiVariant::substitution(
        accession,
        zero_based,
        ref_base.to_string(),
        alt_base.to_string(),
    ))
}

#[test]
fn only_single_base_genomic_substitutions_reach_the_roundtrip_oracle() {
    assert_eq!(
        expected_spdi_from_genomic_substitution("NC_000017.11:g.7674220C>T"),
        Some(SpdiVariant::substitution(
            "NC_000017.11",
            7_674_219,
            "C",
            "T"
        ))
    );
    // The lowest position maps to interbase 0.
    assert_eq!(
        expected_spdi_from_genomic_substitution("NC_000017.11:g.1C>T"),
        Some(SpdiVariant::substitution("NC_000017.11", 0, "C", "T"))
    );
    for input in [
        "NM_000546.6:c.215C>G",        // transcript axis
        "NC_000017.11(TP53):g.100C>T", // gene selector on the accession
        "NC_000017.11:g.100del",       // not a substitution
        "NC_000017.11:g.100AC>T",      // multi-base reference
        "NC_000017.11:g.100C>TT",      // multi-base alternate
        "NC_000017.11:g.100N>T",       // not an ACGT base
        "NC_000017.11:g.100\u{c4}>T",  // a non-ASCII reference must not panic
        "NC_000017.11:g.C>T",          // no position
        "NC_000017.11:g.+100C>T",      // a signed position
        "NC_000017.11:g.0C>T",         // position below 1
        "NC_000017.11:g.100_101C>T",   // a range, not a position
    ] {
        assert_eq!(
            expected_spdi_from_genomic_substitution(input),
            None,
            "{input} must not reach the genomic-substitution oracle"
        );
    }
}

/// Count the conversions that carry a *successful* SPDI result — an
/// `spdi_result` object that is present and reports no error.
///
/// This mirrors exactly how `scripts/fetch_ncbi_variation.py` computes the
/// `hgvs_conversions.successful` header (`spdi_result` present and no `error`
/// key), so the self-consistency check tolerates either way the fetch harness
/// records a failed conversion: as a null `spdi_result` (the committed fixture's
/// convention) or as an `spdi_result` carrying an `error` field (the live
/// fetch's convention — the API's error body is stored inline). Counting merely
/// `spdi_result.is_some()` conflates the two and over-counts whenever a live
/// fetch hits any upstream API errors, which is the failure reported in issue
/// #2017 (`successful` = 19 vs 25 `spdi_result` present).
fn count_successful_conversions(conversions: &[HgvsConversion]) -> usize {
    conversions
        .iter()
        .filter(|conversion| {
            conversion
                .spdi_result
                .as_ref()
                .is_some_and(|result| result.error.is_none())
        })
        .count()
}

/// The `source` value `scripts/fetch_ncbi_variation.py` writes into a live
/// fetch. See [`FixtureKind`].
const LIVE_FETCH_SOURCE: &str = "NCBI Variation Services API";

/// Set to `1` by the weekly External API Validation workflow on the step that
/// runs these tests against a live fetch. See [`FixtureKind`].
const LIVE_FIXTURE_ENV: &str = "FERRO_SPDI_LIVE_FIXTURE";

/// Load the NCBI variation fixture: the committed curated one, or in the weekly
/// External API Validation run the live fetch written over it (see the module
/// docs). Either way the file is expected, so its absence is a real test
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
        "fixture: {:?}, generated {}",
        fixture.source, fixture.generated
    );
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
    // NCBI gives this deletion as NM_000492.4:1589:TCTT:T. SPDI counts along the
    // transcript sequence, not in c. numbering, and NCBI writes a deletion in a
    // repeat over the whole repeat with the kept base on both sides.
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

    // SPDI: NM_000546.6:356:C:G. SPDI counts from 0 along the transcript
    // sequence, and c.1 is n.143 on NM_000546.6, so c.215 is n.357, position 356.
    // HGVS: NM_000546.6:c.215C>G (1-based position 215)

    let result = parse_hgvs("NM_000546.6:c.215C>G");
    assert!(result.is_ok(), "Should parse with 1-based coordinate");
}

#[test]
fn test_spdi_half_open_intervals() {
    // SPDI uses half-open intervals [start, end)
    // HGVS uses closed intervals [start, end]

    // For a 3bp deletion:
    // On a genomic reference, g.1521_1523del is SPDI position 1520 deleting the
    // 3 bases in [1520, 1523). For the c. input below, SPDI counts along the
    // transcript sequence instead (NCBI: NM_000492.4:1589:TCTT:T).

    let result = parse_hgvs("NM_000492.4:c.1521_1523del");
    assert!(result.is_ok());
}

// ============================================================================
// Unresolvable end boundaries on the transcript axes (#1804)
// ============================================================================

/// The refusal added in #1804 must be visible to a **library consumer**, not
/// only to `src/spdi/convert.rs`'s own unit tests.
///
/// This is the integration half of those guards, and it asks a question they
/// cannot: the unit tests call the conversion through the module's own `use`,
/// so they would still pass if the error were swallowed by a wrapper on the way
/// out. `Normalizer::to_spdi` is that wrapper — it is what the PyO3 binding and
/// the `Normalizer` API route through — so it is the surface that has to
/// decline.
///
/// The collapse being pinned: on `main` at `439617c2`,
/// `NM_INTRON.1:c.10_(20_30)delAAAA` and `NM_INTRON.1:c.10_?delAAAA` both
/// converted to `NM_INTRON.1:19:AAAA:`, sharing that triple with each other and
/// with the fully-resolvable `c.10_13delAAAA`.
#[test]
fn an_unresolvable_transcript_end_is_refused_through_the_public_api() {
    use ferro_hgvs::normalize::Normalizer;
    use ferro_hgvs::reference::mock::MockProvider;
    use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
    use ferro_hgvs::spdi::convert::hgvs_to_spdi;

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_INTRON.1".to_string(),
        Some("INTRON".to_string()),
        Strand::Plus,
        "A".repeat(100),
        Some(11),
        Some(90),
        vec![Exon::new(1, 1, 50), Exon::new(2, 51, 100)],
        None,
        None,
        None,
        GenomeBuild::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    let normalizer = Normalizer::new(provider.clone());

    for descriptor in [
        "NM_INTRON.1:c.10_(20_30)delAAAA",
        "NM_INTRON.1:c.10_?delAAAA",
        "NM_INTRON.1:n.10_(20_30)delAAAA",
        "NM_INTRON.1:n.10_?delAAAA",
        "NM_INTRON.1:r.10_(20_30)delaaaa",
        "NM_INTRON.1:r.10_?delaaaa",
    ] {
        let variant = parse_hgvs(descriptor).expect("fixture must parse");

        let err = hgvs_to_spdi(&variant, &provider)
            .expect_err("an end naming no coordinate has no SPDI representation");
        assert!(
            err.to_string().contains("names no single coordinate"),
            "`{descriptor}` was refused for some other reason: {err}"
        );

        // The wrapper the Python binding and the `Normalizer` API use.
        let err = normalizer
            .to_spdi(&variant)
            .expect_err("the refusal must survive the Normalizer wrapper");
        assert!(
            err.to_string().contains("names no single coordinate"),
            "`{descriptor}` lost its reason through `Normalizer::to_spdi`: {err}"
        );
    }

    // Negative control on the same surface: `(13)` is `Mu::Uncertain` — a
    // parenthesised but perfectly numeric position — and must keep converting,
    // so this pins the refusal as "no coordinate" rather than "parenthesised".
    for descriptor in [
        "NM_INTRON.1:c.10_(13)delAAAA",
        "NM_INTRON.1:c.10_13delAAAA",
        "NM_INTRON.1:n.10_(13)delAAAA",
        "NM_INTRON.1:r.10_(13)delaaaa",
    ] {
        let variant = parse_hgvs(descriptor).expect("fixture must parse");
        normalizer
            .to_spdi(&variant)
            .unwrap_or_else(|e| panic!("`{descriptor}` must still convert, got {e}"));
    }
}
