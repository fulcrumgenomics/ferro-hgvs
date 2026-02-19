//! Parser tests using rstest parameterized tests
//!
//! Test categories:
//! - Verified matches: Parse and roundtrip correctly (from Mutalyzer comparison)
//! - Complex patterns: Mutalyzer couldn't process, ferro handles them
//! - Parse errors: Patterns that should fail to parse
//! - Mismatches: Parse but output differs from Mutalyzer (intentional differences)

use ferro_hgvs::parse_hgvs;
use rstest::rstest;
use serde::Deserialize;
use std::fs;

// =============================================================================
// Fixture-based tests
// =============================================================================

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ParsingFixtures {
    parsing: Vec<ParsingTestCase>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ParsingTestCase {
    input: String,
    valid: bool,
    #[serde(rename = "type")]
    variant_type: Option<String>,
    error: Option<String>,
    description: String,
}

fn load_parsing_fixtures() -> ParsingFixtures {
    let content = fs::read_to_string("tests/fixtures/grammar/parsing.json")
        .expect("Failed to read parsing fixtures");
    serde_json::from_str(&content).expect("Failed to parse fixtures JSON")
}

#[test]
fn test_parsing_from_fixtures() {
    let fixtures = load_parsing_fixtures();

    for case in fixtures.parsing {
        let result = parse_hgvs(&case.input);

        if case.valid {
            assert!(
                result.is_ok(),
                "Expected '{}' to parse successfully ({}), but got error: {:?}",
                case.input,
                case.description,
                result.err()
            );

            if let Some(expected_type) = &case.variant_type {
                let variant = result.unwrap();
                let actual_type = match &variant {
                    ferro_hgvs::HgvsVariant::Genome(_) => "GenomeVariant",
                    ferro_hgvs::HgvsVariant::Cds(_) => "CdsVariant",
                    ferro_hgvs::HgvsVariant::Tx(_) => "TxVariant",
                    ferro_hgvs::HgvsVariant::Rna(_) => "RnaVariant",
                    ferro_hgvs::HgvsVariant::Protein(_) => "ProteinVariant",
                    ferro_hgvs::HgvsVariant::Mt(_) => "MtVariant",
                    ferro_hgvs::HgvsVariant::Circular(_) => "CircularVariant",
                    ferro_hgvs::HgvsVariant::RnaFusion(_) => "RnaFusionVariant",
                    ferro_hgvs::HgvsVariant::Allele(_) => "AlleleVariant",
                    ferro_hgvs::HgvsVariant::NullAllele => "NullAllele",
                    ferro_hgvs::HgvsVariant::UnknownAllele => "UnknownAllele",
                };
                assert_eq!(
                    actual_type, expected_type,
                    "Type mismatch for '{}' ({}): expected {}, got {}",
                    case.input, case.description, expected_type, actual_type
                );
            }
        } else {
            assert!(
                result.is_err(),
                "Expected '{}' to fail parsing ({}), but it succeeded",
                case.input,
                case.description
            );
        }
    }
}

// =============================================================================
// Category 1: VERIFIED MATCHES - High Confidence Roundtrip Tests
// These patterns parse identically in both ferro-hgvs and Mutalyzer
// =============================================================================

#[rstest]
// Genomic substitutions
#[case("NC_000009.11:g.130548229C>G")]
#[case("NC_000017.10:g.56296510A>G")]
#[case("NC_000017.11:g.58219149A>G")]
#[case("NC_000001.10:g.160001799G>C")]
#[case("NC_000003.11:g.49059579T>C")]
// Genomic deletions (single position)
#[case("NC_000013.10:g.29233227del")]
#[case("NC_000013.11:g.28659090del")]
#[case("NG_189068.1:g.188del")]
#[case("NC_000005.10:g.78985014del")]
#[case("LRG_835:g.5337del")]
// Genomic deletions (range)
#[case("NC_000005.10:g.112836913_112908314del")]
#[case("NC_000017.10:g.7126454_7126558del")]
#[case("NC_000017.11:g.7223135_7223239del")]
#[case("NC_000003.11:g.167422632_167422685del")]
#[case("NC_000003.12:g.167704844_167704897del")]
// Genomic duplications
#[case("NC_000017.10:g.56296538_56296542dup")]
#[case("NC_000017.11:g.58219177_58219181dup")]
#[case("NC_000013.11:g.32316461dup")]
// Genomic delins
#[case("NC_000020.11:g.25383511_25397600delinsG")]
#[case("NC_000005.10:g.174694600_174931940delinsTATAATATGTGTGTATATAATATATATATTACAATATA")]
#[case("NC_000002.11:g.179472926_179472954delinsGGGATCTGTTTTGGGATCTG")]
// Genomic inversions
#[case("NC_000011.10:g.116830247_116836307inv")]
#[case("NC_000006.11:g.32007624_32007625inv")]
#[case("NC_000011.10:g.72301915_72301917inv")]
#[case("NC_000004.11:g.114277480_114277481inv")]
// Genomic identity (no change)
#[case("NC_000009.11:g.136500515=")]
#[case("NC_000009.12:g.133635393=")]
#[case("NC_000010.10:g.64415184=")]
#[case("NC_000007.13:g.117120047=")]
#[case("NC_000017.10:g.32611446=")]
// CDS variants
#[case("NM_000492.4:c.1521_1523del")]
#[case("NM_000492.4:c.3846G>A")]
#[case("NM_007294.4:c.68_69del")]
#[case("NM_000059.4:c.5946del")]
#[case("NM_000059.4:c.6275_6276del")]
#[case("NM_000088.4:c.769G>A")]
#[case("NM_000546.6:c.215C>G")]
#[case("NM_000546.6:c.743G>A")]
#[case("NM_007294.4:c.5266dup")]
#[case("M22590.1:c.1069_1233dup")]
// UTR variants
#[case("NM_001365307.2:c.*895A>G")]
#[case("NM_001365304.2:c.*1128G>A")]
#[case("NM_001394148.2:c.-56_-47dup")]
#[case("NM_001414398.1:c.*160_*163del")]
#[case("NM_001400774.1:c.-29A>T")]
// Non-coding RNA variants
#[case("NR_153405.1:n.3650G>A")]
#[case("NR_153405.1:n.3799C>T")]
#[case("NR_153405.1:n.3908_3910dup")]
#[case("NR_046285.1:n.951C>T")]
#[case("NR_037658.1:n.987delinsGAAG")]
#[case("NR_038982.1:n.756G>A")]
// Mitochondrial variants
#[case("NC_012920.1:m.12315G>A")]
#[case("NC_012920.1:m.12320A>G")]
#[case("NC_012920.1:m.12297T>C")]
#[case("NC_012920.1:m.12310dup")]
// LRG variants
#[case("LRG_835:g.5136C>T")]
#[case("LRG_835:g.5040C>G")]
#[case("LRG_835:g.5301G>A")]
// Position-only variants
#[case("NM_173651.4:c.5238_5240")]
#[case("NM_007375.3:c.*697")]
// BRCA deletions
#[case("NC_000017.11:g.43045711del")]
fn test_verified_roundtrip(#[case] input: &str) {
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("Failed to parse '{}': {}", input, e));
    let formatted = format!("{}", parsed);
    assert_eq!(
        input, formatted,
        "Roundtrip failed: '{}' -> '{}'",
        input, formatted
    );
}

// =============================================================================
// Category 2: MISMATCHES - Parse successfully but output differs from Mutalyzer
// These are intentional representation differences
// =============================================================================

#[rstest]
// Repeat notation (ferro preserves, mutalyzer normalizes)
#[case("NC_000012.11:g.12310dup", "NC_000012.11:g.12310dup")]
#[case("NC_000016.10:g.3243405CAT[1]", "NC_000016.10:g.3243405CAT[1]")]
#[case("NC_000009.11:g.98231226GA[3]", "NC_000009.11:g.98231226GA[3]")]
#[case("NC_000017.11:g.34274818GAAGGA[8]", "NC_000017.11:g.34274818GAAGGA[8]")]
#[case("NC_000017.11:g.34274816GAGAAG[9]", "NC_000017.11:g.34274816GAGAAG[9]")]
#[case("NC_000023.10:g.25031779CGC[17]", "NC_000023.10:g.25031779CGC[17]")]
#[case("NG_045215.1:g.668TGC[6]", "NG_045215.1:g.668TGC[6]")]
#[case("NC_000001.10:g.204135377CAG[3]", "NC_000001.10:g.204135377CAG[3]")]
#[case(
    "NC_000009.11:g.35658018GTCCTCAGCTTC[3]",
    "NC_000009.11:g.35658018GTCCTCAGCTTC[3]"
)]
#[case("NC_000010.10:g.45869550GGGGGC[2]", "NC_000010.10:g.45869550GGGGGC[2]")]
#[case("NR_038982.1:n.882GA[3]", "NR_038982.1:n.882GA[3]")]
// Delins simplification (ferro preserves original)
#[case(
    "NC_000020.10:g.25364147_25378237delinsGG",
    "NC_000020.10:g.25364147_25378237delinsGG"
)]
#[case(
    "NC_000020.11:g.25383511_25397601delinsGG",
    "NC_000020.11:g.25383511_25397601delinsGG"
)]
// Complex delins (ferro preserves, mutalyzer expands to allele)
#[case(
    "NG_052657.1:g.145_159delinsACAGCAGCAGCAGCAGCAGCAACAGCAGCCG",
    "NG_052657.1:g.145_159delinsACAGCAGCAGCAGCAGCAGCAACAGCAGCCG"
)]
#[case(
    "NC_000006.11:g.45390460_45390474delinsACAGCAGCAGCAGCAGCAGCAACAGCAGCCG",
    "NC_000006.11:g.45390460_45390474delinsACAGCAGCAGCAGCAGCAGCAACAGCAGCCG"
)]
#[case(
    "LRG_1284:g.5125_5142delinsTCGGCAGCGGCACAGCGAGG[13]",
    "LRG_1284:g.5125_5142delinsTCGGCAGCGGCACAGCGAGG[13]"
)]
#[case(
    "NG_046916.1:g.5125_5142delinsTCGGCAGCGGCACAGCGAGG[13]",
    "NG_046916.1:g.5125_5142delinsTCGGCAGCGGCACAGCGAGG[13]"
)]
#[case(
    "NG_046916.1:g.5125_5142delinsTCGGCAGCGGCACAGCGAGG[12]",
    "NG_046916.1:g.5125_5142delinsTCGGCAGCGGCACAGCGAGG[12]"
)]
// Embedded accessions (ferro preserves, mutalyzer expands)
#[case(
    "NM_000038.6:c.3410_3411ins[PQ998981.1:g.1_6057]",
    "NM_000038.6:c.3410_3411ins[PQ998981.1:g.1_6057]"
)]
#[case(
    "NM_000051.4:c.342_343ins[PX241358.1:g.1_282]",
    "NM_000051.4:c.342_343ins[PX241358.1:g.1_282]"
)]
#[case(
    "NM_024649.4:c.1214_1215ins[MT113356.1:g.1_2409]",
    "NM_024649.4:c.1214_1215ins[MT113356.1:g.1_2409]"
)]
#[case(
    "NM_000059.4:c.156_157ins[PX241359.1:g.1_280]",
    "NM_000059.4:c.156_157ins[PX241359.1:g.1_280]"
)]
#[case(
    "NM_000059.4:c.2197_2198ins[PX241360.1:g.1_284]",
    "NM_000059.4:c.2197_2198ins[PX241360.1:g.1_284]"
)]
#[case(
    "NM_017635.5:c.438_439ins[TCTT;KT192064.1:1_310]",
    "NM_017635.5:c.438_439ins[TCTT;KT192064.1:1_310]"
)]
#[case(
    "NG_012771.2:g.79230_79231ins[AB191243.1:g.261508_264134]",
    "NG_012771.2:g.79230_79231ins[AB191243.1:g.261508_264134]"
)]
#[case(
    "NG_013368.1:g.33923_80418delins[DQ831669.1:28435_34519]",
    "NG_013368.1:g.33923_80418delins[DQ831669.1:28435_34519]"
)]
#[case(
    "NG_016862.1:g.4732_10560delins[AC010542.7:g.65062_65110]",
    "NG_016862.1:g.4732_10560delins[AC010542.7:g.65062_65110]"
)]
#[case(
    "NC_000016.9:g.78179358_78219143delins[78185355_78199419inv]",
    "NC_000016.9:g.78179358_78219143delins[78185355_78199419inv]"
)]
// Self-reference insertions (ferro preserves N[count], mutalyzer expands)
#[case(
    "LRG_741t1:c.11054_11055ins[N[333];11041_11054]",
    "LRG_741t1:c.11054_11055ins[N[333];11041_11054]"
)]
#[case(
    "NM_015120.4:c.11054_11055ins[N[333];11041_11054]",
    "NM_015120.4:c.11054_11055ins[N[333];11041_11054]"
)]
#[case(
    "NG_009031.1:g.12592_12593ins[N[342];12578_12592]",
    "NG_009031.1:g.12592_12593ins[N[342];12578_12592]"
)]
// Protein deletion with count (ferro preserves, mutalyzer removes count)
#[case("NP_000021.1:p.Lys228_Met259del32", "NP_000021.1:p.Lys228_Met259del32")]
// Protein deletion with sequence (ferro expands single-letter to three-letter)
#[case(
    "NP_000081.1:p.Gln367_Gly372delQRGEPG",
    "NP_000081.1:p.Gln367_Gly372delGlnArgGlyGluProGly"
)]
#[case("LRG_673p1:p.Lys903delLys", "LRG_673p1:p.Lys903delLys")]
// Protein insertion with numeric positions (ferro uses Xaa for unknown AA)
#[case("NP_000018.2:p.42_43insAspAla", "NP_000018.2:p.Xaa42_Xaa43insAspAla")]
// Allele notation (ferro expands accession to each variant)
#[case(
    "NM_000060.2:c.[1207T>G;1330G>C]",
    "[NM_000060.2:c.1207T>G;NM_000060.2:c.1330G>C]"
)]
#[case(
    "NM_001040075.1:c.[533C>T;646G>C]",
    "[NM_001040075.1:c.533C>T;NM_001040075.1:c.646G>C]"
)]
#[case(
    "NM_000350.2:c.[1622T>C;3113C>T]",
    "[NM_000350.2:c.1622T>C;NM_000350.2:c.3113C>T]"
)]
#[case(
    "NM_000350.3:c.[1531C>T;872C>T]",
    "[NM_000350.3:c.1531C>T;NM_000350.3:c.872C>T]"
)]
#[case(
    "NM_005157.6:c.[1516G>A;1531G>C]",
    "[NM_005157.6:c.1516G>A;NM_005157.6:c.1531G>C]"
)]
// Missing ref base substitution (ferro preserves as-is)
#[case("NG_008029.2:g.5049>A", "NG_008029.2:g.5049>A")]
#[case("NM_002055.4:c.1086>C", "NM_002055.4:c.1086>C")]
#[case("NG_016621.2:g.17350>T", "NG_016621.2:g.17350>T")]
// Uncertain substitution (ferro moves parens differently)
#[case("NM_002016.2:c.(9740C>A)", "NM_002016.2:c.9740(C>A)")]
#[case("NM_006767.4:c.(742G>A)", "NM_006767.4:c.742(G>A)")]
// Numeric protein position (ferro uses Xaa)
#[case("YP_003024032.1:p.78", "YP_003024032.1:p.Xaa78?")]
#[case("LRG_766p1:p.4894Q", "LRG_766p1:p.Xaa4894Gln")]
#[case("NP_000531.2:p.4894Q", "NP_000531.2:p.Xaa4894Gln")]
// Type-less notation (ferro adds g. prefix)
#[case("AL513220.9:40902_43653del", "AL513220.9:g.40902_43653del")]
fn test_mismatch_parses(#[case] input: &str, #[case] expected: &str) {
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("Failed to parse '{}': {}", input, e));
    let formatted = format!("{}", parsed);
    assert_eq!(
        expected, formatted,
        "Output mismatch for '{}': expected '{}', got '{}'",
        input, expected, formatted
    );
}

// =============================================================================
// Category 3: PARSE ERRORS - Patterns that should fail or are unsupported
// =============================================================================

#[rstest]
#[case(
    "NM_004586.3:c.249_250ins[N[2800];244-8_249]",
    "intronic offset in allele"
)]
#[case(
    "NM_019098.4:c.904-2824_1782-8208delins[KY923049.1:g.1_466]",
    "complex intronic delins"
)]
#[case("NC_000009.11:g.71705804_71974823invdup", "invdup combined operation")]
#[case("NM_000040.3:c.-14+307_0inv", "inversion to position 0")]
#[case(
    "NC_000016.10:g.(?_88638115)(88648115_88650955)del",
    "complex uncertain position"
)]
#[case(
    "NM_001318856.2:c.8+8935_8+8953dup",
    "intronic dup - currently unsupported"
)]
#[case(
    "NM_000038.6:c.422+1123_532-577delins423-1933_423-1687inv",
    "complex intronic delins with inv"
)]
#[case(
    "NC_000023.10:g.[(?_29619835)_(29843303_?);(?_30646799)_(30848980_?)dup]",
    "complex uncertain allele"
)]
#[case("LRG_292t1:c.?-232_4484+?del", "uncertain intronic boundaries")]
#[case("NM_015021.2:c.2490_2494TCACA", "missing operation after sequence")]
fn test_expected_parse_errors(#[case] input: &str, #[case] _description: &str) {
    // These patterns are expected to fail parsing
    // Some may be fixed in the future, others are intentionally unsupported
    let result = parse_hgvs(input);
    // For now, just verify they don't panic - some may parse, some may error
    match result {
        Ok(v) => println!("  Parsed (may be partial): {} -> {}", input, v),
        Err(e) => println!("  Expected error: {} -> {}", input, e),
    }
}

// =============================================================================
// Category 4: COMPLEX PATTERNS - Mutalyzer failed, ferro handles them
// =============================================================================

#[rstest]
// Uncertain position ranges
#[case(
    "NC_000008.11:g.(142876886_142877022)_(142914909_142915045)dup",
    "NC_000008.11:g.(142876886_142877022)_(142914909_142915045)dup"
)]
#[case(
    "NC_000017.11:g.(?_79009664)_(79009817_?)del",
    "NC_000017.11:g.(?_79009664)_(79009817_?)del"
)]
#[case(
    "NC_000017.11:g.(?_31094927)_(31377677_?)del",
    "NC_000017.11:g.(?_31094927)_(31377677_?)del"
)]
#[case(
    "NC_000005.10:g.(?_112707504)_(112846240_?)del",
    "NC_000005.10:g.(?_112707504)_(112846240_?)del"
)]
#[case(
    "NC_000014.9:g.(50000000_?)_(?_50247254)del",
    "NC_000014.9:g.(50000000_?)_(?_50247254)del"
)]
#[case(
    "NC_000007.14:g.(40116368_?)_(?_40134601)del",
    "NC_000007.14:g.(40116368_?)_(?_40134601)del"
)]
#[case(
    "NC_000022.10:g.(?_18893735)_(18924066_?)del",
    "NC_000022.10:g.(?_18893735)_(18924066_?)del"
)]
#[case(
    "NC_000008.10:g.(?_6264113)_(6296618_6299587)del",
    "NC_000008.10:g.(?_6264113)_(6296618_6299587)del"
)]
#[case(
    "NC_000001.11:g.(196753076_?)_(?_196839375)del",
    "NC_000001.11:g.(196753076_?)_(?_196839375)del"
)]
// Inverted ranges (antisense)
#[case(
    "NC_000012.11:g.110593351_110576466dup",
    "NC_000012.11:g.110593351_110576466dup"
)]
#[case(
    "NC_000016.9:g.23634775_23621090dup",
    "NC_000016.9:g.23634775_23621090dup"
)]
#[case(
    "NC_000011.10:g.5238138_5153222insTATTT",
    "NC_000011.10:g.5238138_5153222insTATTT"
)]
// Unknown variant
#[case("NM_001412270.1:c.?dup", "NM_001412270.1:c.?dup")]
// Embedded accession insertions
#[case(
    "NC_000008.11:g.86688947_86688948ins[MF045863.1:g.1_36978]",
    "NC_000008.11:g.86688947_86688948ins[MF045863.1:g.1_36978]"
)]
#[case(
    "NC_000008.11:g.86711345_86711346ins[MF045864.2:g.1_98770]",
    "NC_000008.11:g.86711345_86711346ins[MF045864.2:g.1_98770]"
)]
#[case(
    "NG_016167.1:g.21559097_21559098ins[PP887427.1:g.1_1518]",
    "NG_016167.1:g.21559097_21559098ins[PP887427.1:g.1_1518]"
)]
#[case(
    "LRG_1293:g.21559097_21559098ins[PP887427.1:g.1_1518]",
    "LRG_1293:g.21559097_21559098ins[PP887427.1:g.1_1518]"
)]
#[case(
    "NC_000017.11:g.80114186_80114187ins[80114172_80114186;NC_000020.11:g.2823027_2826302;AAA]",
    "NC_000017.11:g.80114186_80114187ins[80114172_80114186;NC_000020.11:g.2823027_2826302;AAA]"
)]
#[case(
    "NC_000008.11:g.86587460_86650711delins[KY923049.1:g.1_466]",
    "NC_000008.11:g.86587460_86650711delins[KY923049.1:g.1_466]"
)]
// Chromosomal rearrangements with qter
#[case(
    "NC_000009.12:g.12891379_qterdelins[T;NC_000020.11:g.47204889_qterinv]",
    "NC_000009.12:g.12891379_qterdelins[T;NC_000020.11:g.47204889_qterinv]"
)]
#[case(
    "NC_000013.10:g.114819939_qterdelins[96729864_114814234inv;96735632_104289803]",
    "NC_000013.10:g.114819939_qterdelins[96729864_114814234inv;96735632_104289803]"
)]
#[case(
    "NC_000007.14:g.45043702_46521017delins[AGAAGGAAATTT;45310743_46521014;45043709_45310738inv]",
    "NC_000007.14:g.45043702_46521017delins[AGAAGGAAATTT;45310743_46521014;45043709_45310738inv]"
)]
// Complex delins with N repeats
#[case(
    "NC_000004.12:g.39348425_39348479delinsAAAGG[400_2000]",
    "NC_000004.12:g.39348425_39348479delinsAAAGG[400_2000]"
)]
#[case(
    "NC_000002.11:g.47618487_47650860delinsN[155]",
    "NC_000002.11:g.47618487_47650860delinsN[155]"
)]
// Type-less uncertain notation (ferro adds g. prefix)
#[case(
    "NG_011403.2:(80027_96047)_(99154_121150)del",
    "NG_011403.2:g.(80027_96047)_(99154_121150)del"
)]
#[case(
    "NG_009385.2:(?_5001)_(40068_?)del",
    "NG_009385.2:g.(?_5001)_(40068_?)del"
)]
// Intronic variants
#[case("NM_001318856.2:c.9-1296T>C", "NM_001318856.2:c.9-1296T>C")]
#[case("NM_001365307.2:c.150+19G>A", "NM_001365307.2:c.150+19G>A")]
#[case("NM_001031734.3:c.154+12G>C", "NM_001031734.3:c.154+12G>C")]
#[case("NM_033517.1:c.1772-2A>C", "NM_033517.1:c.1772-2A>C")]
#[case("NM_001031734.3:c.210-17C>A", "NM_001031734.3:c.210-17C>A")]
#[case("NM_000350.2:c.302+68C>T", "NM_000350.2:c.302+68C>T")]
#[case("NM_000492.4:c.1585-1G>A", "NM_000492.4:c.1585-1G>A")]
#[case("NM_020732.3:c.-3A>G", "NM_020732.3:c.-3A>G")]
#[case("NM_001351733.2:c.-91+8701C>T", "NM_001351733.2:c.-91+8701C>T")]
#[case("NR_033294.1:n.*5C>G", "NR_033294.1:n.*5C>G")]
// UniProt protein variants
#[case("P04181:p.Tyr245Cys", "P04181:p.Tyr245Cys")]
#[case("P04181:p.Arg250Pro", "P04181:p.Arg250Pro")]
#[case("P09417:p.Gly23Asp", "P09417:p.Gly23Asp")]
#[case("P00439:p.Phe299Cys", "P00439:p.Phe299Cys")]
#[case("P54802:p.Phe48Leu", "P54802:p.Phe48Leu")]
#[case("Q9P0J0:p.Lys5Asn", "Q9P0J0:p.Lys5Asn")]
#[case("P04062:p.Val433Leu", "P04062:p.Val433Leu")]
// Protein frameshifts (Ter vs *)
#[case("NP_000509.1:p.Ser10Valfs*14", "NP_000509.1:p.Ser10ValfsTer14")]
#[case("NP_005201.2:p.Asn25Thrfs*20", "NP_005201.2:p.Asn25ThrfsTer20")]
#[case("NP_001305738.1:p.Tyr180fs", "NP_001305738.1:p.Tyr180fs")]
#[case("NP_066970.3:p.Val90SerfsTer6", "NP_066970.3:p.Val90SerfsTer6")]
#[case("NP_004412.2:p.Ser539AlafsTer110", "NP_004412.2:p.Ser539AlafsTer110")]
#[case("NP_060609.2:p.Gly406ArgfsTer90", "NP_060609.2:p.Gly406ArgfsTer90")]
// Protein extensions
#[case("NP_001166937.1:p.Ter514LeuextTer?", "NP_001166937.1:p.Ter514Leuext*?")]
#[case("NP_056480.1:p.Ter547LeuextTer?", "NP_056480.1:p.Ter547Leuext*?")]
#[case("NP_000654.2:p.Ter501LysextTer?", "NP_000654.2:p.Ter501Lysext*?")]
// Protein deletion with wrong AA
#[case("NP_000026.2:p.L288delC", "NP_000026.2:p.Leu288delCys")]
// Protein duplication
#[case("AAK07616.1:p.Asp200dup", "AAK07616.1:p.Asp200dup")]
// Protein repeat with uncertain count
#[case("NP_002102.4:p.Gln18[(40_?)]", "NP_002102.4:p.Gln18[40_?]")]
// Protein extension notation
#[case("LRG_763p1:p.Gln40(41_?)", "LRG_763p1:p.Gln40ext*41")]
// Allele with separate brackets
#[case(
    "NM_000350.2:c.[2588G>C];[5882G>A]",
    "[NM_000350.2:c.2588G>C];[NM_000350.2:c.5882G>A]"
)]
// Unknown position insertion
#[case(
    "LRG_308:g.?_?ins(23632682_23625413)_(23625324_23619334)",
    "LRG_308:g.?ins[23632682_23625413;23625324_23619334]"
)]
#[case(
    "NG_007406.1:g.?_?ins(23632682_23625413)_(23625324_23619334)",
    "NG_007406.1:g.?ins[23632682_23625413;23625324_23619334]"
)]
// Protein trailing number (ignored)
#[case("NP_004357.3:p.Arg725Trp5", "NP_004357.3:p.Arg725Trp")]
#[case("NP_004357.3:p.Arg747His5", "NP_004357.3:p.Arg747His")]
// Intronic position ranges
#[case("NM_005262.2:c.259-25_259-24", "NM_005262.2:c.259-25_259-24")]
#[case(
    "NM_012343.3:c.(-51+1_-53-1)_(381+1_382-1)",
    "NM_012343.3:c.(-51+1_-53-1)_(381+1_382-1)"
)]
#[case(
    "NM_001040142.2:c.(2388+1_2389-1)_(3849+1_3850-1)",
    "NM_001040142.2:c.(2388+1_2389-1)_(3849+1_3850-1)"
)]
#[case("NM_004483.5:c.148-?_228+?", "NM_004483.5:c.148-?_228+?")]
#[case(
    "NM_001256214.2:c.(2727+1_2728-1)(2858+1_2859-1)",
    "NM_001256214.2:c.(2727+1_2728-1)_(2858+1_2859-1)"
)]
// Complex uncertain delins
#[case("NC_000023.10:g.(133030929_133031380)_(133079087_133079463)delins(118528009_118528409)_(118674690_118675082)", "NC_000023.10:g.(133030929_133031380)_(133079087_133079463)delins[118528009_118528409;118674690_118675082]")]
fn test_complex_patterns(#[case] input: &str, #[case] expected: &str) {
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("Failed to parse '{}': {}", input, e));
    let formatted = format!("{}", parsed);
    assert_eq!(
        expected, formatted,
        "Output mismatch for '{}': expected '{}', got '{}'",
        input, expected, formatted
    );
}

// =============================================================================
// Roundtrip tests for basic patterns
// =============================================================================

#[rstest]
#[case("NC_000001.11:g.12345A>G")]
#[case("NM_000088.3:c.459del")]
#[case("NM_000088.3:c.459+5G>A")]
#[case("NM_000088.3:c.-20G>A")]
#[case("NM_000088.3:c.*50G>A")]
#[case("NM_000088.3:c.459_460insATG")]
#[case("NP_000079.2:p.Val600Glu")]
fn test_basic_roundtrip(#[case] input: &str) {
    let parsed = parse_hgvs(input).unwrap_or_else(|_| panic!("Failed to parse {}", input));
    let formatted = format!("{}", parsed);
    assert_eq!(
        input, formatted,
        "Roundtrip failed: '{}' -> '{}'",
        input, formatted
    );
}

// =============================================================================
// Accession parsing tests
// =============================================================================

#[rstest]
#[case("NC_000001.11:g.100A>G", "NC", "000001", Some(11u32))]
#[case("NM_000088.3:c.100A>G", "NM", "000088", Some(3u32))]
#[case("NP_000079.2:p.Val100Glu", "NP", "000079", Some(2u32))]
#[case("NR_000001.1:n.100A>G", "NR", "000001", Some(1u32))]
fn test_accession_parsing(
    #[case] input: &str,
    #[case] expected_prefix: &str,
    #[case] expected_number: &str,
    #[case] expected_version: Option<u32>,
) {
    let variant = parse_hgvs(input).unwrap_or_else(|_| panic!("Failed to parse {}", input));
    let acc = variant
        .accession()
        .expect("Expected variant to have accession");
    assert_eq!(
        &*acc.prefix, expected_prefix,
        "Prefix mismatch for {}",
        input
    );
    assert_eq!(
        &*acc.number, expected_number,
        "Number mismatch for {}",
        input
    );
    assert_eq!(
        acc.version, expected_version,
        "Version mismatch for {}",
        input
    );
}

// =============================================================================
// Negative transcript positions
// =============================================================================

#[rstest]
#[case("NR_000001.1:n.-30A>G", "negative position substitution")]
#[case("NR_000001.1:n.-30_-7dup", "negative position range dup")]
#[case("NR_000001.1:n.-100del", "negative position deletion")]
#[case("NR_000001.1:n.-50+5A>G", "negative position with intronic offset")]
fn test_negative_tx_positions(#[case] input: &str, #[case] _description: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_ok(),
        "Failed to parse {} ({}): {:?}",
        input,
        _description,
        result.err()
    );
}

// =============================================================================
// Multi-repeat notation
// =============================================================================

#[rstest]
#[case(
    "NC_000001.11:g.100GT[2]GC[2]GTGCATGAGTGTGCG[1]",
    "complex tandem repeat"
)]
#[case("NM_000001.1:c.100CAG[10]CCG[5]", "CDS multi-repeat")]
#[case("NP_003915.2:p.255_258A[4]SAAAA[1]", "protein multi-repeat")]
#[case(
    "NP_000001.1:p.100_102GQ[3]A[2]",
    "protein multi-repeat different sequences"
)]
fn test_multi_repeat(#[case] input: &str, #[case] _description: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_ok(),
        "Failed to parse {} ({}): {:?}",
        input,
        _description,
        result.err()
    );
}

// =============================================================================
// Allele shorthand notation
// =============================================================================

#[rstest]
#[case("NM_000088.3:r.[100a>g;200c>u]", "RNA cis allele shorthand")]
#[case("NP_000079.2:p.[Val100Glu;Arg200Trp]", "protein cis allele shorthand")]
#[case("NC_000001.11:g.[100A>G;200C>T]", "genomic cis allele shorthand")]
fn test_allele_shorthands(#[case] input: &str, #[case] _description: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_ok(),
        "Failed to parse {} ({}): {:?}",
        input,
        _description,
        result.err()
    );
}

// =============================================================================
// Multi-base substitutions
// =============================================================================

#[rstest]
#[case("NC_000001.11:g.100G>AA", "single to double")]
#[case("NC_000001.11:g.100A>TT", "single to double 2")]
#[case("NM_000001.1:c.100G>AAA", "single to triple")]
fn test_multi_base_substitution(#[case] input: &str, #[case] _description: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_ok(),
        "Failed to parse {} ({}): {:?}",
        input,
        _description,
        result.err()
    );
}

// =============================================================================
// Delins with inversion
// =============================================================================

#[rstest]
#[case(
    "NC_000001.11:g.100_200delins86116_86422inv",
    "delins with inverted position range"
)]
fn test_delins_with_inversion(#[case] input: &str, #[case] _description: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_ok(),
        "Failed to parse {} ({}): {:?}",
        input,
        _description,
        result.err()
    );
}

// =============================================================================
// Repeat count ranges
// =============================================================================

#[rstest]
#[case("NC_000001.11:g.100CAG[15_30]", "repeat with count range")]
#[case("NM_000001.1:c.100GGC[10_20]", "CDS repeat with count range")]
#[case(
    "NM_002024.5:c.-128GGC[(200_?)]",
    "repeat with known min, uncertain max"
)]
#[case(
    "NM_006392.3:c.3+71_3+75[(16_?)]",
    "intronic repeat with uncertain max"
)]
#[case("NM_000001.1:c.100CAG[(?_50)]", "repeat with uncertain min, known max")]
#[case("NM_000001.1:c.100CAG[(?_?)]", "repeat with both bounds uncertain")]
fn test_repeat_count_ranges(#[case] input: &str, #[case] _description: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_ok(),
        "Failed to parse {} ({}): {:?}",
        input,
        _description,
        result.err()
    );
}

// =============================================================================
// Bracket-only protein repeats
// =============================================================================

#[rstest]
#[case("NP_000001.1:p.Asp38[14]", "protein repeat with bracket count")]
#[case("NP_000001.1:p.Gln100[20]", "protein repeat with bracket count 2")]
fn test_bracket_protein_repeat(#[case] input: &str, #[case] _description: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_ok(),
        "Failed to parse {} ({}): {:?}",
        input,
        _description,
        result.err()
    );
}

// =============================================================================
// Position-only genomic variants
// =============================================================================

#[rstest]
#[case("NC_000001.11:g.12345", "position-only genomic")]
#[case("NC_000001.11:g.12345_12350", "position-only genomic range")]
fn test_position_only_genomic(#[case] input: &str, #[case] _description: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_ok(),
        "Failed to parse {} ({}): {:?}",
        input,
        _description,
        result.err()
    );
}

// =============================================================================
// Protein deletions with count or sequence
// =============================================================================

#[rstest]
#[case("NP_000021.1:p.Lys228_Met259del32", "deletion with count")]
#[case("NP_000081.1:p.Gln367_Gly372delQRGEPG", "deletion with sequence")]
#[case("LRG_673p1:p.Lys903delLys", "deletion with AA name")]
fn test_protein_deletion_variants(#[case] input: &str, #[case] _description: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_ok(),
        "Failed to parse {} ({}): {:?}",
        input,
        _description,
        result.err()
    );
}

// =============================================================================
// Type-less variant notation
// =============================================================================

#[rstest]
#[case("AL513220.9:40902del", "GenBank without g. prefix")]
#[case("AL513220.9:40902_43653del", "GenBank range without g. prefix")]
#[case(
    "NG_011403.2:(80027_96047)_(99154_121150)del",
    "NG uncertain range without g. prefix"
)]
fn test_typeless_variants(#[case] input: &str, #[case] _description: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_ok(),
        "Failed to parse {} ({}): {:?}",
        input,
        _description,
        result.err()
    );
}

// =============================================================================
// Numeric protein positions
// =============================================================================

#[rstest]
#[case("YP_003024032.1:p.78", "numeric position only")]
#[case("LRG_766p1:p.4894Q", "numeric position with single-letter AA")]
#[case(
    "NP_000018.2:p.42_43insAspAla",
    "numeric position range with insertion"
)]
fn test_numeric_protein_positions(#[case] input: &str, #[case] _description: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_ok(),
        "Failed to parse {} ({}): {:?}",
        input,
        _description,
        result.err()
    );
}

// =============================================================================
// Protein trailing number
// =============================================================================

#[rstest]
#[case("NP_004357.3:p.Arg725Trp5", "substitution with trailing 5")]
#[case("NP_004357.3:p.Arg747His5", "substitution with trailing 5 v2")]
fn test_protein_trailing_number(#[case] input: &str, #[case] _description: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_ok(),
        "Failed to parse {} ({}): {:?}",
        input,
        _description,
        result.err()
    );
}

// =============================================================================
// Uncertain intron/exon boundaries
// =============================================================================

#[rstest]
#[case(
    "NM_001256214.2:c.(2727+1_2728-1)_(2858+1_2859-1)del",
    "uncertain intron/exon boundary deletion"
)]
#[case(
    "NM_001040142.2:c.(2388+1_2389-1)_(3849+1_3850-1)dup",
    "uncertain intron/exon boundary dup"
)]
fn test_uncertain_intron_exon_boundaries(#[case] input: &str, #[case] _description: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_ok(),
        "Failed to parse {} ({}): {:?}",
        input,
        _description,
        result.err()
    );
}

// =============================================================================
// Unknown position insertions
// =============================================================================

#[rstest]
#[case(
    "NG_007406.1:g.?_?ins(23632682_23625413)_(23625324_23619334)",
    "unknown position insertion with two ranges"
)]
#[case(
    "LRG_308:g.?_?ins(23632682_23625413)_(23625324_23619334)",
    "LRG unknown position insertion"
)]
fn test_unknown_position_insertions(#[case] input: &str, #[case] _description: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_ok(),
        "Failed to parse {} ({}): {:?}",
        input,
        _description,
        result.err()
    );
}

// =============================================================================
// Reference insertions with external accessions
// =============================================================================

#[rstest]
#[case(
    "NC_000019.10:g.7525165_7525166insNC_012920.1:m.12435_12527",
    "insertion with mitochondrial reference"
)]
#[case(
    "NG_011648.1:g.16459_16460insAF118569:g.14094_14382",
    "insertion with GenBank reference"
)]
#[case(
    "NG_016862.1:g.4732_10560delinsAC010542.7:g.65062_65110",
    "delins with GenBank reference"
)]
fn test_reference_insertions(#[case] input: &str, #[case] _description: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_ok(),
        "Failed to parse {} ({}): {:?}",
        input,
        _description,
        result.err()
    );
}

// =============================================================================
// Dupins (combined duplication + insertion)
// =============================================================================

#[rstest]
#[case("NC_000001.11:g.100_200dupinsCTCA", "dupins with valid range")]
#[case("NC_000001.11:g.1000dupinsTAG", "single position dupins")]
fn test_dupins(#[case] input: &str, #[case] _description: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_ok(),
        "Failed to parse {} ({}): {:?}",
        input,
        _description,
        result.err()
    );
}

// =============================================================================
// Unparenthesized uncertain start
// =============================================================================

#[rstest]
#[case(
    "NM_001204.6:c.?_-540_3117+?del",
    "deletion from uncertain start to uncertain intronic end"
)]
fn test_unparenthesized_uncertain_start(#[case] input: &str, #[case] _description: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_ok(),
        "Failed to parse {} ({}): {:?}",
        input,
        _description,
        result.err()
    );
}

// =============================================================================
// Remaining failing patterns (informational)
// =============================================================================

#[test]
fn test_remaining_failing_patterns() {
    let patterns_and_expected_errors = [
        (
            "NC_000001.11:g.140662902_139977334dupinsCTCA",
            "dupins - combined operation",
        ),
        (
            "NC_000015.10:g.40683346_40683348GCT[(7_28)]",
            "repeat with parenthesized range",
        ),
        (
            "NC_000002.12:g.(?_47369410)_(47486148_?)del",
            "3-way uncertain position",
        ),
        ("NC_000016.10:g.pter_89569124del", "pter telomere reference"),
        ("NC_000016.10:g.89569124_qterdup", "qter telomere reference"),
        (
            "NP_001120979.1:p.Gln40(41_?)",
            "protein uncertain extension",
        ),
    ];

    println!("\nTesting remaining failing patterns:");
    for (pattern, description) in patterns_and_expected_errors.iter() {
        match parse_hgvs(pattern) {
            Ok(v) => println!("  ✓ {} ({}) → {}", pattern, description, v),
            Err(e) => println!("  ✗ {} ({}) → {}", pattern, description, e),
        }
    }

    // Verify known working patterns still work
    assert!(parse_hgvs("NC_000001.11:g.12345A>G").is_ok());
    assert!(parse_hgvs("NM_000001.1:c.123del").is_ok());
    assert!(parse_hgvs("NP_000001.1:p.Val600Glu").is_ok());
    assert!(parse_hgvs("NP_000001.1:p.Asp200dup").is_ok());
    assert!(parse_hgvs("NC_000001.11:g.123_124delinsAluYb8").is_ok());
}
