//! Cross-axis conformance for the HGVS bracket/allele cardinality rule (#493).
//!
//! `[ ]` is allele syntax. The spec admits only two conformant shapes, identically
//! across every coordinate system (`c/g/n/m/o/r/p`): one bracket group with ≥2 cis
//! members (`c.[76A>C;88G>T]`), or ≥2 trans groups (`c.[76A>C];[88G>T]`). A standalone
//! single-member bracket (`c.[76A>C]`, `p.[=]`, `c.[?]`) is non-conformant; the
//! canonical repair drops the redundant brackets. Strict rejects (W3026), lenient
//! unwraps + warns, silent unwraps silently.

use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::HgvsVariant;

const W: &str = "W3026";

/// Every standalone single-member bracket form, paired with its canonical unwrap.
/// Covers all seven coordinate systems plus the whole-entity / marker members.
fn nonconformant_singletons() -> Vec<(&'static str, &'static str)> {
    vec![
        // one per coordinate axis
        ("NC_000001.11:g.[1000G>A]", "NC_000001.11:g.1000G>A"),
        ("NM_000088.3:c.[76A>C]", "NM_000088.3:c.76A>C"),
        ("NR_000088.3:n.[76a>c]", "NR_000088.3:n.76A>C"),
        ("NM_000088.3:r.[76a>c]", "NM_000088.3:r.76a>c"),
        ("NP_000079.2:p.[Val600Glu]", "NP_000079.2:p.Val600Glu"),
        ("NC_012920.1:m.[100A>G]", "NC_012920.1:m.100A>G"),
        ("NC_000001.11:o.[1000G>A]", "NC_000001.11:o.1000G>A"),
        // whole-entity / marker members
        ("NM_000088.3:c.[=]", "NM_000088.3:c.="),
        ("NM_000088.3:c.[?]", "NM_000088.3:c.?"),
        ("NP_000079.2:p.[=]", "NP_000079.2:p.="),
        ("NP_000079.2:p.[(?)]", "NP_000079.2:p.(?)"),
    ]
}

/// Conformant bracket forms that must NOT be flagged: ≥2 cis members or ≥2 trans groups.
fn conformant_brackets() -> Vec<&'static str> {
    vec![
        "NM_000088.3:c.[76A>C;88G>T]",
        "NM_000088.3:c.[76A>C];[88G>T]",
        "NM_000088.3:c.[2376G>C];[=]",
        "NP_000079.2:p.[Ser68Arg];[=]",
    ]
}

#[test]
fn strict_rejects_standalone_singletons() {
    for (input, _) in nonconformant_singletons() {
        let result = parse_hgvs_with_config(input, ErrorConfig::strict());
        assert!(
            result.is_err(),
            "strict mode should reject non-conformant singleton bracket: {input}"
        );
        let msg = format!("{}", result.err().unwrap());
        assert!(
            msg.contains(W),
            "strict rejection of {input} should cite {W}, got: {msg}"
        );
    }
}

#[test]
fn lenient_unwraps_and_warns() {
    for (input, expected) in nonconformant_singletons() {
        let parsed = parse_hgvs_with_config(input, ErrorConfig::lenient())
            .unwrap_or_else(|e| panic!("lenient should accept {input}: {e}"));
        assert!(
            !matches!(parsed.result, HgvsVariant::Allele(_)),
            "lenient should unwrap the singleton wrapper for {input}, still got Allele"
        );
        assert_eq!(
            parsed.result.to_string(),
            expected,
            "lenient canonical form mismatch for {input}"
        );
        assert!(
            parsed.warnings.iter().any(|w| w.error_type.code() == W),
            "lenient should emit {W} for {input}"
        );
    }
}

#[test]
fn silent_unwraps_without_warning() {
    for (input, expected) in nonconformant_singletons() {
        let parsed = parse_hgvs_with_config(input, ErrorConfig::silent())
            .unwrap_or_else(|e| panic!("silent should accept {input}: {e}"));
        assert!(
            !matches!(parsed.result, HgvsVariant::Allele(_)),
            "silent should unwrap the singleton wrapper for {input}"
        );
        assert_eq!(parsed.result.to_string(), expected);
        assert!(
            !parsed.warnings.iter().any(|w| w.error_type.code() == W),
            "silent must not emit {W} for {input}"
        );
    }
}

#[test]
fn conformant_brackets_are_not_flagged() {
    for input in conformant_brackets() {
        // Strict must accept (no cardinality rejection).
        let strict = parse_hgvs_with_config(input, ErrorConfig::strict());
        assert!(
            strict.is_ok(),
            "strict should accept conformant bracket {input}: {:?}",
            strict.err()
        );
        // The conformant bracket wrapper round-trips verbatim.
        assert_eq!(strict.unwrap().result.to_string(), input);
        // Lenient must not raise the cardinality warning.
        let lenient = parse_hgvs_with_config(input, ErrorConfig::lenient()).unwrap();
        assert!(
            !lenient.warnings.iter().any(|w| w.error_type.code() == W),
            "conformant bracket {input} must not raise {W}"
        );
    }
}

#[test]
fn empty_and_zero_brackets_still_error_in_all_modes() {
    // `c.[]` (empty) and `c.[0]` already fail to parse; the cardinality rule
    // does not regress that into silent acceptance.
    for input in ["NM_000088.3:c.[]", "NM_000088.3:c.[0]"] {
        for cfg in [
            ErrorConfig::strict(),
            ErrorConfig::lenient(),
            ErrorConfig::silent(),
        ] {
            assert!(
                parse_hgvs_with_config(input, cfg).is_err(),
                "{input} should remain a parse error"
            );
        }
    }
}
