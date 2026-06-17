//! Issue #278: heteroplasmy prose diagnostic + allele-fraction diagnostic code.
//!
//! Two follow-ups under #133:
//!
//! 1. ClinVar prose `m.<pos><ref>><alt>/<alt2>` is currently rejected
//!    with a generic "Unknown variant type prefix" / "Unexpected trailing
//!    characters". The dispatcher should emit a targeted message pointing
//!    at the three spec-supported alternatives:
//!    a) compound brackets `m.[a;b]`
//!    b) dual fully-qualified slash `var1/var2`
//!    c) spec compact mosaic `m.<pos>=/<alt>`
//!
//! 2. Allele-fraction annotations `[level=70%]`, `(80%)` are already
//!    rejected. The rejection must emit a dedicated SVA code so
//!    downstream tooling can recognize the intent rather than seeing a
//!    generic parse error.

use ferro_hgvs::error::ErrorCode;
use ferro_hgvs::parse_hgvs;

// ---------------------------------------------------------------------------
// Item 1: ClinVar prose multi-allelic shorthand diagnostic.
// ---------------------------------------------------------------------------

mod prose_diagnostic {
    use super::*;

    #[test]
    fn names_compact_mosaic() {
        let input = "NC_012920.1:m.3243A>G/T";
        let err = parse_hgvs(input).expect_err("prose multi-allelic must be rejected");
        let msg = err.to_string();
        assert!(
            msg.contains("compact mosaic") || msg.contains("=/"),
            "diagnostic should mention the spec compact mosaic alternative; got: {msg}"
        );
    }

    #[test]
    fn names_compound_brackets() {
        let input = "NC_012920.1:m.3243A>G/T";
        let err = parse_hgvs(input).expect_err("prose multi-allelic must be rejected");
        let msg = err.to_string();
        // The compound-bracket alternative names the coord type
        // explicitly, so look for `m.[` or the position-keyed form
        // `[3243`. Mere `[`/`]` would match too much (we already include
        // them in the position-marker advice).
        assert!(
            msg.contains("m.[") || msg.contains("[3243"),
            "diagnostic should mention the compound bracket alternative \
             with the m./position scaffold; got: {msg}"
        );
    }

    #[test]
    fn names_dual_qualified_slash() {
        let input = "NC_012920.1:m.3243A>G/C";
        let err = parse_hgvs(input).expect_err("prose multi-allelic must be rejected");
        let msg = err.to_string();
        // The diagnostic must reference a fully qualified RHS in some form.
        assert!(
            msg.contains("fully qualified") || msg.contains("fully-qualified"),
            "diagnostic should mention the dual fully-qualified slash alternative; got: {msg}"
        );
    }

    #[test]
    fn keys_on_alt_base() {
        // Sanity: the diagnostic should fire only for the prose shape
        // where the RHS of the slash is a bare nucleotide letter (no
        // edit, no accession, no `=`). Other malformed inputs should
        // not be mistaken for this case.
        let prose = parse_hgvs("NC_012920.1:m.3243A>G/T").expect_err("prose form must reject");
        assert!(
            prose.to_string().to_ascii_lowercase().contains("clinvar"),
            "the targeted diagnostic should self-identify (e.g. mention \
             ClinVar prose); got: {prose}"
        );
    }

    #[test]
    fn surfaces_structured_w3018_code() {
        let err = parse_hgvs("NC_012920.1:m.3243A>G/T").expect_err("prose form must reject");
        assert_eq!(
            err.code(),
            Some(ErrorCode::ClinVarProseMultiAllelic),
            "structured error code must be W3018 ClinVarProseMultiAllelic; got: {err:?}",
        );
        assert_eq!(
            err.code().map(|c| c.as_str()),
            Some("W3018".to_string()),
            "structured code must render as 'W3018'; got: {err:?}",
        );
    }

    #[test]
    fn legitimate_dual_fully_qualified_slash_parses() {
        // Positive boundary: a real dual fully-qualified slash (each
        // side is a complete `<acc>:m.<pos><edit>`) must NOT trip the
        // prose diagnostic — it is the recommended remediation, not
        // the malformed shape.
        let ok = "NC_012920.1:m.3243A>G/NC_012920.1:m.3243A>T";
        let parsed = parse_hgvs(ok);
        assert!(
            parsed.is_ok(),
            "legitimate dual fully-qualified slash must parse cleanly; got: {parsed:?}",
        );
    }
}

// ---------------------------------------------------------------------------
// Item 2: Allele-fraction diagnostic code (SVA).
// ---------------------------------------------------------------------------

mod allele_fraction {
    use super::*;

    #[test]
    fn bracket_form_emits_sva_code() {
        // [level=70%] and similar "level=XX%" inside square brackets
        // must be rejected with a targeted message that names the
        // offending annotation. The structured `Diagnostic.code`
        // channel carries the W-code so downstream tooling can
        // discriminate this rejection from a generic parse error.
        let input = "NC_012920.1:m.3243A>G[level=70%]";
        let err = parse_hgvs(input).expect_err("[level=70%] must be rejected");
        assert_eq!(
            err.code(),
            Some(ErrorCode::AlleleFractionAnnotation),
            "rejection must carry the dedicated SVA W-code on the structured \
             code channel (W3017 AlleleFractionAnnotation); got: {err:?}"
        );
        assert_eq!(err.code().map(|c| c.as_str()), Some("W3017".to_string()),);
        let msg = err.to_string();
        assert!(
            msg.contains("allele fraction") || msg.contains("heteroplasmy"),
            "rejection must self-identify (allele fraction / heteroplasmy); got: {msg}"
        );
    }

    #[test]
    fn paren_percent_emits_sva_code() {
        let input = "NC_012920.1:m.3243A>G(80%)";
        let err = parse_hgvs(input).expect_err("(80%) must be rejected");
        assert_eq!(
            err.code(),
            Some(ErrorCode::AlleleFractionAnnotation),
            "rejection must carry the dedicated SVA W-code on the structured \
             code channel (W3017 AlleleFractionAnnotation); got: {err:?}"
        );
        let msg = err.to_string();
        assert!(
            msg.contains("allele fraction") || msg.contains("heteroplasmy"),
            "rejection must self-identify (allele fraction / heteroplasmy); got: {msg}"
        );
    }

    #[test]
    fn diagnostic_points_to_metadata() {
        // The diagnostic should tell the user where the annotation
        // *should* live (downstream metadata: VCF FORMAT/AF, ClinVar
        // metadata fields, etc.), not just that it was rejected.
        let inputs = [
            "NC_012920.1:m.3243A>G[level=70%]",
            "NC_012920.1:m.3243A>G(80%)",
        ];
        for input in inputs {
            let err = parse_hgvs(input).expect_err("allele-fraction must reject");
            let msg = err.to_string();
            assert!(
                msg.to_ascii_lowercase().contains("metadata")
                    || msg.to_ascii_lowercase().contains("vcf"),
                "diagnostic should redirect callers to metadata; got `{input}`: {msg}"
            );
        }
    }

    #[test]
    fn w3017_is_registered() {
        // Registry-level pin: the SVA W-code is callable via `ferro explain`
        // and has heteroplasmy-relevant content.
        let info =
            ferro_hgvs::error_handling::get_code_info("W3017").expect("W3017 must be registered");
        assert_eq!(info.name, "AlleleFractionAnnotation");
        assert!(
            info.summary
                .to_ascii_lowercase()
                .contains("allele fraction")
                || info.summary.to_ascii_lowercase().contains("heteroplasmy"),
            "W3017 summary should mention allele fraction / heteroplasmy; got: {}",
            info.summary
        );
    }

    #[test]
    fn w3018_is_registered() {
        // Registry-level pin for the prose-shape SVA W-code.
        let info =
            ferro_hgvs::error_handling::get_code_info("W3018").expect("W3018 must be registered");
        assert_eq!(info.name, "ClinVarProseMultiAllelic");
        assert!(
            info.summary.to_ascii_lowercase().contains("clinvar")
                || info.summary.to_ascii_lowercase().contains("multi-allelic"),
            "W3018 summary should reference ClinVar / multi-allelic shape; got: {}",
            info.summary
        );
    }
}
