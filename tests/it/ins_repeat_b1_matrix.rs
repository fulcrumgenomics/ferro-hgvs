//! B1: ins matching the adjacent repeat unit + c. codon-frame gate.
//!
//! Locks two HGVS spec rules across all coordinate systems:
//!
//! 1. 1 added copy of a unit = `dup`; ≥2 added copies = `[N]` repeat
//!    notation (or, in c. with `unit_len % 3 != 0`, the spec's fallback
//!    forms — see rule 2).
//! 2. In c. context, repeat notation requires `unit_len % 3 == 0`
//!    (codon-frame exception, `docs/recommendations/DNA/repeated.md`).
//!    When the gate fires:
//!      - contraction-with-survivors → plain `del`,
//!      - 1-copy expansion → `dup` (always permitted),
//!      - ≥2-copy expansion / dup → literal `ins<sequence>`.
//!
//! See spec at
//! `docs/superpowers/specs/2026-05-05-b1-ins-repeat-unit-increment-design.md`.

mod genomic {
    use crate::common::synthetic::{normalize_to_string, SyntheticBuilder, PAD_OFFSET};
    use rstest::rstest;

    /// 1-based HGVS position of the first base of the core region.
    const C0: u64 = PAD_OFFSET + 1;

    fn provider(core: &str) -> ferro_hgvs::MockProvider {
        SyntheticBuilder::genomic(core).build()
    }

    fn hgvs(template: &str, args: &[u64]) -> String {
        let mut s = template.to_string();
        for (i, p) in args.iter().enumerate() {
            s = s.replace(&format!("{{{}}}", i), &(C0 + p - 1).to_string());
        }
        s
    }

    #[rstest]
    fn single_copy_ins_in_homopolymer_tract_emits_dup() {
        // Core ACAAAACG: A-tract at core 3-6. insA at 2_3 → dup at the
        // canonical 3' end (1 added copy = dup per spec).
        let p = provider("ACAAAACG");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}_{1}insA", &[2, 3]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}dup", &[6]));
    }

    #[rstest]
    fn multi_copy_ins_homopolymer_unit_emits_repeat_in_genomic() {
        // Core ACAAAACG (4 A's). insAA at 2_3 → A[6] (gate is no-op in g.).
        let p = provider("ACAAAACG");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}_{1}insAA", &[2, 3]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}_{1}A[6]", &[3, 6]));
    }

    #[rstest]
    fn multi_copy_ins_dinucleotide_unit_emits_repeat_in_genomic() {
        // Core CATATATC: 3 ATs at core 2-7. insATAT at 1_2 → AT[5].
        let p = provider("CATATATC");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}_{1}insATAT", &[1, 2]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}_{1}AT[5]", &[2, 7]));
    }

    #[rstest]
    fn multi_copy_ins_codon_unit_emits_repeat_in_genomic() {
        // Core CCAGCAGCAGT: 3 CAGs at core 2-10. insCAGCAG at 1_2 → CAG[5].
        let p = provider("CCAGCAGCAGT");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}_{1}insCAGCAG", &[1, 2]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}_{1}CAG[5]", &[2, 10]));
    }

    #[rstest]
    fn single_copy_ins_into_single_ref_copy_emits_dup() {
        // Core ACGTACGT. insAC at 4_5 → adjacent matches, single copy → dup.
        let p = provider("ACGTACGT");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}_{1}insAC", &[4, 5]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}_{1}dup", &[5, 6]));
    }
}

mod cds_plus {
    use crate::common::synthetic::{hgvs, normalize_to_string, SyntheticBuilder};
    use ferro_hgvs::reference::transcript::Strand;
    use rstest::rstest;

    fn provider(core: &str) -> ferro_hgvs::MockProvider {
        SyntheticBuilder::cds(core, 1, core.len() as u64, Strand::Plus).build()
    }

    #[rstest]
    fn single_copy_ins_in_homopolymer_tract_emits_dup() {
        // 4-A tract; insA → dup at the canonical 3' end. 1 added copy is
        // always `dup` per spec, regardless of coord system.
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insA", &[2, 3]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[6]));
    }

    #[rstest]
    fn multi_copy_ins_homopolymer_unit_blocks_in_cds() {
        // unit_len=1 in c. → codon-frame gate forbids A[N]. Caller emits
        // ins literal at the 3' tract flanking position.
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insAA", &[2, 3]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}_{1}insAA", &[6, 7]));
    }

    #[rstest]
    fn multi_copy_ins_dinucleotide_unit_blocks_in_cds() {
        // unit_len=2 in c. → gate forbids AT[N]. Mirrors the spec's
        // c.1741_1742insTATATATA example.
        let p = provider("CATATATCACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insATAT", &[1, 2]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}_{1}insATAT", &[7, 8]));
    }

    #[rstest]
    fn single_copy_ins_dinucleotide_unit_in_multi_copy_tract_emits_dup() {
        // Single-copy alt of a codon-misaligned unit (len=2) abutting a
        // multi-copy AT tract. The 3'-rule fast path (`ref_count >= 2` in
        // `insertion_to_duplication`) emits `dup` at the tract-aligned 3'-most
        // position; the codon-frame gate does NOT block single-copy
        // additions, since the gate fires only when the alt itself is >=2
        // copies of a non-codon-aligned unit. Regression lock-in for the
        // codon_blocks_dup guard applied to the fast path.
        let p = provider("CATATATCACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insAT", &[1, 2]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}_{1}dup", &[6, 7]));
    }

    #[rstest]
    fn multi_copy_ins_codon_unit_passes_in_cds() {
        // unit_len=3, codon-aligned: gate is no-op, repeat notation fires.
        let p = provider("CCAGCAGCAGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insCAGCAG", &[1, 2]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}_{1}CAG[5]", &[2, 10]));
    }

    #[rstest]
    fn round_trip_dup_form_and_repeat_form_in_cds() {
        // Round-trip lock: dup-form and repeat-form of the same edit must
        // canonicalize identically under the codon-frame gate. 4-A ref,
        // c.3_4dup vs c.3_6A[6] both → `c.6_7insAA`.
        let p1 = provider("ACAAAACGTACGTACGTAC");
        let p2 = provider("ACAAAACGTACGTACGTAC");
        let r1 = normalize_to_string(p1, &hgvs("NM_TEST.1:c.{0}_{1}dup", &[3, 4]));
        let r2 = normalize_to_string(p2, &hgvs("NM_TEST.1:c.{0}_{1}A[6]", &[3, 6]));
        assert_eq!(r1, r2, "dup-form and repeat-form must agree");
        assert_eq!(r1, hgvs("NM_TEST.1:c.{0}_{1}insAA", &[6, 7]));
    }
}

mod cds_minus {
    use crate::common::synthetic::{hgvs, normalize_to_string, SyntheticBuilder};
    use ferro_hgvs::reference::transcript::Strand;
    use rstest::rstest;

    fn provider(core: &str) -> ferro_hgvs::MockProvider {
        SyntheticBuilder::cds(core, 1, core.len() as u64, Strand::Minus).build()
    }

    #[rstest]
    fn single_copy_ins_in_homopolymer_tract_emits_dup() {
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insA", &[2, 3]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[6]));
    }

    #[rstest]
    fn multi_copy_ins_homopolymer_unit_blocks_in_cds_minus() {
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insAA", &[2, 3]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}_{1}insAA", &[6, 7]));
    }

    #[rstest]
    fn multi_copy_ins_dinucleotide_unit_blocks_in_cds_minus() {
        // unit_len=2 in c. → gate forbids AT[N] regardless of strand.
        // Mirrors cds_plus::multi_copy_ins_dinucleotide_unit_blocks_in_cds.
        let p = provider("CATATATCACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insATAT", &[1, 2]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}_{1}insATAT", &[7, 8]));
    }

    #[rstest]
    fn multi_copy_ins_codon_unit_passes_in_cds_minus() {
        let p = provider("CCAGCAGCAGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insCAGCAG", &[1, 2]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}_{1}CAG[5]", &[2, 10]));
    }
}

mod noncoding_plus {
    use crate::common::synthetic::{hgvs, normalize_to_string, SyntheticBuilder};
    use ferro_hgvs::reference::transcript::Strand;
    use rstest::rstest;

    fn provider(core: &str) -> ferro_hgvs::MockProvider {
        SyntheticBuilder::noncoding(core, Strand::Plus).build()
    }

    #[rstest]
    fn multi_copy_ins_homopolymer_unit_emits_repeat_in_noncoding() {
        // n. is non-coding; codon-frame gate does not apply. A[N] permitted.
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}insAA", &[2, 3]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}_{1}A[6]", &[3, 6]));
    }

    #[rstest]
    fn single_copy_ins_in_homopolymer_tract_emits_dup_in_noncoding() {
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}insA", &[2, 3]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[6]));
    }
}

mod noncoding_minus {
    use crate::common::synthetic::{hgvs, normalize_to_string, SyntheticBuilder};
    use ferro_hgvs::reference::transcript::Strand;
    use rstest::rstest;

    fn provider(core: &str) -> ferro_hgvs::MockProvider {
        SyntheticBuilder::noncoding(core, Strand::Minus).build()
    }

    #[rstest]
    fn multi_copy_ins_homopolymer_unit_emits_repeat_in_noncoding_minus() {
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}insAA", &[2, 3]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}_{1}A[6]", &[3, 6]));
    }

    #[rstest]
    fn single_copy_ins_in_homopolymer_tract_emits_dup_in_noncoding_minus() {
        // Mirror of noncoding_plus single-copy case. 1 added copy is always
        // `dup` per spec, regardless of coord system or strand.
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}insA", &[2, 3]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[6]));
    }
}
