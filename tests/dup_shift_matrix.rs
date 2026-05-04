//! A6: 3' shift coverage matrix for dup
//!
//! 7 modules × 8 scenarios. Mirrors `tests/del_shift_matrix.rs` (#81 A5)
//! and `tests/ins_shift_matrix.rs` (#81 A1, A7) for completeness of the
//! ins/del/dup canonical-form trilogy.

mod common;

mod genomic {
    use super::common::synthetic::{normalize_to_string, SyntheticBuilder, PAD_OFFSET};
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
    fn no_shift_dup_base_differs_from_neighbor() {
        // Core: ACGTACGT. Dup 'G' at core pos 3. ref[del_start]='G',
        // ref[del_end]='T' → mismatch, no shift. dup_len=1, stays as dup.
        let p = provider("ACGTACGT");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}dup", &[3]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}dup", &[3]));
    }

    #[rstest]
    fn single_base_dup_in_homopolymer() {
        // Core: ACAAAACG. A-tract at core 3-6 (4 A's). Dup A at core 3.
        // Shifts to 3' end of tract; dup_len=1 stays as dup.
        let p = provider("ACAAAACG");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}dup", &[3]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}dup", &[6]));
    }

    #[rstest]
    // Case 1: AC tandem (ACAC at core 3-6). Dup 1 AC at core 3-4 → dup at core 5-6.
    #[case("ACACACGT", "g.{0}_{1}dup", &[3u64, 4u64], "g.{0}_{1}dup", &[5u64, 6u64])]
    // Case 2: GCA tandem (GCAGCAGCA at core 3-11). Dup 1 GCA at core 3-5 → dup at core 9-11.
    #[case("ACGCAGCAGCATG", "g.{0}_{1}dup", &[3u64, 5u64], "g.{0}_{1}dup", &[9u64, 11u64])]
    fn multi_base_dup_of_one_unit(
        #[case] core: &str,
        #[case] in_template: &str,
        #[case] in_args: &[u64],
        #[case] out_template: &str,
        #[case] out_args: &[u64],
    ) {
        let p = provider(core);
        let input = format!("NC_TEST.1:{}", hgvs(in_template, in_args));
        let expected = format!("NC_TEST.1:{}", hgvs(out_template, out_args));
        let result = normalize_to_string(p, &input);
        assert_eq!(result, expected);
    }

    #[rstest]
    // Case 1: 2 A's dup'd in 5-A homopolymer → A[7]. Pre-tract: 5 A's at 3-7.
    #[case("ACAAAAACG", "g.{0}_{1}dup", &[3u64, 4u64], "g.{0}_{1}A[7]", &[3u64, 7u64])]
    // Case 2: 4 A's dup'd in 4-A homopolymer → A[8] (large k).
    #[case("ACAAAACG", "g.{0}_{1}dup", &[3u64, 6u64], "g.{0}_{1}A[8]", &[3u64, 6u64])]
    // Case 3: 2 GCAs dup'd in 3-GCA tandem → GCA[5].
    #[case("ACGCAGCAGCATG", "g.{0}_{1}dup", &[3u64, 8u64], "g.{0}_{1}GCA[5]", &[3u64, 11u64])]
    // Case 4: cyclic-rotation case (CAGCAG dup'd, finds GCA tract by phase alignment).
    #[case("TTGCAGCAGCATT", "g.{0}_{1}dup", &[4u64, 9u64], "g.{0}_{1}GCA[5]", &[3u64, 11u64])]
    // Case 5: finer-periodicity (ATAT dup'd, smallest_repeat_unit returns AT, 5-AT tract → AT[7]).
    //   Core CCATATATATATCC: AT tract at core 3-12 (5 ATs), CC bookends.
    #[case("CCATATATATATCC", "g.{0}_{1}dup", &[3u64, 6u64], "g.{0}_{1}AT[7]", &[3u64, 12u64])]
    fn multi_base_dup_of_multiple_units(
        #[case] core: &str,
        #[case] in_template: &str,
        #[case] in_args: &[u64],
        #[case] out_template: &str,
        #[case] out_args: &[u64],
    ) {
        let p = provider(core);
        let input = format!("NC_TEST.1:{}", hgvs(in_template, in_args));
        let expected = format!("NC_TEST.1:{}", hgvs(out_template, out_args));
        let result = normalize_to_string(p, &input);
        assert_eq!(result, expected);
    }

    #[rstest]
    // Locks dual-path consistency: dup-form and repeat-form of the same edit
    // must canonicalize identically.
    // 4-A ref, 2 As added: c.3_4dup vs c.3_6A[6] both → c.3_6A[6].
    #[case("ACAAAACG", "g.{0}_{1}dup", &[3u64, 4u64], "g.{0}_{1}A[6]", &[3u64, 6u64])]
    fn round_trip_equivalence(
        #[case] core: &str,
        #[case] dup_template: &str,
        #[case] dup_args: &[u64],
        #[case] repeat_template: &str,
        #[case] repeat_args: &[u64],
    ) {
        let p1 = provider(core);
        let p2 = provider(core);
        let r1 = normalize_to_string(p1, &format!("NC_TEST.1:{}", hgvs(dup_template, dup_args)));
        let r2 = normalize_to_string(
            p2,
            &format!("NC_TEST.1:{}", hgvs(repeat_template, repeat_args)),
        );
        assert_eq!(r1, r2, "dup-form and repeat-form must agree");
    }

    #[rstest]
    fn multi_base_dup_shifts_but_not_tandem_unit_addition() {
        // Scenario 5: 2-unit tract, dup 1 unit → stays as dup shifted.
        // Core ACACGT: AC tandem at core 1-4 (2 ACs). Dup AC at core 1-2.
        // Shifts to core 3-4; tandem path has copies_in_dup=1, doesn't fire.
        let p = provider("ACACGT");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}_{1}dup", &[1, 2]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}_{1}dup", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        // 4-A tract bordered by C on both sides; dup shifts to last A.
        let p = provider("CAAAACG");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}dup", &[2]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}dup", &[5]));
    }

    #[rstest]
    fn dup_at_5prime_edge() {
        // Core AAAACG: A-tract at core 1-4. Dup A at core 1 shifts to core 4.
        let p = provider("AAAACG");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}dup", &[1]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}dup", &[4]));
    }

    #[rstest]
    fn dup_at_3prime_edge() {
        // Core CGAAAAC: A-tract at core 3-6. Dup A at core 6.
        // Trailing C bookend prevents the shift from leaking into PAD.
        let p = provider("CGAAAAC");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}dup", &[6]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}dup", &[6]));
    }
}

mod cds_plus {
    use super::common::synthetic::{normalize_to_string, SyntheticBuilder};
    use ferro_hgvs::reference::transcript::Strand;
    use rstest::rstest;

    fn provider(core: &str) -> ferro_hgvs::MockProvider {
        SyntheticBuilder::cds(core, 1, core.len() as u64, Strand::Plus).build()
    }

    fn hgvs(template: &str, args: &[u64]) -> String {
        let mut s = template.to_string();
        for (i, p) in args.iter().enumerate() {
            s = s.replace(&format!("{{{}}}", i), &p.to_string());
        }
        s
    }

    #[rstest]
    fn no_shift_dup_base_differs_from_neighbor() {
        let p = provider("ACGTACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}dup", &[3]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[3]));
    }

    #[rstest]
    fn single_base_dup_in_homopolymer() {
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}dup", &[3]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[6]));
    }

    #[rstest]
    #[case("ACACACGTACGTACGTACGT", "c.{0}_{1}dup", &[3u64, 4u64], "c.{0}_{1}dup", &[5u64, 6u64])]
    #[case("ACGCAGCAGCATGACGTACG", "c.{0}_{1}dup", &[3u64, 5u64], "c.{0}_{1}dup", &[9u64, 11u64])]
    fn multi_base_dup_of_one_unit(
        #[case] core: &str,
        #[case] in_template: &str,
        #[case] in_args: &[u64],
        #[case] out_template: &str,
        #[case] out_args: &[u64],
    ) {
        let p = provider(core);
        let input = format!("NM_TEST.1:{}", hgvs(in_template, in_args));
        let expected = format!("NM_TEST.1:{}", hgvs(out_template, out_args));
        let result = normalize_to_string(p, &input);
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case("ACAAAAACGTACGTACGTAC", "c.{0}_{1}dup", &[3u64, 4u64], "c.{0}_{1}A[7]", &[3u64, 7u64])]
    #[case("ACAAAACGTACGTACGTAC", "c.{0}_{1}dup", &[3u64, 6u64], "c.{0}_{1}A[8]", &[3u64, 6u64])]
    #[case("ACGCAGCAGCATGACGTACG", "c.{0}_{1}dup", &[3u64, 8u64], "c.{0}_{1}GCA[5]", &[3u64, 11u64])]
    #[case("TTGCAGCAGCATTACGTACGT", "c.{0}_{1}dup", &[4u64, 9u64], "c.{0}_{1}GCA[5]", &[3u64, 11u64])]
    #[case("CCATATATATATCCACGTACGT", "c.{0}_{1}dup", &[3u64, 6u64], "c.{0}_{1}AT[7]", &[3u64, 12u64])]
    fn multi_base_dup_of_multiple_units(
        #[case] core: &str,
        #[case] in_template: &str,
        #[case] in_args: &[u64],
        #[case] out_template: &str,
        #[case] out_args: &[u64],
    ) {
        let p = provider(core);
        let input = format!("NM_TEST.1:{}", hgvs(in_template, in_args));
        let expected = format!("NM_TEST.1:{}", hgvs(out_template, out_args));
        let result = normalize_to_string(p, &input);
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case("ACAAAACGTACGTACGTAC", "c.{0}_{1}dup", &[3u64, 4u64], "c.{0}_{1}A[6]", &[3u64, 6u64])]
    fn round_trip_equivalence(
        #[case] core: &str,
        #[case] dup_template: &str,
        #[case] dup_args: &[u64],
        #[case] repeat_template: &str,
        #[case] repeat_args: &[u64],
    ) {
        let p1 = provider(core);
        let p2 = provider(core);
        let r1 = normalize_to_string(p1, &format!("NM_TEST.1:{}", hgvs(dup_template, dup_args)));
        let r2 = normalize_to_string(
            p2,
            &format!("NM_TEST.1:{}", hgvs(repeat_template, repeat_args)),
        );
        assert_eq!(r1, r2, "dup-form and repeat-form must agree");
    }

    #[rstest]
    fn multi_base_dup_shifts_but_not_tandem_unit_addition() {
        let p = provider("ACACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}dup", &[1, 2]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}_{1}dup", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        let p = provider("CAAAACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}dup", &[2]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[5]));
    }

    #[rstest]
    fn dup_at_5prime_edge() {
        let p = provider("AAAACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}dup", &[1]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[4]));
    }

    #[rstest]
    fn dup_at_3prime_edge() {
        // A-tract at end of transcript; dup clamps at tx_end.
        let p = provider("CGTACGTACGTACGTAAAA");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}dup", &[19]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[19]));
    }
}

mod cds_minus {
    use super::common::synthetic::{normalize_to_string, SyntheticBuilder};
    use ferro_hgvs::reference::transcript::Strand;
    use rstest::rstest;

    fn provider(core: &str) -> ferro_hgvs::MockProvider {
        SyntheticBuilder::cds(core, 1, core.len() as u64, Strand::Minus).build()
    }

    fn hgvs(template: &str, args: &[u64]) -> String {
        let mut s = template.to_string();
        for (i, p) in args.iter().enumerate() {
            s = s.replace(&format!("{{{}}}", i), &p.to_string());
        }
        s
    }

    #[rstest]
    fn no_shift_dup_base_differs_from_neighbor() {
        let p = provider("ACGTACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}dup", &[3]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[3]));
    }

    #[rstest]
    fn single_base_dup_in_homopolymer() {
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}dup", &[3]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[6]));
    }

    #[rstest]
    #[case("ACACACGTACGTACGTACGT", "c.{0}_{1}dup", &[3u64, 4u64], "c.{0}_{1}dup", &[5u64, 6u64])]
    #[case("ACGCAGCAGCATGACGTACG", "c.{0}_{1}dup", &[3u64, 5u64], "c.{0}_{1}dup", &[9u64, 11u64])]
    fn multi_base_dup_of_one_unit(
        #[case] core: &str,
        #[case] in_template: &str,
        #[case] in_args: &[u64],
        #[case] out_template: &str,
        #[case] out_args: &[u64],
    ) {
        let p = provider(core);
        let input = format!("NM_TEST.1:{}", hgvs(in_template, in_args));
        let expected = format!("NM_TEST.1:{}", hgvs(out_template, out_args));
        let result = normalize_to_string(p, &input);
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case("ACAAAAACGTACGTACGTAC", "c.{0}_{1}dup", &[3u64, 4u64], "c.{0}_{1}A[7]", &[3u64, 7u64])]
    #[case("ACAAAACGTACGTACGTAC", "c.{0}_{1}dup", &[3u64, 6u64], "c.{0}_{1}A[8]", &[3u64, 6u64])]
    #[case("ACGCAGCAGCATGACGTACG", "c.{0}_{1}dup", &[3u64, 8u64], "c.{0}_{1}GCA[5]", &[3u64, 11u64])]
    #[case("TTGCAGCAGCATTACGTACGT", "c.{0}_{1}dup", &[4u64, 9u64], "c.{0}_{1}GCA[5]", &[3u64, 11u64])]
    #[case("CCATATATATATCCACGTACGT", "c.{0}_{1}dup", &[3u64, 6u64], "c.{0}_{1}AT[7]", &[3u64, 12u64])]
    fn multi_base_dup_of_multiple_units(
        #[case] core: &str,
        #[case] in_template: &str,
        #[case] in_args: &[u64],
        #[case] out_template: &str,
        #[case] out_args: &[u64],
    ) {
        let p = provider(core);
        let input = format!("NM_TEST.1:{}", hgvs(in_template, in_args));
        let expected = format!("NM_TEST.1:{}", hgvs(out_template, out_args));
        let result = normalize_to_string(p, &input);
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case("ACAAAACGTACGTACGTAC", "c.{0}_{1}dup", &[3u64, 4u64], "c.{0}_{1}A[6]", &[3u64, 6u64])]
    fn round_trip_equivalence(
        #[case] core: &str,
        #[case] dup_template: &str,
        #[case] dup_args: &[u64],
        #[case] repeat_template: &str,
        #[case] repeat_args: &[u64],
    ) {
        let p1 = provider(core);
        let p2 = provider(core);
        let r1 = normalize_to_string(p1, &format!("NM_TEST.1:{}", hgvs(dup_template, dup_args)));
        let r2 = normalize_to_string(
            p2,
            &format!("NM_TEST.1:{}", hgvs(repeat_template, repeat_args)),
        );
        assert_eq!(r1, r2, "dup-form and repeat-form must agree");
    }

    #[rstest]
    fn multi_base_dup_shifts_but_not_tandem_unit_addition() {
        let p = provider("ACACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}dup", &[1, 2]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}_{1}dup", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        let p = provider("CAAAACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}dup", &[2]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[5]));
    }

    #[rstest]
    fn dup_at_5prime_edge() {
        let p = provider("AAAACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}dup", &[1]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[4]));
    }

    #[rstest]
    fn dup_at_3prime_edge() {
        let p = provider("CGTACGTACGTACGTAAAA");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}dup", &[19]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[19]));
    }
}

mod noncoding_plus {
    use super::common::synthetic::{normalize_to_string, SyntheticBuilder};
    use ferro_hgvs::reference::transcript::Strand;
    use rstest::rstest;

    fn provider(core: &str) -> ferro_hgvs::MockProvider {
        SyntheticBuilder::noncoding(core, Strand::Plus).build()
    }

    fn hgvs(template: &str, args: &[u64]) -> String {
        let mut s = template.to_string();
        for (i, p) in args.iter().enumerate() {
            s = s.replace(&format!("{{{}}}", i), &p.to_string());
        }
        s
    }

    #[rstest]
    fn no_shift_dup_base_differs_from_neighbor() {
        let p = provider("ACGTACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}dup", &[3]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[3]));
    }

    #[rstest]
    fn single_base_dup_in_homopolymer() {
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}dup", &[3]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[6]));
    }

    #[rstest]
    #[case("ACACACGTACGTACGTACGT", "n.{0}_{1}dup", &[3u64, 4u64], "n.{0}_{1}dup", &[5u64, 6u64])]
    #[case("ACGCAGCAGCATGACGTACG", "n.{0}_{1}dup", &[3u64, 5u64], "n.{0}_{1}dup", &[9u64, 11u64])]
    fn multi_base_dup_of_one_unit(
        #[case] core: &str,
        #[case] in_template: &str,
        #[case] in_args: &[u64],
        #[case] out_template: &str,
        #[case] out_args: &[u64],
    ) {
        let p = provider(core);
        let input = format!("NR_TEST.1:{}", hgvs(in_template, in_args));
        let expected = format!("NR_TEST.1:{}", hgvs(out_template, out_args));
        let result = normalize_to_string(p, &input);
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case("ACAAAAACGTACGTACGTAC", "n.{0}_{1}dup", &[3u64, 4u64], "n.{0}_{1}A[7]", &[3u64, 7u64])]
    #[case("ACAAAACGTACGTACGTAC", "n.{0}_{1}dup", &[3u64, 6u64], "n.{0}_{1}A[8]", &[3u64, 6u64])]
    #[case("ACGCAGCAGCATGACGTACG", "n.{0}_{1}dup", &[3u64, 8u64], "n.{0}_{1}GCA[5]", &[3u64, 11u64])]
    #[case("TTGCAGCAGCATTACGTACGT", "n.{0}_{1}dup", &[4u64, 9u64], "n.{0}_{1}GCA[5]", &[3u64, 11u64])]
    #[case("CCATATATATATCCACGTACGT", "n.{0}_{1}dup", &[3u64, 6u64], "n.{0}_{1}AT[7]", &[3u64, 12u64])]
    fn multi_base_dup_of_multiple_units(
        #[case] core: &str,
        #[case] in_template: &str,
        #[case] in_args: &[u64],
        #[case] out_template: &str,
        #[case] out_args: &[u64],
    ) {
        let p = provider(core);
        let input = format!("NR_TEST.1:{}", hgvs(in_template, in_args));
        let expected = format!("NR_TEST.1:{}", hgvs(out_template, out_args));
        let result = normalize_to_string(p, &input);
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case("ACAAAACGTACGTACGTAC", "n.{0}_{1}dup", &[3u64, 4u64], "n.{0}_{1}A[6]", &[3u64, 6u64])]
    fn round_trip_equivalence(
        #[case] core: &str,
        #[case] dup_template: &str,
        #[case] dup_args: &[u64],
        #[case] repeat_template: &str,
        #[case] repeat_args: &[u64],
    ) {
        let p1 = provider(core);
        let p2 = provider(core);
        let r1 = normalize_to_string(p1, &format!("NR_TEST.1:{}", hgvs(dup_template, dup_args)));
        let r2 = normalize_to_string(
            p2,
            &format!("NR_TEST.1:{}", hgvs(repeat_template, repeat_args)),
        );
        assert_eq!(r1, r2, "dup-form and repeat-form must agree");
    }

    #[rstest]
    fn multi_base_dup_shifts_but_not_tandem_unit_addition() {
        let p = provider("ACACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}dup", &[1, 2]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}_{1}dup", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        let p = provider("CAAAACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}dup", &[2]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[5]));
    }

    #[rstest]
    fn dup_at_5prime_edge() {
        let p = provider("AAAACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}dup", &[1]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[4]));
    }

    #[rstest]
    fn dup_at_3prime_edge() {
        let p = provider("CGTACGTACGTACGTAAAA");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}dup", &[19]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[19]));
    }
}

mod noncoding_minus {
    use super::common::synthetic::{normalize_to_string, SyntheticBuilder};
    use ferro_hgvs::reference::transcript::Strand;
    use rstest::rstest;

    fn provider(core: &str) -> ferro_hgvs::MockProvider {
        SyntheticBuilder::noncoding(core, Strand::Minus).build()
    }

    fn hgvs(template: &str, args: &[u64]) -> String {
        let mut s = template.to_string();
        for (i, p) in args.iter().enumerate() {
            s = s.replace(&format!("{{{}}}", i), &p.to_string());
        }
        s
    }

    #[rstest]
    fn no_shift_dup_base_differs_from_neighbor() {
        let p = provider("ACGTACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}dup", &[3]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[3]));
    }

    #[rstest]
    fn single_base_dup_in_homopolymer() {
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}dup", &[3]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[6]));
    }

    #[rstest]
    #[case("ACACACGTACGTACGTACGT", "n.{0}_{1}dup", &[3u64, 4u64], "n.{0}_{1}dup", &[5u64, 6u64])]
    #[case("ACGCAGCAGCATGACGTACG", "n.{0}_{1}dup", &[3u64, 5u64], "n.{0}_{1}dup", &[9u64, 11u64])]
    fn multi_base_dup_of_one_unit(
        #[case] core: &str,
        #[case] in_template: &str,
        #[case] in_args: &[u64],
        #[case] out_template: &str,
        #[case] out_args: &[u64],
    ) {
        let p = provider(core);
        let input = format!("NR_TEST.1:{}", hgvs(in_template, in_args));
        let expected = format!("NR_TEST.1:{}", hgvs(out_template, out_args));
        let result = normalize_to_string(p, &input);
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case("ACAAAAACGTACGTACGTAC", "n.{0}_{1}dup", &[3u64, 4u64], "n.{0}_{1}A[7]", &[3u64, 7u64])]
    #[case("ACAAAACGTACGTACGTAC", "n.{0}_{1}dup", &[3u64, 6u64], "n.{0}_{1}A[8]", &[3u64, 6u64])]
    #[case("ACGCAGCAGCATGACGTACG", "n.{0}_{1}dup", &[3u64, 8u64], "n.{0}_{1}GCA[5]", &[3u64, 11u64])]
    #[case("TTGCAGCAGCATTACGTACGT", "n.{0}_{1}dup", &[4u64, 9u64], "n.{0}_{1}GCA[5]", &[3u64, 11u64])]
    #[case("CCATATATATATCCACGTACGT", "n.{0}_{1}dup", &[3u64, 6u64], "n.{0}_{1}AT[7]", &[3u64, 12u64])]
    fn multi_base_dup_of_multiple_units(
        #[case] core: &str,
        #[case] in_template: &str,
        #[case] in_args: &[u64],
        #[case] out_template: &str,
        #[case] out_args: &[u64],
    ) {
        let p = provider(core);
        let input = format!("NR_TEST.1:{}", hgvs(in_template, in_args));
        let expected = format!("NR_TEST.1:{}", hgvs(out_template, out_args));
        let result = normalize_to_string(p, &input);
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case("ACAAAACGTACGTACGTAC", "n.{0}_{1}dup", &[3u64, 4u64], "n.{0}_{1}A[6]", &[3u64, 6u64])]
    fn round_trip_equivalence(
        #[case] core: &str,
        #[case] dup_template: &str,
        #[case] dup_args: &[u64],
        #[case] repeat_template: &str,
        #[case] repeat_args: &[u64],
    ) {
        let p1 = provider(core);
        let p2 = provider(core);
        let r1 = normalize_to_string(p1, &format!("NR_TEST.1:{}", hgvs(dup_template, dup_args)));
        let r2 = normalize_to_string(
            p2,
            &format!("NR_TEST.1:{}", hgvs(repeat_template, repeat_args)),
        );
        assert_eq!(r1, r2, "dup-form and repeat-form must agree");
    }

    #[rstest]
    fn multi_base_dup_shifts_but_not_tandem_unit_addition() {
        let p = provider("ACACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}dup", &[1, 2]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}_{1}dup", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        let p = provider("CAAAACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}dup", &[2]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[5]));
    }

    #[rstest]
    fn dup_at_5prime_edge() {
        let p = provider("AAAACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}dup", &[1]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[4]));
    }

    #[rstest]
    fn dup_at_3prime_edge() {
        let p = provider("CGTACGTACGTACGTAAAA");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}dup", &[19]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[19]));
    }
}

mod rna_plus {
    use super::common::synthetic::{normalize_to_string, SyntheticBuilder};
    use ferro_hgvs::reference::transcript::Strand;
    use rstest::rstest;

    fn provider(core: &str) -> ferro_hgvs::MockProvider {
        SyntheticBuilder::rna(core, Strand::Plus).build()
    }

    fn hgvs(template: &str, args: &[u64]) -> String {
        let mut s = template.to_string();
        for (i, p) in args.iter().enumerate() {
            s = s.replace(&format!("{{{}}}", i), &p.to_string());
        }
        s
    }

    #[rstest]
    fn no_shift_dup_base_differs_from_neighbor() {
        let p = provider("acguacguacguacguacgu");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}dup", &[3]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}dup", &[3]));
    }

    #[rstest]
    fn single_base_dup_in_homopolymer() {
        let p = provider("acaaaacguacguacguac");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}dup", &[3]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}dup", &[6]));
    }

    #[rstest]
    #[case("acacacguacguacguacgu", "r.{0}_{1}dup", &[3u64, 4u64], "r.{0}_{1}dup", &[5u64, 6u64])]
    #[case("acgcagcagcaugacguacg", "r.{0}_{1}dup", &[3u64, 5u64], "r.{0}_{1}dup", &[9u64, 11u64])]
    fn multi_base_dup_of_one_unit(
        #[case] core: &str,
        #[case] in_template: &str,
        #[case] in_args: &[u64],
        #[case] out_template: &str,
        #[case] out_args: &[u64],
    ) {
        let p = provider(core);
        let input = format!("NR_TEST.1:{}", hgvs(in_template, in_args));
        let expected = format!("NR_TEST.1:{}", hgvs(out_template, out_args));
        let result = normalize_to_string(p, &input);
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case("acaaaaacguacguacguac", "r.{0}_{1}dup", &[3u64, 4u64], "r.{0}_{1}a[7]", &[3u64, 7u64])]
    #[case("acaaaacguacguacguac", "r.{0}_{1}dup", &[3u64, 6u64], "r.{0}_{1}a[8]", &[3u64, 6u64])]
    #[case("acgcagcagcaugacguacg", "r.{0}_{1}dup", &[3u64, 8u64], "r.{0}_{1}gca[5]", &[3u64, 11u64])]
    #[case("uugcagcagcauuacguacgu", "r.{0}_{1}dup", &[4u64, 9u64], "r.{0}_{1}gca[5]", &[3u64, 11u64])]
    // Finer-periodicity case uses an AC tandem rather than AT to avoid the
    // pre-existing limitation in the RNA Display path for repeat-form output
    // when the unit contains 'T' (DNA) — see comment in del_shift_matrix.rs.
    #[case("ccacacacacacccguacgu", "r.{0}_{1}dup", &[3u64, 6u64], "r.{0}_{1}ac[7]", &[3u64, 12u64])]
    fn multi_base_dup_of_multiple_units(
        #[case] core: &str,
        #[case] in_template: &str,
        #[case] in_args: &[u64],
        #[case] out_template: &str,
        #[case] out_args: &[u64],
    ) {
        let p = provider(core);
        let input = format!("NR_TEST.1:{}", hgvs(in_template, in_args));
        let expected = format!("NR_TEST.1:{}", hgvs(out_template, out_args));
        let result = normalize_to_string(p, &input);
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case("acaaaacguacguacguac", "r.{0}_{1}dup", &[3u64, 4u64], "r.{0}_{1}a[6]", &[3u64, 6u64])]
    fn round_trip_equivalence(
        #[case] core: &str,
        #[case] dup_template: &str,
        #[case] dup_args: &[u64],
        #[case] repeat_template: &str,
        #[case] repeat_args: &[u64],
    ) {
        let p1 = provider(core);
        let p2 = provider(core);
        let r1 = normalize_to_string(p1, &format!("NR_TEST.1:{}", hgvs(dup_template, dup_args)));
        let r2 = normalize_to_string(
            p2,
            &format!("NR_TEST.1:{}", hgvs(repeat_template, repeat_args)),
        );
        assert_eq!(r1, r2, "dup-form and repeat-form must agree");
    }

    #[rstest]
    fn multi_base_dup_shifts_but_not_tandem_unit_addition() {
        let p = provider("acacguacguacguacguac");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}_{1}dup", &[1, 2]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}_{1}dup", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        let p = provider("caaaacguacguacguacgu");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}dup", &[2]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}dup", &[5]));
    }

    #[rstest]
    fn dup_at_5prime_edge() {
        let p = provider("aaaacguacguacguacgu");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}dup", &[1]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}dup", &[4]));
    }

    #[rstest]
    fn dup_at_3prime_edge() {
        let p = provider("cguacguacguacgaaaa");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}dup", &[18]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}dup", &[18]));
    }
}

mod rna_minus {
    use super::common::synthetic::{normalize_to_string, SyntheticBuilder};
    use ferro_hgvs::reference::transcript::Strand;
    use rstest::rstest;

    fn provider(core: &str) -> ferro_hgvs::MockProvider {
        SyntheticBuilder::rna(core, Strand::Minus).build()
    }

    fn hgvs(template: &str, args: &[u64]) -> String {
        let mut s = template.to_string();
        for (i, p) in args.iter().enumerate() {
            s = s.replace(&format!("{{{}}}", i), &p.to_string());
        }
        s
    }

    #[rstest]
    fn no_shift_dup_base_differs_from_neighbor() {
        let p = provider("acguacguacguacguacgu");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}dup", &[3]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}dup", &[3]));
    }

    #[rstest]
    fn single_base_dup_in_homopolymer() {
        let p = provider("acaaaacguacguacguac");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}dup", &[3]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}dup", &[6]));
    }

    #[rstest]
    #[case("acacacguacguacguacgu", "r.{0}_{1}dup", &[3u64, 4u64], "r.{0}_{1}dup", &[5u64, 6u64])]
    #[case("acgcagcagcaugacguacg", "r.{0}_{1}dup", &[3u64, 5u64], "r.{0}_{1}dup", &[9u64, 11u64])]
    fn multi_base_dup_of_one_unit(
        #[case] core: &str,
        #[case] in_template: &str,
        #[case] in_args: &[u64],
        #[case] out_template: &str,
        #[case] out_args: &[u64],
    ) {
        let p = provider(core);
        let input = format!("NR_TEST.1:{}", hgvs(in_template, in_args));
        let expected = format!("NR_TEST.1:{}", hgvs(out_template, out_args));
        let result = normalize_to_string(p, &input);
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case("acaaaaacguacguacguac", "r.{0}_{1}dup", &[3u64, 4u64], "r.{0}_{1}a[7]", &[3u64, 7u64])]
    #[case("acaaaacguacguacguac", "r.{0}_{1}dup", &[3u64, 6u64], "r.{0}_{1}a[8]", &[3u64, 6u64])]
    #[case("acgcagcagcaugacguacg", "r.{0}_{1}dup", &[3u64, 8u64], "r.{0}_{1}gca[5]", &[3u64, 11u64])]
    #[case("uugcagcagcauuacguacgu", "r.{0}_{1}dup", &[4u64, 9u64], "r.{0}_{1}gca[5]", &[3u64, 11u64])]
    #[case("ccacacacacacccguacgu", "r.{0}_{1}dup", &[3u64, 6u64], "r.{0}_{1}ac[7]", &[3u64, 12u64])]
    fn multi_base_dup_of_multiple_units(
        #[case] core: &str,
        #[case] in_template: &str,
        #[case] in_args: &[u64],
        #[case] out_template: &str,
        #[case] out_args: &[u64],
    ) {
        let p = provider(core);
        let input = format!("NR_TEST.1:{}", hgvs(in_template, in_args));
        let expected = format!("NR_TEST.1:{}", hgvs(out_template, out_args));
        let result = normalize_to_string(p, &input);
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case("acaaaacguacguacguac", "r.{0}_{1}dup", &[3u64, 4u64], "r.{0}_{1}a[6]", &[3u64, 6u64])]
    fn round_trip_equivalence(
        #[case] core: &str,
        #[case] dup_template: &str,
        #[case] dup_args: &[u64],
        #[case] repeat_template: &str,
        #[case] repeat_args: &[u64],
    ) {
        let p1 = provider(core);
        let p2 = provider(core);
        let r1 = normalize_to_string(p1, &format!("NR_TEST.1:{}", hgvs(dup_template, dup_args)));
        let r2 = normalize_to_string(
            p2,
            &format!("NR_TEST.1:{}", hgvs(repeat_template, repeat_args)),
        );
        assert_eq!(r1, r2, "dup-form and repeat-form must agree");
    }

    #[rstest]
    fn multi_base_dup_shifts_but_not_tandem_unit_addition() {
        let p = provider("acacguacguacguacguac");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}_{1}dup", &[1, 2]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}_{1}dup", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        let p = provider("caaaacguacguacguacgu");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}dup", &[2]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}dup", &[5]));
    }

    #[rstest]
    fn dup_at_5prime_edge() {
        let p = provider("aaaacguacguacguacgu");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}dup", &[1]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}dup", &[4]));
    }

    #[rstest]
    fn dup_at_3prime_edge() {
        let p = provider("cguacguacguacgaaaa");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}dup", &[18]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}dup", &[18]));
    }
}
