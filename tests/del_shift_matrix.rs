//! A5: 3' shift coverage matrix for del
//!
//! 7 modules × 8 scenarios. See spec at
//! docs/superpowers/specs/2026-05-02-A5-del-shift-coverage-matrix-design.md.

mod common;

mod genomic {
    use super::common::synthetic::{normalize_to_string, SyntheticBuilder, PAD_OFFSET};
    use rstest::rstest;

    /// 1-based HGVS position of the first base of the core region —
    /// derived from the synthetic builder's padding contract.
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
    fn no_shift_deleted_base_differs_from_neighbor() {
        // Core: ACGTACGT. Delete 'G' at core pos 3. ref[del_start]='G',
        // ref[del_end]='T' → mismatch, no shift.
        let p = provider("ACGTACGT");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}del", &[3]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}del", &[3]));
    }

    #[rstest]
    fn single_base_del_in_homopolymer() {
        // Core: ACAAAACG. A-tract at core pos 3-6 (4 A's). Delete A at core 3.
        // Shifts to 3' end; k=1 stays as del.
        let p = provider("ACAAAACG");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}del", &[3]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}del", &[6]));
    }

    #[rstest]
    // Case 1: AC tandem (ACAC at core 3-6). Delete 1 AC at core 3-4 → del at core 5-6.
    #[case("ACACACGT", "g.{0}_{1}del", &[3u64, 4u64], "g.{0}_{1}del", &[5u64, 6u64])]
    // Case 2: GCA tandem (GCAGCAGCA at core 3-11). Delete 1 GCA at core 3-5 → del at core 9-11.
    #[case("ACGCAGCAGCATG", "g.{0}_{1}del", &[3u64, 5u64], "g.{0}_{1}del", &[9u64, 11u64])]
    fn multi_base_del_in_tandem_one_unit(
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
    // Case 1: 2 A's removed from 5-A homopolymer → A[3].
    #[case("ACAAAAACG", "g.{0}_{1}del", &[3u64, 4u64], "g.{0}_{1}A[3]", &[3u64, 7u64])]
    // Case 2: 2 A's removed from 4-A homopolymer → A[2].
    #[case("ACAAAACG", "g.{0}_{1}del", &[3u64, 4u64], "g.{0}_{1}A[2]", &[3u64, 6u64])]
    // Case 3: 2 GCAs removed from 3-GCA tandem → GCA[1].
    #[case("ACGCAGCAGCATG", "g.{0}_{1}del", &[3u64, 8u64], "g.{0}_{1}GCA[1]", &[3u64, 11u64])]
    // Case 4: full tract removal → falls back to del.
    #[case("ACAAAACG", "g.{0}_{1}del", &[3u64, 6u64], "g.{0}_{1}del", &[3u64, 6u64])]
    // Case 5: cyclic rotation case (delete CAGCAG, find GCA tract by phase alignment).
    #[case("TTGCAGCAGCATT", "g.{0}_{1}del", &[4u64, 9u64], "g.{0}_{1}GCA[1]", &[3u64, 11u64])]
    // Case 6: finer periodicity (delete ATAT, smallest_repeat_unit returns AT).
    //   Core CCATATATATATCC: AT tract at core 3-12 (5 ATs), CC flanks. Bound the AT
    //   tract on both sides so the period doesn't leak into PAD.
    #[case("CCATATATATATCC", "g.{0}_{1}del", &[3u64, 6u64], "g.{0}_{1}AT[3]", &[3u64, 12u64])]
    fn multi_base_del_in_tandem_multiple_units(
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
    // Locks dual-path consistency: del-form and repeat-form of the same edit
    // must canonicalize identically.
    // 4-A ref, 2 As removed: c.3_4del vs c.3_6A[2] both → c.3_6A[2].
    #[case("ACAAAACG", "g.{0}_{1}del", &[3u64, 4u64], "g.{0}_{1}A[2]", &[3u64, 6u64])]
    fn round_trip_equivalence(
        #[case] core: &str,
        #[case] del_template: &str,
        #[case] del_args: &[u64],
        #[case] repeat_template: &str,
        #[case] repeat_args: &[u64],
    ) {
        let p1 = provider(core);
        let p2 = provider(core);
        let r1 = normalize_to_string(p1, &format!("NC_TEST.1:{}", hgvs(del_template, del_args)));
        let r2 = normalize_to_string(
            p2,
            &format!("NC_TEST.1:{}", hgvs(repeat_template, repeat_args)),
        );
        assert_eq!(r1, r2, "del-form and repeat-form must agree");
    }

    #[rstest]
    fn multi_base_del_shifts_but_not_tandem_unit_removal() {
        // Scenario 5: 2-unit tract, k=1 → stays as del shifted (no B2 rewrite).
        // Core ACACGT: AC tandem at core 1-4 (2 ACs). Delete AC at core 1-2.
        // Shifts to core 3-4; ref_count=2, k=1, B2 doesn't fire.
        let p = provider("ACACGT");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}_{1}del", &[1, 2]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}_{1}del", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        // 4-A tract bordered by C on both sides.
        // Core CAAAACG: A-tract at core 2-5 (4 A's). Delete A at core 2 shifts to core 5.
        let p = provider("CAAAACG");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}del", &[2]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}del", &[5]));
    }

    #[rstest]
    fn deletion_at_5prime_edge() {
        // Core AAAACG: A-tract at core 1-4. Delete A at core 1 shifts to core 4.
        let p = provider("AAAACG");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}del", &[1]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}del", &[4]));
    }

    #[rstest]
    fn deletion_at_3prime_edge() {
        // Core CGAAAAC: A-tract at core 3-6. Delete A at core 6.
        // Trailing C bookend prevents the shift from leaking into PAD.
        let p = provider("CGAAAAC");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}del", &[6]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}del", &[6]));
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
    fn no_shift_deleted_base_differs_from_neighbor() {
        // Core: ACGTACGTACGTACGTACGT (20 bp). c.3='G', c.4='T' → no shift.
        let p = provider("ACGTACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}del", &[3]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}del", &[3]));
    }

    #[rstest]
    fn single_base_del_in_homopolymer() {
        // Core: ACAAAACGTACGTACGTAC (19 bp). A-tract at c.3-6.
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}del", &[3]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}del", &[6]));
    }

    #[rstest]
    #[case("ACACACGTACGTACGTACGT", "c.{0}_{1}del", &[3u64, 4u64], "c.{0}_{1}del", &[5u64, 6u64])]
    #[case("ACGCAGCAGCATGACGTACG", "c.{0}_{1}del", &[3u64, 5u64], "c.{0}_{1}del", &[9u64, 11u64])]
    fn multi_base_del_in_tandem_one_unit(
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
    #[case("ACAAAAACGTACGTACGTAC", "c.{0}_{1}del", &[3u64, 4u64], "c.{0}_{1}A[3]", &[3u64, 7u64])]
    #[case("ACAAAACGTACGTACGTAC", "c.{0}_{1}del", &[3u64, 4u64], "c.{0}_{1}A[2]", &[3u64, 6u64])]
    #[case("ACGCAGCAGCATGACGTACG", "c.{0}_{1}del", &[3u64, 8u64], "c.{0}_{1}GCA[1]", &[3u64, 11u64])]
    #[case("ACAAAACGTACGTACGTAC", "c.{0}_{1}del", &[3u64, 6u64], "c.{0}_{1}del", &[3u64, 6u64])]
    #[case("TTGCAGCAGCATTACGTACGT", "c.{0}_{1}del", &[4u64, 9u64], "c.{0}_{1}GCA[1]", &[3u64, 11u64])]
    #[case("CCATATATATATCCACGTACGT", "c.{0}_{1}del", &[3u64, 6u64], "c.{0}_{1}AT[3]", &[3u64, 12u64])]
    fn multi_base_del_in_tandem_multiple_units(
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
    #[case("ACAAAACGTACGTACGTAC", "c.{0}_{1}del", &[3u64, 4u64], "c.{0}_{1}A[2]", &[3u64, 6u64])]
    fn round_trip_equivalence(
        #[case] core: &str,
        #[case] del_template: &str,
        #[case] del_args: &[u64],
        #[case] repeat_template: &str,
        #[case] repeat_args: &[u64],
    ) {
        let p1 = provider(core);
        let p2 = provider(core);
        let r1 = normalize_to_string(p1, &format!("NM_TEST.1:{}", hgvs(del_template, del_args)));
        let r2 = normalize_to_string(
            p2,
            &format!("NM_TEST.1:{}", hgvs(repeat_template, repeat_args)),
        );
        assert_eq!(r1, r2, "del-form and repeat-form must agree");
    }

    #[rstest]
    fn multi_base_del_shifts_but_not_tandem_unit_removal() {
        // 2-unit tract, k=1 → stays as del shifted (no B2).
        let p = provider("ACACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}del", &[1, 2]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}_{1}del", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        // 4-A tract bordered by C on both sides; del shifts to last A.
        let p = provider("CAAAACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}del", &[2]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}del", &[5]));
    }

    #[rstest]
    fn deletion_at_5prime_edge() {
        let p = provider("AAAACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}del", &[1]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}del", &[4]));
    }

    #[rstest]
    fn deletion_at_3prime_edge() {
        // A-tract at end of transcript; del clamps at tx_end.
        let p = provider("CGTACGTACGTACGTAAAA");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}del", &[19]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}del", &[19]));
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
    fn no_shift_deleted_base_differs_from_neighbor() {
        let p = provider("ACGTACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}del", &[3]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}del", &[3]));
    }

    #[rstest]
    fn single_base_del_in_homopolymer() {
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}del", &[3]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}del", &[6]));
    }

    #[rstest]
    #[case("ACACACGTACGTACGTACGT", "c.{0}_{1}del", &[3u64, 4u64], "c.{0}_{1}del", &[5u64, 6u64])]
    #[case("ACGCAGCAGCATGACGTACG", "c.{0}_{1}del", &[3u64, 5u64], "c.{0}_{1}del", &[9u64, 11u64])]
    fn multi_base_del_in_tandem_one_unit(
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
    #[case("ACAAAAACGTACGTACGTAC", "c.{0}_{1}del", &[3u64, 4u64], "c.{0}_{1}A[3]", &[3u64, 7u64])]
    #[case("ACAAAACGTACGTACGTAC", "c.{0}_{1}del", &[3u64, 4u64], "c.{0}_{1}A[2]", &[3u64, 6u64])]
    #[case("ACGCAGCAGCATGACGTACG", "c.{0}_{1}del", &[3u64, 8u64], "c.{0}_{1}GCA[1]", &[3u64, 11u64])]
    #[case("ACAAAACGTACGTACGTAC", "c.{0}_{1}del", &[3u64, 6u64], "c.{0}_{1}del", &[3u64, 6u64])]
    #[case("TTGCAGCAGCATTACGTACGT", "c.{0}_{1}del", &[4u64, 9u64], "c.{0}_{1}GCA[1]", &[3u64, 11u64])]
    #[case("CCATATATATATCCACGTACGT", "c.{0}_{1}del", &[3u64, 6u64], "c.{0}_{1}AT[3]", &[3u64, 12u64])]
    fn multi_base_del_in_tandem_multiple_units(
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
    #[case("ACAAAACGTACGTACGTAC", "c.{0}_{1}del", &[3u64, 4u64], "c.{0}_{1}A[2]", &[3u64, 6u64])]
    fn round_trip_equivalence(
        #[case] core: &str,
        #[case] del_template: &str,
        #[case] del_args: &[u64],
        #[case] repeat_template: &str,
        #[case] repeat_args: &[u64],
    ) {
        let p1 = provider(core);
        let p2 = provider(core);
        let r1 = normalize_to_string(p1, &format!("NM_TEST.1:{}", hgvs(del_template, del_args)));
        let r2 = normalize_to_string(
            p2,
            &format!("NM_TEST.1:{}", hgvs(repeat_template, repeat_args)),
        );
        assert_eq!(r1, r2, "del-form and repeat-form must agree");
    }

    #[rstest]
    fn multi_base_del_shifts_but_not_tandem_unit_removal() {
        let p = provider("ACACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}del", &[1, 2]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}_{1}del", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        let p = provider("CAAAACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}del", &[2]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}del", &[5]));
    }

    #[rstest]
    fn deletion_at_5prime_edge() {
        let p = provider("AAAACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}del", &[1]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}del", &[4]));
    }

    #[rstest]
    fn deletion_at_3prime_edge() {
        let p = provider("CGTACGTACGTACGTAAAA");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}del", &[19]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}del", &[19]));
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
    fn no_shift_deleted_base_differs_from_neighbor() {
        let p = provider("ACGTACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}del", &[3]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}del", &[3]));
    }

    #[rstest]
    fn single_base_del_in_homopolymer() {
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}del", &[3]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}del", &[6]));
    }

    #[rstest]
    #[case("ACACACGTACGTACGTACGT", "n.{0}_{1}del", &[3u64, 4u64], "n.{0}_{1}del", &[5u64, 6u64])]
    #[case("ACGCAGCAGCATGACGTACG", "n.{0}_{1}del", &[3u64, 5u64], "n.{0}_{1}del", &[9u64, 11u64])]
    fn multi_base_del_in_tandem_one_unit(
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
    #[case("ACAAAAACGTACGTACGTAC", "n.{0}_{1}del", &[3u64, 4u64], "n.{0}_{1}A[3]", &[3u64, 7u64])]
    #[case("ACAAAACGTACGTACGTAC", "n.{0}_{1}del", &[3u64, 4u64], "n.{0}_{1}A[2]", &[3u64, 6u64])]
    #[case("ACGCAGCAGCATGACGTACG", "n.{0}_{1}del", &[3u64, 8u64], "n.{0}_{1}GCA[1]", &[3u64, 11u64])]
    #[case("ACAAAACGTACGTACGTAC", "n.{0}_{1}del", &[3u64, 6u64], "n.{0}_{1}del", &[3u64, 6u64])]
    #[case("TTGCAGCAGCATTACGTACGT", "n.{0}_{1}del", &[4u64, 9u64], "n.{0}_{1}GCA[1]", &[3u64, 11u64])]
    #[case("CCATATATATATCCACGTACGT", "n.{0}_{1}del", &[3u64, 6u64], "n.{0}_{1}AT[3]", &[3u64, 12u64])]
    fn multi_base_del_in_tandem_multiple_units(
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
    #[case("ACAAAACGTACGTACGTAC", "n.{0}_{1}del", &[3u64, 4u64], "n.{0}_{1}A[2]", &[3u64, 6u64])]
    fn round_trip_equivalence(
        #[case] core: &str,
        #[case] del_template: &str,
        #[case] del_args: &[u64],
        #[case] repeat_template: &str,
        #[case] repeat_args: &[u64],
    ) {
        let p1 = provider(core);
        let p2 = provider(core);
        let r1 = normalize_to_string(p1, &format!("NR_TEST.1:{}", hgvs(del_template, del_args)));
        let r2 = normalize_to_string(
            p2,
            &format!("NR_TEST.1:{}", hgvs(repeat_template, repeat_args)),
        );
        assert_eq!(r1, r2, "del-form and repeat-form must agree");
    }

    #[rstest]
    fn multi_base_del_shifts_but_not_tandem_unit_removal() {
        let p = provider("ACACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}del", &[1, 2]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}_{1}del", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        let p = provider("CAAAACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}del", &[2]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}del", &[5]));
    }

    #[rstest]
    fn deletion_at_5prime_edge() {
        let p = provider("AAAACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}del", &[1]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}del", &[4]));
    }

    #[rstest]
    fn deletion_at_3prime_edge() {
        let p = provider("CGTACGTACGTACGTAAAA");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}del", &[19]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}del", &[19]));
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
    fn no_shift_deleted_base_differs_from_neighbor() {
        let p = provider("ACGTACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}del", &[3]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}del", &[3]));
    }

    #[rstest]
    fn single_base_del_in_homopolymer() {
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}del", &[3]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}del", &[6]));
    }

    #[rstest]
    #[case("ACACACGTACGTACGTACGT", "n.{0}_{1}del", &[3u64, 4u64], "n.{0}_{1}del", &[5u64, 6u64])]
    #[case("ACGCAGCAGCATGACGTACG", "n.{0}_{1}del", &[3u64, 5u64], "n.{0}_{1}del", &[9u64, 11u64])]
    fn multi_base_del_in_tandem_one_unit(
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
    #[case("ACAAAAACGTACGTACGTAC", "n.{0}_{1}del", &[3u64, 4u64], "n.{0}_{1}A[3]", &[3u64, 7u64])]
    #[case("ACAAAACGTACGTACGTAC", "n.{0}_{1}del", &[3u64, 4u64], "n.{0}_{1}A[2]", &[3u64, 6u64])]
    #[case("ACGCAGCAGCATGACGTACG", "n.{0}_{1}del", &[3u64, 8u64], "n.{0}_{1}GCA[1]", &[3u64, 11u64])]
    #[case("ACAAAACGTACGTACGTAC", "n.{0}_{1}del", &[3u64, 6u64], "n.{0}_{1}del", &[3u64, 6u64])]
    #[case("TTGCAGCAGCATTACGTACGT", "n.{0}_{1}del", &[4u64, 9u64], "n.{0}_{1}GCA[1]", &[3u64, 11u64])]
    #[case("CCATATATATATCCACGTACGT", "n.{0}_{1}del", &[3u64, 6u64], "n.{0}_{1}AT[3]", &[3u64, 12u64])]
    fn multi_base_del_in_tandem_multiple_units(
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
    #[case("ACAAAACGTACGTACGTAC", "n.{0}_{1}del", &[3u64, 4u64], "n.{0}_{1}A[2]", &[3u64, 6u64])]
    fn round_trip_equivalence(
        #[case] core: &str,
        #[case] del_template: &str,
        #[case] del_args: &[u64],
        #[case] repeat_template: &str,
        #[case] repeat_args: &[u64],
    ) {
        let p1 = provider(core);
        let p2 = provider(core);
        let r1 = normalize_to_string(p1, &format!("NR_TEST.1:{}", hgvs(del_template, del_args)));
        let r2 = normalize_to_string(
            p2,
            &format!("NR_TEST.1:{}", hgvs(repeat_template, repeat_args)),
        );
        assert_eq!(r1, r2, "del-form and repeat-form must agree");
    }

    #[rstest]
    fn multi_base_del_shifts_but_not_tandem_unit_removal() {
        let p = provider("ACACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}del", &[1, 2]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}_{1}del", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        let p = provider("CAAAACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}del", &[2]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}del", &[5]));
    }

    #[rstest]
    fn deletion_at_5prime_edge() {
        let p = provider("AAAACGTACGTACGTACGT");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}del", &[1]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}del", &[4]));
    }

    #[rstest]
    fn deletion_at_3prime_edge() {
        let p = provider("CGTACGTACGTACGTAAAA");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}del", &[19]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}del", &[19]));
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
    fn no_shift_deleted_base_differs_from_neighbor() {
        let p = provider("acguacguacguacguacgu");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}del", &[3]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}del", &[3]));
    }

    #[rstest]
    fn single_base_del_in_homopolymer() {
        let p = provider("acaaaacguacguacguac");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}del", &[3]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}del", &[6]));
    }

    #[rstest]
    #[case("acacacguacguacguacgu", "r.{0}_{1}del", &[3u64, 4u64], "r.{0}_{1}del", &[5u64, 6u64])]
    #[case("acgcagcagcaugacguacg", "r.{0}_{1}del", &[3u64, 5u64], "r.{0}_{1}del", &[9u64, 11u64])]
    fn multi_base_del_in_tandem_one_unit(
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
    #[case("acaaaaacguacguacguac", "r.{0}_{1}del", &[3u64, 4u64], "r.{0}_{1}a[3]", &[3u64, 7u64])]
    #[case("acaaaacguacguacguac", "r.{0}_{1}del", &[3u64, 4u64], "r.{0}_{1}a[2]", &[3u64, 6u64])]
    #[case("acgcagcagcaugacguacg", "r.{0}_{1}del", &[3u64, 8u64], "r.{0}_{1}gca[1]", &[3u64, 11u64])]
    #[case("acaaaacguacguacguac", "r.{0}_{1}del", &[3u64, 6u64], "r.{0}_{1}del", &[3u64, 6u64])]
    #[case("uugcagcagcauuacguacgu", "r.{0}_{1}del", &[4u64, 9u64], "r.{0}_{1}gca[1]", &[3u64, 11u64])]
    // Case 6 (finer periodicity) uses an AC tandem rather than AT to avoid
    // a pre-existing limitation in the RNA Display path for repeat-form output:
    // when the unit contains 'T' (DNA), the lowercase RNA emission keeps the
    // 't' instead of converting to 'u'. Filed as a follow-up. The genomic/cds/
    // noncoding modules use the AT tandem which is correct for DNA emission.
    #[case("ccacacacacacccguacgu", "r.{0}_{1}del", &[3u64, 6u64], "r.{0}_{1}ac[3]", &[3u64, 12u64])]
    fn multi_base_del_in_tandem_multiple_units(
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
    #[case("acaaaacguacguacguac", "r.{0}_{1}del", &[3u64, 4u64], "r.{0}_{1}a[2]", &[3u64, 6u64])]
    fn round_trip_equivalence(
        #[case] core: &str,
        #[case] del_template: &str,
        #[case] del_args: &[u64],
        #[case] repeat_template: &str,
        #[case] repeat_args: &[u64],
    ) {
        let p1 = provider(core);
        let p2 = provider(core);
        let r1 = normalize_to_string(p1, &format!("NR_TEST.1:{}", hgvs(del_template, del_args)));
        let r2 = normalize_to_string(
            p2,
            &format!("NR_TEST.1:{}", hgvs(repeat_template, repeat_args)),
        );
        assert_eq!(r1, r2, "del-form and repeat-form must agree");
    }

    #[rstest]
    fn multi_base_del_shifts_but_not_tandem_unit_removal() {
        let p = provider("acacguacguacguacguac");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}_{1}del", &[1, 2]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}_{1}del", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        let p = provider("caaaacguacguacguacgu");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}del", &[2]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}del", &[5]));
    }

    #[rstest]
    fn deletion_at_5prime_edge() {
        let p = provider("aaaacguacguacguacgu");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}del", &[1]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}del", &[4]));
    }

    #[rstest]
    fn deletion_at_3prime_edge() {
        let p = provider("cguacguacguacgaaaa");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}del", &[18]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}del", &[18]));
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
    fn no_shift_deleted_base_differs_from_neighbor() {
        let p = provider("acguacguacguacguacgu");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}del", &[3]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}del", &[3]));
    }

    #[rstest]
    fn single_base_del_in_homopolymer() {
        let p = provider("acaaaacguacguacguac");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}del", &[3]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}del", &[6]));
    }

    #[rstest]
    #[case("acacacguacguacguacgu", "r.{0}_{1}del", &[3u64, 4u64], "r.{0}_{1}del", &[5u64, 6u64])]
    #[case("acgcagcagcaugacguacg", "r.{0}_{1}del", &[3u64, 5u64], "r.{0}_{1}del", &[9u64, 11u64])]
    fn multi_base_del_in_tandem_one_unit(
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
    #[case("acaaaaacguacguacguac", "r.{0}_{1}del", &[3u64, 4u64], "r.{0}_{1}a[3]", &[3u64, 7u64])]
    #[case("acaaaacguacguacguac", "r.{0}_{1}del", &[3u64, 4u64], "r.{0}_{1}a[2]", &[3u64, 6u64])]
    #[case("acgcagcagcaugacguacg", "r.{0}_{1}del", &[3u64, 8u64], "r.{0}_{1}gca[1]", &[3u64, 11u64])]
    #[case("acaaaacguacguacguac", "r.{0}_{1}del", &[3u64, 6u64], "r.{0}_{1}del", &[3u64, 6u64])]
    #[case("uugcagcagcauuacguacgu", "r.{0}_{1}del", &[4u64, 9u64], "r.{0}_{1}gca[1]", &[3u64, 11u64])]
    // Case 6 (finer periodicity): see comment in rna_plus module.
    #[case("ccacacacacacccguacgu", "r.{0}_{1}del", &[3u64, 6u64], "r.{0}_{1}ac[3]", &[3u64, 12u64])]
    fn multi_base_del_in_tandem_multiple_units(
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
    #[case("acaaaacguacguacguac", "r.{0}_{1}del", &[3u64, 4u64], "r.{0}_{1}a[2]", &[3u64, 6u64])]
    fn round_trip_equivalence(
        #[case] core: &str,
        #[case] del_template: &str,
        #[case] del_args: &[u64],
        #[case] repeat_template: &str,
        #[case] repeat_args: &[u64],
    ) {
        let p1 = provider(core);
        let p2 = provider(core);
        let r1 = normalize_to_string(p1, &format!("NR_TEST.1:{}", hgvs(del_template, del_args)));
        let r2 = normalize_to_string(
            p2,
            &format!("NR_TEST.1:{}", hgvs(repeat_template, repeat_args)),
        );
        assert_eq!(r1, r2, "del-form and repeat-form must agree");
    }

    #[rstest]
    fn multi_base_del_shifts_but_not_tandem_unit_removal() {
        let p = provider("acacguacguacguacguac");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}_{1}del", &[1, 2]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}_{1}del", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        let p = provider("caaaacguacguacguacgu");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}del", &[2]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}del", &[5]));
    }

    #[rstest]
    fn deletion_at_5prime_edge() {
        let p = provider("aaaacguacguacguacgu");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}del", &[1]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}del", &[4]));
    }

    #[rstest]
    fn deletion_at_3prime_edge() {
        let p = provider("cguacguacguacgaaaa");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}del", &[18]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}del", &[18]));
    }
}
