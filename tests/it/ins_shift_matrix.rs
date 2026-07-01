//! A7: 3' shift coverage matrix for ins
//!
//! 56 rstest functions across 7 coord-system / strand modules × 8 shuffle
//! scenarios. See spec at:
//! docs/superpowers/specs/2026-05-01-A7-ins-shift-coverage-matrix-design.md.

mod genomic {
    use crate::common::synthetic::{normalize_to_string, SyntheticBuilder, PAD_OFFSET};
    use rstest::rstest;

    /// 1-based HGVS position of the first base of the core region —
    /// derived from the synthetic builder's padding contract.
    const C0: u64 = PAD_OFFSET + 1;

    /// Build a genomic provider over the contig `NC_TEST.1` with the given
    /// core sequence.
    fn provider(core: &str) -> ferro_hgvs::MockProvider {
        SyntheticBuilder::genomic(core).build()
    }

    /// Format an HGVS string by substituting `{n}` placeholders with positions
    /// computed as `C0 + p - 1` (so the caller can reason in core-local coords).
    fn hgvs(template: &str, args: &[u64]) -> String {
        let mut s = template.to_string();
        for (i, p) in args.iter().enumerate() {
            s = s.replace(&format!("{{{}}}", i), &(C0 + p - 1).to_string());
        }
        s
    }

    #[rstest]
    fn no_shift_first_base_mismatch() {
        // Core: ACGTACGT. Insert C between core pos 4 (T) and 5 (A).
        //   ref[5] = 'A' ≠ 'C' (alt[0])    → no 3' shift
        //   ref[4] = 'T' ≠ 'C' (alt[last]) → no preceding-base dup
        // Output: unchanged ins.
        let p = provider("ACGTACGT");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}_{1}insC", &[4, 5]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}_{1}insC", &[4, 5]));
    }

    #[rstest]
    fn single_base_ins_in_homopolymer_one_unit() {
        // Core: ACAAAACG. A-tract at core pos 3-6. Insert A at core pos 2_3.
        // 3'-shifts to end of A-tract (core pos 6), becomes dup.
        let p = provider("ACAAAACG");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}_{1}insA", &[2, 3]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}dup", &[6]));
    }

    #[rstest]
    // Issue #132 / #882: single-copy ins where alt is a non-zero cyclic
    // rotation of the reference repeat unit. Such a rotation is OUT OF PHASE
    // with the tract at the insertion cut: decoding the candidate dup would
    // yield a DIFFERENT sequence than the input insertion, so per HGVS
    // duplication.md:19 ("no evidence the extra copy is in tandem … it should
    // be described as an insertion") the #882 gate rejects the dup and the
    // variant stays a plain `ins` at the input position. (Contrast the
    // in-phase rotation in `multi_base_ins_in_tandem_one_unit`, which is a
    // legitimate dup.)
    //
    // Case 1: GT tract (GTGTGT at core 3..9), insert TG at core 2_3.
    //   alt "TG" is the r=1 (out-of-phase) rotation of canonical "GT"; the
    //   candidate dup GT[core 7_8] would decode differently, so it stays ins.
    #[case("ACGTGTGTAC", "g.{0}_{1}insTG", &[2u64, 3u64], "g.{0}_{1}insTG", &[2u64, 3u64])]
    // Case 2: GCA tract (GCAGCAGCA at core 3..11), insert CAG at core 2_3.
    //   alt "CAG" is the r=1 (out-of-phase) rotation of canonical "GCA"; the
    //   candidate dup would decode differently, so it stays ins.
    #[case("TTGCAGCAGCATT", "g.{0}_{1}insCAG", &[2u64, 3u64], "g.{0}_{1}insCAG", &[2u64, 3u64])]
    fn cyclic_rotation_ins_one_unit(
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
    // Multi-base ins in tandem repeat, ONE unit added → dup.
    // Case 1: AC tandem (ACAC at core pos 3-6), insert AC at core pos 2_3.
    //   Shifts to 3' end; dup at the most-3' AC = core pos 5_6.
    #[case("ACACACGT", "g.{0}_{1}insAC", &[2u64, 3u64], "g.{0}_{1}dup", &[5u64, 6u64])]
    // Case 2: GCA tandem (GCAGCAGCA at core pos 3-11), insert GCA at core pos 2_3.
    //   Shifts to 3' end; dup at the most-3' GCA = core pos 9_11.
    #[case("ACGCAGCAGCATG", "g.{0}_{1}insGCA", &[2u64, 3u64], "g.{0}_{1}dup", &[9u64, 11u64])]
    fn multi_base_ins_in_tandem_one_unit(
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
    // 2+ units added → repeat notation per HGVS spec.
    // Case 1: 2 A's added to AAA tract → A[5].
    //   Core: ACAAACG. ref core[3..6]=AAA. Insert AA at core 2_3.
    //   Shifts to end of A-tract; 3+2=5 As → A[5].
    //   Repeat tract: core[3..5] (1-based, 3 A's at HGVS pos 3..5).
    #[case("ACAAACG", "g.{0}_{1}insAA", &[2u64, 3u64], "g.{0}_{1}A[5]", &[3u64, 5u64])]
    // Case 2: 4 A's added to AAA tract → A[7].
    #[case("ACAAACG", "g.{0}_{1}insAAAA", &[2u64, 3u64], "g.{0}_{1}A[7]", &[3u64, 5u64])]
    // Case 3: 2 GCA units added to GCAGCAGCA tract → GCA[5].
    //   Core: ACGCAGCAGCATG. ref core[3..11]=GCAGCAGCA (3 GCA units).
    //   Insert GCAGCA at core 2_3. Shifts; 3+2=5 GCA → GCA[5].
    //   Repeat tract: core[3..11] (1-based, 9 bases = 3 units).
    #[case("ACGCAGCAGCATG", "g.{0}_{1}insGCAGCA", &[2u64, 3u64], "g.{0}_{1}GCA[5]", &[3u64, 11u64])]
    // Case 4: 2 AC units added to ACACAC tract → AC[5].
    //   Core: GCACACACG. Pairs (3-4)=AC, (5-6)=AC, (7-8)=AC; 3 AC units at core 3..8.
    //   Trailing G at core 9 prevents the AC tract from leaking into padding
    //   (padding starts with 'A', which would otherwise extend the tract).
    //   Insert ACAC at core 2_3 (between C and A). alt[0]=A, ref[3]=A → match → shifts.
    //   3+2=5 AC → AC[5]. Tract: core 3..8 (1-based, 6 bases = 3 units).
    #[case("GCACACACG", "g.{0}_{1}insACAC", &[2u64, 3u64], "g.{0}_{1}AC[5]", &[3u64, 8u64])]
    // Case 5: alt is an OUT-OF-PHASE cyclic rotation of the ref unit (#882).
    //   Core: TTGCAGCAGCATT. Ref tract (GCA)x3 at core 3..11; flanking T's
    //   at core 1-2 and 12-13 prevent the tract from extending into PAD.
    //   Insert CAGCAG at core 2_3 (between T and the start of the tract).
    //   The candidate repeat GCA[5] at core 3..11 would decode to a DIFFERENT
    //   sequence than inserting CAGCAG at the cut (the rotation is out of
    //   phase), so the #882 gate rejects it and the variant stays a plain
    //   `ins` at the input position (duplication.md:19).
    #[case("TTGCAGCAGCATT", "g.{0}_{1}insCAGCAG", &[2u64, 3u64], "g.{0}_{1}insCAGCAG", &[2u64, 3u64])]
    fn multi_base_ins_in_tandem_multiple_units(
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
    fn partial_rotation_no_full_match() {
        // Multi-base ins shifts by 1 (rotation = 1, len = 2). Rotated alt
        // doesn't match adjacent ref, so output is rotated ins, not dup.
        //
        // Layout: core "ACGTAC". Insert GA at core 2_3:
        //   ref[3]=G=alt[0]            → shift, alt rotates to "AG"
        //   ref[4]=T, alt[1]=A         → mismatch, stop
        // Final position: core 3_4, alt rotated to "AG".
        // Preceding ref bytes (core[1..3]="CG") ≠ "AG" → not dup.
        let p = provider("ACGTAC");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}_{1}insGA", &[2, 3]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}_{1}insAG", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        // Genomic-only "boundary" is the end of a homopolymer tract — once
        // the run-of-A's ends, shuffling stops. (Pure end-of-ref clamping
        // is exercised in Task 8 onward via transcript boundaries.)
        //
        // Core: CAAAACG. A-tract at core pos 2..5 (1-based), bordered by
        // C at pos 1 and C at pos 6. Insert A at core 2_3 — shifts to end
        // of A-tract. With a 1-unit insertion (single A), output is dup
        // at core pos 5 (last A before boundary).
        let p = provider("CAAAACG");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}_{1}insA", &[2, 3]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}dup", &[5]));
    }

    #[rstest]
    fn insertion_at_5prime_edge() {
        // Insertion at very start of core: pos 1_2. Verify shuffle behavior
        // doesn't underflow on the 5' boundary.
        //
        // Core: AAAACG. A-tract at core 1..4. Insert A at core 1_2.
        // Shifts to end of A-tract; dup at core pos 4.
        let p = provider("AAAACG");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}_{1}insA", &[1, 2]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}dup", &[4]));
    }

    #[rstest]
    fn insertion_at_3prime_edge() {
        // Insertion near the 3' edge of the A-tract.
        //
        // Core: CGAAAAC. A-tract at core 3..6, bordered by G at pos 2 and
        // C at pos 7. Insert A at core 5_6 (within the tract). Shifts to
        // end (core pos 6); dup at core pos 6.
        let p = provider("CGAAAAC");
        let result = normalize_to_string(p, &hgvs("NC_TEST.1:g.{0}_{1}insA", &[5, 6]));
        assert_eq!(result, hgvs("NC_TEST.1:g.{0}dup", &[6]));
    }
}
mod cds_plus {
    use crate::common::synthetic::{normalize_to_string, SyntheticBuilder};
    use ferro_hgvs::reference::transcript::Strand;
    use rstest::rstest;

    /// Build a plus-strand CDS provider whose CDS spans the entire transcript.
    fn provider(core: &str) -> ferro_hgvs::MockProvider {
        SyntheticBuilder::cds(core, 1, core.len() as u64, Strand::Plus).build()
    }

    /// Format an HGVS string with `{n}` placeholders. CDS positions are
    /// 1-based transcript-relative, no PAD_OFFSET — the caller passes raw
    /// 1-based positions.
    fn hgvs(template: &str, args: &[u64]) -> String {
        let mut s = template.to_string();
        for (i, p) in args.iter().enumerate() {
            s = s.replace(&format!("{{{}}}", i), &p.to_string());
        }
        s
    }

    #[rstest]
    fn no_shift_first_base_mismatch() {
        // Core: ACGTACGTACGTACGTACGTACGTACGTAC (30 bases). c.10=C, c.11=A.
        // Insert C between c.10 and c.11 — alt[0]='C' ≠ ref[c.11]='A' (no shift)
        // and alt[last]='C' ≠ ref[c.10]='C'... wait that's equal, would dup.
        // Pick a different position: c.4_5. c.4=T, c.5=A. Insert C.
        //   alt[0]='C' ≠ ref[c.5]='A' → no shift
        //   alt[last]='C' ≠ ref[c.4]='T' → no preceding-base dup
        let p = provider("ACGTACGTACGTACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insC", &[4, 5]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}_{1}insC", &[4, 5]));
    }

    #[rstest]
    fn single_base_ins_in_homopolymer_one_unit() {
        // Core: ACAAAACGTACGTACGTAC (19 bases). A-tract at c.3..6.
        // Insert A at c.2_3. Shifts to end of A-tract; dup at c.6.
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insA", &[2, 3]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[6]));
    }

    #[rstest]
    // Case 1: AC tandem at c.3-6 (ACAC). Insert AC at c.2_3 → dup at c.5_6.
    #[case("ACACACGTACGTACGTACGTAC", "c.{0}_{1}insAC", &[2u64, 3u64], "c.{0}_{1}dup", &[5u64, 6u64])]
    // Case 2: GCA tandem at c.3-11 (GCAGCAGCA). Insert GCA at c.2_3 → dup at c.9_11.
    #[case("ACGCAGCAGCATGACGTACG", "c.{0}_{1}insGCA", &[2u64, 3u64], "c.{0}_{1}dup", &[9u64, 11u64])]
    fn multi_base_ins_in_tandem_one_unit(
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
    // Issue #132 / #882: single-copy ins where alt is an OUT-OF-PHASE cyclic
    // rotation of the ref repeat unit. The candidate dup would decode to a
    // different sequence than the input insertion, so the #882 gate rejects it
    // and the variant stays a plain `ins` at the input position
    // (duplication.md:19). Mirrors the genomic module case shape.
    //
    // Case 1: GT tract at c.3-8. Insert TG (out-of-phase r=1 of GT) at c.2_3 → stays ins.
    #[case("ACGTGTGTACGTACGTACGT", "c.{0}_{1}insTG", &[2u64, 3u64], "c.{0}_{1}insTG", &[2u64, 3u64])]
    // Case 2: GCA tract at c.3-11. Insert CAG (out-of-phase r=1 of GCA) at c.2_3 → stays ins.
    #[case("TTGCAGCAGCATTGACGTAC", "c.{0}_{1}insCAG", &[2u64, 3u64], "c.{0}_{1}insCAG", &[2u64, 3u64])]
    fn cyclic_rotation_ins_one_unit(
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
    // Case 1: 2 A's added to AAA tract. unit_len=1; codon-frame gate forbids
    // A[N] in c., emit insAA at the 3' tract flanking position per spec.
    #[case("ACAAACGTACGTACGTACGT", "c.{0}_{1}insAA", &[2u64, 3u64], "c.{0}_{1}insAA", &[5u64, 6u64])]
    // Case 2: 4 A's added to AAA tract → insAAAA (gate blocks A[N], unit_len=1).
    #[case("ACAAACGTACGTACGTACGT", "c.{0}_{1}insAAAA", &[2u64, 3u64], "c.{0}_{1}insAAAA", &[5u64, 6u64])]
    // Case 3: 2 GCA units added → GCA[5] (codon-aligned, gate passes).
    #[case("ACGCAGCAGCATGACGTACG", "c.{0}_{1}insGCAGCA", &[2u64, 3u64], "c.{0}_{1}GCA[5]", &[3u64, 11u64])]
    // Case 4: 2 AC units added → insACAC (gate blocks AC[N], unit_len=2).
    #[case("GCACACACGTACGTACGTAC", "c.{0}_{1}insACAC", &[2u64, 3u64], "c.{0}_{1}insACAC", &[8u64, 9u64])]
    fn multi_base_ins_in_tandem_multiple_units(
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
    fn partial_rotation_no_full_match() {
        // Core: ACGTACGTACGTACGTACGTAC. Insert GA at c.2_3.
        //   ref[c.3]=G=alt[0] → shift, alt rotates to "AG"
        //   ref[c.4]=T, alt[1]=A → mismatch, stop.
        // Output: c.3_4insAG.
        let p = provider("ACGTACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insGA", &[2, 3]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}_{1}insAG", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        // Boundary clamp via end-of-transcript: A-tract runs into the last
        // position of the transcript. Shuffling stops at the transcript end.
        //
        // Core: CAAAA (5 bases — entire transcript). A-tract c.2..5; the
        // transcript ENDS at c.5 (end-of-ref boundary). Insert A at c.2_3.
        // Should clamp to c.5 — dup at c.5.
        let p = provider("CAAAA");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insA", &[2, 3]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[5]));
    }

    #[rstest]
    fn insertion_at_5prime_edge() {
        // Insertion at c.1_2. Core: AAAACGTACGTACGTACGTAC.
        // A-tract at c.1..4. Shifts to dup at c.4.
        let p = provider("AAAACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insA", &[1, 2]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[4]));
    }

    #[rstest]
    fn insertion_at_3prime_edge() {
        // Insertion within an A-tract near transcript end.
        // Core: CGTACGTACGTACGAAAAC (19 bases). A-tract at c.15..18.
        // Insert A at c.16_17 (within tract); shifts to c.18, dup at c.18.
        let p = provider("CGTACGTACGTACGAAAAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insA", &[16, 17]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[18]));
    }
}
mod cds_minus {
    use crate::common::synthetic::{normalize_to_string, SyntheticBuilder};
    use ferro_hgvs::reference::transcript::Strand;
    use rstest::rstest;

    /// Build a minus-strand CDS provider whose CDS spans the entire transcript.
    fn provider(core: &str) -> ferro_hgvs::MockProvider {
        SyntheticBuilder::cds(core, 1, core.len() as u64, Strand::Minus).build()
    }

    /// Format an HGVS string with `{n}` placeholders. CDS positions are
    /// 1-based transcript-relative, no PAD_OFFSET — the caller passes raw
    /// 1-based positions.
    fn hgvs(template: &str, args: &[u64]) -> String {
        let mut s = template.to_string();
        for (i, p) in args.iter().enumerate() {
            s = s.replace(&format!("{{{}}}", i), &p.to_string());
        }
        s
    }

    #[rstest]
    fn no_shift_first_base_mismatch() {
        // Core: ACGTACGTACGTACGTACGTACGTACGTAC (30 bases). c.10=C, c.11=A.
        // Insert C between c.10 and c.11 — alt[0]='C' ≠ ref[c.11]='A' (no shift)
        // and alt[last]='C' ≠ ref[c.10]='C'... wait that's equal, would dup.
        // Pick a different position: c.4_5. c.4=T, c.5=A. Insert C.
        //   alt[0]='C' ≠ ref[c.5]='A' → no shift
        //   alt[last]='C' ≠ ref[c.4]='T' → no preceding-base dup
        let p = provider("ACGTACGTACGTACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insC", &[4, 5]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}_{1}insC", &[4, 5]));
    }

    #[rstest]
    fn single_base_ins_in_homopolymer_one_unit() {
        // Core: ACAAAACGTACGTACGTAC (19 bases). A-tract at c.3..6.
        // Insert A at c.2_3. Shifts to end of A-tract; dup at c.6.
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insA", &[2, 3]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[6]));
    }

    #[rstest]
    // Case 1: AC tandem at c.3-6 (ACAC). Insert AC at c.2_3 → dup at c.5_6.
    #[case("ACACACGTACGTACGTACGTAC", "c.{0}_{1}insAC", &[2u64, 3u64], "c.{0}_{1}dup", &[5u64, 6u64])]
    // Case 2: GCA tandem at c.3-11 (GCAGCAGCA). Insert GCA at c.2_3 → dup at c.9_11.
    #[case("ACGCAGCAGCATGACGTACG", "c.{0}_{1}insGCA", &[2u64, 3u64], "c.{0}_{1}dup", &[9u64, 11u64])]
    fn multi_base_ins_in_tandem_one_unit(
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
    // Issue #132 / #882: single-copy ins where alt is an OUT-OF-PHASE cyclic
    // rotation of the ref repeat unit. The candidate dup would decode to a
    // different sequence than the input insertion, so the #882 gate rejects it
    // and the variant stays a plain `ins` at the input position
    // (duplication.md:19). Mirrors the genomic module case shape.
    //
    // Case 1: GT tract at c.3-8. Insert TG (out-of-phase r=1 of GT) at c.2_3 → stays ins.
    #[case("ACGTGTGTACGTACGTACGT", "c.{0}_{1}insTG", &[2u64, 3u64], "c.{0}_{1}insTG", &[2u64, 3u64])]
    // Case 2: GCA tract at c.3-11. Insert CAG (out-of-phase r=1 of GCA) at c.2_3 → stays ins.
    #[case("TTGCAGCAGCATTGACGTAC", "c.{0}_{1}insCAG", &[2u64, 3u64], "c.{0}_{1}insCAG", &[2u64, 3u64])]
    fn cyclic_rotation_ins_one_unit(
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
    // Case 1: 2 A's added to AAA tract. unit_len=1; codon-frame gate forbids
    // A[N] in c., emit insAA at the 3' tract flanking position per spec.
    #[case("ACAAACGTACGTACGTACGT", "c.{0}_{1}insAA", &[2u64, 3u64], "c.{0}_{1}insAA", &[5u64, 6u64])]
    // Case 2: 4 A's added to AAA tract → insAAAA (gate blocks A[N], unit_len=1).
    #[case("ACAAACGTACGTACGTACGT", "c.{0}_{1}insAAAA", &[2u64, 3u64], "c.{0}_{1}insAAAA", &[5u64, 6u64])]
    // Case 3: 2 GCA units added → GCA[5] (codon-aligned, gate passes).
    #[case("ACGCAGCAGCATGACGTACG", "c.{0}_{1}insGCAGCA", &[2u64, 3u64], "c.{0}_{1}GCA[5]", &[3u64, 11u64])]
    // Case 4: 2 AC units added → insACAC (gate blocks AC[N], unit_len=2).
    #[case("GCACACACGTACGTACGTAC", "c.{0}_{1}insACAC", &[2u64, 3u64], "c.{0}_{1}insACAC", &[8u64, 9u64])]
    fn multi_base_ins_in_tandem_multiple_units(
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
    fn partial_rotation_no_full_match() {
        // Core: ACGTACGTACGTACGTACGTAC. Insert GA at c.2_3.
        //   ref[c.3]=G=alt[0] → shift, alt rotates to "AG"
        //   ref[c.4]=T, alt[1]=A → mismatch, stop.
        // Output: c.3_4insAG.
        let p = provider("ACGTACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insGA", &[2, 3]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}_{1}insAG", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        // Boundary clamp via end-of-transcript: A-tract runs into the last
        // position of the transcript. Shuffling stops at the transcript end.
        //
        // Core: CAAAA (5 bases — entire transcript). A-tract c.2..5; the
        // transcript ENDS at c.5 (end-of-ref boundary). Insert A at c.2_3.
        // Should clamp to c.5 — dup at c.5.
        let p = provider("CAAAA");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insA", &[2, 3]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[5]));
    }

    #[rstest]
    fn insertion_at_5prime_edge() {
        // Insertion at c.1_2. Core: AAAACGTACGTACGTACGTAC.
        // A-tract at c.1..4. Shifts to dup at c.4.
        let p = provider("AAAACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insA", &[1, 2]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[4]));
    }

    #[rstest]
    fn insertion_at_3prime_edge() {
        // Insertion within an A-tract near transcript end.
        // Core: CGTACGTACGTACGAAAAC (19 bases). A-tract at c.15..18.
        // Insert A at c.16_17 (within tract); shifts to c.18, dup at c.18.
        let p = provider("CGTACGTACGTACGAAAAC");
        let result = normalize_to_string(p, &hgvs("NM_TEST.1:c.{0}_{1}insA", &[16, 17]));
        assert_eq!(result, hgvs("NM_TEST.1:c.{0}dup", &[18]));
    }
}
mod noncoding_plus {
    use crate::common::synthetic::{normalize_to_string, SyntheticBuilder};
    use ferro_hgvs::reference::transcript::Strand;
    use rstest::rstest;

    /// Build a plus-strand noncoding provider whose transcript spans the entire
    /// core sequence.
    fn provider(core: &str) -> ferro_hgvs::MockProvider {
        SyntheticBuilder::noncoding(core, Strand::Plus).build()
    }

    /// Format an HGVS string with `{n}` placeholders. Noncoding positions are
    /// 1-based transcript-relative, no PAD_OFFSET — the caller passes raw
    /// 1-based positions.
    fn hgvs(template: &str, args: &[u64]) -> String {
        let mut s = template.to_string();
        for (i, p) in args.iter().enumerate() {
            s = s.replace(&format!("{{{}}}", i), &p.to_string());
        }
        s
    }

    #[rstest]
    fn no_shift_first_base_mismatch() {
        // Core: ACGTACGTACGTACGTACGTACGTACGTAC (30 bases). n.10=C, n.11=A.
        // Insert C between n.10 and n.11 — alt[0]='C' ≠ ref[n.11]='A' (no shift)
        // and alt[last]='C' ≠ ref[n.10]='C'... wait that's equal, would dup.
        // Pick a different position: n.4_5. n.4=T, n.5=A. Insert C.
        //   alt[0]='C' ≠ ref[n.5]='A' → no shift
        //   alt[last]='C' ≠ ref[n.4]='T' → no preceding-base dup
        let p = provider("ACGTACGTACGTACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}insC", &[4, 5]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}_{1}insC", &[4, 5]));
    }

    #[rstest]
    fn single_base_ins_in_homopolymer_one_unit() {
        // Core: ACAAAACGTACGTACGTAC (19 bases). A-tract at n.3..6.
        // Insert A at n.2_3. Shifts to end of A-tract; dup at n.6.
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}insA", &[2, 3]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[6]));
    }

    #[rstest]
    // Case 1: AC tandem at n.3-6 (ACAC). Insert AC at n.2_3 → dup at n.5_6.
    #[case("ACACACGTACGTACGTACGTAC", "n.{0}_{1}insAC", &[2u64, 3u64], "n.{0}_{1}dup", &[5u64, 6u64])]
    // Case 2: GCA tandem at n.3-11 (GCAGCAGCA). Insert GCA at n.2_3 → dup at n.9_11.
    #[case("ACGCAGCAGCATGACGTACG", "n.{0}_{1}insGCA", &[2u64, 3u64], "n.{0}_{1}dup", &[9u64, 11u64])]
    fn multi_base_ins_in_tandem_one_unit(
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
    // Issue #132 / #882: single-copy ins where alt is an OUT-OF-PHASE cyclic
    // rotation of the ref repeat unit. The candidate dup would decode to a
    // different sequence than the input insertion, so the #882 gate rejects it
    // and the variant stays a plain `ins` at the input position
    // (duplication.md:19). Mirrors the genomic module case shape.
    //
    // Case 1: GT tract at n.3-8. Insert TG (out-of-phase r=1 of GT) at n.2_3 → stays ins.
    #[case("ACGTGTGTACGTACGTACGT", "n.{0}_{1}insTG", &[2u64, 3u64], "n.{0}_{1}insTG", &[2u64, 3u64])]
    // Case 2: GCA tract at n.3-11. Insert CAG (out-of-phase r=1 of GCA) at n.2_3 → stays ins.
    #[case("TTGCAGCAGCATTGACGTAC", "n.{0}_{1}insCAG", &[2u64, 3u64], "n.{0}_{1}insCAG", &[2u64, 3u64])]
    fn cyclic_rotation_ins_one_unit(
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
    // Case 1: 2 A's added to AAA tract → A[5]. Tract: n.3..5.
    #[case("ACAAACGTACGTACGTACGT", "n.{0}_{1}insAA", &[2u64, 3u64], "n.{0}_{1}A[5]", &[3u64, 5u64])]
    // Case 2: 4 A's added to AAA tract → A[7]. Tract: n.3..5.
    #[case("ACAAACGTACGTACGTACGT", "n.{0}_{1}insAAAA", &[2u64, 3u64], "n.{0}_{1}A[7]", &[3u64, 5u64])]
    // Case 3: 2 GCA units added → GCA[5]. Tract: n.3..11.
    #[case("ACGCAGCAGCATGACGTACG", "n.{0}_{1}insGCAGCA", &[2u64, 3u64], "n.{0}_{1}GCA[5]", &[3u64, 11u64])]
    // Case 4: 2 AC units added to ACACAC tract → AC[5]. Tract: n.3..8.
    #[case("GCACACACGTACGTACGTAC", "n.{0}_{1}insACAC", &[2u64, 3u64], "n.{0}_{1}AC[5]", &[3u64, 8u64])]
    fn multi_base_ins_in_tandem_multiple_units(
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
    fn partial_rotation_no_full_match() {
        // Core: ACGTACGTACGTACGTACGTAC. Insert GA at n.2_3.
        //   ref[n.3]=G=alt[0] → shift, alt rotates to "AG"
        //   ref[n.4]=T, alt[1]=A → mismatch, stop.
        // Output: n.3_4insAG.
        let p = provider("ACGTACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}insGA", &[2, 3]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}_{1}insAG", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        // Boundary clamp via end-of-transcript: A-tract runs into the last
        // position of the transcript. Shuffling stops at the transcript end.
        //
        // Core: CAAAA (5 bases — entire transcript). A-tract n.2..5; the
        // transcript ENDS at n.5 (end-of-ref boundary). Insert A at n.2_3.
        // Should clamp to n.5 — dup at n.5.
        let p = provider("CAAAA");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}insA", &[2, 3]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[5]));
    }

    #[rstest]
    fn insertion_at_5prime_edge() {
        // Insertion at n.1_2. Core: AAAACGTACGTACGTACGTAC.
        // A-tract at n.1..4. Shifts to dup at n.4.
        let p = provider("AAAACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}insA", &[1, 2]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[4]));
    }

    #[rstest]
    fn insertion_at_3prime_edge() {
        // Insertion within an A-tract near transcript end.
        // Core: CGTACGTACGTACGAAAAC (19 bases). A-tract at n.15..18.
        // Insert A at n.16_17 (within tract); shifts to n.18, dup at n.18.
        let p = provider("CGTACGTACGTACGAAAAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}insA", &[16, 17]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[18]));
    }
}
mod noncoding_minus {
    use crate::common::synthetic::{normalize_to_string, SyntheticBuilder};
    use ferro_hgvs::reference::transcript::Strand;
    use rstest::rstest;

    /// Build a minus-strand noncoding provider whose transcript spans the entire
    /// core sequence.
    fn provider(core: &str) -> ferro_hgvs::MockProvider {
        SyntheticBuilder::noncoding(core, Strand::Minus).build()
    }

    /// Format an HGVS string with `{n}` placeholders. Noncoding positions are
    /// 1-based transcript-relative, no PAD_OFFSET — the caller passes raw
    /// 1-based positions.
    fn hgvs(template: &str, args: &[u64]) -> String {
        let mut s = template.to_string();
        for (i, p) in args.iter().enumerate() {
            s = s.replace(&format!("{{{}}}", i), &p.to_string());
        }
        s
    }

    #[rstest]
    fn no_shift_first_base_mismatch() {
        // Core: ACGTACGTACGTACGTACGTACGTACGTAC (30 bases). n.10=C, n.11=A.
        // Insert C between n.10 and n.11 — alt[0]='C' ≠ ref[n.11]='A' (no shift)
        // and alt[last]='C' ≠ ref[n.10]='C'... wait that's equal, would dup.
        // Pick a different position: n.4_5. n.4=T, n.5=A. Insert C.
        //   alt[0]='C' ≠ ref[n.5]='A' → no shift
        //   alt[last]='C' ≠ ref[n.4]='T' → no preceding-base dup
        let p = provider("ACGTACGTACGTACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}insC", &[4, 5]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}_{1}insC", &[4, 5]));
    }

    #[rstest]
    fn single_base_ins_in_homopolymer_one_unit() {
        // Core: ACAAAACGTACGTACGTAC (19 bases). A-tract at n.3..6.
        // Insert A at n.2_3. Shifts to end of A-tract; dup at n.6.
        let p = provider("ACAAAACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}insA", &[2, 3]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[6]));
    }

    #[rstest]
    // Case 1: AC tandem at n.3-6 (ACAC). Insert AC at n.2_3 → dup at n.5_6.
    #[case("ACACACGTACGTACGTACGTAC", "n.{0}_{1}insAC", &[2u64, 3u64], "n.{0}_{1}dup", &[5u64, 6u64])]
    // Case 2: GCA tandem at n.3-11 (GCAGCAGCA). Insert GCA at n.2_3 → dup at n.9_11.
    #[case("ACGCAGCAGCATGACGTACG", "n.{0}_{1}insGCA", &[2u64, 3u64], "n.{0}_{1}dup", &[9u64, 11u64])]
    fn multi_base_ins_in_tandem_one_unit(
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
    // Issue #132 / #882: single-copy ins where alt is an OUT-OF-PHASE cyclic
    // rotation of the ref repeat unit. The candidate dup would decode to a
    // different sequence than the input insertion, so the #882 gate rejects it
    // and the variant stays a plain `ins` at the input position
    // (duplication.md:19). Mirrors the genomic module case shape.
    //
    // Case 1: GT tract at n.3-8. Insert TG (out-of-phase r=1 of GT) at n.2_3 → stays ins.
    #[case("ACGTGTGTACGTACGTACGT", "n.{0}_{1}insTG", &[2u64, 3u64], "n.{0}_{1}insTG", &[2u64, 3u64])]
    // Case 2: GCA tract at n.3-11. Insert CAG (out-of-phase r=1 of GCA) at n.2_3 → stays ins.
    #[case("TTGCAGCAGCATTGACGTAC", "n.{0}_{1}insCAG", &[2u64, 3u64], "n.{0}_{1}insCAG", &[2u64, 3u64])]
    fn cyclic_rotation_ins_one_unit(
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
    // Case 1: 2 A's added to AAA tract → A[5]. Tract: n.3..5.
    #[case("ACAAACGTACGTACGTACGT", "n.{0}_{1}insAA", &[2u64, 3u64], "n.{0}_{1}A[5]", &[3u64, 5u64])]
    // Case 2: 4 A's added to AAA tract → A[7]. Tract: n.3..5.
    #[case("ACAAACGTACGTACGTACGT", "n.{0}_{1}insAAAA", &[2u64, 3u64], "n.{0}_{1}A[7]", &[3u64, 5u64])]
    // Case 3: 2 GCA units added → GCA[5]. Tract: n.3..11.
    #[case("ACGCAGCAGCATGACGTACG", "n.{0}_{1}insGCAGCA", &[2u64, 3u64], "n.{0}_{1}GCA[5]", &[3u64, 11u64])]
    // Case 4: 2 AC units added to ACACAC tract → AC[5]. Tract: n.3..8.
    #[case("GCACACACGTACGTACGTAC", "n.{0}_{1}insACAC", &[2u64, 3u64], "n.{0}_{1}AC[5]", &[3u64, 8u64])]
    fn multi_base_ins_in_tandem_multiple_units(
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
    fn partial_rotation_no_full_match() {
        // Core: ACGTACGTACGTACGTACGTAC. Insert GA at n.2_3.
        //   ref[n.3]=G=alt[0] → shift, alt rotates to "AG"
        //   ref[n.4]=T, alt[1]=A → mismatch, stop.
        // Output: n.3_4insAG.
        let p = provider("ACGTACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}insGA", &[2, 3]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}_{1}insAG", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        // Boundary clamp via end-of-transcript: A-tract runs into the last
        // position of the transcript. Shuffling stops at the transcript end.
        //
        // Core: CAAAA (5 bases — entire transcript). A-tract n.2..5; the
        // transcript ENDS at n.5 (end-of-ref boundary). Insert A at n.2_3.
        // Should clamp to n.5 — dup at n.5.
        let p = provider("CAAAA");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}insA", &[2, 3]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[5]));
    }

    #[rstest]
    fn insertion_at_5prime_edge() {
        // Insertion at n.1_2. Core: AAAACGTACGTACGTACGTAC.
        // A-tract at n.1..4. Shifts to dup at n.4.
        let p = provider("AAAACGTACGTACGTACGTAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}insA", &[1, 2]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[4]));
    }

    #[rstest]
    fn insertion_at_3prime_edge() {
        // Insertion within an A-tract near transcript end.
        // Core: CGTACGTACGTACGAAAAC (19 bases). A-tract at n.15..18.
        // Insert A at n.16_17 (within tract); shifts to n.18, dup at n.18.
        let p = provider("CGTACGTACGTACGAAAAC");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:n.{0}_{1}insA", &[16, 17]));
        assert_eq!(result, hgvs("NR_TEST.1:n.{0}dup", &[18]));
    }
}
mod rna_plus {
    use crate::common::synthetic::{normalize_to_string, SyntheticBuilder};
    use ferro_hgvs::reference::transcript::Strand;
    use rstest::rstest;

    /// Build a plus-strand RNA provider whose transcript spans the entire
    /// core sequence.
    fn provider(core: &str) -> ferro_hgvs::MockProvider {
        SyntheticBuilder::rna(core, Strand::Plus).build()
    }

    /// Format an HGVS string with `{n}` placeholders. RNA positions are
    /// 1-based transcript-relative, no PAD_OFFSET — the caller passes raw
    /// 1-based positions.
    fn hgvs(template: &str, args: &[u64]) -> String {
        let mut s = template.to_string();
        for (i, p) in args.iter().enumerate() {
            s = s.replace(&format!("{{{}}}", i), &p.to_string());
        }
        s
    }

    #[rstest]
    fn no_shift_first_base_mismatch() {
        // Core: acguacguacguacguacguacguacguac (30 bases, RNA). r.4=u, r.5=a.
        // Insert c between r.4 and r.5 — alt[0]='c' ≠ ref[r.5]='a' (no shift)
        // and alt[last]='c' ≠ ref[r.4]='u' → no preceding-base dup.
        let p = provider("acguacguacguacguacguacguacguac");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}_{1}insc", &[4, 5]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}_{1}insc", &[4, 5]));
    }

    #[rstest]
    fn single_base_ins_in_homopolymer_one_unit() {
        // Core: acaaaacguacguacguac (19 bases, RNA). A-tract at r.3..6.
        // Insert a at r.2_3. Shifts to end of A-tract; dup at r.6.
        let p = provider("acaaaacguacguacguac");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}_{1}insa", &[2, 3]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}dup", &[6]));
    }

    #[rstest]
    // Case 1: ac tandem at r.3-6 (acac). Insert ac at r.2_3 → dup at r.5_6.
    #[case("acacacguacguacguacguac", "r.{0}_{1}insac", &[2u64, 3u64], "r.{0}_{1}dup", &[5u64, 6u64])]
    // Case 2: gca tandem at r.3-11 (gcagcagca). Insert gca at r.2_3 → dup at r.9_11.
    #[case("acgcagcagcaugacguacg", "r.{0}_{1}insgca", &[2u64, 3u64], "r.{0}_{1}dup", &[9u64, 11u64])]
    fn multi_base_ins_in_tandem_one_unit(
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
    // Issue #132 / #882: single-copy ins where alt is an OUT-OF-PHASE cyclic
    // rotation of the ref repeat unit. The candidate dup would decode to a
    // different sequence than the input insertion, so the #882 gate rejects it
    // and the variant stays a plain `ins` at the input position
    // (duplication.md:19). RNA equivalent of the genomic case shape.
    //
    // Note: avoid units containing 'u' — ferro stores RNA refs as DNA
    // (U→T). A `ug` rotation alt against a `gu` ref tract becomes byte
    // 'U' vs 'T' and the rotation iterator won't match. The case below
    // uses only A/C/G, which are unambiguous between DNA and RNA byte
    // representations.
    //
    // Case 1: gca tract at r.3-11. Insert cag (out-of-phase r=1 of gca) at r.2_3 → stays ins.
    #[case("uugcagcagcauugacguac", "r.{0}_{1}inscag", &[2u64, 3u64], "r.{0}_{1}inscag", &[2u64, 3u64])]
    fn cyclic_rotation_ins_one_unit(
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
    // Case 1: 2 a's added to aaa tract → a[5]. Tract: r.3..5.
    #[case("acaaacguacguacguacgu", "r.{0}_{1}insaa", &[2u64, 3u64], "r.{0}_{1}a[5]", &[3u64, 5u64])]
    // Case 2: 4 a's added to aaa tract → a[7]. Tract: r.3..5.
    #[case("acaaacguacguacguacgu", "r.{0}_{1}insaaaa", &[2u64, 3u64], "r.{0}_{1}a[7]", &[3u64, 5u64])]
    // Case 3: 2 gca units added → gca[5]. Tract: r.3..11.
    #[case("acgcagcagcaugacguacg", "r.{0}_{1}insgcagca", &[2u64, 3u64], "r.{0}_{1}gca[5]", &[3u64, 11u64])]
    // Case 4: 2 ac units added to acacac tract → ac[5]. Tract: r.3..8.
    #[case("gcacacacguacguacguac", "r.{0}_{1}insacac", &[2u64, 3u64], "r.{0}_{1}ac[5]", &[3u64, 8u64])]
    fn multi_base_ins_in_tandem_multiple_units(
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
    fn partial_rotation_no_full_match() {
        // Core: acguacguacguacguacguac (RNA). Insert ga at r.2_3.
        //   ref[r.3]=g=alt[0] → shift, alt rotates to "ag"
        //   ref[r.4]=u, alt[1]=a → mismatch, stop.
        // Output: r.3_4insag.
        let p = provider("acguacguacguacguacguac");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}_{1}insga", &[2, 3]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}_{1}insag", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        // Boundary clamp via end-of-transcript: a-tract runs into the last
        // position of the transcript. Shuffling stops at the transcript end.
        //
        // Core: caaaa (5 bases, RNA). a-tract r.2..5; the
        // transcript ENDS at r.5 (end-of-ref boundary). Insert a at r.2_3.
        // Should clamp to r.5 — dup at r.5.
        let p = provider("caaaa");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}_{1}insa", &[2, 3]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}dup", &[5]));
    }

    #[rstest]
    fn insertion_at_5prime_edge() {
        // Insertion at r.1_2. Core: aaaacguacguacguacguac (RNA).
        // a-tract at r.1..4. Shifts to dup at r.4.
        let p = provider("aaaacguacguacguacguac");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}_{1}insa", &[1, 2]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}dup", &[4]));
    }

    #[rstest]
    fn insertion_at_3prime_edge() {
        // Insertion within an a-tract near transcript end.
        // Core: cguacguacguacgaaaac (19 bases, RNA). a-tract at r.15..18.
        // Insert a at r.16_17 (within tract); shifts to r.18, dup at r.18.
        let p = provider("cguacguacguacgaaaac");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}_{1}insa", &[16, 17]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}dup", &[18]));
    }
}
mod rna_minus {
    use crate::common::synthetic::{normalize_to_string, SyntheticBuilder};
    use ferro_hgvs::reference::transcript::Strand;
    use rstest::rstest;

    /// Build a minus-strand RNA provider whose transcript spans the entire
    /// core sequence.
    fn provider(core: &str) -> ferro_hgvs::MockProvider {
        SyntheticBuilder::rna(core, Strand::Minus).build()
    }

    /// Format an HGVS string with `{n}` placeholders. RNA positions are
    /// 1-based transcript-relative, no PAD_OFFSET — the caller passes raw
    /// 1-based positions.
    fn hgvs(template: &str, args: &[u64]) -> String {
        let mut s = template.to_string();
        for (i, p) in args.iter().enumerate() {
            s = s.replace(&format!("{{{}}}", i), &p.to_string());
        }
        s
    }

    #[rstest]
    fn no_shift_first_base_mismatch() {
        // Core: acguacguacguacguacguacguacguac (30 bases, RNA). r.4=u, r.5=a.
        // Insert c between r.4 and r.5 — alt[0]='c' ≠ ref[r.5]='a' (no shift)
        // and alt[last]='c' ≠ ref[r.4]='u' → no preceding-base dup.
        let p = provider("acguacguacguacguacguacguacguac");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}_{1}insc", &[4, 5]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}_{1}insc", &[4, 5]));
    }

    #[rstest]
    fn single_base_ins_in_homopolymer_one_unit() {
        // Core: acaaaacguacguacguac (19 bases, RNA). a-tract at r.3..6.
        // Insert a at r.2_3. Shifts to end of a-tract; dup at r.6.
        let p = provider("acaaaacguacguacguac");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}_{1}insa", &[2, 3]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}dup", &[6]));
    }

    #[rstest]
    // Case 1: ac tandem at r.3-6 (acac). Insert ac at r.2_3 → dup at r.5_6.
    #[case("acacacguacguacguacguac", "r.{0}_{1}insac", &[2u64, 3u64], "r.{0}_{1}dup", &[5u64, 6u64])]
    // Case 2: gca tandem at r.3-11 (gcagcagca). Insert gca at r.2_3 → dup at r.9_11.
    #[case("acgcagcagcaugacguacg", "r.{0}_{1}insgca", &[2u64, 3u64], "r.{0}_{1}dup", &[9u64, 11u64])]
    fn multi_base_ins_in_tandem_one_unit(
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
    // Issue #132 / #882: single-copy ins where alt is an OUT-OF-PHASE cyclic
    // rotation of the ref repeat unit. The candidate dup would decode to a
    // different sequence than the input insertion, so the #882 gate rejects it
    // and the variant stays a plain `ins` at the input position
    // (duplication.md:19). RNA equivalent of the genomic case shape.
    //
    // Note: avoid units containing 'u' — ferro stores RNA refs as DNA
    // (U→T). A `ug` rotation alt against a `gu` ref tract becomes byte
    // 'U' vs 'T' and the rotation iterator won't match. The case below
    // uses only A/C/G, which are unambiguous between DNA and RNA byte
    // representations.
    //
    // Case 1: gca tract at r.3-11. Insert cag (out-of-phase r=1 of gca) at r.2_3 → stays ins.
    #[case("uugcagcagcauugacguac", "r.{0}_{1}inscag", &[2u64, 3u64], "r.{0}_{1}inscag", &[2u64, 3u64])]
    fn cyclic_rotation_ins_one_unit(
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
    // Case 1: 2 a's added to aaa tract → a[5]. Tract: r.3..5.
    #[case("acaaacguacguacguacgu", "r.{0}_{1}insaa", &[2u64, 3u64], "r.{0}_{1}a[5]", &[3u64, 5u64])]
    // Case 2: 4 a's added to aaa tract → a[7]. Tract: r.3..5.
    #[case("acaaacguacguacguacgu", "r.{0}_{1}insaaaa", &[2u64, 3u64], "r.{0}_{1}a[7]", &[3u64, 5u64])]
    // Case 3: 2 gca units added → gca[5]. Tract: r.3..11.
    #[case("acgcagcagcaugacguacg", "r.{0}_{1}insgcagca", &[2u64, 3u64], "r.{0}_{1}gca[5]", &[3u64, 11u64])]
    // Case 4: 2 ac units added to acacac tract → ac[5]. Tract: r.3..8.
    #[case("gcacacacguacguacguac", "r.{0}_{1}insacac", &[2u64, 3u64], "r.{0}_{1}ac[5]", &[3u64, 8u64])]
    fn multi_base_ins_in_tandem_multiple_units(
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
    fn partial_rotation_no_full_match() {
        // Core: acguacguacguacguacguac (RNA). Insert ga at r.2_3.
        //   ref[r.3]=g=alt[0] → shift, alt rotates to "ag"
        //   ref[r.4]=u, alt[1]=a → mismatch, stop.
        // Output: r.3_4insag.
        let p = provider("acguacguacguacguacguac");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}_{1}insga", &[2, 3]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}_{1}insag", &[3, 4]));
    }

    #[rstest]
    fn shift_clamped_by_boundary() {
        // Boundary clamp via end-of-transcript: a-tract runs into the last
        // position of the transcript. Shuffling stops at the transcript end.
        //
        // Core: caaaa (5 bases, RNA). a-tract r.2..5; the
        // transcript ENDS at r.5 (end-of-ref boundary). Insert a at r.2_3.
        // Should clamp to r.5 — dup at r.5.
        let p = provider("caaaa");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}_{1}insa", &[2, 3]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}dup", &[5]));
    }

    #[rstest]
    fn insertion_at_5prime_edge() {
        // Insertion at r.1_2. Core: aaaacguacguacguacguac (RNA).
        // a-tract at r.1..4. Shifts to dup at r.4.
        let p = provider("aaaacguacguacguacguac");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}_{1}insa", &[1, 2]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}dup", &[4]));
    }

    #[rstest]
    fn insertion_at_3prime_edge() {
        // Insertion within an a-tract near transcript end.
        // Core: cguacguacguacgaaaac (19 bases, RNA). a-tract at r.15..18.
        // Insert a at r.16_17 (within tract); shifts to r.18, dup at r.18.
        let p = provider("cguacguacguacgaaaac");
        let result = normalize_to_string(p, &hgvs("NR_TEST.1:r.{0}_{1}insa", &[16, 17]));
        assert_eq!(result, hgvs("NR_TEST.1:r.{0}dup", &[18]));
    }
}
