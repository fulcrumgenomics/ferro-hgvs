//! Boundary detection tests
//!
//! Tests for normalization boundary detection including:
//! - Exon/intron boundaries
//! - UTR/CDS boundaries
//! - Transcript start/end boundaries

use ferro_hgvs::normalize::boundary::Boundaries;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};

fn make_multi_exon_transcript() -> Transcript {
    // Transcript with 3 exons:
    // Exon 1: 1-100 (5' UTR: 1-20, CDS: 21-100)
    // Exon 2: 201-300 (CDS: 201-300)
    // Exon 3: 401-500 (CDS: 401-450, 3' UTR: 451-500)
    Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        "A".repeat(500),
        Some(21),
        Some(450),
        vec![
            Exon::new(1, 1, 100),
            Exon::new(2, 201, 300),
            Exon::new(3, 401, 500),
        ],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    )
}

fn make_single_exon_transcript() -> Transcript {
    // Simple single-exon transcript
    // 5' UTR: 1-10, CDS: 11-90, 3' UTR: 91-100
    Transcript::new(
        "NM_SIMPLE.1".to_string(),
        Some("SIMPLE".to_string()),
        Strand::Plus,
        "A".repeat(100),
        Some(11),
        Some(90),
        vec![Exon::new(1, 1, 100)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    )
}

mod boundary_struct {
    use super::*;

    #[test]
    fn test_boundaries_contains() {
        let bounds = Boundaries::new(10, 50);

        // Inside bounds
        assert!(bounds.contains(10)); // Left edge
        assert!(bounds.contains(30)); // Middle
        assert!(bounds.contains(50)); // Right edge

        // Outside bounds
        assert!(!bounds.contains(9));
        assert!(!bounds.contains(51));
    }

    #[test]
    fn test_boundaries_width() {
        let bounds = Boundaries::new(10, 50);

        // Test width calculation
        assert_eq!(bounds.right - bounds.left, 40);
    }
}

mod exon_boundaries {
    use super::*;

    #[test]
    fn test_exon_contains_position() {
        let exon = Exon::new(1, 100, 200);

        assert!(exon.contains(100)); // Start
        assert!(exon.contains(150)); // Middle
        assert!(exon.contains(200)); // End
        assert!(!exon.contains(99));
        assert!(!exon.contains(201));
    }

    #[test]
    fn test_find_exon_at_position() {
        let tx = make_multi_exon_transcript();

        // Position in exon 1
        let exon1 = tx.exon_at(50);
        assert!(exon1.is_some());
        assert_eq!(exon1.unwrap().number, 1);

        // Position in exon 2
        let exon2 = tx.exon_at(250);
        assert!(exon2.is_some());
        assert_eq!(exon2.unwrap().number, 2);

        // Position in exon 3
        let exon3 = tx.exon_at(450);
        assert!(exon3.is_some());
        assert_eq!(exon3.unwrap().number, 3);

        // Position in intron (between exons)
        let intron = tx.exon_at(150);
        assert!(intron.is_none());
    }
}

mod transcript_regions {
    use super::*;

    #[test]
    fn test_is_in_5utr() {
        let tx = make_single_exon_transcript();

        // 5' UTR is positions 1-10
        assert!(is_in_5utr(&tx, 1));
        assert!(is_in_5utr(&tx, 10));
        assert!(!is_in_5utr(&tx, 11)); // CDS start
        assert!(!is_in_5utr(&tx, 50)); // Middle of CDS
    }

    #[test]
    fn test_is_in_cds() {
        let tx = make_single_exon_transcript();

        // CDS is positions 11-90
        assert!(!is_in_cds(&tx, 10)); // 5' UTR
        assert!(is_in_cds(&tx, 11)); // CDS start
        assert!(is_in_cds(&tx, 50)); // Middle
        assert!(is_in_cds(&tx, 90)); // CDS end
        assert!(!is_in_cds(&tx, 91)); // 3' UTR
    }

    #[test]
    fn test_is_in_3utr() {
        let tx = make_single_exon_transcript();

        // 3' UTR is positions 91-100
        assert!(!is_in_3utr(&tx, 90)); // CDS end
        assert!(is_in_3utr(&tx, 91));
        assert!(is_in_3utr(&tx, 100));
    }

    fn is_in_5utr(tx: &Transcript, pos: u64) -> bool {
        if let Some(cds_start) = tx.cds_start {
            pos < cds_start
        } else {
            false
        }
    }

    fn is_in_cds(tx: &Transcript, pos: u64) -> bool {
        match (tx.cds_start, tx.cds_end) {
            (Some(start), Some(end)) => pos >= start && pos <= end,
            _ => false,
        }
    }

    fn is_in_3utr(tx: &Transcript, pos: u64) -> bool {
        if let Some(cds_end) = tx.cds_end {
            pos > cds_end
        } else {
            false
        }
    }
}

mod boundary_edge_cases {
    use super::*;

    #[test]
    fn test_variant_at_exon_boundary() {
        let tx = make_multi_exon_transcript();

        // Position at exon 1 end (100)
        let exon = tx.exon_at(100);
        assert!(exon.is_some());
        assert_eq!(exon.unwrap().end, 100);

        // Position at exon 2 start (201)
        let exon = tx.exon_at(201);
        assert!(exon.is_some());
        assert_eq!(exon.unwrap().start, 201);
    }

    #[test]
    fn test_variant_at_cds_boundary() {
        let tx = make_single_exon_transcript();

        // CDS start boundary (position 11)
        assert!(!is_in_5utr(&tx, 11));
        assert!(is_in_cds(&tx, 11));

        // CDS end boundary (position 90)
        assert!(is_in_cds(&tx, 90));
        assert!(!is_in_3utr(&tx, 90));
    }

    #[test]
    fn test_non_coding_transcript() {
        // Transcript without CDS (e.g., lncRNA)
        let tx = Transcript::new(
            "NR_TEST.1".to_string(),
            Some("NCRNA".to_string()),
            Strand::Plus,
            "A".repeat(100),
            None,
            None,
            vec![Exon::new(1, 1, 100)],
            None,
            None,
            None,
            Default::default(),
            ManeStatus::default(),
            None,
            None,
        );

        assert!(!tx.is_coding());
        assert!(tx.exon_at(50).is_some());
    }

    fn is_in_5utr(tx: &Transcript, pos: u64) -> bool {
        if let Some(cds_start) = tx.cds_start {
            pos < cds_start
        } else {
            false
        }
    }

    fn is_in_cds(tx: &Transcript, pos: u64) -> bool {
        match (tx.cds_start, tx.cds_end) {
            (Some(start), Some(end)) => pos >= start && pos <= end,
            _ => false,
        }
    }

    fn is_in_3utr(tx: &Transcript, pos: u64) -> bool {
        if let Some(cds_end) = tx.cds_end {
            pos > cds_end
        } else {
            false
        }
    }
}

mod normalization_boundaries {
    use super::*;
    use ferro_hgvs::normalize::boundary::get_cds_boundaries;
    use ferro_hgvs::normalize::NormalizeConfig;

    #[test]
    fn test_cds_boundaries_single_exon() {
        let tx = make_single_exon_transcript();
        let config = NormalizeConfig::default();

        // Position in middle of CDS
        let bounds = get_cds_boundaries(&tx, 50, &config);
        assert!(bounds.is_ok());
        let bounds = bounds.unwrap();

        // Boundaries should be within the exon (may include UTRs depending on config)
        assert!(bounds.left >= 1);
        assert!(bounds.right <= 100);
    }

    #[test]
    fn test_cds_boundaries_cross_allowed() {
        let tx = make_single_exon_transcript();
        let config = NormalizeConfig::default().allow_crossing_boundaries();

        let bounds = get_cds_boundaries(&tx, 50, &config);
        assert!(bounds.is_ok());
        let bounds = bounds.unwrap();

        // With cross_boundaries=true, should allow full transcript range
        assert_eq!(bounds.left, 1);
        assert_eq!(bounds.right, tx.sequence.len() as u64);
    }
}
