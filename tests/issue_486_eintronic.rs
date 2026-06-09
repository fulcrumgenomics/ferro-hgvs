//! #486 EINTRONIC — reject intronic offsets on a bare transcript reference.

use ferro_hgvs::error::FerroError;
use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer};

const PAD: &str = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
const PAD_OFFSET: u64 = 256;

fn padded(core: &str) -> String {
    format!("{}{}{}", PAD, core, PAD)
}

/// Plus-strand 2-exon coding transcript with a 10-bp intron seeded
/// `ATATATATGG` at intronic positions c.30+1..=c.30+10, with full genomic
/// sequence so indels resolve (mirrors issue_214's provider_with_intronic_at4).
fn provider_intronic_nm() -> MockProvider {
    let exon1 = "ATGCATGCATGCATGCATGCATGCATGCAT"; // 30 bp
    let intron = "ATATATATGG"; // 10 bp
    let exon2 = "AAACAACATGGAAAAAAAAAAAAAAAAAAA"; // 30 bp
    let core = format!("{}{}{}", exon1, intron, exon2);
    let tx_seq = format!("{}{}", exon1, exon2);
    let p_off = PAD_OFFSET;
    let transcript = Transcript::new(
        "NM_INTRON.1".to_string(),
        Some("INTRONGENE".to_string()),
        Strand::Plus,
        tx_seq,
        Some(1),
        Some(60),
        vec![
            Exon::with_genomic(1, 1, 30, p_off, p_off + 29),
            Exon::with_genomic(2, 31, 60, p_off + 40, p_off + 69),
        ],
        Some("chr_intron".to_string()),
        Some(p_off),
        Some(p_off + 69),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );
    let mut p = MockProvider::new();
    p.add_genomic_sequence("chr_intron", padded(&core));
    p.add_transcript(transcript);
    p
}

fn normalize(
    provider: MockProvider,
    cfg: NormalizeConfig,
    input: &str,
) -> Result<String, FerroError> {
    let normalizer = Normalizer::with_config(provider, cfg);
    let variant = parse_hgvs(input).expect("parse failed");
    normalizer.normalize(&variant).map(|v| format!("{}", v))
}

#[test]
fn strict_rejects_bare_nm_intronic_substitution() {
    let err = normalize(
        provider_intronic_nm(),
        NormalizeConfig::strict(),
        "NM_INTRON.1:c.30+1A>G",
    )
    .expect_err("strict mode must reject a bare-NM intronic substitution as EINTRONIC");
    assert!(
        matches!(err, FerroError::IntronicVariant { .. }),
        "expected IntronicVariant (EINTRONIC), got {:?}",
        err
    );
}
