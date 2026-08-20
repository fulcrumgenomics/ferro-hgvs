//! #2161 — the lone-member derivation fast-path.
//!
//! `rederive(recommended_form = true)` short-circuits a lone `g.`/`m.`
//! del/ins/dup straight to its final normalize step, skipping the confluence
//! loop and `from_sequences`. Such a member has no sibling to re-partition
//! against and cannot denote an identity, so its sequence-first derivation lands
//! on the exact canonical form the final normalize produces positionally — this
//! file pins that the fast-path is **byte-identical** to the full derivation and
//! that its boundary (del/ins/dup in; `delins`/`inv`/alleles out) is load-bearing.
//!
//! # How the two paths are reconstructed
//!
//! `rederive(.., false)` returns the *derived* form and is **not** gated on
//! `recommended_form`, so it never takes the fast-path — it is the full
//! sequence-first derivation. Re-applying the final normalize to it
//! (`normalize_gated_repartition`, exactly what `rederive(.., true)` runs at its
//! settled exit) reconstructs the pre-fast-path `rederive(true)` output. The
//! shipped fast-path is `rederive(.., true)` itself. The guard asserts the two
//! agree on every input.

use ferro_hgvs::{parse_hgvs, FerroError, FromSequencesOptions, JsonProvider, Normalizer};
use std::io::Write;

/// A genome-capable provider over four structured contigs (the same shapes the
/// inversion-geometry guard uses): a general mix, an inverted repeat, a tandem
/// tract, and a GC palindrome — so the corpus exercises homopolymer shifting,
/// repeat tracts, and reverse-complement coincidence, not one flat sequence.
fn provider() -> JsonProvider {
    let doc = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [],
        "genomic_sequences": {
            "NC_TEST.1": "GCTAGCATGCATGCGTACAGTCGATCGATCAAAAAATGCAGTCAGTGGATCCGATTACGATCAGCT",
            "NC_TEST.2": "AACCGGTTAATCGATCGATTGCACGTACGTGCAATCGATCGATTAACCGGTTAACCGGTTAACCGG",
            "NC_TEST.3": "GGCCAATTGGCCATATATATATATATGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATT",
            "NC_TEST.4": "ATATATATATGGGCGCGCCCATATATATATGGGCGCGCCCATATATATATGGGCGCGCCCATATAT",
        },
    });
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(doc.to_string().as_bytes()).unwrap();
    JsonProvider::from_json(f.path()).unwrap()
}

/// `rederive(recommended_form = true)` — the shipped path, which takes the
/// fast-path for a lone del/ins/dup.
fn rederive_shipped(nz: &Normalizer<JsonProvider>, desc: &str) -> Result<String, FerroError> {
    let v = parse_hgvs(desc)?;
    let opts = FromSequencesOptions::default();
    Ok(nz.rederive(&v, &opts, true)?.to_string())
}

/// The full sequence-first path, reconstructed without the fast-path:
/// derive (`recommended_form = false`, never fast-pathed), then re-apply the
/// same final normalize the shipped path runs at its settled exit.
fn rederive_full(nz: &Normalizer<JsonProvider>, desc: &str) -> Result<String, FerroError> {
    let v = parse_hgvs(desc)?;
    let opts = FromSequencesOptions::default();
    let derived = nz.rederive(&v, &opts, false)?;
    Ok(nz.normalize_gated_repartition(&derived)?.to_string())
}

/// The fast-path is byte-identical to the full derivation across a diverse
/// corpus — every lone del/ins/dup on four contigs, plus controls (lone
/// `delins`/`inv` and multi-member alleles) the fast-path must **not** take.
#[test]
fn the_lone_member_fast_path_is_byte_identical() {
    let nz = Normalizer::new(provider());
    let opts = FromSequencesOptions::default();
    let contigs = ["NC_TEST.1", "NC_TEST.2", "NC_TEST.3", "NC_TEST.4"];

    let mut compared = 0usize;
    let mut fired = 0usize;
    for ctg in contigs {
        for a in 5u64..=55 {
            // Fast-path members (must fire) and controls (must not).
            let cases: [(String, bool); 8] = [
                (format!("{ctg}:g.{a}del"), true),
                (format!("{ctg}:g.{a}_{}del", a + 2), true),
                (format!("{ctg}:g.{a}dup"), true),
                (format!("{ctg}:g.{a}_{}dup", a + 3), true),
                (format!("{ctg}:g.{a}_{}insAT", a + 1), true),
                (format!("{ctg}:g.{a}_{}inv", a + 3), false),
                (format!("{ctg}:g.{a}_{}delinsGC", a + 3), false),
                (format!("{ctg}:g.[{a}del;{}del]", a + 5), false),
            ];
            for (desc, expect_fast) in &cases {
                let Ok(v) = parse_hgvs(desc) else { continue };
                assert_eq!(
                    nz.lone_member_fast_path_applies(&v),
                    *expect_fast,
                    "fast-path boundary wrong for {desc}"
                );
                if *expect_fast {
                    fired += 1;
                }
                // Compare the two paths, and never let an *asymmetric* outcome
                // slip through: if one path yields an output and the other
                // declines, the fast-path and the full derivation disagree on
                // whether this input even has a result — exactly the divergence
                // this guard exists to catch. Only a *symmetric* decline is
                // legitimate (some inputs simply do not derive, which is not a
                // fast-path question).
                match (rederive_shipped(&nz, desc), rederive_full(&nz, desc)) {
                    (Ok(shipped), Ok(full)) => {
                        assert_eq!(shipped, full, "fast-path changed the output of {desc}");
                        compared += 1;
                    }
                    (Err(_), Err(_)) => {}
                    (shipped, full) => panic!(
                        "asymmetric derivation for {desc}: shipped={shipped:?} full={full:?}"
                    ),
                }
                let _ = &opts;
            }
        }
    }

    // Non-vacuity: the fast-path must actually fire, or the equality above is
    // proving nothing. Four contigs x 51 positions x 5 fast-path members is well
    // over a thousand; a floor far below that still catches a predicate that
    // stopped firing.
    assert!(
        fired > 500,
        "fast-path fired on only {fired} inputs — guard is vacuous"
    );
    assert!(
        compared > 500,
        "compared only {compared} inputs — guard is vacuous"
    );
}

/// The identity boundary is load-bearing: a palindromic `inv` and a self-`delins`
/// denote no change, and the derivation spells that as whole-reference `g.=`
/// while the positional normalize keeps `span=`. So these two classes genuinely
/// diverge — which is why the fast-path excludes `inv`/`delins`. If the predicate
/// were widened to include them, `the_lone_member_fast_path_is_byte_identical`
/// would go red; this test pins the divergence that makes the exclusion necessary.
#[test]
fn the_identity_boundary_is_load_bearing() {
    let nz = Normalizer::new(provider());

    // NC_TEST.1[6..9] = CATG, whose reverse complement is itself (a palindrome),
    // and a delins replacing it with the same bases — both denote no change.
    for desc in ["NC_TEST.1:g.6_9inv", "NC_TEST.1:g.6_9delinsCATG"] {
        let v = parse_hgvs(desc).unwrap();
        assert!(
            !nz.lone_member_fast_path_applies(&v),
            "{desc} must be excluded from the fast-path"
        );
        let shipped = rederive_shipped(&nz, desc).expect("rederive");
        let full = rederive_full(&nz, desc).expect("full");
        // The two paths agree (both go through the derivation, since the
        // fast-path is excluded) …
        assert_eq!(shipped, full, "{desc}");
        // … and the derivation collapses the no-change edit to whole-reference
        // identity, which a positional normalize of the *input* does not — the
        // exact divergence the exclusion avoids.
        let positional = nz.normalize_gated_repartition(&v).unwrap().to_string();
        assert_ne!(
            shipped, positional,
            "{desc}: derivation and positional normalize must differ here, or the \
             fast-path exclusion would be unnecessary"
        );
    }
}
