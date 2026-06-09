//! Measurement (ignored): how much does routing the default through the fast
//! path save on a representative corpus? Times the forced-generic
//! `parse_variant` against the default `parse_hgvs` (now fast-path) over
//! `FERRO_DIFF_CORPUS_TXT` (newline-delimited, e.g. the 500k ClinVar inputs),
//! min-of-reps.
//!
//!   cargo test --release --test fast_path_flip_bench --features dev -- \
//!     --nocapture --ignored

use ferro_hgvs::hgvs::parser::variant::parse_variant;
use ferro_hgvs::parse_hgvs;
use std::hint::black_box;
use std::time::Instant;

#[test]
#[ignore = "measurement, run explicitly with --release --ignored --nocapture; needs FERRO_DIFF_CORPUS_TXT"]
fn measure_flip_net_effect() {
    let Ok(path) = std::env::var("FERRO_DIFF_CORPUS_TXT") else {
        println!("\nSKIP: set FERRO_DIFF_CORPUS_TXT to a newline-delimited HGVS file\n");
        return;
    };
    let inputs: Vec<String> = std::fs::read_to_string(&path)
        .expect("read corpus")
        .lines()
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty())
        .collect();
    if inputs.is_empty() {
        println!("\nSKIP: corpus file {path} has no non-empty lines\n");
        return;
    }
    println!("\ncorpus: {} inputs", inputs.len());

    // Count how many the fast path actually claims (the rest fall back), so we
    // can interpret the net number.
    let mut fast_claimed = 0usize;
    for s in &inputs {
        if matches!(
            ferro_hgvs::hgvs::parser::fast_path::try_fast_path(s),
            ferro_hgvs::hgvs::parser::fast_path::FastPathResult::Success(_)
        ) {
            fast_claimed += 1;
        }
    }
    println!(
        "fast-path claims {fast_claimed}/{} ({:.1}%); rest fall back to generic",
        inputs.len(),
        100.0 * fast_claimed as f64 / inputs.len() as f64
    );

    let reps = 5;
    let mut best_generic = f64::MAX;
    let mut best_fast = f64::MAX;
    // Interleave reps so transient system load hits both equally.
    for _ in 0..reps {
        let t = Instant::now();
        let mut ok = 0usize;
        for s in &inputs {
            ok += black_box(parse_variant(s)).is_ok() as usize;
        }
        black_box(ok);
        best_generic = best_generic.min(t.elapsed().as_secs_f64());

        let t = Instant::now();
        let mut ok = 0usize;
        for s in &inputs {
            ok += black_box(parse_hgvs(s)).is_ok() as usize;
        }
        black_box(ok);
        best_fast = best_fast.min(t.elapsed().as_secs_f64());
    }

    let n = inputs.len() as f64;
    println!("\n=== flip net effect (min of {reps} reps) ===");
    println!(
        "  forced-generic parse_variant: {:.3} s   ({:.2} M/s)",
        best_generic,
        n / best_generic / 1e6
    );
    println!(
        "  default parse_hgvs (fast):    {:.3} s   ({:.2} M/s)",
        best_fast,
        n / best_fast / 1e6
    );
    println!(
        "  speedup (generic/fast):  {:.3}x   ({:+.1}%)\n",
        best_generic / best_fast,
        100.0 * (best_generic / best_fast - 1.0)
    );
}
