//! Golden tests for the perf-table renderer.
#![cfg(feature = "dev")]

use ferro_hgvs::perf_table::PerfResults;

fn fixture() -> PerfResults {
    PerfResults::load("tests/fixtures/perf/sample_perf_results.json").expect("load fixture")
}

#[test]
fn fixture_loads_and_has_both_ops() {
    let r = fixture();
    assert!(r.operations.contains_key("parse"));
    assert!(r.operations.contains_key("normalize"));
}

#[test]
fn parse_cross_tool_renders_expected() {
    let r = fixture();
    let out = r.render_cross_tool("parse").unwrap();
    let expected = "\
| Tool | Throughput @ 1 worker | Throughput @ 8 workers | ferro speedup @ 8w |
|------|----------------------:|-----------------------:|-------------------:|
| ferro | 1.2M/s | 8.0M/s | — |
| mutalyzer | 21/s | single-threaded | 53,000× |
| biocommons | 20/s | single-threaded | 57,000× |
| hgvs-rs | 2/s | single-threaded | 570,000× |
";
    assert_eq!(out, expected);
}

#[test]
fn normalize_marks_not_run() {
    let r = fixture();
    let out = r.render_cross_tool("normalize").unwrap();
    // biocommons + hgvs-rs are not_run -> em-dash in both worker columns.
    assert!(out.contains("| biocommons | — | — | — |"));
    assert!(out.contains("| hgvs-rs | — | — | — |"));
    assert!(out.contains("| mutalyzer | 19/s |"));
}

#[test]
fn ferro_scaling_renders_expected() {
    let r = fixture();
    let out = r.render_ferro_scaling("parse").unwrap();
    assert_eq!(
        out,
        "\
| Threads | 1 | 2 | 4 | 8 |
|---------|--:|--:|--:|--:|
| ferro parse | 1.2M/s | 2.3M/s | 4.4M/s | 8.0M/s |
"
    );
}
