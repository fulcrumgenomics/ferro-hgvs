# Pickup Notes: In-Memory hgvs-rs Provider

## Status

Tasks 1-5 complete and committed. Task 6 (smoke test) blocked on slow UTA view queries.

## What's done

- `InMemoryProvider` struct with full `Provider` trait implementation (`src/benchmark/inmemory_provider.rs`)
- Pre-created normalizer pool to avoid per-chunk Provider init (`src/benchmark/hgvs_rs.rs`)
- `HgvsRsNormalizer::with_provider()` constructor accepting `Arc<dyn Provider>`
- `--in-memory` CLI flag wired through `ferro-benchmark normalize hgvs-rs`
- All compiles clean, clippy clean

## What's next

### Immediate: Serialize-to-disk cache (Option A+B)

The current `--in-memory` flag loads all data live from UTA PostgreSQL + SeqRepo at startup.
This takes 12-35 minutes because of three slow PostgreSQL views (`tx_similarity_v`,
`tx_def_summary_v`, `tx_exon_aln_v`) running under Docker/Rosetta emulation.

Plan:
1. **Add `ferro-benchmark prepare in-memory-cache` command** (or extend existing `prepare hgvs-rs`)
   - Replace slow `_v` view queries with direct table JOINs (Option B) to speed up loading
   - Skip `tx_similarity_v` entirely — `get_similar_transcripts` is never called during normalization
   - Serialize `InMemoryProvider` to a binary cache file via `bincode`
   - Estimated prepare time: ~1-2 min (down from 12-35 min)

2. **Update `--in-memory` to deserialize from cache file**
   - Deserialize at benchmark time: ~2-5s
   - Zero DB I/O during benchmarking

3. **Prerequisite: hgvs-rs serde derives**
   - Draft PR open: varfish-org/hgvs-rs#266 (branch `nh/serde-provider-records` on nh13 fork)
   - Adds `Serialize, Deserialize` to 7 Provider record types + enables chrono/serde
   - If not merged in time, can use a git dependency or wrapper types as workaround

4. **Run smoke test** (Task 6) once cache is working

## Key files

- `src/benchmark/inmemory_provider.rs` — InMemoryProvider struct + Provider trait impl
- `src/benchmark/hgvs_rs.rs` — normalizer pool + with_provider + in_memory flag
- `src/bin/benchmark.rs` — CLI --in-memory flag
- `docs/superpowers/plans/2026-03-27-in-memory-hgvs-rs-provider.md` — original plan

## Related

- ferro-hgvs issue #18 — hgvs-rs benchmark investigation
- varfish-org/hgvs-rs#265 — filed issue about production provider recommendations
- varfish-org/hgvs-rs#266 — draft PR for serde derives
- hgvs-rs local branch: `~/work/git/hgvs-rs` on `nh/serde-provider-records`
