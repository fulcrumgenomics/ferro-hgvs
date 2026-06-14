# Scoping: pre-stored transcript sequences (eliminate the runtime FASTA fetch+strip)

## Problem

In a profile of `ferro normalize` over a 500k-variant ClinVar corpus (after the
teardown / preprocessor / FxHashMap / reserve-maps wins), the single largest
remaining steady-state hot spot is `MultiFastaProvider::get_sequence_from_index`
(~14% self / ~20% inclusive), dominated by the **full-transcript fetch on a
transcript-cache miss**.

### Current path (per transcript-cache miss)

`get_transcript_on_build_inner` → `get_sequence_from_index(entry, 0, length)`:

1. `pread` the **line-wrapped** transcript bytes from the prepared, uncompressed
   `transcripts/*.fna` (indexed by `.fai`: offset, `line_bases`, `line_bytes`).
2. Allocate a read buffer (`vec![0u8; n]` — alloc + zero).
3. Walk every byte: strip `\n`/`\r`, `to_ascii_uppercase`, push into a second
   `Vec` (alloc).
4. `String::from_utf8(..)` — O(n) validation.

So each cold transcript pays: a buffer alloc+zero, a full byte-walk that
strips + uppercases into a second alloc, and a UTF-8 validation — over the whole
transcript (RefSeq RNA averages ~3 kb; 273,423 transcripts total). The
`transcript_cache` LRU means this only happens on a **miss**; an attempted
buffer-reuse micro-opt measured ~0% because the buffer was never the cost — the
strip + uppercase + validate are, and they are largely irreducible *given the
current on-disk format* (newline-wrapped, mixed-case `.fna`).

Note: the cost is workload-dependent. A transcript-**sorted** input is ~all cache
hits (≈0 fetches); an unsorted / random-access annotation stream (the realistic
case) misses often, so this lever's value scales with transcript locality.

## Idea

Do the strip + uppercase + validate **once, at `ferro prepare` time**, and store
the clean sequence so the runtime fetch is a near-free slice instead of a
read-strip-uppercase-validate.

## Constraints (from the codebase)

- `#![forbid(unsafe_code)]` style — **no `str::from_utf8_unchecked`**. A runtime
  `&str` view of mmap'd bytes still costs an `O(n)` `from_utf8` validation (fast,
  SIMD, on clean ASCII — but not free). A zero-copy `&str` with *no* validation
  is only available through rkyv's `ArchivedString` (valid-by-construction) or
  `unsafe` (excluded).
- `memmap2 = "0.9"` is already an **unconditional** dependency (the cdot `.rkyv`
  loader mmaps), so an mmap'd sequence store adds no new dependency and no
  resident-RAM baseline (OS pages on demand).
- The cdot `.rkyv` archive currently stores transcript **metadata only**
  (`RkyvTx`: exons, CDS bounds, strand, contig — no sequence). Sequence lives
  exclusively in the FASTA. The runtime `Transcript { sequence: Option<String> }`
  is filled from the FASTA on miss.
- Genome (`g.`) and protein (`p.`) sequence fetches also go through
  `get_sequence_from_index`; **only the per-transcript full fetch** is in scope —
  genome windows stay on the FASTA path.

## Options

### Option A — Sidecar pre-cleaned sequence blob + offset index (recommended)

`ferro prepare` writes one extra artifact, e.g. `transcripts/transcripts.seqblob`:
all transcript sequences, **newline-free + uppercased + UTF-8-validated**,
concatenated; plus a per-transcript `(offset: u64, len: u32)` index. Store the
index either as new fields on `RkyvTx` (bumps the cdot rkyv `format_version`) or
as a small sidecar (`transcript_count × 12 B ≈ 3 MB`).

Runtime: mmap the blob once (held by `MultiFastaProvider`); the transcript fetch
becomes `&blob[offset .. offset+len]` → `from_utf8` (validate) → done. No pread of
wrapped bytes, no strip, no uppercase, no read-buffer alloc.

- **Win:** removes the strip loop, the uppercase, the read-buffer alloc, and the
  newline bytes from the read. `from_utf8` remains. Estimate **~5–8%** of normalize
  wall on a cache-miss-heavy workload (roughly the strip/uppercase/alloc share of
  the ~14%); near-0 on a fully-sorted workload.
- **Memory:** mmap → on-demand paging, ~0 resident baseline; index ~3 MB.
- **Disk:** the clean blob is ~the size of the uncompressed transcript sequence
  (~0.8–1 GB). It can largely *replace* the uncompressed `.fna` for transcripts
  (those exist only to be read+stripped at runtime), so net disk may be flat.
- **Risk:** moderate. New artifact + format bump + a runtime fetch branch.
  Backward-compatible: if the blob is absent (older prepared reference), fall back
  to the existing FASTA path — graceful degradation.

### Option B — Sequences in the cdot rkyv as `ArchivedString`, queried zero-copy

Add the sequence to `RkyvTx` and query the archived map **in place** (rkyv
zero-copy), so the runtime sequence is `archived_tx.sequence.as_str()` — a true
zero-copy `&str` with **no `from_utf8`** (ArchivedString is valid by
construction).

- **Win:** larger — eliminates fetch + strip + validate entirely.
- **Cost:** this *requires* the rkyv zero-copy refactor (stop rebuilding the
  runtime `CdotMaps` in `maps_from_archived`; query `ArchivedHashMap`/archived
  transcripts directly), which touches the whole reference query layer (archived
  vs native types, mmap lifetimes threaded through every `get_*`). It also moves
  ~1 GB of sequence into the archive (on-demand paged via the existing `.rkyv`
  mmap, so RAM is OK, but the archive ~doubles on disk).
- **Risk:** high / large. This is the "big architectural change" flagged earlier;
  Option A is a decoupled subset of its benefit.

### Option C — Pre-unwrap + uppercase the on-disk `.fna` (minimal change)

`prepare` writes the transcript `.fna` unwrapped (one line per sequence) and
uppercased; the `.fai` then has `line_bases == length`. Runtime can skip the
newline-strip (no `\n` in range) and the uppercase.

- **Win:** removes the strip + uppercase, smaller than A (still allocates the read
  buffer + `from_utf8`, still preads from the per-FASTA files via the mutex'd
  handle cache).
- **Risk:** low code, but changes the on-disk `.fna`/`.fai` format and any other
  consumer that assumes wrapped FASTA; muddies "standard FASTA" expectations.

## Recommendation

**Option A.** It captures most of the realizable win (strip + uppercase + alloc +
newline-read removal) as a self-contained, backward-compatible change, decoupled
from the large rkyv zero-copy refactor (Option B). Reach for B only as part of a
broader reference-layer rework, and treat Option A's blob as a stepping stone.

## Work breakdown (Option A)

1. **`prepare`**: after decompressing transcript FASTAs, stream each sequence →
   strip newlines, uppercase, validate UTF-8 → append to `transcripts.seqblob`,
   recording `(id → offset, len)`. Add the blob to the manifest. (~½ day)
2. **Index storage**: add `seq_offset: u64`, `seq_len: u32` to `RkyvTx`/`CdotTx`
   runtime maps; bump cdot `format_version` (old caches self-heal/regenerate).
   Alternatively a sidecar index file. (~½ day)
3. **Runtime load**: `MultiFastaProvider` mmaps the blob (memmap2) when present;
   store `Option<Mmap>`. (~¼ day)
4. **Runtime fetch**: in `get_transcript_on_build_inner`, when the blob + offset
   are present, build the sequence from `&blob[off..off+len]` (`from_utf8`);
   else fall back to `get_sequence_from_index`. Keep windowed `get_sequence` for
   transcripts routed through the blob too (slice sub-ranges). (~½ day)
5. **Manifest / `ferro check`**: surface the blob; validate presence/coherence.
   (~¼ day)
6. **Tests + parity**: unit tests (blob round-trips a few transcripts); the
   **500k normalize golden must stay byte-identical**; a differential check that
   blob-sourced sequences equal FASTA-sourced for a transcript sample;
   `prepare` regenerates the test reference. (~½ day)
7. **Benchmark**: cache-miss-heavy (unsorted) A/B + reprofile to confirm the
   `get_sequence_from_index` share drops. (~¼ day)

**Estimate:** ~2.5–3 focused days for a reviewable PR.

## Parity & risk notes

- **Byte-identical sequences** vs the FASTA strip path are the invariant
  (soft-masked lowercase → uppercased identically; ambiguity codes preserved;
  exact transcript bounds). The 500k golden + a FASTA-vs-blob differential are the
  oracles.
- **Backward compatibility**: absent blob ⇒ fall back to FASTA. No hard break.
- **mmap lifetime**: the `Mmap` is owned by the provider; sequence `&str`/`&[u8]`
  borrow from it — lifetimes already managed the same way for the cdot `.rkyv`
  mmap, so this is established practice here.
- **Out of scope**: genome/protein fetches; the rkyv zero-copy refactor (Option B).
