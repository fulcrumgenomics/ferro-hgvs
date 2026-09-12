//! A 2-bit packed representation of decoded reference bases.
//!
//! The reference genome is stored on disk as FASTA text, and every base access
//! today pays a `decode_range` pass — newline-strip + ASCII-uppercase over the
//! raw bytes — which profiling puts at ~76% of `ferro normalize` runtime. That
//! decode is pure function of the FASTA and can be done **once**, offline, and
//! stored ready-to-use, exactly as bwa's `.pac` does.
//!
//! This module is that store's core: pack a record's already-decoded bases
//! (uppercased, newline-free — i.e. `decode_range` output) into 2 bits per base,
//! and read any range back **byte-identically**. Bases outside `A/C/G/T` (`N`
//! runs, rare IUPAC ambiguity codes) cannot be 2-bit encoded, so they pack as a
//! placeholder and their literal bytes are kept in a run-length **ambiguity
//! side-table** that the reader overlays — so `N` and IUPAC reproduce verbatim.
//!
//! `decode_range` (in `multi_fasta`) remains the correctness oracle: a read from
//! this store must equal a `decode_range` read of the same range, for every range.

use memmap2::Mmap;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

/// File magic for the on-disk store (`.pac`-style sidecar).
const MAGIC: &[u8; 8] = b"FERROSEQ";
/// On-disk format version; bumped on any layout change so a stale sidecar is
/// rejected rather than misread. v2 changed the fingerprint semantics (the
/// provider now folds source-file size+mtime into it, see
/// [`crate::reference::multi_fasta`]) and hardened `open` to fail closed on a
/// malformed directory, so a v1 sidecar is rejected and rebuilt.
const FORMAT_VERSION: u32 = 2;

/// 2-bit code for a base, `A=0 C=1 G=2 T=3` (bwa's mapping). Anything else maps
/// to `0` as a placeholder and is recorded as a hole (see [`RecordSpan::holes`]).
#[inline]
fn code_of(base: u8) -> u8 {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0,
    }
}

/// The base a 2-bit code decodes to. Inverse of [`code_of`] on `A/C/G/T`.
const CODE_TO_BASE: [u8; 4] = *b"ACGT";

/// For each packed byte, its four decoded bases in stream order — base `p` of the
/// byte lives at bits `p*2`, exactly as the per-base unpack reads them. Lets
/// [`read_range`] decode a whole aligned byte with one table lookup and a 4-byte
/// copy instead of four shift/mask/index/push steps.
const PACKED_TO_BASES: [[u8; 4]; 256] = {
    let mut table = [[0u8; 4]; 256];
    let mut b = 0usize;
    while b < 256 {
        table[b] = [
            CODE_TO_BASE[b & 0b11],
            CODE_TO_BASE[(b >> 2) & 0b11],
            CODE_TO_BASE[(b >> 4) & 0b11],
            CODE_TO_BASE[(b >> 6) & 0b11],
        ];
        b += 1;
    }
    table
};

/// A fingerprint over a set of `(record name, length)` pairs, tying a written
/// store to the *record set* of the source index it was built from. FNV-1a per
/// record, combined by wrapping addition — which is commutative and associative,
/// so the result is genuinely independent of record iteration order (a property
/// [`fingerprint_is_order_independent`] pins).
///
/// This detects a record being added, removed, renamed, or changed in length —
/// but **not** a same-length change to a record's bases (a corrected base, a
/// soft-mask flip). Detecting that cheaply at open time is impossible without
/// re-reading the FASTA (the decode this store exists to avoid), so the provider
/// (`MultiFastaProvider::sequence_fingerprint`) additionally folds each source
/// FASTA file's size and mtime into the stamp it writes and checks — any in-place
/// rewrite bumps the mtime, so the sidecar is rejected and the text path used.
/// A stale fingerprint at open time is a fall-back to the text path, never a wrong
/// answer.
pub(crate) fn fingerprint(records: &[(&str, u64)]) -> u64 {
    let mut acc: u64 = 0;
    for (name, length) in records {
        acc = acc.wrapping_add(record_hash(name, *length));
    }
    acc
}

/// Per-record FNV-1a over `name`, a NUL separator, and the little-endian length.
/// The NUL keeps `("AB", 1)` from colliding with `("A", …)` on a shared tail.
fn record_hash(name: &str, length: u64) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in name.as_bytes() {
        h = (h ^ b as u64).wrapping_mul(0x100000001b3);
    }
    // FNV step for the NUL separator (XOR with 0 is the identity, so only the
    // multiply remains).
    h = h.wrapping_mul(0x100000001b3);
    for &b in &length.to_le_bytes() {
        h = (h ^ b as u64).wrapping_mul(0x100000001b3);
    }
    h
}

/// Fold one more `u64` into a fingerprint commutatively, for the provider to mix
/// in per-source-file freshness stamps (size, mtime) on top of the record set.
pub(crate) fn fingerprint_mix(acc: u64, value: u64) -> u64 {
    acc.wrapping_add(value.wrapping_mul(0x100000001b3))
}

/// A maximal run of a single repeated non-`ACGT` byte, so it reproduces verbatim.
/// Covers record-relative positions `[pos, pos + len)`, every one equal to `byte`.
///
/// Run-length (one byte per run, not per position) keeps GRCh38's multi-megabase
/// `N` runs to a few bytes each; a maximal run is split at byte-value changes, so
/// a rare mixed IUPAC stretch (`…NRY…`) becomes several single-byte runs.
#[derive(Debug, Clone)]
struct Hole {
    pos: u64,
    len: u64,
    byte: u8,
}

/// One record's placement in the packed stream, plus its ambiguity holes.
#[derive(Debug, Clone)]
struct RecordSpan {
    /// The record's name (FASTA accession), used to look it up by the provider.
    name: String,
    /// Starting base index of this record within the concatenated packed stream.
    base_offset: u64,
    /// Number of bases in the record.
    length: u64,
    /// Non-`ACGT` runs, record-relative, sorted by `pos`, non-overlapping.
    holes: Vec<Hole>,
}

/// Read record `span`'s bases over `[start, end)` out of the 2-bit `packed`
/// stream, byte-identical to the decoded input. Shared by the in-memory
/// [`PackedSequence`] and the mmap-backed [`MappedSequence`].
fn read_range(packed: &[u8], span: &RecordSpan, start: u64, end: u64) -> Vec<u8> {
    debug_assert!(
        end <= span.length,
        "range end {end} past record length {}",
        span.length
    );
    if end <= start {
        return Vec::new();
    }
    // Unpack the 2-bit codes for [start, end) into A/C/G/T. Decode the aligned
    // interior a whole byte (4 bases) at a time via the lookup table; the leading
    // and trailing partial bytes fall back to the per-base read. Byte-identical to
    // a pure per-base loop (same codes, same CODE_TO_BASE) — the round-trip test
    // exercises every start/end offset, including mid-byte edges.
    let mut out = Vec::with_capacity((end - start) as usize);
    let gend = span.base_offset + end;
    let mut g = span.base_offset + start;
    let base = |g: u64| CODE_TO_BASE[((packed[(g / 4) as usize] >> ((g % 4) * 2)) & 0b11) as usize];
    // Leading partial byte, up to the next 4-base boundary.
    while g < gend && !g.is_multiple_of(4) {
        out.push(base(g));
        g += 1;
    }
    // Whole aligned bytes.
    while g + 4 <= gend {
        out.extend_from_slice(&PACKED_TO_BASES[packed[(g / 4) as usize] as usize]);
        g += 4;
    }
    // Trailing partial byte.
    while g < gend {
        out.push(base(g));
        g += 1;
    }
    // Overlay any ambiguity run intersecting [start, end) with its literal byte,
    // so non-ACGT bases reproduce exactly rather than as the A placeholder they
    // packed to.
    for hole in &span.holes {
        let hole_end = hole.pos + hole.len;
        if hole.pos >= end || hole_end <= start {
            continue;
        }
        let lo = hole.pos.max(start);
        let hi = hole_end.min(end);
        for p in lo..hi {
            out[(p - start) as usize] = hole.byte;
        }
    }
    out
}

/// A 2-bit packed set of decoded records with byte-identical range read-back.
#[derive(Debug, Clone)]
pub(crate) struct PackedSequence {
    /// 2-bit codes, 4 bases/byte, records concatenated in build order. Base `k`
    /// of the stream lives in `packed[k / 4]` at bit `(k % 4) * 2`.
    packed: Vec<u8>,
    records: Vec<RecordSpan>,
}

impl PackedSequence {
    /// Build from each record's fully-decoded bases (uppercased, newline-free),
    /// synthesizing record names `rec{i}`. For the round-trip test; callers with
    /// real accessions use [`Self::from_named_records`].
    #[cfg(test)]
    pub(crate) fn from_records(records: &[&[u8]]) -> Self {
        let named: Vec<(String, &[u8])> = records
            .iter()
            .enumerate()
            .map(|(i, &b)| (format!("rec{i}"), b))
            .collect();
        let refs: Vec<(&str, &[u8])> = named.iter().map(|(n, b)| (n.as_str(), *b)).collect();
        Self::from_named_records(&refs)
    }

    /// Build from each record's name and fully-decoded bases (uppercased,
    /// newline-free — i.e. `decode_range` output over the whole record).
    pub(crate) fn from_named_records(records: &[(&str, &[u8])]) -> Self {
        let total_bases: u64 = records.iter().map(|(_, r)| r.len() as u64).sum();
        let mut packed = vec![0u8; total_bases.div_ceil(4) as usize];
        let mut spans = Vec::with_capacity(records.len());
        let mut base_offset: u64 = 0;

        for (name, bases) in records {
            let mut holes: Vec<Hole> = Vec::new();
            // The run in progress, extended only while the byte value repeats; a
            // different non-ACGT byte closes it and opens a new one.
            let mut run: Option<Hole> = None;
            for (i, &b) in bases.iter().enumerate() {
                let global = base_offset + i as u64;
                packed[(global / 4) as usize] |= code_of(b) << (((global % 4) * 2) as u8);
                let is_acgt = matches!(b, b'A' | b'C' | b'G' | b'T');
                if is_acgt {
                    if let Some(h) = run.take() {
                        holes.push(h);
                    }
                    continue;
                }
                match run.as_mut() {
                    Some(h) if h.byte == b => h.len += 1,
                    _ => {
                        if let Some(h) = run.take() {
                            holes.push(h);
                        }
                        run = Some(Hole {
                            pos: i as u64,
                            len: 1,
                            byte: b,
                        });
                    }
                }
            }
            if let Some(h) = run.take() {
                holes.push(h);
            }
            spans.push(RecordSpan {
                name: (*name).to_string(),
                base_offset,
                length: bases.len() as u64,
                holes,
            });
            base_offset += bases.len() as u64;
        }

        PackedSequence {
            packed,
            records: spans,
        }
    }

    /// Bases for record `r` over `[start, end)`, byte-identical to the decoded
    /// input the store was built from. `end` must be `<= record length`. The
    /// production read path uses [`MappedSequence::range`]; this in-memory reader
    /// exists for the round-trip test.
    #[cfg(test)]
    pub(crate) fn range(&self, r: usize, start: u64, end: u64) -> Vec<u8> {
        read_range(&self.packed, &self.records[r], start, end)
    }

    /// Serialize to `path` as the on-disk sidecar. `fingerprint` ties the store
    /// to the FASTA/`.fai` it was built from, so a later mismatch is rejected.
    ///
    /// Written atomically: the bytes go to a `<path>.tmp` sibling and are renamed
    /// into place only after a successful flush, so a crash or disk-full mid-write
    /// never leaves a partial `sequence_store.pac` at the final name (which `open`
    /// would then have to reject rather than a reader silently trusting it).
    ///
    /// Layout (all little-endian): magic, version, fingerprint, record_count,
    /// packed_len; then per record `name_len,name, base_offset,length,
    /// hole_count, (pos,len,byte)×hole_count`; then the packed bytes last.
    pub(crate) fn write(&self, path: &Path, fingerprint: u64) -> io::Result<()> {
        let tmp = {
            let mut s = path.as_os_str().to_owned();
            s.push(".tmp");
            PathBuf::from(s)
        };
        // Scope the writer so the file is closed (and flushed) before the rename.
        let result = (|| {
            let mut w = BufWriter::new(File::create(&tmp)?);
            w.write_all(MAGIC)?;
            w.write_all(&FORMAT_VERSION.to_le_bytes())?;
            w.write_all(&fingerprint.to_le_bytes())?;
            w.write_all(&(self.records.len() as u32).to_le_bytes())?;
            w.write_all(&(self.packed.len() as u64).to_le_bytes())?;
            for span in &self.records {
                let name = span.name.as_bytes();
                w.write_all(&(name.len() as u16).to_le_bytes())?;
                w.write_all(name)?;
                w.write_all(&span.base_offset.to_le_bytes())?;
                w.write_all(&span.length.to_le_bytes())?;
                w.write_all(&(span.holes.len() as u32).to_le_bytes())?;
                for h in &span.holes {
                    w.write_all(&h.pos.to_le_bytes())?;
                    w.write_all(&h.len.to_le_bytes())?;
                    w.write_all(&[h.byte])?;
                }
            }
            w.write_all(&self.packed)?;
            w.flush()
        })();
        match result {
            Ok(()) => std::fs::rename(&tmp, path),
            Err(e) => {
                // Best-effort cleanup of the partial temp file; the write error is
                // what the caller needs to see.
                let _ = std::fs::remove_file(&tmp);
                Err(e)
            }
        }
    }
}

/// A cursor over an in-memory byte slice, for parsing the sidecar directory.
struct Cursor<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn need(&self, n: usize) -> io::Result<()> {
        if self.pos + n > self.buf.len() {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "sequence store truncated",
            ));
        }
        Ok(())
    }
    fn take(&mut self, n: usize) -> io::Result<&'a [u8]> {
        self.need(n)?;
        let s = &self.buf[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }
    fn u16(&mut self) -> io::Result<u16> {
        Ok(u16::from_le_bytes(self.take(2)?.try_into().unwrap()))
    }
    fn u32(&mut self) -> io::Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn u64(&mut self) -> io::Result<u64> {
        Ok(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
    fn u8(&mut self) -> io::Result<u8> {
        Ok(self.take(1)?[0])
    }
}

/// An mmap-backed store opened from an on-disk sidecar. The directory (records +
/// holes) is parsed into memory; the packed 2-bit stream stays in the mmap and is
/// read as a zero-copy slice.
pub(crate) struct MappedSequence {
    _mmap: Mmap,
    /// Byte offset of the packed stream within the mmap.
    packed_offset: usize,
    packed_len: usize,
    records: Vec<RecordSpan>,
    by_name: HashMap<String, usize>,
}

impl MappedSequence {
    /// Open and validate the sidecar at `path`.
    ///
    /// Returns `Ok(None)` — *fall back to `decode_range`* — when the file is
    /// absent, its magic/version is unrecognized, or its `fingerprint` does not
    /// match `expected_fingerprint` (the FASTA changed under it). `Err` is a real
    /// I/O or corruption error.
    pub(crate) fn open(path: &Path, expected_fingerprint: u64) -> io::Result<Option<Self>> {
        let file = match File::open(path) {
            Ok(f) => f,
            Err(e) if e.kind() == io::ErrorKind::NotFound => return Ok(None),
            Err(e) => return Err(e),
        };
        // SAFETY: the sidecar is a ferro-prepared artifact; concurrent truncation
        // by another writer is out of scope, as it is for the FASTA it mirrors.
        let mmap = unsafe { Mmap::map(&file)? };
        let mut c = Cursor { buf: &mmap, pos: 0 };

        if c.take(8)? != MAGIC || c.u32()? != FORMAT_VERSION {
            return Ok(None);
        }
        let fingerprint = c.u64()?;
        if fingerprint != expected_fingerprint {
            return Ok(None);
        }
        let record_count = c.u32()? as usize;
        let packed_len = c.u64()? as usize;

        // Fail closed on a malformed-but-fingerprint-matching directory. The
        // counts and spans below are untrusted on-disk `u32`/`u64`s; without
        // these guards a garbage `record_count`/`hole_count` drives a
        // multi-gigabyte `Vec`/`HashMap` allocation before the cursor's `need()`
        // bound ever runs, and a span past the packed stream makes `read_range`
        // index out of bounds — a release panic (the only end guard there is a
        // `debug_assert`, compiled out) rather than the intended text-path
        // fallback. Every record needs ≥ MIN_RECORD_BYTES and every hole exactly
        // HOLE_BYTES on disk, so the mmap length bounds both counts.
        const MIN_RECORD_BYTES: usize = 2 + 8 + 8 + 4; // name_len, base_offset, length, hole_count (name may be empty)
        const HOLE_BYTES: usize = 8 + 8 + 1; // pos, len, byte
        let malformed = |what: &str| io::Error::new(io::ErrorKind::InvalidData, what);
        if record_count > mmap.len() / MIN_RECORD_BYTES {
            return Err(malformed("sequence store record count implausible"));
        }

        let mut records = Vec::with_capacity(record_count);
        let mut by_name = HashMap::with_capacity(record_count);
        for idx in 0..record_count {
            let name_len = c.u16()? as usize;
            let name = String::from_utf8(c.take(name_len)?.to_vec())
                .map_err(|_| malformed("sequence store holds a non-UTF8 record name"))?;
            let base_offset = c.u64()?;
            let length = c.u64()?;
            let hole_count = c.u32()? as usize;
            if hole_count > (mmap.len() - c.pos) / HOLE_BYTES {
                return Err(malformed("sequence store hole count implausible"));
            }
            let mut holes = Vec::with_capacity(hole_count);
            for _ in 0..hole_count {
                let pos = c.u64()?;
                let len = c.u64()?;
                let byte = c.u8()?;
                holes.push(Hole { pos, len, byte });
            }
            by_name.insert(name.clone(), idx);
            records.push(RecordSpan {
                name,
                base_offset,
                length,
                holes,
            });
        }

        let packed_offset = c.pos;
        if packed_offset + packed_len > mmap.len() {
            return Err(malformed("sequence store packed stream truncated"));
        }
        // Every record's bases must lie within the packed stream (4 bases/byte),
        // or read_range indexes past `packed`. Checked with overflow-safe
        // arithmetic since these are untrusted u64s.
        let capacity_bases = (packed_len as u64).checked_mul(4);
        for span in &records {
            let end = span.base_offset.checked_add(span.length);
            match (end, capacity_bases) {
                (Some(end), Some(cap)) if end <= cap => {}
                _ => {
                    return Err(malformed(
                        "sequence store record span exceeds packed stream",
                    ))
                }
            }
        }
        Ok(Some(MappedSequence {
            _mmap: mmap,
            packed_offset,
            packed_len,
            records,
            by_name,
        }))
    }

    /// Index of the record named `name`, if present.
    pub(crate) fn record_index(&self, name: &str) -> Option<usize> {
        self.by_name.get(name).copied()
    }

    /// Length in bases of record `r`.
    #[cfg(test)]
    pub(crate) fn record_len(&self, r: usize) -> u64 {
        self.records[r].length
    }

    /// Bases for record `r` over `[start, end)`, byte-identical to `decode_range`.
    pub(crate) fn range(&self, r: usize, start: u64, end: u64) -> Vec<u8> {
        let packed = &self._mmap[self.packed_offset..self.packed_offset + self.packed_len];
        read_range(packed, &self.records[r], start, end)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::ErrorKind;
    use tempfile::tempdir;

    /// The fingerprint is documented order-independent; the fold must actually be
    /// commutative, or a hash-map iteration-order difference between build and
    /// load spuriously invalidates a valid sidecar. Any permutation of the same
    /// record set folds to the same value, and a real change (a length edit) does
    /// not.
    #[test]
    fn fingerprint_is_order_independent() {
        let a = fingerprint(&[("chr1", 100), ("chr2", 50), ("chrM", 16569)]);
        let b = fingerprint(&[("chrM", 16569), ("chr1", 100), ("chr2", 50)]);
        let c = fingerprint(&[("chr2", 50), ("chrM", 16569), ("chr1", 100)]);
        assert_eq!(a, b, "permuted record order must fold identically");
        assert_eq!(a, c, "permuted record order must fold identically");
        // A same-name length change is detected.
        assert_ne!(
            a,
            fingerprint(&[("chr1", 101), ("chr2", 50), ("chrM", 16569)]),
            "a length change must change the fingerprint"
        );
        // A rename is detected.
        assert_ne!(
            a,
            fingerprint(&[("chr1", 100), ("chr2", 50), ("chrMT", 16569)]),
            "a rename must change the fingerprint"
        );
    }

    /// Assemble a minimal on-disk directory by hand so a malformed one can be
    /// tested without a real pack. `packed_len` and the single record's `length`
    /// are the knobs the fail-closed checks key on.
    fn write_raw_store(
        path: &Path,
        fp: u64,
        record_count: u32,
        packed_len: u64,
        record: Option<(&str, u64, u64)>, // (name, base_offset, length)
    ) {
        let mut b: Vec<u8> = Vec::new();
        b.extend_from_slice(MAGIC);
        b.extend_from_slice(&FORMAT_VERSION.to_le_bytes());
        b.extend_from_slice(&fp.to_le_bytes());
        b.extend_from_slice(&record_count.to_le_bytes());
        b.extend_from_slice(&packed_len.to_le_bytes());
        if let Some((name, base_offset, length)) = record {
            b.extend_from_slice(&(name.len() as u16).to_le_bytes());
            b.extend_from_slice(name.as_bytes());
            b.extend_from_slice(&base_offset.to_le_bytes());
            b.extend_from_slice(&length.to_le_bytes());
            b.extend_from_slice(&0u32.to_le_bytes()); // hole_count
        }
        std::fs::write(path, &b).unwrap();
    }

    /// A sidecar whose header/fingerprint match but whose directory is malformed
    /// must fail closed (an `Err`, which the provider turns into a text-path
    /// fallback) rather than panic or read out of bounds in a release build.
    #[test]
    fn open_fails_closed_on_a_malformed_directory() {
        let dir = tempdir().unwrap();
        let fp: u64 = 0x1234_5678_9ABC_DEF0;
        // MappedSequence is not Debug (it owns an Mmap), so assert on the result
        // without unwrap_err.
        let assert_invalid = |path: &Path, what: &str| match MappedSequence::open(path, fp) {
            Err(e) => assert_eq!(e.kind(), ErrorKind::InvalidData, "{what}"),
            Ok(_) => panic!("{what}: expected InvalidData, got Ok"),
        };

        // (a) A record span past the packed stream (length 100, capacity 0).
        let p1 = dir.path().join("span.pac");
        write_raw_store(&p1, fp, 1, 0, Some(("chr1", 0, 100)));
        assert_invalid(&p1, "span overflow must be rejected");

        // (b) An implausible record count (huge) against a tiny file — must be
        // rejected before any capacity allocation.
        let p2 = dir.path().join("count.pac");
        write_raw_store(&p2, fp, u32::MAX, 0, None);
        assert_invalid(&p2, "implausible record count must be rejected");

        // (c) A truncated packed stream (the record fits the declared base
        // capacity, but the packed bytes are not present) is still rejected.
        let p3 = dir.path().join("trunc.pac");
        write_raw_store(&p3, fp, 1, 1_000_000, Some(("chr1", 0, 4)));
        assert_invalid(&p3, "truncated packed stream must be rejected");
    }

    /// The persisted sidecar reopens byte-identically, resolves records by name,
    /// and rejects a stale fingerprint or a missing file by falling back (`None`)
    /// rather than erroring.
    #[test]
    fn persisted_store_reopens_byte_identically_and_rejects_wrong_fingerprint() {
        let rec0: &[u8] = b"ACGTACGTNNNNNRYACGTGGCC";
        let rec1: &[u8] = b"NNACGTACGTNN";
        let named: &[(&str, &[u8])] = &[("chr1", rec0), ("chr2", rec1)];
        let store = PackedSequence::from_named_records(named);

        let dir = tempdir().unwrap();
        let path = dir.path().join("ref.pac");
        let fp: u64 = 0xABCD_1234_5678_9ABC;
        store.write(&path, fp).unwrap();

        // Correct fingerprint => opens; every range matches the source bytes.
        let mapped = MappedSequence::open(&path, fp)
            .unwrap()
            .expect("store should open");
        for (r, (name, bases)) in named.iter().enumerate() {
            assert_eq!(mapped.record_index(name), Some(r), "record {name} index");
            assert_eq!(mapped.record_len(r), bases.len() as u64);
            let n = bases.len() as u64;
            for start in 0..=n {
                for end in start..=n {
                    assert_eq!(
                        mapped.range(r, start, end),
                        &bases[start as usize..end as usize],
                        "record {name} range [{start},{end})",
                    );
                }
            }
        }
        assert_eq!(mapped.record_index("absent"), None);

        // Stale fingerprint => fall back (None), not an error.
        assert!(
            MappedSequence::open(&path, fp ^ 1).unwrap().is_none(),
            "stale fingerprint"
        );
        // Absent sidecar => fall back (None).
        assert!(MappedSequence::open(&dir.path().join("missing.pac"), fp)
            .unwrap()
            .is_none());
    }

    /// The oracle: a read from the packed store must equal the decoded bytes it
    /// was built from, for every range — including ranges that start/end inside
    /// an `N` run or an IUPAC ambiguity code, empty ranges, and whole records.
    #[test]
    fn range_round_trips_decoded_bases_including_holes() {
        // ACGT, an interior N run, a stray IUPAC code, then ACGT again — the
        // shapes a real contig carries (telomere/centromere N, rare ambiguity).
        let rec0: &[u8] = b"ACGTACGTNNNNNRYACGTGGCC";
        let rec1: &[u8] = b"NNACGTACGTNN"; // leading and trailing holes
        let records: &[&[u8]] = &[rec0, rec1];
        let store = PackedSequence::from_records(records);

        for (r, bases) in records.iter().enumerate() {
            let n = bases.len() as u64;
            // Every (start, end) with start <= end <= n.
            for start in 0..=n {
                for end in start..=n {
                    let got = store.range(r, start, end);
                    let want = &bases[start as usize..end as usize];
                    assert_eq!(
                        got,
                        want,
                        "record {r} range [{start},{end}) mismatch: got {:?} want {:?}",
                        String::from_utf8_lossy(&got),
                        String::from_utf8_lossy(want),
                    );
                }
            }
        }
    }
}
