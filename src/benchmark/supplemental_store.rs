//! Upsert-by-accession provisioning for a supplemental transcript FASTA.
//!
//! # Why this module exists
//!
//! "Do we have sequence for this accession?" used to be answered by the
//! `*.metadata.json` sidecar, while the sequence itself lived in the FASTA and was
//! written in **append** mode. Two artifacts, updated non-atomically, answering one
//! question — so any run that left them disagreeing caused the next run to re-fetch an
//! accession whose sequence was already on disk and append a second copy of it. A FASTA
//! has no uniqueness constraint, so a duplicate header is legal FASTA and silently wrong
//! here: the file is being used as a keyed store with no key.
//!
//! # The rule this module establishes
//!
//! **The FASTA is the single authority for presence and sequence.** Everything else is
//! derived from it:
//!
//! | artifact | status | recovered how |
//! |---|---|---|
//! | `<name>.fna` | **authoritative** | — |
//! | `<name>.fna.fai` | derived | rebuilt from the committed FASTA on every commit |
//! | `<name>.metadata.json` | derived *key set*, carried-forward CDS annotation | keys reconciled against the FASTA on every commit; a missing row is synthesized from the FASTA record |
//!
//! The sidecar is not *fully* derived and this module does not pretend it is: CDS
//! coordinates come from a GenBank record and cannot be recovered from a FASTA. What is
//! derived is the sidecar's **key set** and each row's `sequence_length` — which is
//! precisely the part that was being asked the presence question. A row whose CDS
//! annotation already exists is carried forward untouched.
//!
//! # Why the FASTA and not the `.fai`
//!
//! **This repository already decided this question, for the same defect on the neighbouring
//! artifact.** `prepare::run_backfill` reads its present set from the backfill FASTA's
//! headers and says why in so many words (`src/prepare/mod.rs`, at the
//! `load_fasta_header_accessions` call): "not its derived `.fai`: a prior run interrupted
//! after the append but before reindex — or a deleted `.fai` — would otherwise look empty
//! and re-fetch and re-append the same accession, duplicating its `>accession` record".
//! That is this defect exactly, and the rule it settled on is the one applied here rather
//! than a second rule invented alongside it.
//!
//! The reasoning holds independently: the `.fai` is a *second* artifact that can drift from
//! the FASTA, which is what the sidecar was. On the operator's shipped reference it is seven
//! minutes newer than the `.fna` it indexes. Nominating it as the authority inverts the
//! dependency and reintroduces the defect class one level over.
//!
//! Scanning the FASTA's headers instead costs one sequential read — measured at **2.08 s
//! over the 1.3 GB shipped artifact**, against a provisioning run whose fetch half is a
//! networked, rate-limited operation measured in hours.
//!
//! Writing the index here is therefore load-bearing rather than tidy: the merge re-emits
//! every record at [`FASTA_LINE_BASES`], so a commit invalidates any externally-built index
//! it does not replace.
//!
//! # Atomicity boundary
//!
//! Fetched records are appended to a **staging** file (`<name>.fna.incoming`), which is
//! durable and resumable: an interrupted run's staged records count as present on the next
//! run, so no accession is fetched twice. [`SupplementalStore::commit`] then streams the
//! committed FASTA and the staging file into `<name>.fna.tmp`, keyed by accession — one
//! record per accession, staged records replacing committed ones — and `fsync`s and
//! renames it into place.
//!
//! **That rename is the commit point, and it is the only artifact that needs one.** The
//! index and the sidecar are rewritten after it, and a crash between the rename and either
//! of those is self-healing: the next `open` re-derives both from the FASTA. Making the
//! derived artifacts genuinely derived is what buys atomicity across three files without
//! needing a three-file atomic write, which is not available.
//!
//! # Idempotence
//!
//! [`SupplementalStore::is_dirty`] is the fixed-point test: no staged records, no
//! duplicate accession in the FASTA, and a sidecar whose key set already equals the
//! FASTA's. When it is false the store writes **nothing**, so a second provisioning run
//! over an artifact the first one produced leaves every byte alone. When it is true — the
//! shipped duplicate-bearing artifact is the case that matters — the commit is also the
//! repair: dedup-by-accession is not a special mode, it is what the merge does.
//!
//! **Start the fixed-point test from a DAMAGED artifact, not a clean one.** Measured by
//! reconstructing the previous behaviour and re-running these tests: run-twice-from-empty
//! **passes** on the old code, because on an artifact the old code itself produced, the
//! FASTA and the sidecar agree and its sidecar-keyed resume filter therefore skips
//! correctly. That agreement is the coincidence this module exists to stop relying on, so a
//! fixed-point test seeded from a clean state is satisfied *by* the bug. The tests that
//! fail on the old behaviour are the ones seeded from a state where the two artifacts
//! disagree — a FASTA record with no sidecar row, and a duplicate-bearing FASTA.

use crate::benchmark::cache::{SupplementalMetadata, SupplementalTranscript};
use crate::FerroError;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::ffi::OsString;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

/// Bases per sequence line in every FASTA record this module writes.
///
/// Matches the width the previous append-based writer used, so a committed record is
/// byte-identical to the one the old path would have written.
const FASTA_LINE_BASES: usize = 70;

/// Suffix of the durable, resumable file that fetched records are appended to.
const STAGING_SUFFIX: &str = ".incoming";

/// Suffix of the transient file that a commit writes before renaming it into place.
const TEMP_SUFFIX: &str = ".tmp";

/// Return `path` with `suffix` appended to its file name.
///
/// Deliberately not `Path::with_extension`, which would *replace* `.fna` rather than
/// extend it — `patterns_transcripts.fna.with_extension("tmp")` is
/// `patterns_transcripts.tmp`, which collides across sibling artifacts.
fn sibling_with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut name = OsString::from(path.as_os_str());
    name.push(suffix);
    PathBuf::from(name)
}

/// The sidecar path for a supplemental FASTA.
///
/// `patterns_transcripts.fna` -> `patterns_transcripts.metadata.json`. This *does* replace
/// the extension, matching the reader in [`crate::reference::multi_fasta`] and the
/// historical writer, so the name is not this module's to choose.
fn metadata_path_for(fasta_path: &Path) -> PathBuf {
    fasta_path.with_extension("metadata.json")
}

/// One record of a supplemental FASTA.
///
/// Deliberately not `pub`: nothing outside this module constructs or inspects one, and
/// the store's public surface is the accession-keyed questions, not the record type.
#[derive(Debug, Clone, PartialEq, Eq)]
struct FastaRecord {
    /// Accession as it appears first on the header line, without the leading `>`.
    pub accession: String,
    /// Remainder of the header line, possibly empty.
    pub description: String,
    /// Sequence with all line breaks removed.
    pub sequence: String,
}

impl FastaRecord {
    /// Gene symbol as the sidecar records it for a record with no GenBank annotation.
    ///
    /// Mirrors the CON/FASTA-fallback insert exactly: the description up to its first
    /// `(`, trimmed, and `None` when there is no description at all.
    fn gene_symbol(&self) -> Option<String> {
        if self.description.is_empty() {
            None
        } else {
            self.description
                .split('(')
                .next()
                .map(|s| s.trim().to_string())
        }
    }

    /// The sidecar row this record implies on its own, carrying no CDS annotation.
    fn synthesized_metadata(&self) -> SupplementalTranscript {
        SupplementalTranscript {
            id: self.accession.clone(),
            gene_symbol: self.gene_symbol(),
            cds_start: None,
            cds_end: None,
            sequence_length: self.sequence.len(),
        }
    }
}

/// Streaming FASTA reader yielding one [`FastaRecord`] at a time.
///
/// Streaming rather than slurping because the artifact this runs against is 1.3 GB.
struct FastaReader<R: BufRead> {
    reader: R,
    /// Header line already read but belonging to the *next* record.
    pending_header: Option<String>,
    /// Set once the underlying reader is exhausted.
    exhausted: bool,
}

impl<R: BufRead> FastaReader<R> {
    fn new(reader: R) -> Self {
        Self {
            reader,
            pending_header: None,
            exhausted: false,
        }
    }

    /// Read the next record, or `None` at end of input.
    ///
    /// A header with no sequence is skipped rather than returned: this module's whole
    /// premise is that a record in the FASTA means bases on disk, and an empty record
    /// would be a presence claim with nothing behind it.
    fn next_record(&mut self) -> Result<Option<FastaRecord>, FerroError> {
        loop {
            let header = match self.take_header()? {
                Some(header) => header,
                None => return Ok(None),
            };

            let mut sequence = String::new();
            loop {
                let mut line = String::new();
                let read = self
                    .reader
                    .read_line(&mut line)
                    .map_err(|e| FerroError::Io {
                        msg: format!("Failed to read FASTA: {}", e),
                    })?;
                if read == 0 {
                    self.exhausted = true;
                    break;
                }
                let trimmed = line.trim_end_matches(['\n', '\r']);
                if let Some(header) = trimmed.strip_prefix('>') {
                    self.pending_header = Some(header.to_string());
                    break;
                }
                sequence.push_str(trimmed.trim());
            }

            let (accession, description) = split_header(&header);
            if accession.is_empty() || sequence.is_empty() {
                continue;
            }
            return Ok(Some(FastaRecord {
                accession,
                description,
                sequence,
            }));
        }
    }

    /// Return the next header line's contents (without `>`), or `None` at end of input.
    fn take_header(&mut self) -> Result<Option<String>, FerroError> {
        if let Some(header) = self.pending_header.take() {
            return Ok(Some(header));
        }
        loop {
            if self.exhausted {
                return Ok(None);
            }
            let mut line = String::new();
            let read = self
                .reader
                .read_line(&mut line)
                .map_err(|e| FerroError::Io {
                    msg: format!("Failed to read FASTA: {}", e),
                })?;
            if read == 0 {
                self.exhausted = true;
                return Ok(None);
            }
            let trimmed = line.trim_end_matches(['\n', '\r']);
            if let Some(header) = trimmed.strip_prefix('>') {
                return Ok(Some(header.to_string()));
            }
        }
    }
}

/// Split a FASTA header body into its accession and the rest.
fn split_header(header: &str) -> (String, String) {
    match header.split_once(char::is_whitespace) {
        Some((accession, rest)) => (accession.to_string(), rest.trim().to_string()),
        None => (header.trim().to_string(), String::new()),
    }
}

/// Open a FASTA for streaming, treating an absent file as empty input.
fn open_fasta(path: &Path) -> Result<Option<FastaReader<BufReader<File>>>, FerroError> {
    match File::open(path) {
        Ok(file) => Ok(Some(FastaReader::new(BufReader::new(file)))),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
        Err(e) => Err(FerroError::Io {
            msg: format!("Failed to open {}: {}", path.display(), e),
        }),
    }
}

/// Line terminating a completely-written record in the staging file.
///
/// The staging file is this module's own format, not a FASTA anyone else reads, so it can
/// carry a completeness marker that a FASTA cannot. A `;`-prefixed line is the classic
/// FASTA comment and cannot collide with a header (`>`) or with a sequence line, which
/// [`write_record`] emits as bases only.
const STAGED_RECORD_TERMINATOR: &str = ";end";

/// Streaming reader over the staging file that yields only completely-written records.
///
/// # Why this exists
///
/// Without it, a crash part-way through [`SupplementalStore::stage_record`] leaves a
/// truncated trailing record, and the next run would count its accession as present —
/// so the accession is never re-fetched and the *truncated* sequence is merged into the
/// authoritative FASTA. A torn write silently becoming the authority is the exact defect
/// class this module exists to close, so the staging file must not reintroduce it one
/// level down.
///
/// A terminator can only be complete if everything before it was written, so its presence
/// is a sound completeness witness for an append-only sequential writer.
struct StagedReader<R: BufRead> {
    reader: R,
    /// Header of the record currently being accumulated.
    header: Option<String>,
    sequence: String,
    /// Records dropped because they were never terminated.
    torn: usize,
}

impl<R: BufRead> StagedReader<R> {
    fn new(reader: R) -> Self {
        Self {
            reader,
            header: None,
            sequence: String::new(),
            torn: 0,
        }
    }

    /// Read the next completely-written record, or `None` at end of input.
    fn next_record(&mut self) -> Result<Option<FastaRecord>, FerroError> {
        loop {
            let mut line = String::new();
            let read = self
                .reader
                .read_line(&mut line)
                .map_err(|e| FerroError::Io {
                    msg: format!("Failed to read staging file: {}", e),
                })?;
            if read == 0 {
                // An in-progress record at end of input was never terminated.
                if self.header.take().is_some() {
                    self.torn += 1;
                }
                self.sequence.clear();
                return Ok(None);
            }
            let trimmed = line.trim_end_matches(['\n', '\r']);
            if let Some(header) = trimmed.strip_prefix('>') {
                // A pending header here means the previous record never terminated.
                if self.header.replace(header.to_string()).is_some() {
                    self.torn += 1;
                }
                self.sequence.clear();
            } else if trimmed == STAGED_RECORD_TERMINATOR {
                if let Some(header) = self.header.take() {
                    let (accession, description) = split_header(&header);
                    let sequence = std::mem::take(&mut self.sequence);
                    if !accession.is_empty() && !sequence.is_empty() {
                        return Ok(Some(FastaRecord {
                            accession,
                            description,
                            sequence,
                        }));
                    }
                }
            } else {
                self.sequence.push_str(trimmed.trim());
            }
        }
    }
}

/// Open the staging file for streaming, treating an absent file as empty input.
fn open_staged(path: &Path) -> Result<Option<StagedReader<BufReader<File>>>, FerroError> {
    match File::open(path) {
        Ok(file) => Ok(Some(StagedReader::new(BufReader::new(file)))),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
        Err(e) => Err(FerroError::Io {
            msg: format!("Failed to open {}: {}", path.display(), e),
        }),
    }
}

/// Accessions of every completely-written record in the staging file.
fn scan_staged_accessions(path: &Path) -> Result<HashSet<String>, FerroError> {
    let mut accessions = HashSet::new();
    let Some(mut reader) = open_staged(path)? else {
        return Ok(accessions);
    };
    while let Some(record) = reader.next_record()? {
        accessions.insert(record.accession);
    }
    if reader.torn > 0 {
        eprintln!(
            "  Discarding {} incompletely-written record(s) from {}; they will be re-fetched",
            reader.torn,
            path.display()
        );
    }
    Ok(accessions)
}

/// Count records per accession in a FASTA, treating an absent file as empty.
///
/// A count above one is a duplicated accession — the damage this module repairs.
fn scan_accession_counts(path: &Path) -> Result<BTreeMap<String, usize>, FerroError> {
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    let Some(mut reader) = open_fasta(path)? else {
        return Ok(counts);
    };
    while let Some(record) = reader.next_record()? {
        *counts.entry(record.accession).or_insert(0) += 1;
    }
    Ok(counts)
}

/// Write one record and return the `.fai` line describing where it landed.
///
/// `offset` is advanced by the bytes written, so a caller streaming records into a fresh
/// file can build an exact index without a second pass over it.
fn write_record(
    writer: &mut impl Write,
    record: &FastaRecord,
    offset: &mut u64,
) -> Result<String, FerroError> {
    let io_err = |e: std::io::Error| FerroError::Io {
        msg: format!("Failed to write FASTA record {}: {}", record.accession, e),
    };

    let header = if record.description.is_empty() {
        format!(">{}\n", record.accession)
    } else {
        format!(">{} {}\n", record.accession, record.description)
    };
    writer.write_all(header.as_bytes()).map_err(io_err)?;
    *offset += header.len() as u64;

    let sequence_offset = *offset;
    for chunk in record.sequence.as_bytes().chunks(FASTA_LINE_BASES) {
        writer.write_all(chunk).map_err(io_err)?;
        writer.write_all(b"\n").map_err(io_err)?;
        *offset += chunk.len() as u64 + 1;
    }

    // samtools' `.fai` describes the geometry of the *first* sequence line, which for a
    // record shorter than one full line is the record itself.
    let line_bases = record.sequence.len().min(FASTA_LINE_BASES);
    Ok(format!(
        "{}\t{}\t{}\t{}\t{}\n",
        record.accession,
        record.sequence.len(),
        sequence_offset,
        line_bases,
        line_bases + 1
    ))
}

/// What a commit did, for the caller to report.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct CommitReport {
    /// Records in the committed FASTA afterwards — one per accession, by construction.
    pub records: usize,
    /// Duplicate records dropped by this commit.
    pub duplicates_removed: usize,
    /// Committed records replaced by a freshly staged one for the same accession.
    pub records_replaced: usize,
    /// Sidecar rows dropped because the FASTA has no record for their accession.
    pub metadata_rows_dropped: usize,
    /// Sidecar rows synthesized from a FASTA record that had none.
    pub metadata_rows_synthesized: usize,
}

/// A supplemental FASTA plus the derived artifacts kept consistent with it.
///
/// See the [module documentation](self) for the authority rule and the atomicity boundary.
pub struct SupplementalStore {
    fasta_path: PathBuf,
    staging_path: PathBuf,
    metadata_path: PathBuf,
    /// Records per accession in the committed FASTA.
    committed: BTreeMap<String, usize>,
    /// Accessions already written to the staging file, committed or not.
    staged: HashSet<String>,
    /// Opened lazily, so a run that stages nothing creates no staging file.
    staging_writer: Option<BufWriter<File>>,
}

impl SupplementalStore {
    /// Open the store for `fasta_path`, scanning it and any staging file left behind.
    ///
    /// Neither file needs to exist; both are treated as empty when absent, which is what
    /// a genuine first run looks like.
    pub fn open(fasta_path: &Path) -> Result<Self, FerroError> {
        let staging_path = sibling_with_suffix(fasta_path, STAGING_SUFFIX);
        let committed = scan_accession_counts(fasta_path)?;
        let staged = scan_staged_accessions(&staging_path)?;
        Ok(Self {
            fasta_path: fasta_path.to_path_buf(),
            staging_path,
            metadata_path: metadata_path_for(fasta_path),
            committed,
            staged,
            staging_writer: None,
        })
    }

    /// Path of the sidecar this store keeps consistent with its FASTA.
    pub fn metadata_path(&self) -> &Path {
        &self.metadata_path
    }

    /// Whether bases for `accession` are already on disk, committed or staged.
    ///
    /// **This is the presence question, and the FASTA is the only artifact that answers
    /// it.** The sidecar is never consulted: a row in it is an annotation claim, not a
    /// sequence.
    pub fn has_sequence(&self, accession: &str) -> bool {
        self.committed.contains_key(accession) || self.staged.contains(accession)
    }

    /// Number of records the committed FASTA holds beyond one per accession.
    pub fn duplicate_record_count(&self) -> usize {
        self.committed.values().map(|n| n.saturating_sub(1)).sum()
    }

    /// Accessions the committed FASTA holds more than one record for.
    pub fn duplicated_accessions(&self) -> Vec<&str> {
        self.committed
            .iter()
            .filter(|(_, n)| **n > 1)
            .map(|(acc, _)| acc.as_str())
            .collect()
    }

    /// Whether the on-disk artifacts are anything other than the fixed point.
    ///
    /// False means a commit would change no byte, so the caller must not perform one —
    /// that is what makes a second provisioning run a genuine no-op rather than a rewrite
    /// that happens to produce the same content.
    ///
    /// **Both derived artifacts are terms, not just the sidecar.** Declaring the `.fai`
    /// derived and then not noticing that it is missing would leave it derived in the
    /// documentation only — and a deleted index makes this artifact invisible to
    /// `prepare::load_existing_accessions`, which enumerates references by their `.fai`.
    pub fn is_dirty(&self, metadata: &SupplementalMetadata) -> bool {
        !self.staged.is_empty()
            || self.duplicate_record_count() > 0
            || self.metadata_keys_disagree(metadata)
            || self.index_disagrees()
    }

    /// Whether the sidecar's key set differs from the committed FASTA's accession set.
    fn metadata_keys_disagree(&self, metadata: &SupplementalMetadata) -> bool {
        if metadata.transcripts.len() != self.committed.len() {
            return true;
        }
        !self
            .committed
            .keys()
            .all(|acc| metadata.transcripts.contains_key(acc))
    }

    /// Whether the `.fai` is absent or names a different accession set than the FASTA.
    ///
    /// Compares names only. Detecting a *stale offset* would mean re-deriving the whole
    /// index, which costs what rebuilding it costs — so this catches a deleted, truncated
    /// or accession-drifted index and deliberately does not claim to catch an index whose
    /// names are right and whose offsets are wrong. Only this module writes the pair, and
    /// it writes them in one commit, so that state is not reachable through this code.
    fn index_disagrees(&self) -> bool {
        let fai_path = sibling_with_suffix(&self.fasta_path, ".fai");
        // An absent FASTA needs no index; an absent index for a populated FASTA is drift.
        if self.committed.is_empty() {
            return fai_path.exists();
        }
        let Ok(text) = std::fs::read_to_string(&fai_path) else {
            return true;
        };
        let indexed: BTreeMap<&str, usize> =
            text.lines()
                .filter(|l| !l.is_empty())
                .fold(BTreeMap::new(), |mut acc, line| {
                    *acc.entry(line.split('\t').next().unwrap_or(""))
                        .or_insert(0) += 1;
                    acc
                });
        indexed.len() != self.committed.len()
            || !self
                .committed
                .keys()
                .all(|acc| indexed.contains_key(acc.as_str()))
    }

    /// Append a fetched record to the staging file.
    ///
    /// Nothing is written to the FASTA itself until [`Self::commit`]. Staging an
    /// accession makes [`Self::has_sequence`] true for it immediately, so the same run
    /// cannot fetch it twice.
    pub fn stage_record(
        &mut self,
        accession: &str,
        description: &str,
        sequence: &str,
    ) -> Result<(), FerroError> {
        if sequence.is_empty() {
            return Ok(());
        }
        if self.staging_writer.is_none() {
            let file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(&self.staging_path)
                .map_err(|e| FerroError::Io {
                    msg: format!("Failed to open {}: {}", self.staging_path.display(), e),
                })?;
            self.staging_writer = Some(BufWriter::new(file));
        }

        let record = FastaRecord {
            accession: accession.to_string(),
            description: description.to_string(),
            sequence: sequence.to_string(),
        };
        let writer = self
            .staging_writer
            .as_mut()
            .expect("staging writer was just opened");
        let mut ignored_offset = 0u64;
        write_record(writer, &record, &mut ignored_offset)?;
        // Marks the record complete. A crash before this line lands leaves the record
        // unterminated, and the next run discards it rather than treating a truncated
        // sequence as present. See `StagedReader`.
        writeln!(writer, "{}", STAGED_RECORD_TERMINATOR).map_err(|e| FerroError::Io {
            msg: format!("Failed to write {}: {}", self.staging_path.display(), e),
        })?;
        self.staged.insert(record.accession);
        Ok(())
    }

    /// Merge staged records into the FASTA keyed by accession, then rebuild the derived
    /// artifacts.
    ///
    /// One record survives per accession: a staged record replaces the committed one for
    /// its accession, and every other committed record is carried through in its original
    /// order. "Carried through" is about content — accession, description and bases — not
    /// about bytes: every record is re-emitted at [`FASTA_LINE_BASES`], so a file wrapped
    /// at some other width is normalized. That is why the `.fai` is rewritten here and not
    /// left to a later `index_fasta`.
    ///
    /// `metadata` is reconciled against the result — orphan rows dropped, missing rows
    /// synthesized, existing CDS annotation preserved — and written atomically, as is the
    /// `.fai`.
    ///
    /// The caller must have checked [`Self::is_dirty`]; committing a clean store would
    /// rewrite three files to their existing content.
    pub fn commit(
        mut self,
        metadata: &mut SupplementalMetadata,
    ) -> Result<CommitReport, FerroError> {
        // Flush staged bytes to the OS before reading the staging file back.
        if let Some(writer) = self.staging_writer.take() {
            let file = writer.into_inner().map_err(|e| FerroError::Io {
                msg: format!("Failed to flush {}: {}", self.staging_path.display(), e),
            })?;
            file.sync_all().map_err(|e| FerroError::Io {
                msg: format!("Failed to sync {}: {}", self.staging_path.display(), e),
            })?;
        }

        let fasta_tmp = sibling_with_suffix(&self.fasta_path, TEMP_SUFFIX);
        let fai_path = sibling_with_suffix(&self.fasta_path, ".fai");
        let fai_tmp = sibling_with_suffix(&fai_path, TEMP_SUFFIX);

        let mut report = CommitReport::default();
        let mut offset = 0u64;
        // Accessions already written to the merged output, so the *second* record for an
        // accession is dropped whichever file it came from. This is the uniqueness
        // constraint the FASTA format does not supply.
        let mut emitted: HashSet<String> = HashSet::new();
        // Length and description of every surviving record, for sidecar reconciliation.
        let mut surviving: HashMap<String, SupplementalTranscript> = HashMap::new();

        {
            let out = File::create(&fasta_tmp).map_err(|e| FerroError::Io {
                msg: format!("Failed to create {}: {}", fasta_tmp.display(), e),
            })?;
            let mut writer = BufWriter::new(out);
            // Index lines are streamed alongside the records rather than accumulated:
            // the shipped artifact's index is 4.9 MB, and a merge must not hold a
            // multi-megabyte buffer it has somewhere to put.
            let index_out = File::create(&fai_tmp).map_err(|e| FerroError::Io {
                msg: format!("Failed to create {}: {}", fai_tmp.display(), e),
            })?;
            let mut index_writer = BufWriter::new(index_out);
            let index_err = |e: std::io::Error| FerroError::Io {
                msg: format!("Failed to write {}: {}", fai_tmp.display(), e),
            };

            // Committed records first, in their existing order, minus the ones a staged
            // record replaces and minus every duplicate.
            if let Some(mut reader) = open_fasta(&self.fasta_path)? {
                while let Some(record) = reader.next_record()? {
                    if self.staged.contains(&record.accession) {
                        report.records_replaced += 1;
                        continue;
                    }
                    if !emitted.insert(record.accession.clone()) {
                        report.duplicates_removed += 1;
                        continue;
                    }
                    surviving.insert(record.accession.clone(), record.synthesized_metadata());
                    let fai_line = write_record(&mut writer, &record, &mut offset)?;
                    index_writer
                        .write_all(fai_line.as_bytes())
                        .map_err(index_err)?;
                }
            }

            // Then the staged records, appended in the order they were fetched. A repeat
            // within one staging file keeps the first, which is unreachable in practice
            // because staging an accession makes `has_sequence` true for it.
            //
            // Read through `StagedReader`, so an incompletely-written trailing record is
            // dropped instead of being merged as a truncated sequence.
            if let Some(mut reader) = open_staged(&self.staging_path)? {
                while let Some(record) = reader.next_record()? {
                    if !emitted.insert(record.accession.clone()) {
                        report.duplicates_removed += 1;
                        continue;
                    }
                    surviving.insert(record.accession.clone(), record.synthesized_metadata());
                    let fai_line = write_record(&mut writer, &record, &mut offset)?;
                    index_writer
                        .write_all(fai_line.as_bytes())
                        .map_err(index_err)?;
                }
            }

            let file = writer.into_inner().map_err(|e| FerroError::Io {
                msg: format!("Failed to flush {}: {}", fasta_tmp.display(), e),
            })?;
            file.sync_all().map_err(|e| FerroError::Io {
                msg: format!("Failed to sync {}: {}", fasta_tmp.display(), e),
            })?;

            let index_file = index_writer.into_inner().map_err(|e| FerroError::Io {
                msg: format!("Failed to flush {}: {}", fai_tmp.display(), e),
            })?;
            index_file.sync_all().map_err(|e| FerroError::Io {
                msg: format!("Failed to sync {}: {}", fai_tmp.display(), e),
            })?;
        }

        // ---- Commit point. Everything below is derived and re-derivable. ----
        std::fs::rename(&fasta_tmp, &self.fasta_path).map_err(|e| FerroError::Io {
            msg: format!(
                "Failed to install {} over {}: {}",
                fasta_tmp.display(),
                self.fasta_path.display(),
                e
            ),
        })?;
        report.records = emitted.len();

        // The index describes the file that was just installed, so its own rename comes
        // after that one and never before.
        std::fs::rename(&fai_tmp, &fai_path).map_err(|e| FerroError::Io {
            msg: format!(
                "Failed to install {} over {}: {}",
                fai_tmp.display(),
                fai_path.display(),
                e
            ),
        })?;

        // Reconcile the sidecar against what the FASTA actually holds. A row with no
        // record is an orphan and goes; a record with no row gains one; a row that
        // already carries CDS annotation keeps it, because a FASTA cannot re-supply it.
        let before = metadata.transcripts.len();
        metadata
            .transcripts
            .retain(|acc, _| surviving.contains_key(acc));
        report.metadata_rows_dropped = before - metadata.transcripts.len();
        for (accession, synthesized) in surviving {
            metadata.transcripts.entry(accession).or_insert_with(|| {
                report.metadata_rows_synthesized += 1;
                synthesized
            });
        }

        let json = serde_json::to_vec_pretty(&metadata).map_err(|e| FerroError::Io {
            msg: format!(
                "Failed to serialize {}: {}",
                self.metadata_path.display(),
                e
            ),
        })?;
        let metadata_tmp = sibling_with_suffix(&self.metadata_path, TEMP_SUFFIX);
        write_atomically(&metadata_tmp, &self.metadata_path, &json)?;

        // Last, so a crash before this point simply replays a merge that is deterministic
        // and produces the identical FASTA.
        if self.staging_path.exists() {
            std::fs::remove_file(&self.staging_path).map_err(|e| FerroError::Io {
                msg: format!("Failed to remove {}: {}", self.staging_path.display(), e),
            })?;
        }

        Ok(report)
    }
}

/// Write `bytes` to `tmp`, `fsync` it, then rename it onto `final_path`.
fn write_atomically(tmp: &Path, final_path: &Path, bytes: &[u8]) -> Result<(), FerroError> {
    {
        let mut file = File::create(tmp).map_err(|e| FerroError::Io {
            msg: format!("Failed to create {}: {}", tmp.display(), e),
        })?;
        file.write_all(bytes).map_err(|e| FerroError::Io {
            msg: format!("Failed to write {}: {}", tmp.display(), e),
        })?;
        file.sync_all().map_err(|e| FerroError::Io {
            msg: format!("Failed to sync {}: {}", tmp.display(), e),
        })?;
    }
    std::fs::rename(tmp, final_path).map_err(|e| FerroError::Io {
        msg: format!(
            "Failed to install {} over {}: {}",
            tmp.display(),
            final_path.display(),
            e
        ),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::benchmark::cache::{fetch_fasta_to_file_with, EfetchSource};
    use std::cell::RefCell;
    use std::time::Duration;

    /// Build a sidecar with the given rows.
    fn metadata(rows: Vec<SupplementalTranscript>) -> SupplementalMetadata {
        SupplementalMetadata {
            generated_at: "2026-01-06T00:00:00Z".to_string(),
            transcripts: rows.into_iter().map(|t| (t.id.clone(), t)).collect(),
        }
    }

    /// A sidecar row carrying neither a CDS nor a length — the shape 116,572 of the
    /// shipped artifact's 140,027 rows have.
    fn hollow_row(accession: &str) -> SupplementalTranscript {
        SupplementalTranscript {
            id: accession.to_string(),
            gene_symbol: None,
            cds_start: None,
            cds_end: None,
            sequence_length: 0,
        }
    }

    fn annotated_row(accession: &str, cds: (u64, u64), length: usize) -> SupplementalTranscript {
        SupplementalTranscript {
            id: accession.to_string(),
            gene_symbol: Some("GENEA".to_string()),
            cds_start: Some(cds.0),
            cds_end: Some(cds.1),
            sequence_length: length,
        }
    }

    /// Write a FASTA holding exactly the given `(accession, description, sequence)` records,
    /// in order, with duplicates preserved.
    fn write_fasta(path: &Path, records: &[(&str, &str, &str)]) {
        let mut body = String::new();
        for (accession, description, sequence) in records {
            if description.is_empty() {
                body.push_str(&format!(">{}\n", accession));
            } else {
                body.push_str(&format!(">{} {}\n", accession, description));
            }
            for chunk in sequence.as_bytes().chunks(FASTA_LINE_BASES) {
                body.push_str(std::str::from_utf8(chunk).unwrap());
                body.push('\n');
            }
        }
        std::fs::write(path, body).unwrap();
    }

    /// Write a staging file: complete records are terminated, and `torn` optionally
    /// appends an unterminated trailing record — what a crash mid-write leaves behind.
    fn write_staged(path: &Path, records: &[(&str, &str, &str)], torn: Option<(&str, &str)>) {
        let mut body = String::new();
        for (accession, description, sequence) in records {
            body.push_str(&format!(">{} {}\n", accession, description));
            for chunk in sequence.as_bytes().chunks(FASTA_LINE_BASES) {
                body.push_str(std::str::from_utf8(chunk).unwrap());
                body.push('\n');
            }
            body.push_str(STAGED_RECORD_TERMINATOR);
            body.push('\n');
        }
        if let Some((accession, partial)) = torn {
            body.push_str(&format!(">{}\n{}\n", accession, partial));
        }
        std::fs::write(path, body).unwrap();
    }

    /// Accessions of every record in a FASTA, in file order, duplicates included.
    fn headers(path: &Path) -> Vec<String> {
        let text = std::fs::read_to_string(path).unwrap();
        text.lines()
            .filter_map(|l| l.strip_prefix('>'))
            .map(|h| split_header(h).0)
            .collect()
    }

    /// Sequence recorded for `accession`, read back through the record stream.
    fn sequence_of(path: &Path, accession: &str) -> Option<String> {
        let mut reader = open_fasta(path).unwrap()?;
        while let Some(record) = reader.next_record().unwrap() {
            if record.accession == accession {
                return Some(record.sequence);
            }
        }
        None
    }

    // ---------------------------------------------------------------- store layer

    #[test]
    fn presence_is_answered_by_the_fasta_and_never_by_the_sidecar() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        write_fasta(&fasta, &[("NM_000001.1", "GENEA", "ACGTACGT")]);

        let store = SupplementalStore::open(&fasta).unwrap();

        // Bases on disk, sidecar row hollow -> still present. This is the accession the
        // issue's suggested acceptance test is about.
        assert!(store.has_sequence("NM_000001.1"));
        // A sidecar row with no FASTA record is not presence: the sidecar is not asked.
        assert!(!store.has_sequence("NM_000002.1"));
    }

    #[test]
    fn a_staged_record_replaces_the_committed_one_rather_than_appending() {
        // This is the property that makes PR #1732's usability-keyed resume predicate
        // safe: even when an accession whose bases are already on disk is re-fetched, the
        // write path upserts rather than appends.
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        write_fasta(
            &fasta,
            &[
                ("NM_000001.1", "GENEA", "AAAA"),
                ("NM_000002.1", "", "CCCC"),
            ],
        );
        let mut meta = metadata(vec![hollow_row("NM_000001.1"), hollow_row("NM_000002.1")]);

        let mut store = SupplementalStore::open(&fasta).unwrap();
        store
            .stage_record("NM_000001.1", "GENEA", "GGGGGG")
            .unwrap();
        let report = store.commit(&mut meta).unwrap();

        assert_eq!(headers(&fasta), vec!["NM_000002.1", "NM_000001.1"]);
        assert_eq!(sequence_of(&fasta, "NM_000001.1").unwrap(), "GGGGGG");
        assert_eq!(sequence_of(&fasta, "NM_000002.1").unwrap(), "CCCC");
        assert_eq!(report.records, 2);
        assert_eq!(report.records_replaced, 1);
        assert_eq!(report.duplicates_removed, 0);
    }

    #[test]
    fn a_duplicated_accession_is_collapsed_to_one_record() {
        // The shipped `patterns_transcripts.fna.bak` is in this state: 159,456 records for
        // 140,026 unique accessions.
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        write_fasta(
            &fasta,
            &[
                ("NM_000001.1", "GENEA", "AAAA"),
                ("NM_000002.1", "", "CCCC"),
                ("NM_000001.1", "GENEA", "AAAA"),
                ("NM_000001.1", "GENEA", "AAAA"),
            ],
        );
        let mut meta = metadata(vec![hollow_row("NM_000001.1"), hollow_row("NM_000002.1")]);

        let store = SupplementalStore::open(&fasta).unwrap();
        assert_eq!(store.duplicate_record_count(), 2);
        assert_eq!(store.duplicated_accessions(), vec!["NM_000001.1"]);
        assert!(store.is_dirty(&meta));

        let report = store.commit(&mut meta).unwrap();

        assert_eq!(headers(&fasta), vec!["NM_000001.1", "NM_000002.1"]);
        assert_eq!(report.duplicates_removed, 2);
        assert_eq!(report.records, 2);
        // The first occurrence is kept, so a repair preserves the record that was already
        // being served rather than picking an arbitrary copy.
        assert_eq!(sequence_of(&fasta, "NM_000001.1").unwrap(), "AAAA");
    }

    #[test]
    fn an_interrupted_runs_staging_file_counts_as_already_fetched() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        write_fasta(&fasta, &[("NM_000001.1", "GENEA", "AAAA")]);
        // What an interrupted run leaves behind: one complete staged record.
        write_staged(
            &sibling_with_suffix(&fasta, STAGING_SUFFIX),
            &[("NM_000002.1", "GENEB", "CCCC")],
            None,
        );

        let store = SupplementalStore::open(&fasta).unwrap();
        assert!(store.has_sequence("NM_000001.1"));
        assert!(
            store.has_sequence("NM_000002.1"),
            "staged bases are on disk"
        );

        let mut meta = metadata(vec![hollow_row("NM_000001.1")]);
        assert!(store.is_dirty(&meta));
        store.commit(&mut meta).unwrap();

        assert_eq!(headers(&fasta), vec!["NM_000001.1", "NM_000002.1"]);
        assert!(!sibling_with_suffix(&fasta, STAGING_SUFFIX).exists());
    }

    #[test]
    fn the_sidecar_key_set_is_reconciled_against_the_fasta() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        write_fasta(&fasta, &[("NM_000001.1", "GENEA (variant 1)", "AAAACCCC")]);
        // One orphan row (no FASTA record) and no row at all for the record that exists.
        let mut meta = metadata(vec![hollow_row("NM_999999.9")]);

        let store = SupplementalStore::open(&fasta).unwrap();
        assert!(store.is_dirty(&meta), "key sets disagree");
        let report = store.commit(&mut meta).unwrap();

        assert_eq!(report.metadata_rows_dropped, 1);
        assert_eq!(report.metadata_rows_synthesized, 1);
        assert_eq!(meta.transcripts.len(), 1);
        let row = &meta.transcripts["NM_000001.1"];
        assert_eq!(row.sequence_length, 8);
        assert_eq!(row.gene_symbol.as_deref(), Some("GENEA"));
        assert_eq!(row.cds_start, None, "a FASTA cannot supply a CDS");
    }

    #[test]
    fn existing_cds_annotation_is_carried_forward_not_overwritten() {
        // The sidecar's key set is derived; its CDS payload is not, and a commit must not
        // reset an annotated row to the FASTA-only shape.
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        write_fasta(&fasta, &[("NM_000001.1", "GENEA", "AAAACCCC")]);
        let mut meta = metadata(vec![
            annotated_row("NM_000001.1", (2, 7), 8),
            hollow_row("NM_999999.9"),
        ]);

        SupplementalStore::open(&fasta)
            .unwrap()
            .commit(&mut meta)
            .unwrap();

        let row = &meta.transcripts["NM_000001.1"];
        assert_eq!(row.cds_start, Some(2));
        assert_eq!(row.cds_end, Some(7));
    }

    #[test]
    fn a_reconciled_artifact_is_not_dirty() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        write_fasta(&fasta, &[("NM_000001.1", "GENEA", "AAAA")]);
        // A hollow row is still a row, so the key sets agree and this predicate is
        // satisfied. Whether that row's DERIVED `sequence_length` is stale is PR #1732's
        // question, and it is deliberately not a term here: re-deriving it from the FASTA
        // rewrites only the sidecar, so it must not make a 1.3 GB artifact dirty.
        let meta = metadata(vec![hollow_row("NM_000001.1")]);
        // A complete artifact carries its index; without one the store is dirty, which
        // `a_missing_or_drifted_index_makes_the_store_dirty` pins separately.
        std::fs::write(
            sibling_with_suffix(&fasta, ".fai"),
            "NM_000001.1\t4\t19\t4\t5\n",
        )
        .unwrap();

        assert!(!SupplementalStore::open(&fasta).unwrap().is_dirty(&meta));
    }

    #[test]
    fn the_rebuilt_index_locates_every_record_exactly() {
        use std::io::{Read, Seek, SeekFrom};

        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        // One record shorter than a line, one spanning three lines, so both the
        // `line_bases` branch and multi-line offset arithmetic are exercised.
        let long: String = std::iter::repeat_n('G', 155).collect();
        write_fasta(&fasta, &[("NM_000001.1", "GENEA", "ACGT")]);
        let mut store = SupplementalStore::open(&fasta).unwrap();
        store.stage_record("NM_000002.1", "", &long).unwrap();
        let mut meta = metadata(vec![]);
        store.commit(&mut meta).unwrap();

        let fai = std::fs::read_to_string(sibling_with_suffix(&fasta, ".fai")).unwrap();
        let mut file = File::open(&fasta).unwrap();
        let mut seen = 0;
        for line in fai.lines() {
            let fields: Vec<&str> = line.split('\t').collect();
            let (length, offset): (usize, u64) =
                (fields[1].parse().unwrap(), fields[2].parse().unwrap());
            let line_bases: usize = fields[3].parse().unwrap();
            let line_bytes: usize = fields[4].parse().unwrap();
            assert_eq!(line_bytes, line_bases + 1);
            assert_eq!(line_bases, length.min(FASTA_LINE_BASES));

            // Seek to the recorded offset and read the record's bases back. `take` +
            // `read_to_end` rather than a bare `read`, which may return a short buffer
            // and produce a spurious failure rather than a real one.
            file.seek(SeekFrom::Start(offset)).unwrap();
            let mut raw = Vec::new();
            Read::by_ref(&mut file)
                .take((length + length / FASTA_LINE_BASES + 1) as u64)
                .read_to_end(&mut raw)
                .unwrap();
            let bases: String = String::from_utf8_lossy(&raw)
                .chars()
                .filter(|c| c.is_ascii_alphabetic())
                .take(length)
                .collect();
            let expected = sequence_of(&fasta, fields[0]).unwrap();
            assert_eq!(bases, expected, "index misplaces {}", fields[0]);
            seen += 1;
        }
        assert_eq!(seen, 2);
    }

    #[test]
    fn a_header_with_no_sequence_is_not_a_presence_claim() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        std::fs::write(&fasta, ">NM_000001.1 GENEA\n>NM_000002.1\nACGT\n").unwrap();

        let store = SupplementalStore::open(&fasta).unwrap();
        assert!(!store.has_sequence("NM_000001.1"), "no bases behind it");
        assert!(store.has_sequence("NM_000002.1"));
    }

    #[test]
    fn a_temp_path_extends_the_name_rather_than_replacing_its_extension() {
        // `with_extension` would turn `supp.fna` into `supp.tmp`, which collides with the
        // temp path of every sibling artifact in the directory.
        let path = Path::new("/x/patterns_transcripts.fna");
        assert_eq!(
            sibling_with_suffix(path, TEMP_SUFFIX),
            Path::new("/x/patterns_transcripts.fna.tmp")
        );
        assert_eq!(
            metadata_path_for(path),
            Path::new("/x/patterns_transcripts.metadata.json")
        );
    }

    // ------------------------------------------------- end-to-end provisioning path

    /// An [`EfetchSource`] serving canned payloads and counting what it was asked for.
    struct FakeEfetch {
        /// `accession -> (gene, sequence)`.
        records: HashMap<String, (String, String)>,
        genbank_calls: RefCell<Vec<String>>,
        fasta_calls: RefCell<Vec<String>>,
    }

    impl FakeEfetch {
        fn new(records: &[(&str, &str, &str)]) -> Self {
            Self {
                records: records
                    .iter()
                    .map(|(a, g, s)| (a.to_string(), (g.to_string(), s.to_string())))
                    .collect(),
                genbank_calls: RefCell::new(Vec::new()),
                fasta_calls: RefCell::new(Vec::new()),
            }
        }

        fn calls(&self) -> usize {
            self.genbank_calls.borrow().len() + self.fasta_calls.borrow().len()
        }
    }

    impl EfetchSource for FakeEfetch {
        fn genbank(&self, ids: &str) -> Result<Option<String>, FerroError> {
            self.genbank_calls.borrow_mut().push(ids.to_string());
            let mut text = String::new();
            for id in ids.split(',') {
                let Some((gene, sequence)) = self.records.get(id) else {
                    continue;
                };
                text.push_str(&format!(
                    "LOCUS       {id}    {len} bp    mRNA\n\
                     VERSION     {id}\n\
                     FEATURES             Location/Qualifiers\n     \
                     gene            1..{len}\n                     \
                     /gene=\"{gene}\"\n     \
                     CDS             1..{len}\n\
                     ORIGIN\n        1 {sequence}\n//\n",
                    len = sequence.len()
                ));
            }
            Ok((!text.is_empty()).then_some(text))
        }

        fn fasta(&self, ids: &str) -> Option<String> {
            self.fasta_calls.borrow_mut().push(ids.to_string());
            None
        }

        fn rate_delay(&self) -> Duration {
            Duration::ZERO
        }
    }

    /// Content **and modification time** of the three artifacts.
    ///
    /// The mtime is not redundant with the bytes: the merge is deterministic, so a run
    /// that needlessly rewrote all three would produce byte-identical content and satisfy
    /// a content-only comparison. "The second run is a no-op" means it performed no write,
    /// which only the mtime can witness.
    fn artifact_state(fasta: &Path) -> Vec<(String, Vec<u8>, std::time::SystemTime)> {
        [
            fasta.to_path_buf(),
            sibling_with_suffix(fasta, ".fai"),
            metadata_path_for(fasta),
        ]
        .iter()
        .map(|p| {
            let bytes = std::fs::read(p).unwrap_or_else(|e| panic!("read {}: {}", p.display(), e));
            let mtime = std::fs::metadata(p).unwrap().modified().unwrap();
            (p.display().to_string(), bytes, mtime)
        })
        .collect()
    }

    #[test]
    fn provisioning_is_idempotent_and_the_second_run_writes_nothing() {
        // The acceptance property. It is deliberately not "a hollow row does not
        // double-append": stated as a fixed point it holds every write path to account
        // for itself, including ones nobody anticipated.
        //
        // NOTE, because this test reads like the guard and is not the whole of it: seeded
        // from an EMPTY directory it also passes on the pre-fix behaviour, since an
        // artifact that behaviour produced has a FASTA and a sidecar that agree, so its
        // sidecar-keyed resume filter skips correctly. The tests that are red on the old
        // behaviour are the ones seeded from a state where the two disagree —
        // `an_accession_in_the_fasta_but_missing_from_the_sidecar_is_not_re_fetched` and
        // `a_duplicate_bearing_artifact_is_repaired_and_then_stable`. Do not delete those
        // as redundant with this one.
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        let source = FakeEfetch::new(&[
            ("NM_000001.1", "GENEA", "ACGTACGTAC"),
            ("NM_000002.1", "GENEB", "TTTTGGGGCCCC"),
        ]);
        let requested = vec!["NM_000001.1".to_string(), "NM_000002.1".to_string()];

        let first = fetch_fasta_to_file_with(&requested, &fasta, false, &source).unwrap();
        assert_eq!(first, 2);
        let after_first = artifact_state(&fasta);
        let calls_after_first = source.calls();
        assert!(calls_after_first > 0, "the first run must actually fetch");

        let second = fetch_fasta_to_file_with(&requested, &fasta, false, &source).unwrap();

        assert_eq!(second, 0);
        assert_eq!(
            source.calls(),
            calls_after_first,
            "the second run must not fetch anything"
        );
        assert_eq!(
            after_first,
            artifact_state(&fasta),
            "the second run must leave every artifact byte-identical"
        );
        assert_eq!(headers(&fasta), vec!["NM_000001.1", "NM_000002.1"]);
        assert!(
            !sibling_with_suffix(&fasta, STAGING_SUFFIX).exists(),
            "no staging file survives a clean run"
        );
    }

    #[test]
    fn an_accession_in_the_fasta_but_missing_from_the_sidecar_is_not_re_fetched() {
        // The route that produced the shipped `.bak`: the sidecar is written once, after
        // the final batch, so an interrupted run leaves the FASTA ahead of it. Keyed on
        // the sidecar this re-fetches and appends a second record; keyed on the FASTA it
        // is already present.
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        write_fasta(&fasta, &[("NM_000001.1", "GENEA", "ACGTACGTAC")]);
        // Sidecar knows nothing about it.
        std::fs::write(
            metadata_path_for(&fasta),
            serde_json::to_vec_pretty(&metadata(vec![])).unwrap(),
        )
        .unwrap();

        let source = FakeEfetch::new(&[("NM_000001.1", "GENEA", "ACGTACGTAC")]);
        let fetched =
            fetch_fasta_to_file_with(&["NM_000001.1".to_string()], &fasta, false, &source).unwrap();

        assert_eq!(fetched, 0);
        assert_eq!(source.calls(), 0, "nothing was fetched");
        assert_eq!(
            headers(&fasta),
            vec!["NM_000001.1"],
            "the FASTA gained no second record"
        );
        // The run still repaired the sidecar, because its key set had drifted.
        let repaired: SupplementalMetadata =
            serde_json::from_slice(&std::fs::read(metadata_path_for(&fasta)).unwrap()).unwrap();
        assert_eq!(repaired.transcripts.len(), 1);
        assert_eq!(repaired.transcripts["NM_000001.1"].sequence_length, 10);
    }

    #[test]
    fn a_hollow_sidecar_row_does_not_produce_a_second_fasta_record() {
        // The issue's own suggested acceptance test: run the resume filter over a cohort
        // containing a hollow row whose sequence is already in the FASTA, and assert the
        // FASTA gains no second header for it.
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        write_fasta(&fasta, &[("NM_000001.1", "GENEA", "ACGTACGTAC")]);
        std::fs::write(
            metadata_path_for(&fasta),
            serde_json::to_vec_pretty(&metadata(vec![hollow_row("NM_000001.1")])).unwrap(),
        )
        .unwrap();
        // A complete artifact has an index too; it is part of what must not move.
        std::fs::write(
            sibling_with_suffix(&fasta, ".fai"),
            "NM_000001.1\t10\t19\t10\t11\n",
        )
        .unwrap();

        let source = FakeEfetch::new(&[("NM_000001.1", "GENEA", "ACGTACGTAC")]);
        let before = artifact_state(&fasta);

        fetch_fasta_to_file_with(&["NM_000001.1".to_string()], &fasta, false, &source).unwrap();

        assert_eq!(headers(&fasta), vec!["NM_000001.1"]);
        assert_eq!(source.calls(), 0);
        // Key sets already agreed, so the FASTA and its index were already the fixed point
        // and neither is rewritten. The sidecar is: #1732 re-derives the hollow row's
        // `sequence_length` from the FASTA, which is a **derived** field and therefore
        // repairable with no source call — which is why `source.calls()` above is still 0.
        // Pinned by `a_hollow_sidecar_row_is_repaired_from_the_fasta_without_a_fetch`.
        assert_eq!(
            before[..2],
            artifact_state(&fasta)[..2],
            "repairing a derived field must not rewrite the authoritative FASTA or its index"
        );
    }

    #[test]
    fn a_hollow_sidecar_row_is_repaired_from_the_fasta_without_a_fetch() {
        // #1722's producer half, and the ruling it was fixed under: **never fetch to
        // recover a value you can compute**. `sequence_length` is derived — every site
        // that builds a row sets it to the length of the sequence it just wrote — so a
        // hollow row whose accession is already in the FASTA is a stale derived field, not
        // a missing fetch. All 116,572 hollow rows of the shipped artifact are in its
        // FASTA, which is why #1791's presence-keyed resume filter admits none of them and
        // why the repair has to be local.
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        write_fasta(&fasta, &[("NM_000001.1", "GENEA", "ACGTACGTAC")]);
        std::fs::write(
            metadata_path_for(&fasta),
            serde_json::to_vec_pretty(&metadata(vec![hollow_row("NM_000001.1")])).unwrap(),
        )
        .unwrap();
        std::fs::write(
            sibling_with_suffix(&fasta, ".fai"),
            "NM_000001.1\t10\t19\t10\t11\n",
        )
        .unwrap();

        let source = FakeEfetch::new(&[("NM_000001.1", "GENEA", "ACGTACGTAC")]);
        fetch_fasta_to_file_with(&["NM_000001.1".to_string()], &fasta, false, &source).unwrap();

        assert_eq!(
            source.calls(),
            0,
            "a derived field is computed from the FASTA, never fetched"
        );
        let repaired: SupplementalMetadata =
            serde_json::from_slice(&std::fs::read(metadata_path_for(&fasta)).unwrap()).unwrap();
        let row = &repaired.transcripts["NM_000001.1"];
        assert_eq!(
            row.sequence_length, 10,
            "the hollow row is no longer hollow"
        );
        assert_eq!(
            row.cds_start, None,
            "a CDS is sourced from GenBank and stays absent"
        );
    }

    #[test]
    fn a_repaired_sidecar_is_the_fixed_point_on_the_next_run() {
        // The repair must converge, or every provisioning run rewrites a 24 MB sidecar for
        // nothing. Seeded from the damaged state, as this module's fixed-point tests are.
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        write_fasta(&fasta, &[("NM_000001.1", "GENEA", "ACGTACGTAC")]);
        std::fs::write(
            metadata_path_for(&fasta),
            serde_json::to_vec_pretty(&metadata(vec![hollow_row("NM_000001.1")])).unwrap(),
        )
        .unwrap();
        std::fs::write(
            sibling_with_suffix(&fasta, ".fai"),
            "NM_000001.1\t10\t19\t10\t11\n",
        )
        .unwrap();
        let source = FakeEfetch::new(&[("NM_000001.1", "GENEA", "ACGTACGTAC")]);
        let requested = vec!["NM_000001.1".to_string()];

        fetch_fasta_to_file_with(&requested, &fasta, false, &source).unwrap();
        let repaired = artifact_state(&fasta);

        fetch_fasta_to_file_with(&requested, &fasta, false, &source).unwrap();

        assert_eq!(
            repaired,
            artifact_state(&fasta),
            "a repaired artifact is the fixed point; the second run writes no byte"
        );
        assert_eq!(source.calls(), 0);
    }

    #[test]
    fn a_repair_survives_the_commit_that_runs_alongside_it() {
        // The two repairs compose or they do not, and `commit` fills only *missing* rows
        // (`or_insert_with`), so a row it already has is carried through as-is. That is
        // what makes the ordering load-bearing: the derived lengths are re-derived BEFORE
        // the commit, so the commit carries them out rather than preserving the hollow
        // shape. Seeded from an artifact that is damaged in both ways at once — duplicate
        // records AND hollow rows — because that is the only state where the interaction
        // is observable.
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        write_fasta(
            &fasta,
            &[
                ("NM_000001.1", "GENEA", "ACGTACGTAC"),
                ("NM_000002.1", "GENEB", "TTTT"),
                ("NM_000001.1", "GENEA", "ACGTACGTAC"),
            ],
        );
        std::fs::write(
            metadata_path_for(&fasta),
            serde_json::to_vec_pretty(&metadata(vec![
                hollow_row("NM_000001.1"),
                hollow_row("NM_000002.1"),
            ]))
            .unwrap(),
        )
        .unwrap();

        let source = FakeEfetch::new(&[]);
        fetch_fasta_to_file_with(
            &["NM_000001.1".to_string(), "NM_000002.1".to_string()],
            &fasta,
            false,
            &source,
        )
        .unwrap();

        let committed: SupplementalMetadata =
            serde_json::from_slice(&std::fs::read(metadata_path_for(&fasta)).unwrap()).unwrap();
        assert_eq!(committed.transcripts["NM_000001.1"].sequence_length, 10);
        assert_eq!(committed.transcripts["NM_000002.1"].sequence_length, 4);
        assert_eq!(source.calls(), 0);
    }

    #[test]
    fn a_duplicate_bearing_artifact_is_repaired_and_then_stable() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        write_fasta(
            &fasta,
            &[
                ("NM_000001.1", "GENEA", "ACGTACGTAC"),
                ("NM_000002.1", "GENEB", "TTTT"),
                ("NM_000001.1", "GENEA", "ACGTACGTAC"),
            ],
        );
        std::fs::write(
            metadata_path_for(&fasta),
            serde_json::to_vec_pretty(&metadata(vec![
                hollow_row("NM_000001.1"),
                hollow_row("NM_000002.1"),
            ]))
            .unwrap(),
        )
        .unwrap();

        let source = FakeEfetch::new(&[]);
        let requested = vec!["NM_000001.1".to_string(), "NM_000002.1".to_string()];

        fetch_fasta_to_file_with(&requested, &fasta, false, &source).unwrap();
        assert_eq!(headers(&fasta), vec!["NM_000001.1", "NM_000002.1"]);
        let repaired = artifact_state(&fasta);

        // Repair reaches a fixed point rather than churning.
        fetch_fasta_to_file_with(&requested, &fasta, false, &source).unwrap();
        assert_eq!(repaired, artifact_state(&fasta));
        assert_eq!(source.calls(), 0);
    }

    #[test]
    fn an_incompletely_written_staged_record_is_discarded_not_merged() {
        // A crash part-way through `stage_record` leaves a truncated trailing record. If
        // it counted as present, the accession would never be re-fetched and the
        // TRUNCATED sequence would be merged into the authoritative FASTA — a torn write
        // silently becoming the authority, which is the defect class this module closes.
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        write_staged(
            &sibling_with_suffix(&fasta, STAGING_SUFFIX),
            &[("NM_000001.1", "GENEA", "ACGTACGT")],
            Some(("NM_000002.1", "TTTT")), // header written, terminator never reached
        );

        let store = SupplementalStore::open(&fasta).unwrap();
        assert!(store.has_sequence("NM_000001.1"), "complete record counts");
        assert!(
            !store.has_sequence("NM_000002.1"),
            "a torn record is not a presence claim, so it is re-fetched"
        );

        let mut meta = metadata(vec![]);
        store.commit(&mut meta).unwrap();
        assert_eq!(headers(&fasta), vec!["NM_000001.1"]);
        assert_eq!(
            sequence_of(&fasta, "NM_000002.1"),
            None,
            "the truncated sequence must not reach the committed FASTA"
        );
    }

    #[test]
    fn a_missing_or_drifted_index_makes_the_store_dirty() {
        // The `.fai` is declared derived, so noticing it is gone is part of that claim
        // being true rather than documentary. A deleted index also makes this artifact
        // invisible to `prepare::load_existing_accessions`, which enumerates by `.fai`.
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        write_fasta(&fasta, &[("NM_000001.1", "GENEA", "ACGT")]);
        let meta = metadata(vec![hollow_row("NM_000001.1")]);

        // No index at all.
        assert!(SupplementalStore::open(&fasta).unwrap().is_dirty(&meta));

        // An index naming a different accession set.
        std::fs::write(
            sibling_with_suffix(&fasta, ".fai"),
            "NM_999999.9\t4\t19\t4\t5\n",
        )
        .unwrap();
        assert!(SupplementalStore::open(&fasta).unwrap().is_dirty(&meta));

        // A commit reconciles it, and the result is the fixed point.
        let mut meta2 = meta.clone();
        SupplementalStore::open(&fasta)
            .unwrap()
            .commit(&mut meta2)
            .unwrap();
        assert!(!SupplementalStore::open(&fasta).unwrap().is_dirty(&meta2));
    }

    #[test]
    fn a_run_whose_every_fetch_failed_rewrites_nothing() {
        // Reaching the commit having staged nothing is a real outcome — a down or
        // rate-limiting endpoint fails every batch — and a failed run must not rewrite
        // the artifact it could not add to.
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        write_fasta(&fasta, &[("NM_000001.1", "GENEA", "ACGTACGTAC")]);
        // Seeded with a HEALTHY row, whose length already agrees with the FASTA, so the
        // only thing this test can observe is the failed fetch. A hollow row here would
        // also make #1732's derived-length repair fire, and the sidecar would move for a
        // reason that has nothing to do with the property being pinned — that repair is
        // orthogonal and has its own tests.
        std::fs::write(
            metadata_path_for(&fasta),
            serde_json::to_vec_pretty(&metadata(vec![annotated_row("NM_000001.1", (1, 9), 10)]))
                .unwrap(),
        )
        .unwrap();
        std::fs::write(
            sibling_with_suffix(&fasta, ".fai"),
            "NM_000001.1\t10\t19\t10\t11\n",
        )
        .unwrap();
        let before = artifact_state(&fasta);

        // `NM_000002.1` is requested and genuinely absent, so the fetch loop runs — and
        // the fake has no record for it, so every batch comes back empty.
        let source = FakeEfetch::new(&[]);
        let fetched = fetch_fasta_to_file_with(
            &["NM_000001.1".to_string(), "NM_000002.1".to_string()],
            &fasta,
            false,
            &source,
        )
        .unwrap();

        assert_eq!(fetched, 0);
        assert!(source.calls() > 0, "the fetch loop must actually have run");
        assert_eq!(
            before,
            artifact_state(&fasta),
            "a run that staged nothing must not rewrite the artifact"
        );
    }

    #[test]
    fn a_dry_run_writes_nothing() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("supp.fna");
        let source = FakeEfetch::new(&[("NM_000001.1", "GENEA", "ACGT")]);

        fetch_fasta_to_file_with(&["NM_000001.1".to_string()], &fasta, true, &source).unwrap();

        assert!(!fasta.exists());
        assert!(!metadata_path_for(&fasta).exists());
        assert_eq!(source.calls(), 0);
    }
}
