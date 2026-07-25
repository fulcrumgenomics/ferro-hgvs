//! Construction of the FASTA index a [`MultiFastaProvider`] serves from.
//!
//! This module owns index construction and nothing else. Every provider
//! constructor will route through [`PreparedIndex::from_dirs`] (the routing
//! lands in the next commit), so there will be exactly one producer of this
//! state — the property a future on-disk index artifact depends on, since an
//! artifact loader would otherwise be a second, silently divergent producer.
//!
//! [`MultiFastaProvider`]: crate::reference::multi_fasta::MultiFastaProvider

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use rustc_hash::FxHashMap;

use crate::error::FerroError;

/// Index entry for a sequence in a FASTA file.
#[derive(Debug, Clone, Copy)]
pub(crate) struct FastaIndexEntry {
    /// Index into the file table (see [`PreparedIndex::files`]) this entry's
    /// index was built against.
    pub(crate) file_id: u32,
    /// Length of the sequence.
    pub(crate) length: u64,
    /// Byte offset to the start of sequence data.
    pub(crate) offset: u64,
    /// Number of bases per line.
    pub(crate) line_bases: u64,
    /// Number of bytes per line (including newline).
    pub(crate) line_bytes: u64,
}

/// Parse the integer version suffix of an accession (`"NM_003002.4"` -> `Some(4)`),
/// or `None` when there is no trailing numeric `.<n>`.
fn accession_version(accession: &str) -> Option<u32> {
    accession
        .rsplit_once('.')
        .and_then(|(_, v)| v.parse::<u32>().ok())
}

/// Build the `base -> versioned` fallback map from a fully populated FASTA
/// index, deterministically choosing the **highest** version per base
/// accession.
///
/// Building this by iterating the index in hash order with last-writer-wins
/// made version fallback nondeterministic across runs (the index uses a
/// fixed-seed hasher now, but relying on iteration order for *which* version
/// wins is fragile regardless). Keying on `(version, accession)` gives a stable
/// total order: the newest version always wins, independent of insertion or
/// iteration order.
pub(super) fn build_base_to_versioned(
    index: &FxHashMap<String, FastaIndexEntry>,
) -> FxHashMap<String, String> {
    let mut map: FxHashMap<String, String> = FxHashMap::default();
    for name in index.keys() {
        let Some(base) = name.split('.').next() else {
            continue;
        };
        let replace = match map.get(base) {
            None => true,
            Some(existing) => {
                (accession_version(name), name.as_str())
                    > (accession_version(existing), existing.as_str())
            }
        };
        if replace {
            map.insert(base.to_string(), name.clone());
        }
    }
    map
}

/// Whether `index` holds any genomic (chromosome) sequence, named either
/// NCBI RefSeq style (`NC_*`) or UCSC style (`chr*`). Computed once at
/// construction and cached in `MultiFastaProvider::has_genomic_data` — the
/// index can hold hundreds of thousands of keys and this predicate is hit on
/// the per-variant normalization path.
pub(super) fn index_has_genomic_data(index: &FxHashMap<String, FastaIndexEntry>) -> bool {
    index
        .keys()
        .any(|k| k.starts_with("NC_") || k.starts_with("chr"))
}

/// The FASTA index and everything derived from it.
// Temporary: nothing constructs this until the provider is routed through it
// in the next commit; this allow is removed there.
#[allow(dead_code)]
#[derive(Debug, Default)]
pub(crate) struct PreparedIndex {
    /// FASTA files backing the index, addressed by [`FastaIndexEntry::file_id`].
    ///
    /// APPEND-ONLY: `file_id` is a positional index into this `Vec`, so an
    /// entry must never be reordered or removed once assigned. A future
    /// `sort()` or `retain()` here would silently repoint every existing
    /// `FastaIndexEntry` at the wrong FASTA — a correctness bug with no
    /// compiler signal and, absent a targeted test, no test signal either.
    pub files: Vec<PathBuf>,
    /// Sequence accession -> location.
    pub index: FxHashMap<String, FastaIndexEntry>,
    /// Base accession (no version) -> highest-versioned accession present.
    pub base_to_versioned: FxHashMap<String, String>,
    /// Whether the index holds any genomic (chromosome) sequence.
    pub has_genomic_data: bool,
    /// Protein FASTA index (NP_/XP_), addressed into the same `files` table.
    pub protein_index: FxHashMap<String, FastaIndexEntry>,
}

// Temporary: nothing constructs a `PreparedIndex` until the provider is
// routed through it in the next commit; this allow is removed there.
#[allow(dead_code)]
impl PreparedIndex {
    /// Build from `dirs`, in order. A later directory wins for an accession
    /// present in more than one.
    pub(crate) fn from_dirs(dirs: &[&Path]) -> Result<Self, FerroError> {
        let (index, files) = load_index_from_dirs(dirs)?;
        Ok(Self {
            base_to_versioned: build_base_to_versioned(&index),
            has_genomic_data: index_has_genomic_data(&index),
            index,
            files,
            protein_index: FxHashMap::default(),
        })
    }

    /// Index a protein FASTA (`NP_`/`XP_` accessions) into the same file-id
    /// space as the rest of the reference, appending its entries to
    /// [`Self::protein_index`]. Returns the number of entries parsed from
    /// `fai` — NOT net growth in `protein_index`: `extend` overwrites an
    /// accession already present (e.g. from an earlier protein FASTA), so a
    /// caller that sums this return value across multiple calls can overcount
    /// relative to `protein_index.len()`.
    ///
    /// Mirrors the file-table bookkeeping [`load_index_from_dirs`] does for
    /// directory-scanned FASTAs: the canonical path is looked up in the
    /// existing `files` table first, so a protein FASTA that is also reachable
    /// from one of the scanned directories reuses its `file_id` rather than
    /// getting a second table slot and a second cached handle. Unlike that
    /// loader, a parse failure here does not unwind the whole load — the
    /// caller is expected to warn and continue past a single bad protein
    /// FASTA — so the `files` push is deferred until `load_fai_index` actually
    /// succeeds; pushing first (as `load_index_from_dirs` safely does, because
    /// its own failures propagate with `?` and unwind) would otherwise leave
    /// an orphan table slot referenced by no index entry.
    pub(crate) fn add_protein_fasta(
        &mut self,
        fasta: &Path,
        fai: &Path,
    ) -> Result<usize, FerroError> {
        let canonical = fasta.canonicalize().map_err(|e| FerroError::Io {
            msg: format!(
                "Failed to canonicalize protein FASTA '{}': {}",
                fasta.display(),
                e
            ),
        })?;
        let existing_id =
            self.files.iter().position(|f| f == &canonical).map(|pos| {
                u32::try_from(pos).expect("table index already bounded by u32 fits u32")
            });
        let file_id = match existing_id {
            Some(id) => id,
            None => u32::try_from(self.files.len()).map_err(|_| FerroError::Io {
                msg: "too many FASTA files in the reference (max 4294967296)".to_string(),
            })?,
        };
        // Computed before the push below, and the push itself happens only
        // once `load_fai_index` has actually succeeded (see the doc comment
        // above) — a failure here propagates via `?` before any push occurs,
        // so `file_id` staying valid does not depend on the push happening.
        let idx = load_fai_index(fai, file_id)?;
        let added = idx.len();
        if existing_id.is_none() {
            self.files.push(canonical);
        }
        self.protein_index.extend(idx);
        Ok(added)
    }

    /// Look up a sequence entry by exact accession.
    pub(crate) fn entry(&self, name: &str) -> Option<FastaIndexEntry> {
        self.index.get(name).copied()
    }

    /// Resolve a `file_id` to its FASTA path.
    pub(crate) fn path_for(&self, file_id: u32) -> Option<&Path> {
        self.files.get(file_id as usize).map(PathBuf::as_path)
    }
}

/// Load a FASTA index (.fai) file, tagging every entry with `file_id`.
///
/// The caller owns path canonicalization (see [`load_index_from_dirs`]) — it
/// happens once per FASTA file rather than once per entry.
pub(super) fn load_fai_index(
    fai_path: &Path,
    file_id: u32,
) -> Result<FxHashMap<String, FastaIndexEntry>, FerroError> {
    let file = File::open(fai_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open FAI file: {}", e),
    })?;
    let reader = BufReader::new(file);

    let mut index: FxHashMap<String, FastaIndexEntry> = FxHashMap::default();

    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Failed to read FAI line: {}", e),
        })?;

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }

        // Parse index values with proper error handling instead of silent defaults
        let name = fields[0].to_string();
        let length: u64 = fields[1].parse().map_err(|_| FerroError::Io {
            msg: format!(
                "Invalid length '{}' in FAI for sequence '{}'",
                fields[1], name
            ),
        })?;
        let offset: u64 = fields[2].parse().map_err(|_| FerroError::Io {
            msg: format!(
                "Invalid offset '{}' in FAI for sequence '{}'",
                fields[2], name
            ),
        })?;
        let line_bases: u64 = fields[3].parse().map_err(|_| FerroError::Io {
            msg: format!(
                "Invalid line_bases '{}' in FAI for sequence '{}'",
                fields[3], name
            ),
        })?;
        let line_bytes: u64 = fields[4].parse().map_err(|_| FerroError::Io {
            msg: format!(
                "Invalid line_bytes '{}' in FAI for sequence '{}'",
                fields[4], name
            ),
        })?;

        // Validate that critical fields are non-zero to prevent divide-by-zero
        if line_bases == 0 || line_bytes == 0 {
            return Err(FerroError::Io {
                msg: format!(
                    "Invalid FAI entry for '{}': line_bases={}, line_bytes={} (must be > 0)",
                    name, line_bases, line_bytes
                ),
            });
        }

        let entry = FastaIndexEntry {
            file_id,
            length,
            offset,
            line_bases,
            line_bytes,
        };

        index.insert(name, entry);
    }

    Ok(index)
}

/// Load and combine the FASTA indexes for `dirs`, in order.
///
/// Directories are processed in the order given and entries are inserted into a
/// single map, so a later directory wins for an accession present in more than
/// one — the same last-wins precedence the previous per-directory-then-`extend`
/// construction had, and which real references depend on (`supplemental/`
/// overlaps `transcripts/` on ~80k accessions).
///
/// Returns the combined index and the file table (`files.len()` is the number
/// of distinct physical FASTA files backing the index — canonical-path dedupe
/// means a FASTA reachable by two on-disk paths has its `.fai` parsed twice
/// but occupies one table slot); entries are tagged with a `file_id` indexing
/// into this table rather than each carrying a cloned path.
pub(super) fn load_index_from_dirs(
    dirs: &[&Path],
) -> Result<(FxHashMap<String, FastaIndexEntry>, Vec<PathBuf>), FerroError> {
    let mut index: FxHashMap<String, FastaIndexEntry> = FxHashMap::default();
    let mut files: Vec<PathBuf> = Vec::new();
    // Reuses a `file_id` when the same physical file is reached by more than
    // one on-disk path (e.g. a symlink, or the same FASTA listed under two
    // manifest-scanned directories) — without this, two table slots holding
    // one canonical path would get two distinct ids and therefore two cached
    // handles for one physical file, breaking the invariant the canonicalize
    // call below documents.
    let mut file_ids_by_canonical_path: FxHashMap<PathBuf, u32> = FxHashMap::default();

    for dir in dirs {
        if !dir.exists() {
            continue;
        }

        let entries = std::fs::read_dir(dir).map_err(|e| FerroError::Io {
            msg: format!("Failed to read directory {}: {}", dir.display(), e),
        })?;

        // Sort so file ordering — and therefore which file wins a within-directory
        // duplicate — is deterministic. `read_dir` order is unspecified, and 65
        // accessions appear in both genome/GRCh38.fna.fai and GRCh37.fna.fai.
        let mut paths: Vec<PathBuf> = Vec::new();
        for entry in entries {
            let entry = entry.map_err(|e| FerroError::Io {
                msg: format!("Failed to read directory entry: {}", e),
            })?;
            paths.push(entry.path());
        }
        paths.sort();

        let mut entries_from_dir: usize = 0;

        for path in paths {
            // Skip hidden files and macOS AppleDouble sidecars (names beginning
            // with '.', e.g. `._LRG_430.fasta`). These are not real FASTAs — an
            // accompanying `._*.fasta.fai` is a binary extended-attribute blob,
            // and parsing it as a FAI index would abort the entire reference
            // load with a "stream did not contain valid UTF-8" error. Such files
            // are commonly introduced when a reference directory is copied
            // through macOS tooling.
            let is_hidden = path
                .file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.starts_with('.'));
            if is_hidden {
                continue;
            }

            let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
            if ext != "fna" && ext != "fa" && ext != "fasta" {
                continue;
            }

            let fai_path = PathBuf::from(format!("{}.fai", path.display()));
            if !fai_path.exists() {
                continue;
            }

            // Canonicalized once per FASTA rather than cloned into every entry.
            // Resolving symlinks here also keeps the open-file cache keyed on
            // one identity per physical file, and — since the result is
            // trusted as an absolute lookup path — resolves `..`/symlink
            // components before that trust is extended, helping prevent path
            // traversal.
            let canonical = path.canonicalize().map_err(|e| FerroError::Io {
                msg: format!(
                    "Failed to canonicalize FASTA path '{}': {}",
                    path.display(),
                    e
                ),
            })?;
            // Pushed into `files` here, before `load_fai_index` runs — unlike
            // the protein-FASTA loader below, which defers its push until
            // parsing succeeds. The two are a deliberate pair, not a
            // contradiction: a parse failure here is propagated with `?`,
            // which unwinds the whole load, so no orphan table slot can ever
            // be observed; the protein loader instead warns and `continue`s
            // past a parse failure, so it MUST defer the push or it would
            // leave a slot referenced by nothing.
            let file_id = match file_ids_by_canonical_path.get(&canonical) {
                Some(&id) => id,
                None => {
                    let id = u32::try_from(files.len()).map_err(|_| FerroError::Io {
                        msg: "too many FASTA files in the reference (max 4294967296)".to_string(),
                    })?;
                    file_ids_by_canonical_path.insert(canonical.clone(), id);
                    files.push(canonical);
                    id
                }
            };

            let file_index = load_fai_index(&fai_path, file_id)?;
            for (name, entry) in file_index {
                index.insert(name, entry);
                entries_from_dir += 1;
            }
        }

        // A directory that exists but yields no indexed FASTA is an error, not a
        // silent skip: `from_directory` errored on this before the refactor, and
        // silently skipping it here would let a later directory's overlapping
        // accessions be served from the wrong place (e.g. `transcripts/` winning
        // accessions `supplemental/` was supposed to). Contrast with a
        // non-existent directory above, which is legitimately optional and
        // skipped without error.
        if entries_from_dir == 0 {
            return Err(FerroError::Io {
                msg: format!(
                    "No indexed FASTA files found in {}. Run 'ferro-benchmark prepare' first.",
                    dir.display()
                ),
            });
        }
    }

    Ok((index, files))
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::io::Write;

    use tempfile::tempdir;

    #[test]
    fn from_dirs_builds_index_and_derived_maps() {
        let dir = tempdir().unwrap();
        std::fs::write(
            dir.path().join("t.fna"),
            ">NM_000088.3\nAAAA\n>NM_000088.4\nCCCC\n",
        )
        .unwrap();
        std::fs::write(
            dir.path().join("t.fna.fai"),
            "NM_000088.3\t4\t13\t4\t5\nNM_000088.4\t4\t31\t4\t5\n",
        )
        .unwrap();

        let prepared = PreparedIndex::from_dirs(&[dir.path()]).unwrap();

        assert_eq!(prepared.index.len(), 2);
        assert_eq!(prepared.files.len(), 1);
        // Highest version wins the base-accession fallback.
        assert_eq!(
            prepared
                .base_to_versioned
                .get("NM_000088")
                .map(String::as_str),
            Some("NM_000088.4")
        );
        assert!(!prepared.has_genomic_data);
        assert_eq!(prepared.entry("NM_000088.3").unwrap().file_id, 0);
        assert!(prepared.path_for(0).is_some());
        assert!(prepared.path_for(1).is_none());
    }

    #[test]
    fn from_dirs_detects_genomic_data() {
        let dir = tempdir().unwrap();
        std::fs::write(dir.path().join("g.fna"), ">NC_000017.11\nAAAA\n").unwrap();
        std::fs::write(dir.path().join("g.fna.fai"), "NC_000017.11\t4\t14\t4\t5\n").unwrap();

        let prepared = PreparedIndex::from_dirs(&[dir.path()]).unwrap();
        assert!(prepared.has_genomic_data);
    }

    /// #520 edge 3 / add_protein_fasta: exercises the three properties the
    /// method must preserve from the inline `from_manifest_inner` logic it
    /// replaces — canonical-path dedupe against the existing file table, the
    /// `file_id` computed before any push, and a failed parse leaving no
    /// orphan file-table slot.
    #[test]
    fn add_protein_fasta_dedupes_and_defers_push_on_failure() {
        let dir = tempdir().unwrap();
        std::fs::write(dir.path().join("t.fna"), ">NM_1.1\nAAAA\n").unwrap();
        std::fs::write(dir.path().join("t.fna.fai"), "NM_1.1\t4\t13\t4\t5\n").unwrap();
        let mut prepared = PreparedIndex::from_dirs(&[dir.path()]).unwrap();
        assert_eq!(prepared.files.len(), 1);

        let protein_fasta = dir.path().join("p.fasta");
        std::fs::write(&protein_fasta, ">NP_1.1\nMM\n").unwrap();
        let protein_fai = dir.path().join("p.fasta.fai");
        std::fs::write(&protein_fai, "NP_1.1\t2\t12\t2\t3\n").unwrap();

        // A new protein FASTA gets the next file_id and its entries land in
        // protein_index.
        let added = prepared
            .add_protein_fasta(&protein_fasta, &protein_fai)
            .unwrap();
        assert_eq!(added, 1);
        assert_eq!(prepared.files.len(), 2);
        assert_eq!(prepared.protein_index.get("NP_1.1").unwrap().file_id, 1);

        // Re-adding the SAME canonical path reuses the existing file_id rather
        // than growing the file table.
        let added_again = prepared
            .add_protein_fasta(&protein_fasta, &protein_fai)
            .unwrap();
        assert_eq!(added_again, 1);
        assert_eq!(
            prepared.files.len(),
            2,
            "the canonical path is deduped against the existing table"
        );

        // A parse failure (missing .fai) must not leave an orphan file-table slot.
        let other_protein_fasta = dir.path().join("p2.fasta");
        std::fs::write(&other_protein_fasta, ">NP_2.1\nMM\n").unwrap();
        let missing_fai = dir.path().join("does-not-exist.fai");
        assert!(prepared
            .add_protein_fasta(&other_protein_fasta, &missing_fai)
            .is_err());
        assert_eq!(
            prepared.files.len(),
            2,
            "a parse failure must not push an orphan file-table slot"
        );
    }

    #[test]
    fn test_load_fai_index() {
        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("test.fna");
        let fai_path = dir.path().join("test.fna.fai");

        // Create FASTA file
        let mut fasta_file = File::create(&fasta_path).unwrap();
        writeln!(fasta_file, ">NM_000001.1").unwrap();
        writeln!(fasta_file, "ATGCATGCAT").unwrap();

        // Create FAI file
        let mut fai_file = File::create(&fai_path).unwrap();
        writeln!(fai_file, "NM_000001.1\t10\t13\t10\t11").unwrap();

        let index = load_fai_index(&fai_path, 0).unwrap();

        assert!(index.contains_key("NM_000001.1"));
        let entry = index.get("NM_000001.1").unwrap();
        assert_eq!(entry.file_id, 0);
        assert_eq!(entry.length, 10);
        assert_eq!(entry.offset, 13);
    }

    #[test]
    fn build_base_to_versioned_picks_highest_version_deterministically() {
        let entry = |_name: &str| FastaIndexEntry {
            file_id: 0,
            length: 1,
            offset: 0,
            line_bases: 1,
            line_bytes: 2,
        };
        let mut index: FxHashMap<String, FastaIndexEntry> = FxHashMap::default();
        for v in ["NM_000088.3", "NM_000088.4", "NM_000088.10", "NM_000088.2"] {
            index.insert(v.to_string(), entry(v));
        }
        let map = build_base_to_versioned(&index);
        // Highest *numeric* version wins (10 > 4 > 3 > 2) — not lexical, which
        // would pick ".4" over ".10" — and the result is independent of the
        // (hash) iteration order of the index.
        assert_eq!(
            map.get("NM_000088").map(String::as_str),
            Some("NM_000088.10")
        );
    }
}
