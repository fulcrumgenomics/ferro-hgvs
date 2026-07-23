//! Reference content fingerprint (#1001).
//!
//! A path-independent FNV-1a signature over the prepared reference's content:
//! `transcript_count`, sorted `available_prefixes`, artifact BASENAMES (whose
//! names encode source build/version), and FNV content-stamps of the small
//! derived artifacts ([`CONTENT_STAMPED_ARTIFACTS`]). Excludes `prepared_at`
//! (#933) and the bulk cdot/genome/FASTA artifacts (prohibitive to hash; guarded
//! by version-tagged basenames). Promoted here from a test helper so
//! `ferro prepare` can stamp it and the loader can verify it — this is now the
//! single source of truth for the identity (previously duplicated in tests).

use std::fs;
use std::path::Path;

/// FNV-1a (64-bit) hex digest of `bytes` — small, dependency-free, and stable
/// across platforms/toolchains. Shared by the reference-identity signature and
/// the #905 content stamps so they compose with one identical hash.
pub fn fnv1a_hex(bytes: &[u8]) -> String {
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        hash ^= b as u64;
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    format!("{hash:016x}")
}

/// The tracked artifacts whose *content* is stamped into the reference identity
/// (#905 Part 2): every manifest artifact small enough to FNV-hash on each
/// identity computation. A stamp catches an *in-place* content change — a rewrite
/// that keeps the same filename (a standalone derive/backfill example that
/// bypasses `manifest.save()`, or a hand-edit) — which the Part-1 basename set
/// cannot, and which matters most for the artifacts with *unversioned* basenames
/// (`LRG_RefSeqGene`, `lrg_refseq_mapping.txt`, the derived/legacy JSONs), where
/// Part 1 gives no content protection at all.
///
/// The cut is size, not artifact class, so there are no arbitrary sibling gaps:
/// deliberately excluded are the bulk artifacts (all ≥100 MB on a live reference)
/// — the cdot JSONs (~200–500 MB), the genome FASTAs (~3 GB each),
/// `supplemental_fasta` (~1.3 GB), and the transcript/protein/refseqgene/LRG/
/// ensembl sequence-FASTA lists (~100–600 MB). Hashing any of those on every
/// `reference_identity()` call would be prohibitive I/O; they carry
/// build/version-tagged basenames (Part 1), and a real data update bumps that
/// version — so a *same-basename* in-place rewrite of a bulk artifact is the one
/// content change the identity does not catch, an accepted gap now that
/// `prepared_at` is no longer hashed (#933). Absent fields contribute no stamp.
///
/// NOTE: `refseqgene_alignments`, `refseqgene_alignments_grch37`, and
/// `refseqgene_summary` (LRG_RefSeqGene) are ~2.5 MB each on a live reference —
/// not tiny; the ~8 MB total is acceptable because the loader already reads
/// most of these during provider construction (a page-cache-warm re-read). The
/// bulk cdot/genome/transcript FASTAs are deliberately excluded.
pub const CONTENT_STAMPED_ARTIFACTS: &[&str] = &[
    "derived_transcript_placements",
    "derived_refseqgene_placements",
    "ng_hosted_transcripts",
    "backfill_transcripts_fasta",
    "refseqgene_summary",
    "refseqgene_alignments",
    "refseqgene_alignments_grch37",
    "assembly_report",
    "assembly_report_grch37",
    "lrg_refseq_mapping",
    "legacy_transcripts_fasta",
    "legacy_transcripts_metadata",
    "legacy_genbank_fasta",
    "legacy_genbank_metadata",
    "canonical_overrides",
];

/// Read each present content-stamped artifact (resolving its manifest-relative
/// path against `base_dir`) and inject a `derived_artifact_stamps` object
/// (field name → FNV-1a hex of the bytes) into `manifest`. A missing or
/// unreadable artifact contributes no stamp (so a removed artifact also changes
/// the signature); a no-op when none are present.
pub fn inject_content_stamps(manifest: &mut serde_json::Value, base_dir: &Path) {
    let mut stamps = serde_json::Map::new();
    for &field in CONTENT_STAMPED_ARTIFACTS {
        let Some(rel) = manifest.get(field).and_then(serde_json::Value::as_str) else {
            continue;
        };
        let rel_path = Path::new(rel);
        let path = if rel_path.is_absolute() {
            rel_path.to_path_buf()
        } else {
            base_dir.join(rel_path)
        };
        if let Ok(bytes) = fs::read(&path) {
            stamps.insert(
                field.to_string(),
                serde_json::Value::String(fnv1a_hex(&bytes)),
            );
        }
    }
    if !stamps.is_empty() {
        if let Some(obj) = manifest.as_object_mut() {
            obj.insert(
                "derived_artifact_stamps".to_string(),
                serde_json::Value::Object(stamps),
            );
        }
    }
}

/// Compute the reference identity (#764) from a parsed manifest value.
///
/// Builds a canonical, path-independent content signature and FNV-1a hashes it.
/// The signature draws only on fields that move when the prepared-reference
/// *content* changes and that are invariant to machine-local paths (notably
/// EXCLUDING `prepared_at`, a wall-clock stamp — see the inline note in the body
/// and #933):
/// - `transcript_count` and the sorted `available_prefixes` — content-derived;
/// - the *basenames* (not full paths) of **every** content-bearing artifact the
///   manifest tracks (cdot / genome / transcript / derived-placement / NG-hosted
///   / RefSeqGene / supplemental / backfill / legacy / LRG / ensembl / canonical),
///   whose names encode the source build/version — Part 1 of #905, catching an
///   artifact being added/removed or version-renamed;
/// - the per-derived-artifact **content stamps** (`derived_artifact_stamps`) —
///   Part 2 of #905 (the root fix): an FNV-1a of each out-of-band-regenerated
///   artifact's bytes, so an **in-place** content change bumps the identity even
///   when `prepared_at` and every basename are unchanged (the gap that let the
///   #790/#795 `derived_transcript_placements` and #859/#871 `ng_hosted_transcripts`
///   additions leave the identity stale).
///
/// Split out from [`reference_identity`] so it can be unit-tested without a live
/// prepared reference.
pub fn reference_identity_from_manifest(manifest: &serde_json::Value) -> String {
    // Basename of a string-valued manifest field, if present.
    let basename = |key: &str| -> Option<String> {
        manifest
            .get(key)
            .and_then(serde_json::Value::as_str)
            .map(|p| {
                Path::new(p)
                    .file_name()
                    .map(|n| n.to_string_lossy().into_owned())
                    .unwrap_or_else(|| p.to_string())
            })
    };
    // Basenames of an array-valued manifest field, if present.
    let basenames = |key: &str| -> Vec<String> {
        manifest
            .get(key)
            .and_then(serde_json::Value::as_array)
            .map(|arr| {
                arr.iter()
                    .filter_map(serde_json::Value::as_str)
                    .map(|p| {
                        Path::new(p)
                            .file_name()
                            .map(|n| n.to_string_lossy().into_owned())
                            .unwrap_or_else(|| p.to_string())
                    })
                    .collect()
            })
            .unwrap_or_default()
    };

    // Collect the content signature components, then canonicalize: every entry is
    // `key=value`, the multi-valued artifact lists are sorted, and the whole set
    // is joined deterministically so the same content always yields the same
    // signature regardless of manifest key ordering or host paths.
    let mut parts: Vec<String> = Vec::new();
    // `prepared_at` is deliberately EXCLUDED from the signature (#933): it is a
    // wall-clock stamp refreshed on every `save()`, so hashing it made the
    // identity non-reproducible — a re-bless of byte-identical data produced a
    // new identity and forced a fixture re-record for no content change. The
    // identity is now a pure content signature (transcript_count + artifact
    // basenames, whose names encode source build/version, + FNV content-stamps
    // of the derived artifacts). Tradeoff: a same-basename, same-count artifact
    // whose bytes changed in place is caught only for the content-stamped
    // derived artifacts, not the big cdot/genome/transcript FASTAs — but a real
    // data update bumps the version in the filename (a basename change), so this
    // is a negligible gap versus the reproducibility win.
    if let Some(count) = manifest
        .get("transcript_count")
        .and_then(serde_json::Value::as_u64)
    {
        parts.push(format!("transcript_count={count}"));
    }
    // Part 1 (#905): every single-path content-bearing artifact's basename.
    // Expanded from the original cdot/genome-only set to cover the derived,
    // NG-hosted, RefSeqGene, supplemental, backfill, legacy, LRG, ensembl, and
    // canonical artifacts — so adding/removing an artifact, or a version-name
    // change on any of them, bumps the identity (previously only the cdot/genome/
    // transcript basenames did). Mirrors `ReferenceManifest::for_each_path`.
    for key in [
        "cdot_json",
        "cdot_grch37_json",
        "ensembl_cdot_json",
        "ensembl_cdot_grch37_json",
        "genome_fasta",
        "genome_grch37_fasta",
        "refseqgene_alignments",
        "refseqgene_alignments_grch37",
        "assembly_report",
        "assembly_report_grch37",
        "derived_refseqgene_placements",
        "ng_hosted_transcripts",
        "derived_transcript_placements",
        "refseqgene_summary",
        "lrg_refseq_mapping",
        "supplemental_fasta",
        "backfill_transcripts_fasta",
        "legacy_transcripts_fasta",
        "legacy_transcripts_metadata",
        "legacy_genbank_fasta",
        "legacy_genbank_metadata",
        "canonical_overrides",
    ] {
        if let Some(name) = basename(key) {
            parts.push(format!("{key}={name}"));
        }
    }
    for key in [
        "transcript_fastas",
        "protein_fastas",
        "refseqgene_fastas",
        "lrg_fastas",
        "lrg_xmls",
        "ensembl_transcript_fastas",
    ] {
        let mut names = basenames(key);
        names.sort();
        if !names.is_empty() {
            parts.push(format!("{key}=[{}]", names.join(",")));
        }
    }
    // `available_prefixes` is content-derived (which accession namespaces the
    // reference serves); include it, sorted for order-independence.
    if let Some(arr) = manifest
        .get("available_prefixes")
        .and_then(serde_json::Value::as_array)
    {
        let mut prefixes: Vec<String> = arr
            .iter()
            .filter_map(|v| v.as_str().map(str::to_owned))
            .collect();
        prefixes.sort();
        if !prefixes.is_empty() {
            parts.push(format!("available_prefixes=[{}]", prefixes.join(",")));
        }
    }
    // Part 2 (#905, root fix): per-derived-artifact content stamps folded in from
    // the `derived_artifact_stamps` object — an FNV-1a of each out-of-band-
    // regenerated artifact's bytes, computed from the files at identity time by
    // `inject_content_stamps` (the live [`reference_identity`] path) or supplied
    // directly (unit tests). An **in-place** content regeneration of a
    // constant-named artifact (without moving `prepared_at` or a basename)
    // changes its stamp here, so the identity bumps — the case the basename part
    // cannot catch.
    if let Some(stamps) = manifest
        .get("derived_artifact_stamps")
        .and_then(serde_json::Value::as_object)
    {
        let mut stamp_parts: Vec<String> = stamps
            .iter()
            .filter_map(|(k, v)| v.as_str().map(|s| format!("stamp:{k}={s}")))
            .collect();
        stamp_parts.sort();
        parts.extend(stamp_parts);
    }
    parts.sort();
    let signature = parts.join("\n");
    fnv1a_hex(signature.as_bytes())
}

/// Compute the identity from a manifest value plus the on-disk artifacts under
/// `reference_dir`: clones `manifest`, injects the derived-artifact content
/// stamps read from disk, then hashes. This is the entry point used by
/// `ferro prepare` (stamp) and the loader (verify).
pub fn reference_identity(reference_dir: &Path, manifest: &serde_json::Value) -> String {
    let mut with_stamps = manifest.clone();
    inject_content_stamps(&mut with_stamps, reference_dir);
    reference_identity_from_manifest(&with_stamps)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fnv1a_known_vector() {
        // FNV-1a of the empty string is the offset basis.
        assert_eq!(fnv1a_hex(b""), "cbf29ce484222325");
    }

    #[test]
    fn schema_version_and_identity_field_are_excluded_from_the_hash() {
        let base = serde_json::json!({
            "transcript_count": 3,
            "cdot_json": "cdot-0.2.32.refseq.GRCh38.json",
            "available_prefixes": ["NM_", "NC_"],
        });
        let id = reference_identity_from_manifest(&base);
        let mut with_extras = base.clone();
        with_extras["manifest_schema_version"] = serde_json::json!(2);
        with_extras["reference_identity"] = serde_json::json!("deadbeefdeadbeef");
        // Neither field is in a hashed key list, so the digest is unchanged.
        assert_eq!(reference_identity_from_manifest(&with_extras), id);
    }

    // Reference identity from a parsed manifest (#764).

    #[test]
    fn reference_identity_is_deterministic_and_path_independent() {
        // Two manifests with identical *content* but different absolute paths and
        // a different manifest key order must yield the same identity.
        let a = serde_json::json!({
            "prepared_at": "2026-06-17T09:44:18.428902+00:00",
            "transcript_count": 273423,
            "cdot_json": "/host-a/ferro/cdot/cdot-0.2.32.refseq.GRCh38.json",
            "genome_fasta": "/host-a/ferro/genome/GRCh38.fna",
            "transcript_fastas": [
                "/host-a/ferro/transcripts/human.2.rna.fna.gz",
                "/host-a/ferro/transcripts/human.1.rna.fna.gz",
            ],
        });
        let b = serde_json::json!({
            // Different host paths, and the transcript list in a different order.
            "transcript_fastas": [
                "/elsewhere/data/transcripts/human.1.rna.fna.gz",
                "/elsewhere/data/transcripts/human.2.rna.fna.gz",
            ],
            "genome_fasta": "/elsewhere/data/genome/GRCh38.fna",
            "cdot_json": "/elsewhere/data/cdot/cdot-0.2.32.refseq.GRCh38.json",
            "transcript_count": 273423,
            "prepared_at": "2026-06-17T09:44:18.428902+00:00",
        });
        let id_a = reference_identity_from_manifest(&a);
        // Deterministic: same input twice is the same identity.
        assert_eq!(id_a, reference_identity_from_manifest(&a));
        // Path-independent and order-independent: identical content matches.
        assert_eq!(id_a, reference_identity_from_manifest(&b));
        // Looks like a 16-hex-digit FNV-1a digest.
        assert_eq!(id_a.len(), 16, "identity is a 64-bit hex digest: {id_a}");
        assert!(
            id_a.chars().all(|c| c.is_ascii_hexdigit()),
            "hex digest: {id_a}"
        );
    }

    #[test]
    fn reference_identity_changes_when_content_changes() {
        let base = serde_json::json!({
            "prepared_at": "2026-06-17T09:44:18.428902+00:00",
            "transcript_count": 273423,
            "cdot_json": "cdot/cdot-0.2.32.refseq.GRCh38.json",
            "genome_fasta": "genome/GRCh38.fna",
        });
        let id_base = reference_identity_from_manifest(&base);

        // A re-bless bumps `prepared_at` but leaves every artifact byte-identical:
        // the identity must be UNCHANGED (#933). `prepared_at` is not hashed, so
        // re-blessing the same data is idempotent and forces no fixture re-record.
        let mut reprepared = base.clone();
        reprepared["prepared_at"] = serde_json::json!("2026-06-18T00:00:00.000000+00:00");
        assert_eq!(
            id_base,
            reference_identity_from_manifest(&reprepared),
            "a re-bless of identical data (new prepared_at only) must NOT change the identity"
        );

        // A different cdot version (encoded in the filename) → different identity.
        let mut new_cdot = base.clone();
        new_cdot["cdot_json"] = serde_json::json!("cdot/cdot-0.2.33.refseq.GRCh38.json");
        assert_ne!(
            id_base,
            reference_identity_from_manifest(&new_cdot),
            "a different cdot version must change the identity"
        );

        // A different transcript count → different identity.
        let mut new_count = base.clone();
        new_count["transcript_count"] = serde_json::json!(273424);
        assert_ne!(
            id_base,
            reference_identity_from_manifest(&new_count),
            "a different transcript_count must change the identity"
        );
    }

    #[test]
    fn reference_identity_changes_when_a_derived_artifact_is_added_or_renamed() {
        // Part 1 (#905): a content-bearing artifact outside the original
        // cdot/genome/transcript set now contributes — adding it, and
        // version-renaming it, each change the identity.
        let base = serde_json::json!({
            "prepared_at": "2026-06-17T09:44:18.428902+00:00",
            "transcript_count": 273423,
            "cdot_json": "cdot/cdot-0.2.32.refseq.GRCh38.json",
            "genome_fasta": "genome/GRCh38.fna",
        });
        let id_base = reference_identity_from_manifest(&base);

        let mut with_derived = base.clone();
        with_derived["derived_transcript_placements"] =
            serde_json::json!("derived/derived_transcript_placements.json");
        let id_added = reference_identity_from_manifest(&with_derived);
        assert_ne!(
            id_base, id_added,
            "adding a derived-placements artifact must change the identity"
        );

        let mut renamed = with_derived.clone();
        renamed["derived_transcript_placements"] =
            serde_json::json!("derived/derived_transcript_placements.v2.json");
        assert_ne!(
            id_added,
            reference_identity_from_manifest(&renamed),
            "a version-name change on the artifact must change the identity"
        );
    }

    #[test]
    fn reference_identity_changes_when_a_content_stamp_changes_in_place() {
        // Part 2 (#905) acceptance: an in-place content change of a
        // constant-named derived artifact — surfaced only via its
        // `derived_artifact_stamps` entry — bumps the identity even though every
        // basename, `prepared_at`, and `transcript_count` are byte-identical.
        let base = serde_json::json!({
            "prepared_at": "2026-06-17T09:44:18.428902+00:00",
            "transcript_count": 273423,
            "cdot_json": "cdot/cdot-0.2.32.refseq.GRCh38.json",
            "genome_fasta": "genome/GRCh38.fna",
            "derived_transcript_placements": "derived/derived_transcript_placements.json",
            "derived_artifact_stamps": { "derived_transcript_placements": "0123456789abcdef" },
        });
        let id_base = reference_identity_from_manifest(&base);

        let mut restamped = base.clone();
        restamped["derived_artifact_stamps"] =
            serde_json::json!({ "derived_transcript_placements": "fedcba9876543210" });
        assert_ne!(
            id_base,
            reference_identity_from_manifest(&restamped),
            "an in-place content change (new stamp, unchanged basename) must bump the identity"
        );

        // A manifest with no stamps map (prepared before #905) is still
        // deterministic and distinct from a stamped one.
        let mut no_stamps = base.clone();
        no_stamps
            .as_object_mut()
            .unwrap()
            .remove("derived_artifact_stamps");
        let id_no_stamps = reference_identity_from_manifest(&no_stamps);
        assert_ne!(id_base, id_no_stamps);
        assert_eq!(id_no_stamps, reference_identity_from_manifest(&no_stamps));
    }

    #[test]
    fn inject_content_stamps_hashes_files_and_tracks_in_place_changes() {
        // #905 Part 2 computation: stamps are the FNV-1a of the artifact bytes,
        // read fresh from disk (resolving the manifest-relative path), and an
        // in-place rewrite changes the stamp — hence the identity — with the
        // manifest value otherwise untouched.
        let dir = tempfile::tempdir().unwrap();
        let base = dir.path();
        let artifact = base.join("derived_transcript_placements.json");
        std::fs::write(&artifact, br#"{"v":1}"#).unwrap();
        let manifest = || {
            serde_json::json!({
                "derived_transcript_placements": "derived_transcript_placements.json",
            })
        };

        let mut m1 = manifest();
        inject_content_stamps(&mut m1, base);
        let s1 = m1["derived_artifact_stamps"]["derived_transcript_placements"]
            .as_str()
            .unwrap()
            .to_owned();
        assert_eq!(s1, fnv1a_hex(br#"{"v":1}"#));

        std::fs::write(&artifact, br#"{"v":2}"#).unwrap();
        let mut m2 = manifest();
        inject_content_stamps(&mut m2, base);
        let s2 = m2["derived_artifact_stamps"]["derived_transcript_placements"]
            .as_str()
            .unwrap();
        assert_ne!(s1, s2, "an in-place content change must change the stamp");
        assert_ne!(
            reference_identity_from_manifest(&m1),
            reference_identity_from_manifest(&m2),
            "the in-place change must bump the identity"
        );

        // A wired-but-absent artifact contributes no stamp.
        let mut m3 = serde_json::json!({ "derived_transcript_placements": "missing.json" });
        inject_content_stamps(&mut m3, base);
        assert!(
            m3.get("derived_artifact_stamps").is_none(),
            "an absent artifact must not add a stamps map"
        );
    }
}
