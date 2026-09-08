//! Strategy output types: segments (L1) and typed members / partitions (L2).

use serde::{Deserialize, Serialize};

/// A segment (L1 output): an aligned region bounded by robustly-unchanged walls.
/// Half-open ranges into the reference and resulting sequences respectively.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Segment {
    pub ref_start: usize,
    pub ref_end: usize,
    pub res_start: usize,
    pub res_end: usize,
}

/// The HGVS edit kind a member is typed as. `Ord` is derived so members and
/// edit-mix maps have a stable key ordering.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum EditKind {
    Sub,
    Del,
    Ins,
    Delins,
    Dup,
    Inv,
    Identity,
}

/// One HGVS member (L2 output). `ref_start..ref_end` is the half-open reference
/// span consumed; `inserted` carries new bases for Ins/Delins/Dup/Sub and the
/// revcomp payload is implied for Inv (so `inserted` is empty for Inv and Del).
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Member {
    pub kind: EditKind,
    pub ref_start: usize,
    pub ref_end: usize,
    #[serde(with = "serde_bytes_lossy")]
    pub inserted: Vec<u8>,
}

/// An ordered partition of a variant into members. Structural equality is the
/// basis for the cross-arm agreement metric.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Partition {
    pub members: Vec<Member>,
}

/// Serialize `Vec<u8>` of ACGT bytes as a UTF-8 string for readable JSONL.
mod serde_bytes_lossy {
    use serde::{Deserialize, Deserializer, Serializer};
    pub fn serialize<S: Serializer>(v: &[u8], s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str(&String::from_utf8_lossy(v))
    }
    pub fn deserialize<'de, D: Deserializer<'de>>(d: D) -> Result<Vec<u8>, D::Error> {
        Ok(String::deserialize(d)?.into_bytes())
    }
}
