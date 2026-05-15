//! HGVS variant types
//!
//! The main variant enum and supporting types for representing
//! complete HGVS variant descriptions.

use super::edit::{NaEdit, ProteinEdit};
use super::interval::{CdsInterval, GenomeInterval, ProtInterval, RnaInterval, TxInterval};
use super::uncertainty::Mu;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::sync::Arc;

/// Interned accession prefixes for common RefSeq and Ensembl types.
/// Using static `Arc<str>` avoids allocating new strings for common prefixes.
mod interned {
    use std::sync::{Arc, OnceLock};

    macro_rules! interned_prefix {
        ($name:ident, $value:literal) => {
            pub fn $name() -> Arc<str> {
                static INSTANCE: OnceLock<Arc<str>> = OnceLock::new();
                INSTANCE.get_or_init(|| Arc::from($value)).clone()
            }
        };
    }

    // RefSeq prefixes
    interned_prefix!(nc, "NC");
    interned_prefix!(ng, "NG");
    interned_prefix!(nt, "NT");
    interned_prefix!(nw, "NW");
    interned_prefix!(nm, "NM");
    interned_prefix!(nr, "NR");
    interned_prefix!(np, "NP");
    interned_prefix!(xm, "XM");
    interned_prefix!(xr, "XR");
    interned_prefix!(xp, "XP");

    // Ensembl prefixes
    interned_prefix!(enst, "ENST");
    interned_prefix!(ensg, "ENSG");
    interned_prefix!(ensp, "ENSP");
    interned_prefix!(ense, "ENSE");
    interned_prefix!(ensr, "ENSR");

    // LRG prefix
    interned_prefix!(lrg, "LRG");

    // Empty string for assembly refs
    interned_prefix!(empty, "");

    /// Get an interned prefix if available, otherwise return None
    #[inline]
    pub fn get_prefix(s: &str) -> Option<Arc<str>> {
        match s {
            "NC" => Some(nc()),
            "NG" => Some(ng()),
            "NT" => Some(nt()),
            "NW" => Some(nw()),
            "NM" => Some(nm()),
            "NR" => Some(nr()),
            "NP" => Some(np()),
            "XM" => Some(xm()),
            "XR" => Some(xr()),
            "XP" => Some(xp()),
            "ENST" => Some(enst()),
            "ENSG" => Some(ensg()),
            "ENSP" => Some(ensp()),
            "ENSE" => Some(ense()),
            "ENSR" => Some(ensr()),
            "LRG" => Some(lrg()),
            "" => Some(empty()),
            _ => None,
        }
    }
}

/// Accession number with optional version
///
/// Uses `Arc<str>` internally for cheap cloning - accession strings are frequently
/// cloned during parsing and normalization, and `Arc<str>` reduces this to a
/// simple reference count increment.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Accession {
    /// Accession prefix (e.g., "NC", "NM", "NP", "ENST")
    pub prefix: Arc<str>,
    /// Accession number
    pub number: Arc<str>,
    /// Version number (e.g., .3 in NM_000088.3)
    pub version: Option<u32>,
    /// Whether this is an Ensembl-style accession (no underscore separator)
    #[serde(default)]
    pub ensembl_style: bool,
    /// Assembly name (e.g., "GRCh37", "GRCh38") for assembly/chromosome notation
    #[serde(default)]
    pub assembly: Option<Arc<str>>,
    /// Chromosome name (e.g., "chr1", "chr23") for assembly/chromosome notation
    #[serde(default)]
    pub chromosome: Option<Arc<str>>,
    /// Genomic context accession for compound reference syntax
    /// e.g., NC_000013.11 in NC_000013.11(NM_004119.3):c.…
    #[serde(default)]
    pub genomic_context: Option<Box<Accession>>,
}

impl Accession {
    /// Intern a prefix string, using a static Arc if available
    #[inline]
    fn intern_prefix(prefix: impl Into<Arc<str>>) -> Arc<str> {
        let prefix: Arc<str> = prefix.into();
        // Try to use an interned version to avoid allocation
        interned::get_prefix(&prefix).unwrap_or(prefix)
    }

    pub fn new(
        prefix: impl Into<Arc<str>>,
        number: impl Into<Arc<str>>,
        version: Option<u32>,
    ) -> Self {
        let prefix = Self::intern_prefix(prefix);
        let ensembl_style = Self::is_ensembl_prefix(&prefix);
        // LRG accessions never carry a version per HGVS spec (SVD-WG008).
        // Drop any provided version so Display canonicalizes to the spec form.
        let version = if Self::is_lrg_prefix(&prefix) {
            None
        } else {
            version
        };
        Self {
            prefix,
            number: number.into(),
            version,
            ensembl_style,
            assembly: None,
            chromosome: None,
            genomic_context: None,
        }
    }

    /// Create with explicit Ensembl style setting
    pub fn with_style(
        prefix: impl Into<Arc<str>>,
        number: impl Into<Arc<str>>,
        version: Option<u32>,
        ensembl_style: bool,
    ) -> Self {
        let prefix = Self::intern_prefix(prefix);
        // LRG accessions never carry a version per HGVS spec (SVD-WG008).
        let version = if Self::is_lrg_prefix(&prefix) {
            None
        } else {
            version
        };
        Self {
            prefix,
            number: number.into(),
            version,
            ensembl_style,
            assembly: None,
            chromosome: None,
            genomic_context: None,
        }
    }

    /// Create an assembly/chromosome reference (e.g., GRCh37(chr23))
    pub fn from_assembly(assembly: impl Into<Arc<str>>, chromosome: impl Into<Arc<str>>) -> Self {
        Self {
            prefix: interned::empty(),
            number: interned::empty(),
            version: None,
            ensembl_style: false,
            assembly: Some(assembly.into()),
            chromosome: Some(chromosome.into()),
            genomic_context: None,
        }
    }

    /// Create a compound reference with genomic context
    /// e.g., NC_000013.11(NM_004119.3) where `context` is NC_000013.11
    /// and the primary accession is NM_004119.3
    pub fn with_genomic_context(mut self, context: Accession) -> Self {
        self.genomic_context = Some(Box::new(context));
        self
    }

    /// Check if this is an assembly/chromosome style reference
    pub fn is_assembly_ref(&self) -> bool {
        self.assembly.is_some() && self.chromosome.is_some()
    }

    /// Check if a prefix is an Ensembl-style prefix
    pub fn is_ensembl_prefix(prefix: &str) -> bool {
        matches!(prefix, "ENST" | "ENSG" | "ENSP" | "ENSE" | "ENSR")
    }

    /// Check if this accession is an Ensembl accession
    pub fn is_ensembl(&self) -> bool {
        self.ensembl_style || Self::is_ensembl_prefix(&self.prefix)
    }

    /// Check if a prefix is the canonical LRG prefix (`"LRG"`, uppercase).
    ///
    /// Per HGVS spec the LRG prefix is case-sensitive (lowercase `lrg_` is not
    /// recognized as an LRG accession; it parses as a generic chromosome
    /// name). See `assets/hgvs-nomenclature/docs/background/refseq.md`.
    pub fn is_lrg_prefix(prefix: &str) -> bool {
        prefix == "LRG"
    }

    /// Check if this accession is an LRG accession.
    pub fn is_lrg(&self) -> bool {
        Self::is_lrg_prefix(&self.prefix)
    }

    /// Validate Ensembl accession format
    /// Returns true if valid, false if invalid
    pub fn validate_ensembl(&self) -> bool {
        if !self.is_ensembl() {
            return true; // Not an Ensembl ID, skip validation
        }
        // Ensembl IDs should have 11-15 digit numbers
        let digit_count = self.number.len();
        (11..=15).contains(&digit_count) && self.number.chars().all(|c| c.is_ascii_digit())
    }

    /// Infer the *primary* expected variant coordinate type from the accession
    /// prefix and (for LRG) the trailing discriminator inside `number`.
    ///
    /// This is informational, not a validator: it returns the canonical
    /// coordinate type for the accession class. The accession may legally
    /// appear with other coordinate types in the variant block — for
    /// example, `NM_*` always returns `"c"` even though `NM_*:n.…` and
    /// `NM_*:r.…` parse correctly, and `LRG_<N>t<M>` always returns `"c"`
    /// even though `LRG_<N>t<M>:n.…` (non-coding) and
    /// `LRG_<N>t<M>:r.…` (RNA) are valid per the HGVS spec.
    ///
    /// LRG accessions take three forms (see https://www.lrg-sequence.org/faq/
    /// and `assets/hgvs-nomenclature/docs/background/refseq.md`):
    /// - `LRG_<N>`        — the LRG genomic record itself; uses `g.` coordinates.
    /// - `LRG_<N>t<M>`    — transcript M of LRG_N; primary `c.` (coding).
    /// - `LRG_<N>p<M>`    — protein M of LRG_N; uses `p.` coordinates.
    ///
    /// The `t<M>` / `p<M>` discriminator is captured inside `self.number`
    /// (e.g. `prefix="LRG"`, `number="1t1"`), so the dispatch must inspect the
    /// trailing discriminator rather than only the prefix.
    pub fn inferred_variant_type(&self) -> Option<&'static str> {
        match &*self.prefix {
            // RefSeq genomic
            "NC" | "NG" | "NT" | "NW" => Some("g"),
            // RefSeq transcript/mRNA
            "NM" => Some("c"),
            // RefSeq non-coding RNA
            "NR" => Some("n"),
            // RefSeq protein
            "NP" => Some("p"),
            // Ensembl transcript
            "ENST" => Some("c"),
            // Ensembl gene (genomic)
            "ENSG" => Some("g"),
            // Ensembl protein
            "ENSP" => Some("p"),
            // LRG: dispatch on the t<M> / p<M> discriminator captured in `number`.
            // See `is_lrg_prefix` / `is_lrg` for the prefix-class predicate.
            p if Self::is_lrg_prefix(p) => Some(Self::lrg_inferred_variant_type(&self.number)),
            // UniProt protein (single letter prefix: O, P, Q, A-N, R-Z)
            p if p.len() == 1 && p.chars().next().is_some_and(|c| c.is_ascii_uppercase()) => {
                Some("p")
            }
            _ => None,
        }
    }

    /// Map an LRG `number` field (e.g. `"1"`, `"1t1"`, `"1p1"`) to its variant type.
    ///
    /// The discriminator is `t<digits>` for transcripts (coding, `c.`) or
    /// `p<digits>` for protein products (`p.`). A bare numeric `number`
    /// (the LRG record itself) maps to genomic (`g.`). Any malformed input
    /// falls back to `g`, matching the prior behavior for unrecognized shapes.
    fn lrg_inferred_variant_type(number: &str) -> &'static str {
        let bytes = number.as_bytes();
        // Find the discriminator letter, if any. A valid LRG number is one or
        // more digits, optionally followed by `t` or `p` and one or more digits.
        let disc_pos = bytes.iter().position(|&b| b == b't' || b == b'p');
        match disc_pos {
            None => "g",
            Some(pos) => {
                // Validate: digits before, digits after, no further non-digits.
                let before_ok = pos > 0 && bytes[..pos].iter().all(|b| b.is_ascii_digit());
                let after = &bytes[pos + 1..];
                let after_ok = !after.is_empty() && after.iter().all(|b| b.is_ascii_digit());
                if before_ok && after_ok {
                    if bytes[pos] == b't' {
                        "c"
                    } else {
                        "p"
                    }
                } else {
                    // Malformed discriminator — be conservative and treat as genomic.
                    "g"
                }
            }
        }
    }

    /// Check if this is a UniProt accession
    pub fn is_uniprot(&self) -> bool {
        // UniProt: single uppercase letter prefix, 5-character number, no version
        self.prefix.len() == 1
            && self
                .prefix
                .chars()
                .next()
                .is_some_and(|c| c.is_ascii_uppercase())
            && self.number.len() == 5
            && self.number.chars().all(|c| c.is_ascii_alphanumeric())
    }

    /// Full accession string without version
    pub fn base(&self) -> String {
        // Assembly/chromosome notation: GRCh37(chr23)
        if let (Some(assembly), Some(chromosome)) = (&self.assembly, &self.chromosome) {
            return format!("{}({})", assembly, chromosome);
        }
        let base = if self.ensembl_style {
            format!("{}{}", self.prefix, self.number)
        } else {
            format!("{}_{}", self.prefix, self.number)
        };
        // Compound reference: context(base)
        if let Some(ctx) = &self.genomic_context {
            format!("{}({})", ctx.full(), base)
        } else {
            base
        }
    }

    /// Full accession string with version if present
    pub fn full(&self) -> String {
        // Assembly/chromosome notation doesn't have versions
        if self.is_assembly_ref() {
            return self.base();
        }
        let full = match self.version {
            Some(v) => {
                if self.ensembl_style {
                    format!("{}{}.{}", self.prefix, self.number, v)
                } else {
                    format!("{}_{}.{}", self.prefix, self.number, v)
                }
            }
            None => {
                if self.ensembl_style {
                    format!("{}{}", self.prefix, self.number)
                } else {
                    format!("{}_{}", self.prefix, self.number)
                }
            }
        };
        // Compound reference: context(full)
        if let Some(ctx) = &self.genomic_context {
            format!("{}({})", ctx.full(), full)
        } else {
            full
        }
    }

    /// Returns the bare accession string for provider/transcript lookups,
    /// stripping any genomic context wrapper. For compound references like
    /// `NC_000013.11(NM_004119.3)`, this returns just `"NM_004119.3"`.
    pub fn transcript_accession(&self) -> String {
        // Assembly/chromosome notation (e.g. GRCh38(chr1)) uses base() directly
        if self.is_assembly_ref() {
            return self.base();
        }

        match self.version {
            Some(v) => {
                if self.ensembl_style {
                    format!("{}{}.{}", self.prefix, self.number, v)
                } else {
                    format!("{}_{}.{}", self.prefix, self.number, v)
                }
            }
            None => {
                if self.ensembl_style {
                    format!("{}{}", self.prefix, self.number)
                } else {
                    format!("{}_{}", self.prefix, self.number)
                }
            }
        }
    }
}

impl fmt::Display for Accession {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.full())
    }
}

/// Location and edit combined
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct LocEdit<L, E> {
    pub location: L,
    pub edit: Mu<E>,
}

impl<L, E> LocEdit<L, E> {
    pub fn new(location: L, edit: E) -> Self {
        Self {
            location,
            edit: Mu::Certain(edit),
        }
    }

    /// Create a predicted/uncertain loc_edit (edit wrapped in parentheses)
    /// Used for patterns like c.(9740C>A) where the variant is predicted
    pub fn new_predicted(location: L, edit: E) -> Self {
        Self {
            location,
            edit: Mu::Uncertain(edit),
        }
    }

    pub fn with_uncertainty(location: L, edit: Mu<E>) -> Self {
        Self { location, edit }
    }
}

impl<L: fmt::Display, E: fmt::Display> fmt::Display for LocEdit<L, E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.location, self.edit)
    }
}

/// Main HGVS variant enum
///
/// Phase of allele variants (cis vs trans vs mosaic vs chimeric)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AllelePhase {
    /// Same chromosome (cis): [var1;var2]
    Cis,
    /// Different chromosomes (trans): [var1];[var2]
    Trans,
    /// Phase unknown
    Unknown,
    /// Mosaic - variant present in some cells: var1/var2
    /// (somatic variation within an individual)
    Mosaic,
    /// Chimeric - different cell lines: var1//var2
    /// (two distinct cell populations in one individual)
    Chimeric,
}

impl fmt::Display for AllelePhase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AllelePhase::Cis => write!(f, "cis"),
            AllelePhase::Trans => write!(f, "trans"),
            AllelePhase::Unknown => write!(f, "unknown"),
            AllelePhase::Mosaic => write!(f, "mosaic"),
            AllelePhase::Chimeric => write!(f, "chimeric"),
        }
    }
}

impl AllelePhase {
    /// Returns true if this is a mosaic phase
    pub fn is_mosaic(&self) -> bool {
        matches!(self, AllelePhase::Mosaic)
    }

    /// Returns true if this is a chimeric phase
    pub fn is_chimeric(&self) -> bool {
        matches!(self, AllelePhase::Chimeric)
    }
}

/// Allele variant containing multiple variants
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct AlleleVariant {
    /// The variants in this allele
    pub variants: Vec<HgvsVariant>,
    /// Phase relationship
    pub phase: AllelePhase,
}

impl AlleleVariant {
    /// Create a new allele variant
    pub fn new(variants: Vec<HgvsVariant>, phase: AllelePhase) -> Self {
        Self { variants, phase }
    }

    /// Create a cis allele (same chromosome)
    pub fn cis(variants: Vec<HgvsVariant>) -> Self {
        Self::new(variants, AllelePhase::Cis)
    }

    /// Create a trans allele (different chromosomes)
    pub fn trans(variants: Vec<HgvsVariant>) -> Self {
        Self::new(variants, AllelePhase::Trans)
    }

    /// Create a mosaic (somatic variation, some cells have variant)
    pub fn mosaic(variants: Vec<HgvsVariant>) -> Self {
        Self::new(variants, AllelePhase::Mosaic)
    }

    /// Create a chimeric (two distinct cell populations)
    pub fn chimeric(variants: Vec<HgvsVariant>) -> Self {
        Self::new(variants, AllelePhase::Chimeric)
    }

    /// Create an unknown phase allele
    pub fn unknown_phase(variants: Vec<HgvsVariant>) -> Self {
        Self::new(variants, AllelePhase::Unknown)
    }

    /// Detect whether the allele contains a self-cancelling (`del` + `dup`)
    /// pair that the HGVS spec disallows ("descriptions removing part of a
    /// reference sequence replacing it with part of the same sequence are
    /// not allowed", `recommendations/general.md` line 47).
    ///
    /// Returns `Some((idx_del, idx_dup))` for the first offending pair found
    /// (deterministic by ascending pair index), or `None` otherwise.
    ///
    /// Detection is conservative: both variants must reference the same
    /// accession AND have simple integer position ranges (no intronic
    /// offsets, no UTR markers, no uncertain boundaries). Variants outside
    /// that shape are skipped — the spec example
    /// `c.[762_768del;767_774dup]` falls inside it.
    pub fn detect_self_cancelling_pair(variants: &[HgvsVariant]) -> Option<(usize, usize)> {
        // Indexed loops are intentional: we need both `i` and `j` to return
        // the pair indices to the caller.
        #[allow(clippy::needless_range_loop)]
        for i in 0..variants.len() {
            let Some((kind_i, range_i, acc_i)) = self_cancelling_descriptor(&variants[i]) else {
                continue;
            };
            for j in (i + 1)..variants.len() {
                let Some((kind_j, range_j, acc_j)) = self_cancelling_descriptor(&variants[j])
                else {
                    continue;
                };
                if acc_i.full() != acc_j.full() {
                    continue;
                }
                let (del_idx, dup_idx) = match (kind_i, kind_j) {
                    (SelfCancellingEditKind::Del, SelfCancellingEditKind::Dup) => (i, j),
                    (SelfCancellingEditKind::Dup, SelfCancellingEditKind::Del) => (j, i),
                    _ => continue,
                };
                if ranges_overlap(range_i, range_j) {
                    return Some((del_idx, dup_idx));
                }
            }
        }
        None
    }

    /// Build a new allele after running the spec-mandated self-cancelling
    /// check. Returns `FerroError::Parse` with
    /// `ErrorCode::SelfCancellingAllele` if the allele violates HGVS
    /// `recommendations/general.md` line 47.
    pub fn try_new_validated(
        variants: Vec<HgvsVariant>,
        phase: AllelePhase,
    ) -> Result<Self, crate::error::FerroError> {
        if let Some((dl, du)) = Self::detect_self_cancelling_pair(&variants) {
            return Err(crate::error::FerroError::Parse {
                pos: 0,
                msg: format!(
                    "Self-cancelling allele: variants at index {} (del) and {} (dup) describe \
                     overlapping reference positions",
                    dl, du
                ),
                diagnostic: Some(Box::new(
                    crate::error::Diagnostic::new()
                        .with_code(crate::error::ErrorCode::SelfCancellingAllele)
                        .with_hint(
                            "HGVS does not allow describing both a deletion and a duplication \
                             of overlapping reference positions in the same allele \
                             (recommendations/general.md line 47); drop one edit or rewrite \
                             as a single delins",
                        ),
                )),
            });
        }
        Ok(Self::new(variants, phase))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SelfCancellingEditKind {
    Del,
    Dup,
}

/// Extract `(edit-kind, [start, end] base range, accession)` if the variant
/// is a candidate for self-cancelling-pair detection. Returns `None` for any
/// variant that is not a simple `del` / `dup` on bare integer positions.
fn self_cancelling_descriptor(
    v: &HgvsVariant,
) -> Option<(SelfCancellingEditKind, (i64, i64), &Accession)> {
    use crate::hgvs::edit::NaEdit;
    let (edit, range, accession) = match v {
        HgvsVariant::Cds(inner) => (
            inner.loc_edit.edit.inner()?,
            cds_simple_range(&inner.loc_edit.location)?,
            &inner.accession,
        ),
        HgvsVariant::Genome(inner) => (
            inner.loc_edit.edit.inner()?,
            genome_simple_range(&inner.loc_edit.location)?,
            &inner.accession,
        ),
        HgvsVariant::Tx(inner) => (
            inner.loc_edit.edit.inner()?,
            tx_simple_range(&inner.loc_edit.location)?,
            &inner.accession,
        ),
        HgvsVariant::Rna(inner) => (
            inner.loc_edit.edit.inner()?,
            rna_simple_range(&inner.loc_edit.location)?,
            &inner.accession,
        ),
        HgvsVariant::Mt(inner) => (
            inner.loc_edit.edit.inner()?,
            genome_simple_range(&inner.loc_edit.location)?,
            &inner.accession,
        ),
        HgvsVariant::Circular(inner) => (
            inner.loc_edit.edit.inner()?,
            genome_simple_range(&inner.loc_edit.location)?,
            &inner.accession,
        ),
        _ => return None,
    };
    let kind = match edit {
        NaEdit::Deletion { .. } => SelfCancellingEditKind::Del,
        NaEdit::Duplication { .. } => SelfCancellingEditKind::Dup,
        _ => return None,
    };
    Some((kind, range, accession))
}

fn cds_simple_range(interval: &crate::hgvs::interval::CdsInterval) -> Option<(i64, i64)> {
    let s = interval.start.inner()?;
    let e = interval.end.inner()?;
    if s.offset.is_some() || e.offset.is_some() || s.utr3 || e.utr3 {
        return None;
    }
    Some((s.base, e.base))
}

fn genome_simple_range(interval: &crate::hgvs::interval::GenomeInterval) -> Option<(i64, i64)> {
    let s = interval.start.inner()?;
    let e = interval.end.inner()?;
    Some((s.base as i64, e.base as i64))
}

fn tx_simple_range(interval: &crate::hgvs::interval::TxInterval) -> Option<(i64, i64)> {
    let s = interval.start.inner()?;
    let e = interval.end.inner()?;
    if s.offset.is_some() || e.offset.is_some() || s.downstream || e.downstream {
        return None;
    }
    Some((s.base, e.base))
}

fn rna_simple_range(interval: &crate::hgvs::interval::RnaInterval) -> Option<(i64, i64)> {
    let s = interval.start.inner()?;
    let e = interval.end.inner()?;
    if s.offset.is_some() || e.offset.is_some() || s.utr3 || e.utr3 {
        return None;
    }
    Some((s.base, e.base))
}

fn ranges_overlap(a: (i64, i64), b: (i64, i64)) -> bool {
    let lo = a.0.max(b.0);
    let hi = a.1.min(b.1);
    lo <= hi
}

/// Write the `ACC:<type>.` prefix for compact allele form.
///
/// Per HGVS spec, bracketed `?` is only valid as a whole-allele marker (`[?]`),
/// never mixed with concrete edits inside the same bracket. Compact form is
/// suppressed when any sub-variant carries the per-variant `?` (e.g. `c.?`).
fn use_compact_form(variants: &[HgvsVariant]) -> bool {
    !variants.is_empty()
        && HgvsVariant::all_share_accession_and_type(variants)
        && !variants.iter().any(HgvsVariant::is_loc_edit_unknown)
}

fn write_compact_prefix(f: &mut fmt::Formatter<'_>, first: &HgvsVariant) -> fmt::Result {
    let accession = first
        .accession()
        .expect("compact form requires an accession; guarded by all_share_accession_and_type");
    write!(f, "{}", accession)?;
    if let Some(gene) = first.gene_symbol() {
        write!(f, "({})", gene)?;
    }
    write!(f, ":{}.", first.variant_type())
}

/// Write `accession(gene):` when `gene_symbol` is set, otherwise `accession:`.
///
/// Per HGVS Nomenclature, the gene-symbol selector is informational disambiguation.
/// ferro's round-trip policy is "preserve when present in input; do not synthesize
/// when absent" (#121). This helper centralizes the selector emission so all
/// per-kind `Display` impls follow the same rule.
fn write_accession_with_optional_gene(
    f: &mut fmt::Formatter<'_>,
    accession: &Accession,
    gene_symbol: Option<&str>,
) -> fmt::Result {
    write!(f, "{}", accession)?;
    if let Some(gene) = gene_symbol {
        write!(f, "({})", gene)?;
    }
    Ok(())
}

/// True iff `variant` carries a position-bound (not whole-entity)
/// identity edit (`<pos>=` or `<pos><ref>=`). This is the LHS shape
/// of the HGVS spec compact mosaic / chimeric form.
fn is_position_identity_lhs(variant: &HgvsVariant) -> bool {
    let edit_is_pos_identity = |edit: &Mu<NaEdit>| {
        matches!(
            edit.inner(),
            Some(NaEdit::Identity {
                whole_entity: false,
                ..
            })
        )
    };
    match variant {
        HgvsVariant::Genome(v) => edit_is_pos_identity(&v.loc_edit.edit),
        HgvsVariant::Cds(v) => edit_is_pos_identity(&v.loc_edit.edit),
        HgvsVariant::Tx(v) => edit_is_pos_identity(&v.loc_edit.edit),
        HgvsVariant::Rna(v) => edit_is_pos_identity(&v.loc_edit.edit),
        HgvsVariant::Mt(v) => edit_is_pos_identity(&v.loc_edit.edit),
        HgvsVariant::Circular(v) => edit_is_pos_identity(&v.loc_edit.edit),
        _ => false,
    }
}

/// True iff `variant` carries a whole-entity identity edit (the
/// canonical shape produced by `create_identity_variant_from`, used
/// in the `var/=` shorthand). On Display this should render as bare
/// `=` rather than the full synthetic-position form.
fn is_whole_entity_identity_rhs(variant: &HgvsVariant) -> bool {
    let edit_is_whole = |edit: &Mu<NaEdit>| {
        matches!(
            edit.inner(),
            Some(NaEdit::Identity {
                whole_entity: true,
                ..
            })
        )
    };
    match variant {
        HgvsVariant::Genome(v) => edit_is_whole(&v.loc_edit.edit),
        HgvsVariant::Cds(v) => edit_is_whole(&v.loc_edit.edit),
        HgvsVariant::Tx(v) => edit_is_whole(&v.loc_edit.edit),
        HgvsVariant::Rna(v) => edit_is_whole(&v.loc_edit.edit),
        HgvsVariant::Mt(v) => edit_is_whole(&v.loc_edit.edit),
        HgvsVariant::Circular(v) => edit_is_whole(&v.loc_edit.edit),
        _ => false,
    }
}

/// True iff `a` and `b` carry the same coord-system arm AND their
/// `loc_edit.location` values compare equal. Used to decide whether
/// the spec compact mosaic Display form applies to a 2-variant Allele.
fn intervals_match_for_compact_mosaic(a: &HgvsVariant, b: &HgvsVariant) -> bool {
    match (a, b) {
        (HgvsVariant::Genome(x), HgvsVariant::Genome(y)) => {
            x.loc_edit.location == y.loc_edit.location
        }
        (HgvsVariant::Cds(x), HgvsVariant::Cds(y)) => x.loc_edit.location == y.loc_edit.location,
        (HgvsVariant::Tx(x), HgvsVariant::Tx(y)) => x.loc_edit.location == y.loc_edit.location,
        (HgvsVariant::Rna(x), HgvsVariant::Rna(y)) => x.loc_edit.location == y.loc_edit.location,
        (HgvsVariant::Mt(x), HgvsVariant::Mt(y)) => x.loc_edit.location == y.loc_edit.location,
        (HgvsVariant::Circular(x), HgvsVariant::Circular(y)) => {
            x.loc_edit.location == y.loc_edit.location
        }
        _ => false,
    }
}

/// Shared Display path for `Mosaic` (`/`) and `Chimeric` (`//`) phases.
///
/// Three branches in priority order (mutually exclusive on each input):
///
/// 1. **Spec compact mosaic / chimeric** — LHS is a position-bound `=`
///    identity edit, RHS shares accession + coord type + interval.
///    Emits `<lhs-full>=/<rhs-bare-edit>`. Per
///    `recommendations/DNA/{substitution,deletion,duplication}.md`.
///
/// 2. **`var/=` shorthand cleanup** — LHS is a non-identity variant,
///    RHS is a whole-entity identity edit (the synthetic shape
///    produced by `create_identity_variant_from`), and the two share
///    accession + coord type. Emits `<lhs-full>/=` instead of
///    expanding the synthetic identity to its `<acc>:<type>.1=` form.
///
/// 3. **Long form** — fallback per-variant join.
fn write_mosaic_or_chimeric(
    f: &mut fmt::Formatter<'_>,
    variants: &[HgvsVariant],
    sep: &str,
) -> fmt::Result {
    if variants.len() == 2 && use_compact_form(variants) {
        let lhs = &variants[0];
        let rhs = &variants[1];

        // Branch 1: spec compact (LHS pos-identity, intervals match).
        if is_position_identity_lhs(lhs) && intervals_match_for_compact_mosaic(lhs, rhs) {
            write!(f, "{}", lhs)?;
            write!(f, "{}", sep)?;
            // Emit the RHS edit alone (no accession, no type prefix, no
            // position) — that is the spec compact RHS shape.
            return write_bare_edit(f, rhs);
        }

        // Branch 2: var/= cleanup. LHS must NOT be position-identity
        // (otherwise branch 1 owns it); RHS must be whole-entity Identity.
        if !is_position_identity_lhs(lhs) && is_whole_entity_identity_rhs(rhs) {
            write!(f, "{}", lhs)?;
            write!(f, "{}", sep)?;
            write!(f, "=")?;
            return Ok(());
        }
    }

    // Branch 3: long-form join.
    for (i, v) in variants.iter().enumerate() {
        if i > 0 {
            write!(f, "{}", sep)?;
        }
        write!(f, "{}", v)?;
    }
    Ok(())
}

/// Write the `NaEdit` portion of `variant` (no accession, no type
/// prefix, no position). Used by `write_mosaic_or_chimeric` for the
/// spec compact RHS.
///
/// Variants that are not NA-coord arms are written with their full
/// `Display` (callers should not invoke this on those — guarded by
/// `use_compact_form`).
fn write_bare_edit(f: &mut fmt::Formatter<'_>, variant: &HgvsVariant) -> fmt::Result {
    match variant {
        HgvsVariant::Genome(v) => write!(f, "{}", v.loc_edit.edit),
        HgvsVariant::Cds(v) => write!(f, "{}", v.loc_edit.edit),
        HgvsVariant::Tx(v) => write!(f, "{}", v.loc_edit.edit),
        HgvsVariant::Rna(v) => write!(f, "{}", v.loc_edit.edit),
        HgvsVariant::Mt(v) => write!(f, "{}", v.loc_edit.edit),
        HgvsVariant::Circular(v) => write!(f, "{}", v.loc_edit.edit),
        _ => write!(f, "{}", variant),
    }
}

impl fmt::Display for AlleleVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Per HGVS spec, `[var]` denotes a multi-variant allele; a singleton
        // is just the inner variant. (Mosaic/chimeric already emit bare via
        // their separator-only loops; this unifies cis/trans/unknown.)
        if self.variants.len() == 1 {
            return self.variants[0].fmt(f);
        }
        match self.phase {
            AllelePhase::Cis => {
                if use_compact_form(&self.variants) {
                    // Compact form: ACC:g.[edit1;edit2]
                    write_compact_prefix(f, &self.variants[0])?;
                    write!(f, "[")?;
                    for (i, v) in self.variants.iter().enumerate() {
                        if i > 0 {
                            write!(f, ";")?;
                        }
                        v.fmt_loc_edit(f)?;
                    }
                    write!(f, "]")
                } else {
                    // Expanded form: [ACC:g.edit1;ACC:g.edit2]
                    write!(f, "[")?;
                    for (i, v) in self.variants.iter().enumerate() {
                        if i > 0 {
                            write!(f, ";")?;
                        }
                        write!(f, "{}", v)?;
                    }
                    write!(f, "]")
                }
            }
            AllelePhase::Trans => {
                if use_compact_form(&self.variants) {
                    // Compact form: ACC:g.[edit1];[edit2]
                    write_compact_prefix(f, &self.variants[0])?;
                    for (i, v) in self.variants.iter().enumerate() {
                        if i > 0 {
                            write!(f, ";")?;
                        }
                        write!(f, "[")?;
                        v.fmt_loc_edit(f)?;
                        write!(f, "]")?;
                    }
                    Ok(())
                } else {
                    // Expanded form: [ACC:g.edit1];[ACC:g.edit2]
                    for (i, v) in self.variants.iter().enumerate() {
                        if i > 0 {
                            write!(f, ";")?;
                        }
                        write!(f, "[{}]", v)?;
                    }
                    Ok(())
                }
            }
            AllelePhase::Unknown => {
                // Unknown-phase compact form has no surrounding brackets, so a
                // singleton would print as `ACC:type.edit` and lose the allele
                // wrapper. Only use compact form when there are 2+ sub-variants.
                if self.variants.len() > 1 && use_compact_form(&self.variants) {
                    // Compact form: ACC:g.edit1(;)edit2
                    write_compact_prefix(f, &self.variants[0])?;
                    for (i, v) in self.variants.iter().enumerate() {
                        if i > 0 {
                            write!(f, "(;)")?;
                        }
                        v.fmt_loc_edit(f)?;
                    }
                    Ok(())
                } else {
                    // Expanded form: [ACC:g.edit1(;)ACC:g.edit2]
                    write!(f, "[")?;
                    for (i, v) in self.variants.iter().enumerate() {
                        if i > 0 {
                            write!(f, "(;)")?;
                        }
                        write!(f, "{}", v)?;
                    }
                    write!(f, "]")
                }
            }
            AllelePhase::Mosaic => write_mosaic_or_chimeric(f, &self.variants, "/"),
            AllelePhase::Chimeric => write_mosaic_or_chimeric(f, &self.variants, "//"),
        }
    }
}

/// Represents all types of HGVS variants:
/// - g. (genomic)
/// - c. (coding DNA)
/// - n. (non-coding transcript)
/// - r. (RNA)
/// - p. (protein)
/// - m. (mitochondrial)
/// - o. (circular DNA) - SVD-WG006
/// - RNA fusion (::) - SVD-WG007
/// - Allele (compound variants)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum HgvsVariant {
    /// Genomic variant (g.)
    Genome(GenomeVariant),
    /// Coding DNA variant (c.)
    Cds(CdsVariant),
    /// Non-coding transcript variant (n.)
    Tx(TxVariant),
    /// RNA variant (r.)
    Rna(RnaVariant),
    /// Protein variant (p.)
    Protein(ProteinVariant),
    /// Mitochondrial variant (m.)
    Mt(MtVariant),
    /// Circular DNA variant (o.) - SVD-WG006
    Circular(CircularVariant),
    /// RNA Fusion variant (::) - SVD-WG007
    RnaFusion(RnaFusionVariant),
    /// Allele/compound variant
    Allele(AlleleVariant),
    /// Null allele marker [0] - indicates no variant on this chromosome (hemizygous)
    NullAllele,
    /// Unknown allele marker [?] - indicates unknown variant on second chromosome
    UnknownAllele,
}

impl HgvsVariant {
    /// Get the accession for this variant, if available.
    ///
    /// For allele variants, returns the accession of the first variant.
    /// For RNA fusion variants, returns the 5' partner accession.
    /// Returns `None` for NullAllele, UnknownAllele, or empty allele variants.
    pub fn accession(&self) -> Option<&Accession> {
        match self {
            HgvsVariant::Genome(v) => Some(&v.accession),
            HgvsVariant::Cds(v) => Some(&v.accession),
            HgvsVariant::Tx(v) => Some(&v.accession),
            HgvsVariant::Rna(v) => Some(&v.accession),
            HgvsVariant::Protein(v) => Some(&v.accession),
            HgvsVariant::Mt(v) => Some(&v.accession),
            HgvsVariant::Circular(v) => Some(&v.accession),
            HgvsVariant::RnaFusion(v) => Some(&v.five_prime.accession),
            HgvsVariant::Allele(a) => a.variants.first().and_then(|v| v.accession()),
            HgvsVariant::NullAllele | HgvsVariant::UnknownAllele => None,
        }
    }

    /// Get the gene-symbol selector for this variant, if one was present in the
    /// original input.
    ///
    /// Mirrors `accession()` and is used by `Display` to preserve the selector
    /// on round-trip (#121). Returns `None` for variant kinds that do not carry
    /// a single top-level selector (`Allele`, `RnaFusion`, `NullAllele`,
    /// `UnknownAllele`); selectors on those nested kinds are emitted by their
    /// own `Display` impls.
    pub fn gene_symbol(&self) -> Option<&str> {
        match self {
            HgvsVariant::Genome(v) => v.gene_symbol.as_deref(),
            HgvsVariant::Cds(v) => v.gene_symbol.as_deref(),
            HgvsVariant::Tx(v) => v.gene_symbol.as_deref(),
            HgvsVariant::Rna(v) => v.gene_symbol.as_deref(),
            HgvsVariant::Protein(v) => v.gene_symbol.as_deref(),
            HgvsVariant::Mt(v) => v.gene_symbol.as_deref(),
            HgvsVariant::Circular(v) => v.gene_symbol.as_deref(),
            HgvsVariant::RnaFusion(_)
            | HgvsVariant::Allele(_)
            | HgvsVariant::NullAllele
            | HgvsVariant::UnknownAllele => None,
        }
    }

    /// Get the variant type string
    pub fn variant_type(&self) -> &'static str {
        match self {
            HgvsVariant::Genome(_) => "g",
            HgvsVariant::Cds(_) => "c",
            HgvsVariant::Tx(_) => "n",
            HgvsVariant::Rna(_) => "r",
            HgvsVariant::Protein(_) => "p",
            HgvsVariant::Mt(_) => "m",
            HgvsVariant::Circular(_) => "o",
            HgvsVariant::RnaFusion(_) => "r::r",
            HgvsVariant::Allele(_) => "allele",
            HgvsVariant::NullAllele => "null",
            HgvsVariant::UnknownAllele => "unknown",
        }
    }

    /// Check if this is an allele (compound) variant
    pub fn is_allele(&self) -> bool {
        matches!(self, HgvsVariant::Allele(_))
    }

    /// Check if this is a null allele marker
    pub fn is_null_allele(&self) -> bool {
        matches!(self, HgvsVariant::NullAllele)
    }

    /// Check if this is an unknown allele marker
    pub fn is_unknown_allele(&self) -> bool {
        matches!(self, HgvsVariant::UnknownAllele)
    }

    /// Format just the position+edit portion (without accession and coordinate prefix).
    ///
    /// For `NM_000088.3:c.459A>G`, this writes `459A>G`. Used by `AlleleVariant::Display`
    /// to emit the spec-correct compact form (`ACC:c.[edit1;edit2]`).
    ///
    /// Variant types that have no simple loc/edit form (`Allele`, `RnaFusion`,
    /// `NullAllele`, `UnknownAllele`) fall back to their full `Display` output. Callers
    /// guard against this via `all_share_accession_and_type`, which excludes those
    /// types from the compact branch.
    fn fmt_loc_edit(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            HgvsVariant::Genome(v) => v.fmt_loc_edit(f),
            HgvsVariant::Cds(v) => v.fmt_loc_edit(f),
            HgvsVariant::Tx(v) => v.fmt_loc_edit(f),
            HgvsVariant::Rna(v) => v.fmt_loc_edit(f),
            HgvsVariant::Protein(v) => v.fmt_loc_edit(f),
            HgvsVariant::Mt(v) => v.fmt_loc_edit(f),
            HgvsVariant::Circular(v) => v.fmt_loc_edit(f),
            // These types have no compact form — fall back to full Display.
            HgvsVariant::RnaFusion(v) => write!(f, "{}", v),
            HgvsVariant::Allele(a) => write!(f, "{}", a),
            HgvsVariant::NullAllele => write!(f, "0"),
            HgvsVariant::UnknownAllele => write!(f, "?"),
        }
    }

    /// True if this variant's loc/edit portion is the per-variant unknown form
    /// (`g.?`, `c.?`, `n.?`, `r.?`, `p.?`, `m.?`, `o.?`).
    ///
    /// The HGVS spec uses bracketed `?` only as a whole-allele marker (`[?]`),
    /// never mixed with concrete edits inside the same bracket. Compact form is
    /// suppressed when any sub-variant matches this so the output isn't visually
    /// ambiguous with the spec-sanctioned `[?]` form.
    pub(crate) fn is_loc_edit_unknown(&self) -> bool {
        match self {
            HgvsVariant::Genome(v) => is_na_edit_unknown(&v.loc_edit.edit),
            HgvsVariant::Cds(v) => is_na_edit_unknown(&v.loc_edit.edit),
            HgvsVariant::Tx(v) => is_na_edit_unknown(&v.loc_edit.edit),
            HgvsVariant::Rna(v) => is_na_edit_unknown(&v.loc_edit.edit),
            HgvsVariant::Mt(v) => is_na_edit_unknown(&v.loc_edit.edit),
            HgvsVariant::Circular(v) => is_na_edit_unknown(&v.loc_edit.edit),
            HgvsVariant::Protein(v) => match &v.loc_edit.edit {
                Mu::Unknown => true,
                Mu::Certain(e) | Mu::Uncertain(e) => e.is_whole_protein_unknown(),
            },
            HgvsVariant::RnaFusion(_)
            | HgvsVariant::Allele(_)
            | HgvsVariant::NullAllele
            | HgvsVariant::UnknownAllele => false,
        }
    }

    /// Check if all variants in a slice share the same accession, coordinate type, and
    /// gene-symbol selector. Used to determine whether the compact allele form can be
    /// used without losing information.
    ///
    /// The compact form (`ACC(GENE):c.[edit1;edit2]`) can carry exactly one
    /// `(GENE)` selector on the prefix. Per the #121 round-trip policy
    /// ("preserve when present; do not synthesize"), the compact form is only
    /// safe when every sub-variant agrees on `gene_symbol`: both `Some(g)`
    /// (same `g`) and all-`None` are safe; any disagreement (different `g`,
    /// or `Some` vs `None`) falls back to the expanded form so each
    /// sub-variant emits its own selector via per-variant `Display`.
    pub(crate) fn all_share_accession_and_type(variants: &[HgvsVariant]) -> bool {
        let Some(first) = variants.first() else {
            return true;
        };
        let first_acc = first.accession();
        let first_type = first.variant_type();
        let first_gene = first.gene_symbol();

        // Don't use compact form for types that aren't simple coordinate-based variants,
        // or when the first variant has no accession (e.g. NullAllele, UnknownAllele)
        if first_acc.is_none() || matches!(first_type, "allele" | "null" | "unknown" | "r::r") {
            return false;
        }

        variants[1..].iter().all(|v| {
            v.variant_type() == first_type
                && v.accession() == first_acc
                && v.gene_symbol() == first_gene
        })
    }
}

/// True if a nucleotide-edit `Mu` wrapper represents the per-variant unknown form
/// (e.g. `c.?`, `r.?`).
fn is_na_edit_unknown(edit: &Mu<NaEdit>) -> bool {
    match edit {
        Mu::Unknown => true,
        Mu::Certain(e) | Mu::Uncertain(e) => e.is_whole_entity_unknown(),
    }
}

impl fmt::Display for HgvsVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            HgvsVariant::Genome(v) => write!(f, "{}", v),
            HgvsVariant::Cds(v) => write!(f, "{}", v),
            HgvsVariant::Tx(v) => write!(f, "{}", v),
            HgvsVariant::Rna(v) => write!(f, "{}", v),
            HgvsVariant::Protein(v) => write!(f, "{}", v),
            HgvsVariant::Mt(v) => write!(f, "{}", v),
            HgvsVariant::Circular(v) => write!(f, "{}", v),
            HgvsVariant::RnaFusion(v) => write!(f, "{}", v),
            HgvsVariant::Allele(a) => write!(f, "{}", a),
            HgvsVariant::NullAllele => write!(f, "0"),
            HgvsVariant::UnknownAllele => write!(f, "?"),
        }
    }
}

/// Genomic variant (g.)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GenomeVariant {
    pub accession: Accession,
    pub gene_symbol: Option<String>,
    pub loc_edit: LocEdit<GenomeInterval, NaEdit>,
}

impl GenomeVariant {
    /// Format just the position+edit portion (without `accession:g.` prefix).
    fn fmt_loc_edit(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // For whole-entity identity (g.=) or unknown (g.?), skip the position
        if let Some(edit) = self.loc_edit.edit.inner() {
            if edit.is_whole_entity() {
                return write!(f, "{}", self.loc_edit.edit);
            }
        }
        write!(f, "{}", self.loc_edit)
    }
}

impl fmt::Display for GenomeVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write_accession_with_optional_gene(f, &self.accession, self.gene_symbol.as_deref())?;
        write!(f, ":g.")?;
        self.fmt_loc_edit(f)
    }
}

/// Coding DNA variant (c.)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CdsVariant {
    pub accession: Accession,
    pub gene_symbol: Option<String>,
    pub loc_edit: LocEdit<CdsInterval, NaEdit>,
}

impl CdsVariant {
    /// Format just the position+edit portion (without `accession:c.` prefix).
    fn fmt_loc_edit(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // For whole-entity identity (c.=) or unknown (c.?), skip the position
        if let Some(edit) = self.loc_edit.edit.inner() {
            if edit.is_whole_entity_identity() || edit.is_whole_entity_unknown() {
                return write!(f, "{}", self.loc_edit.edit);
            }
        }
        write!(f, "{}", self.loc_edit)
    }
}

impl fmt::Display for CdsVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write_accession_with_optional_gene(f, &self.accession, self.gene_symbol.as_deref())?;
        write!(f, ":c.")?;
        self.fmt_loc_edit(f)
    }
}

/// Non-coding transcript variant (n.)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct TxVariant {
    pub accession: Accession,
    pub gene_symbol: Option<String>,
    pub loc_edit: LocEdit<TxInterval, NaEdit>,
}

impl TxVariant {
    fn fmt_loc_edit(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.loc_edit)
    }
}

impl fmt::Display for TxVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write_accession_with_optional_gene(f, &self.accession, self.gene_symbol.as_deref())?;
        write!(f, ":n.")?;
        self.fmt_loc_edit(f)
    }
}

/// RNA variant (r.)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct RnaVariant {
    pub accession: Accession,
    pub gene_symbol: Option<String>,
    pub loc_edit: LocEdit<RnaInterval, NaEdit>,
}

impl RnaVariant {
    /// Format just the position+edit portion (without `accession:r.` prefix).
    ///
    /// RNA uses lowercase nucleotides per HGVS spec.
    fn fmt_loc_edit(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.loc_edit.edit {
            Mu::Certain(edit) => {
                // For whole-entity patterns (r.=, r.?, r.spl, r.0), skip the position
                if edit.is_whole_entity() {
                    write!(f, "{}", edit.to_rna_string())
                } else {
                    write!(f, "{}{}", self.loc_edit.location, edit.to_rna_string())
                }
            }
            Mu::Uncertain(edit) => {
                // For whole-entity patterns, skip the position
                if edit.is_whole_entity() {
                    write!(f, "({})", edit.to_rna_string())
                } else {
                    write!(f, "{}({})", self.loc_edit.location, edit.to_rna_string())
                }
            }
            Mu::Unknown => write!(f, "?"),
        }
    }
}

impl fmt::Display for RnaVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write_accession_with_optional_gene(f, &self.accession, self.gene_symbol.as_deref())?;
        write!(f, ":r.")?;
        self.fmt_loc_edit(f)
    }
}

/// Protein variant (p.)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ProteinVariant {
    pub accession: Accession,
    pub gene_symbol: Option<String>,
    pub loc_edit: LocEdit<ProtInterval, ProteinEdit>,
}

impl ProteinVariant {
    /// Format just the position+edit portion (without `accession:p.` prefix).
    fn fmt_loc_edit(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // For whole-protein identity (p.= or p.(=)), no-protein (p.0), or whole-protein unknown (p.?), skip the position
        if let Some(edit) = self.loc_edit.edit.inner() {
            if edit.is_whole_protein_identity()
                || edit.is_no_protein()
                || edit.is_whole_protein_unknown()
            {
                return write!(f, "{}", self.loc_edit.edit);
            }
        }
        // For predicted protein changes (uncertain edit), wrap position+edit in parentheses
        // e.g., p.(Arg248Gln) instead of p.Arg248(Gln)
        if self.loc_edit.edit.is_uncertain() {
            if let Some(edit) = self.loc_edit.edit.inner() {
                return write!(f, "({}{})", self.loc_edit.location, edit);
            }
        }
        write!(f, "{}", self.loc_edit)
    }
}

impl fmt::Display for ProteinVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write_accession_with_optional_gene(f, &self.accession, self.gene_symbol.as_deref())?;
        write!(f, ":p.")?;
        self.fmt_loc_edit(f)
    }
}

/// Mitochondrial variant (m.)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct MtVariant {
    pub accession: Accession,
    pub gene_symbol: Option<String>,
    pub loc_edit: LocEdit<GenomeInterval, NaEdit>,
}

impl MtVariant {
    fn fmt_loc_edit(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.loc_edit)
    }
}

impl fmt::Display for MtVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write_accession_with_optional_gene(f, &self.accession, self.gene_symbol.as_deref())?;
        write!(f, ":m.")?;
        self.fmt_loc_edit(f)
    }
}

/// Circular DNA variant (o.) - SVD-WG006
///
/// Used for circular DNA sequences like plasmids, bacterial chromosomes,
/// and mitochondrial DNA where variants can span the origin.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CircularVariant {
    pub accession: Accession,
    pub gene_symbol: Option<String>,
    pub loc_edit: LocEdit<GenomeInterval, NaEdit>,
}

impl CircularVariant {
    fn fmt_loc_edit(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.loc_edit)
    }
}

impl fmt::Display for CircularVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write_accession_with_optional_gene(f, &self.accession, self.gene_symbol.as_deref())?;
        write!(f, ":o.")?;
        self.fmt_loc_edit(f)
    }
}

/// RNA Fusion breakpoint - SVD-WG007
///
/// Represents one side of an RNA fusion with accession, gene symbol, and interval
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct RnaFusionBreakpoint {
    pub accession: Accession,
    pub gene_symbol: Option<String>,
    pub interval: RnaInterval,
}

impl fmt::Display for RnaFusionBreakpoint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(ref gene) = self.gene_symbol {
            write!(f, "{}({}):r.{}", self.accession, gene, self.interval)
        } else {
            write!(f, "{}:r.{}", self.accession, self.interval)
        }
    }
}

/// RNA Fusion variant (::) - SVD-WG007
///
/// Represents a fusion transcript where two RNA sequences are joined.
/// Common in cancer genomics (BCR-ABL1, EML4-ALK, etc.)
///
/// Format: `accession:r.interval::accession:r.interval`
/// Example: `NM_152263.2:r.-115_775::NM_002609.3:r.1580_*1924`
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct RnaFusionVariant {
    /// 5' partner (upstream)
    pub five_prime: RnaFusionBreakpoint,
    /// 3' partner (downstream)
    pub three_prime: RnaFusionBreakpoint,
}

impl RnaFusionVariant {
    pub fn new(five_prime: RnaFusionBreakpoint, three_prime: RnaFusionBreakpoint) -> Self {
        Self {
            five_prime,
            three_prime,
        }
    }
}

impl fmt::Display for RnaFusionVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}::{}", self.five_prime, self.three_prime)
    }
}

/// Check if a variant causes a frameshift (net indel length not divisible by 3).
///
/// Returns false if the indel length is 0 or cannot be determined.
pub fn is_frameshift(variant: &HgvsVariant) -> bool {
    match crate::python_helpers::get_indel_length(variant) {
        Some(len) if len != 0 => len % 3 != 0,
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::edit::Base;
    use crate::hgvs::location::GenomePos;

    #[test]
    fn test_accession_display() {
        let acc = Accession::new("NM", "000088", Some(3));
        assert_eq!(format!("{}", acc), "NM_000088.3");

        let acc_no_version = Accession::new("NM", "000088", None);
        assert_eq!(format!("{}", acc_no_version), "NM_000088");
    }

    #[test]
    fn test_genome_variant_display() {
        let variant = GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(12345)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        };
        assert_eq!(format!("{}", variant), "NC_000001.11:g.12345A>G");
    }

    #[test]
    fn test_cds_variant_display() {
        use crate::hgvs::location::CdsPos;

        let variant = CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(459)),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        };
        assert_eq!(format!("{}", variant), "NM_000088.3:c.459del");
    }

    #[test]
    fn test_allele_cis_display() {
        use crate::hgvs::location::CdsPos;

        let var1 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let var2 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(200)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let allele = AlleleVariant::cis(vec![var1, var2]);
        let allele_variant = HgvsVariant::Allele(allele);

        assert_eq!(
            format!("{}", allele_variant),
            "NM_000088.3:c.[100A>G;200C>T]"
        );
    }

    #[test]
    fn test_allele_trans_display() {
        use crate::hgvs::location::CdsPos;

        let var1 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let var2 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(200)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let allele = AlleleVariant::trans(vec![var1, var2]);
        let allele_variant = HgvsVariant::Allele(allele);

        assert_eq!(
            format!("{}", allele_variant),
            "NM_000088.3:c.[100A>G];[200C>T]"
        );
    }

    #[test]
    fn test_allele_cis_three_variants() {
        use crate::hgvs::location::CdsPos;

        let make_var = |pos, alt| {
            HgvsVariant::Cds(CdsVariant {
                accession: Accession::new("NM", "000088", Some(3)),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(pos)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: alt,
                    },
                ),
            })
        };

        let allele = AlleleVariant::cis(vec![
            make_var(100, Base::G),
            make_var(200, Base::C),
            make_var(300, Base::T),
        ]);

        assert_eq!(
            format!("{}", HgvsVariant::Allele(allele)),
            "NM_000088.3:c.[100A>G;200A>C;300A>T]"
        );
    }

    #[test]
    fn test_allele_trans_three_variants() {
        use crate::hgvs::location::CdsPos;

        let make_var = |pos, alt| {
            HgvsVariant::Cds(CdsVariant {
                accession: Accession::new("NM", "000088", Some(3)),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(pos)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: alt,
                    },
                ),
            })
        };

        let allele = AlleleVariant::trans(vec![
            make_var(100, Base::G),
            make_var(200, Base::C),
            make_var(300, Base::T),
        ]);

        assert_eq!(
            format!("{}", HgvsVariant::Allele(allele)),
            "NM_000088.3:c.[100A>G];[200A>C];[300A>T]"
        );
    }

    #[test]
    fn test_allele_unknown_phase_compact_form() {
        // 2+ unknown-phase variants sharing accession + coordinate type emit the
        // bracket-less compact form `ACC:c.edit1(;)edit2`.
        use crate::hgvs::location::CdsPos;

        let make_var = |pos, alt| {
            HgvsVariant::Cds(CdsVariant {
                accession: Accession::new("NM", "000088", Some(3)),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(pos)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: alt,
                    },
                ),
            })
        };

        let pair =
            AlleleVariant::unknown_phase(vec![make_var(100, Base::G), make_var(200, Base::C)]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(pair)),
            "NM_000088.3:c.100A>G(;)200A>C"
        );

        let triple = AlleleVariant::unknown_phase(vec![
            make_var(100, Base::G),
            make_var(200, Base::C),
            make_var(300, Base::T),
        ]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(triple)),
            "NM_000088.3:c.100A>G(;)200A>C(;)300A>T"
        );
    }

    #[test]
    fn test_allele_mixed_coordinate_types_uses_expanded_form() {
        // Same accession, different coordinate types (c. and n.) → expanded form.
        // `all_share_accession_and_type` requires both fields to match.
        use crate::hgvs::location::{CdsPos, TxPos};

        let coding = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let noncoding = HgvsVariant::Tx(TxVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(200)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let cis = AlleleVariant::cis(vec![coding.clone(), noncoding.clone()]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(cis)),
            "[NM_000088.3:c.100A>G;NM_000088.3:n.200C>T]"
        );

        let trans = AlleleVariant::trans(vec![coding, noncoding]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(trans)),
            "[NM_000088.3:c.100A>G];[NM_000088.3:n.200C>T]"
        );
    }

    #[test]
    fn test_allele_singleton_emits_bare_form() {
        use crate::hgvs::location::CdsPos;

        let make_var = || {
            HgvsVariant::Cds(CdsVariant {
                accession: Accession::new("NM", "000088", Some(3)),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(100)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::G,
                    },
                ),
            })
        };
        let bare = "NM_000088.3:c.100A>G";

        for allele in [
            AlleleVariant::cis(vec![make_var()]),
            AlleleVariant::trans(vec![make_var()]),
            AlleleVariant::unknown_phase(vec![make_var()]),
            AlleleVariant::mosaic(vec![make_var()]),
            AlleleVariant::chimeric(vec![make_var()]),
        ] {
            assert_eq!(format!("{}", HgvsVariant::Allele(allele)), bare);
        }
    }

    #[test]
    fn test_allele_mixed_accession_uses_expanded_form() {
        use crate::hgvs::location::CdsPos;

        // Different accessions → expanded form (no compact shorthand)
        let var1 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let var2 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000099", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(200)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let cis = AlleleVariant::cis(vec![var1.clone(), var2.clone()]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(cis)),
            "[NM_000088.3:c.100A>G;NM_000099.1:c.200C>T]"
        );

        let trans = AlleleVariant::trans(vec![var1, var2]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(trans)),
            "[NM_000088.3:c.100A>G];[NM_000099.1:c.200C>T]"
        );
    }

    #[test]
    fn test_allele_with_unknown_sub_variant_uses_expanded_form() {
        // Per HGVS spec, bracketed `?` is reserved for the whole-allele marker (`[?]`).
        // An allele containing a `c.?` sub-variant must NOT collapse to `[?;100A>G]` —
        // that would visually conflict with the spec-sanctioned form.
        use crate::hgvs::location::CdsPos;

        let unknown = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::whole_entity_unknown(),
            ),
        });
        let concrete = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let trans = AlleleVariant::trans(vec![concrete, unknown]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(trans)),
            "[NM_000088.3:c.100A>G];[NM_000088.3:c.?]"
        );
    }

    #[test]
    fn test_allele_compact_form_requires_matching_gene_symbol() {
        // Hand-built cis allele: same accession, same coordinate type, but the
        // sub-variants disagree on gene_symbol. The compact form
        // (`ACC(GENE):c.[edit1;edit2]`) can only carry one selector on the
        // prefix, so collapsing two distinct selectors would silently drop
        // the second. Per the #121 round-trip policy ("preserve when present;
        // do not synthesize"), we fall back to the expanded form which
        // emits each sub-variant's selector verbatim.
        use crate::hgvs::location::CdsPos;

        let acc = Accession::new("NM", "000088", Some(3));
        let v_a = HgvsVariant::Cds(CdsVariant {
            accession: acc.clone(),
            gene_symbol: Some("COL1A1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let v_b = HgvsVariant::Cds(CdsVariant {
            accession: acc.clone(),
            gene_symbol: Some("COL1A2".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(200)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let cis = AlleleVariant::cis(vec![v_a.clone(), v_b.clone()]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(cis)),
            "[NM_000088.3(COL1A1):c.100A>G;NM_000088.3(COL1A2):c.200C>T]",
            "mismatched gene symbols must use expanded form to preserve both"
        );

        let trans = AlleleVariant::trans(vec![v_a, v_b]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(trans)),
            "[NM_000088.3(COL1A1):c.100A>G];[NM_000088.3(COL1A2):c.200C>T]",
            "trans expanded form already carries per-variant selectors"
        );
    }

    #[test]
    fn test_allele_compact_form_requires_matching_gene_symbol_some_vs_none() {
        // First sub-variant has Some(gene), second has None. Compact form
        // would either lose the gene symbol (emit bare prefix) or synthesize
        // one for the bare sub-variant (emit prefix-with-gene). Both violate
        // the #121 policy; fall back to expanded form.
        use crate::hgvs::location::CdsPos;

        let acc = Accession::new("NM", "000088", Some(3));
        let v_with = HgvsVariant::Cds(CdsVariant {
            accession: acc.clone(),
            gene_symbol: Some("COL1A1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let v_bare = HgvsVariant::Cds(CdsVariant {
            accession: acc.clone(),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(200)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let cis = AlleleVariant::cis(vec![v_with, v_bare]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(cis)),
            "[NM_000088.3(COL1A1):c.100A>G;NM_000088.3:c.200C>T]",
            "Some/None gene_symbol mismatch must use expanded form"
        );
    }

    #[test]
    fn test_allele_compact_form_keeps_compact_when_all_match() {
        // All sub-variants share Some(same gene): compact form is safe and
        // emits the selector once on the prefix (regression test for the
        // happy path so the new compact-form rule isn't over-constrained).
        use crate::hgvs::location::CdsPos;

        let acc = Accession::new("NM", "000088", Some(3));
        let mk = |pos: i64, alt: Base| {
            HgvsVariant::Cds(CdsVariant {
                accession: acc.clone(),
                gene_symbol: Some("COL1A1".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(pos)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: alt,
                    },
                ),
            })
        };
        let cis = AlleleVariant::cis(vec![mk(100, Base::G), mk(200, Base::T)]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(cis)),
            "NM_000088.3(COL1A1):c.[100A>G;200A>T]",
            "matching gene symbols must keep compact form"
        );

        // All-None: existing compact behavior must be preserved.
        let mk_none = |pos: i64, alt: Base| {
            HgvsVariant::Cds(CdsVariant {
                accession: acc.clone(),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(pos)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: alt,
                    },
                ),
            })
        };
        let cis_none = AlleleVariant::cis(vec![mk_none(100, Base::G), mk_none(200, Base::T)]);
        assert_eq!(
            format!("{}", HgvsVariant::Allele(cis_none)),
            "NM_000088.3:c.[100A>G;200A>T]",
            "all-None gene_symbol must use compact form unchanged"
        );
    }

    #[test]
    fn test_allele_null_variant_singleton_delegates() {
        // NullAllele has no accession; the singleton path must not invoke
        // `use_compact_form` (which would deref the missing accession).
        let null = HgvsVariant::NullAllele;
        let cis = AlleleVariant::cis(vec![null]);
        assert_eq!(format!("{}", HgvsVariant::Allele(cis)), "0");
    }

    #[test]
    fn test_allele_accession() {
        use crate::hgvs::location::CdsPos;

        let var1 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });

        let allele = AlleleVariant::cis(vec![var1]);
        let allele_variant = HgvsVariant::Allele(allele);

        assert_eq!(
            &*allele_variant
                .accession()
                .expect("Expected accession")
                .prefix,
            "NM"
        );
        assert!(allele_variant.is_allele());
    }

    #[test]
    fn test_allele_mosaic_display() {
        use crate::hgvs::location::CdsPos;

        let var1 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let var2 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(200)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let allele = AlleleVariant::mosaic(vec![var1, var2]);
        let allele_variant = HgvsVariant::Allele(allele);

        // Mosaic uses single slash separator: var1/var2
        assert_eq!(
            format!("{}", allele_variant),
            "NM_000088.3:c.100A>G/NM_000088.3:c.200C>T"
        );
    }

    #[test]
    fn test_allele_chimeric_display() {
        use crate::hgvs::location::CdsPos;

        let var1 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let var2 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(200)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let allele = AlleleVariant::chimeric(vec![var1, var2]);
        let allele_variant = HgvsVariant::Allele(allele);

        // Chimerism uses double slash separator: var1//var2
        assert_eq!(
            format!("{}", allele_variant),
            "NM_000088.3:c.100A>G//NM_000088.3:c.200C>T"
        );
    }

    #[test]
    fn test_allele_phase_variants() {
        // Verify all phase types are correctly identified
        assert!(!AllelePhase::Cis.is_mosaic());
        assert!(!AllelePhase::Trans.is_mosaic());
        assert!(!AllelePhase::Unknown.is_mosaic());
        assert!(AllelePhase::Mosaic.is_mosaic());
        assert!(!AllelePhase::Chimeric.is_mosaic());

        assert!(!AllelePhase::Cis.is_chimeric());
        assert!(!AllelePhase::Trans.is_chimeric());
        assert!(!AllelePhase::Unknown.is_chimeric());
        assert!(!AllelePhase::Mosaic.is_chimeric());
        assert!(AllelePhase::Chimeric.is_chimeric());
    }

    // ============== Accession Tests ==============

    #[test]
    fn test_accession_new() {
        let acc = Accession::new("NM", "000088", Some(3));
        assert_eq!(&*acc.prefix, "NM");
        assert_eq!(&*acc.number, "000088");
        assert_eq!(acc.version, Some(3));
        assert!(!acc.ensembl_style);
    }

    #[test]
    fn test_accession_with_style() {
        let acc = Accession::with_style("ENST", "00000012345", Some(1), true);
        assert!(acc.ensembl_style);
        assert_eq!(format!("{}", acc), "ENST00000012345.1");
    }

    #[test]
    fn test_accession_from_assembly() {
        let acc = Accession::from_assembly("GRCh37", "chr1");
        assert!(acc.is_assembly_ref());
        assert_eq!(format!("{}", acc), "GRCh37(chr1)");
    }

    #[test]
    fn test_accession_is_ensembl_prefix() {
        assert!(Accession::is_ensembl_prefix("ENST"));
        assert!(Accession::is_ensembl_prefix("ENSG"));
        assert!(Accession::is_ensembl_prefix("ENSP"));
        assert!(Accession::is_ensembl_prefix("ENSE"));
        assert!(Accession::is_ensembl_prefix("ENSR"));
        assert!(!Accession::is_ensembl_prefix("NM"));
        assert!(!Accession::is_ensembl_prefix("NC"));
    }

    #[test]
    fn test_accession_is_ensembl() {
        let ensembl = Accession::new("ENST", "00000012345", Some(1));
        assert!(ensembl.is_ensembl());

        let refseq = Accession::new("NM", "000088", Some(3));
        assert!(!refseq.is_ensembl());
    }

    #[test]
    fn test_accession_validate_ensembl() {
        // Valid Ensembl ID (11 digits)
        let valid = Accession::new("ENST", "00000012345", Some(1));
        assert!(valid.validate_ensembl());

        // Non-Ensembl IDs should always validate
        let refseq = Accession::new("NM", "000088", Some(3));
        assert!(refseq.validate_ensembl());
    }

    #[test]
    fn test_accession_inferred_variant_type() {
        assert_eq!(
            Accession::new("NC", "000001", Some(11)).inferred_variant_type(),
            Some("g")
        );
        assert_eq!(
            Accession::new("NG", "012345", Some(1)).inferred_variant_type(),
            Some("g")
        );
        assert_eq!(
            Accession::new("NM", "000088", Some(3)).inferred_variant_type(),
            Some("c")
        );
        assert_eq!(
            Accession::new("NR", "123456", Some(1)).inferred_variant_type(),
            Some("n")
        );
        assert_eq!(
            Accession::new("NP", "000079", Some(2)).inferred_variant_type(),
            Some("p")
        );
        assert_eq!(
            Accession::new("ENST", "00000012345", Some(1)).inferred_variant_type(),
            Some("c")
        );
        assert_eq!(
            Accession::new("ENSG", "00000012345", Some(1)).inferred_variant_type(),
            Some("g")
        );
        assert_eq!(
            Accession::new("ENSP", "00000012345", Some(1)).inferred_variant_type(),
            Some("p")
        );
        assert_eq!(
            Accession::new("LRG", "1", None).inferred_variant_type(),
            Some("g")
        );
        // LRG transcript discriminator (`t<M>`) → coding (c.)
        // The transcript discriminator is captured in `number` (e.g. "1t1"),
        // so the dispatch must inspect it rather than only the prefix.
        // See https://www.lrg-sequence.org/faq/ — `t<N>` is the transcript.
        assert_eq!(
            Accession::new("LRG", "1t1", None).inferred_variant_type(),
            Some("c")
        );
        assert_eq!(
            Accession::new("LRG", "292t1", None).inferred_variant_type(),
            Some("c")
        );
        // LRG protein discriminator (`p<M>`) → protein (p.)
        assert_eq!(
            Accession::new("LRG", "1p1", None).inferred_variant_type(),
            Some("p")
        );
        assert_eq!(
            Accession::new("LRG", "673p1", None).inferred_variant_type(),
            Some("p")
        );
        // UniProt single-letter prefix
        assert_eq!(
            Accession::new("P", "12345", None).inferred_variant_type(),
            Some("p")
        );
        // Unknown prefix
        assert_eq!(
            Accession::new("XX", "12345", None).inferred_variant_type(),
            None
        );
    }

    #[test]
    fn test_accession_is_uniprot() {
        // Valid UniProt
        let uniprot = Accession::new("P", "12345", None);
        assert!(uniprot.is_uniprot());

        // Not UniProt - prefix too long
        let not_uniprot = Accession::new("NM", "000088", Some(3));
        assert!(!not_uniprot.is_uniprot());

        // Not UniProt - number wrong length
        let wrong_length = Accession::new("P", "123", None);
        assert!(!wrong_length.is_uniprot());
    }

    #[test]
    fn test_accession_is_lrg_prefix() {
        assert!(Accession::is_lrg_prefix("LRG"));
        assert!(!Accession::is_lrg_prefix("lrg"));
        assert!(!Accession::is_lrg_prefix("LRG_"));
        assert!(!Accession::is_lrg_prefix("NM"));
        assert!(!Accession::is_lrg_prefix(""));
    }

    #[test]
    fn test_accession_is_lrg() {
        assert!(Accession::new("LRG", "1", None).is_lrg());
        assert!(Accession::new("LRG", "1t1", None).is_lrg());
        assert!(Accession::new("LRG", "1p1", None).is_lrg());
        assert!(!Accession::new("NM", "000088", Some(3)).is_lrg());
        assert!(!Accession::new("NC", "000001", Some(11)).is_lrg());
    }

    #[test]
    fn test_accession_base() {
        let acc = Accession::new("NM", "000088", Some(3));
        assert_eq!(acc.base(), "NM_000088");

        let ensembl = Accession::with_style("ENST", "00000012345", Some(1), true);
        assert_eq!(ensembl.base(), "ENST00000012345");

        let assembly = Accession::from_assembly("GRCh38", "chrX");
        assert_eq!(assembly.base(), "GRCh38(chrX)");
    }

    #[test]
    fn test_accession_full() {
        let acc = Accession::new("NM", "000088", Some(3));
        assert_eq!(acc.full(), "NM_000088.3");

        let no_version = Accession::new("NM", "000088", None);
        assert_eq!(no_version.full(), "NM_000088");

        let ensembl = Accession::with_style("ENST", "00000012345", Some(1), true);
        assert_eq!(ensembl.full(), "ENST00000012345.1");
    }

    // ============== LocEdit Tests ==============

    #[test]
    fn test_loc_edit_new() {
        use crate::hgvs::location::GenomePos;

        let loc_edit: LocEdit<GenomeInterval, NaEdit> = LocEdit::new(
            GenomeInterval::point(GenomePos::new(100)),
            NaEdit::Deletion {
                sequence: None,
                length: None,
            },
        );
        assert!(matches!(loc_edit.edit, Mu::Certain(_)));
    }

    #[test]
    fn test_loc_edit_with_uncertainty() {
        use crate::hgvs::location::GenomePos;

        let loc_edit: LocEdit<GenomeInterval, NaEdit> = LocEdit::with_uncertainty(
            GenomeInterval::point(GenomePos::new(100)),
            Mu::Uncertain(NaEdit::Deletion {
                sequence: None,
                length: None,
            }),
        );
        assert!(matches!(loc_edit.edit, Mu::Uncertain(_)));
    }

    // ============== AllelePhase Tests ==============

    #[test]
    fn test_allele_phase_display() {
        assert_eq!(format!("{}", AllelePhase::Cis), "cis");
        assert_eq!(format!("{}", AllelePhase::Trans), "trans");
        assert_eq!(format!("{}", AllelePhase::Unknown), "unknown");
        assert_eq!(format!("{}", AllelePhase::Mosaic), "mosaic");
        assert_eq!(format!("{}", AllelePhase::Chimeric), "chimeric");
    }

    // ============== AlleleVariant Tests ==============

    #[test]
    fn test_allele_variant_unknown_phase() {
        use crate::hgvs::location::CdsPos;

        let var1 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let var2 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(200)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let allele = AlleleVariant::unknown_phase(vec![var1, var2]);
        assert_eq!(allele.phase, AllelePhase::Unknown);
    }

    // ============== HgvsVariant Helper Tests ==============

    #[test]
    fn test_hgvs_variant_accession() {
        let variant = HgvsVariant::Genome(GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(12345)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let acc = variant.accession().expect("Expected accession");
        assert_eq!(&*acc.prefix, "NC");
        assert_eq!(&*acc.number, "000001");
    }

    #[test]
    fn test_hgvs_variant_is_allele() {
        use crate::hgvs::location::CdsPos;

        let var1 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        assert!(!var1.is_allele());

        let allele = AlleleVariant::cis(vec![var1.clone()]);
        let allele_variant = HgvsVariant::Allele(allele);
        assert!(allele_variant.is_allele());
    }

    #[test]
    fn test_hgvs_variant_protein_type() {
        use crate::hgvs::location::{AminoAcid, ProtPos};

        let protein = HgvsVariant::Protein(ProteinVariant {
            accession: Accession::new("NP", "000079", Some(2)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                ProtInterval::point(ProtPos::new(AminoAcid::Val, 600)),
                ProteinEdit::Substitution {
                    reference: AminoAcid::Val,
                    alternative: AminoAcid::Glu,
                },
            ),
        });
        assert_eq!(protein.variant_type(), "p");
        assert!(matches!(protein, HgvsVariant::Protein(_)));

        let genomic = HgvsVariant::Genome(GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(12345)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        assert_eq!(genomic.variant_type(), "g");
        assert!(!matches!(genomic, HgvsVariant::Protein(_)));
    }

    #[test]
    fn test_hgvs_variant_variant_type() {
        let genomic = HgvsVariant::Genome(GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(12345)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        assert_eq!(genomic.variant_type(), "g");

        use crate::hgvs::location::CdsPos;
        let cds = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        assert_eq!(cds.variant_type(), "c");
    }

    #[test]
    fn test_cds_variant_with_gene_symbol() {
        use crate::hgvs::location::CdsPos;

        let variant = CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: Some("COL1A1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        };

        // Per #121: Display preserves the gene-symbol selector when set.
        assert_eq!(variant.gene_symbol, Some("COL1A1".to_string()));
        let display = format!("{}", variant);
        assert_eq!(display, "NM_000088.3(COL1A1):c.100A>G");
    }

    #[test]
    fn test_genome_variant_with_gene_symbol() {
        let variant = GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: Some("BRCA1".to_string()),
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(12345)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        };

        // Per #121: Display preserves the gene-symbol selector when set.
        assert_eq!(variant.gene_symbol, Some("BRCA1".to_string()));
        let display = format!("{}", variant);
        assert_eq!(display, "NC_000001.11(BRCA1):g.12345A>G");
    }

    #[test]
    fn test_transcript_accession_bare() {
        let acc = Accession::new("NM", "004119", Some(3));
        assert_eq!(acc.transcript_accession(), "NM_004119.3");
        assert_eq!(acc.full(), "NM_004119.3");
    }

    #[test]
    fn test_transcript_accession_with_genomic_context() {
        let inner = Accession::new("NM", "004119", Some(3));
        let outer = Accession::new("NC", "000013", Some(11));
        let compound = inner.with_genomic_context(outer);

        // full() includes the compound wrapper
        assert_eq!(compound.full(), "NC_000013.11(NM_004119.3)");
        // transcript_accession() returns just the bare inner accession
        assert_eq!(compound.transcript_accession(), "NM_004119.3");
    }

    #[test]
    fn test_transcript_accession_assembly_ref() {
        let acc = Accession::from_assembly("GRCh38", "chr1");
        assert_eq!(acc.transcript_accession(), "GRCh38(chr1)");
    }

    #[test]
    fn test_transcript_accession_ensembl_with_context() {
        let inner = Accession::with_style("ENST", "00000241453", Some(7), true);
        let outer = Accession::new("NC", "000013", Some(11));
        let compound = inner.with_genomic_context(outer);

        assert_eq!(compound.full(), "NC_000013.11(ENST00000241453.7)");
        assert_eq!(compound.transcript_accession(), "ENST00000241453.7");
    }

    // ----- Issue #115: self-cancelling allele detection (E3006) -----

    use crate::hgvs::parser::variant::parse_variant;

    #[test]
    fn test_self_cancelling_detect_spec_example() {
        // Spec example: c.[762_768del;767_774dup]
        let del = parse_variant("NM_004006.2:c.762_768del").unwrap();
        let dup = parse_variant("NM_004006.2:c.767_774dup").unwrap();
        let pair = AlleleVariant::detect_self_cancelling_pair(&[del, dup]);
        assert!(
            pair.is_some(),
            "spec-example overlapping del+dup must be detected"
        );
        let (dl, du) = pair.unwrap();
        assert_eq!(dl, 0, "first index is the del");
        assert_eq!(du, 1, "second index is the dup");
    }

    #[test]
    fn test_self_cancelling_no_overlap() {
        let del = parse_variant("NM_004006.2:c.100_110del").unwrap();
        let dup = parse_variant("NM_004006.2:c.500_510dup").unwrap();
        assert!(AlleleVariant::detect_self_cancelling_pair(&[del, dup]).is_none());
    }

    #[test]
    fn test_self_cancelling_two_dels() {
        let a = parse_variant("NM_004006.2:c.100_110del").unwrap();
        let b = parse_variant("NM_004006.2:c.105_115del").unwrap();
        assert!(AlleleVariant::detect_self_cancelling_pair(&[a, b]).is_none());
    }

    #[test]
    fn test_self_cancelling_different_accessions() {
        let del = parse_variant("NM_004006.2:c.100_110del").unwrap();
        let dup = parse_variant("NM_000088.3:c.105_115dup").unwrap();
        assert!(AlleleVariant::detect_self_cancelling_pair(&[del, dup]).is_none());
    }

    #[test]
    fn test_self_cancelling_dup_first_then_del() {
        // Order in the allele shouldn't matter; the dup-first variant is
        // still a self-cancelling pair.
        let dup = parse_variant("NM_004006.2:c.767_774dup").unwrap();
        let del = parse_variant("NM_004006.2:c.762_768del").unwrap();
        let pair = AlleleVariant::detect_self_cancelling_pair(&[dup, del]);
        assert!(pair.is_some());
        let (dl, du) = pair.unwrap();
        assert_eq!(dl, 1, "del is at index 1");
        assert_eq!(du, 0, "dup is at index 0");
    }

    #[test]
    fn test_self_cancelling_genomic_overlap() {
        let del = parse_variant("NC_000017.11:g.43044295_43044300del").unwrap();
        let dup = parse_variant("NC_000017.11:g.43044298_43044305dup").unwrap();
        assert!(AlleleVariant::detect_self_cancelling_pair(&[del, dup]).is_some());
    }
}
