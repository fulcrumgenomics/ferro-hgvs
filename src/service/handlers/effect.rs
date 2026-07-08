//! Effect prediction endpoint - predict protein consequences using Sequence Ontology terms

use axum::{extract::State, http::StatusCode, response::Json};
use std::time::Instant;

use crate::hgvs::location::CdsPos;
use crate::service::{
    server::AppState,
    types::{
        EffectRequest, EffectResponse, ErrorResponse, NmdPrediction, ProteinConsequence,
        SequenceEffect, ServiceError,
    },
    validation::validate_hgvs,
};

/// Predict the effect of an HGVS variant
///
/// This endpoint analyzes an HGVS variant and predicts its effect using
/// Sequence Ontology (SO) terms. It can optionally include NMD prediction.
pub async fn predict_effect(
    State(state): State<AppState>,
    Json(request): Json<EffectRequest>,
) -> Result<Json<EffectResponse>, (StatusCode, Json<ErrorResponse>)> {
    let start = Instant::now();

    // Validate input HGVS
    if let Err(validation_error) = validate_hgvs(&request.hgvs) {
        let error = ServiceError::InvalidHgvs(validation_error.to_string());
        return Err((
            StatusCode::from_u16(error.status_code()).unwrap_or(StatusCode::INTERNAL_SERVER_ERROR),
            Json(error.to_response()),
        ));
    }

    // Parse the input HGVS
    let hgvs_str = request.hgvs.clone();
    let parse_result =
        tokio::task::spawn_blocking(move || crate::hgvs::parser::parse_hgvs_lenient(&hgvs_str))
            .await
            .map_err(|e| {
                let error = ServiceError::InternalError(format!("Task error: {}", e));
                (StatusCode::INTERNAL_SERVER_ERROR, Json(error.to_response()))
            })?;

    let elapsed_ms = start.elapsed().as_millis() as u64;

    match parse_result {
        Ok(result) => {
            // Predict effect based on variant type and location
            let (effect, protein_consequence, nmd_prediction) =
                predict_from_variant(&state, &result.result, request.include_nmd);

            Ok(Json(EffectResponse {
                input: request.hgvs,
                effect,
                protein_consequence,
                nmd_prediction,
                error: None,
                processing_time_ms: elapsed_ms,
            }))
        }
        Err(e) => Ok(Json(EffectResponse {
            input: request.hgvs,
            effect: None,
            protein_consequence: None,
            nmd_prediction: None,
            error: Some(format!("Failed to parse input: {}", e)),
            processing_time_ms: elapsed_ms,
        })),
    }
}

/// Predict effect from parsed variant
fn predict_from_variant(
    state: &AppState,
    variant: &crate::hgvs::variant::HgvsVariant,
    include_nmd: bool,
) -> (
    Option<SequenceEffect>,
    Option<ProteinConsequence>,
    Option<NmdPrediction>,
) {
    use crate::hgvs::variant::HgvsVariant;

    match variant {
        HgvsVariant::Cds(v) => {
            // Analyze the edit type to predict effect
            let accession = v.accession.to_string();
            let cds_pos = extract_cds_position(&v.loc_edit);

            // Span length derived from the position interval; required for
            // accurate delins frameshift classification (closes #394 item 1).
            let span_len = span_len_from_cds_interval(&v.loc_edit.location);

            // Get the edit info
            let edit_opt = v.loc_edit.edit.inner();
            let (edit_type, is_frameshift, ref_len, alt_len) = match edit_opt {
                Some(edit) => analyze_na_edit(edit, span_len),
                None => ("unknown", false, 0, 0),
            };

            // Determine if intronic
            let is_intronic = cds_pos.is_intronic();

            // `analyze_na_edit` conservative-skips with `(0, 0, false)`
            // for a deletion / duplication / insertion whose span is
            // undecidable (no explicit sequence/length and no
            // position-interval span — e.g. `c.?_123del`). For those,
            // `is_frameshift = false` does NOT mean "in-frame"; the frame
            // delta is simply unknown. Emitting `inframe_deletion` /
            // `inframe_insertion` / `duplication` would assert an
            // in-frame call we can't make, so surface a neutral
            // `coding_sequence_variant` and skip the protein / NMD
            // predictions that consume the (unreliable) frameshift signal.
            //
            // Scoped to del/dup/ins only: `delins` already reports the
            // length-independent `indel`, and `inversion` is HIGH-impact
            // regardless of span — neither makes an inframe claim, so
            // neither needs the guard.
            let coding = !is_intronic && !cds_pos.utr3 && cds_pos.base > 0;
            let span_undecidable = matches!(edit_type, "deletion" | "duplication" | "insertion")
                && ref_len == 0
                && alt_len == 0;

            // Predict SO effect
            let effect = if coding && span_undecidable {
                SequenceEffect {
                    so_term: "SO:0001580".to_string(),
                    name: "coding_sequence_variant".to_string(),
                    description: "A coding-sequence change whose span/length is undecidable"
                        .to_string(),
                    impact: "MODERATE".to_string(),
                }
            } else {
                predict_cds_effect(edit_type, is_intronic, is_frameshift, &cds_pos)
            };

            // Try to get protein consequence if we have cdot data.
            // Pass `is_frameshift` through from `analyze_na_edit` rather
            // than letting `predict_protein_consequence` recompute it
            // from `(alt_len - ref_len)`. The classifier is the
            // authoritative source — it returns `false` for shapes
            // where the net delta can't be determined (intronic /
            // mixed-axis / non-literal insert with unknown length),
            // and `ref_len`/`alt_len` are reported as `0` in those
            // cases. A naive recompute here would mis-flag those as
            // frameshift whenever `alt_len % 3 != 0` (e.g. an
            // unknown-span delins with a single-nt literal insert).
            let protein_consequence = if coding && !span_undecidable {
                predict_protein_consequence(
                    state,
                    &accession,
                    &cds_pos,
                    edit_opt,
                    edit_type,
                    is_frameshift,
                    ref_len,
                )
            } else {
                None
            };

            // A delins whose net frame delta is undeterminable reports
            // `is_frameshift = false` conservatively (not "proven in-frame").
            // The effect/protein paths above are frame-agnostic for delins (they
            // report the length-independent `indel`), but NMD must NOT consume
            // this signal — it would emit a definitive `predicted: false` we
            // cannot justify. `predict_nmd_for_cds` declines (returns `None`)
            // when this flag is set.
            let frame_undecidable =
                edit_opt.is_some_and(|edit| delins_frame_undecidable(edit, span_len));

            // NMD prediction if requested. Skip when the span is
            // undecidable — an NMD call rests on the frameshift signal,
            // which is unknown for a conservative-skip del/dup/ins.
            let nmd_prediction = if include_nmd && !span_undecidable {
                predict_nmd_for_cds(
                    state,
                    &accession,
                    &cds_pos,
                    edit_opt,
                    is_frameshift,
                    edit_type,
                    frame_undecidable,
                )
            } else {
                None
            };

            (Some(effect), protein_consequence, nmd_prediction)
        }
        HgvsVariant::Genome(v) => {
            if let Some(edit) = v.loc_edit.edit.inner() {
                // span_len matters only when downstream consumers
                // (protein-consequence / NMD) use the frameshift signal;
                // the genome path keeps only `edit_type`, so passing
                // `None` skips the unused interval-length computation.
                let (edit_type, _, _, _) = analyze_na_edit(edit, None);
                let effect = predict_genomic_effect(edit_type);
                (Some(effect), None, None)
            } else {
                (None, None, None)
            }
        }
        HgvsVariant::Protein(v) => {
            // For protein variants, we can directly describe the effect
            let effect = SequenceEffect {
                so_term: "SO:0001583".to_string(),
                name: "missense_variant".to_string(),
                description: "A sequence variant that changes one or more bases".to_string(),
                impact: "MODERATE".to_string(),
            };
            let protein = ProteinConsequence {
                hgvs_p: format!("{}", v),
                ref_aa: "".to_string(),
                alt_aa: "".to_string(),
                position: 0,
                is_frameshift: false,
            };
            (Some(effect), Some(protein), None)
        }
        HgvsVariant::Tx(v) => {
            if let Some(edit) = v.loc_edit.edit.inner() {
                // Same rationale as the Genome arm above — the Tx path
                // keeps only `edit_type` for the effect description.
                let (edit_type, _, _, _) = analyze_na_edit(edit, None);
                let effect = SequenceEffect {
                    so_term: "SO:0001619".to_string(),
                    name: "non_coding_transcript_variant".to_string(),
                    description: format!("A {} in a non-coding transcript", edit_type),
                    impact: "MODIFIER".to_string(),
                };
                (Some(effect), None, None)
            } else {
                (None, None, None)
            }
        }
        _ => (None, None, None),
    }
}

/// Extract CdsPos from a LocEdit
fn extract_cds_position(
    loc_edit: &crate::hgvs::variant::LocEdit<
        crate::hgvs::interval::CdsInterval,
        crate::hgvs::edit::NaEdit,
    >,
) -> CdsPos {
    let interval = &loc_edit.location;
    if let Some(start) = interval.start.inner() {
        CdsPos {
            base: start.base,
            offset: start.offset,
            utr3: start.utr3,
            special: start.special,
        }
    } else {
        CdsPos::new(1)
    }
}

/// Analyze nucleic acid edit to determine type and properties.
///
/// Returns `(edit_type, is_frameshift, ref_len, alt_len)`. The
/// `span_len` parameter is the length of the position interval (`end -
/// start + 1`) when computable from base coordinates; pass `None` when
/// the interval has intronic offsets, decorated positions (pter/qter/
/// cen), or unknown endpoints. Used by the `Delins`, `Deletion`, and
/// `Duplication` arms to derive `ref_len` (and `is_frameshift` =
/// `net_delta % 3 != 0`) per the HGVS frameshift spec
/// (`recommendations/protein/frameshift.md`).
///
/// **Conservative-skip contract:** when neither an explicit `sequence`
/// / `length` nor a computable `span_len` is available, `is_frameshift`
/// is reported as `false` and `ref_len` / `alt_len` are reported as
/// `0`. This is **undecidable, not proven in-frame** — downstream
/// callers must not treat `is_frameshift = false` as a positive
/// declaration of in-frame status when the variant's lengths are
/// genuinely unknown. The same contract applies to the `Insertion`
/// arm when `InsertedSequence::len()` returns `None`. (`Substitution`
/// is always 1↔1 by construction, so its `is_frameshift = false` IS a
/// positive declaration.)
///
/// Per-arm `alt_len` semantics under the conservative-skip branch:
/// - `Deletion`: `0` (semantically correct — nothing inserted).
/// - `Duplication`: `0` (placeholder — a real duplication inserts at
///   least one copy, so this `0` is a sentinel for "unknown",
///   matching the `ref_len = 0` shape, NOT a claim that the alt is
///   empty).
/// - `Delins`: `0` if `InsertedSequence::len()` is also `None`.
pub fn analyze_na_edit(
    edit: &crate::hgvs::edit::NaEdit,
    span_len: Option<usize>,
) -> (&'static str, bool, usize, usize) {
    use crate::hgvs::edit::NaEdit;

    match edit {
        NaEdit::Substitution {
            reference,
            alternative,
        } => (
            "substitution",
            false,
            reference.to_string().len(),
            alternative.to_string().len(),
        ),
        NaEdit::SubstitutionNoRef { alternative } => {
            ("substitution", false, 1, alternative.to_string().len())
        }
        NaEdit::Deletion { sequence, length } => {
            // Closes #427 (follow-up to #394 item 1). The previous
            // `unwrap_or(1)` fallback over-predicted frameshift on the
            // canonical short-form `c.100_102del` (where `sequence:
            // None, length: None` is set by the parser and the actual
            // span is 3 bp, in-frame). Mirror the Delins arm: prefer
            // the explicit `sequence` / `length`, fall back to the
            // position-interval `span_len`, and conservative-skip
            // (`is_frameshift = false`, `ref_len = 0`) when none of
            // those is available — better than emitting a frameshift
            // on every short-form deletion.
            let len_opt = sequence
                .as_ref()
                .map(|s| s.to_string().len())
                .or(length.map(|l| l as usize))
                .or(span_len);
            match len_opt {
                Some(len) => ("deletion", len % 3 != 0, len, 0),
                None => ("deletion", false, 0, 0),
            }
        }
        NaEdit::Insertion { sequence } => {
            // `InsertedSequence::len()` is the spec-aligned primitive:
            // `Some(n)` only when the inserted nt count is known
            // (`Literal`, `Count`, `PositionRange`, `Repeat` with
            // `Exact` count, ...). For non-literal shapes whose nt
            // count is unknown (`Range`, `Reference`, `Named`,
            // `Uncertain`, ...), it returns `None`. The previous
            // `to_string().len()` was the *rendered notation* length,
            // which is unrelated to nt count for those variants
            // (e.g. `ins[NC_...]` rendered to ~25 chars). Conservative
            // fallback `is_frameshift = false` when length is unknown.
            match sequence.len() {
                Some(len) => ("insertion", len % 3 != 0, 0, len),
                None => ("insertion", false, 0, 0),
            }
        }
        NaEdit::Delins { sequence, .. } => {
            // Closes #394 item 1. Previously hardcoded `is_frameshift =
            // true` with "Can't determine ref length easily"; the
            // information IS available via `span_len` (ref) and
            // `InsertedSequence::len()` (alt). Conservative fallback
            // (`is_frameshift = false`) when either ref_len or alt_len
            // is unknown — better than the previous always-true default
            // which over-predicted frameshift on every non-trivial
            // delins, and better than the intermediate fix that used
            // `to_string().len()` (rendered notation length, not nt
            // count) for non-literal inserted sequences.
            let alt_len_opt = sequence.len();
            let (ref_len, alt_len, is_frameshift) = match (span_len, alt_len_opt) {
                (Some(rl), Some(al)) => {
                    let net_delta = al as isize - rl as isize;
                    (rl, al, net_delta.rem_euclid(3) != 0)
                }
                // Either endpoint of the (ref, alt) net-delta is
                // unknown — frameshift undecidable; conservative false.
                (Some(rl), None) => (rl, 0, false),
                (None, Some(al)) => (0, al, false),
                (None, None) => (0, 0, false),
            };
            ("delins", is_frameshift, ref_len, alt_len)
        }
        NaEdit::Duplication {
            sequence, length, ..
        } => {
            // Closes #427 (follow-up to #394 item 1). Same fallback
            // chain as the Deletion arm above: `sequence` → `length` →
            // `span_len` → conservative-skip. The `alt_len = len * 2`
            // pattern (one copy of the ref bases is inserted before
            // the duplicated span) is preserved.
            let len_opt = sequence
                .as_ref()
                .map(|s| s.to_string().len())
                .or(length.map(|l| l as usize))
                .or(span_len);
            match len_opt {
                Some(len) => ("duplication", len % 3 != 0, len, len * 2),
                None => ("duplication", false, 0, 0),
            }
        }
        NaEdit::Inversion { sequence, length } => {
            // Closes #438. Per
            // `assets/hgvs-nomenclature/docs/recommendations/DNA/inversion.md`,
            // an inversion reverse-complements its declared range in
            // place. The DNA-level net length change is zero by
            // construction (`alt_len == ref_len`), so the frameshift
            // check `(alt_len - ref_len) % 3 != 0` is always false.
            // The previous unconditional `is_frameshift = true`
            // over-predicted frameshift on every inversion and fed
            // the wrong signal into `predict_protein_consequence` /
            // NMD logic. Inversion may still disrupt protein function
            // (the inverted bases code for different amino acids),
            // but that's a missense / in-frame change, not a frameshift.
            //
            // Fallback chain mirrors the post-#427 Deletion/Duplication
            // arms: prefer explicit `sequence` (e.g. `c.100_102invATG`,
            // 3 bp) over explicit `length` (`c.100_102inv3`) over the
            // position-interval `span_len`. Conservative-skip
            // (`ref_len = alt_len = 0`) when none of those is available.
            let len_opt = sequence
                .as_ref()
                .map(|s| s.to_string().len())
                .or(length.map(|l| l as usize))
                .or(span_len);
            match len_opt {
                Some(len) => ("inversion", false, len, len),
                None => ("inversion", false, 0, 0),
            }
        }
        NaEdit::Repeat { .. } => ("repeat", false, 0, 0),
        NaEdit::Identity { .. } => ("identity", false, 0, 0),
        NaEdit::Unknown { .. } => ("unknown", false, 0, 0),
        _ => ("other", false, 0, 0),
    }
}

/// Whether a `delins`'s net frame delta is **undeterminable** — the deleted
/// span length (from the position interval) or the inserted nucleotide count is
/// unknown. In that state [`analyze_na_edit`] reports `is_frameshift = false`
/// *conservatively* (not "proven in-frame"), so a frame-dependent consumer must
/// decline rather than assert a definitive result.
///
/// This mirrors the `(span_len, alt_len)` match in [`analyze_na_edit`]'s
/// `Delins` arm: the net delta is undecidable exactly when either endpoint of
/// the `(ref, alt)` length pair is `None`. It is consumed only by the NMD path
/// — an undecidable delins must yield *no* NMD prediction rather than a
/// definitive `predicted: false` (issue #806 review). The effect/protein paths
/// are frame-agnostic for delins (they report the length-independent `indel`),
/// so they are unaffected.
///
/// Non-`delins` edits return `false`: the undecidable del/dup/ins case is the
/// `span_undecidable` skip applied separately by the caller, and every other
/// edit either carries a determinable frame delta or is not frame-relevant.
fn delins_frame_undecidable(edit: &crate::hgvs::edit::NaEdit, span_len: Option<usize>) -> bool {
    match edit {
        crate::hgvs::edit::NaEdit::Delins { sequence, .. } => {
            span_len.is_none() || sequence.len().is_none()
        }
        _ => false,
    }
}

/// Compute the position-interval span length for a c. (CDS) interval.
///
/// Returns `Some(end.base - start.base + 1)` only when both endpoints are
/// present, integer-only (no offset), non-unknown (`c.?` placeholder),
/// and on the same coordinate sub-axis (both 5'UTR, both CDS, or both
/// 3'UTR).
///
/// Mixed `c.-N_M` style ranges (5'UTR-to-CDS or CDS-to-3'UTR) have no
/// well-defined integer span in tx-frame from the c. values alone
/// because HGVS skips `c.0`: `c.-1_1` covers 2 bases, not the naive
/// `1 - (-1) + 1 = 3`. The same-sign + same-utr3 guards collapse those
/// shapes to `None` per the conservative-skip contract.
///
/// Returns `None` for any case where the span can't be computed from
/// base coordinates alone — the conservative-skip contract used by the
/// effect classifier for delins frameshift inference.
pub fn span_len_from_cds_interval(interval: &crate::hgvs::interval::CdsInterval) -> Option<usize> {
    let start = interval.start.inner()?;
    let end = interval.end.inner()?;
    // Unknown placeholder (`c.?`): `CdsPos::is_unknown()` keys off
    // `CDS_BASE_UNKNOWN == 0 && !utr3`. The parser rejects `c.0` (#269)
    // but programmatic construction can produce this shape, so guard
    // explicitly to avoid emitting `Some(1)` for a placeholder.
    // Special positions (pter/qter/cen) also carry `base == 0` but with
    // `special.is_some()`; they must be skipped too — mirroring the
    // `span_len_from_genome_interval` guard below.
    if start.is_unknown() || start.is_special() || end.is_unknown() || end.is_special() {
        return None;
    }
    // Intronic-offset positions have no base-only span.
    if start.offset.is_some() || end.offset.is_some() {
        return None;
    }
    // Mixed axis sides (5'UTR ↔ CDS, CDS ↔ 3'UTR, etc.) — base values
    // are in different numbering frames, so `end - start + 1` does not
    // count tx-frame nucleotides.
    if start.utr3 != end.utr3 {
        return None;
    }
    // 5'UTR vs CDS distinction (signed base): same `utr3` flag but
    // negative-vs-positive bases also cross the CDS-start boundary,
    // and HGVS skips `c.0` (`c.-1_1` spans 2 bases, not 3).
    if (start.base < 0) != (end.base < 0) {
        return None;
    }
    if end.base < start.base {
        return None;
    }
    Some(((end.base - start.base) + 1) as usize)
}

/// Compute the position-interval span length for a g. (genome) interval.
///
/// Returns `Some(end.base - start.base + 1)` only when both endpoints are
/// present, integer-only, and non-special. See
/// [`span_len_from_cds_interval`] for the conservative-skip contract.
pub fn span_len_from_genome_interval(
    interval: &crate::hgvs::interval::GenomeInterval,
) -> Option<usize> {
    let start = interval.start.inner()?;
    let end = interval.end.inner()?;
    if start.offset.is_some() || end.offset.is_some() {
        return None;
    }
    if start.is_special() || end.is_special() {
        return None;
    }
    if end.base < start.base {
        return None;
    }
    Some(((end.base - start.base) + 1) as usize)
}

/// Compute the position-interval span length for an n. (tx) interval.
///
/// Returns `Some(end.base - start.base + 1)` only when both endpoints are
/// present, integer-only (no offset), and non-downstream (no `n.*N`).
/// See [`span_len_from_cds_interval`] for the conservative-skip
/// contract.
pub fn span_len_from_tx_interval(interval: &crate::hgvs::interval::TxInterval) -> Option<usize> {
    let start = interval.start.inner()?;
    let end = interval.end.inner()?;
    if start.offset.is_some() || end.offset.is_some() {
        return None;
    }
    // Per the helper's contract, any downstream (`n.*N`) endpoint is a
    // conservative skip — both pure-downstream spans (`n.*5_*10`) and
    // mixed downstream/non-downstream spans cross axes that the
    // `end.base - start.base + 1` arithmetic does not count correctly.
    if start.downstream || end.downstream {
        return None;
    }
    // The parser rejects `n.0` (`src/hgvs/parser/position.rs:268-275`)
    // but programmatic construction can produce that shape; guard
    // explicitly for symmetry with `span_len_from_cds_interval`'s
    // `is_unknown()` guard.
    if start.base == 0 || end.base == 0 {
        return None;
    }
    if (start.base < 0) != (end.base < 0) {
        return None;
    }
    if end.base < start.base {
        return None;
    }
    Some(((end.base - start.base) + 1) as usize)
}

/// Predict effect for CDS variant
fn predict_cds_effect(
    edit_type: &str,
    is_intronic: bool,
    is_frameshift: bool,
    cds_pos: &CdsPos,
) -> SequenceEffect {
    if is_intronic {
        // Check for splice site impact
        if let Some(offset) = cds_pos.offset {
            if offset.abs() <= 2 {
                return SequenceEffect {
                    so_term: "SO:0001629".to_string(),
                    name: "splice_site_variant".to_string(),
                    description: "A sequence variant that changes the first two or last two bases of an intron".to_string(),
                    impact: "HIGH".to_string(),
                };
            } else if offset.abs() <= 8 {
                return SequenceEffect {
                    so_term: "SO:0001630".to_string(),
                    name: "splice_region_variant".to_string(),
                    description: "A sequence variant in which a change has occurred within the region of the splice site".to_string(),
                    impact: "LOW".to_string(),
                };
            }
        }
        return SequenceEffect {
            so_term: "SO:0001627".to_string(),
            name: "intron_variant".to_string(),
            description: "A transcript variant occurring within an intron".to_string(),
            impact: "MODIFIER".to_string(),
        };
    }

    // 5' UTR
    if cds_pos.base < 0 {
        return SequenceEffect {
            so_term: "SO:0001623".to_string(),
            name: "5_prime_UTR_variant".to_string(),
            description: "A UTR variant of the 5' UTR".to_string(),
            impact: "MODIFIER".to_string(),
        };
    }

    // 3' UTR
    if cds_pos.utr3 {
        return SequenceEffect {
            so_term: "SO:0001624".to_string(),
            name: "3_prime_UTR_variant".to_string(),
            description: "A UTR variant of the 3' UTR".to_string(),
            impact: "MODIFIER".to_string(),
        };
    }

    // Within CDS
    match edit_type {
        "substitution" => SequenceEffect {
            so_term: "SO:0001583".to_string(),
            name: "missense_variant".to_string(),
            description: "A codon change resulting in a different amino acid (predicted)"
                .to_string(),
            impact: "MODERATE".to_string(),
        },
        "deletion" => {
            if is_frameshift {
                SequenceEffect {
                    so_term: "SO:0001589".to_string(),
                    name: "frameshift_variant".to_string(),
                    description:
                        "A sequence variant which causes a disruption of the translational reading frame"
                            .to_string(),
                    impact: "HIGH".to_string(),
                }
            } else {
                SequenceEffect {
                    so_term: "SO:0001822".to_string(),
                    name: "inframe_deletion".to_string(),
                    description:
                        "An inframe non-synonymous variant that deletes bases from the coding sequence"
                            .to_string(),
                    impact: "MODERATE".to_string(),
                }
            }
        }
        "insertion" => {
            if is_frameshift {
                SequenceEffect {
                    so_term: "SO:0001589".to_string(),
                    name: "frameshift_variant".to_string(),
                    description:
                        "A sequence variant which causes a disruption of the translational reading frame"
                            .to_string(),
                    impact: "HIGH".to_string(),
                }
            } else {
                SequenceEffect {
                    so_term: "SO:0001821".to_string(),
                    name: "inframe_insertion".to_string(),
                    description:
                        "An inframe non-synonymous variant that inserts bases into the coding sequence"
                            .to_string(),
                    impact: "MODERATE".to_string(),
                }
            }
        }
        "delins" => SequenceEffect {
            so_term: "SO:1000032".to_string(),
            name: "indel".to_string(),
            description: "A sequence alteration which includes both deletion and insertion"
                .to_string(),
            impact: "MODERATE".to_string(),
        },
        "duplication" => {
            if is_frameshift {
                SequenceEffect {
                    so_term: "SO:0001589".to_string(),
                    name: "frameshift_variant".to_string(),
                    description:
                        "A sequence variant which causes a disruption of the translational reading frame"
                            .to_string(),
                    impact: "HIGH".to_string(),
                }
            } else {
                SequenceEffect {
                    so_term: "SO:1000035".to_string(),
                    name: "duplication".to_string(),
                    description:
                        "An insertion which derives from a copy of a sequence immediately adjacent"
                            .to_string(),
                    impact: "MODERATE".to_string(),
                }
            }
        }
        "inversion" => SequenceEffect {
            so_term: "SO:1000036".to_string(),
            name: "inversion".to_string(),
            description: "A continuous nucleotide sequence is inverted in the same position"
                .to_string(),
            impact: "HIGH".to_string(),
        },
        _ => SequenceEffect {
            so_term: "SO:0001580".to_string(),
            name: "coding_sequence_variant".to_string(),
            description: "A sequence variant that changes the coding sequence".to_string(),
            impact: "MODERATE".to_string(),
        },
    }
}

/// Predict effect for genomic variant
fn predict_genomic_effect(edit_type: &str) -> SequenceEffect {
    match edit_type {
        "substitution" => SequenceEffect {
            so_term: "SO:0001483".to_string(),
            name: "SNV".to_string(),
            description: "A single nucleotide variant".to_string(),
            impact: "MODIFIER".to_string(),
        },
        "deletion" => SequenceEffect {
            so_term: "SO:0000159".to_string(),
            name: "deletion".to_string(),
            description: "A sequence alteration where nucleotides are removed".to_string(),
            impact: "MODIFIER".to_string(),
        },
        "insertion" => SequenceEffect {
            so_term: "SO:0000667".to_string(),
            name: "insertion".to_string(),
            description: "The sequence of one or more nucleotides added".to_string(),
            impact: "MODIFIER".to_string(),
        },
        _ => SequenceEffect {
            so_term: "SO:0001060".to_string(),
            name: "sequence_variant".to_string(),
            description: "A sequence variant".to_string(),
            impact: "MODIFIER".to_string(),
        },
    }
}

/// Predict protein consequence using transcript data.
///
/// `is_frameshift` must be the value returned by [`analyze_na_edit`] —
/// the classifier already accounts for unknown ref/alt lengths and
/// reports `false` for those cases. Recomputing from `alt_len - ref_len`
/// here would silently overturn that decision whenever
/// `alt_len.rem_euclid(3) != 0` (e.g. unknown-span delins where
/// `ref_len = 0` is a placeholder, not a true zero-length ref).
fn predict_protein_consequence(
    state: &AppState,
    accession: &str,
    cds_pos: &CdsPos,
    edit: Option<&crate::hgvs::edit::NaEdit>,
    edit_type: &str,
    is_frameshift: bool,
    ref_len: usize,
) -> Option<ProteinConsequence> {
    // Need cdot data for protein prediction
    let cdot = state.cdot.as_ref()?;

    // Get transcript (verify it exists)
    let cdot_tx = cdot.get_transcript(accession)?;

    // Calculate protein position: (cds_pos - 1) / 3 + 1
    let prot_position = if cds_pos.base > 0 {
        ((cds_pos.base - 1) / 3 + 1) as u64
    } else {
        return None; // UTR position
    };

    // Calculate codon phase (position within codon: 0, 1, or 2)
    let codon_phase = ((cds_pos.base - 1) % 3) as u8;

    // Get the protein accession: the authoritative cdot value if present, else
    // the transcript accession itself. We compute it here (before resolving
    // residues) and thread it through `resolve_residues*` so the sequence-backed
    // path uses the SAME accession as the no-reference fallback below. The
    // sequence-bearing reference `Transcript.protein_id` can be absent even when
    // cdot carries the protein accession, so resolving from the transcript there
    // would flip `NP_*:p...` to `NM_*:p...` purely by enabling `state.reference`
    // (#806 review). We do NOT infer `NP_*`/`XP_*` from `NM_*`/`XM_*` by
    // preserving the number — RefSeq does not guarantee the NM and NP numbers
    // match, so that inference is frequently wrong (#808).
    let prot_acc = cdot_tx
        .protein
        .clone()
        .unwrap_or_else(|| accession.to_string());

    // Try to resolve real amino-acid residues from a sequence-bearing
    // reference (issue #806). When `state.reference` is available we always
    // resolve the *reference* residue at the variant's first CDS codon (cheap,
    // all edit classes); for substitutions we additionally resolve the real
    // alternate residue and full HGVS-p notation via the engine's
    // `predict_substitution_protein`. See `resolve_residues` for the honest
    // `?` contract when sequence is unavailable or the edit class is deferred
    // to #498.
    if let Some(resolved) = resolve_residues(
        state,
        accession,
        &prot_acc,
        cds_pos,
        edit,
        edit_type,
        is_frameshift,
        ref_len,
    ) {
        return Some(resolved);
    }

    // `is_frameshift` is authoritative from the classifier; the
    // `edit_type != "substitution"` guard is belt-and-braces (the
    // classifier already returns `false` for the substitution arms).
    let is_frameshift = is_frameshift && edit_type != "substitution";

    // Build HGVS protein notation
    // Without sequence data, we use position-based notation with uncertainty markers
    let hgvs_p = if is_frameshift {
        // Frameshift: p.(Xxx123fs)
        format!("{}:p.(?{}fs)", prot_acc, prot_position)
    } else if edit_type == "substitution" {
        // Substitution: p.(Xxx123?) - unknown AA change at position
        format!("{}:p.(?{}?)", prot_acc, prot_position)
    } else if edit_type == "deletion" {
        // In-frame deletion
        let end_pos = prot_position + (ref_len / 3).max(1) as u64 - 1;
        if prot_position == end_pos {
            format!("{}:p.(?{}del)", prot_acc, prot_position)
        } else {
            format!("{}:p.(?{}_?{}del)", prot_acc, prot_position, end_pos)
        }
    } else if edit_type == "insertion" {
        // In-frame insertion
        format!(
            "{}:p.(?{}_?{}ins?)",
            prot_acc,
            prot_position,
            prot_position + 1
        )
    } else {
        // Other: just show position
        format!("{}:p.(?{}?)", prot_acc, prot_position)
    };

    // Honest fallback (issue #806, state 2 = "data unavailable"). We reach
    // this branch only when `resolve_residues` returned `None` — i.e. no
    // sequence-bearing reference is configured (`state.reference` is `None`) or
    // the transcript carries no CDS sequence. The `?` markers therefore mean
    // "the service cannot obtain CDS bases", NOT "not implemented": when a
    // reference IS available, `resolve_residues` always returns real residues
    // (and never reaches here). The position-based HGVS-p notation above is the
    // best the service can emit without bases.
    Some(ProteinConsequence {
        hgvs_p,
        ref_aa: format!("?(pos{})", codon_phase + 1), // Show codon position (1-3)
        alt_aa: "?".to_string(),
        position: prot_position,
        is_frameshift,
    })
}

/// Resolve real amino-acid residues for a CDS variant using a sequence-bearing
/// reference provider (issue #806).
///
/// Returns `None` when no real sequence is available — `state.reference` is
/// unconfigured, the transcript is missing or carries no sequence, or the CDS
/// position is not a resolvable in-CDS position — so the caller can fall back to
/// the honest "data unavailable" (`?`) signal.
///
/// When sequence IS available it returns a [`ProteinConsequence`] with the real
/// **reference** residue at the variant's first CDS codon for every edit class.
/// For substitutions it additionally resolves the real **alternate** residue
/// and full HGVS-p notation via the engine's
/// [`predict_substitution_protein`](crate::project::protein::predict_substitution_protein),
/// inheriting its initiation-codon (`p.(Met1?)`) and stop-loss-extension
/// handling. For non-substitution indels the alternate residue and full
/// protein-level notation are deferred to #498, so `alt_aa` carries the `?`
/// sentinel meaning "not yet implemented for this edit class" — distinct from
/// the no-reference "data unavailable" `?`.
///
/// This helper is the shared seam #498 (full c.→p.) inherits: it already fetches
/// a sequence-bearing [`Transcript`](crate::reference::transcript::Transcript)
/// and calls into the same protein machinery the engine uses.
// Threads both the transcript `accession` (to fetch sequence from the provider)
// and the authoritative protein `prot_acc` (for HGVS-p), plus the edit context;
// splitting these into a struct would obscure the simple pass-through.
#[allow(clippy::too_many_arguments)]
fn resolve_residues(
    state: &AppState,
    accession: &str,
    prot_acc: &str,
    cds_pos: &CdsPos,
    edit: Option<&crate::hgvs::edit::NaEdit>,
    edit_type: &str,
    is_frameshift: bool,
    ref_len: usize,
) -> Option<ProteinConsequence> {
    // A reference residue is only defined for an in-CDS position.
    if cds_pos.base <= 0 || cds_pos.utr3 || cds_pos.offset.is_some() {
        return None;
    }
    let provider = state.reference.as_ref()?;
    let transcript =
        crate::reference::provider::ReferenceProvider::get_transcript(provider.as_ref(), accession)
            .ok()?;
    resolve_residues_from_transcript(
        &transcript,
        prot_acc,
        cds_pos,
        edit,
        edit_type,
        is_frameshift,
        ref_len,
    )
}

/// Pure residue-resolution core, separated from `AppState` so it is unit
/// testable against a synthetic sequence-bearing
/// [`Transcript`](crate::reference::transcript::Transcript) (issue #806).
///
/// See [`resolve_residues`] for the honest `?` contract. Returns `None` when the
/// reference residue cannot be read from `transcript.sequence` (e.g. the
/// transcript carries no sequence) so the caller falls back to the
/// "data unavailable" signal.
fn resolve_residues_from_transcript(
    transcript: &crate::reference::transcript::Transcript,
    prot_acc: &str,
    cds_pos: &CdsPos,
    edit: Option<&crate::hgvs::edit::NaEdit>,
    edit_type: &str,
    is_frameshift: bool,
    ref_len: usize,
) -> Option<ProteinConsequence> {
    use crate::project::protein::{predict_substitution_protein, read_ref_codon, translate};

    let cds_base = cds_pos.base; // i64, 1-based CDS coordinate
    let prot_position = ((cds_base - 1) / 3 + 1) as u64;

    // Reference residue at the variant's first CDS codon — resolved for every
    // edit class (this is the "no silent ? where data exists" guarantee).
    let (ref_codon, frame) = read_ref_codon(transcript, cds_base).ok()?;
    let ref_aa = translate(&ref_codon)?;

    // Substitutions get full alt-residue + HGVS-p resolution via the engine.
    // This arm is *terminal* for every supported substitution form: once we
    // have sequence we either resolve the real alternate residue or emit an
    // EXPLICIT unresolvable `p.(?N?)` state — we never fall through to the
    // indel formatter, which would silently relabel a substitution as a
    // deferred indel (#806 review).
    if edit_type == "substitution" {
        // Normalize the supported substitution forms to a `Substitution` with a
        // concrete reference base. `SubstitutionNoRef` (`c.N>Alt`) carries only
        // the alternate allele, so derive the reference base from the sequence
        // codon at this position (the engine requires a `Substitution`); this
        // lets the no-ref form resolve the real alternate residue when sequence
        // exists rather than degrading to `p.(?N?)`.
        let normalized_sub = match edit {
            Some(sub @ crate::hgvs::edit::NaEdit::Substitution { .. }) => Some(sub.clone()),
            Some(crate::hgvs::edit::NaEdit::SubstitutionNoRef { alternative }) => {
                let ref_base = crate::hgvs::edit::Base::from_char(
                    ref_codon.as_bytes()[frame as usize] as char,
                )?;
                Some(crate::hgvs::edit::NaEdit::Substitution {
                    reference: ref_base,
                    alternative: *alternative,
                })
            }
            _ => None,
        };
        if let Some(sub) = normalized_sub.as_ref() {
            match predict_substitution_protein(transcript, cds_base, sub, prot_acc) {
                Ok(variant) => {
                    // The residue fields must agree with the engine's `hgvs_p`.
                    // For a substitution affecting the translation initiation
                    // codon the engine reports the uncertain initiator form
                    // (`p.(Met1?)`), where the protein consequence — and thus
                    // the alternate residue — is unpredictable and residue 1 is
                    // the initiator Met regardless of the reference codon.
                    // Anchoring `ref_aa`/`alt_aa` to the raw codon translation
                    // here would contradict that contract (e.g. `ref_aa = Leu`
                    // for a non-ATG start while `hgvs_p` says `Met1?`), so
                    // derive both from the predicted variant instead. See
                    // `substitution_residues_from_variant`.
                    let (ref_aa_str, alt_aa, position) =
                        substitution_residues_from_variant(&variant, ref_aa, prot_position);
                    return Some(ProteinConsequence {
                        hgvs_p: format!("{}", variant),
                        ref_aa: ref_aa_str,
                        alt_aa,
                        position,
                        is_frameshift: false,
                    });
                }
                Err(_) => {
                    // Prediction failed despite available sequence (e.g. an
                    // explicit reference base that mismatches the codon, or an
                    // unsupported codon shape). Represent this as an EXPLICIT
                    // unresolvable protein state (`p.(?N?)`) for the variant's
                    // own position — NOT a deferred indel notation. We still
                    // report the real reference residue we read from sequence.
                    return Some(ProteinConsequence {
                        hgvs_p: format!("{}:p.(?{}?)", prot_acc, prot_position),
                        ref_aa: ref_aa.to_three_letter().to_string(),
                        alt_aa: "?".to_string(),
                        position: prot_position,
                        is_frameshift: false,
                    });
                }
            }
        }
    }

    // Non-substitution edit classes: real reference residue, but the alternate
    // residue / full protein notation is #498's remit. The `?` here means
    // "not yet implemented for this edit class", not "data unavailable".
    let hgvs_p = build_indel_hgvs_p(prot_acc, edit_type, prot_position, is_frameshift, ref_len);
    Some(ProteinConsequence {
        hgvs_p,
        ref_aa: ref_aa.to_three_letter().to_string(),
        alt_aa: "?".to_string(),
        position: prot_position,
        is_frameshift,
    })
}

/// Derive the `(ref_aa, alt_aa, position)` residue fields for a substitution
/// `ProteinConsequence` from the engine's predicted protein `variant`, so they
/// never contradict the variant's `hgvs_p`.
///
/// `raw_ref_aa` (the literal codon translation) and `raw_position` are the
/// caller's fallbacks for the synonymous-identity shape, where `ref_aa` and
/// `alt_aa` are legitimately equal.
///
/// The cases that need the variant's own contract rather than the raw codon:
/// * **Plain residue substitution** (`p.(RefNAlt)`): take ref/alt/position from
///   the variant — the engine already resolved them.
/// * **Uncertain initiator / whole-protein unknown** (`p.(Met1?)`, `p.?`): the
///   protein consequence is unpredictable, so the alternate residue is `"?"`;
///   the reference residue and position come from the variant's location
///   (residue 1 is the initiator `Met` on an ATG start regardless of the raw
///   codon — anchoring `ref_aa` to the raw codon here would print e.g. `Leu`
///   for a non-ATG start while `hgvs_p` says `Met1?`).
/// * **C-terminal extension** (stop-loss, `p.(Ter3GlnextTer2)`): `ref_aa` is the
///   terminator at the location and `alt_aa` is the read-through residue — never
///   `Ter==Ter`, which the raw `ref/ref` fallback would wrongly produce.
/// * **Synonymous identity** (`p.(=)`): `ref_aa == alt_aa` is correct; use the
///   raw codon residue for both.
fn substitution_residues_from_variant(
    variant: &crate::hgvs::variant::HgvsVariant,
    raw_ref_aa: crate::hgvs::location::AminoAcid,
    raw_position: u64,
) -> (String, String, u64) {
    use crate::hgvs::edit::ProteinEdit;
    use crate::hgvs::variant::HgvsVariant;

    let raw_ref = raw_ref_aa.to_three_letter().to_string();
    let HgvsVariant::Protein(p) = variant else {
        return (raw_ref.clone(), raw_ref, raw_position);
    };
    // The variant's location residue/position — the engine's contract for the
    // resolved fields. Falls back to the raw codon position if the boundary is
    // unknown (should not happen for the shapes handled below).
    let (loc_ref_aa, loc_position) = match p.loc_edit.location.start.inner() {
        Some(pos) => (pos.aa.to_three_letter().to_string(), pos.number),
        None => (raw_ref.clone(), raw_position),
    };
    match p.loc_edit.edit.inner() {
        Some(ProteinEdit::Substitution { alternative, .. }) => (
            loc_ref_aa,
            alternative.to_three_letter().to_string(),
            loc_position,
        ),
        // Uncertain initiator (`p.(Met1?)`) or whole-protein unknown (`p.?`):
        // the alternate residue is unpredictable. Use the variant's location
        // residue/position so the fields agree with `hgvs_p`.
        Some(ProteinEdit::Unknown { .. }) => (loc_ref_aa, "?".to_string(), loc_position),
        // Identity (synonymous): no residue change, so `ref_aa == alt_aa` is the
        // honest report. The raw codon translation is the residue.
        Some(ProteinEdit::Identity { .. }) => (raw_ref.clone(), raw_ref, raw_position),
        // C-terminal extension (stop-loss, e.g. `p.(Ter3GlnextTer2)`): the
        // reference residue is the terminator at the variant's location and the
        // alternate is the read-through residue (`new_aa`). Reporting the raw
        // codon for both would print `ref_aa == alt_aa == Ter`, which is wrong —
        // a stop-loss is not synonymous. Anchor `ref_aa` to the location (the
        // `Ter`) and take the new residue from the edit when the engine resolved
        // it, falling back to `"?"` when it did not.
        Some(ProteinEdit::Extension { new_aa, .. }) => {
            let alt = new_aa
                .map(|aa| aa.to_three_letter().to_string())
                .unwrap_or_else(|| "?".to_string());
            (loc_ref_aa, alt, loc_position)
        }
        // Any other structured edit a substitution might predict (frameshift,
        // delins, …) is not representable as a single alternate residue here, so
        // don't claim the reference residue is the alternate. Anchor to the
        // location residue and mark the alternate unknown.
        _ => (loc_ref_aa, "?".to_string(), loc_position),
    }
}

/// Build position-based HGVS-p notation for a non-substitution indel when only
/// the reference residue is resolved (full notation deferred to #498).
fn build_indel_hgvs_p(
    prot_acc: &str,
    edit_type: &str,
    prot_position: u64,
    is_frameshift: bool,
    ref_len: usize,
) -> String {
    if is_frameshift {
        format!("{}:p.(?{}fs)", prot_acc, prot_position)
    } else {
        match edit_type {
            "deletion" => {
                // Preserve the ranged deletion form for in-frame multi-codon
                // deletions (e.g. `p.(?2_?3del)`); collapse to single-position
                // form only when the deletion spans one codon. Mirrors the
                // no-reference fallback in `predict_protein_consequence`.
                let end_pos = prot_position + (ref_len / 3).max(1) as u64 - 1;
                if prot_position == end_pos {
                    format!("{}:p.(?{}del)", prot_acc, prot_position)
                } else {
                    format!("{}:p.(?{}_?{}del)", prot_acc, prot_position, end_pos)
                }
            }
            "insertion" => format!(
                "{}:p.(?{}_?{}ins?)",
                prot_acc,
                prot_position,
                prot_position + 1
            ),
            _ => format!("{}:p.(?{}?)", prot_acc, prot_position),
        }
    }
}

/// The "50–55 nt" NMD-escape boundary: a premature termination codon (PTC)
/// within this many nucleotides upstream of the last exon-exon junction escapes
/// nonsense-mediated decay (Nagy & Maquat 1998; the canonical rule is commonly
/// cited as 50–55 nt). We use 55 nt — the conservative upper bound widely
/// adopted by annotation tools (e.g. Ensembl VEP's NMD plugin) — so that
/// borderline PTCs are classified as escaping rather than triggering.
const NMD_LAST_JUNCTION_ESCAPE_NT: u64 = 55;

/// Confidence for the NMD call when the variant introduces no PTC, or for a
/// single-exon transcript — both are structural certainties, not predictions.
const NMD_CONFIDENCE_STRUCTURAL: f64 = 0.9;

/// Confidence for a nonsense-substitution NMD call: the PTC location is exact
/// (the substituted codon), so the only residual uncertainty is the EJC model.
const NMD_CONFIDENCE_EXACT_PTC: f64 = 0.85;

/// Confidence for a frameshift NMD call: the reported PTC position is a 5' proxy
/// for the true (downstream) PTC, so these calls are less certain.
const NMD_CONFIDENCE_FRAMESHIFT_PTC: f64 = 0.7;

/// Predict NMD for a CDS variant using junction-based PTC logic (issue #806).
///
/// Replaces the previous "last 10% of CDS" fraction heuristic with the
/// biological rule: a PTC escapes NMD if it lies in the **last exon** or within
/// [`NMD_LAST_JUNCTION_ESCAPE_NT`] nucleotides upstream of the **last exon-exon
/// junction**; otherwise it triggers NMD. A single-exon transcript has no
/// exon-exon junction and therefore can never trigger NMD.
///
/// The exon structure and CDS bounds come from the coordinate-only
/// [`CdotTranscript`](crate::data::cdot::CdotTranscript) (always available when
/// cdot is loaded — NMD needs only coordinates, not sequence). The PTC location
/// is approximated by the variant's CDS position mapped to a transcript
/// position: for a nonsense substitution this is exact; for a frameshift the
/// true PTC is at or 3' of the variant, so using the variant position is the
/// standard conservative proxy (reflected in the confidence score).
///
/// Returns `None` (no NMD prediction) when the variant's PTC status cannot be
/// determined — either a substitution whose nonsense state is unknown (no
/// sequence-bearing reference is configured) or a delins whose net frame delta
/// is undeterminable (`frame_undecidable`) — rather than asserting a definitive
/// `predicted: false` we cannot justify.
fn predict_nmd_for_cds(
    state: &AppState,
    accession: &str,
    cds_pos: &CdsPos,
    edit: Option<&crate::hgvs::edit::NaEdit>,
    is_frameshift: bool,
    edit_type: &str,
    frame_undecidable: bool,
) -> Option<NmdPrediction> {
    let cdot = state.cdot.as_ref()?;
    let cdot_tx = cdot.get_transcript(accession)?;

    // Only an in-CDS, non-intronic position can introduce a PTC we can place.
    if cds_pos.base <= 0 || cds_pos.utr3 || cds_pos.offset.is_some() {
        return None;
    }

    // Does this variant plausibly introduce a PTC? Frameshifts do; a
    // substitution does iff it creates a stop codon (nonsense) — detected
    // exactly via the sequence-bearing reference. In-frame indels and
    // missense/synonymous substitutions do not introduce a PTC.
    //
    // Two "unknown" states must NOT collapse into "no PTC" — each would emit a
    // definitive `predicted: false` NMD call we cannot justify. Instead they
    // propagate `None` so this function returns `None` (no NMD prediction), the
    // honest signal that the call cannot be made:
    //   - a substitution whose nonsense state is undeterminable because no
    //     sequence-bearing reference is configured; and
    //   - a delins whose net frame delta is undeterminable (`frame_undecidable`):
    //     `is_frameshift = false` is conservative here, not "proven in-frame".
    let introduces_ptc = if is_frameshift {
        Some(true)
    } else if edit_type == "substitution" {
        is_nonsense_substitution(state, accession, cds_pos, edit)
    } else if frame_undecidable {
        None
    } else {
        Some(false)
    }?;

    nmd_from_junction(cdot_tx, cds_pos, is_frameshift, introduces_ptc)
}

/// Pure junction-based NMD decision, separated from `AppState` so it is unit
/// testable against a synthetic [`CdotTranscript`] (issue #806).
///
/// `introduces_ptc` is the caller's verdict on whether the variant introduces a
/// premature termination codon (frameshift, or a nonsense substitution). All
/// coordinates are taken from `cdot_tx` in its 0-based transcript coordinate
/// system; `cds_pos` is the variant's 1-based CDS position.
fn nmd_from_junction(
    cdot_tx: &crate::data::cdot::CdotTranscript,
    cds_pos: &CdsPos,
    is_frameshift: bool,
    introduces_ptc: bool,
) -> Option<NmdPrediction> {
    if !introduces_ptc {
        return Some(NmdPrediction {
            predicted: false,
            confidence: NMD_CONFIDENCE_STRUCTURAL,
            reason: "Variant does not introduce a premature termination codon".to_string(),
        });
    }

    // Map the PTC's 1-based CDS position to a 0-based transcript position.
    //
    // For a nonsense substitution the PTC is the stop codon the substitution
    // creates, so the distance-to-junction must be measured from that codon's
    // anchor (its first base), not the edited nucleotide. Normalizing to the
    // codon start makes every substitution that yields the same stop codon
    // classify identically against the 55-nt boundary; without it a hit on
    // codon base 2 or 3 shifts the distance by 1-2 nt and can flip a borderline
    // call. A frameshift's PTC position is a deliberate 5' proxy at the edited
    // base (see the confidence-constant docs), so it is left un-normalized.
    let ptc_cds_base = if is_frameshift {
        cds_pos.base
    } else {
        // Codon start (1-based): for base b, b - ((b - 1) mod 3).
        cds_pos.base - (cds_pos.base - 1).rem_euclid(3)
    };
    let ptc_tx = cdot_tx.cds_to_tx(ptc_cds_base)?;

    // No exon structure means we cannot reason about junctions at all — stay
    // undecided rather than emit a definitive "no NMD" call we cannot justify.
    if cdot_tx.exons.is_empty() {
        return None;
    }

    // Single-exon transcript: no exon-exon junction exists, so the EJC-based
    // NMD machinery has nothing to act on — a PTC here cannot trigger NMD.
    if cdot_tx.exons.len() == 1 {
        return Some(NmdPrediction {
            predicted: false,
            confidence: NMD_CONFIDENCE_STRUCTURAL,
            reason: "Single-exon transcript has no exon-exon junction; PTC cannot trigger NMD"
                .to_string(),
        });
    }

    // Last exon-exon junction = first base (tx-start, column index 2) of the
    // 3'-most exon in transcript order, in 0-based tx coordinates. Derive it as
    // the maximum tx-start across exons rather than trusting `.last()`, so a
    // missing or unsorted exon ordering invariant cannot silently flip a call.
    let last_junction_tx = cdot_tx.exons.iter().map(|exon| exon[2]).max()?;

    // A frameshift's PTC position is a 5' proxy (see the constant docs), so it
    // carries lower confidence than an exact nonsense-substitution PTC.
    let confidence = if is_frameshift {
        NMD_CONFIDENCE_FRAMESHIFT_PTC
    } else {
        NMD_CONFIDENCE_EXACT_PTC
    };

    let (predicted, reason) = if ptc_tx >= last_junction_tx {
        (false, "PTC in the last exon - escapes NMD".to_string())
    } else if last_junction_tx - ptc_tx <= NMD_LAST_JUNCTION_ESCAPE_NT {
        (
            false,
            format!(
                "PTC within {} nt upstream of the last exon-exon junction - escapes NMD",
                NMD_LAST_JUNCTION_ESCAPE_NT
            ),
        )
    } else {
        (
            true,
            format!(
                "PTC more than {} nt upstream of the last exon-exon junction - triggers NMD",
                NMD_LAST_JUNCTION_ESCAPE_NT
            ),
        )
    };

    Some(NmdPrediction {
        predicted,
        confidence,
        reason,
    })
}

/// Determine whether a substitution at `cds_pos` creates a stop codon (nonsense)
/// using the sequence-bearing reference.
///
/// Translates the reference and alternate codons (the latter via the engine's
/// [`predict_substitution_protein`](crate::project::protein::predict_substitution_protein))
/// and reports `true` iff the alternate residue is the terminator and the
/// reference residue was not — i.e. a *new* premature stop.
///
/// Both `c.Ref>Alt` and the no-reference `c.N>Alt` form are accepted; the latter
/// has its reference base resolved from sequence in the core (mirroring the
/// protein-resolution path) so it is not silently excluded from NMD prediction.
///
/// Returns `None` when sequence is unavailable (reference unconfigured or
/// transcript missing sequence) or the edit is not a substitution, so the caller
/// can treat PTC-introduction as unknown rather than asserting a nonsense call
/// without bases.
fn is_nonsense_substitution(
    state: &AppState,
    accession: &str,
    cds_pos: &CdsPos,
    edit: Option<&crate::hgvs::edit::NaEdit>,
) -> Option<bool> {
    use crate::hgvs::edit::NaEdit;

    // Both `c.Ref>Alt` (Substitution) and the no-reference `c.N>Alt`
    // (SubstitutionNoRef) are nonsense candidates — the core resolves the real
    // reference base from sequence for the latter. Reject only non-substitution
    // edits here so the no-ref form is not dropped before reaching the core.
    match edit? {
        NaEdit::Substitution { .. } | NaEdit::SubstitutionNoRef { .. } => {}
        _ => return None,
    }
    let provider = state.reference.as_ref()?;
    let transcript =
        crate::reference::provider::ReferenceProvider::get_transcript(provider.as_ref(), accession)
            .ok()?;
    is_nonsense_substitution_from_transcript(&transcript, accession, cds_pos, edit)
}

/// Pure nonsense-substitution detection core, separated from `AppState` for unit
/// testing against a synthetic sequence-bearing `Transcript` (issue #806). See
/// [`is_nonsense_substitution`] for the contract.
fn is_nonsense_substitution_from_transcript(
    transcript: &crate::reference::transcript::Transcript,
    accession: &str,
    cds_pos: &CdsPos,
    edit: Option<&crate::hgvs::edit::NaEdit>,
) -> Option<bool> {
    use crate::hgvs::edit::{Base, NaEdit, ProteinEdit};
    use crate::hgvs::variant::HgvsVariant;
    use crate::project::protein::{predict_substitution_protein, read_ref_codon};

    // Normalize the supported substitution forms to a `Substitution` with a
    // concrete reference base before calling the engine, which requires one.
    // `SubstitutionNoRef` (`c.N>Alt`) carries only the alternate allele, so
    // derive the reference base from the sequence codon at this position —
    // mirroring the protein-resolution path so sequence-backed no-ref nonsense
    // variants flow through the same NMD logic rather than being dropped.
    let normalized_sub;
    let edit = match edit? {
        sub @ NaEdit::Substitution { .. } => sub,
        NaEdit::SubstitutionNoRef { alternative } => {
            let (ref_codon, frame) = read_ref_codon(transcript, cds_pos.base).ok()?;
            let reference = Base::from_char(ref_codon.as_bytes()[frame as usize] as char)?;
            normalized_sub = NaEdit::Substitution {
                reference,
                alternative: *alternative,
            };
            &normalized_sub
        }
        _ => return None,
    };
    let prot_acc = transcript
        .protein_id
        .clone()
        .unwrap_or_else(|| accession.to_string());
    let variant = predict_substitution_protein(transcript, cds_pos.base, edit, &prot_acc).ok()?;

    // A substitution introduces a new PTC iff its predicted protein edit is a
    // residue substitution whose alternate residue is the terminator. A
    // reference stop (stop-loss) yields an extension, not a Substitution, so it
    // correctly returns `false` here.
    let HgvsVariant::Protein(p) = variant else {
        return Some(false);
    };
    match p.loc_edit.edit.inner() {
        Some(ProteinEdit::Substitution {
            reference,
            alternative,
        }) => Some(
            *alternative == crate::hgvs::location::AminoAcid::Ter
                && *reference != crate::hgvs::location::AminoAcid::Ter,
        ),
        _ => Some(false),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // -------------------------------------------------------------------------
    // Issue #806: real amino-acid resolution + junction-based NMD
    // -------------------------------------------------------------------------

    use crate::data::cdot::CdotTranscript;
    use crate::hgvs::edit::{Base, NaEdit};
    use crate::reference::transcript::{Exon as TxExon, Strand, Transcript};

    /// Build a sequence-bearing single-exon `Transcript` whose CDS starts at
    /// tx position 1 (1-based inclusive). `seq` is the full transcript = CDS.
    fn seq_transcript(id: &str, seq: &str, protein_id: Option<&str>) -> Transcript {
        Transcript {
            id: id.to_string(),
            strand: Strand::Plus,
            sequence: Some(seq.to_string()),
            cds_start: Some(1),
            cds_end: Some(seq.len() as u64),
            exons: vec![TxExon::new(1, 1, seq.len() as u64)],
            protein_id: protein_id.map(|s| s.to_string()),
            ..Default::default()
        }
    }

    /// Like [`seq_transcript`] but with an explicit 1-based inclusive CDS end
    /// shorter than the sequence — so the bases past `cds_end` act as a 3'UTR
    /// (needed to exercise the stop-loss read-through / extension path).
    fn seq_transcript_with_cds(
        id: &str,
        seq: &str,
        cds_end: u64,
        protein_id: Option<&str>,
    ) -> Transcript {
        Transcript {
            id: id.to_string(),
            strand: Strand::Plus,
            sequence: Some(seq.to_string()),
            cds_start: Some(1),
            cds_end: Some(cds_end),
            exons: vec![TxExon::new(1, 1, seq.len() as u64)],
            protein_id: protein_id.map(|s| s.to_string()),
            ..Default::default()
        }
    }

    /// Build a `CdotTranscript` with the given exon tx-spans (each `(tx_start,
    /// tx_end_excl)`, 0-based) and `cds_start = 0` (CDS begins at tx 0).
    fn cdot_with_exons(exon_tx_spans: &[(u64, u64)]) -> CdotTranscript {
        let exons: Vec<[u64; 4]> = exon_tx_spans
            .iter()
            .enumerate()
            .map(|(i, &(ts, te))| {
                // Arbitrary but monotonic genomic coords; NMD math uses only the
                // tx columns (indices 2 and 3).
                let g = 1000 + (i as u64) * 1000;
                [g, g + (te - ts), ts, te]
            })
            .collect();
        let tx_len = exon_tx_spans.last().map(|&(_, te)| te).unwrap_or(0);
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: None,
            contig: "NC_000001.11".to_string(),
            strand: Strand::Plus,
            exons,
            cds_start: Some(0),
            cds_end: Some(tx_len),
            exon_cigars: Vec::new(),
            gene_id: None,
            protein: None,
        }
    }

    fn cds(base: i64) -> CdsPos {
        CdsPos::new(base)
    }

    // ---- Amino-acid resolution (real residues, not "?") ----

    #[test]
    fn resolve_real_missense_residues() {
        // CDS: ATG (Met) CGT (Arg) ... ; substitute CGT -> CAT = His at codon 2.
        // c.5G>A changes the middle base of codon 2 (CGT) to CAT.
        let seq = "ATGCGTGGGTAA"; // Met Arg Gly Ter
        let tx = seq_transcript("NM_TEST.1", seq, Some("NP_TEST.1"));
        let edit = NaEdit::Substitution {
            reference: Base::G,
            alternative: Base::A,
        };
        let pc = resolve_residues_from_transcript(
            &tx,
            "NP_TEST.1",
            &cds(5),
            Some(&edit),
            "substitution",
            false,
            1,
        )
        .expect("residues should resolve");
        assert_eq!(pc.ref_aa, "Arg", "real reference residue, not '?'");
        assert_eq!(pc.alt_aa, "His", "real alternate residue, not '?'");
        assert_eq!(pc.position, 2);
        assert!(!pc.ref_aa.contains('?') && !pc.alt_aa.contains('?'));
    }

    #[test]
    fn resolve_indel_with_reference_gives_real_ref_residue_not_silent_question() {
        // Finding #4 guard: when a reference IS available, an indel still
        // resolves the real *reference* residue (alt_aa is the documented
        // not-yet-implemented "?", never a silent "?" for ref_aa).
        let seq = "ATGCGTGGGTAA"; // Met Arg Gly Ter
        let tx = seq_transcript("NM_TEST.1", seq, Some("NP_TEST.1"));
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let pc = resolve_residues_from_transcript(
            &tx,
            "NP_TEST.1",
            &cds(4), // codon 2 = CGT = Arg
            Some(&edit),
            "deletion",
            true,
            3,
        )
        .expect("ref residue should resolve");
        assert_eq!(
            pc.ref_aa, "Arg",
            "indel-with-reference resolves real ref residue"
        );
        assert_eq!(
            pc.alt_aa, "?",
            "indel alt residue is deferred to #498 (not-yet-implemented)"
        );
    }

    #[test]
    fn resolve_inframe_multicodon_deletion_keeps_ranged_form() {
        // A sequence-backed in-frame deletion spanning >1 codon must keep the
        // ranged HGVS-p form `p.(?N_?Mdel)` rather than collapsing to the
        // single-position `p.(?Ndel)`. `ref_len` flows through
        // `resolve_residues_from_transcript` into `build_indel_hgvs_p`.
        let seq = "ATGCGTGGGGAATAA"; // Met Arg Gly Glu Ter
        let tx = seq_transcript("NM_TEST.1", seq, Some("NP_TEST.1"));
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let pc = resolve_residues_from_transcript(
            &tx,
            "NP_TEST.1",
            &cds(4), // codon 2
            Some(&edit),
            "deletion",
            false,
            6, // two codons deleted -> aa range of length 2
        )
        .expect("ref residue should resolve");
        assert_eq!(pc.ref_aa, "Arg", "real reference residue at codon 2");
        assert_eq!(
            pc.hgvs_p, "NP_TEST.1:p.(?2_?3del)",
            "multi-codon deletion keeps the ranged form, not collapsed p.(?2del)"
        );
    }

    #[test]
    fn resolve_returns_none_without_sequence() {
        // No sequence on the transcript => cannot read a codon => None, so the
        // caller emits the honest "data unavailable" "?" fallback.
        let mut tx = seq_transcript("NM_TEST.1", "ATGCGTGGGTAA", Some("NP_TEST.1"));
        tx.sequence = None;
        let edit = NaEdit::Substitution {
            reference: Base::G,
            alternative: Base::A,
        };
        let pc = resolve_residues_from_transcript(
            &tx,
            "NP_TEST.1",
            &cds(5),
            Some(&edit),
            "substitution",
            false,
            1,
        );
        assert!(
            pc.is_none(),
            "no sequence => None => honest 'data unavailable'"
        );
    }

    #[test]
    fn resolve_initiation_codon_substitution_residues_match_hgvs_p() {
        // A substitution anywhere in the initiation codon (CDS 1-3) has an
        // unpredictable protein consequence: the engine emits `p.(Met1?)` on an
        // ATG start. The residue fields must agree with that contract -- ref_aa
        // = Met (residue 1 is the initiator Met) and alt_aa = "?" -- NOT the raw
        // codon translation, which would otherwise contradict `hgvs_p` for a
        // non-ATG-changing edit and (here) yield a bogus concrete alt residue.
        let seq = "ATGCGTGGGTAA"; // Met Arg Gly Ter
        let tx = seq_transcript("NM_TEST.1", seq, Some("NP_TEST.1"));
        // c.2T>A hits codon 1 (ATG); engine returns p.(Met1?).
        let edit = NaEdit::Substitution {
            reference: Base::T,
            alternative: Base::A,
        };
        let pc = resolve_residues_from_transcript(
            &tx,
            "NP_TEST.1",
            &cds(2),
            Some(&edit),
            "substitution",
            false,
            1,
        )
        .expect("initiation-codon substitution should resolve");
        assert_eq!(pc.hgvs_p, "NP_TEST.1:p.(Met1?)");
        assert_eq!(pc.ref_aa, "Met", "residue 1 is the initiator Met");
        assert_eq!(
            pc.alt_aa, "?",
            "initiation-codon alt residue is unpredictable, must not contradict Met1?"
        );
        assert_eq!(pc.position, 1);
    }

    #[test]
    fn resolve_non_atg_initiation_codon_substitution_residues_match_hgvs_p() {
        // On a non-ATG start codon a raw translation of codon 1 would give a
        // concrete non-Met residue (e.g. Leu for CTG), which contradicts the
        // engine's initiator output. The substitution path reports `p.(Met1?)`
        // for any initiation-codon substitution; the residue fields must agree
        // (ref_aa = Met, alt_aa = "?"), never the raw codon residue.
        let seq = "CTGCGTGGGTAA"; // CTG (Leu by raw translation) start ...
        let tx = seq_transcript("NM_TEST.1", seq, Some("NP_TEST.1"));
        // c.1C>A hits the non-ATG start codon.
        let edit = NaEdit::Substitution {
            reference: Base::C,
            alternative: Base::A,
        };
        let pc = resolve_residues_from_transcript(
            &tx,
            "NP_TEST.1",
            &cds(1),
            Some(&edit),
            "substitution",
            false,
            1,
        )
        .expect("non-ATG initiation-codon substitution should resolve");
        assert_eq!(pc.hgvs_p, "NP_TEST.1:p.(Met1?)");
        assert_ne!(
            pc.ref_aa, "Leu",
            "must not report the raw CTG residue that contradicts the initiator output"
        );
        assert_eq!(
            pc.ref_aa, "Met",
            "initiator residue 1 is Met per the engine"
        );
        assert_eq!(pc.alt_aa, "?");
    }

    #[test]
    fn resolve_uses_cdot_protein_accession_when_transcript_lacks_protein_id() {
        // Finding 1 (#806 review): the sequence-backed path must use the
        // authoritative cdot protein accession (threaded as `prot_acc`), NOT
        // `transcript.protein_id`. A reference `Transcript` can lack
        // `protein_id` even when cdot carries `NP_*`; resolving from the
        // transcript there would flip `NP_*:p...` to `NM_*:p...` purely by
        // enabling `state.reference`. Here the transcript has NO protein_id but
        // we pass the cdot accession `NP_TEST.1` — the HGVS-p must carry it.
        let seq = "ATGCGTGGGTAA"; // Met Arg Gly Ter
        let tx = seq_transcript("NM_TEST.1", seq, None); // no protein_id
        let edit = NaEdit::Substitution {
            reference: Base::G,
            alternative: Base::A,
        };
        let pc = resolve_residues_from_transcript(
            &tx,
            "NP_TEST.1", // cdot-derived protein accession
            &cds(5),
            Some(&edit),
            "substitution",
            false,
            1,
        )
        .expect("residues should resolve");
        assert_eq!(
            pc.hgvs_p, "NP_TEST.1:p.(Arg2His)",
            "uses cdot NP_ accession even though transcript.protein_id is None"
        );
        assert!(
            !pc.hgvs_p.starts_with("NM_"),
            "must not fall back to the NM_ transcript accession"
        );

        // Same guarantee for the deferred indel path (build_indel_hgvs_p).
        let del = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let pc_del = resolve_residues_from_transcript(
            &tx,
            "NP_TEST.1",
            &cds(4),
            Some(&del),
            "deletion",
            true,
            1,
        )
        .expect("ref residue should resolve");
        assert!(
            pc_del.hgvs_p.starts_with("NP_TEST.1:"),
            "indel path also uses the cdot NP_ accession, got {:?}",
            pc_del.hgvs_p
        );
    }

    #[test]
    fn resolve_substitution_no_ref_resolves_real_alt_residue() {
        // `c.N>Alt` (SubstitutionNoRef) carries only the alternate allele. When
        // sequence exists the reference base is derived from the codon, so the
        // form must resolve the REAL alternate residue — never degrade to the
        // deferred `p.(?N?)` / indel notation (#806 review, finding 2a).
        // CDS: ATG (Met) CGT (Arg) ...; c.5N>A rewrites the middle base of codon
        // 2 (CGT -> CAT = His), matching `resolve_real_missense_residues`.
        let seq = "ATGCGTGGGTAA"; // Met Arg Gly Ter
        let tx = seq_transcript("NM_TEST.1", seq, Some("NP_TEST.1"));
        let edit = NaEdit::SubstitutionNoRef {
            alternative: Base::A,
        };
        let pc = resolve_residues_from_transcript(
            &tx,
            "NP_TEST.1",
            &cds(5),
            Some(&edit),
            "substitution",
            false,
            1,
        )
        .expect("no-ref substitution should resolve against sequence");
        assert_eq!(pc.ref_aa, "Arg", "ref residue read from the codon");
        assert_eq!(pc.alt_aa, "His", "no-ref alt residue resolved, not '?'");
        assert_eq!(pc.position, 2);
        assert_eq!(
            pc.hgvs_p, "NP_TEST.1:p.(Arg2His)",
            "no-ref form yields fully-resolved HGVS-p, not p.(?2?)"
        );
        assert!(
            !pc.hgvs_p.contains("?"),
            "no deferred '?' when sequence exists"
        );
    }

    #[test]
    fn resolve_substitution_explicit_ref_mismatch_still_resolves() {
        // An explicit reference base that disagrees with the sequence must NOT
        // collapse to a sequence-backed `p.(?N?)`: the engine reads the real
        // codon from sequence and applies the alt regardless of the stated ref,
        // so the prediction still resolves to real residues (#806 review,
        // finding 2b). Codon 2 is CGT (Arg); we state ref=A (wrong) but alt=A,
        // giving CGT -> CAT = His just as the correct-ref case would.
        let seq = "ATGCGTGGGTAA"; // Met Arg Gly Ter
        let tx = seq_transcript("NM_TEST.1", seq, Some("NP_TEST.1"));
        let edit = NaEdit::Substitution {
            reference: Base::A, // deliberately wrong vs. the codon's G
            alternative: Base::A,
        };
        let pc = resolve_residues_from_transcript(
            &tx,
            "NP_TEST.1",
            &cds(5),
            Some(&edit),
            "substitution",
            false,
            1,
        )
        .expect("explicit-ref-mismatch substitution should still resolve");
        assert_eq!(
            pc.hgvs_p, "NP_TEST.1:p.(Arg2His)",
            "explicit-ref mismatch resolves from sequence, not p.(?2?)"
        );
        assert!(
            !pc.hgvs_p.contains("?"),
            "explicit-ref mismatch must not yield a sequence-backed p.(?N?)"
        );
        assert_eq!(pc.ref_aa, "Arg");
        assert_eq!(pc.alt_aa, "His");
    }

    #[test]
    fn resolve_stop_loss_extension_reports_readthrough_residue_not_ter_ter() {
        // A stop-loss substitution predicts a C-terminal extension
        // (`p.(Ter3GlnextTer2)`), NOT a residue substitution. The residue fields
        // must report the read-through residue as the alternate — never the raw
        // `Ter==Ter` the generic ref/ref fallback would produce (#806 review:
        // "stop-loss substitutions report alt_aa as the reference residue").
        // CDS "ATGAAATAA" (Met-Lys-Ter; stop at c.7-9) + 3'UTR "GGGTAA".
        // c.7T>C turns TAA -> CAA (Gln), reading through to a downstream Ter.
        let seq = "ATGAAATAAGGGTAA";
        let tx = seq_transcript_with_cds("NM_TEST.1", seq, 9, Some("NP_TEST.1"));
        let edit = NaEdit::Substitution {
            reference: Base::T,
            alternative: Base::C,
        };
        let pc = resolve_residues_from_transcript(
            &tx,
            "NP_TEST.1",
            &cds(7),
            Some(&edit),
            "substitution",
            false,
            1,
        )
        .expect("stop-loss substitution should resolve");
        assert_eq!(
            pc.hgvs_p, "NP_TEST.1:p.(Ter3GlnextTer2)",
            "stop-loss yields a C-terminal extension"
        );
        assert_eq!(pc.ref_aa, "Ter", "reference residue is the terminator");
        assert_eq!(
            pc.alt_aa, "Gln",
            "alternate is the read-through residue, NOT Ter (no ref==alt)"
        );
        assert_ne!(
            pc.ref_aa, pc.alt_aa,
            "a stop-loss is not synonymous: ref_aa must not equal alt_aa"
        );
        assert_eq!(pc.position, 3, "extension is at the original stop codon");
    }

    #[test]
    fn resolve_synonymous_substitution_reports_equal_ref_and_alt() {
        // A synonymous substitution predicts identity (`p.(=)`), where reporting
        // `ref_aa == alt_aa` is the honest, correct result. CDS codon 2 is CGT
        // (Arg); c.6T>A -> CGA, still Arg.
        let seq = "ATGCGTGGGTAA"; // Met Arg Gly Ter
        let tx = seq_transcript("NM_TEST.1", seq, Some("NP_TEST.1"));
        let edit = NaEdit::Substitution {
            reference: Base::T,
            alternative: Base::A,
        };
        let pc = resolve_residues_from_transcript(
            &tx,
            "NP_TEST.1",
            &cds(6),
            Some(&edit),
            "substitution",
            false,
            1,
        )
        .expect("synonymous substitution should resolve");
        assert_eq!(
            pc.ref_aa, "Arg",
            "synonymous change keeps the reference residue"
        );
        assert_eq!(pc.alt_aa, "Arg", "synonymous: alt residue equals ref");
        assert_eq!(pc.position, 2);
    }

    // ---- Nonsense detection ----

    #[test]
    fn detects_nonsense_substitution() {
        // Codon 2 CGT (Arg). C>T at c.4 makes TGT (Cys) — not nonsense.
        // Use CAA (Gln) at codon 2 and C>T at c.4 -> TAA (Ter): nonsense.
        let seq = "ATGCAAGGGTAA"; // Met Gln Gly Ter
        let tx = seq_transcript("NM_TEST.1", seq, Some("NP_TEST.1"));
        let edit = NaEdit::Substitution {
            reference: Base::C,
            alternative: Base::T,
        };
        let is_nonsense =
            is_nonsense_substitution_from_transcript(&tx, "NM_TEST.1", &cds(4), Some(&edit));
        assert_eq!(is_nonsense, Some(true), "CAA->TAA is a new stop (nonsense)");
    }

    #[test]
    fn missense_is_not_nonsense() {
        let seq = "ATGCGTGGGTAA"; // Met Arg Gly Ter
        let tx = seq_transcript("NM_TEST.1", seq, Some("NP_TEST.1"));
        let edit = NaEdit::Substitution {
            reference: Base::G,
            alternative: Base::A,
        };
        let is_nonsense =
            is_nonsense_substitution_from_transcript(&tx, "NM_TEST.1", &cds(5), Some(&edit));
        assert_eq!(
            is_nonsense,
            Some(false),
            "CGT->CAT (His) is missense, not nonsense"
        );
    }

    #[test]
    fn nonsense_unknown_without_sequence_stays_none() {
        // Without a sequence the nonsense state is undeterminable; the detector
        // must report `None` (unknown), NOT `Some(false)`. `predict_nmd_for_cds`
        // relies on this to decline an NMD call rather than emit a false-negative
        // `predicted: false`.
        let mut tx = seq_transcript("NM_TEST.1", "ATGCAAGGGTAA", Some("NP_TEST.1"));
        tx.sequence = None;
        let edit = NaEdit::Substitution {
            reference: Base::C,
            alternative: Base::T,
        };
        let is_nonsense =
            is_nonsense_substitution_from_transcript(&tx, "NM_TEST.1", &cds(4), Some(&edit));
        assert_eq!(
            is_nonsense, None,
            "no sequence => unknown nonsense state, never a coerced false"
        );
    }

    #[test]
    fn detects_nonsense_no_ref_substitution() {
        // `c.N>Alt` (SubstitutionNoRef) must flow through nonsense detection just
        // like `c.Ref>Alt`: the reference base is derived from sequence, so a
        // no-ref nonsense variant is not silently excluded from NMD (#806 review:
        // "Normalize no-ref substitutions before NMD detection"). Codon 2 CAA
        // (Gln); c.4N>T derives ref C from the codon -> TAA (Ter): nonsense.
        let seq = "ATGCAAGGGTAA"; // Met Gln Gly Ter
        let tx = seq_transcript("NM_TEST.1", seq, Some("NP_TEST.1"));
        let edit = NaEdit::SubstitutionNoRef {
            alternative: Base::T,
        };
        let is_nonsense =
            is_nonsense_substitution_from_transcript(&tx, "NM_TEST.1", &cds(4), Some(&edit));
        assert_eq!(
            is_nonsense,
            Some(true),
            "no-ref c.4N>T (CAA->TAA) is a new stop (nonsense), not dropped"
        );
    }

    #[test]
    fn no_ref_missense_substitution_is_not_nonsense() {
        // No-ref counterpart of `missense_is_not_nonsense`: a resolved missense
        // must report `Some(false)`, not `None`. Codon 2 CGT (Arg); c.5N>A
        // derives ref G -> CAT (His): missense.
        let seq = "ATGCGTGGGTAA"; // Met Arg Gly Ter
        let tx = seq_transcript("NM_TEST.1", seq, Some("NP_TEST.1"));
        let edit = NaEdit::SubstitutionNoRef {
            alternative: Base::A,
        };
        let is_nonsense =
            is_nonsense_substitution_from_transcript(&tx, "NM_TEST.1", &cds(5), Some(&edit));
        assert_eq!(
            is_nonsense,
            Some(false),
            "no-ref CGT->CAT (His) is missense, resolved (not unknown)"
        );
    }

    // ---- Junction-based NMD ----

    #[test]
    fn nmd_ptc_in_last_exon_escapes() {
        // Exons (tx 0-based): [0,100),[100,200),[200,250). Last junction tx=200.
        let cdot_tx = cdot_with_exons(&[(0, 100), (100, 200), (200, 250)]);
        // c.201 -> tx 200 == last junction -> last exon -> escape.
        let pred = nmd_from_junction(&cdot_tx, &cds(201), true, true).unwrap();
        assert!(!pred.predicted, "PTC in last exon escapes NMD");
        assert!(pred.reason.contains("last exon"));
    }

    #[test]
    fn nmd_ptc_within_55nt_of_last_junction_escapes() {
        let cdot_tx = cdot_with_exons(&[(0, 100), (100, 200), (200, 250)]);
        // c.180 -> tx 179 -> 200-179 = 21 nt upstream -> escape.
        let pred = nmd_from_junction(&cdot_tx, &cds(180), true, true).unwrap();
        assert!(!pred.predicted, "PTC within 55 nt of last junction escapes");
    }

    #[test]
    fn nmd_55nt_boundary_is_inclusive() {
        let cdot_tx = cdot_with_exons(&[(0, 100), (100, 200), (200, 250)]);
        // c.146 -> tx 145 -> 200-145 = 55 -> exactly at boundary -> escape.
        let escape = nmd_from_junction(&cdot_tx, &cds(146), true, true).unwrap();
        assert!(
            !escape.predicted,
            "exactly 55 nt upstream escapes (inclusive)"
        );
        // c.145 -> tx 144 -> 200-144 = 56 -> beyond boundary -> trigger.
        let trigger = nmd_from_junction(&cdot_tx, &cds(145), true, true).unwrap();
        assert!(trigger.predicted, "56 nt upstream triggers NMD");
    }

    #[test]
    fn nmd_ptc_mid_cds_triggers() {
        let cdot_tx = cdot_with_exons(&[(0, 100), (100, 200), (200, 250)]);
        // c.10 -> tx 9 -> 200-9 = 191 nt upstream -> trigger.
        let pred = nmd_from_junction(&cdot_tx, &cds(10), true, true).unwrap();
        assert!(pred.predicted, "PTC far 5' of last junction triggers NMD");
        assert!(pred.reason.contains("triggers NMD"));
    }

    #[test]
    fn nmd_single_exon_never_triggers() {
        let cdot_tx = cdot_with_exons(&[(0, 250)]);
        // Even a frameshift PTC near the 5' end cannot trigger NMD: no junction.
        let pred = nmd_from_junction(&cdot_tx, &cds(10), true, true).unwrap();
        assert!(!pred.predicted, "single-exon transcript never triggers NMD");
        assert!(pred.reason.contains("Single-exon"));
    }

    #[test]
    fn nmd_empty_exons_is_undecided() {
        // Missing exon structure must stay undecided (`None`), not collapse into
        // a definitive "no NMD" call: `cds_to_tx` only needs `cds_start`, so the
        // junction logic is reached even with no exons (#806 review).
        let cdot_tx = cdot_with_exons(&[]);
        assert!(
            nmd_from_junction(&cdot_tx, &cds(4), false, true).is_none(),
            "empty exon metadata yields no NMD prediction"
        );
    }

    #[test]
    fn nmd_last_junction_uses_max_tx_start_not_order() {
        // The last junction is the maximum exon tx-start, independent of the
        // order exons appear in: an unsorted list must classify identically to
        // the sorted one (#806 review: derive the junction by max, not `.last()`).
        let sorted = cdot_with_exons(&[(0, 100), (100, 200), (200, 250)]);
        let unsorted = cdot_with_exons(&[(200, 250), (0, 100), (100, 200)]);
        // c.50 -> tx 49. Against the true last junction (tx 200) this is 151 nt
        // upstream -> triggers NMD. A naive `.last()` on the unsorted list would
        // pick tx-start 100 (51 nt upstream) and wrongly escape, so this pins the
        // order-independent max-based derivation.
        let from_sorted = nmd_from_junction(&sorted, &cds(50), true, true).unwrap();
        let from_unsorted = nmd_from_junction(&unsorted, &cds(50), true, true).unwrap();
        assert!(from_sorted.predicted);
        assert_eq!(
            from_sorted.predicted, from_unsorted.predicted,
            "junction derived by max tx-start is order-independent"
        );
    }

    #[test]
    fn nmd_nonsense_substitution_mid_cds_triggers_end_to_end() {
        // Ties nonsense detection (sequence-bearing Transcript) to the
        // junction decision (CdotTranscript): a CAA->TAA nonsense at codon 2
        // (c.4), far 5' of the last junction, triggers NMD.
        let seq = "ATGCAAGGGTAA"; // Met Gln Gly Ter
        let tx = seq_transcript("NM_TEST.1", seq, Some("NP_TEST.1"));
        let edit = NaEdit::Substitution {
            reference: Base::C,
            alternative: Base::T,
        };
        let introduces_ptc =
            is_nonsense_substitution_from_transcript(&tx, "NM_TEST.1", &cds(4), Some(&edit))
                .unwrap_or(false);
        assert!(introduces_ptc, "nonsense substitution introduces a PTC");

        // Multi-exon transcript with the last junction well 3' of codon 2.
        let cdot_tx = cdot_with_exons(&[(0, 100), (100, 200), (200, 250)]);
        let pred = nmd_from_junction(&cdot_tx, &cds(4), false, introduces_ptc).unwrap();
        assert!(pred.predicted, "nonsense PTC mid-CDS triggers NMD");
    }

    #[test]
    fn nmd_nonsense_substitution_is_phase_independent_at_boundary() {
        // A nonsense substitution must classify by the stop CODON's anchor, not
        // the edited nucleotide: all three bases of the same codon yield the
        // same stop and must classify identically against the 55-nt boundary.
        //
        // Exons [0,100),[100,200),[200,250): last junction tx=200. The codon at
        // c.145/146/147 (codon start c.145 -> tx 144 -> distance 56 -> trigger)
        // straddles the boundary by edited-base: c.146 -> tx 145 -> distance 55
        // would naively escape. Normalizing to codon start makes all three
        // trigger.
        let cdot_tx = cdot_with_exons(&[(0, 100), (100, 200), (200, 250)]);
        for base in [145, 146, 147] {
            let pred = nmd_from_junction(&cdot_tx, &cds(base), false, true).unwrap();
            assert!(
                pred.predicted,
                "nonsense at c.{base} (codon start c.145, 56 nt upstream) must trigger NMD"
            );
        }
    }

    #[test]
    fn nmd_frameshift_uses_edited_base_not_codon_start() {
        // A frameshift's PTC position is a 5' proxy at the edited base -- it must
        // NOT be normalized to a codon start (that would shift the proxy and the
        // classification). c.146 -> tx 145 -> distance 55 -> escape for a
        // frameshift, distinct from the nonsense-substitution case above.
        let cdot_tx = cdot_with_exons(&[(0, 100), (100, 200), (200, 250)]);
        let pred = nmd_from_junction(&cdot_tx, &cds(146), true, true).unwrap();
        assert!(
            !pred.predicted,
            "frameshift proxy at tx 145 (55 nt upstream) escapes; not codon-normalized"
        );
    }

    #[test]
    fn nmd_no_ptc_does_not_trigger() {
        let cdot_tx = cdot_with_exons(&[(0, 100), (100, 200), (200, 250)]);
        // In-frame / missense (introduces_ptc = false), even mid-CDS.
        let pred = nmd_from_junction(&cdot_tx, &cds(10), false, false).unwrap();
        assert!(
            !pred.predicted,
            "variant without a PTC does not trigger NMD"
        );
        assert!(pred.reason.contains("does not introduce"));
    }

    #[test]
    fn test_analyze_substitution() {
        let result = crate::hgvs::parser::parse_hgvs_lenient("NM_000249.4:c.350C>T").unwrap();
        if let crate::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
            if let Some(edit) = v.loc_edit.edit.inner() {
                let (edit_type, is_frameshift, _, _) = analyze_na_edit(edit, None);
                assert_eq!(edit_type, "substitution");
                assert!(!is_frameshift);
            }
        }
    }

    #[test]
    fn test_analyze_frameshift_deletion() {
        // Updated for #427: pass `span_len` derived from the interval
        // (mirrors the production caller in `analyze_variant_effect`).
        // Single-position `c.350del` → span_len = 1 → 1 bp deletion =
        // frameshift.
        let result = crate::hgvs::parser::parse_hgvs_lenient("NM_000249.4:c.350del").unwrap();
        if let crate::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
            if let Some(edit) = v.loc_edit.edit.inner() {
                let span_len = span_len_from_cds_interval(&v.loc_edit.location);
                let (edit_type, is_frameshift, _, _) = analyze_na_edit(edit, span_len);
                assert_eq!(edit_type, "deletion");
                assert!(is_frameshift); // Single base deletion causes frameshift
            }
        }
    }

    #[test]
    fn test_analyze_inframe_deletion() {
        // Closes the "TODO: calculate actual length" gap from the
        // pre-#427 code. The 3-bp range deletion `c.350_352del` is now
        // correctly identified as in-frame via `span_len = 3`.
        let result = crate::hgvs::parser::parse_hgvs_lenient("NM_000249.4:c.350_352del").unwrap();
        if let crate::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
            if let Some(edit) = v.loc_edit.edit.inner() {
                let span_len = span_len_from_cds_interval(&v.loc_edit.location);
                let (edit_type, is_frameshift, ref_len, alt_len) = analyze_na_edit(edit, span_len);
                assert_eq!(edit_type, "deletion");
                assert_eq!(ref_len, 3);
                assert_eq!(alt_len, 0);
                assert!(!is_frameshift, "3 bp deletion is in-frame");
            }
        }
    }

    // ---- delins frame-undecidability → NMD declines (issue #806 review) ----
    //
    // A delins whose net frame delta is undeterminable reports `is_frameshift =
    // false` *conservatively*. `predict_nmd_for_cds` keys off this flag to
    // return `None` (no NMD prediction) instead of a definitive `predicted:
    // false` it cannot justify — these tests pin the flag for the shapes that
    // feed that decision.

    #[test]
    fn delins_unknown_alt_len_is_frame_undecidable() {
        use crate::hgvs::edit::InsertedSequence;
        // Ref span known (3 bp) but the inserted nt count is unknown
        // (`Range`/`Named`/`Reference`-style insert): net delta undeterminable.
        let edit = NaEdit::Delins {
            sequence: InsertedSequence::Range(5, 10),
            deleted: None,
            deleted_length: None,
        };
        assert!(
            delins_frame_undecidable(&edit, Some(3)),
            "unknown inserted length leaves the net frame delta undecidable"
        );
    }

    #[test]
    fn delins_unknown_span_with_literal_is_frame_undecidable() {
        use crate::hgvs::edit::{InsertedSequence, Sequence};
        // Inserted length known (single literal base) but the deleted span is
        // unknown (no position-interval span): still undecidable. This is the
        // "unknown-span delins with a single-nt literal insert" shape.
        let edit = NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::new(vec![Base::A])),
            deleted: None,
            deleted_length: None,
        };
        assert!(
            delins_frame_undecidable(&edit, None),
            "unknown deleted span leaves the net frame delta undecidable"
        );
    }

    #[test]
    fn delins_known_span_and_alt_is_frame_decidable() {
        use crate::hgvs::edit::{InsertedSequence, Sequence};
        // Both endpoints known (3 bp deleted, 2 bp inserted): the net delta is
        // fully determined, so this is NOT undecidable — NMD may proceed.
        let edit = NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::new(vec![Base::A, Base::T])),
            deleted: None,
            deleted_length: None,
        };
        assert!(
            !delins_frame_undecidable(&edit, Some(3)),
            "both ref span and alt length known => frame delta is decidable"
        );
    }

    #[test]
    fn non_delins_is_never_frame_undecidable() {
        // A del/dup/ins with an unknown span is handled by the caller's separate
        // `span_undecidable` skip, not this delins-specific flag — so the flag
        // must stay `false` for non-delins edits regardless of `span_len`.
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        assert!(!delins_frame_undecidable(&edit, None));
        assert!(!delins_frame_undecidable(&edit, Some(2)));
    }

    #[test]
    fn test_predict_splice_site_effect() {
        let cds_pos = CdsPos {
            base: 117,
            offset: Some(-2),
            utr3: false,
            special: None,
        };
        let effect = predict_cds_effect("deletion", true, false, &cds_pos);
        assert_eq!(effect.name, "splice_site_variant");
        assert_eq!(effect.impact, "HIGH");
    }

    #[test]
    fn test_predict_splice_region_effect() {
        let cds_pos = CdsPos {
            base: 117,
            offset: Some(-5),
            utr3: false,
            special: None,
        };
        let effect = predict_cds_effect("substitution", true, false, &cds_pos);
        assert_eq!(effect.name, "splice_region_variant");
        assert_eq!(effect.impact, "LOW");
    }

    #[test]
    fn test_predict_intron_effect() {
        let cds_pos = CdsPos {
            base: 117,
            offset: Some(-50),
            utr3: false,
            special: None,
        };
        let effect = predict_cds_effect("substitution", true, false, &cds_pos);
        assert_eq!(effect.name, "intron_variant");
        assert_eq!(effect.impact, "MODIFIER");
    }

    #[test]
    fn test_predict_5_utr_effect() {
        let cds_pos = CdsPos {
            base: -10,
            offset: None,
            utr3: false,
            special: None,
        };
        let effect = predict_cds_effect("substitution", false, false, &cds_pos);
        assert_eq!(effect.name, "5_prime_UTR_variant");
    }

    #[test]
    fn test_predict_3_utr_effect() {
        let cds_pos = CdsPos {
            base: 10,
            offset: None,
            utr3: true,
            special: None,
        };
        let effect = predict_cds_effect("substitution", false, false, &cds_pos);
        assert_eq!(effect.name, "3_prime_UTR_variant");
    }

    #[test]
    fn test_predict_frameshift_effect() {
        let cds_pos = CdsPos::new(350);
        let effect = predict_cds_effect("deletion", false, true, &cds_pos);
        assert_eq!(effect.name, "frameshift_variant");
        assert_eq!(effect.impact, "HIGH");
    }

    // -------------------------------------------------------------------------
    // span_len_from_cds_interval: special position guards (pter/qter/cen)
    // -------------------------------------------------------------------------

    /// A `pter` start endpoint makes the span uncomputable: `is_special()` must
    /// cause `span_len_from_cds_interval` to return `None`, not `Some(1)`.
    #[test]
    fn span_len_cds_interval_pter_start_returns_none() {
        use crate::hgvs::interval::CdsInterval;
        let interval = CdsInterval::new(CdsPos::pter(), CdsPos::new(50));
        assert_eq!(
            span_len_from_cds_interval(&interval),
            None,
            "pter start is special — span must be None"
        );
    }

    /// A `qter` end endpoint must also return `None`.
    #[test]
    fn span_len_cds_interval_qter_end_returns_none() {
        use crate::hgvs::interval::CdsInterval;
        let interval = CdsInterval::new(CdsPos::new(1), CdsPos::qter());
        assert_eq!(
            span_len_from_cds_interval(&interval),
            None,
            "qter end is special — span must be None"
        );
    }

    /// A `cen` position at either end must return `None`.
    #[test]
    fn span_len_cds_interval_cen_returns_none() {
        use crate::hgvs::interval::CdsInterval;
        let interval = CdsInterval::new(CdsPos::cen(), CdsPos::new(10));
        assert_eq!(
            span_len_from_cds_interval(&interval),
            None,
            "cen start is special — span must be None"
        );
    }

    /// Sanity check: a normal interval still computes correctly after the new
    /// guard is added (regression guard for the fix).
    #[test]
    fn span_len_cds_interval_normal_interval_unaffected() {
        use crate::hgvs::interval::CdsInterval;
        let interval = CdsInterval::new(CdsPos::new(10), CdsPos::new(19));
        assert_eq!(
            span_len_from_cds_interval(&interval),
            Some(10),
            "normal CDS interval should still compute span correctly"
        );
    }
}
