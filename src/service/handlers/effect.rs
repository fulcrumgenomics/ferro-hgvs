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
            let (edit_type, is_frameshift, ref_len, alt_len) =
                if let Some(edit) = v.loc_edit.edit.inner() {
                    analyze_na_edit(edit, span_len)
                } else {
                    ("unknown", false, 0, 0)
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
                    edit_type,
                    is_frameshift,
                    ref_len,
                    alt_len,
                )
            } else {
                None
            };

            // NMD prediction if requested. Skip when the span is
            // undecidable — an NMD call rests on the frameshift signal,
            // which is unknown for a conservative-skip del/dup/ins.
            let nmd_prediction = if include_nmd && !span_undecidable {
                predict_nmd_for_cds(state, &accession, &cds_pos, is_frameshift, edit_type)
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
    edit_type: &str,
    is_frameshift: bool,
    ref_len: usize,
    _alt_len: usize,
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

    // `is_frameshift` is authoritative from the classifier; the
    // `edit_type != "substitution"` guard is belt-and-braces (the
    // classifier already returns `false` for the substitution arms).
    // `_alt_len` is unused for the frameshift decision after this fix
    // — the classifier handles unknown-length insert shapes correctly;
    // kept in the signature for symmetry with `ref_len`.
    let is_frameshift = is_frameshift && edit_type != "substitution";

    // Get the protein accession: the authoritative cdot value if present, else
    // the transcript accession itself. We do NOT infer `NP_*`/`XP_*` from
    // `NM_*`/`XM_*` by preserving the number — RefSeq does not guarantee the NM
    // and NP numbers match, so that inference is frequently wrong (#808).
    let prot_acc = cdot_tx
        .protein
        .clone()
        .unwrap_or_else(|| accession.to_string());

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

    // Note: Full amino acid lookup requires CDS sequence data
    // which is not currently available in the service.
    // The "?" markers indicate uncertain/unknown amino acids.
    Some(ProteinConsequence {
        hgvs_p,
        ref_aa: format!("?(pos{})", codon_phase + 1), // Show codon position (1-3)
        alt_aa: "?".to_string(),
        position: prot_position,
        is_frameshift,
    })
}

/// Predict NMD for CDS variant
fn predict_nmd_for_cds(
    state: &AppState,
    accession: &str,
    cds_pos: &CdsPos,
    is_frameshift: bool,
    edit_type: &str,
) -> Option<NmdPrediction> {
    // NMD prediction requires:
    // 1. Knowing if the variant introduces a premature termination codon (PTC)
    // 2. The position of the PTC relative to the last exon-exon junction

    // For now, we use simplified rules:
    // - Frameshift variants early in the CDS are more likely to trigger NMD
    // - Variants in the last exon typically escape NMD

    let cdot = state.cdot.as_ref()?;
    let cdot_tx = cdot.get_transcript(accession)?;

    let cds_start = cdot_tx.cds_start?;
    let cds_end = cdot_tx.cds_end?;
    let cds_length = cds_end.saturating_sub(cds_start);

    // Calculate relative position in CDS
    let relative_pos = if cds_pos.base > 0 && cds_length > 0 {
        cds_pos.base as f64 / cds_length as f64
    } else {
        0.0
    };

    // Check if in last exon (simplified - last 10% of CDS rule)
    let near_3_end = relative_pos > 0.9;

    // Predict NMD likelihood
    let (predicted, confidence, reason) = if !is_frameshift && edit_type != "deletion" {
        (
            false,
            0.9,
            "Non-frameshift variant unlikely to trigger NMD".to_string(),
        )
    } else if near_3_end {
        (
            false,
            0.8,
            "Variant in last exon region - likely escapes NMD".to_string(),
        )
    } else if relative_pos < 0.5 && is_frameshift {
        (
            true,
            0.7,
            "Early frameshift likely to introduce PTC and trigger NMD".to_string(),
        )
    } else if is_frameshift {
        (
            true,
            0.5,
            "Frameshift may introduce PTC, NMD possible".to_string(),
        )
    } else {
        (
            false,
            0.3,
            "Insufficient information for NMD prediction".to_string(),
        )
    };

    Some(NmdPrediction {
        predicted,
        confidence,
        reason,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

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
