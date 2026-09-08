//! The partitioner registry: an open `name -> arm` map that replaces the closed
//! `FERRO_PARTITION` enum (design §8).
//!
//! Every stable legacy name (`live`, `shadow`, `canonical`, `canonical-coalesced`)
//! is registered here, each backed by the existing `merge::partition_block*`
//! dispatch wrapped as a [`Strategy::Mono`] — verified byte-identical to today's
//! in-module dispatch by `tests/it/partition_registry_equivalence.rs`, not
//! assumed (design §9 tripwire). Task 4 adds the bake-off arms alongside, and the
//! ruled driver is registered under `ruled` ([`DERIVED_BLOCK_PARTITION_KEY`]) as
//! [`RuledArm`].
//!
//! # Dependency direction
//!
//! This module reaches back into `crate::normalize::merge` for one reason: the
//! four legacy wrappers ([`LegacyRuleArm`]) and the ruled arm ([`RuledArm`]) need
//! `merge`-private `partition_block*`/`Piece` internals, so they are defined there
//! and merely *registered* here. It is one of the confined `partition → normalize` seams the
//! [`crate::partition`] module doc enumerates.

use std::collections::HashMap;
use std::sync::OnceLock;

use crate::normalize::merge::{legacy_rule_arms, LegacyRuleArm, RuledArm};
use crate::partition::arm::DialBConfig;
use crate::partition::arms::arms::{AlwaysOneDelinsL2, TrimOnlyL1};
use crate::partition::arms::r2_winner;
use crate::partition::partitioner::Partitioner;
use crate::partition::strategy::Strategy;

/// The registry-key counterpart of `merge::DERIVED_BLOCK_PARTITION_RULE` — the
/// name the derivation surface pins to (design §8).
///
/// The two-surface pin (#1834/#2155) is stated once more in the arm-*name*
/// vocabulary a registry selects from, so the property survives the move from a
/// closed enum to an open map: `canonicalize_from_sequence` resolves its arm at
/// runtime while `derive_block_members` stays pinned, and the two must agree.
/// Kept in lockstep with the enum pin by
/// `merge::partition_rule_knob::the_two_surfaces_agree_on_a_pinned_registry_key`,
/// and registered: [`registry`] maps it to [`RuledArm`].
pub const DERIVED_BLOCK_PARTITION_KEY: &str = "ruled";

/// Wrap a legacy-rule arm as a named [`Strategy::Mono`].
///
/// `Strategy::Mono` already dispatches to `Monolithic::partition` (Fable review
/// #9), so registration is mechanical: only `Strategy` carries the blanket
/// `Partitioner` impl, so the arm is boxed inside a `Strategy`, never registered
/// as a bare `Monolithic`.
fn mono(arm: LegacyRuleArm) -> Strategy {
    Strategy::Mono {
        name: crate::partition::arm::Monolithic::name(&arm).to_owned(),
        arm: Box::new(arm),
    }
}

/// The process-wide registry: every stable arm name mapped to its partitioner.
///
/// `&'static` and built once, like the arm selection itself. The values are
/// `Box<dyn Partitioner + Send + Sync>` (a `Strategy`, or [`RuledArm`]): design §8
/// sketches `+ Sync`, but a value held in a `static OnceLock` must be `Send + Sync`
/// (a `static` must be `Sync`, and `OnceLock<T>: Sync` needs `T: Send + Sync`), so
/// the bound is tightened. Every arm is a plain data struct, so `Send` is free.
pub fn registry() -> &'static HashMap<&'static str, Box<dyn Partitioner + Send + Sync>> {
    static REGISTRY: OnceLock<HashMap<&'static str, Box<dyn Partitioner + Send + Sync>>> =
        OnceLock::new();
    REGISTRY.get_or_init(|| {
        // The keys are the arms' own stable `&'static str` names, taken from the
        // same source of truth `merge` builds the arms from, so a name cannot
        // drift between the key and the arm.
        let names = ["live", "shadow", "canonical", "canonical-coalesced"];
        let mut map: HashMap<&'static str, Box<dyn Partitioner + Send + Sync>> = HashMap::new();
        for (name, arm) in names.into_iter().zip(legacy_rule_arms()) {
            debug_assert_eq!(
                crate::partition::arm::Monolithic::name(&arm),
                name,
                "the registry key must match the arm's own name",
            );
            map.insert(name, Box::new(mono(arm)));
        }
        // Task 4: the bake-off arms, NON-default — reachable only by name, never
        // selected unless a caller asks for one. Each is a `Strategy` (the blanket
        // `Partitioner` impl), so registration is just a boxed insert. Unlike the
        // four legacy wrappers, these do NOT decline-to-`live`: an unsound partition
        // surfaces as `Err(PartitionError::Unsound)` through `Partitioner::partition`
        // (the hard-error rule for the bake-off arms).
        for arm in bakeoff_arms() {
            // The map key is the arm's own stable `&'static str` name, so a name
            // cannot drift between key and value. `Strategy::name` returns a `&str`
            // into the arm's `String`; the literal key is asserted to match it.
            let name = bakeoff_arm_key(arm.name());
            debug_assert_eq!(
                arm.name(),
                name,
                "the registry key must match the bake-off arm's own name",
            );
            map.insert(name, Box::new(arm));
        }
        // The ruled driver, under the name the shipped default and the derivation
        // surface pin to, so `arm(DERIVED_BLOCK_PARTITION_KEY)` resolves.
        debug_assert_eq!(RuledArm::NAME, DERIVED_BLOCK_PARTITION_KEY);
        map.insert(RuledArm::NAME, Box::new(RuledArm));
        map
    })
}

/// The `&'static str` key for a bake-off arm whose `Strategy` name is `name`, so the
/// map can be keyed by `&'static str` while the `Strategy` owns a `String`. Every
/// registered bake-off name is a literal here; an unrecognised name panics rather
/// than leak a `String::leak`, keeping the key set closed and auditable.
fn bakeoff_arm_key(name: &str) -> &'static str {
    match name {
        "trim-only/op-extract-v2/ledger-r5" => "trim-only/op-extract-v2/ledger-r5",
        "axis-peel/op-extract-v2/ledger-r5" => "axis-peel/op-extract-v2/ledger-r5",
        "peel-comp/op-extract-v2/ledger-r5" => "peel-comp/op-extract-v2/ledger-r5",
        "peel-strict/op-extract-v2/ledger-r5" => "peel-strict/op-extract-v2/ledger-r5",
        "some-l1/always-one-delins/off" => "some-l1/always-one-delins/off",
        other => panic!("unregistered bake-off arm name: {other:?}"),
    }
}

/// The bake-off arms registered alongside the four legacy names (Task 4).
///
/// The four R2-winner arms (`{trim-only,axis-peel,peel-comp,peel-strict}/
/// op-extract-v2/ledger-r5`) are the graded tournament survivors, taken straight
/// from [`r2_winner`] so the registered arm is the exact one the confluence harness
/// measured. The fifth, `some-l1/always-one-delins/off`, is the L2 control (always
/// one `delins` over the trimmed block) — the arm that proves the Direction-2
/// adapter honors a DECLARED kind: on a reverse-complement-shaped block the geometry
/// ladder would promote to `inv`, and this arm's `Delins` choice must survive.
fn bakeoff_arms() -> Vec<Strategy> {
    let mut arms = r2_winner::all();
    arms.push(Strategy::Grid {
        name: "some-l1/always-one-delins/off".into(),
        l1: Box::new(TrimOnlyL1),
        l2: Box::new(AlwaysOneDelinsL2),
        dial_b: DialBConfig::all_off(),
    });
    arms
}

/// Look up a registered arm by name.
pub fn arm(name: &str) -> Option<&'static (dyn Partitioner + Send + Sync)> {
    registry().get(name).map(|arm| arm.as_ref())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::partition::block_ctx::{BlockCtx, FrameContext, Molecule, Provenance};

    /// The four stable legacy names are all registered and all resolve.
    #[test]
    fn the_four_legacy_names_are_registered() {
        for name in ["live", "shadow", "canonical", "canonical-coalesced"] {
            assert!(arm(name).is_some(), "missing `{name}`");
        }
        assert!(arm("nope").is_none());
        assert!(registry().len() >= 4);
    }

    /// The Task 4 bake-off arms all resolve by name and are NON-default (reachable
    /// only via an explicit lookup).
    #[test]
    fn the_bakeoff_arms_are_registered_and_non_default() {
        for name in [
            "trim-only/op-extract-v2/ledger-r5",
            "axis-peel/op-extract-v2/ledger-r5",
            "peel-comp/op-extract-v2/ledger-r5",
            "peel-strict/op-extract-v2/ledger-r5",
            "some-l1/always-one-delins/off",
        ] {
            let a = arm(name).unwrap_or_else(|| panic!("missing bake-off arm `{name}`"));
            assert_eq!(a.name(), name, "the registry key must match the arm's name");
        }
        // The four legacy names plus the five bake-off arms.
        assert!(registry().len() >= 9);
    }

    /// The derivation surface's pinned key names a registered arm, and that arm is
    /// the ruled driver: its members cover exactly the runs `partition_ruled` cuts
    /// on the 3' placement, and it declares the canonical family. Before `ruled` was
    /// registered, `normalize_block_via_arm` on the pinned key returned
    /// `Unrenderable` for an unknown arm name.
    #[test]
    fn the_derived_block_partition_key_resolves_to_the_ruled_driver() {
        use crate::normalize::{normalize_block_via_arm, ShuffleDirection};
        use crate::partition::driver::partition_ruled;

        let a = arm(DERIVED_BLOCK_PARTITION_KEY)
            .unwrap_or_else(|| panic!("registry lacks `{DERIVED_BLOCK_PARTITION_KEY}`"));
        assert_eq!(a.name(), DERIVED_BLOCK_PARTITION_KEY);
        assert!(a.cuts_with_canonical(), "the ruled arm is canonical-family");

        let frame = FrameContext::NonCoding;
        let provenance = Provenance::none();
        for (reference, resulting) in [
            (&b"ACGT"[..], &b"AGT"[..]),
            (&b"ACGTACGT"[..], &b"ACCTTCGT"[..]),
            (&b"AACCGGTTA"[..], &b"AAACCGGTA"[..]),
        ] {
            let ctx = BlockCtx {
                reference,
                resulting,
                frame: &frame,
                molecule: Molecule::Dna,
                provenance: &provenance,
            };
            let cut = partition_ruled(&ctx, ShuffleDirection::ThreePrime)
                .expect("the ruled driver partitions this block");
            let partition = a.partition(&ctx).expect("the ruled arm is Ok");
            let member_spans: Vec<(usize, usize)> = partition
                .members
                .iter()
                .map(|m| (m.ref_start, m.ref_end))
                .collect();
            let run_spans: Vec<(usize, usize)> = cut
                .runs()
                .iter()
                .map(|r| (r.ref_start, r.ref_end))
                .collect();
            assert_eq!(member_spans, run_spans);
            assert!(
                normalize_block_via_arm(DERIVED_BLOCK_PARTITION_KEY, reference, resulting).is_ok(),
                "`normalize_block_via_arm` must resolve and render the pinned key"
            );
        }
        // A lone deletion renders as one, not as a payload-less `delins` (the
        // adapter brackets every partition).
        assert_eq!(
            normalize_block_via_arm(DERIVED_BLOCK_PARTITION_KEY, b"ACGT", b"AGT").as_deref(),
            Ok("g.[2del]")
        );
    }

    /// Each registered arm's declared coalesce-eligibility matches its name, and
    /// a legacy arm is always `Ok` on a well-formed block (Q1).
    #[test]
    fn legacy_arms_are_sound_and_declare_their_family() {
        let frame = FrameContext::NonCoding;
        let provenance = Provenance::none();
        let ctx = BlockCtx {
            reference: b"ACGT",
            resulting: b"AGT",
            frame: &frame,
            molecule: Molecule::Dna,
            provenance: &provenance,
        };
        for (name, cuts) in [
            ("live", false),
            ("shadow", false),
            ("canonical", true),
            ("canonical-coalesced", true),
        ] {
            let a = arm(name).unwrap();
            assert_eq!(a.cuts_with_canonical(), cuts, "`{name}` family");
            assert!(a.partition(&ctx).is_ok(), "`{name}` must be Ok");
        }
    }
}
