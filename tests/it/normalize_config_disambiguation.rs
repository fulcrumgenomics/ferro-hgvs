//! The batch-run config and the normalizer config must stay distinguishable.
//!
//! Two public types were both called `NormalizeConfig`:
//! [`ferro_hgvs::commands::NormalizeCommandConfig`] (reference directory, progress,
//! workers — the shape of a batch *run*) and [`ferro_hgvs::normalize::NormalizeConfig`]
//! (shuffle direction, window, and the **error mode**). With one name between them,
//! `NormalizeConfig::default()` resolved to whichever a single `use` line happened to
//! bring into scope, and the two defaults differ in what they even contain — so the
//! mistake produced a config silently missing an error mode rather than a compile error.
//!
//! This file is the executable half of that rename. Both types are imported here **by
//! their own names, in one scope**, so the rename cannot be undone silently: reverting it
//! breaks this file's compilation — as `E0432: unresolved import` when the new name simply
//! disappears (verified), or as a duplicate-import error if both types end up sharing a
//! name again. Either way the build stops. *That* is the guard; the assertions below only
//! document what each type carries.
//!
//! The `const` below pins `normalize_batch`'s signature for the same reason. A test that
//! merely *called* the function would need a temp file and would exercise batch I/O this
//! test has no interest in; coercing the function to a concrete `fn` pointer checks the
//! parameter type at compile time and runs nothing.

use ferro_hgvs::commands::{BatchResults, NormalizeCommandConfig};
use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::normalize::NormalizeConfig;
use ferro_hgvs::FerroError;

/// `normalize_batch` must keep taking the *command* config, not the normalizer one.
///
/// Instantiated at the owned `PathBuf` rather than `&Path` deliberately: naming a
/// reference type here would pin `P` to one concrete lifetime, and the resulting `fn`
/// item is then less general than the `for<'a>` pointer this constant declares.
const _NORMALIZE_BATCH_TAKES_THE_COMMAND_CONFIG: fn(
    std::path::PathBuf,
    &NormalizeCommandConfig,
) -> Result<BatchResults, FerroError> = ferro_hgvs::commands::normalize_batch::<std::path::PathBuf>;

#[test]
fn the_batch_config_and_the_normalizer_config_are_distinct_types() {
    let batch = NormalizeCommandConfig::default();
    let normalizer = NormalizeConfig::default();

    // Naming each config's own fields is what separates the two types here. Their
    // default *values* are deliberately not pinned: nothing in this change claims a
    // particular worker count, and a test called "are distinct types" going red
    // because someone retuned a default would point at the wrong thing.
    let _run_shape = (batch.workers, batch.show_progress, &batch.reference_dir);

    // The normalizer config is the one that carries the error mode — the field the
    // old shared name made it possible to lose by accident. This value *is* pinned,
    // because `src/normalize/config.rs` states it: `Default` and
    // `NormalizeConfig::lenient` both give it `ErrorConfig::lenient()`.
    assert_eq!(normalizer.error_config.mode, ErrorMode::Lenient);
}

/// `NormalizeConfig::lenient()` and `NormalizeConfig::default()` agree on the error mode.
///
/// Named in the doc comment the rename added to `src/normalize/config.rs`, so it is worth
/// having something check it rather than leaving it as prose.
#[test]
fn the_normalizer_config_is_lenient_by_default_and_by_constructor() {
    assert_eq!(
        NormalizeConfig::default().error_config.mode,
        NormalizeConfig::lenient().error_config.mode
    );
    assert_eq!(
        NormalizeConfig::lenient().error_config.mode,
        ErrorMode::Lenient
    );
}
