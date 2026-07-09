//! HGVS arbitration: explain and adjudicate parse/normalize/projection
//! disagreements between ferro and another tool. See the design spec
//! `specs/2026-07-09-hgvs-arbitration-skill-design.md`.

pub mod category;

pub use category::ArbitrationCategory;
