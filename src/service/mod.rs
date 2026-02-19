//! Web service module for multi-tool HGVS variant normalization
//!
//! This module provides a web API that can run multiple HGVS normalization tools:
//! - ferro (native Rust, fastest)
//! - mutalyzer (HTTP API)
//! - biocommons/hgvs (Python subprocess)
//! - hgvs-rs (native Rust)
//!
//! The service provides both single variant and batch processing capabilities
//! with unified output formats and agreement analysis.

pub mod config;
pub mod handlers;
pub mod server;
pub mod tools;
pub mod types;
pub mod validation;

pub use config::ServiceConfig;
pub use server::{create_app, spawn_health_check_task};
pub use types::*;
