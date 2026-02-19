//! Ferro HGVS Web Service
//!
//! Multi-tool HGVS variant normalization web service that provides REST API
//! endpoints for parsing and normalizing HGVS variants using multiple tools:
//! - ferro (native Rust)
//! - mutalyzer (HTTP API)
//! - biocommons/hgvs (Python subprocess)
//! - hgvs-rs (native Rust)

use std::net::SocketAddr;
use std::path::{Path, PathBuf};

use clap::{Parser, Subcommand};
use tracing::{error, info};
use tracing_subscriber::util::SubscriberInitExt;

use ferro_hgvs::service::{create_app, spawn_health_check_task, ServiceConfig};

#[derive(Parser)]
#[command(name = "ferro-web")]
#[command(about = "Multi-tool HGVS variant normalization web service")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Start the web service
    Serve {
        /// Configuration file path
        #[arg(short, long, default_value = "config/service.toml")]
        config: PathBuf,

        /// Override host address
        #[arg(long)]
        host: Option<String>,

        /// Override port
        #[arg(short, long)]
        port: Option<u16>,

        /// Log level (trace, debug, info, warn, error)
        #[arg(long, default_value = "info")]
        log_level: String,

        /// Enable JSON logging
        #[arg(long)]
        json_logs: bool,

        /// Open browser automatically
        #[arg(long)]
        open: bool,
    },

    /// Generate a sample configuration file
    Config {
        /// Output path for configuration file
        #[arg(short, long, default_value = "config/service.toml")]
        output: PathBuf,

        /// Overwrite existing file
        #[arg(long)]
        force: bool,
    },

    /// Check configuration and tool availability
    Check {
        /// Configuration file path
        #[arg(short, long, default_value = "config/service.toml")]
        config: PathBuf,
    },
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Serve {
            config,
            host,
            port,
            log_level,
            json_logs,
            open,
        } => serve_command(config, host, port, log_level, json_logs, open).await,
        Commands::Config { output, force } => config_command(output, force).await,
        Commands::Check { config } => check_command(config).await,
    }
}

async fn serve_command(
    config_path: PathBuf,
    host_override: Option<String>,
    port_override: Option<u16>,
    log_level: String,
    json_logs: bool,
    open_browser: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    // Initialize tracing
    init_tracing(&log_level, json_logs)?;

    info!("Starting ferro-hgvs web service");

    // Load configuration
    let mut config = load_or_create_config(&config_path).await?;

    // Apply overrides
    if let Some(host) = host_override {
        config.server.host = host;
    }
    if let Some(port) = port_override {
        config.server.port = port;
    }

    // Validate configuration
    if let Err(e) = config.validate() {
        error!("Configuration validation failed: {}", e);
        return Err(e.into());
    }

    info!("Configuration loaded successfully");
    info!("Enabled tools: {:?}", config.enabled_tools());

    // Create the application
    let (app, state) = create_app(config.clone())?;

    // Start background health check task
    spawn_health_check_task(state);

    // Create server address
    let addr = SocketAddr::new(config.server.host.parse()?, config.server.port);

    info!("Starting server on http://{}", addr);

    // Start the server
    let listener = tokio::net::TcpListener::bind(&addr).await?;

    // Determine the URL to use for browser (use localhost if bound to 0.0.0.0)
    let browser_host = if config.server.host == "0.0.0.0" {
        "localhost"
    } else {
        &config.server.host
    };
    let url = format!("http://{}:{}", browser_host, config.server.port);

    // Use emoji only in debug builds
    #[cfg(debug_assertions)]
    {
        info!("Ferro HGVS web service running on {}", url);
        info!("API documentation available at {}/", url);
        info!("Health check available at {}/health", url);
    }
    #[cfg(not(debug_assertions))]
    {
        info!("Ferro HGVS web service running on {}", url);
        info!("API documentation available at {}/", url);
        info!("Health check available at {}/health", url);
    }

    // Open browser if requested
    if open_browser {
        if let Err(e) = open::that(&url) {
            error!("Failed to open browser: {}", e);
        }
    }

    axum::serve(listener, app).await?;

    Ok(())
}

async fn config_command(
    output_path: PathBuf,
    force: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    if output_path.exists() && !force {
        error!(
            "Configuration file already exists: {}",
            output_path.display()
        );
        error!("Use --force to overwrite");
        std::process::exit(1);
    }

    // Create output directory if it doesn't exist
    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    // Generate default configuration
    let config = ServiceConfig::default();

    config.to_file(&output_path)?;

    println!(
        "Sample configuration file created: {}",
        output_path.display()
    );
    println!("Edit the file to configure your tools and settings");

    Ok(())
}

async fn check_command(config_path: PathBuf) -> Result<(), Box<dyn std::error::Error>> {
    println!("Checking configuration and tool availability...");

    // Load configuration
    let config = load_or_create_config(&config_path).await?;

    // Validate configuration
    match config.validate() {
        Ok(()) => println!("Configuration is valid"),
        Err(e) => {
            println!("Configuration validation failed: {}", e);
            return Err(e.into());
        }
    }

    // Check tool availability
    println!("\nTool Status:");

    let enabled_tools = config.enabled_tools();
    if enabled_tools.is_empty() {
        println!("WARNING: No tools are enabled");
        return Ok(());
    }

    // Try to create the application (this will initialize tools)
    match create_app(config) {
        Ok((_, state)) => {
            println!("Application created successfully");

            // Run health checks
            let health_results = state.tool_manager.health_check_all().await;

            for (tool_name, result) in health_results {
                use ferro_hgvs::service::types::health_check::HealthCheckResult;
                match result {
                    HealthCheckResult::Healthy => println!("  OK {}: Available", tool_name),
                    HealthCheckResult::Degraded { reason } => {
                        println!("  WARN {}: Degraded - {}", tool_name, reason)
                    }
                    HealthCheckResult::Unhealthy { reason } => {
                        println!("  ERROR {}: {}", tool_name, reason)
                    }
                }
            }
        }
        Err(e) => {
            println!("Failed to initialize application: {}", e);
            return Err(e.into());
        }
    }

    println!("\nHealth check completed");
    Ok(())
}

async fn load_or_create_config(
    config_path: &Path,
) -> Result<ServiceConfig, Box<dyn std::error::Error>> {
    if config_path.exists() {
        info!("Loading configuration from {}", config_path.display());
        Ok(ServiceConfig::from_file(config_path)?)
    } else {
        info!("Configuration file not found, using defaults");
        println!(
            "WARNING: Configuration file not found: {}",
            config_path.display()
        );
        println!("TIP: Run 'ferro-web config' to generate a sample configuration file");

        // Return default config but warn user
        Ok(ServiceConfig::default())
    }
}

fn init_tracing(level: &str, _json_logs: bool) -> Result<(), Box<dyn std::error::Error>> {
    use tracing_subscriber::{fmt, layer::SubscriberExt, EnvFilter};

    // Parse and validate log level
    let filter =
        EnvFilter::try_new(level).map_err(|e| format!("Invalid log level '{}': {}", level, e))?;

    // Initialize with standard formatting
    // Note: JSON logging requires the 'json' feature in tracing-subscriber
    tracing_subscriber::registry()
        .with(filter)
        .with(fmt::layer())
        .init();

    info!("Tracing initialized with level: {}", level);

    Ok(())
}
