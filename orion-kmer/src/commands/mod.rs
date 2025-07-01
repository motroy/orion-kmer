pub mod build;
pub mod compare;
pub mod count;
pub mod query;

use crate::cli::Commands; // Keep only Commands
use anyhow::Result;

pub fn dispatch_command(command: Commands, threads: usize, verbose: u8) -> Result<()> {
    // Setup logging based on verbosity
    let log_level = match verbose {
        0 => log::LevelFilter::Warn,
        1 => log::LevelFilter::Info,
        2 => log::LevelFilter::Debug,
        _ => log::LevelFilter::Trace,
    };
    env_logger::Builder::new().filter_level(log_level).init();

    // Initialize rayon thread pool
    crate::utils::initialize_rayon_pool(crate::utils::get_num_threads(threads))?;

    match command {
        Commands::Count(args) => count::run_count(args),
        Commands::Build(args) => build::run_build(args),
        Commands::Compare(args) => compare::run_compare(args),
        Commands::Query(args) => query::run_query(args),
    }
}
