use anyhow::Result;
use log::debug;

/// Determines the number of threads to use.
/// If `cli_threads` is 0, it uses all available logical cores.
/// Otherwise, it uses the number specified in `cli_threads`.
pub fn get_num_threads(cli_threads: usize) -> usize {
    let num_threads = if cli_threads == 0 {
        num_cpus::get()
    } else {
        cli_threads
    };
    debug!("Using {} threads for processing.", num_threads);
    num_threads
}

/// Initializes the Rayon global thread pool with the specified number of threads.
pub fn initialize_rayon_pool(num_threads: usize) -> Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()?;
    Ok(())
}
