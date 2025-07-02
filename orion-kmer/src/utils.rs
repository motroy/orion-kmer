use anyhow::{Context, Result};
use log::{debug, info}; // Added info
use std::{fs::File, io::BufReader, path::Path}; // Added Path and File, BufReader

use crate::db_types::KmerDbV2; // Import KmerDbV2

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

/// Loads a KmerDbV2 from the specified file path.
pub fn load_kmer_db_v2(path: &Path) -> Result<KmerDbV2> {
    info!("Loading k-mer database (KmerDbV2) from: {:?}", path);
    let file = File::open(path)
        .with_context(|| format!("Failed to open k-mer database file: {:?}", path))?;
    let reader = BufReader::new(file);
    let kmer_db: KmerDbV2 = bincode::deserialize_from(reader)
        .with_context(|| format!("Failed to deserialize KmerDbV2 from {:?}", path))?;

    info!(
        "Successfully loaded KmerDbV2 from {:?} (k={}, {} references, {} total unique k-mers)",
        path,
        kmer_db.k,
        kmer_db.num_references(),
        kmer_db.total_unique_kmers()
    );
    Ok(kmer_db)
}

use indicatif::{ProgressBar, ProgressStyle};
use psutil::process::Process;
use std::time::Instant;

/// Wraps a function to provide progress tracking, execution time, and max RAM usage.
pub fn track_progress_and_resources<F, R>(
    task_description: &str,
    total_items: u64,
    func: F,
) -> Result<R>
where
    F: FnOnce(&ProgressBar) -> Result<R>,
{
    info!("Starting task: {}", task_description);
    let start_time = Instant::now();

    let pb = ProgressBar::new(total_items);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
            )
            .unwrap_or_else(|e| {
                debug!("Error setting progress bar style: {}", e);
                ProgressStyle::default_bar()
            })
            .progress_chars("#>-"),
    );

    let result = func(&pb);

    pb.finish_with_message(format!("{} completed.", task_description));

    let duration = start_time.elapsed();
    info!("Task '{}' finished in {:.2?}", task_description, duration);

    match Process::current() {
        Ok(process) => match process.memory_info() {
            Ok(mem_info) => {
                info!(
                    "Max RAM usage for task '{}': {} MB",
                    task_description,
                    mem_info.rss() / 1024 / 1024
                );
            }
            Err(e) => {
                debug!("Failed to get memory info: {}", e);
            }
        },
        Err(e) => {
            debug!("Failed to get current process: {}", e);
        }
    }

    result
}
