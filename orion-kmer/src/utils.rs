use anyhow::{Context, Result};
use flate2::{read::MultiGzDecoder, write::GzEncoder, Compression as GzCompression};
use log::{debug, info}; // Added info
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read, Write}, // Added BufRead, Write, BufWriter, Read
    path::Path,
};
use xz2::{read::XzDecoder, write::XzEncoder};
use zstd::{stream::read::Decoder as ZstdDecoder, stream::write::Encoder as ZstdEncoder};

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
/// Handles decompression automatically based on file extension.
pub fn load_kmer_db_v2(path: &Path) -> Result<KmerDbV2> {
    info!("Loading k-mer database (KmerDbV2) from: {:?}", path);
    // Use get_input_reader to handle potential compression
    let mut reader = get_input_reader(path)
        .with_context(|| format!("Failed to get input reader for k-mer database: {:?}", path))?;

    // bincode::deserialize_from directly takes a Read, which Box<dyn BufRead> implements.
    let kmer_db: KmerDbV2 = bincode::deserialize_from(&mut reader)
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

// Helper function to get file extension as lowercase string
fn get_extension(path: &Path) -> Option<String> {
    path.extension()
        .and_then(|ext| ext.to_str())
        .map(|s| s.to_lowercase())
}

/// Opens a file for reading, handling decompression based on file extension.
/// Supported extensions: .gz, .xz, .zst.
/// Returns a `Box<dyn BufRead>` for generic reading.
pub fn get_input_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open input file: {:?}", path))?;
    let extension = get_extension(path);

    match extension.as_deref() {
        Some("gz") => {
            info!("Reading GZipped file: {:?}", path);
            let decoder = MultiGzDecoder::new(file);
            Ok(Box::new(BufReader::new(decoder)))
        }
        Some("xz") => {
            info!("Reading XZ compressed file: {:?}", path);
            let decoder = XzDecoder::new(file);
            Ok(Box::new(BufReader::new(decoder)))
        }
        Some("zst") | Some("zstd") => {
            info!("Reading Zstandard compressed file: {:?}", path);
            let decoder = ZstdDecoder::new(file)
                .with_context(|| format!("Failed to create ZstdDecoder for {:?}", path))?;
            Ok(Box::new(BufReader::new(decoder)))
        }
        _ => {
            info!("Reading uncompressed file: {:?}", path);
            Ok(Box::new(BufReader::new(file)))
        }
    }
}

/// Opens a file for writing, handling compression based on file extension.
/// Supported extensions: .gz, .xz, .zst.
/// Returns a `Box<dyn Write>` for generic writing.
/// Note: The writers are typically `BufWriter`s wrapping compressing encoders.
pub fn get_output_writer(path: &Path) -> Result<Box<dyn Write>> {
    let file = File::create(path)
        .with_context(|| format!("Failed to create output file: {:?}", path))?;
    let extension = get_extension(path);

    match extension.as_deref() {
        Some("gz") => {
            info!("Writing GZipped file: {:?}", path);
            // BufWriter is recommended by flate2 for performance.
            // The GzEncoder itself is not necessarily buffered internally in the way BufWriter is.
            let encoder = GzEncoder::new(file, GzCompression::default());
            Ok(Box::new(BufWriter::new(encoder)))
        }
        Some("xz") => {
            info!("Writing XZ compressed file: {:?}", path);
            // XzEncoder is buffered, but wrapping in BufWriter is harmless and consistent.
            let encoder = XzEncoder::new(file, 6); // Compression level 6 is a good default
            Ok(Box::new(BufWriter::new(encoder)))
        }
        Some("zst") | Some("zstd") => {
            info!("Writing Zstandard compressed file: {:?}", path);
            // ZstdEncoder benefits from a BufWriter.
            let encoder = ZstdEncoder::new(file, 0) // 0 is default compression level for zstd crate
                .with_context(|| format!("Failed to create ZstdEncoder for {:?}", path))?
                .auto_finish(); // Ensures finish is called on drop
            Ok(Box::new(BufWriter::new(encoder)))
        }
        _ => {
            info!("Writing uncompressed file: {:?}", path);
            Ok(Box::new(BufWriter::new(file)))
        }
    }
}
