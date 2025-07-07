use anyhow::{Context, Result, bail};
use flate2::{read::MultiGzDecoder, write::GzEncoder, Compression as GzCompression};
use log::{debug, info}; // Added info
use sevenz_rust2::{Password, SevenZReader, SevenZWriter, Aes256Sha256Algorithm, CompressionMethod, FileInfo};
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Cursor, Read, Write}, // Added BufRead, Write, BufWriter, Read, Cursor
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
        Some("7z") => {
            info!("Reading 7z compressed file: {:?}", path);
            // sevenz-rust2 reads the whole entry into memory.
            // We find the first non-directory entry and read it.
            let mut sz_reader = SevenZReader::open(path, Password::empty())
                .with_context(|| format!("Failed to open 7z archive: {:?}", path))?;

            let mut entry_data: Option<Vec<u8>> = None;
            // sz_reader.for_each_entries(|entry, reader| { // API seems to have changed or my memory was off
            // Let's try to find the first file entry directly
            for (idx, entry) in sz_reader.archive().files.iter().enumerate() {
                if !entry.is_directory() {
                    let mut buffer = Vec::new();
                    sz_reader.read_entry_to_end(idx, &mut buffer)
                        .with_context(|| format!("Failed to read entry {} from 7z archive: {:?}", entry.name(), path))?;
                    entry_data = Some(buffer);
                    break; // Found the first file, stop.
                }
            }

            if let Some(data) = entry_data {
                Ok(Box::new(BufReader::new(Cursor::new(data))))
            } else {
                bail!("No suitable file entry found in 7z archive: {:?}", path)
            }
        }
        _ => {
            info!("Reading uncompressed file: {:?}", path);
            Ok(Box::new(BufReader::new(file)))
        }
    }
}

/// Opens a file for writing, handling compression based on file extension.
/// Supported extensions: .gz, .xz, .zst, .7z
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
        Some("7z") => {
            info!("Writing 7z compressed file: {:?}", path);
            // sevenz-rust2 write_entry takes a Read source.
            // So, we create a temporary buffer that implements Write.
            // On flush/drop, this buffer's content will be compressed into the 7z file.
            Ok(Box::new(SevenZCompressingWriter::new(file)?))
        }
        _ => {
            info!("Writing uncompressed file: {:?}", path);
            Ok(Box::new(BufWriter::new(file)))
        }
    }
}

struct SevenZCompressingWriter {
    target_file: Option<File>, // The actual output file for the .7z archive
    buffer: Cursor<Vec<u8>>, // In-memory buffer to store data before compression
    // We need path to give a name to the entry in 7z archive
    // entry_name: String, // Name for the entry within the 7z archive
}

impl SevenZCompressingWriter {
    fn new(target_file: File) -> Result<Self> {
        // let entry_name = path
        //     .file_name()
        //     .and_then(|n| n.to_str())
        //     .map(|s| s.trim_end_matches(".7z").to_string())
        //     .unwrap_or_else(|| "data.bin".to_string());

        Ok(Self {
            target_file: Some(target_file),
            buffer: Cursor::new(Vec::new()),
            // entry_name,
        })
    }

    fn Suffix(mut self) -> Result<()> { // Changed name to Suffix to avoid conflict with Write::flush
        if let Some(target_file) = self.target_file.take() {
            // Reset cursor position to the beginning to read its content
            self.buffer.set_position(0);

            // Use SevenZWriter to write the buffer to the target_file
            let mut sz_writer = SevenZWriter::new(target_file)
                .with_context(|| "Failed to create SevenZWriter")?;

            // Configure the entry - use a generic name or derive from path
            let entry_name = "content"; // Or derive from self.path if stored
            let entry_file_info = FileInfo::new(
                entry_name.into(),
                false, // is_directory
                None, // last_modified_date
                None, // attributes
                None, // created_date
                None, // last_accessed_date
                None, // crc (will be calculated)
                None, // size (will be calculated)
            );


            // CompressionMethod::Lzma2 is a good default
            sz_writer.push_entry(entry_file_info, Some(CompressionMethod::Lzma2(None)), &mut self.buffer)
                .with_context(|| "Failed to write entry to 7z archive")?;

            sz_writer.finish().with_context(|| "Failed to finish writing 7z archive")?;
            info!("Successfully wrote and finalized 7z archive.");
        }
        Ok(())
    }
}

impl Write for SevenZCompressingWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.buffer.write(buf)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.buffer.flush()
        // Actual compression happens on drop (or explicit finish)
    }
}

impl Drop for SevenZCompressingWriter {
    fn drop(&mut self) {
        // Ensure data is written on drop
        if self.target_file.is_some() { // Only if not already finalized
            if let Err(e) = self.Suffix() { // Renamed from flush to Suffix
                // Log error, but can't propagate it from drop.
                // Consider adding an explicit close/finish method to the public API
                // that users *should* call to handle errors properly.
                log::error!("Failed to write 7z archive on drop: {:?}", e);
            }
        }
    }
}
