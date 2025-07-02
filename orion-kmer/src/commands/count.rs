use anyhow::{Context, Result};
use dashmap::DashMap;
// use flate2::read::MultiGzDecoder; // Removed
use log::{debug, info};
use needletail::{Sequence, parse_fastx_file};
// use rayon::prelude::*;
use std::{
    fs::File,
    io::{BufWriter, Write}, // Removed BufRead, BufReader
    // path::Path, // Removed
    sync::atomic::{AtomicUsize, Ordering},
};

use crate::{
    cli::CountArgs,
    errors::OrionKmerError,
    kmer::{canonical_u64, seq_to_u64, u64_to_seq},
    utils::track_progress_and_resources, // Import the wrapper function
};
// use indicatif::ProgressBar; // Not needed if progress is per file

fn process_sequence_chunk(seq_chunk: &[u8], k: u8, kmer_counts: &DashMap<u64, AtomicUsize>) {
    if seq_chunk.len() < k as usize {
        return;
    }

    for window in seq_chunk.windows(k as usize) {
        if let Some(kmer_val) = seq_to_u64(window, k) {
            let canonical_kmer = canonical_u64(kmer_val, k);
            kmer_counts
                .entry(canonical_kmer)
                .or_insert_with(|| AtomicUsize::new(0))
                .fetch_add(1, Ordering::Relaxed);
        }
        // else: sequence contained 'N' or other non-ACGT char, skip this k-mer
    }
}

pub fn run_count(args: CountArgs) -> Result<()> {
    info!("Starting count command with args: {:?}", args);

    if args.kmer_size == 0 || args.kmer_size > 32 {
        return Err(OrionKmerError::InvalidKmerSize(args.kmer_size).into());
    }
    let k = args.kmer_size;

    let kmer_counts: DashMap<u64, AtomicUsize> = DashMap::new();
    let num_files = args.input_files.len() as u64;

    track_progress_and_resources("Counting k-mers from input files", num_files, |pb_files| {
        for input_path in &args.input_files {
            let path_str = input_path.to_string_lossy();
            info!("Processing file: {}", path_str);
            // Update progress bar message for the current file
            pb_files.set_message(format!("Processing: {}", path_str));

            let mut reader = parse_fastx_file(input_path)
                .with_context(|| format!("Failed to open or parse file: {}", path_str))?;

            info!("Processing records from {}...", path_str);
            let mut record_count = 0;
            while let Some(record) = reader.next() {
                let record =
                    record.with_context(|| format!("Error reading record from {}", path_str))?;
                let norm_seq = record.normalize(false);
                process_sequence_chunk(&norm_seq, k, &kmer_counts);
                record_count += 1;
                if record_count % 100_000 == 0 {
                    debug!("Processed {} records from {}", record_count, path_str);
                }
                // If we wanted a per-file progress bar based on records, we'd need total records first.
                // For now, the main progress bar is for files.
            }
            info!(
                "Finished processing {} records from {}. Unique k-mers so far: {}",
                record_count,
                path_str,
                kmer_counts.len()
            );
            pb_files.inc(1); // Increment file progress bar
        }
        Ok(())
    })?;

    info!(
        "Finished processing all input files. Found {} unique canonical k-mers.",
        kmer_counts.len()
    );

    // Outputting results
    debug!("Opening output file: {:?}", args.output_file);
    let output_file = File::create(&args.output_file)
        .with_context(|| format!("Failed to create output file: {:?}", args.output_file))?;
    let mut writer = BufWriter::new(output_file);

    let mut kmer_vec: Vec<(u64, usize)> = kmer_counts
        .into_iter()
        .filter_map(|(kmer_val, count_atomic)| {
            let count = count_atomic.into_inner();
            if count >= args.min_count {
                Some((kmer_val, count))
            } else {
                None
            }
        })
        .collect();

    // Sort for consistent output (optional, but good for testing)
    kmer_vec.sort_by_key(|item| item.0);

    info!(
        "Writing {} k-mers (count >= {}) to output file...",
        kmer_vec.len(),
        args.min_count
    );

    for (kmer_val, count) in kmer_vec {
        let kmer_seq_bytes = u64_to_seq(kmer_val, k);
        // This allocation to String can be slow for many k-mers.
        // Consider writing bytes directly if performance becomes an issue.
        let kmer_str = String::from_utf8(kmer_seq_bytes)
            .context("Failed to convert k-mer bytes to string (should not happen)")?;
        writeln!(writer, "{}\t{}", kmer_str, count)
            .context("Failed to write k-mer count to output file")?;
    }

    writer.flush().context("Failed to flush output writer")?;
    info!("Successfully wrote k-mer counts to {:?}", args.output_file);

    Ok(())
}
