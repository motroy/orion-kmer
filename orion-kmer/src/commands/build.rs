use anyhow::{Context, Result};
use dashmap::DashSet; // Using DashSet for concurrent k-mer collection per file
use log::{debug, info};
use needletail::{parse_fastx_reader, Sequence}; // Corrected import order & added parse_fastx_reader
use std::{
    collections::HashSet, // Keep HashSet for final storage in KmerDbV2
    // fs::File, // No longer directly used for output
    // io::{BufWriter, Write}, // No longer directly used for output
    path::PathBuf, // For getting filename
};

use crate::{
    cli::BuildArgs,
    db_types::KmerDbV2, // Import the new database structure
    errors::OrionKmerError,
    kmer::{canonical_u64, seq_to_u64},
    utils::{get_buffered_file_reader, get_output_writer, track_progress_and_resources}, // Import the wrapper function and I/O helpers
};
// use indicatif::ProgressBar; // Required for the closure signature - actually not needed

// This function processes sequences for a single file and populates a DashSet for that file.
// It now accepts a ProgressBar to update progress within the file processing.
fn process_sequences_for_file(
    file_path: &PathBuf,
    k: u8,
    file_kmer_set: &DashSet<u64>,
    // pb: &ProgressBar, // Optional: if we want fine-grained progress per file
) -> Result<()> {
    let path_str = file_path.to_string_lossy();
    info!("Processing genome file: {}", path_str);

    // Estimate total records for a more accurate inner progress bar if needed
    // This is a placeholder; actual record count might require a preliminary scan or be omitted
    // let total_records = ...;
    // pb.set_length(total_records); // If using a per-file progress bar

    // Use get_buffered_file_reader, needletail will handle decompression
    let input_buf_reader = get_buffered_file_reader(file_path)
        .with_context(|| format!("Failed to get buffered file reader for file: {}", path_str))?;

    // Pass the BufRead to parse_fastx_reader instead of a path to parse_fastx_file
    let mut reader = parse_fastx_reader(input_buf_reader)
        .with_context(|| format!("Failed to parse FASTA/Q content from: {}", path_str))?;

    let mut record_count = 0;
    while let Some(record) = reader.next() {
        let record = record.with_context(|| format!("Error reading record from {}", path_str))?;
        let norm_seq = record.normalize(false); // Ensure uppercase, no N conversion yet

        if norm_seq.len() >= k as usize {
            for window in norm_seq.windows(k as usize) {
                // seq_to_u64 handles 'N' by returning None
                if let Some(kmer_val) = seq_to_u64(window, k) {
                    let canonical_kmer = canonical_u64(kmer_val, k);
                    file_kmer_set.insert(canonical_kmer);
                }
            }
        }
        record_count += 1;
        if record_count % 100_000 == 0 {
            // Keep debug logging for detailed progress
            debug!(
                "Processed {} records from {}. Current unique k-mers for this file: {}",
                record_count,
                path_str,
                file_kmer_set.len()
            );
        }
        // pb.inc(1); // Increment per-file progress bar if used
    }
    info!(
        "Finished processing {} records from {}. Found {} unique k-mers for this file.",
        record_count,
        path_str,
        file_kmer_set.len()
    );
    Ok(())
}

pub fn run_build(args: BuildArgs) -> Result<()> {
    info!("Starting build command with args: {:?}", args);

    if args.kmer_size == 0 || args.kmer_size > 32 {
        return Err(OrionKmerError::InvalidKmerSize(args.kmer_size).into());
    }
    let k = args.kmer_size;

    let mut kmer_db_v2 = KmerDbV2::new(k);
    let num_files = args.genome_files.len() as u64;

    // Wrap the main file processing loop
    track_progress_and_resources("Building k-mer database", num_files, |pb_files| {
        for input_path in &args.genome_files {
            // For each file, create a new DashSet to collect its k-mers.
            let file_kmer_set: DashSet<u64> = DashSet::new();

            // Pass the main progress bar `pb_files` if process_sequences_for_file
            // is to update it directly (e.g. if it was for sequences, not files).
            // Here, we are processing file by file, so `pb_files.inc(1)` is done after each file.
            process_sequences_for_file(input_path, k, &file_kmer_set)?;
            // Consider adding a nested progress bar inside process_sequences_for_file
            // if individual file processing is very long and has measurable units (e.g. sequences).

            let final_file_kmers: HashSet<u64> = file_kmer_set.into_iter().collect();

            let reference_name = input_path.file_name().map_or_else(
                || input_path.to_string_lossy().into_owned(),
                |os_str| os_str.to_string_lossy().into_owned(),
            );

            info!(
                "Adding {} unique k-mers from reference '{}' to the database.",
                final_file_kmers.len(),
                reference_name
            );
            kmer_db_v2.add_reference(reference_name.clone(), final_file_kmers); // Use clone if reference_name is used after
            pb_files.set_message(format!("Processed: {}", reference_name));
            pb_files.inc(1);
        }
        Ok(()) // Return Ok from the closure
    })?;

    info!(
        "Finished processing all input files. Database contains {} references and a total of {} unique canonical k-mers across all references.",
        kmer_db_v2.num_references(),
        kmer_db_v2.total_unique_kmers()
    );

    debug!(
        "Opening output database file for writing: {:?}",
        args.output_file
    );
    // Use get_output_writer to handle potential compression
    let mut writer = get_output_writer(&args.output_file).with_context(|| {
        format!(
            "Failed to get output writer for database file: {:?}",
            args.output_file
        )
    })?;

    bincode::serialize_into(&mut writer, &kmer_db_v2).with_context(|| {
        format!(
            "Failed to serialize k-mer database (KmerDbV2) to {:?}",
            args.output_file
        )
    })?;

    // The writer from get_output_writer is already buffered (e.g. BufWriter wrapping an encoder)
    // and will flush on drop (especially ZstdEncoder with auto_finish, GzEncoder, XzEncoder).
    // Explicit flush can still be called if desired, but might be redundant.
    writer
        .flush()
        .context("Failed to flush output database writer")?;
    info!(
        "Successfully wrote k-mer database (KmerDbV2) to {:?}",
        args.output_file
    );

    Ok(())
}
