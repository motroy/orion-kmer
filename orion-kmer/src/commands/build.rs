use anyhow::{Context, Result};
use dashmap::DashSet; // Using DashSet for concurrent k-mer collection per file
use log::{debug, info};
use needletail::{parse_fastx_file, Sequence}; // Corrected import order
use std::{
    collections::HashSet, // Keep HashSet for final storage in KmerDbV2
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf, // For getting filename
};

use crate::{
    cli::BuildArgs,
    db_types::KmerDbV2, // Import the new database structure
    errors::OrionKmerError,
    kmer::{canonical_u64, seq_to_u64},
};

// This function processes sequences for a single file and populates a DashSet for that file.
fn process_sequences_for_file(
    file_path: &PathBuf,
    k: u8,
    file_kmer_set: &DashSet<u64>,
) -> Result<()> {
    let path_str = file_path.to_string_lossy();
    info!("Processing genome file: {}", path_str);

    let mut reader = parse_fastx_file(file_path)
        .with_context(|| format!("Failed to open or parse FASTA/Q file: {}", path_str))?;

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
            debug!(
                "Processed {} records from {}. Current unique k-mers for this file: {}",
                record_count,
                path_str,
                file_kmer_set.len()
            );
        }
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

    for input_path in &args.genome_files {
        // For each file, create a new DashSet to collect its k-mers.
        let file_kmer_set: DashSet<u64> = DashSet::new();

        process_sequences_for_file(input_path, k, &file_kmer_set)?;

        // Convert the DashSet for this file to a HashSet
        let final_file_kmers: HashSet<u64> = file_kmer_set.into_iter().collect();

        // Use the filename as the reference identifier.
        // PathBuf.file_name() returns Option<&OsStr>.
        let reference_name = input_path
            .file_name()
            .map_or_else(
                || input_path.to_string_lossy().into_owned(), // Fallback to full path if no filename
                |os_str| os_str.to_string_lossy().into_owned(),
            );

        info!(
            "Adding {} unique k-mers from reference '{}' to the database.",
            final_file_kmers.len(),
            reference_name
        );
        kmer_db_v2.add_reference(reference_name, final_file_kmers);
    }

    info!(
        "Finished processing all input files. Database contains {} references and a total of {} unique canonical k-mers across all references.",
        kmer_db_v2.num_references(),
        kmer_db_v2.total_unique_kmers()
    );

    debug!(
        "Opening output database file for writing: {:?}",
        args.output_file
    );
    let output_file = File::create(&args.output_file).with_context(|| {
        format!(
            "Failed to create output database file: {:?}",
            args.output_file
        )
    })?;
    let mut writer = BufWriter::new(output_file);

    bincode::serialize_into(&mut writer, &kmer_db_v2).with_context(|| {
        format!(
            "Failed to serialize k-mer database (KmerDbV2) to {:?}",
            args.output_file
        )
    })?;

    writer
        .flush()
        .context("Failed to flush output database writer")?;
    info!(
        "Successfully wrote k-mer database (KmerDbV2) to {:?}",
        args.output_file
    );

    Ok(())
}
