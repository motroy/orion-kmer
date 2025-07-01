use anyhow::{Context, Result};
use dashmap::DashSet; // Using DashSet to store unique k-mers concurrently
use log::{debug, info};
use needletail::{Sequence, parse_fastx_file};
use serde::{Deserialize, Serialize};
use std::{
    collections::HashSet,
    fs::File,
    io::{BufWriter, Write},
};

use crate::{
    cli::BuildArgs,
    errors::OrionKmerError,
    kmer::{canonical_u64, seq_to_u64},
};

#[derive(Serialize, Deserialize, Debug)]
pub struct KmerDb {
    pub k: u8,
    pub kmers: HashSet<u64>,
}

fn process_sequence_for_build(seq_chunk: &[u8], k: u8, kmer_set: &DashSet<u64>) {
    if seq_chunk.len() < k as usize {
        return;
    }

    for window in seq_chunk.windows(k as usize) {
        if let Some(kmer_val) = seq_to_u64(window, k) {
            let canonical_kmer = canonical_u64(kmer_val, k);
            kmer_set.insert(canonical_kmer);
        }
        // else: sequence contained 'N' or other non-ACGT char, skip this k-mer
    }
}

pub fn run_build(args: BuildArgs) -> Result<()> {
    info!("Starting build command with args: {:?}", args);

    if args.kmer_size == 0 || args.kmer_size > 32 {
        return Err(OrionKmerError::InvalidKmerSize(args.kmer_size).into());
    }
    let k = args.kmer_size;

    // DashSet is suitable for concurrent insertion of unique items.
    let kmer_set: DashSet<u64> = DashSet::new();

    for input_path in &args.genome_files {
        info!("Processing genome file: {:?}", input_path);
        let path_str = input_path.to_string_lossy();
        let mut reader = parse_fastx_file(input_path) // Removed &
            .with_context(|| format!("Failed to open or parse FASTA file: {}", path_str))?;

        info!("Processing records from {}...", path_str);
        let mut record_count = 0;
        while let Some(record) = reader.next() {
            let record =
                record.with_context(|| format!("Error reading record from {}", path_str))?;
            let norm_seq = record.normalize(false); // Ensure uppercase
            process_sequence_for_build(&norm_seq, k, &kmer_set);
            record_count += 1;
            if record_count % 100_000 == 0 {
                debug!(
                    "Processed {} records from {}. Current unique k-mers: {}",
                    record_count,
                    path_str,
                    kmer_set.len()
                );
            }
        }
        info!(
            "Finished processing {} records from {}. Current unique k-mers: {}",
            record_count,
            path_str,
            kmer_set.len()
        );
    }

    info!(
        "Finished processing all input files. Found {} unique canonical k-mers.",
        kmer_set.len()
    );

    // Convert DashSet to HashSet for serialization as KmerDb
    let final_kmers: HashSet<u64> = kmer_set.into_iter().collect();

    let kmer_db = KmerDb {
        k,
        kmers: final_kmers,
    };

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

    bincode::serialize_into(&mut writer, &kmer_db).with_context(|| {
        format!(
            "Failed to serialize k-mer database to {:?}",
            args.output_file
        )
    })?;

    writer
        .flush()
        .context("Failed to flush output database writer")?;
    info!(
        "Successfully wrote k-mer database to {:?}",
        args.output_file
    );

    Ok(())
}
