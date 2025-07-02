use anyhow::{Context, Result};
use log::info;
use needletail::{parse_fastx_file, Sequence}; // Corrected import order
use rayon::prelude::*;
use std::{
    collections::HashSet, // Required for the unified k-mer set
    fs::File,
    io::{BufWriter, Write}, // Removed BufReader
    sync::Mutex,
};

use crate::{
    cli::QueryArgs,
    // KmerDbV2 is loaded by load_kmer_db_v2, direct import not needed here
    // db_types::KmerDbV2,
    errors::OrionKmerError,
    kmer::{canonical_u64, seq_to_u64},
    utils::{load_kmer_db_v2, track_progress_and_resources}, // Import the wrapper
};
// use indicatif::ProgressBar; // For passing to the closure - actually not needed

// Removed local load_kmer_db function

pub fn run_query(args: QueryArgs) -> Result<()> {
    info!("Starting query command with args: {:?}", args);

    // Load the KmerDbV2 database
    let kmer_db_v2 = load_kmer_db_v2(&args.database_file)?;
    let k = kmer_db_v2.k;

    if k == 0 || k > 32 {
        return Err(OrionKmerError::InvalidKmerSize(k).into());
    }

    // Get the unified set of all k-mers from the database for querying
    let db_all_kmers: HashSet<u64> = kmer_db_v2.get_all_kmers_unified();
    info!(
        "Querying reads from {:?} against database with k={} ({} unique k-mers in DB)",
        args.reads_file,
        k,
        db_all_kmers.len()
    );

    let mut reader = parse_fastx_file(&args.reads_file)
        .with_context(|| format!("Failed to open or parse FASTQ file: {:?}", args.reads_file))?;

    let output_file = File::create(&args.output_file).with_context(|| {
        format!(
            "Failed to create output file for matching reads: {:?}",
            args.output_file
        )
    })?;
    let writer = Mutex::new(BufWriter::new(output_file));

    let mut records = Vec::new();
    while let Some(record) = reader.next() {
        let record =
            record.with_context(|| format!("Error reading record from {:?}", args.reads_file))?;
        records.push((record.id().to_vec(), record.sequence().to_owned()));
    }

    info!(
        "Collected {} reads. Starting parallel query...",
        records.len()
    );

    let num_records = records.len() as u64;
    let matching_read_ids: Vec<Vec<u8>> = track_progress_and_resources(
        "Querying reads against database",
        num_records,
        |pb_query| {
            let result: Vec<Vec<u8>> = records
                .par_iter()
                .filter_map(|(read_id_bytes, read_seq_vec)| {
                    let mut kmer_hits = 0;
                    let norm_seq: &[u8] = read_seq_vec;

                    if norm_seq.len() < k as usize {
                        return None;
                    }

                    for window in norm_seq.windows(k as usize) {
                        if let Some(kmer_val) = seq_to_u64(window, k) {
                            let canonical_kmer = canonical_u64(kmer_val, k);
                            if db_all_kmers.contains(&canonical_kmer) {
                                kmer_hits += 1;
                            }
                        }
                    }

                    // Increment progress bar after processing each read
                    // Note: In a parallel iterator, direct incrementing like this might lead to
                    // frequent updates. For very large datasets, consider custom logic or
                    // relying on the automatic updates if the overhead is too high.
                    // However, indicatif is generally efficient.
                    pb_query.inc(1);


                    if kmer_hits >= args.min_hits {
                        Some(read_id_bytes.clone())
                    } else {
                        None
                    }
                })
                .collect();
            Ok(result)
        },
    )?;

    info!(
        "Found {} reads matching criteria (min_hits: {}). Writing to output...",
        matching_read_ids.len(),
        args.min_hits
    );

    let mut locked_writer = writer.lock().unwrap();
    for read_id_vec in &matching_read_ids {
        locked_writer.write_all(read_id_vec)?;
        locked_writer.write_all(b"\n")?;
    }
    locked_writer
        .flush()
        .context("Failed to flush output writer for query results")?;

    info!(
        "Successfully wrote matching read IDs to {:?}",
        args.output_file
    );

    Ok(())
}
