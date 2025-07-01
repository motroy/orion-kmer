use anyhow::{Context, Result};
use log::info;
use needletail::{Sequence, parse_fastx_file};
use rayon::prelude::*;
use std::{
    fs::File,
    io::{BufReader, BufWriter, Write},
    sync::Mutex, // For collecting read IDs in parallel
}; // debug was unused

use crate::{
    cli::QueryArgs,
    commands::build::KmerDb, // To load the database
    errors::OrionKmerError,
    kmer::{canonical_u64, seq_to_u64},
};

// Function to load KmerDb, similar to the one in compare.rs
// Consider moving to a shared utils module if used by more commands.
fn load_kmer_db(path: &std::path::PathBuf) -> Result<KmerDb> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open k-mer database file: {:?}", path))?;
    let reader = BufReader::new(file);
    let kmer_db: KmerDb = bincode::deserialize_from(reader)
        .with_context(|| format!("Failed to deserialize k-mer database from {:?}", path))?;
    info!(
        "Successfully loaded k-mer database from {:?} (k={}, {} kmers)",
        path,
        kmer_db.k,
        kmer_db.kmers.len()
    );
    Ok(kmer_db)
}

pub fn run_query(args: QueryArgs) -> Result<()> {
    info!("Starting query command with args: {:?}", args);

    // Load the k-mer database
    let kmer_db = load_kmer_db(&args.database_file)?;
    let k = kmer_db.k;

    if k == 0 || k > 32 {
        // This check should ideally be redundant if DB was created by this tool,
        // but good for defense if DB is from elsewhere or corrupted.
        return Err(OrionKmerError::InvalidKmerSize(k).into());
    }

    info!(
        "Querying reads from {:?} against database with k={}",
        args.reads_file, k
    );

    let mut reader = parse_fastx_file(&args.reads_file)
        .with_context(|| format!("Failed to open or parse FASTQ file: {:?}", args.reads_file))?;

    let output_file = File::create(&args.output_file).with_context(|| {
        format!(
            "Failed to create output file for matching reads: {:?}",
            args.output_file
        )
    })?;
    let writer = Mutex::new(BufWriter::new(output_file)); // Mutex for parallel writes

    // Collect records to enable parallel processing with Rayon
    // This can be memory intensive for very large FASTQ files.
    // Alternatives:
    // 1. Chunked reading and parallel processing of chunks.
    // 2. Using a library that provides a parallel FASTQ iterator.
    let mut records = Vec::new();
    while let Some(record) = reader.next() {
        let record =
            record.with_context(|| format!("Error reading record from {:?}", args.reads_file))?;
        // Store owned data (ID and sequence)
        records.push((record.id().to_vec(), record.sequence().to_owned())); // to_owned()
    }

    info!(
        "Collected {} reads. Starting parallel query...",
        records.len()
    );

    let matching_read_ids: Vec<Vec<u8>> = records
        .par_iter() // Parallel iteration over reads
        .filter_map(|(read_id_bytes, read_seq_vec)| {
            // read_seq_vec is Vec<u8>
            let mut kmer_hits = 0;
            // Use needletail::sequence::normalize which takes &[u8] and returns Vec<u8>
            // However, it's better to normalize on the fly or use an iterator if possible
            // For now, let's assume read_seq_vec is already what we want to iterate over,
            // or normalize it. record.normalize(iupac_to_n) returns an iterator.
            // The current read_seq_vec is Vec<u8> from record.sequence().to_owned().
            // seq_to_u64 expects &[u8].

            // Normalizing the sequence directly if it contains non-ACGT.
            // seq_to_u64 handles case, but not 'N'.
            // We need to ensure only valid DNA chars go to seq_to_u64, or it returns None.
            // Iterating over a normalized version is safer.
            let norm_seq: &[u8] = read_seq_vec; // It's already Vec<u8>, so slice is fine. Renamed variable.

            if norm_seq.len() < k as usize {
                // Use norm_seq here
                return None; // Read shorter than k-mer size
            }

            for window in norm_seq.windows(k as usize) {
                if let Some(kmer_val) = seq_to_u64(window, k) {
                    let canonical_kmer = canonical_u64(kmer_val, k);
                    if kmer_db.kmers.contains(&canonical_kmer) {
                        kmer_hits += 1;
                    }
                }
            }

            if kmer_hits >= args.min_hits {
                Some(read_id_bytes.clone()) // Clone read_id_bytes for collection
            } else {
                None
            }
        })
        .collect();

    info!(
        "Found {} reads matching criteria (min_hits: {}). Writing to output...",
        matching_read_ids.len(),
        args.min_hits
    );

    let mut locked_writer = writer.lock().unwrap();
    for read_id_vec in &matching_read_ids {
        // Iterate by reference
        // Read IDs in FASTQ often start with '@'. The spec doesn't say if it should be trimmed.
        // Outputting the full ID as present in the file.
        locked_writer.write_all(read_id_vec)?; // Pass &Vec<u8> which coerces to &[u8]
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
